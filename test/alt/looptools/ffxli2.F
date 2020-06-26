#include "externals.h"
#include "types.h"


*###[ ffxli2:
	subroutine ffxli2(xdilog,xlog,x,ier)
***#[*comment:***********************************************************
*									*
*	Computes the dilogarithm (Li2, Sp) for (real) x	to precision	*
*	precx. It is assumed that -1<=x<=1/2.  As it is available anyway*
*	log(1-x) = -Li1(x) is also passed.				*
*									*
*	Input:	x	(real)						*
*									*
*	Output: xdilog	(real)    Li2(x)				*
*		xlog	(real)    log(1-x) = -Li1(x)			*
*		ier	(integer) 0=OK, 1=num prob, 2=error		*
*									*
*	Calls:	log,dfflo1						*
*									*
***#]*comment:***********************************************************
*  #[ declarations:
	implicit none
*
*	arguments
*
	integer ier
	RealType xdilog,xlog,x
*
*	local variables
*
	integer ipi12
	RealType dfflo1,u,u2,a,ffbnd,
     +		xprec,bdn02,bdn05,bdn10,bdn15
	ComplexType zxdilo,zlog
	external ffbnd,dfflo1
	save xprec,bdn02,bdn05,bdn10,bdn15
*
*	common blocks
*
#include "ff.h"
*  #] declarations:
*  #[ initialisations:
	data xprec /-1D0/
	if ( xprec .ne. precx ) then
	    xprec = precx
	    bdn02 = ffbnd(1,2,bf)
	    bdn05 = ffbnd(1,5,bf)
	    bdn10 = ffbnd(1,10,bf)
	    bdn15 = ffbnd(1,15,bf)
	endif
*  #] initialisations:
*  #[ if the argument is too large...
	if ( x .lt. -1.5 .or. x .gt. .75 ) then
	    call ffzxdl(zxdilo,ipi12,zlog,x,0,ier)
	    if ( Im(zxdilo) .ne. 0 ) then
	    	call fferr(52,ier)
	    endif
	    xdilog = Re(zxdilo) + ipi12*pi12
	    xlog = Re(zlog)
	    return
	endif
*  #] if the argument is too large...
*  #[ exceptional cases:
	if ( x .eq. -1 ) then
	    xdilog = -pi12
	    xlog = log(2D0)
	    return
	elseif ( x .eq. .5D0 ) then
	    xdilog = - xlg2**2/2 + pi12
	    xlog = - xlg2
	    return
	elseif ( abs(x) .lt. precx ) then
	    xdilog = x
	    xlog = -x
	    return
	endif
*  #] exceptional cases:
*  #[ calculate dilog:
	if ( abs(x) .lt. xloss ) then
	    xlog = dfflo1(x,ier)
	else
	    xlog = log(1-x)
	endif
	u = -xlog
	u2 = u*u
	a = abs(u2)
	if ( a .gt. bdn15 ) then
	    xdilog = u2*(bf(16) + u2*(bf(17) + u2*(bf(18) +
     +		    u2*(bf(19) + u2*bf(20) ))))
	else
	    xdilog = 0
	endif
	if ( a .gt. bdn10 ) then
	    xdilog = u2*(bf(11) + u2*(bf(12) + u2*(bf(13) +
     +		    u2*(bf(14) + u2*(bf(15) + xdilog)))))
	endif
	if ( a .gt. bdn05 ) then
	    xdilog = u2*(bf(6) + u2*(bf(7) + u2*(bf(8) +
     +		    u2*(bf(9) + u2*(bf(10) + xdilog)))))
	endif
	if ( a .gt. bdn02 ) then
	    xdilog = u2*(bf(3) + u2*(bf(4) + u2*(bf(5) + xdilog)))
	endif
*	watch the powers of u.
	xdilog = u + u2*(bf(1) + u*(bf(2) + xdilog))
*  #] calculate dilog:
*###] ffxli2:
	end


*###[ ffzxdl:
	subroutine ffzxdl(zxdilo,ipi12,zlog,x,ieps,ier)
***#[*comment:***************************************************
*	Computes the dilogarithm (Li2, Sp) for any (real) x	*
*	to precision precx. If an error message is given add	*
*	more bf's. For x > 1 the imaginary part is		*
*	 -/+i*pi*log(x), corresponding to x+ieps.		*
*	The number of factors pi^2/12 is passed separately in	*
*	ipi12 for accuracy.  We also calculate log(1-x)		*
*	which is likely to be needed.				*
*								*
*	Input:	x	(real)					*
*		ieps	(integer,+/-1)				*
*								*
*	Output: zxdilo	(complex) the dilog mod factors pi2/12	*
*		ipi12	(integer) these factors			*
*		zlog	(complex) log(1-x)			*
*								*
*	Calls:	log,dfflo1					*
*								*
***#]*comment:***************************************************
*  #[ declarations:
	implicit none
*
*	arguments
*
	integer ipi12,ieps,ier
	RealType x
	ComplexType zxdilo,zlog
*
*	local variables
*
	integer jsgn
	RealType fact,u,u2,dfflo1,ffbnd,a,xdilo,
     +		xprec,bdn02,bdn05,bdn10,bdn15
	ComplexType cy,cfact
	external ffbnd,dfflo1
	save xprec,bdn02,bdn05,bdn10,bdn15
*
*	common blocks
*
#include "ff.h"
*  #] declarations:
*  #[ initialisations:
	data xprec /-1D0/
	if ( xprec .ne. precx ) then
	    xprec = precx
	    bdn02 = ffbnd(1,2,bf)
	    bdn05 = ffbnd(1,5,bf)
	    bdn10 = ffbnd(1,10,bf)
	    bdn15 = ffbnd(1,15,bf)
	endif
*  #] initialisations:
*  #[ exceptional cases:
	if ( x .eq. 1) then
	    zxdilo = 0
	    zlog = -99999
	    ipi12 = 2
	    return
	elseif (x .eq. -1) then
	    zxdilo = 0
	    zlog = xlg2
	    ipi12 = -1
	    return
	elseif (x .eq. .5D0) then
	    zxdilo = - xlg2**2/2
	    zlog = -xlg2
	    ipi12 = 1
	    return
	elseif ( abs(x) .lt. precx ) then
	    zxdilo = x
	    zlog = -x
	    ipi12 = 0
	    return
	endif
*  #] exceptional cases:
*  #[ transform to (-1,.5):
	if (x .lt. -1) then
	    fact = log(-x)
	    cy = - fact**2/2
	    ipi12 = -2
	    if ( -x*xloss .gt. 1 ) then
		u = -dfflo1(1/x,ier)
	    else
		u = -log(1-1/x)
	    endif
	    zlog = log(1-x)
	    jsgn = -1
	elseif ( x .lt. .5D0) then
	    cy = 0
	    ipi12 = 0
	    if ( abs(x) .lt. xloss ) then
		zlog = dfflo1(x,ier)
	    else
		zlog = log(1-x)
	    endif
	    u = -Re(zlog)
	    jsgn = 1
	elseif ( x .le. 2 ) then
	    u = -log(x)
	    if ( abs(1-x) .lt. xalogm ) then
		cy = 0
	    elseif ( x .lt. 1 ) then
		zlog = log(1-x)
		cy = Re(u)*zlog
	    elseif ( ieps .gt. 0 ) then
		zlog = ToComplex(log(x-1),-pi)
		cy = Re(u)*zlog
	    else
		zlog = ToComplex(log(x-1),+pi)
		cy = Re(u)*zlog
	    endif
	    ipi12 = 2
	    jsgn = -1
	else
	    if ( ieps .gt. 0 ) then
		cfact = ToComplex(log(x),-pi)
		zlog = ToComplex(log(x-1),-pi)
	    else
		cfact = ToComplex(log(x),+pi)
		zlog = ToComplex(log(x-1),+pi)
	    endif
	    cy = - cfact**2/2
	    ipi12 = -2
	    if ( x*xloss .gt. 1 ) then
		u = -dfflo1(1/x,ier)
	    else
		u = -log(1-1/x)
	    endif
	    jsgn = -1
	endif
*  #] transform to (-1,.5):
*  #[ calculate dilog:
	if ( abs(u) .lt. xalog2 ) then
	    xdilo = u
	else
	u2 = u**2
	a = abs(u2)
	if ( a .gt. bdn15 ) then
	    xdilo = u2*(bf(16) + u2*(bf(17) + u2*(bf(18) +
     +		    u2*(bf(19) + u2*bf(20) ))))
	else
	    xdilo = 0
	endif
	if ( a .gt. bdn10 ) then
	    xdilo = u2*(bf(11) + u2*(bf(12) + u2*(bf(13) +
     +		    u2*(bf(14) + u2*(bf(15) + xdilo)))))
	endif
	if ( a .gt. bdn05 ) then
	    xdilo = u2*(bf(6) + u2*(bf(7) + u2*(bf(8) +
     +		    u2*(bf(9) + u2*(bf(10) + xdilo)))))
	endif
	if ( a .gt. bdn02 ) then
	    xdilo = u2*(bf(3) + u2*(bf(4) + u2*(bf(5) + xdilo)))
	endif
*	watch the powers of u.
	xdilo = u + u2*(bf(1) + u*(bf(2) + xdilo))
	endif
	if(jsgn.eq.1)then
	    zxdilo =  Re(xdilo) + cy
	else
	    zxdilo = -Re(xdilo) + cy
	endif
*  #] calculate dilog:
*###] ffzxdl:
	end


*###[ zxfflg:
	ComplexType function zxfflg(x,ieps,y,ier)
***#[*comment:***********************************************************
*									*
*	Calculate the complex logarithm of x.  The following cases	*
*	are treted separately:						*
*		|x| too small:		give warning and return 0	*
*					(for Absoft, Apollo DN300)	*
*		|x| < 0:		take sign according to ieps	*
*									*
***#]*comment:***********************************************************
*  #[ declarations:
*
*	arguments
*
	implicit none
	integer ieps,ier
	RealType x,y
*
*	local variables
*
	RealType xlog
*
*	common blocks
*
#include "ff.h"
*  #] declarations:
*  #[ calculations:
	if ( abs(x) .lt. xalogm ) then
	    zxfflg = 0
	elseif ( x .gt. 0 ) then
	    zxfflg = log(x)
	else
	    xlog = log(-x)
*	    checked imaginary parts 19-May-1988
	    if ( abs(ieps) .eq. 1 ) then
		if ( y*ieps .lt. 0 ) then
		    zxfflg = ToComplex(xlog,-pi)
		else
		    zxfflg = ToComplex(xlog,pi)
		endif
	    elseif ( ieps .eq. 2 ) then
		zxfflg = ToComplex(xlog,-pi)
	    elseif ( ieps .eq. -2 ) then
		zxfflg = ToComplex(xlog,+pi)
	    else
		call fferr(52,ier)
		zxfflg = ToComplex(xlog,pi)
	    endif
	endif
*  #] calculations:
*###] zxfflg:
	end


*###[ dfflo1:
	RealType function dfflo1(x,ier)
***#[*comment:***************************************************
*	calculates log(1-x) for |x|<.14 in a faster way to ~15	*
*	significant figures.					*
***#]*comment:***************************************************
*  #[ declarations:
	implicit none
	integer ier
	RealType x,bdn01,bdn05,bdn10,bdn15,bdn19,xprec,
     +		xa,ffbnd
	ComplexType zxfflg
	external ffbnd,zxfflg
	save xprec,bdn01,bdn05,bdn10,bdn15,bdn19
#include "ff.h"
*  #] declarations:
*  #[ initialisation:
	data xprec /-1D0/
	if ( xprec .ne. precx ) then
	    xprec = precx
*	    determine the boundaries for 1,5,10,15 terms
	    bdn01 = ffbnd(1,1,xninv)
	    bdn05 = ffbnd(1,5,xninv)
	    bdn10 = ffbnd(1,10,xninv)
	    bdn15 = ffbnd(1,15,xninv)
	    bdn19 = ffbnd(1,19,xninv)
	endif
*  #] initialisation:
*  #[ calculations:
	xa = abs(x)
	if ( xa .gt. bdn19 ) then
	    dfflo1 = Re(zxfflg(1-x,0,0D0,ier))
	    return
	endif
	if ( xa .gt. bdn15 ) then
	    dfflo1 = x*( xninv(16) + x*( xninv(17) + x*( xninv(18) +
     +		x*( xninv(19) + x*xninv(20) ))))
	else
	    dfflo1 = 0
	endif
	if ( xa .gt. bdn10 ) then
	    dfflo1 = x*( xninv(11) + x*( xninv(12) + x*( xninv(13) +
     +		x*( xninv(14) + x*( xninv(15) + dfflo1 )))))
	endif
	if ( xa .gt. bdn05 ) then
	    dfflo1 = x*( xninv(6) + x*( xninv(7) + x*( xninv(8) +
     +		x*( xninv(9) + x*( xninv(10) + dfflo1 )))))
	endif
	if ( xa .gt. bdn01 ) then
	    dfflo1 = x*( xninv(2) + x*( xninv(3) + x*( xninv(4) +
     +		x*( xninv(5) + dfflo1 ))))
	endif
	dfflo1 = - x*( xninv(1) + dfflo1 )
*  #] calculations:
*###] dfflo1:
	end


*###[ dfflo2:
	RealType function dfflo2(x,ier)
***#[*comment:***************************************************
*	calculates log(1-x)+x for |x|<.14 in a faster way to	*
*	~15 significant figures.				*
***#]*comment:***************************************************
*  #[ declarations:
	implicit none
	integer ier
	RealType x,bdn01,bdn05,bdn10,bdn15,bdn18,xprec,
     +		xa,ffbnd,dfflo1
	external ffbnd,dfflo1
	save xprec,bdn01,bdn05,bdn10,bdn15,bdn18
#include "ff.h"
*  #] declarations:
*  #[ initialisation:
	data xprec /-1D0/
	if ( xprec .ne. precx ) then
	    xprec = precx
*	    determine the boundaries for 1,5,10,15 terms
	    bdn01 = ffbnd(1,1,xninv(2))
	    bdn05 = ffbnd(1,5,xninv(2))
	    bdn10 = ffbnd(1,10,xninv(2))
	    bdn15 = ffbnd(1,15,xninv(2))
	    bdn18 = ffbnd(1,18,xninv(2))
	endif
*  #] initialisation:
*  #[ calculations:
	xa = abs(x)
	if ( xa .gt. bdn18 ) then
	    dfflo2 = dfflo1(x,ier) + x
	    return
	endif
	if ( xa .gt. bdn15 ) then
	    dfflo2 = x*( xninv(17) + x*( xninv(18) + x*( xninv(19) +
     +		x*xninv(20) )))
	else
	    dfflo2 = 0
	endif
	if ( xa .gt. bdn10 ) then
	    dfflo2 = x*( xninv(12) + x*( xninv(13) + x*( xninv(14) +
     +		x*( xninv(15) + x*( xninv(16) + dfflo2 )))))
	endif
	if ( xa .gt. bdn05 ) then
	    dfflo2 = x*( xninv(7) + x*( xninv(8) + x*( xninv(9) +
     +		x*( xninv(10) + x*( xninv(11) + dfflo2 )))))
	endif
	if ( xa .gt. bdn01 ) then
	    dfflo2 = x*( xninv(3) + x*( xninv(4) + x*( xninv(5) +
     +		x*( xninv(6) + dfflo2 ))))
	endif
	dfflo2 = - x**2*( xninv(2) + dfflo2 )
*  #] calculations:
*###] dfflo2:
	end


*###[ dfflo3:
	RealType function dfflo3(x,ier)
***#[*comment:***************************************************
*	calculates log(1-x)+x+x^2/2 for |x|<.14 in a faster 	*
*	way to ~15 significant figures.				*
***#]*comment:***************************************************
*  #[ declarations:
	implicit none
	integer ier
	RealType x,bdn01,bdn05,bdn10,bdn15,xprec,
     +		xa,ffbnd,dfflo2
	external ffbnd,dfflo2
	save xprec,bdn01,bdn05,bdn10,bdn15
#include "ff.h"
*  #] declarations:
*  #[ initialisation:
	data xprec /-1D0/
	if ( xprec .ne. precx ) then
	    xprec = precx
*	    determine the boundaries for 1,5,10,15 terms
	    bdn01 = ffbnd(1,1,xninv(3))
	    bdn05 = ffbnd(1,5,xninv(3))
	    bdn10 = ffbnd(1,10,xninv(3))
	    bdn15 = ffbnd(1,15,xninv(3))
	endif
*  #] initialisation:
*  #[ calculations:
	xa = abs(x)
	if ( xa .gt. bdn15 ) then
	    dfflo3 = dfflo2(x,ier) + x**2/2
	    return
	endif
	if ( xa .gt. bdn10 ) then
	    dfflo3 = x*( xninv(13) + x*( xninv(14) + x*( xninv(15) +
     +		x*( xninv(16) + x*xninv(17) ))))
	else
	    dfflo3 = 0
	endif
	if ( xa .gt. bdn05 ) then
	    dfflo3 = x*( xninv(8) + x*( xninv(9) + x*( xninv(10) +
     +		x*( xninv(11) + x*( xninv(12) + dfflo3 )))))
	endif
	if ( xa .gt. bdn01 ) then
	    dfflo3 = x*( xninv(4) + x*( xninv(5) + x*( xninv(6) +
     +		x*( xninv(7) + dfflo3 ))))
	endif
	dfflo3 = - x**3*( xninv(3) + dfflo3 )
*  #] calculations:
*###] dfflo3:
	end


*###[ ffxl22:
	subroutine ffxl22(xl22,x,ier)
***#[*comment:***************************************************
*	calculates Li2(2-x) for |x|<.14 in a faster way to ~15	*
*	significant figures.					*
***#]*comment:***************************************************
*  #[ declarations:
	implicit none
	integer ier,init
	RealType xl22,x,bdn01,bdn05,bdn10,bdn15,bdn20,bdn25,
     +		xprec,xa,ffbnd,dilog2(29)
	external ffbnd
	save xprec,bdn01,bdn05,bdn10,bdn15,bdn20,bdn25,init,dilog2
#include "ff.h"
	data xprec /-1D0/
	data init /0/
	if ( init .eq. 0 ) then
	    init = 1
* taylor(dilog(x-1),x,30);
	    dilog2( 1) = 0.d0
	    dilog2( 2) = 1/4.d0
	    dilog2( 3) = 1/6.d0
	    dilog2( 4) = 5/48.d0
	    dilog2( 5) = 1/15.d0
	    dilog2( 6) = 2/45.d0
	    dilog2( 7) = 13/420.d0
	    dilog2( 8) = 151/6720.d0
	    dilog2( 9) = 16/945.d0
	    dilog2(10) = 83/6300.d0
	    dilog2(11) = 73/6930.d0
	    dilog2(12) = 1433/166320.d0
	    dilog2(13) = 647/90090.d0
	    dilog2(14) = 15341/2522520.d0
	    dilog2(15) = 28211/5405400.d0
	    dilog2(16) = 10447/2306304.d0
	    dilog2(17) = 608/153153.d0
	    dilog2(18) = 19345/5513508.d0
	    dilog2(19) = 18181/5819814.d0
	    dilog2(20) = 130349/46558512.d0
	    dilog2(21) = 771079/305540235.d0
	    dilog2(22) = 731957/320089770.d0
	    dilog2(23) = 2786599/1338557220.d0
	    dilog2(24) = 122289917/64250746560.d0
	    dilog2(25) = 14614772/8365982625.d0
	    dilog2(26) = 140001721/87006219300.d0
	    dilog2(27) = 134354573/90352612350.d0
	    dilog2(28) = 774885169/562194032400.d0
	    dilog2(29) = 745984697/582272390700.d0
	endif
*  #] declarations:
*  #[ initialisation:
	if ( xprec .ne. precx ) then
	    xprec = precx
*	    determine the boundaries for 1,5,10,15,20 terms
	    bdn01 = ffbnd(2,1,dilog2)
	    bdn05 = ffbnd(2,5,dilog2)
	    bdn10 = ffbnd(2,10,dilog2)
	    bdn15 = ffbnd(2,15,dilog2)
	    bdn20 = ffbnd(2,20,dilog2)
	    bdn25 = ffbnd(2,25,dilog2)
	endif
*  #] initialisation:
*  #[ calculations:
	xa = abs(x)
	if ( xa .gt. bdn25 ) then
	    call ffwarn(230,ier,precx,dilog2(27)*xa**25)
	endif
	if ( xa .gt. bdn20 ) then
	    xl22 = x*( dilog2(22) + x*( dilog2(23) + x*( dilog2(24) +
     +		x*( dilog2(25) + x*dilog2(26) ))))
	else
	    xl22 = 0
	endif
	if ( xa .gt. bdn15 ) then
	    xl22 = x*( dilog2(17) + x*( dilog2(18) + x*( dilog2(19) +
     +		x*( dilog2(20) + x*dilog2(21) ))))
	endif
	if ( xa .gt. bdn10 ) then
	    xl22 = x*( dilog2(12) + x*( dilog2(13) + x*( dilog2(14) +
     +		x*( dilog2(15) + x*dilog2(16) ))))
	endif
	if ( xa .gt. bdn05 ) then
	    xl22 = x*( dilog2(7) + x*( dilog2(8) + x*( dilog2(9) +
     +		x*( dilog2(10) + x*( dilog2(11) + xl22 )))))
	endif
	if ( xa .gt. bdn01 ) then
	    xl22 = x*( dilog2(3) + x*( dilog2(4) + x*( dilog2(5) +
     +		x*( dilog2(6) + xl22 ))))
	endif
	xl22 = - x**2*( dilog2(2) + xl22 )
*  #] calculations:
*###] ffxl22:
	end
