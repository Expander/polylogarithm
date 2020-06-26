* ini.F
* routines for initializing and setting some parameters
* this file is part of LoopTools
* last modified 14 Apr 18 th

#include "externals.h"
#include "types.h"
#include "defs.h"


	subroutine clearcache
	implicit none

#include "lt.h"

	integer i

	do i = 1, ncaches
	  cacheptr(1,0,i) = 0
	  cacheptr(2,0,i) = 0
	  savedptr(1,i) = 0
	  savedptr(2,i) = 0
	enddo
	end

************************************************************************

	subroutine markcache
	implicit none

#include "lt.h"

	integer i

	do i = 1, ncaches
	  savedptr(1,i) = cacheptr(1,0,i)
	  savedptr(2,i) = cacheptr(2,0,i)
	enddo
	end

************************************************************************

	subroutine restorecache
	implicit none

#include "lt.h"

	integer i

	do i = 1, ncaches
	  cacheptr(1,0,i) = savedptr(1,i)
	  cacheptr(2,0,i) = savedptr(2,i)
	enddo
	end

************************************************************************
* Legacy function, provided for compatibility only.
* Works only approximately as before!

	subroutine setcachelast(base, offset)
	implicit none
	ComplexType base(*)
	integer offset

	logical ini
	data ini /.TRUE./

	if( ini ) then
	  print *, "setcachelast is deprecated"
	  print *, "use clearcache or restorecache instead"
	  ini = .FALSE.
	endif

	if( offset .eq. 0 ) then
	  call clearcache
	else
	  call restorecache
	endif
	end

************************************************************************
* Legacy function, provided for compatibility only.
* Works only approximately as before!

	integer function getcachelast(base)
	implicit none
	ComplexType base(*)

	logical ini
	data ini /.TRUE./

	if( ini ) then
	  print *, "getcachelast is deprecated"
	  print *, "use markcache instead"
	  ini = .FALSE.
	endif

	getcachelast = 1
	call markcache
	end

************************************************************************

	subroutine setmudim(mudim_)
	implicit none
	RealType mudim_

#include "lt.h"

	if( abs(mudim - mudim_) .gt. diffeps ) call clearcache
	mudim = mudim_
	end

************************************************************************

	RealType function getmudim()
	implicit none

#include "lt.h"

	getmudim = mudim
	end

************************************************************************

	subroutine setdelta(delta_)
	implicit none
	RealType delta_

#include "lt.h"

	if( abs(delta - delta_) .gt. diffeps ) call clearcache
	delta = delta_
	end

************************************************************************

	RealType function getdelta()
	implicit none

#include "lt.h"

	getdelta = delta
	end

************************************************************************

	subroutine setuvdiv(uvdiv_)
	implicit none
	RealType uvdiv_

#include "lt.h"

	if( abs(uvdiv - uvdiv_) .gt. diffeps ) call clearcache
	uvdiv = uvdiv_
	end

************************************************************************

	RealType function getuvdiv()
	implicit none

#include "lt.h"

	getuvdiv = uvdiv
	end

************************************************************************

	subroutine setlambda(lam_)
	implicit none
	RealType lam_

#include "lt.h"

	RealType lambda_

	if( lam_ .ge. 0 .or. lam_ .eq. -1 .or. lam_ .eq. -2 ) then
	  lambda_ = dim(lam_, 0D0)
	  epsi = int(dim(0D0, lam_))
	else
	  print *, "illegal value for lambda"
	  lambda_ = 0
	  epsi = 0
	endif
	if( abs(lambda - lambda_) .gt. diffeps ) call clearcache
	lambda = lambda_
	end

************************************************************************

	RealType function getlambda()
	implicit none

#include "lt.h"

	getlambda = lambda
	end

************************************************************************

	integer function getepsi()
	implicit none

#include "lt.h"

	getepsi = epsi
	end

************************************************************************

	subroutine setminmass(minmass_)
	implicit none
	RealType minmass_

#include "lt.h"

	if( abs(minmass - minmass_) .gt. diffeps ) call clearcache
	minmass = minmass_
	end

************************************************************************

	RealType function getminmass()
	implicit none

#include "lt.h"

	getminmass = minmass
	end

************************************************************************

	subroutine setmaxdev(maxdev_)
	implicit none
	RealType maxdev_

#include "lt.h"

	maxdev = maxdev_
	end

************************************************************************

	RealType function getmaxdev()
	implicit none

#include "lt.h"

	getmaxdev = maxdev
	end

************************************************************************

	subroutine setdiffeps(diffeps_)
	implicit none
	RealType diffeps_

#include "lt.h"

	diffeps = diffeps_
	end

************************************************************************

	RealType function getdiffeps()
	implicit none

#include "lt.h"

	getdiffeps = diffeps
	end

************************************************************************

	subroutine setzeroeps(zeroeps_)
	implicit none
	RealType zeroeps_

#include "lt.h"

	zeroeps = zeroeps_
	end

************************************************************************

	RealType function getzeroeps()
	implicit none

#include "lt.h"

	getzeroeps = zeroeps
	end

************************************************************************

	subroutine setwarndigits(warndigits_)
	implicit none
	integer warndigits_

#include "lt.h"

	warndigits = warndigits_
	end

************************************************************************

	integer function getwarndigits()
	implicit none

#include "lt.h"

	getwarndigits = warndigits
	end

************************************************************************

	subroutine seterrdigits(errdigits_)
	implicit none
	integer errdigits_

#include "lt.h"

	errdigits = errdigits_
	end

************************************************************************

	integer function geterrdigits()
	implicit none

#include "lt.h"

	geterrdigits = errdigits
	end

************************************************************************

	subroutine setversionkey(versionkey_)
	implicit none
	integer versionkey_

#include "lt.h"

	versionkey = versionkey_
	call clearcache
	end

************************************************************************

	integer function getversionkey()
	implicit none

#include "lt.h"

	getversionkey = versionkey
	end

************************************************************************

	subroutine setdebugkey(debugkey_)
	implicit none
	integer debugkey_

#include "lt.h"

	debugkey = debugkey_
	end

************************************************************************

	integer function getdebugkey()
	implicit none

#include "lt.h"

	getdebugkey = debugkey
	end

************************************************************************

	subroutine setdebugrange(debugfrom_, debugto_)
	implicit none
	integer debugfrom_, debugto_

#include "lt.h"

	debugfrom = debugfrom_
	debugto = debugto_
	end

************************************************************************

	subroutine setcmpbits(cmpbits_)
	implicit none
	integer cmpbits_

#include "lt.h"

	cmpbits = max(cmpbits_, 12 + QUAD*4)
	end

************************************************************************

	integer function getcmpbits()
	implicit none

#include "lt.h"

	getcmpbits = cmpbits
	end

************************************************************************

* This silly subroutine is called from ffini while determining
* the working precision of the machine we're running on.
* It works around the optimizer to guarantee that we're not in
* fact determining the precision of the FPU registers.

	subroutine ffset(res, x)
	implicit none
	RealType res, x

	res = x
	end

************************************************************************

	block data namedata
	implicit none

#include "ltnames.h"

	data paraname(1:Paa,1) /
     &    "m" /

	data paraname(1:Pbb,2) /
     &    "m1", "m2", "p" /

	data paraname(1:Pcc,3) /
     &    "m1", "m2", "m3", "p1", "p2", "p1p2" /

	data paraname(1:Pdd,4) /
     &    "m1", "m2", "m3", "m4",
     &    "p1", "p2", "p3", "p4", "p1p2", "p2p3" /

	data paraname(1:Pee,5) /
     &    "m1", "m2", "m3", "m4", "m5",
     &    "p1", "p2", "p3", "p4", "p5",
     &    "p1p2", "p2p3", "p3p4", "p4p5", "p5p1" /

	data coeffname(1:Naa,1) /
     &    "aa0", "aa0:1", "aa0:2",
     &    "aa00", "aa00:1", "aa00:2" /

	data coeffname(1:Nbb,2) /
     &    "bb0", "bb0:1", "bb0:2",
     &    "bb1", "bb1:1", "bb1:2",
     &    "bb00", "bb00:1", "bb00:2",
     &    "bb11", "bb11:1", "bb11:2",
     &    "bb001", "bb001:1", "bb001:2",
     &    "bb111", "bb111:1", "bb111:2",
     &    "dbb0", "dbb0:1", "dbb0:2",
     &    "dbb1", "dbb1:1", "dbb1:2",
     &    "dbb00", "dbb00:1", "dbb00:2",
     &    "dbb11", "dbb11:1", "dbb11:2",
     &    "dbb001", "dbb001:1", "dbb001:2" /

	data coeffname(1:Ncc,3) /
     &    "cc0", "cc0:1", "cc0:2",
     &    "cc1", "cc1:1", "cc1:2",
     &    "cc2", "cc2:1", "cc2:2",
     &    "cc00", "cc00:1", "cc00:2",
     &    "cc11", "cc11:1", "cc11:2",
     &    "cc12", "cc12:1", "cc12:2",
     &    "cc22", "cc22:1", "cc22:2",
     &    "cc001", "cc001:1", "cc001:2",
     &    "cc002", "cc002:1", "cc002:2",
     &    "cc111", "cc111:1", "cc111:2",
     &    "cc112", "cc112:1", "cc112:2",
     &    "cc122", "cc122:1", "cc122:2",
     &    "cc222", "cc222:1", "cc222:2",
     &    "cc0000", "cc0000:1", "cc0000:2",
     &    "cc0011", "cc0011:1", "cc0011:2",
     &    "cc0012", "cc0012:1", "cc0012:2",
     &    "cc0022", "cc0022:1", "cc0022:2",
     &    "cc1111", "cc1111:1", "cc1111:2",
     &    "cc1112", "cc1112:1", "cc1112:2",
     &    "cc1122", "cc1122:1", "cc1122:2",
     &    "cc1222", "cc1222:1", "cc1222:2",
     &    "cc2222", "cc2222:1", "cc2222:2" /

	data coeffname(1:Ndd,4) /
     &    "dd0", "dd0:1", "dd0:2",
     &    "dd1", "dd1:1", "dd1:2",
     &    "dd2", "dd2:1", "dd2:2",
     &    "dd3", "dd3:1", "dd3:2",
     &    "dd00", "dd00:1", "dd00:2",
     &    "dd11", "dd11:1", "dd11:2",
     &    "dd12", "dd12:1", "dd12:2",
     &    "dd13", "dd13:1", "dd13:2",
     &    "dd22", "dd22:1", "dd22:2",
     &    "dd23", "dd23:1", "dd23:2",
     &    "dd33", "dd33:1", "dd33:2",
     &    "dd001", "dd001:1", "dd001:2",
     &    "dd002", "dd002:1", "dd002:2",
     &    "dd003", "dd003:1", "dd003:2",
     &    "dd111", "dd111:1", "dd111:2", 
     &    "dd112", "dd112:1", "dd112:2",
     &    "dd113", "dd113:1", "dd113:2",
     &    "dd122", "dd122:1", "dd122:2",
     &    "dd123", "dd123:1", "dd123:2",
     &    "dd133", "dd133:1", "dd133:2",
     &    "dd222", "dd222:1", "dd222:2",
     &    "dd223", "dd223:1", "dd223:2",
     &    "dd233", "dd233:1", "dd233:2",
     &    "dd333", "dd333:1", "dd333:2",
     &    "dd0000", "dd0000:1", "dd0000:2",
     &    "dd0011", "dd0011:1", "dd0011:2",
     &    "dd0012", "dd0012:1", "dd0012:2",
     &    "dd0013", "dd0013:1", "dd0013:2",
     &    "dd0022", "dd0022:1", "dd0022:2",
     &    "dd0023", "dd0023:1", "dd0023:2",
     &    "dd0033", "dd0033:1", "dd0033:2",
     &    "dd1111", "dd1111:1", "dd1111:2",
     &    "dd1112", "dd1112:1", "dd1112:2",
     &    "dd1113", "dd1113:1", "dd1113:2",
     &    "dd1122", "dd1122:1", "dd1122:2",
     &    "dd1123", "dd1123:1", "dd1123:2",
     &    "dd1133", "dd1133:1", "dd1133:2",
     &    "dd1222", "dd1222:1", "dd1222:2",
     &    "dd1223", "dd1223:1", "dd1223:2",
     &    "dd1233", "dd1233:1", "dd1233:2",
     &    "dd1333", "dd1333:1", "dd1333:2",
     &    "dd2222", "dd2222:1", "dd2222:2",
     &    "dd2223", "dd2223:1", "dd2223:2",
     &    "dd2233", "dd2233:1", "dd2233:2",
     &    "dd2333", "dd2333:1", "dd2333:2",
     &    "dd3333", "dd3333:1", "dd3333:2",
     &    "dd00001", "dd00001:1", "dd00001:2",
     &    "dd00002", "dd00002:1", "dd00002:2",
     &    "dd00003", "dd00003:1", "dd00003:2",
     &    "dd00111", "dd00111:1", "dd00111:2",
     &    "dd00112", "dd00112:1", "dd00112:2",
     &    "dd00113", "dd00113:1", "dd00113:2",
     &    "dd00122", "dd00122:1", "dd00122:2",
     &    "dd00123", "dd00123:1", "dd00123:2",
     &    "dd00133", "dd00133:1", "dd00133:2",
     &    "dd00222", "dd00222:1", "dd00222:2",
     &    "dd00223", "dd00223:1", "dd00223:2",
     &    "dd00233", "dd00233:1", "dd00233:2",
     &    "dd00333", "dd00333:1", "dd00333:2",
     &    "dd11111", "dd11111:1", "dd11111:2", 
     &    "dd11112", "dd11112:1", "dd11112:2",
     &    "dd11113", "dd11113:1", "dd11113:2",
     &    "dd11122", "dd11122:1", "dd11122:2",
     &    "dd11123", "dd11123:1", "dd11123:2",
     &    "dd11133", "dd11133:1", "dd11133:2",
     &    "dd11222", "dd11222:1", "dd11222:2",
     &    "dd11223", "dd11223:1", "dd11223:2",
     &    "dd11233", "dd11233:1", "dd11233:2",
     &    "dd11333", "dd11333:1", "dd11333:2",
     &    "dd12222", "dd12222:1", "dd12222:2",
     &    "dd12223", "dd12223:1", "dd12223:2",
     &    "dd12233", "dd12233:1", "dd12233:2",
     &    "dd12333", "dd12333:1", "dd12333:2",
     &    "dd13333", "dd13333:1", "dd13333:2",
     &    "dd22222", "dd22222:1", "dd22222:2",
     &    "dd22223", "dd22223:1", "dd22223:2",
     &    "dd22233", "dd22233:1", "dd22233:2",
     &    "dd22333", "dd22333:1", "dd22333:2",
     &    "dd23333", "dd23333:1", "dd23333:2",
     &    "dd33333", "dd33333:1", "dd33333:2" /

	data coeffname(1:Nee,5) /
     &    "ee0", "ee0:1", "ee0:2",
     &    "ee1", "ee1:1", "ee1:2",
     &    "ee2", "ee2:1", "ee2:2",
     &    "ee3", "ee3:1", "ee3:2",
     &    "ee4", "ee4:1", "ee4:2",
     &    "ee00", "ee00:1", "ee00:2",
     &    "ee11", "ee11:1", "ee11:2",
     &    "ee12", "ee12:1", "ee12:2",
     &    "ee13", "ee13:1", "ee13:2",
     &    "ee14", "ee14:1", "ee14:2",
     &    "ee22", "ee22:1", "ee22:2",
     &    "ee23", "ee23:1", "ee23:2",
     &    "ee24", "ee24:1", "ee24:2",
     &    "ee33", "ee33:1", "ee33:2",
     &    "ee34", "ee34:1", "ee34:2",
     &    "ee44", "ee44:1", "ee44:2",
     &    "ee001", "ee001:1", "ee001:2",
     &    "ee002", "ee002:1", "ee002:2",
     &    "ee003", "ee003:1", "ee003:2",
     &    "ee004", "ee004:1", "ee004:2",
     &    "ee111", "ee111:1", "ee111:2",
     &    "ee112", "ee112:1", "ee112:2",
     &    "ee113", "ee113:1", "ee113:2",
     &    "ee114", "ee114:1", "ee114:2",
     &    "ee122", "ee122:1", "ee122:2",
     &    "ee123", "ee123:1", "ee123:2",
     &    "ee124", "ee124:1", "ee124:2",
     &    "ee133", "ee133:1", "ee133:2",
     &    "ee134", "ee134:1", "ee134:2",
     &    "ee144", "ee144:1", "ee144:2",
     &    "ee222", "ee222:1", "ee222:2",
     &    "ee223", "ee223:1", "ee223:2",
     &    "ee224", "ee224:1", "ee224:2",
     &    "ee233", "ee233:1", "ee233:2",
     &    "ee234", "ee234:1", "ee234:2",
     &    "ee244", "ee244:1", "ee244:2",
     &    "ee333", "ee333:1", "ee333:2",
     &    "ee334", "ee334:1", "ee334:2",
     &    "ee344", "ee344:1", "ee344:2",
     &    "ee444", "ee444:1", "ee444:2",
     &    "ee0000", "ee0000:1", "ee0000:2",
     &    "ee0011", "ee0011:1", "ee0011:2",
     &    "ee0012", "ee0012:1", "ee0012:2",
     &    "ee0013", "ee0013:1", "ee0013:2",
     &    "ee0014", "ee0014:1", "ee0014:2",
     &    "ee0022", "ee0022:1", "ee0022:2",
     &    "ee0023", "ee0023:1", "ee0023:2",
     &    "ee0024", "ee0024:1", "ee0024:2",
     &    "ee0033", "ee0033:1", "ee0033:2",
     &    "ee0034", "ee0034:1", "ee0034:2",
     &    "ee0044", "ee0044:1", "ee0044:2",
     &    "ee1111", "ee1111:1", "ee1111:2",
     &    "ee1112", "ee1112:1", "ee1112:2",
     &    "ee1113", "ee1113:1", "ee1113:2",
     &    "ee1114", "ee1114:1", "ee1114:2",
     &    "ee1122", "ee1122:1", "ee1122:2",
     &    "ee1123", "ee1123:1", "ee1123:2",
     &    "ee1124", "ee1124:1", "ee1124:2",
     &    "ee1133", "ee1133:1", "ee1133:2",
     &    "ee1134", "ee1134:1", "ee1134:2",
     &    "ee1144", "ee1144:1", "ee1144:2",
     &    "ee1222", "ee1222:1", "ee1222:2",
     &    "ee1223", "ee1223:1", "ee1223:2",
     &    "ee1224", "ee1224:1", "ee1224:2",
     &    "ee1233", "ee1233:1", "ee1233:2",
     &    "ee1234", "ee1234:1", "ee1234:2",
     &    "ee1244", "ee1244:1", "ee1244:2",
     &    "ee1333", "ee1333:1", "ee1333:2",
     &    "ee1334", "ee1334:1", "ee1334:2",
     &    "ee1344", "ee1344:1", "ee1344:2",
     &    "ee1444", "ee1444:1", "ee1444:2",
     &    "ee2222", "ee2222:1", "ee2222:2",
     &    "ee2223", "ee2223:1", "ee2223:2",
     &    "ee2224", "ee2224:1", "ee2224:2",
     &    "ee2233", "ee2233:1", "ee2233:2",
     &    "ee2234", "ee2234:1", "ee2234:2",
     &    "ee2244", "ee2244:1", "ee2244:2",
     &    "ee2333", "ee2333:1", "ee2333:2",
     &    "ee2334", "ee2334:1", "ee2334:2", 
     &    "ee2344", "ee2344:1", "ee2344:2",
     &    "ee2444", "ee2444:1", "ee2444:2",
     &    "ee3333", "ee3333:1", "ee3333:2",
     &    "ee3334", "ee3334:1", "ee3334:2",
     &    "ee3344", "ee3344:1", "ee3344:2",
     &    "ee3444", "ee3444:1", "ee3444:2",
     &    "ee4444", "ee4444:1", "ee4444:2" /
	end

