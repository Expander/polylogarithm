* Li2.F
* the dilogarithm function
* this file is part of LoopTools
* last modified 7 Apr 18 th

#include "externals.h"
#include "types.h"
#include "defs.h"


	ComplexType function XLi2(x)
	implicit none
	ArgType x

	RealType pi12
	parameter (pi12 = .822467033424113218236207583323D0)

	ComplexType res, dummy
	integer ier, ipi12

	ier = 0
#ifdef COMPLEXPARA
	call ffzzdl(res, ipi12, dummy, x, ier)
#else
	call ffzxdl(res, ipi12, dummy, x, -1, ier)
#endif
	XLi2 = res + ipi12*pi12
	end

************************************************************************
* adapter code for C++

	subroutine XLi2sub(res, x)
	implicit none
	ComplexType res
	ArgType x

	ComplexType XLi2
	external XLi2

	res = XLi2(x)
        end

