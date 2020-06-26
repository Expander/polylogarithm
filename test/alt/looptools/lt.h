* lt.h
* internal common blocks for the LoopTools routines
* this file is part of LoopTools
* last modified 14 Apr 18 th


#include "ff.h"

* the cache-pointer structure is (see cache.c):
* 1. int valid
* 2. Node *last
* 3. Node *first
* 4. (not used)

	integer ncaches
	parameter (ncaches = 10)

	integer*8 cacheptr(4,0:QUAD,ncaches)
	integer*8 savedptr(2,ncaches)
	RealType maxdev
	integer epsi, warndigits, errdigits
	integer serial, versionkey
	integer debugkey, debugfrom, debugto

	common /ltvars/
     &    cacheptr, savedptr,
     &    maxdev,
     &    epsi, warndigits, errdigits,
     &    serial, versionkey,
     &    debugkey, debugfrom, debugto

	integer cmpbits

	common /ltcache/ cmpbits

	ComplexType cache(2,ncaches)
	equivalence (cacheptr, cache)

#ifndef DEBUGLEVEL
#define DEBUGLEVEL ibits(debugkey,8,2)
#endif

