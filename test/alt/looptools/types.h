* types.h
* real-based type declarations
* this file is part of LoopTools
* last modified 9 Apr 18 th


#ifndef TYPES_H
#define TYPES_H

#define RealType double precision
#define ComplexType double complex
#define Re DBLE
#define Im DIMAG
#define Conjugate DCONJG
#define ToComplex DCMPLX

#if QUADSIZE == 16
#define RealQuad real*16
#define ComplexQuad complex*32
#elif QUADSIZE == 10
#define RealQuad real*10
#define ComplexQuad complex*20
#else
#define RealQuad RealType
#define ComplexQuad ComplexType
#endif

#define Sq(c) Re((c)*Conjugate(c))
#define Sqrtc(c) sqrt(ToComplex(c))

#endif

