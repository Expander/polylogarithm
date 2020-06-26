	character*80 e1
	parameter (e1="ffca0: minimum value complex logarit"//
     +   "hm gives problem, change mu.")
	character*80 e2
	parameter (e2="ffxa0: minimum value real logarithm "//
     +   "gives problem, change mu.")
	character*80 e3
	parameter (e3="ffcb0: minimum value complex logarit"//
     +   "hm gives problem, change mu.")
	character*80 e4
	parameter (e4="ffxb0: minimum value real logarithm "//
     +   "gives problem, change mu.")
	character*80 e5
	parameter (e5="ffcb0p: cannot handle complex k^2 yet")
	character*80 e6
	parameter (e6="ffcb0p: minimum value complex log giv"//
     +   "es problem in unequal masses.")
	character*80 e7
	parameter (e7="ffxb0p: divergence for k->0, m1=m2=0.")
	character*80 e8
	parameter (e8="ffxb0p: minimum value real log gives "//
     +   "problem in equal masses.")
	character*80 e9
	parameter (e9="ffxb0p: minimum value real log gives "//
     +   "problem in unequal masses.")
	character*80 e10
	parameter (e10="ffcc0p: cannot handle two spacelike m"//
     +   "omenta and one zero.")
	character*80 e11
	parameter (e11="ffxc0p: cannot handle two spacelike m"//
     +   "omenta and one zero.")
	character*80 e12
	parameter (e12="ffcs3: illegal code for isoort(1) (s"//
     +   "hould not occur)")
	character*80 e13
	parameter (e13="ffcs3: illegal code for isoort(2) (s"//
     +   "hould not occur)")
	character*80 e14
	parameter (e14="ffcs3: imaginary part wrong, will be"//
     +   " improved later")
	character*80 e15
	parameter (e15="ffcs3: isoort = -1,0 not yet ready")
	character*80 e16
	parameter (e16="ffcs3: illegal combination in isoort"//
     +   " (should not occur)")
	character*80 e17
	parameter (e17="ffcxs3: illegal code for isoort(1) (s"//
     +   "hould not occur)")
	character*80 e18
	parameter (e18="ffcxs3: illegal code for isoort(2) (s"//
     +   "hould not occur)")
	character*80 e19
	parameter (e19="ffcs4: imaginary part is wrong (shou"//
     +   "ld be updated)")
	character*80 e20
	parameter (e20="ffdcrr: Taylor expansion in 1/x not y"//
     +   "et ready")
	character*80 e21
	parameter (e21="ffdcxr: imaginary part is wrong")
	character*80 e22
	parameter (e22="ffdcxr: Taylor expansion in 1/x not y"//
     +   "et ready")
	character*80 e23
	parameter (e23="ffcrr: minimum value complex log cau"//
     +   "ses correction term to be wrong.")
	character*80 e24
	parameter (e24="ffcxr: minimum value real log causes"//
     +   " correction term to be wrong.")
	character*80 e25
	parameter (e25="ffcrr: illegal code for iclas1 (shou"//
     +   "ld not occur)")
	character*80 e26
	parameter (e26="ffcxr: illegal code for iclas1 (shou"//
     +   "ld not occur)")
	character*80 e27
	parameter (e27="ffcrr: illegal code for iclas2 (shou"//
     +   "ld not occur)")
	character*80 e28
	parameter (e28="ffcxr: illegal code for iclas2 (shou"//
     +   "ld not occur)")
	character*80 e29
	parameter (e29="ffxli2: argument too large (should no"//
     +   "t occur)")
	character*80 e30
	parameter (e30="ffzli2: argument too large (should no"//
     +   "t occur)")
	character*80 e31
	parameter (e31="ffzzdl: imaginary part dilog is undef"//
     +   "ined for real x > 1.")
	character*80 e32
	parameter (e32="nffeta: eta is not defined for real n"//
     +   "egative numbers a,b, ab.")
	character*80 e33
	parameter (e33="nffet1: eta is not defined for real n"//
     +   "egative numbers a,b, ab.")
	character*80 e34
	parameter (e34="ffcota: illegal flag (should not occu"//
     +   "r)")
	character*80 e35
	parameter (e35="ffrota: illegal flag (should not occu"//
     +   "r)")
	character*80 e36
	parameter (e36="ffccyz: I took the wrong value for ca"//
     +   "lpha... (should not occur)")
	character*80 e37
	parameter (e37="ffxxyz: I took the wrong value for al"//
     +   "pha... (should not occur)")
	character*80 e38
	parameter (e38="ffcoot: a=0, trying to find two roots"//
     +   " of a linear equation ...")
	character*80 e39
	parameter (e39="ffroot: a=0, trying to find two roots"//
     +   " of a linear equation ...")
	character*80 e40
	parameter (e40="ffrot3: all three external masses zer"//
     +   "o !")
	character*80 e41
	parameter (e41="ffxc0: lambda(p1,p2,p3) < 0, unphysi"//
     +   "cal configuration")
	character*80 e42
	parameter (e42="ffxc0: cannot handle this case (p1,p"//
     +   "2,p3 dependent, on threshold)")
	character*80 e43
	parameter (e43="ffcxs3: illegal code for isoort(1) (s"//
     +   "hould not occur)")
	character*80 e44
	parameter (e44="ffxd0: lambda(p1,p2,p3,p4) < 0, unph"//
     +   "ysical configuration")
	character*80 e45
	parameter (e45="ffxd0: cannot handle this case (p1,p"//
     +   "2,p3 dependent, on threshold)")
	character*80 e46
	parameter (e46="ffxd0p: correction terms for Ai <0 in"//
     +   "finite (mass zero?)")
	character*80 e47
	parameter (e47="ffcxyz: p_i^2 = 0 (should not occur)")
	character*80 e48
	parameter (e48="ffeta: answer not consistent with no"//
     +   "rmal result (old)")
	character*80 e49
	parameter (e49="ffcc0: cannot handle complex externa"//
     +   "l momenta or im > 0")
	character*80 e50
	parameter (e50="ffcd0: cannot handle complex externa"//
     +   "l momenta.")
	character*80 e51
	parameter (e51="zfflog: imaginary part undefined for "//
     +   "real z < 0.")
	character*80 e52
	parameter (e52="zxfflg: imaginary part undefined for "//
     +   "x < 0.")
	character*80 e53
	parameter (e53="ffcs3: eta changes within (0,1), add"//
     +   " sophisticated terms...")
	character*80 e54
	parameter (e54="ffrot4: cannot find any physical vert"//
     +   "ex to apply transformation.")
	character*80 e55
	parameter (e55="fftra0: too many vectors parallel, p_"//
     +   "1.p_7 or p_2.p_7 is zero.")
	character*80 e56
	parameter (e56="zfflog: tiny imaginary part in confli"//
     +   "ct with ieps prescription.")
	character*80 e57
	parameter (e57="ffxe0: lambda(p1,p2,p3,p4,p5) < 0, u"//
     +   "nphysical")
	character*80 e58
	parameter (e58="ffxc0j: IR divergent C0 with lambda(p"//
     +   "1,p2,p3)=0.")
	character*80 e59
	parameter (e59="ffxc0i: IR divergent C0 with lambda2=0.")
	character*80 e60
	parameter (e60="ffxc0j: IR divergent C0 obtained from"//
     +   " D0 is singular. Contact author.")
	character*80 e61
	parameter (e61="ffxd0p: IR divergent D0 with lambda2=0.")
	character*80 e62
	parameter (e62="ffxc0p: I never expected complex root"//
     +   "s in an IR divergent diagram.")
	character*80 e63
	parameter (e63="ffxd0p: can only handle one IR diverg"//
     +   "ence per 3point function")
	character*80 e64
	parameter (e64="ffxd0p: cannot handle a threshold in"//
     +   " (3,4), rotated wrongly.")
	character*80 e65
	parameter (e65="ffcxr: IR divergence but iclass!=3. "//
     +   " should not occur.")
	character*80 e66
	parameter (e66="ffcxs3: different imaginary signs sho"//
     +   "uld not occur for ipole=3.")
	character*80 e67
	parameter (e67="ffxdbd: I cannot use this algorithm f"//
     +   "or a linear IR divergence")
	character*80 e68
	parameter (e68="ffxd0: cannot find a proj. transform"//
     +   "ation; try another permutation.")
	character*80 e69
	parameter (e69="ff5ind: could not find independent mo"//
     +   "menta (should not occur).")
	character*80 e70
	parameter (e70="ffxdna: lambda(pi,pj,pk) < 0, unphysi"//
     +   "cal configuration")
	character*80 e71
	parameter (e71="ffxdna: cannot handle lambda(pi,pj,pk"//
     +   ") = 0, dependent momenta.")
	character*80 e72
	parameter (e72="ffxd0e: could not find a stable root;"//
     +   " please try another permutation")
	character*80 e73
	parameter (e73="ffxdir: cannot handle a linearly dive"//
     +   "rgent four point function (yet)")
	character*80 e74
	parameter (e74="ffxdbd: IR divergent B0' without cuto"//
     +   "ff in /ffregul/")
	character*80 e75
	parameter (e75="ffdcxr: dyz=0, should not occur")
	character*80 e76
	parameter (e76="ffdcrr: cdwz=0, but iepsz!=iepsz and "//
     +   "significant")
	character*80 e77
	parameter (e77="ffdcrr: cdyz=0, should not occur")
	character*80 e78
	parameter (e78="ffdcc0: imaginary part wrong")
	character*80 e79
	parameter (e79="ffdcs: cannot handle isoort=0")
	character*80 e80
	parameter (e80="ffdcs: mixed up iep's, 2*pi^2 wrong "//
     +   "somewhere")
	character*80 e81
	parameter (e81="ffdcs: wrong value for isoort")
	character*80 e82
	parameter (e82="ffdxc0: imaginary part Ai < 0 terms unc"//
     +   "ertain")
	character*80 e83
	parameter (e83="ffxc0j: sorry, complex roots not yet "//
     +   "supported here")
	character*80 e84
	parameter (e84="ffxc0p: imaginary part Ai < 0 terms unc"//
     +   "ertain")
	character*80 e85
	parameter (e85="ffxd0a: t3 = t4, don''t know what to do")
	character*80 e86
	parameter (e86="ffxdbp: cannot compute derivative, la"//
     +   "m=0")
	character*80 e87
	parameter (e87="ffxdi: dependent momenta not yet sup"//
     +   "ported (boundary of phase space)")
	character*80 e88
	parameter (e88="ffxxyz: xk = 0 not yet implemented")
	character*80 e92
	parameter (e92="ffxc1: cannot invert matrix with zer"//
     +   "o determinant.")
	character*80 e93
	parameter (e93="ffze0: Im(m^2) > 0")
	character*80 e94
	parameter (e94="ffze0: Im(p^2) != 0")
	character*80 e95
	parameter (e95="ffzf0: Im(m^2) > 0")
	character*80 e96
	parameter (e96="ffzf0: Im(p^2) != 0")
	character*80 e97
	parameter (e97="ffxc0j: ill-defined IR-divergent C0 "//
     +   "for massless charged particles.")
	character*80 e98
	parameter (e98="ffxdbd: ill-defined IR-divergent D0 "//
     +   "for massless charged particles.")
	character*80 e100
	parameter (e100="ffrcvr: probably underflow, I do"//
     +   " not know where or how severe.")
	character*80 e101
	parameter (e101="ffxdb1: case not defined")
	character*80 e102
	parameter (e102="ffxdb11: case not defined")
	character*80 e103
	parameter (e103="ffd0c: cannot handle this case")
	character*80 e104
	parameter (e104="ffwbeta: prefactor 1/(SV-TU) = 1/0 "//
     +    "for all y")
	character*80 e105
	parameter (e105="ffT_lin: prefactor 1/(SV-TU) = 1/0 "//
     +    "for all y")
	character*80 e99
	parameter (e99="ffT13: prefactor 1/(SV-TU) = 1/0 "//
     +    "for all y")
	character*80 e89
	parameter (e89="ffS2: log(0) singularity")
	character*80 e90
	parameter (e90="ffS3n: end-point singularity")
	character*80 e91
	parameter (e91="ffS3n: log(0) singularity")
	character*80 error(105)
	data error / e1,e2,e3,e4,e5,e6,e7,e8,e9,
     +   e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,
     +   e20,e21,e22,e23,e24,e25,e26,e27,e28,e29,
     +   e30,e31,e32,e33,e34,e35,e36,e37,e38,e39,
     +   e40,e41,e42,e43,e44,e45,e46,e47,e48,e49,
     +   e50,e51,e52,e53,e54,e55,e56,e57,e58,e59,
     +   e60,e61,e62,e63,e64,e65,e66,e67,e68,e69,
     +   e70,e71,e72,e73,e74,e75,e76,e77,e78,e79,
     +   e80,e81,e82,e83,e84,e85,e86,e87,e88,e89,
     +   e90,e91,e92,e93,e94,e95,e96,e97,e98,e99,
     +   e100,e101,e102,e103,e104,e105 /
