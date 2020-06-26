	character*80 w1
	parameter (w1="ffcb0p: warning: instability in case one mas"//
     +   "s zero, may be solved later.")
	character*80 w2
	parameter (w2="ffcb0p: warning: not enough terms in Taylor "//
     +   "expansion ma=mb. May be serious!")
	character*80 w3
	parameter (w3="ffcb0p: warning: minimum value complex logar"//
     +   "ithm gives problem in equal masses.")
	character*80 w4
	parameter (w4="ffcb0p: warning: cancellations in equal mass"//
     +   "es (should not occur).")
	character*80 w5
	parameter (w5="ffcb0p: warning: not enough terms in expansi"//
     +   "on1 k2 zero. May be serious!")
	character*80 w6
	parameter (w6="ffcb0p: warning: not enough terms in expansi"//
     +   "on2 k2 zero, May be serious!")
	character*80 w7
	parameter (w7="ffcb0p: warning: cancellations in final addi"//
     +   "ng up, contact author if serious.")
	character*80 w8
	parameter (w8="ffc1lg: warning: the combination 1-z*log(1-1"//
     +   "/z) id unstable.")
	character*80 w9
	parameter (w9="ffcayl: warning: not enough terms in Taylor "//
     +   "expansion, may be serious.")
	character*80 w10
	parameter (w10="ffcb0p: warning: cancellation in dotproduct "//
     +   "s1.s2")
	character*80 w11
	parameter (w11="ffcot2: warning: cancellation in dotproduct "//
     +   "p.si ")
	character*80 w12
	parameter (w12="ffcdbp: warning: not enough terms in Taylor "//
     +   "expansion, may be serious")
	character*80 w13
	parameter (w13="ffcdbp: warning: cancellations in case one m"//
     +   "ass equal to zero")
	character*80 w14
	parameter (w14="ffxb0p: warning: instability in case one mas"//
     +   "s zero, may be solved later.")
	character*80 w15
	parameter (w15="ffxb0p: warning: not enough terms in Taylor "//
     +   "expansion ma=mb. May be serious!")
	character*80 w16
	parameter (w16="ffxb0p: warning: minimum value real logarith"//
     +   "m gives problem in equal masses.")
	character*80 w17
	parameter (w17="ffxb0p: warning: cancellations in equal mass"//
     +   "es (should not occur).")
	character*80 w18
	parameter (w18="ffxb0p: warning: cancellations in equal mass"//
     +   "es, complex roots, can be avoided.")
	character*80 w19
	parameter (w19="ffxb0p: warning: not enough terms in expansi"//
     +   "on1 k2 zero, may be serious!")
	character*80 w20
	parameter (w20="ffxb0p: warning: not enough terms in expansi"//
     +   "on2 k2 zero, may be serious!")
	character*80 w21
	parameter (w21="ffxb0p: warning: cancellations between s2 an"//
     +   "d alpha, should not be serious")
	character*80 w22
	parameter (w22="ffd1lg: warning: the combination 1-z*log(1-1"//
     +   "/z) id unstable.")
	character*80 w23
	parameter (w23="ffxb0p: warning: cancellations in lambda equ"//
     +   "al masses.")
	character*80 w24
	parameter (w24="ffxb0p: warning: cancellation in dotproduct "//
     +   "s1.s2")
	character*80 w25
	parameter (w25="ffdot2: warning: cancellation in dotproduct "//
     +   "p.si")
	character*80 w26
	parameter (w26="ffcc0:  warning: cancellation between the tw"//
     +   "o twopoint functions.")
	character*80 w27
	parameter (w27="ffcc0:  warning: cancellation in final summi"//
     +   "ng up.")
	character*80 w28
	parameter (w28="ffxc0:  warning: cancellation between the tw"//
     +   "o twopoint functions.")
	character*80 w29
	parameter (w29="ffxc0:  warning: cancellation in final summi"//
     +   "ng up.")
	character*80 w30
	parameter (w30="ffcc0p: warning: numerical problems in cw(j+"//
     +   "2,1), not used")
	character*80 w31
	parameter (w31="ffcc0p: warning: cancellations in cdwz(j,i,1"//
     +   "), not used")
	character*80 w32
	parameter (w32="ffcc0p: warning: numerical problems in cw(j+"//
     +   "2,3), not used")
	character*80 w33
	parameter (w33="ffcc0p: warning: cancellations in cdwz(j,i,3"//
     +   "), not used")
	character*80 w34
	parameter (w34="ffxc0p: warning: numerical problems in w(j+2"//
     +   ",1), not used")
	character*80 w35
	parameter (w35="ffxc0p: warning: cancellations in dwz(j,i,1)"//
     +   ", not used")
	character*80 w36
	parameter (w36="ffxc0p: warning: numerical problems in cw(j+"//
     +   "2,1), not used")
	character*80 w37
	parameter (w37="ffxc0p: warning: cancellations in cdwz(j,i,1"//
     +   "), not used")
	character*80 w38
	parameter (w38="ffxc0p: warning: numerical problems in w(j+2"//
     +   ",3), not used")
	character*80 w39
	parameter (w39="ffxc0p: warning: cancellations in dwz(j,i,3)"//
     +   ", not used")
	character*80 w40
	parameter (w40="ffxc0p: warning: numerical problems in cw(j+"//
     +   "2,3), not used")
	character*80 w41
	parameter (w41="ffxc0p: warning: cancellations in cdwz(j,i,3"//
     +   "), not used")
	character*80 w42
	parameter (w42="ffcs3:  warning: problems with range complex"//
     +   " numbers")
	character*80 w43
	parameter (w43="ffcs3:  warning: cancellations in czz1 in sp"//
     +   "ecial case")
	character*80 w44
	parameter (w44="ffcxs3: warning: cancellations in zz1 in spe"//
     +   "cial case")
	character*80 w45
	parameter (w45="ffdcrr: warning: not enough terms in Taylor "//
     +   "series (may be serious)")
	character*80 w46
	parameter (w46="ffdcxr: warning: not enough terms in Taylor "//
     +   "series (may be serious)")
	character*80 w47
	parameter (w47="ffcrr:  warning: problems with dynamical ran"//
     +   "ge complex numbers")
	character*80 w48
	parameter (w48="ffcrr:  warning: y0 = y1, so R has been take"//
     +   "n zero")
	character*80 w49
	parameter (w49="ffcrr:  warning: very large correction terms.")
	character*80 w50
	parameter (w50="ffcrr:  warning: minimum value complex log c"//
     +   "auses loss of precision.")
	character*80 w51
	parameter (w51="ffcxr:  warning: y0 = y1, so R has been take"//
     +   "n zero")
	character*80 w52
	parameter (w52="ffcxr:  warning: very large correction terms.")
	character*80 w53
	parameter (w53="ffcxr:  warning: minimum value real log caus"//
     +   "es loss of precision.")
	character*80 w54
	parameter (w54="ffcrr:  warning: not enough terms in Taylor "//
     +   "series (may be serious)")
	character*80 w55
	parameter (w55="ffcxr:  warning: not enough terms in Taylor "//
     +   "series (may be serious)")
	character*80 w56
	parameter (w56="ffcrr:  warning: cancellations in cd2yzz + c"//
     +   "zz")
	character*80 w57
	parameter (w57="ffcrr:  warning: cancellations in cd2yzz - c"//
     +   "zz1")
	character*80 w58
	parameter (w58="ffcxr:  warning: cancellations in d2yzz + zz")
	character*80 w59
	parameter (w59="ffcxr:  warning: cancellations in d2yzz - zz1")
	character*80 w60
	parameter (w60="ffxli2: warning: not enough terms in expansi"//
     +   "on (may be serious)")
	character*80 w61
	parameter (w61="ffzli2: warning: not enough terms in expansi"//
     +   "on (may be serious)")
	character*80 w62
	parameter (w62="dfflo1: warning: not enough terms in expansi"//
     +   "on. calling log.")
	character*80 w63
	parameter (w63="zfflo1: warning: not enough terms in expansi"//
     +   "on. calling log.")
	character*80 w64
	parameter (w64="ffzxdl: warning: minimum value real log give"//
     +   "s problems.")
	character*80 w65
	parameter (w65="ffzzdl: warning: minimum value complex log g"//
     +   "ives problems.")
	character*80 w66
	parameter (w66="ffzxdl: warning: not enough terms in expansi"//
     +   "on (may be serious)")
	character*80 w67
	parameter (w67="ffzzdl: warning: not enough terms in expansi"//
     +   "on (may be serious)")
	character*80 w68
	parameter (w68="ffclmb: warning: cancellation in calculation"//
     +   " lambda.")
	character*80 w69
	parameter (w69="ffxlmb: warning: cancellation in calculation"//
     +   " lambda.")
	character*80 w70
	parameter (w70="ffcel2: warning: cancellation in calculation"//
     +   " delta_{pi pj}^{pi pj}")
	character*80 w71
	parameter (w71="ffdel2: warning: cancellation in calculation"//
     +   " delta_{pi pj}^{pi pj}")
	character*80 w72
	parameter (w72="ffcel3: warning: cancellation in calculation"//
     +   " delta_{s1 s2 s3}^{s1 s2 s3}")
	character*80 w73
	parameter (w73="ffdel3: warning: cancellation in calculation"//
     +   " delta_{s1 s2 s3}^{s1 s2 s3}")
	character*80 w74
	parameter (w74="ffcl3m: warning: cancellation in (delta_{sj"//
     +   " sk}^{si mu})^2")
	character*80 w75
	parameter (w75="ffdl3m: warning: cancellation in (delta_{sj"//
     +   " sk}^{si mu})^2")
	character*80 w76
	parameter (w76="ffeta:  warning: still cancellations. (not u"//
     +   "sed)")
	character*80 w77
	parameter (w77="ffceta: warning: still cancellations. (not u"//
     +   "sed)")
	character*80 w78
	parameter (w78="ffcdwz: warning: still cancelations in cw3pm"//
     +   " - cz3mp (not used)")
	character*80 w79
	parameter (w79="ffdwz:  warning: still cancelations in w3pm "//
     +   "- z3mp (not used)")
	character*80 w80
	parameter (w80="ffdcxr: warning: minimum value real log caus"//
     +   "es problems.")
	character*80 w81
	parameter (w81="ffdcxr: warning: ieps <> iepsz, imaginary pa"//
     +   "rt will be wrong")
	character*80 w82
	parameter (w82="ffdcrr: warning: minimum value complex log c"//
     +   "auses problems.")
	character*80 w83
	parameter (w83="ffdl2s: warning: cancellations in delta_{s1'"//
     +   "s2'}^{s1 s2}")
	character*80 w84
	parameter (w84="ffxd0:  warning: cancellation in final summi"//
     +   "ng up.")
	character*80 w85
	parameter (w85="ffdl3s: warning: cancellation in calculation"//
     +   " delta^(si sj sk)_(sl sm sn)")
	character*80 w86
	parameter (w86="ffcc0:  warning: cancellations among input p"//
     +   "arameters")
	character*80 w87
	parameter (w87="ffxc0:  warning: cancellations among input p"//
     +   "arameters (import difference)")
	character*80 w88
	parameter (w88="ffabcd: warning: cancellations in (2*s3.s4^2"//
     +   " - s3^2*s4^2), try with del2")
	character*80 w89
	parameter (w89="ffabcd: warning: cancellations in somb")
	character*80 w90
	parameter (w90="ffabcd: warning: cancellations in d")
	character*80 w91
	parameter (w91="ffabcd: warning: xc not yet accurate (can be"//
     +   " improved)")
	character*80 w92
	parameter (w92="ffdl2p: warning: cancellations in delta_{p1"//
     +   " s2}^{p1 p2}")
	character*80 w93
	parameter (w93="ffdl2t: warning: cancellations in delta_{p1"//
     +   " s4}^{s3 s4}")
	character*80 w94
	parameter (w94="ffcb0:  warning: cancellations between cma a"//
     +   "nd cmb (add input parameters)")
	character*80 w95
	parameter (w95="ffcb0:  warning: cancellations between ck an"//
     +   "d cma (add input parameters)")
	character*80 w96
	parameter (w96="ffcb0:  warning: cancellations between ck an"//
     +   "d cmb (add input parameters)")
	character*80 w97
	parameter (w97="ffxb0:  warning: cancellations between xma a"//
     +   "nd xmb (add input parameters)")
	character*80 w98
	parameter (w98="ffxb0:  warning: cancellations between xk an"//
     +   "d xma (add input parameters)")
	character*80 w99
	parameter (w99="ffxb0:  warning: cancellations between xk an"//
     +   "d xmb (add input parameters)")
	character*80 w100
	parameter (w100="ffdot3: warning: cancellations in dotproduct"//
     +   " s_i.s_{i+1}")
	character*80 w101
	parameter (w101="ffdot3: warning: cancellations in dotproduct"//
     +   " p_i.s_i")
	character*80 w102
	parameter (w102="ffdot3: warning: cancellations in dotproduct"//
     +   " p_i.s_{i+1}")
	character*80 w103
	parameter (w103="ffdot3: warning: cancellations in dotproduct"//
     +   " p_i.s_{i+2}")
	character*80 w104
	parameter (w104="ffdot3: warning: cancellations in dotproduct"//
     +   " p_i.p_{i+1}")
	character*80 w105
	parameter (w105="ffdot4: warning: cancellations in dotproduct"//
     +   " s_i.s_{i+1}")
	character*80 w106
	parameter (w106="ffdot4: warning: cancellations in dotproduct"//
     +   " s_i.s_{i-1}")
	character*80 w107
	parameter (w107="ffdot4: warning: cancellations in dotproduct"//
     +   " p_i.s_i")
	character*80 w108
	parameter (w108="ffdot4: warning: cancellations in dotproduct"//
     +   " p_i.s_{i+1}")
	character*80 w109
	parameter (w109="ffdot4: warning: cancellations in dotproduct"//
     +   " p_{i-1}.s_i")
	character*80 w110
	parameter (w110="ffdot4: warning: cancellations in dotproduct"//
     +   " p_i.s_{i+2}")
	character*80 w111
	parameter (w111="ffdot4: warning: cancellations in dotproduct"//
     +   " p_{i+1}.s_i")
	character*80 w112
	parameter (w112="ffdot4: warning: cancellations in dotproduct"//
     +   " p_{i+2}.s_{i+1}")
	character*80 w113
	parameter (w113="ffdot4: warning: cancellations in dotproduct"//
     +   " p_i.p_{i+1}")
	character*80 w114
	parameter (w114="ffdot4: warning: cancellations in dotproduct"//
     +   " p_{i+1}.p_{i+2}")
	character*80 w115
	parameter (w115="ffdot4: warning: cancellations in dotproduct"//
     +   " p_{i+2}.p_i")
	character*80 w116
	parameter (w116="ffdot4: warning: cancellations in dotproduct"//
     +   " p_5.p_7")
	character*80 w117
	parameter (w117="ffdot4: warning: cancellations in dotproduct"//
     +   " p_6.p_8")
	character*80 w118
	parameter (w118="ffdot4: warning: cancellations in dotproduct"//
     +   " p_9.p_10")
	character*80 w119
	parameter (w119="ffxd0:  warning: sum is close to the minimum"//
     +   " of the range.")
	character*80 w120
	parameter (w120="ffxc0:  warning: sum is close to the minimum"//
     +   " of the range.")
	character*80 w121
	parameter (w121="ffxd0:  warning: cancellations among input p"//
     +   "arameters (import difference)")
	character*80 w122
	parameter (w122="ff2d22: warning: cancellations (delta_{sjsk"//
     +   "}_{si mu} delta_{smsn}^{mu nu})^2")
	character*80 w123
	parameter (w123="ff2dl2: warning: cancellations delta^{si mu"//
     +   "}_{sj sk} delta^{mu sl}_{sm sn}")
	character*80 w124
	parameter (w124="ff3dl2: warning: cancellations d^{i mu}_{jl"//
     +   "} d^{mu nu}_{lm} d^{nu n}_{op}")
	character*80 w125
	parameter (w125="fftran: warning: cancellations in s'_i^2 - s"//
     +   "'_j^2")
	character*80 w126
	parameter (w126="fftran: warning: cancellations in p'_i^2 - s"//
     +   "'_j^2")
	character*80 w127
	parameter (w127="fftran: warning: cancellations in p'_i^2 - p"//
     +   "'_j^2")
	character*80 w128
	parameter (w128="zfflog: warning: taking log of number close "//
     +   "to 1, must be cured.")
	character*80 w129
	parameter (w129="zxfflg: warning: taking log of number close "//
     +   "to 1, must be cured.")
	character*80 w130
	parameter (w130="ffcrr:  warning: cancellations in calculatin"//
     +   "g 2y-1-z...")
	character*80 w131
	parameter (w131="ffxtra: warning: cancellations in extra term"//
     +   "s, working on it")
	character*80 w132
	parameter (w132="dfflo1: warning: cancellations because of wr"//
     +   "ong call, should not occur")
	character*80 w133
	parameter (w133="zfflo1: warning: cancellations because of wr"//
     +   "ong call, should not occur")
	character*80 w134
	parameter (w134="ffcs4:  warning: cancellations in cd2yzz + c"//
     +   "zz")
	character*80 w135
	parameter (w135="ffcd0:  warning: cancellations among input p"//
     +   "arameters (import difference)")
	character*80 w136
	parameter (w136="ffcd0:  warning: cancellation in final summi"//
     +   "ng up.")
	character*80 w137
	parameter (w137="ffcd0:  warning: sum is close to the minimum"//
     +   " of the range.")
	character*80 w138
	parameter (w138="ffdl3p: warning: cancellations in delta_{p1"//
     +   " p2 p3}^{p1 p2 p3}")
	character*80 w139
	parameter (w139="ffxd0p: warning: problems calculating sqrt(d"//
     +   "elta(si,s3)) - sqrt(delta(si,s4))")
	character*80 w140
	parameter (w140="ffdxc0: warning: problems calculating yzzy ="//
     +   " y(4)z(3) - y(3)z(4)")
	character*80 w141
	parameter (w141="ffcd0p: warning: problems calculating sqrt(d"//
     +   "elta(si,s3)) - sqrt(delta(si,s4))")
	character*80 w142
	parameter (w142="ffdcc0: warning: problems calculating yzzy ="//
     +   " y(4)z(3) - y(3)z(4)")
	character*80 w143
	parameter (w143="ffdel4: warning: cancellation in calculation"//
     +   " delta_{s1 s2 s3 s4}^{s1 s2 s3 s4}")
	character*80 w144
	parameter (w144="fftran: warning: cancellation in calculation"//
     +   " s_i'.p_{jk}'")
	character*80 w145
	parameter (w145="fftran: warning: cancellation in calculation"//
     +   " p_{ji}'.p_{lk}'")
	character*80 w146
	parameter (w146="fftran: warning: cancellation in calculation"//
     +   " Ai - Aj")
	character*80 w147
	parameter (w147="ffdxc0: warning: problems calculating yyzz ="//
     +   " y(4) - y(3) - z(3) + z(4)")
	character*80 w148
	parameter (w148="ffdxc0: warning: problems calculating cancel"//
     +   "lations extra terms")
	character*80 w149
	parameter (w149="ffcb0:  warning: cancellations between Delta"//
     +   ", B0' and log(m1*m2/mu^2)/2")
	character*80 w150
	parameter (w150="ffxb0:  warning: cancellations between Delta"//
     +   ", B0' and log(m1*m2/mu^2)/2")
	character*80 w151
	parameter (w151="ffzli2: warning: real part complex dilog ver"//
     +   "y small and not stable")
	character*80 w152
	parameter (w152="ffxxyz: warning: cancellations in y - 2*z (w"//
     +   "ill be solved)")
	character*80 w153
	parameter (w153="ffxd0:  warning: cancellation in u=+p5^2+p6^"//
     +   "2+p7^2+p8^2-p9^2-p10^2, import it!")
	character*80 w154
	parameter (w154="ffxd0:  warning: cancellation in v=-p5^2+p6^"//
     +   "2-p7^2+p8^2+p9^2+p10^2, import it!")
	character*80 w155
	parameter (w155="ffxd0:  warning: cancellation in w=+p5^2-p6^"//
     +   "2+p7^2-p8^2+p9^2+p10^2, import it!")
	character*80 w156
	parameter (w156="ffxc0i: warning: cancellations in dotproduct"//
     +   " p_i.s_j")
	character*80 w157
	parameter (w157="ffxc0i: warning: cancellations in final summ"//
     +   "ing up")
	character*80 w158
	parameter (w158="ffxe0:  warning: cancellations among input p"//
     +   "arameters (import difference)")
	character*80 w159
	parameter (w159="ffdl4p: warning: cancellations in delta_{p1"//
     +   " p2 p3 p4}^{p1 p2 p3 p4}")
	character*80 w160
	parameter (w160="ffdel5: warning: cancellation in calculation"//
     +   " delta_{s1s2s3s4s5}^{s1s2s3s4s5}")
	character*80 w161
	parameter (w161="ffxe0a: warning: cancellation in final summi"//
     +   "ng up.")
	character*80 w162
	parameter (w162="ffxe0a: warning: sum is close to the minimum"//
     +   " of the range.")
	character*80 w163
	parameter (w163="ffxc1:  warning: cancellations in cc1.")
	character*80 w164
	parameter (w164="ffxd1:  warning: cancellations in cd1.")
	character*80 w165
	parameter (w165="ffdl2i: warning: cancellations in delta_{p1"//
     +   " p2}^{p3 p4}")
	character*80 w166
	parameter (w166="ffdl3q: warning: cancellations in delta_{p5"//
     +   " p6 p7}^{p(i1) p(i2) p(i3)}")
	character*80 w167
	parameter (w167="ffxb1:  warning: cancellations in cb1.")
	character*80 w168
	parameter (w168="ffxe0:  warning: cancellations in (p_i+p_{i+"//
     +   "2})^2 (may not be serious)")
	character*80 w169
	parameter (w169="ffdl4r: warning: cancellations in delta_{p1"//
     +   " p2 p3 p4}^{s1 s2 s3 s4}")
	character*80 w170
	parameter (w170="ffdl4s: warning: cancellations in delta_{p1"//
     +   "p2p3p4}^{si pj pk pl}, to be improved")
	character*80 w171
	parameter (w171="ffxe1:  warning: cancellations in ce1")
	character*80 w172
	parameter (w172="ffceta: warning: cancellations in extra term"//
     +   "s for 4point function")
	character*80 w173
	parameter (w173="ffceta: warning: cancellations between alpha"//
     +   " and w-")
	character*80 w174
	parameter (w174="ffceta: warning: cancellations between alpha"//
     +   " and w+")
	character*80 w175
	parameter (w175="ffceta: warning: cancellations between a and"//
     +   " z")
	character*80 w176
	parameter (w176="ffceta: warning: cancellations between a and"//
     +   " y")
	character*80 w177
	parameter (w177="ffcdbd: warning: cancellations in summing up")
	character*80 w178
	parameter (w178="ffkfun: warning: cancellations between z and"//
     +   " (m-mp)^2")
	character*80 w179
	parameter (w179="ffkfun: warning: 4*m*mp/(z-(m-mp)^2) ~ 1, ca"//
     +   "n be solved")
	character*80 w180
	parameter (w180="ffxc0p: warning: delta^{s1,s2,s3}_{s1,s2,s3"//
     +   "} not stable, can be solved.")
	character*80 w181
	parameter (w181="ffxc0p: warning: cancellations in complex di"//
     +   "scriminant, can be solved")
	character*80 w182
	parameter (w182="ffcd0e: warning: still cancellations in del4"//
     +   " with only complex in poles")
	character*80 w183
	parameter (w183="ffcc0a: warning: cannot deal properly with t"//
     +   "hreshold of this type")
	character*80 w184
	parameter (w184="ffcran: warning: cancellations in s'(i).p'(k"//
     +   "j)")
	character*80 w185
	parameter (w185="ffcran: warning: cancellations in p'(ji).p'("//
     +   "lk)")
	character*80 w186
	parameter (w186="ffcd0p: warning: cancellations in cel2")
	character*80 w187
	parameter (w187="ffdel6: warning: cancellations in coefficien"//
     +   "t F0, can be improved")
	character*80 w188
	parameter (w188="ffdl5r: warning: cancellations in coefficien"//
     +   "t E0, can be improved")
	character*80 w189
	parameter (w189="ffxdi:  warning: cancellations in cd2del")
	character*80 w190
	parameter (w190="ffxdi:  warning: cancellations in cd2pp")
	character*80 w191
	parameter (w191="ffxf0a: warning: cancellations in F0 as sum "//
     +   "of 6 E0's - near threshold?")
	character*80 w192
	parameter (w192="ffxf0a: warning: sum is close to minimum of "//
     +   "range")
	character*80 w193
	parameter (w193="ffxf0:  warning: cancellations among input p"//
     +   "arameters (import difference)")
	character*80 w194
	parameter (w194="ffxdbd: warning: cancellations in summing up")
	character*80 w195
	parameter (w195="ffdot6: warning: cancellations in dotproduct"//
     +   " s_i.s_{i+1}")
	character*80 w196
	parameter (w196="ffdot6: warning: cancellations in dotproduct"//
     +   " s_i.s_{i-1}")
	character*80 w197
	parameter (w197="ffdot6: warning: cancellations in dotproduct"//
     +   " p_i.s_i")
	character*80 w198
	parameter (w198="ffdot6: warning: cancellations in dotproduct"//
     +   " p_i.s_{i+1}")
	character*80 w199
	parameter (w199="ffdot6: warning: cancellations in dotproduct"//
     +   " p_{i-1}.s_i")
	character*80 w200
	parameter (w200="ffdot6: warning: cancellations in dotproduct"//
     +   " p_i.s_{i+2}")
	character*80 w201
	parameter (w201="ffdot6: warning: cancellations in dotproduct"//
     +   " p_{i+1}.s_i")
	character*80 w202
	parameter (w202="ffdot6: warning: cancellations in dotproduct"//
     +   " p_{i+2}.s_{i+1}")
	character*80 w203
	parameter (w203="ffdot6: warning: cancellations in dotproduct"//
     +   " p_i.p_{i+1}")
	character*80 w204
	parameter (w204="ffdot6: warning: cancellations in dotproduct"//
     +   " p_{i+1}.p_{i+2}")
	character*80 w205
	parameter (w205="ffdot6: warning: cancellations in dotproduct"//
     +   " p_{i+2}.p_i")
	character*80 w206
	parameter (w206="ffdot6: warning: cancellations in dotproduct"//
     +   " p_{i+2}.s_{i+2}")
	character*80 w207
	parameter (w207="ffdot6: warning: cancellations in dotproduct"//
     +   " s_i.s{i+3}")
	character*80 w208
	parameter (w208="ffdot6: warning: cancellations in dotproduct"//
     +   " pi.pj")
	character*80 w209
	parameter (w209="ffxdna: warning: cancellations in 1+/-a, une"//
     +   "xpected...")
	character*80 w210
	parameter (w210="ffxdna: warning: cancellations in b-a, unexp"//
     +   "ected...")
	character*80 w211
	parameter (w211="ffcd0c: warning: cancellations in subtractio"//
     +   "n of IR pole (to be expected)")
	character*80 w212
	parameter (w212="ffcd0c: warning: cancellations in computatio"//
     +   "n prop1 for threshold")
	character*80 w213
	parameter (w213="ffcd0c: warning: cancellations in computatio"//
     +   "n prop2 for threshold")
	character*80 w214
	parameter (w214="ffxb2a: warning: cancellations in B2d")
	character*80 w215
	parameter (w215="ffxd0p: warning: cancellations in complex de"//
     +   "l3mi")
	character*80 w216
	parameter (w216="ffzcnp: warning: cancellations in y (can be "//
     +   "fixed, contact author)")
	character*80 w217
	parameter (w217="ffzdnp: warning: cancellations in delta^(pi "//
     +   "si+1)_(pi pi+1)")
	character*80 w218
	parameter (w218="ffzdnp: warning: cancellations in (delta^(m"//
     +   "u si+1)_(pi pi+1))^2")
	character*80 w219
	parameter (w219="ffzcnp: warning: cancellations in z (can be "//
     +   "fixed, contact author)")
	character*80 w220
	parameter (w220="ffxb1:  warning: not enough terms in Taylor "//
     +   "expansion, may be serious")
	character*80 w221
	parameter (w221="ffxdb0: warning: cancellations in computatio"//
     +   "n 'diff'")
	character*80 w222
	parameter (w222="ffxdb0: warning: still cancellations is spli"//
     +   "t-up 1")
	character*80 w223
	parameter (w223="ffxdb0: warning: still cancellations is s1")
	character*80 w224
	parameter (w224="ffxdb0: warning: cancellations in B0', compl"//
     +   "ex args (can be improved)")
	character*80 w225
	parameter (w225="ffxb2p: warning: cancellations in B21 (after"//
     +   " a lot of effort)")
	character*80 w226
	parameter (w226="ffxb2p: warning: cancellations in B22")
	character*80 w227
	parameter (w227="ffxb2a: warning: cancellations in B21")
	character*80 w228
	parameter (w228="ffxbdp: warning: cancellations in case p^2=0")
	character*80 w229
	parameter (w229="ffxdpv: warning: cancellations in going from"//
     +   " delta- to PV-scheme")
	character*80 w230
	parameter (w230="ffxl22: warning: not enough terms in Taylor "//
     +   "expansion Li2(2-x)")
	character*80 w231
	parameter (w231="dfflo2: warning: not enough terms in taylor "//
     +   "expansion, using log(1-x)+x")
	character*80 w232
	parameter (w232="dfflo3: warning: not enough terms in taylor "//
     +   "expansion, using log(1-x)+x+x^2/2")
	character*80 w233
	parameter (w233="ffcdbp: warning: cancellations in equal mass"//
     +   "es case")
	character*80 w234
	parameter (w234="ffcbdp: warning: cancellations in case p^2=0")
	character*80 w235
	parameter (w235="ffcbdp: warning: cancellations in small diff.")
	character*80 w236
	parameter (w236="ffcbdp: warning: cancellations in 1-alpha")
	character*80 w237
	parameter (w237="ffcbdp: warning: cancellations in s2-alpha, "//
     +   "may not be serious")
	character*80 w238
	parameter (w238="ffcbdp: warning: not enough terms in Taylor "//
     +   "expansion, may be serious")
	character*80 w239
	parameter (w239="ffcbdp: warning: cancellations in s1-(1-alph"//
     +   "a), may not be serious")
	character*80 w240
	parameter (w240="ffcbdp: warning: cancellations in final resu"//
     +   "lt")
	character*80 w241
	parameter (w241="ffxe2:  warning: cancellations in E2 (can ma"//
     +   "ybe be done better)")
	character*80 w242
	parameter (w242="ffxe3:  warning: cancellations in E3 (can ma"//
     +   "ybe be done better)")
	character*80 w243
	parameter (w243="ffxe3:  warning: cancellations in adding det"//
     +   "erminants (may not be serious)")
	character*80 w244
	parameter (w244="ffcdna: warning: cancellations in del45")
	character*80 w245
	parameter (w245="ffcdna: warning: cancellations in del543m")
	character*80 w246
	parameter (w246="ffcdna: warning: cancellations in B")
	character*80 w247
	parameter (w247="ffcdna: warning: cancellations in C")
	character*80 w248
	parameter (w248="ffcdna: warning: cancellations between z1 an"//
     +   "d alpha")
	character*80 w249
	parameter (w249="ffcdna: warning: cancellations between z2 an"//
     +   "d alpha")
	character*80 w250
	parameter (w250="ffcdna: warning: cancellations in 1 + r*x1 ")
	character*80 w251
	parameter (w251="ffcdna: warning: cancellations in 1 + r*x2")
	character*80 w252
	parameter (w252="ffcdna: warning: cancellations between r*x1 "//
     +   "and r*x2")
	character*80 w253
	parameter (w253="ffd0c: warning: something wrong with the "//
     +   "rotation")
	character*80 w254
	parameter (w254="ffTn: warning: numerical cancellation "//
     +   "in in-triangle check")
	character*80 w255
	parameter (w255="ffRn: warning: 3-point Landau singularity")
	character*80 w256
	parameter (w256="ffRn: warning: Im(a.b) in the 1st theta "//
     +   "function is zero")
	character*80 w257
	parameter (w257="ffRn: warning: Im(a.b) in the 2nd theta "//
     +   "function is zero")
	character*80 w258
	parameter (w258="ffint3: cannot handle complex x yet")
	character*80 warn(258)
	data warn / w1,w2,w3,w4,w5,w6,w7,w8,w9,
     +   w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,
     +   w20,w21,w22,w23,w24,w25,w26,w27,w28,w29,
     +   w30,w31,w32,w33,w34,w35,w36,w37,w38,w39,
     +   w40,w41,w42,w43,w44,w45,w46,w47,w48,w49,
     +   w50,w51,w52,w53,w54,w55,w56,w57,w58,w59,
     +   w60,w61,w62,w63,w64,w65,w66,w67,w68,w69,
     +   w70,w71,w72,w73,w74,w75,w76,w77,w78,w79,
     +   w80,w81,w82,w83,w84,w85,w86,w87,w88,w89,
     +   w90,w91,w92,w93,w94,w95,w96,w97,w98,w99,
     +   w100,w101,w102,w103,w104,w105,w106,w107,w108,w109,
     +   w110,w111,w112,w113,w114,w115,w116,w117,w118,w119,
     +   w120,w121,w122,w123,w124,w125,w126,w127,w128,w129,
     +   w130,w131,w132,w133,w134,w135,w136,w137,w138,w139,
     +   w140,w141,w142,w143,w144,w145,w146,w147,w148,w149,
     +   w150,w151,w152,w153,w154,w155,w156,w157,w158,w159,
     +   w160,w161,w162,w163,w164,w165,w166,w167,w168,w169,
     +   w170,w171,w172,w173,w174,w175,w176,w177,w178,w179,
     +   w180,w181,w182,w183,w184,w185,w186,w187,w188,w189,
     +   w190,w191,w192,w193,w194,w195,w196,w197,w198,w199,
     +   w200,w201,w202,w203,w204,w205,w206,w207,w208,w209,
     +   w210,w211,w212,w213,w214,w215,w216,w217,w218,w219,
     +   w220,w221,w222,w223,w224,w225,w226,w227,w228,w229,
     +   w230,w231,w232,w233,w234,w235,w236,w237,w238,w239,
     +   w240,w241,w242,w243,w244,w245,w246,w247,w248,w249,
     +   w250,w251,w252,w253,w254,w255,w256,w257,w258 /
