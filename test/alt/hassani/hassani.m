R[x_, N_] := Sum[(1/n^2 + 1/n Log[x])/x^n, {n,1,N}]

Spence[x_, N_:100] := -1/2 Log[x]^2 - Pi^2/6 + R[x, N]

r1 = { Li2[x_] :> -Li2[1/x] - Pi^2/6 - 1/2 Log[-x]^2 }

r2 = { Li2[x_] :> -Li2[1-x] + Pi^2/6 - Log[x] Log[1-x] }
