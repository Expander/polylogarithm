u[x_] := -Log[1 - x];
X[0, n_] := BernoulliB[n];
X[p_, n_] := Sum[Binomial[n, k] BernoulliB[n - k]/(k + 1) X[p - 1, k], {k, 0, n}];
Li[p_, x_, n_:100] := Sum[X[p - 2, n] u[x]^(n + 1)/Factorial[n + 1], {n, 0, n}];
Rem[p_,z_] := -(2 Pi I)^p/Factorial[p] BernoulliB[p, 1/2 + Log[-z]/(2 Pi I)];
