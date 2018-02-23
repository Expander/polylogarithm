u[z_] := -Log[1 - z];

X[0, n_] := BernoulliB[n];

X[p_, n_] :=
    Sum[Binomial[n, k] BernoulliB[n - k]/(k + 1) X[p - 1, k], {k, 0, n}];

(* Series expansion of polylogarithm Li_p(z) up to given order *)
Li[p_, z_, order_:50] :=
    Sum[X[p - 2, n] u[z]^(n + 1)/Factorial[n + 1], {n, 0, order}];

(* Remainder function on the r.h.s. of
   Li[p,z] + (-1)^p Li[p,1/z] = Rem[p,z] *)
Rem[p_, z_] :=
    -(2 Pi I)^p/Factorial[p] BernoulliB[p, 1/2 + Log[-z]/(2 Pi I)];
