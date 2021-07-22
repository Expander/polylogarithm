(*

  This Mathematica script calculates the coefficients of the expansion
  of Clausen[2,x] in the intervals [-Pi/2,Pi/2] and [Pi/2,3Pi/2] in
  terms of Chebyshev polynomials.  It reproduces the results found in
  [K.S. Kölbig, Journal of Computational and Applied Mathematics 64
  (1995) 295-297].
 *)

ChebyshevCoefficients[n_Integer?Positive, f_Function] :=
    Module[{c, xk, coeffs},
           (* zeros of Chebyshev polynomials *)
           xk = Pi (Range[n] - 1/2)/n;
           (* Chebyshev coefficients *)
           c[j_] := 2/n Total[Cos[j xk] (f /@ Cos[xk])];
           (* coefficients in expansion of even Chebyshev polynomials *)
           coeffs = Table[c[k], {k, 0, n-1, 2}];
           coeffs[[1]] *= 1/2;
           coeffs
    ]

(* interval mapping y(x): [a,b] -> [-1,1] *)
y[x_, a_, b_] := (a + b - 2*x)/(a - b)

(* inverse interval mapping x(y): [-1,1] -> [a,b] *)
x[y_, a_, b_] := y * (b - a)/2 + (b + a)/2

(* Clausen function, n = 2 *)
Cl2[x_] := I/2 (PolyLog[2, E^(-I*x)] - PolyLog[2, E^(I*x)])

(* expansion interval [a,b] *)
a = -Pi/2;
b = Pi/2;

(* function to approximate on [a,b] *)
f1 = 2/x[#,a,b]^2 (Cl2[x[#,a,b]]/x[#,a,b] - 1 + Log[Abs[x[#,a,b]]])&

coeffs1 = ChebyshevCoefficients[38, f1]

(* expansion interval [a,b] *)
a = Pi/2;
b = 3 Pi/2;

(* function to approximate on [a,b] *)
f2 = Cl2[x[#,a,b]]/(Pi - x[#,a,b])&

coeffs2 = ChebyshevCoefficients[60, f2]

(* precision *)
prec = 37

Print @ DecimalForm[#, {prec, prec}]& @ Re @ N[#, prec]& @ coeffs1

Print @ DecimalForm[#, {prec, prec}]& @ Re @ N[#, prec]& @ coeffs2
