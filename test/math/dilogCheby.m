(*

  This Mathematica script calculates the coefficients of the expansion
  of (-1)PolyLog[2,T] with T = -x in the interval [-1,0] in terms of
  Chebyshev polynomials.  It reproduces the results found in [Yudell
  L. Luke: Mathematical functions and their approximations, Academic
  Press Inc., New York 1975, p.67, Table 3.12].

 *)

ChebyshevCoefficients[n_Integer?Positive, f_Function] :=
    Module[{c, xk, coeffs},
           (* zeros of Chebyshev polynomials *)
           xk = Pi (Range[n] - 1/2)/n;
           (* Chebyshev coefficients *)
           c[j_] := 2/n Total[Cos[j xk] (f /@ Cos[xk])];
           (* coefficients in expansion of Chebyshev polynomials *)
           coeffs = Table[c[k], {k, 0, n-1}];
           coeffs[[1]] *= 1/2;
           coeffs
    ]

(* expansion interval [a,b] = [-1,0] *)
a = -1;
b = 0;

(* interval mapping y(x): [a,b] -> [-1,1] *)
y[x_] := (a + b - 2*x)/(a - b)

(* inverse interval mapping x(y): [-1,1] -> [a,b] *)
x[y_] := y * (b - a)/2 + (b + a)/2

(* function to approximate, # is in [-1,1] *)
f = -PolyLog[2, x[-#]]&

coeffs = ChebyshevCoefficients[45, f]

(* precision *)
prec = 37

Print @ DecimalForm[#, {prec, prec}]& @ Re @ N[#, prec]& @ coeffs
