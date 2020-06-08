(*
  Approximation of PolyLog[2,x]/x by an economized Pade approximation
  on the interval [0, 1/2).
 *)

Needs["FunctionApproximations`"];

(* maximum number of terms in Li2 series,
  so that [0,0.5]^nMax << machine precision *)
nMax = 200;

(* output precision *)
outPrec = 20;

interval = {0, 1/2};

half = (interval[[2]] - interval[[1]])/2;

(* series for PolyLog[2,x]/x *)
Li2x[x_, n_:nMax] := Sum[x^(i-1)/i^2, {i,1,n}];

(* bring rational function to standard form *)
PolynomialStandardForm[expr_, x_] :=
    Module[{n = Numerator[expr], d = Denominator[expr], c},
           c = d /. x -> 0;
           Expand[n/c] / Expand[d/c]
    ];

approx = EconomizedRationalApproximation[Li2x[x], {x, interval, 7, 7}] /. (-half + x) -> y;

approx = PolynomialStandardForm[approx, y];

FormatCoeffs[expr_, x_, prec_] :=
    N[#, prec]& /@ CoefficientList[expr, x]

Print["Numerator coefficients: ",
      FormatCoeffs[#,y,outPrec]& @ Numerator[approx]];

Print["Denominator coefficients: ",
      FormatCoeffs[#,y,outPrec]& @ Denominator[approx]];

(* calculate maximum error *)
diff = (approx - Li2x[x]) /. y -> x - half;

maxErr = Max @ Abs @ Table[N[diff, outPrec] /. x -> 1/2*n/10, {n, 0, 10}];

Print["max error: ", InputForm @ maxErr];
