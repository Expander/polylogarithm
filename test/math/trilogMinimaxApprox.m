(*
  Approximation of PolyLog[3,x]/x by a rational function approximation
  on the intervals [0, 1/2] and [-1, 0].
 *)

Needs["FunctionApproximations`"];

(* maximum number of terms in Li3 series,
  so that [0,0.5]^nMax << machine precision *)
nMax = 200;

(* output precision *)
outPrec = 20;

interval = {0, 1/2};

(* series for PolyLog[3,x]/x *)
Li3x[x_, n_:nMax] := Sum[x^(i-1)/i^3, {i,1,n}];

(* bring rational function to standard form *)
PolynomialStandardForm[expr_, x_] :=
    Module[{n = Numerator[expr], d = Denominator[expr], c},
           c = d /. x -> 0;
           Expand[n/c] / Expand[d/c]
    ];

approx = MiniMaxApproximation[Li3x[x], {x, interval, 5, 6}, WorkingPrecision -> 100];

maxErr = approx[[2,2]]

Print["max rel. error: ", InputForm @ maxErr];

approx = PolynomialStandardForm[approx[[2,1]], x];

FormatCoeffs[expr_, x_, prec_] :=
    N[#, prec]& /@ CoefficientList[expr, x]

Print["Numerator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Numerator[approx]];

Print["Denominator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Denominator[approx]];

(* interval [-1,0] *)

interval = {-1, 0};

li3 = Normal@Series[PolyLog[3,x]/x, {x,-1,100}]

approx = MiniMaxApproximation[li3, {x, interval, 5, 6}, WorkingPrecision -> 100];

maxErr = approx[[2,2]]

Print["max rel. error: ", InputForm @ maxErr];

approx = PolynomialStandardForm[approx[[2,1]], x];

FormatCoeffs[expr_, x_, prec_] :=
    N[#, prec]& /@ CoefficientList[expr, x]

Print["Numerator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Numerator[approx]];

Print["Denominator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Denominator[approx]];
