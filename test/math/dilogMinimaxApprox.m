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

approx = MiniMaxApproximation[Li2x[x], {x, interval, 5, 6}, WorkingPrecision -> 100];

maxErr = Max @ Cases[Abs[Li2x[#] - approx[[2,1]] /. x -> #]& /@ approx[[1]], Except[Indeterminate]];

Print["max error: ", InputForm @ maxErr];

approx = PolynomialStandardForm[approx[[2,1]], x];

FormatCoeffs[expr_, x_, prec_] :=
    N[#, prec]& /@ CoefficientList[expr, x]

Print["Numerator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Numerator[approx]];

Print["Denominator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Denominator[approx]];
