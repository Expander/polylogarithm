(*
  Approximation of PolyLog[4,x] by a rational function approximations.
 *)

Needs["FunctionApproximations`"];

(* maximum number of terms in the series *)
nMax = 200;

(* output precision *)
outPrec = 20;

(* series for PolyLog[4,x]/x *)
Li4x[x_, n_:nMax] := Sum[x^(i-1)/i^4, {i,1,n}];

(* bring rational function to standard form *)
PolynomialStandardForm[expr_, x_] :=
    Module[{n = Numerator[expr], d = Denominator[expr], c},
           c = d /. x -> 0;
           Expand[n/c] / Expand[d/c]
    ];

FormatCoeffs[expr_, x_, prec_] :=
    N[#, prec]& /@ CoefficientList[expr, x]

(* inversion formula rest term *)

(* -(2 Pi I)^n/Factorial[n] BernoulliB[n, 1/2 + Log[-x]/(2 Pi I)] /. n -> 4 *)

(* [-1,0] *)

Print["Approximation of Li4[x]/x on [-1,0]"];

li4 = Normal@Series[PolyLog[4,x]/x, {x,-1,100}]

approx = MiniMaxApproximation[li4, {x, {-1, 0}, 5, 6}, WorkingPrecision -> nMax];

maxErr = approx[[2,2]]

Print["max rel. error: ", InputForm @ maxErr];

approx = PolynomialStandardForm[approx[[2,1]], x];

Print["Numerator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Numerator[approx]];

Print["Denominator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Denominator[approx]];

(* [0,1/2] *)

Print["Approximation of Li4[x]/x on [0,1/2]"];

approx = MiniMaxApproximation[Li4x[x], {x, {0, 1/2}, 5, 5}, WorkingPrecision -> 100];

maxErr = approx[[2,2]]

Print["max rel. error: ", InputForm @ maxErr];

approx = PolynomialStandardForm[approx[[2,1]], x];

Print["Numerator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Numerator[approx]];

Print["Denominator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Denominator[approx]];

(* [1/2,8/10] *)

Print["Approximation of Li4[x] on [1/2,8/10]"];

approx = MiniMaxApproximation[PolyLog[4,x], {x, {1/2, 8/10}, 6, 6}, WorkingPrecision -> 100];

maxErr = approx[[2,2]]

Print["max rel. error: ", InputForm @ maxErr];

approx = PolynomialStandardForm[approx[[2,1]], x];

Print["Numerator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Numerator[approx]];

Print["Denominator coefficients: ",
      FormatCoeffs[#,x,outPrec]& @ Denominator[approx]];
