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

appox = EconomizedRationalApproximation[Li2x[x], {x, interval, 7, 7}] /. (-half + x) -> y;

FormatCoeffs[expr_, x_, prec_] :=
    N[#, prec]& /@ CoefficientList[expr, x]

Print["Numerator coefficients: ",
      FormatCoeffs[#,y,outPrec]& @ Numerator[appox]];

Print["Denominator coefficients: ",
      FormatCoeffs[#,y,outPrec]& @ Denominator[appox]];

diff = (appox - Li2x[x]) /. y -> x - half;

maxErr = Max @ Abs @ Table[N[diff, outPrec] /. x -> 1/2*n/10, {n, 0, 10}];

Print["max error: ", InputForm @ maxErr];
