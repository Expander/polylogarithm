Needs["FunctionApproximations`"];

(* maximum number of terms in Li2 series,
  so that [0,0.5]^nMax << machine precision *)
nMax = 200;

(* series for PolyLog[2,x]/x *)
Li2x[x_, n_:nMax] := Sum[x^(i-1)/i^2, {i,1,n}];

{abscissa, {appox, maxErr}} =
MiniMaxApproximation[Li2x[x], {x, {0,1/2}, 6, 6},
                     WorkingPrecision -> 20,
                     MaxIterations -> 1000];

Print["Numerator coefficients: ",
      CForm @ CoefficientList[#,x]& @ Numerator[appox]];

Print["Denominator coefficients: ",
      CForm @ CoefficientList[#,x]& @ Denominator[appox]];

Print["max error: ", InputForm @ maxErr];
