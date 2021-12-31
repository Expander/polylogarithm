(*
  Approximations of Cl2[2,x] by rational function approximations
  on the intervals [0, Pi/2) and [Pi/2, Pi).
 *)

Needs["FunctionApproximations`"];

Get["../test/math/clausenBernoulli.m"];
Get["../test/math/clausenWu.m"];

(* maximum number of terms in Cl2 series *)
nMax = 100;

(* output precision *)
outPrec = 17;

(* bring rational function to standard form *)
PolynomialStandardForm[expr_, x_] :=
    Module[{n = Numerator[expr], d = Denominator[expr], c},
           c = d /. x -> 0;
           Expand[n/c] / Expand[d/c]
    ]

FormatCoeffs[expr_, x_, prec_] :=
    ScientificForm @ CForm @ N[#, prec]& /@ CoefficientList[expr, x]

CalcPade[fn_, interval_, {numTerms_, denTerms_}] :=
    Module[{half = interval[[1]] + (interval[[2]] - interval[[1]])/2,
            approx, diff, maxErr, y, range},
           approx = EconomizedRationalApproximation[fn[x], {x, interval, numTerms, denTerms}] /. x -> (y + half);
           approx = PolynomialStandardForm[approx, y];
           Print["Pade approximant for ", InputForm[fn], " on the interval ", InputForm[interval]];
           Print["Numerator coefficients: ",
                 FormatCoeffs[#,y,outPrec]& @ Numerator[approx]];
           Print["Denominator coefficients: ",
                 FormatCoeffs[#,y,outPrec]& @ Denominator[approx]];
           (* calculate maximum deviation *)
           Print["half interval: ", InputForm[half]];
           diff = (approx - fn[x]) /. y -> x - half;
           range = Range[interval[[1]] + 10^(-outPrec), interval[[2]], (interval[[2]] - interval[[1]])/100];
           maxErr = Max @ Abs @ N[diff /. x -> #, outPrec]& @ range;
           Print["max error: ", InputForm @ maxErr];
    ]

CalcMinimax[fn_, interval_, {numTerms_, denTerms_}] :=
    Module[{approx, maxErr},
           approx = MiniMaxApproximation[fn[x], {x, interval, numTerms, denTerms}, WorkingPrecision -> 1000];
           maxErr = approx[[2,2]];
           approx = PolynomialStandardForm[approx[[2,1]], x];
           Print["Approximant for ", InputForm[fn], " on the interval ", InputForm[interval]];
           Print["Numerator coefficients: ",
                 FormatCoeffs[#,x,outPrec]& @ Numerator[approx]];
           Print["Denominator coefficients: ",
                 FormatCoeffs[#,x,outPrec]& @ Denominator[approx]];
           (* calculate maximum deviation *)
           Print["max error: ", InputForm @ N[maxErr,6]];
           approx
    ]

Print["======================================================="];
Print[" Cl2 "];
Print["======================================================="];

(* interval = {0, Pi/2}; *)
fLo[x_] := Module[{y}, Expand[1/y^2 (Cl2Lo[y, nMax]/y - 1 + Log[y])] /. y -> x]

CalcMinimax[N[fLo[Sqrt[#]], 10*outPrec]&, {0, (Pi/2)^2}, {3, 3}];

(* interval = {Pi/2, Pi}; *)
fHi[x_] := Cl2Hi[x, nMax]/(Pi - x)

(* transformation *)
trans[x_] := Sqrt[x] + Pi

(* inverse transformation *)
itrans[x_] := (Pi - x)^2

CalcPade[N[fHi[trans[#]], 10*outPrec]&, itrans /@ {Pi/2, Pi}, {5,5}];

Print["======================================================="];
Print[" Cl3 "];
Print["======================================================="];

(* interval = {0, Pi/2}; *)
fLo[x_] := Module[{y}, Normal[Series[Expand[(Cl[3, y, nMax] - Zeta[3])/y^2 - Log[y]/2], {y,0,nMax}]] /. y -> x]

CalcMinimax[N[fLo[Sqrt[#]], 10*outPrec]&, {0, (Pi/2)^2}, {3,3}];

(* interval = {Pi/2, Pi}; *)
fHi[x_] := Cl[3, x, nMax]

CalcPade[N[fHi[trans[#]], 10*outPrec]&, itrans /@ {Pi/2, Pi}, {5,5}];

Print["======================================================="];
Print[" Cl4 "];
Print["======================================================="];

(* interval = {0, Pi/2}; *)
fLo[x_] := Module[{y}, Normal[Series[Expand[(Cl[4, y, nMax]/y - Zeta[3])/y^2 - Log[y]/6], {y,0,nMax}]] /. y -> x]

CalcMinimax[N[fLo[Sqrt[#]], 10*outPrec]&, {0, (Pi/2)^2}, {3,3}];

(* interval = {Pi/2, Pi}; *)
fHi[x_] := Cl[4, x, nMax]/(Pi - x)

CalcPade[N[fHi[trans[#]], 10*outPrec]&, itrans /@ {Pi/2, Pi}, {5,5}];

Print["======================================================="];
Print[" Cl5 "];
Print["======================================================="];

(* interval = {0, Pi/2}; *)
fLo[x_] := Module[{y}, Normal[Series[Expand[Cl[5, y, nMax] + 12 y^4 Log[y]/288], {y,0,nMax}]] /. y -> x]

CalcMinimax[N[fLo[Sqrt[#]], 10*outPrec]&, {0, (Pi/2)^2}, {3,4}];

(* interval = {Pi/2, Pi}; *)
fHi[x_] := Cl[5, x, nMax]

CalcPade[N[fHi[trans[#]], 10*outPrec]&, itrans /@ {Pi/2, Pi}, {5,5}];

Print["======================================================="];
Print[" Cl6 "];
Print["======================================================="];

(* interval = {0, Pi/2}; *)
fLo[x_] := Module[{y}, Normal[Series[Expand[Cl[6, y, nMax]/y + y^4 Log[y]/120], {y,0,nMax}]] /. y -> x]

CalcMinimax[N[fLo[Sqrt[#]], 10*outPrec]&, {0, (Pi/2)^2}, {3,3}];

(* interval = {Pi/2, Pi}; *)
fHi[x_] := Cl[6, x, nMax]/(Pi - x)

CalcPade[N[fHi[trans[#]], 10*outPrec]&, itrans /@ {Pi/2, Pi}, {4,5}];
