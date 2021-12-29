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

CalcPade[fn_, interval_, nTerms_] :=
    Module[{half = interval[[1]] + (interval[[2]] - interval[[1]])/2,
            approx, diff, maxErr, y, range},
           approx = EconomizedRationalApproximation[fn[x], {x, interval, nTerms, nTerms}] /. x -> (y + half);
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

CalcMinimax[fn_, interval_, nTerms_] :=
    Module[{approx, maxErr},
           approx = MiniMaxApproximation[fn[x], {x, interval, nTerms, nTerms}, WorkingPrecision -> 1000];
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

(* interval = {0, Pi/2}; *)
fLo[x_] := Module[{y}, Expand[1/y^2 (Cl2Lo[y, nMax]/y - 1 + Log[y])] /. y -> x]

CalcMinimax[N[fLo[Sqrt[#]], 10*outPrec]&, {0, (Pi/2)^2}, 3];

(* interval = {Pi/2, Pi}; *)
fHi[x_] := Cl2Hi[x, nMax]/(Pi - x)

(* transformation *)
trans[x_] := Sqrt[x] + Pi

(* inverse transformation *)
itrans[x_] := (Pi - x)^2

CalcPade[N[fHi[trans[#]], 10*outPrec]&, itrans /@ {Pi/2, Pi}, 5];

(* Cl3[x] *)

(* interval = {0, Pi/2}; *)
fLo[x_] := (Cl[3, x, 100] - Zeta[3])/x^2 - Log[x]/2

CalcPade[N[fLo[Sqrt[#]], 10*outPrec]&, {0, (Pi/2)^2}, 3];

(* interval = {Pi/2, Pi}; *)
fHi[x_] := Cl[3, x, 100]

CalcPade[N[fHi[trans[#]], 10*outPrec]&, itrans /@ {Pi/2, Pi}, 5];

(* Cl4[x] *)

(* interval = {0, Pi/2}; *)
fLo[x_] := (Cl[4, x, 100]/x - Zeta[3])/x^2 - Log[x]/6

CalcPade[N[fLo[Sqrt[#]], 10*outPrec]&, {0, (Pi/2)^2}, 3];

(* interval = {Pi/2, Pi}; *)
fHi[x_] := Cl[4, x, 100]/(Pi - x)

CalcPade[N[fHi[trans[#]], 10*outPrec]&, itrans /@ {Pi/2, Pi}, 5];
