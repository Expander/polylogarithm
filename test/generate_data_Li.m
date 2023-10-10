NEval[x_, prec_] := N[x, 3 prec]

(* calculate polylogarithm of order s with precision prec *)
GeneratePoint[s_, prec_, re_, im_] :=
    Module[{val = PolyLog[s, re + I im]},
           {
               NEval[#, prec]& @ {re, im},
               {Re[NEval[val, prec]], Im[NEval[val, prec]]}
           }
    ]

(* calculate polylogarithms of order s with precision prec on a
  grid *)
GenerateGridData[s_, prec_, {min_, max_, step_}] :=
    Flatten /@ Join @@ Table[
        GeneratePoint[s, prec, re, im],
        {re, min, max, step},
        {im, min, max, step}
    ]

(* calculate polylogarithms of order s with precision prec on a
  grid close to unity *)
GenerateUnitData[s_, prec_, n_, frac_, dir_] :=
    Flatten /@ Join @@ Table[
        GeneratePoint[s, prec, 1 + dir frac^k, im],
        {k, 1, n, 1},
        {im, 0, 0, 1}
    ]

(* data arount {re, im} from given direction {reDir, imDir} *)
GenerateLimitData[s_, prec_, n_, frac_, {re_, im_}, {reDir_, imDir_}] :=
    Flatten /@ Join @@ Table[
        GeneratePoint[s, prec, re + reDir frac^k, im + imDir frac^m],
        {k, 1, n, 1},
        {m, 1, n, 1}
    ]

(* data arount {re, im} from all directions *)
GenerateLimitData[s_, prec_, n_, frac_, {re_, im_}] :=
    Join[
        GenerateLimitData[s, prec, n, frac, {re, im}, {+1, +1}],
        GenerateLimitData[s, prec, n, frac, {re, im}, {+1, -1}],
        GenerateLimitData[s, prec, n, frac, {re, im}, {-1, +1}],
        GenerateLimitData[s, prec, n, frac, {re, im}, {-1, -1}]
    ]

ExportData[s_, prec_] :=
    Module[{data, omega = 1/2 + Sqrt[3]/2, filename},
           Print["Generating data for Li" <> ToString[s]];
           data = Join[
               GenerateGridData[s, prec, {-5, 5, 1/10}],
               GenerateGridData[s, prec, {-10, 30, 4/10}],
               GenerateGridData[s, prec, {-1000, 1000, 100}],
               GenerateUnitData[s, prec, prec, 1/10, -1],
               GenerateUnitData[s, prec, prec, 1/10, +1],
               GenerateLimitData[s, prec, prec, 1/10, {+1, 0}],
               GenerateLimitData[s, prec, prec, 1/10, { 0, 0}],
               GenerateLimitData[s, prec, prec, 1/10, {-1, 0}]
               (* Join @@@ { *)
               (*     GeneratePoint[s, prec,  0  , 0], *)
               (*     GeneratePoint[s, prec,  1/2, 0], *)
               (*     GeneratePoint[s, prec,  1  , 0], *)
               (*     GeneratePoint[s, prec,  3/2, 0], *)
               (*     GeneratePoint[s, prec, -0  , 0], *)
               (*     GeneratePoint[s, prec, -1/2, 0], *)
               (*     GeneratePoint[s, prec, -1  , 0], *)
               (*     GeneratePoint[s, prec, -3/2, 0], *)
               (*     GeneratePoint[s, prec, -(Sqrt[5] - 1)/2, 0], *)
               (*     GeneratePoint[s, prec, -(Sqrt[5] + 1)/2, 0], *)
               (*     GeneratePoint[s, prec,  (Sqrt[5] + 1)/2, 0], *)
               (*     GeneratePoint[s, prec,  (Sqrt[5] + 3)/2, 0], *)
               (*     GeneratePoint[s, prec, omega, 0], *)
               (*     GeneratePoint[s, prec, omega^2, 0], *)
               (*     GeneratePoint[s, prec, 1 + omega, 0], *)
               (*     GeneratePoint[s, prec, 1/(1 + omega), 0] *)
               (* } *)
           ];
           filename = "Li" <> ToString[s] <> ".txt";
           Print["Writing data to ", filename];
           Export[filename, data, "Table"];
    ]

ExportData[2, 40];
ExportData[3, 40];
ExportData[4, 40];
ExportData[5, 40];
ExportData[6, 40];

ExportData[1/2, 40];
ExportData[3/2, 40];
ExportData[5/2, 40];
ExportData[21/2, 40];
ExportData[41/2, 40];
