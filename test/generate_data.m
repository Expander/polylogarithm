(* calculate polylogarithm of order s with precision prec *)
GeneratePoint[s_, prec_, re_, im_] :=
    Module[{val = PolyLog[s, re + I im]},
           N[#, prec]& @ {
               {re, im},
               {Re[val], Im[val]}
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

ExportData[s_, prec_] :=
    Module[{data, omega = 1/2 + Sqrt[3]/2, filename},
           Print["Generating data for Li" <> ToString[s]];
           data = Join[
               GenerateGridData[s, prec, {-5, 5, 1/10}],
               GenerateUnitData[s, prec, prec, 1/10, -1],
               GenerateUnitData[s, prec, prec, 1/10, +1],
               Join @@@ {
                   GeneratePoint[s, prec,  0  , 0],
                   GeneratePoint[s, prec,  1/2, 0],
                   GeneratePoint[s, prec,  1  , 0],
                   GeneratePoint[s, prec,  3/2, 0],
                   GeneratePoint[s, prec, -0  , 0],
                   GeneratePoint[s, prec, -1/2, 0],
                   GeneratePoint[s, prec, -1  , 0],
                   GeneratePoint[s, prec, -3/2, 0],
                   GeneratePoint[s, prec, -(Sqrt[5] - 1)/2, 0],
                   GeneratePoint[s, prec, -(Sqrt[5] + 1)/2, 0],
                   GeneratePoint[s, prec,  (Sqrt[5] + 1)/2, 0],
                   GeneratePoint[s, prec,  (Sqrt[5] + 3)/2, 0],
                   GeneratePoint[s, prec, omega, 0],
                   GeneratePoint[s, prec, omega^2, 0],
                   GeneratePoint[s, prec, 1 + omega, 0],
                   GeneratePoint[s, prec, 1/(1 + omega), 0]
               }
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
