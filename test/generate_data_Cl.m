(* calculate function with precision prec *)
GeneratePoint[f_, prec_, x_] :=
    { N[x, prec], Re[N[f[x], prec]] }

(* calculate functions with precision prec on a grid *)
GenerateGridData[f_, prec_, {min_, max_, step_}] :=
    Table[GeneratePoint[f, prec, x], {x, min, max, step}]

(* data arount {re, im} from given direction {reDir, imDir} *)
GenerateLimitData[f_, prec_, n_, frac_, x_, xDir_] :=
    Table[GeneratePoint[f, prec, x + xDir frac^k], {k, 1, n, 1}]

(* data arount {re, im} from all directions *)
GenerateLimitData[f_, prec_, n_, frac_, x_] :=
    Join[
        GenerateLimitData[f, prec, n, frac, x, +1],
        GenerateLimitData[f, prec, n, frac, x, -1]
    ]

ExportData[f_, prec_] :=
    Module[{data, filename},
           Print["Generating data for " <> ToString[f]];
           data = Join[
               GenerateGridData[f, prec, {0, 2 Pi, 1/100}],
               GenerateLimitData[f, prec, prec, 1/10, 0,    +1],
               GenerateLimitData[f, prec, prec, 1/10, Pi/2    ],
               GenerateLimitData[f, prec, prec, 1/10, Pi      ],
               GenerateLimitData[f, prec, prec, 1/10, 3Pi/2   ],
               GenerateLimitData[f, prec, prec, 1/10, 2Pi,  -1]
           ];
           filename = ToString[f] <> ".txt";
           Print["Writing data to ", filename];
           Export[filename, data, "Table"];
    ]


Cl2[x_] := ResourceFunction["ClausenCl"][2, x]
Cl3[x_] := ResourceFunction["ClausenCl"][3, x]
Cl4[x_] := ResourceFunction["ClausenCl"][4, x]
Cl5[x_] := ResourceFunction["ClausenCl"][5, x]
Cl6[x_] := ResourceFunction["ClausenCl"][6, x]

ExportData[Cl2, 40];
ExportData[Cl3, 40];
ExportData[Cl4, 40];
ExportData[Cl5, 40];
ExportData[Cl6, 40];
