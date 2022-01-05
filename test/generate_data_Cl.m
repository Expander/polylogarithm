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
               GenerateLimitData[f, prec, prec  , 1/10, 0,    +1],
               GenerateLimitData[f, prec, prec  , 1/10, Pi/2    ],
               GenerateLimitData[f, prec, prec  , 1/10, Pi      ],
               GenerateLimitData[f, prec, prec  , 1/10, 3Pi/2   ],
               GenerateLimitData[f, prec, prec-3, 1/10, 2Pi,  -1]
           ];
           filename = ToString[f] <> ".txt";
           Print["Writing data to ", filename];
           Export[filename, data, "Table"];
    ]


Cl2[x_] := ResourceFunction["ClausenCl"][2, x]
Cl3[x_] := ResourceFunction["ClausenCl"][3, x]
Cl4[x_] := ResourceFunction["ClausenCl"][4, x]
Cl5[x_] := ResourceFunction["ClausenCl"][5, N[x,40]]
Cl6[x_] := ResourceFunction["ClausenCl"][6, N[x,40]]
Cl7[x_] := ResourceFunction["ClausenCl"][7, N[x,40]]
Cl8[x_] := ResourceFunction["ClausenCl"][8, N[x,40]]
Cl9[x_] := ResourceFunction["ClausenCl"][9, N[x,40]]
Cl10[x_] := ResourceFunction["ClausenCl"][10, N[x,40]]
Cl11[x_] := ResourceFunction["ClausenCl"][11, N[x,40]]
Cl12[x_] := ResourceFunction["ClausenCl"][12, N[x,40]]
Cl13[x_] := ResourceFunction["ClausenCl"][13, N[x,40]]
Cl14[x_] := ResourceFunction["ClausenCl"][14, N[x,40]]
Cl15[x_] := ResourceFunction["ClausenCl"][15, N[x,40]]
Cl16[x_] := ResourceFunction["ClausenCl"][16, N[x,40]]
Cl1000[x_] := ResourceFunction["ClausenCl"][1000, N[x,40]]
Cl1001[x_] := ResourceFunction["ClausenCl"][1001, N[x,40]]
Cl1000000[x_] := ResourceFunction["ClausenCl"][1000000, N[x,40]]

ExportData[Cl2, 40];
ExportData[Cl3, 40];
ExportData[Cl4, 40];
ExportData[Cl5, 40];
ExportData[Cl6, 40];
ExportData[Cl7, 40];
ExportData[Cl8, 40];
ExportData[Cl9, 40];
ExportData[Cl10, 40];
ExportData[Cl11, 40];
ExportData[Cl12, 40];
ExportData[Cl13, 40];
ExportData[Cl14, 40];
ExportData[Cl15, 40];
ExportData[Cl16, 40];
ExportData[Cl1000, 40];
ExportData[Cl1001, 40];
ExportData[Cl1000000, 40];
