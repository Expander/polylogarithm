(* Functional properties of the dilogarithm *)

r1 = { Li2[z_] :> -Li2[1/z] - Pi^2/6 - 1/2 Log[-z]^2 }

r2 = { Li2[z_] :> -Li2[1-z] + Pi^2/6 - Log[z] Log[1-z] }
