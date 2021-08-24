(*

  This Mathematica script contains implementations of the Clausen[2,x]
  function in the intervals [0,Pi) and (0,2Pi) in terms of Bernoulli
  numbers [Abramowitz and Stegun, "Handbook of Mathematical Functions
  with Formulas, Graphs, and Mathematical Tables", 27.8.2-3].

 *)

(* x in [0, Pi) *)
Cl2Lo[x_, nMax_] :=
    x - x Log[x] + x^3/2 Sum[Abs[BernoulliB[2n + 2]] x^(2n)/((n + 1) Factorial[2n + 3]), {n,0,nMax}]

(* x in (0, 2Pi) *)
Cl2Hi[x_, nMax_] :=
    (Pi - x) (Log[2] - Sum[(2^(2n) - 1) Abs[BernoulliB[2n]] (Pi - x)^(2n)/(2n Factorial[2n + 1]), {n,1,nMax}])

Cl2[0, _] := 0

Cl2[Pi, _] := 0

Cl2[x_, nMax_] := -Cl2[-x, nMax] /; x < 0

Cl2[x_, nMax_] := Cl2[Mod[x, 2Pi], nMax] /; x >= 2Pi

Cl2[x_, nMax_] := -Cl2[2Pi - x, nMax] /; x > Pi

Cl2[x_, nMax_] := Cl2Lo[x, nMax] /; x > 0 && x < Pi/2

Cl2[x_, nMax_] := Cl2Hi[x, nMax] /; x >= Pi/2 && x < Pi
