(*

  Implementation of the Clausen function Cl_n(theta) for n >= 2 from

  Jiming Wu, Xiaoping Zhang, Dongjie Liu: "An efficient calculation of
  the Clausen functions Cl_n(0)(n >= 2)"

 *)

(* Eq.(2.11) *)
Ncal[n_, theta_, nMax_] := 1/(n+1) (
    theta^(n+1)/(n+1)
    + Sum[(-1)^k BernoulliB[2k] theta^(2k+n+1)/((2k+n+1)Factorial[2k]), {k,1,nMax}]
)

Pcal[n_, theta_] :=
    Sum[(-1)^(Floor[(n-1)/2] + Floor[(i-1)/2]) theta^(n-i)/Factorial[n-i] Cl[i,0], {i,2,n}]

(* Eq.(1.4) *)
Cl[n_, 0, ___] := 0 /; EvenQ[n]

(* Eq.(1.4) *)
Cl[n_, 0, ___] := Zeta[n] /; OddQ[n]

(* Eq.(1.7), restrict theta to [0, Infinity) *)
Cl[n_, theta_, nMax_] :=
    (-1)^(n+1) Cl[n, -theta, nMax] /; theta < 0

(* Eq.(1.6), restrict theta to [0, Pi] *)
Cl[n_, theta_, nMax_] :=
    Cl[n, theta - Ceiling[theta/(2Pi)] 2Pi, nMax] /; theta > Pi

(* Eq.(2.8), (2.13) *)
Cl[n_, theta_, nMax_] := (
    (-1)^(Floor[(n+1)/2]) theta^(n-1)/Factorial[n-1] Log[2 Sin[theta/2]]
    + (-1)^(Floor[n/2] + 1)/Factorial[n-2] *
      Sum[(-1)^i Binomial[n-2,i] theta^i Ncal[n-2-i, theta, nMax], {i,0,n-2}]
    + Pcal[n, theta]
)
