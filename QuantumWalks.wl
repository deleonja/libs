(* ::Package:: *)

(* If ForScience paclet not installed, install it. See https://github.com/MMA-ForScience/ForScience *)
If[Length[PacletFind["ForScience"]]==0, PacletInstall[FileNameJoin[{DirectoryName[$InputFileName], "ForScience-0.88.45.paclet"}]]];


BeginPackage["QuantumWalks`"]


(* For nice formatting of usage messages, see https://github.com/MMA-ForScience/ForScience *)
<<ForScience`;


ClearAll[
  Shift, Coin, DTQWStep, DTQW, PositionProbabilityDistribution
]


Shift::usage = FormatUsage["Shift[t] yields a sparse array of the Shift operator for a 1D DTQW in an infinite line at time ```t```."];
Coin::usage = FormatUsage["Coin[t] yields a sparse array of the Haddamard Coin operator for a 1D DTQW in an infinite line at time ```t```.
Coin[t,C] yields a sparse array of the Coin operator ```C``` for a 1D DTQW in an infinite line at time ```t```."];
DTQWStep::usage = FormatUsage["DTQWStep[t] yields the unitary matrix of a Haddamard 1D DTQW in an infinite line at time ```t```.
DTQWStep[t,C] yields the unitary matrix of a 1D DTQW in an infinite line, using coin ```C```, at time ```t```."];
DTQW::usage = FormatUsage["DTQW[\[Psi]_0,t] yields the state at time ```t``` of a 1D Haddamard DTQW in an infinite line with initial state ```\[Psi]_0```.
DTQW[\[Psi]_0,t] yields the state at time ```t``` of a 1D DTQW in an infinite line, with coin ```C```, with initial state ```\[Psi]_0```."];
PositionProbabilityDistribution::usage = FormatUsage["PositionProbabilityDistribution[\[Psi],t] yields the position probability distribution of the state ```\[Psi]``` of a 1D DTQW at time ```t```."];
ExpValPosition::usage = FormatUsage["ExpValPosition[\[Psi],t] returns the expected value of position for the state \[Psi] of a 1D DTQW at time ```t```."];


0Begin["`Private`"]


Shift[t_] := Module[{},
  (* Check if t is an integer *)
  If[! IntegerQ[t], 
   Return[Message[Shift::intarg, t]]];
  
  (* Proceed with the original implementation *)
  SparseArray[
    Join[
      Table[{i, i + 2}, {i, Range[2, # - 2, 2]}], 
      Table[{i, i - 2}, {i, Range[3, #, 2]}]
    ] -> 1., {#, #}] &[2*(2*t + 1)]
]

(* Define a message for non-integer argument *)
Shift::intarg = "Argument `1` must be an integer (time).";


Coin::intarg = "Argument `1` must be an integer (time).";
Coin::matarg = "Argument `1` must be a matrix (coin argument).";

Coin[t_] := Module[{},
  (* Check if t is an integer *)
  If[! IntegerQ[t], 
   Return[Message[Coin::intarg, t]]];
  
  (* Proceed with the original implementation *)
  KroneckerProduct[IdentityMatrix[2 t + 1, SparseArray], SparseArray[FourierMatrix[2]]]
]

Coin[t_, c_] := Module[{},
  (* Check if t is an integer *)
  If[! IntegerQ[t], 
   Return[Message[Coin::intarg, t]]];
  
  (* Check if c is a matrix *)
  If[! MatrixQ[c], 
   Return[Message[Coin::matarg, c]]];
  
  (* Proceed with the original implementation *)
  KroneckerProduct[IdentityMatrix[2 t + 1, SparseArray], c]
]


DTQWStep[t_] := Shift[t] . Coin[t]
DTQWStep[t_, c_] := Shift[t] . Coin[t, c]
DTQWStep[t_, c_, psi_] := Chop[DTQWStep[t, c] . ArrayPad[psi, 2]]


DTQW[psi0_, t_] := Module[{psi},
  psi = ArrayPad[psi0, 2];
  Do[psi = ArrayPad[Chop[DTQWStep[i] . psi], 2], {i, t}];
  psi
]

DTQW[psi0_, t_, c_] := Module[{psi},
  psi = ArrayPad[psi0, 2];
  Do[psi = ArrayPad[Chop[DTQWStep[i, c] . psi], 2], {i, t}];
  psi
]


PositionProbabilityDistribution[psi_, tmax_] := Chop[
  Total[Abs[psi[[# ;; # + 1]]]^2] & /@ Range[1, 2*(2*tmax + 3), 2]
]


ExpValPosition[\[Psi]_,t_]:=PositionProbabilityDistribution[\[Psi],t] . Range[-t-1,t+1]


End[]


EndPackage[]
