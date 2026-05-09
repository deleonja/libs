(* ::Package:: *)

BeginPackage["QMB`"];


Get["ForScience`"];


(* ::Section:: *)
(*Public definitions*)


Quiet[
SpectralFormFactorRMT::usage = FormatUsage[
"ConnectedSpectralFormFactor[\[Tau], \[Beta]] returns the connected spectral \
form factor K_c(\[Tau]) at the scaled time \[Tau]. \
The Dyson index \[Beta] determines the ensemble: 0 (Poisson), \
1 (GOE), 2 (GUE), or 4 (GSE)."
];
, {FrontEndObject::notavail, First::normal}];


(* ::Section:: *)
(*Private definitions*)


Begin["`RMT`Private`"];


ConnectedSpectralFormFactor::invBeta = "Dyson index \[Beta]=`1` is not supported. \
Valid values are 0 (Poisson), 1 (GOE), 2 (GUE), and 4 (GSE).";

(* Enables automatic vectorization over lists of scaled times *)
SetAttributes[SpectralFormFactorRMT, Listable];

SpectralFormFactorRMT[ScaledTime_, DysonIndex_Integer] := 
 Module[{AbsTau = Abs[ScaledTime]},
  Switch[DysonIndex,
   0, 
   1,
   
   1, 
   Piecewise[{
     {2 * AbsTau - AbsTau * Log[1 + 2 * AbsTau], AbsTau <= 1},
     {2 - AbsTau * Log[(2 * AbsTau + 1) / (2 * AbsTau - 1)], AbsTau > 1}
   }],
   
   2, 
   Piecewise[{
     {AbsTau, AbsTau <= 1},
     {1, AbsTau > 1}
   }],
   
   4, 
   Piecewise[{
     {AbsTau / 2 - (AbsTau / 4) * Log[Abs[1 - AbsTau]], AbsTau <= 2},
     {1, AbsTau > 2}
   }],
   
   _, 
   Message[SpectralFormFactorRMT::invBeta, DysonIndex]; 
   $Failed
  ]
 ]


End[];


EndPackage[];
