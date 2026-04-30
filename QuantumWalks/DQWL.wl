(* ::Package:: *)

BeginPackage["QuantumWalks`"]


(* For nice formatting of usage messages, see https://github.com/MMA-ForScience/ForScience *)
<<ForScience`;


(* ::Section::Closed:: *)
(*Public definitions*)


(* ::Subsection::Closed:: *)
(*DTQW*)


Shift::usage = FormatUsage[
"Shift[t] yields a sparse array of the Shift operator for a 1D DTQW in an infinite \
	line at time ```t```."
];

Coin::usage = FormatUsage[
"Coin[t] yields a sparse array of the Haddamard Coin operator for a 1D DTQW in an \
	infinite line at time ```t```.
Coin[t,coinOperator] yields a sparse array of the Coin operator ```coinOperator``` for a 1D DTQW in an \
	infinite line at time ```t```.
Coin[\[Alpha], dim] yields the coin operator of Alonso-Lobo (2025) for a position space \
	of size ```dim``` (Idenity \[CircleTimes] C).
Coin[matrix, dim] yields a global coin operator using ```matrix``` for a position \
	space of size ```dim``` (Identity \[CircleTimes] ```matrix```)."
];

Coin::intarg = "Argument `1` must be an integer (time or dimension).";
Coin::matarg = "Argument `1` must be a matrix (coin argument).";
Coin::invArg = "Invalid arguments for Coin.";

DTQWStep::usage = FormatUsage[
"DTQWStep[t] yields the unitary matrix of a Haddamard 1D DTQW in an infinite line at \
	time ```t```.
DTQWStep[t,coinOperator] yields the unitary matrix of a 1D DTQW in an infinite line, using coin \
	```coinOperator```, at time ```t```."];

DTQW::usage = FormatUsage[
"DTQW[\[Psi]_0,t] yields the state at time ```t``` of a 1D Haddamard DTQW in an infinite \
	line with initial state ```\[Psi]_0```.
DTQW[\[Psi]_0,t,coinOperator] yields the state at time ```t``` of a 1D DTQW in an infinite line, with \
	coin ```coinOperator```, with initial state ```\[Psi]_0```."];

PositionProbabilityDistribution::usage = FormatUsage[
"PositionProbabilityDistribution[\[Psi],t] yields the position probability distribution \
	of the state ```\[Psi]``` of a 1D DTQW at time ```t```."
];

ExpValPosition::usage = FormatUsage["ExpValPosition[\[Psi],t] returns the expected value \
	of position for the state \[Psi] of a 1D DTQW at time ```t```."
];


TransportVector::usage = FormatUsage[
"TransportVector[CoinMatrix] gives the transport vector for a DTQW given \
	the ```CoinMatrix```.
TransportVector[C1,C2] gives the transport vector for a \
	DTQW with step operator U = S.```C2```.S.```C1```
TransportVector[\[CurlyTheta],\[Phi],\[Theta]] gives the transport vector for a DTQW with a SU(2) \
	coin operator with axis {```\[CurlyTheta]```,```\[Phi]```} and angle of rotation ```\[Theta]```.
"];


(* ::Subsection::Closed:: *)
(*Parrondo's paradox*)


L::usage = FormatUsage[
	"LoosingStrategy[\[Theta], \[Theta]_a, \[Theta]_b] returns the loosing inequality."
];


W::usage = FormatUsage[
	"WinningStrategy[\[Theta], \[Theta]_a, \[Theta]_b] returns the winning inequality."
];


CriticalAngle::usage = FormatUsage[
	"CriticalAngle[avgPos] takes a list of sublists ```avgPos```, where each subslit \
	is of the form '''{\[Theta],\[LeftAngleBracket]x(\[Theta])\[RightAngleBracket]}''', and returns the value '''\[Theta]''' such that \
	'''\[LeftAngleBracket]x(\[Theta])\[RightAngleBracket]''' is the closest to zero of all sublsists."
];


(* ::Section::Closed:: *)
(*Private definitions*)


Begin["`DQWL`Private`"];


(* ::Subsection::Closed:: *)
(*DTQW*)


Shift[t_] := Shift[t] = Module[{},
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


Coin[t_Integer] := 
  KroneckerProduct[IdentityMatrix[2 t + 1, SparseArray], SparseArray[FourierMatrix[2]]]

(* Caso 1D: Coin Custom dependiente del tiempo *)
Coin[t_Integer, c_?MatrixQ] := 
  KroneckerProduct[IdentityMatrix[2 t + 1, SparseArray], c]

(* Caso de dimensi\[OAcute]n arbitraria: *)
Coin[\[Alpha]_?NumericQ, dim_Integer] := 
    KroneckerProduct[
        IdentityMatrix[dim, SparseArray], 
        SparseArray[{{Cos[\[Alpha]], Sin[\[Alpha]]}, {-Exp[I*Pi/4.]*Sin[\[Alpha]], Exp[I*Pi/4.]*Cos[\[Alpha]]}}]
    ];

(* Caso de dimensi\[OAcute]n arbitraria: *)
Coin[coinMatrix_?MatrixQ, dim_Integer] := 
    KroneckerProduct[
        IdentityMatrix[dim, SparseArray], 
        SparseArray[coinMatrix]
    ];
    
(* === VALIDACI\[CapitalOAcute]N Y MENSAJES DE ERROR (Fallbacks) === *)

(* Error: El argumento temporal no es entero *)
Coin[t_] /; !IntegerQ[t] := 
    (Message[Coin::intarg, t]; $Failed);

(* Error: En el caso 1D con moneda custom, el segundo argumento no es matriz *)
Coin[t_Integer, c_] /; !MatrixQ[c] := 
    (Message[Coin::matarg, c]; $Failed);

(* Error: En el caso 2D (Rotaci\[OAcute]n o Matriz), la dimensi\[OAcute]n no es entera *)
Coin[arg_, dim_] /; (NumericQ[arg] || MatrixQ[arg]) && !IntegerQ[dim] := 
    (Message[Coin::intarg, dim]; $Failed);

(* Catch-all: Cualquier otra combinaci\[OAcute]n inv\[AAcute]lida *)
Coin[___] := (Message[Coin::invArg, "Arguments not recognized. Check the usage \
	function with ?Coin"]; $Failed);


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


ExpValPosition[\[Psi]_, t_]:=PositionProbabilityDistribution[\[Psi],t] . Range[-t-1,t+1]


TransportVector[c_?MatrixQ] := Module[
  {k, s, uk, trU, omegaK, vgK, nVectorK, integrand, sigmas},
  
  (* Matrices de Pauli *)
  sigmas = {PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]};
  
  (* Operador de desplazamiento est\[AAcute]ndar en el espacio de momentos *)
  s = DiagonalMatrix[{Exp[-I*k], Exp[I*k]}];
  
  (* Operador de evoluci\[OAcute]n *)
  uk = s . c;
  
  (* Cuasi-energ\[IAcute]a *)
  trU = ComplexExpand[Re[Tr[uk]]];
  omegaK = ArcCos[trU / 2];
  
  (* Velocidad de grupo *)
  vgK = D[omegaK, k];
  
  (* Vector espectral unitario *)
  nVectorK = Table[
    (I / (2 * Sin[omegaK])) * Tr[uk . sigmas[[j]]],
    {j, 1, 3}
  ];
  
  (* Integrando *)
  integrand = vgK * nVectorK;
  
  (* Integraci\[OAcute]n num\[EAcute]rica *)
  Re[ 1/(2 Pi) * NIntegrate[integrand, {k, -Pi, Pi}, 
      Method -> "LocalAdaptive", 
      Exclusions -> {Sin[omegaK] == 0}] 
  ]
]

TransportVector[c1_?MatrixQ, c2_?MatrixQ] := Module[
  {k, s, uk, trU, omegaK, vgK, nVectorK, integrand, sigmas},
  
  (* Matrices de Pauli *)
  sigmas = {PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]};
  
  (* Operadores de desplazamiento en el espacio de Fourier (simb\[OAcute]licos en k) *)
  s = DiagonalMatrix[{Exp[-I*k], Exp[I*k]}];
  
  (* Operador de evoluci\[OAcute]n en el espacio de momentos *)
  uk = s . c2 . s . c1;
  
  (* Cuasi-energ\[IAcute]a (banda positiva) *)
  (* ComplexExpand y Re aseguran que no queden residuos imaginarios que confundan a ArcCos *)
  trU = ComplexExpand[Re[Tr[uk]]];
  omegaK = ArcCos[trU / 2];
  
  (* Velocidad de grupo obtenida por derivaci\[OAcute]n anal\[IAcute]tica exacta de omega *)
  vgK = D[omegaK, k];
  
  (* Vector espectral unitario *)
  (* Nota: Se usa 'nVectorK' para evitar el conflicto con el s\[IAcute]mbolo protegido 'N' *)
  nVectorK = Table[
    (I / (2 * Sin[omegaK])) * Tr[uk . sigmas[[j]]],
    {j, 1, 3}
  ];
  
  (* Integrando del vector de transporte *)
  integrand = vgK * nVectorK;
  
  (* Integraci\[OAcute]n num\[EAcute]rica sobre la zona de Brillouin *)
  (* Re limpia cualquier ruido num\[EAcute]rico imaginario del orden de $10^{-16}i$ *)
  Re[ 1/(2 Pi) * NIntegrate[integrand, {k, -Pi, Pi}, 
      Method -> "LocalAdaptive", 
      Exclusions -> {Sin[omegaK] == 0}] 
  ]
]

TransportVector[alpha_, phi_, theta_] := Module[
    {nVector, zVector, cosHalf, sinHalf, prefactorA, term1, term2, term3},
    
    (* 1. Definici\[OAcute]n de vectores base y eje de rotaci\[OAcute]n *)
    nVector = {Sin[alpha] * Cos[phi], Sin[alpha] * Sin[phi], Cos[alpha]};
    zVector = {0, 0, 1};
    
    (* 2. Pre-c\[AAcute]lculo de t\[EAcute]rminos trigonom\[EAcute]tricos de la moneda *)
    cosHalf = Cos[theta / 2];
    sinHalf = Sin[theta / 2];
    
    (* 3. C\[AAcute]lculo del prefactor escalar A *)
    prefactorA = (1 - sinHalf * Sin[alpha]) / (cosHalf^2 + sinHalf^2 * Cos[alpha]^2);
    
    (* 4. Construcci\[OAcute]n geom\[EAcute]trica del vector *)
    term1 = sinHalf^2 * Dot[nVector, zVector] * nVector;
    term2 = cosHalf^2 * zVector;
    term3 = sinHalf * cosHalf * Cross[zVector, nVector];
    
    (* Retorno simplificado *)
    Simplify[prefactorA * (term1 + term2 + term3)]
]


(* ::Subsection::Closed:: *)
(*Parrondo's paradox*)


L[\[Theta]_,\[Theta]a_,\[Theta]b_]:=\[Theta]a<=Mod[\[Theta],2.Pi]<=\[Theta]b


W[\[Theta]_,\[Theta]a_,\[Theta]b_]:=Mod[\[Theta],2.Pi]<=\[Theta]a||Mod[\[Theta],2.Pi]>=\[Theta]b


CriticalAngle[avgPos_] := 
	SortBy[
		Discard[avgPos, Round[#[[1]], 10.^-6] == Round[2Pi, 10.^-6] || Round[#[[1]], 10.^-6] == 0. &],
		Abs[#[[2]]]&
	][[1, 1]]


End[]


EndPackage[]
