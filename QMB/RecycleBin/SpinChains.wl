(* ::Package:: *)

(* ::Section:: *)
(*Public definitions*)


(* ::Subsection:: *)
(*Spin chains*)


(* ::Subsubsection::Closed:: *)
(*Symmetries*)


SpinParityEigenvectors::usage = "SpinParityEigenvectors[L] gives a list of {even, odd} eigenvectors of the L-spin system parity operator P; P\!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"1\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"L\"]}]},\n\"Ket\"]\) = \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"L\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"1\"]}]},\n\"Ket\"]\), \!\(\*SubscriptBox[\(k\), \(i\)]\)=0,1.";


TranslationEigenvectorRepresentatives::usage = FormatUsage[
  "TranslationEigenvectorRepresentatives[L] returns a list of sublists, each containing:\
  the decimal representation of a bit-string representative eigenvector,\
  its pseudomomentum k, and the length of its translation orbit,\
  for a system of ```L``` qubits."
];


BlockDiagonalize::usage = FormatUsage[
	"BlockDiagonalize[matrix,opts] returns ```matrix``` in block-diagonal form. Only \
	option is '''Symmetry'''."
];


Symmetry::usage = FormatUsage[
	"Symmetry is an option for '''BlockDiagonalize''' to specify the symmetry according to which \
	```matrix``` in '''BlockDiagonalize'''[```matrix```] is block diagonalized. It takes the values: \
	\"Translation\". Default option is \"Translation\"."
];


(* ::Subsubsection::Closed:: *)
(*Hamiltonians*)


IsingHamiltonian::usage = FormatUsage[
	"IsingHamiltonian[h_x,h_z,J,L,opts] returns the Hamiltonian \ 
	H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - ```J``` \[Sum]_{*i=1*}^{*L-1*} \[Sigma]^z_i \[Sigma]^z_{*i+1*} \
	with boundary conditions specified by option BoundaryConditions (default is \"Open\")."
];


BoundaryConditions::usage = FormatUsage[
	"BoundaryConditions is an option for '''IsingHamiltonian''' to specify the boundary conditions. It \
	takes the values \"Open\" or \"Periodic\". Default option is \"Open\"."
];


IsingNNOpenHamiltonian::replaced = "Function `1` has been replaced by `2`.";


(*Quiet[
IsingNNOpenHamiltonian::usage = FormatUsage["IsingNNOpenHamiltonian[h_x,h_z,J,L] returns the Hamiltonian H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - ```J```\[Sum]_{*i=1*}^{*L-1*} \[Sigma]^z_i \[Sigma]^z_{*i+1*}.
IsingNNOpenHamiltonian[h_x,h_z,{J_1,...,J_L},L] returns the Hamiltonian H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - \[Sum]_{*i=1*}^{*L-1*} ```J_i``` \[Sigma]^z_i \[Sigma]^z_{*i+1*}."];
, {FrontEndObject::notavail, First::normal}];*)


IsingNNClosedHamiltonian::replaced = "Function `1` has been replaced by `2`.";


(*IsingNNClosedHamiltonian::usage = "IsingNNClosedHamiltonian[\!\(\*
StyleBox[SubscriptBox[\"h\", \"x\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"h\", \"z\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"J\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)] returns the Hamiltonian \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\)(\!\(\*SubscriptBox[\(h\), \(x\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(x\)]\) + \!\(\*SubscriptBox[\(h\), \(z\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\)) + \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), L]\) \!\(\*SubscriptBox[\(J\), \(i\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\)\!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(z\)]\) with \!\(\*SubscriptBox[\(\[Sigma]\), \(L + 1\)]\) = \!\(\*SubscriptBox[\(\[Sigma]\), \(1\)]\).";*)


ClosedXXZHamiltonian::usage = "ClosedXXZHamiltonian[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[CapitalDelta]\",\nFontSlant->\"Italic\"]\)] returns the closed XXZ 1/2-spin chain as in appendix A.1 of Quantum 8, 1510 (2024).";


OpenXXZHamiltonian::usage= "OpenXXZHamiltonian[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[CapitalDelta]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"h1\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"h2\",\nFontSlant->\"Italic\"]\)] returns the open XXZ 1/2-spin chain as in appendix A.2 of Quantum 8, 1510 (2024).";


Quiet[
LeaSpinChainHamiltonian::usage = FormatUsage["LeaSpinChainHamiltonian[J_{*xy*},J_z,\[Omega],\[Epsilon]_d,L,d] returns the spin-1/2 chain H = \[Sum]_{*i=1*}^{*L-1*} ```J_{*xy*}```(S^x_i S^x_{*i+1*} + S^y_i S^y_{*i+1*}) + ```J_z```S^z_i S^z_{*i+1*} + \[Sum]_{*i=1*}^{*L*} ```\[Omega]``` S^z_i + \[Epsilon]_d S^z_d. [Eq. (1) in Am. J. Phys. 80, 246\[Dash]251 (2012)]."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
XXZOpenHamiltonian::usage = FormatUsage["XXZOpenHamiltonian[J_{*xy*},J_z,\[Omega],\[Epsilon]_d,L,d] returns the spin-1/2 chain \n H = \[Sum]_{*i=1*}^{*L-1*} ```J_{*xy*}```(S^x_i S^x_{*i+1*} + S^y_i S^y_{*i+1*}) + ```J_z```S^z_i S^z_{*i+1*} + \[Sum]_{*i=1*}^{*L*} ```\[Omega]``` S^z_i + \[Epsilon]_d S^z_d. \n [Eq. (1) in Am. J. Phys. 80, 246\[Dash]251 (2012)]."];
, {FrontEndObject::notavail, First::normal}];


HeisenbergXXXwNoise::usage="HeisenbergXXXwNoise[hz,L] returns the Heisenberg XXX spin 1/2 chain with noise: \!\(\*FormBox[\(H\\\  = \\\ \*FractionBox[\(1\), \(4\)]\\\ \(\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(L - 1\)]\\\ \((\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(x\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(x\)]\\\  + \\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(y\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(y\)]\\\  + \\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(z\)])\)\)\\\  + \\\ \*FractionBox[\(1\), \(2\)]\\\ \(\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\*SubsuperscriptBox[\(h\), \(i\), \(z\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\\\ \(\((open\\\ boundaries)\)\(.\)\)\)\),
TraditionalForm]\)";


QuantumGameOfLifeHamiltonian::usage = FormatUsage[
  "QuantumGameOfLifeHamiltonian[L, opts] returns the Hamiltonian for the Quantum Game of Life (QGL) \
  model for a spin-1/2 chain of length ```L``` (arXiv:2510.16570):
  H = \[Sum] \!\(\[Sigma]\_i\^x\) (\!\(\[ScriptN]\_i\^2\) + \!\(\[ScriptN]\_i\^3\)).
  
  Options:
  Momentum -> \"All\" (default) | Integer k (0 <= k < L).
  Parity -> \"None\" (default) | 1 | -1.
  
  Note: Parity projection (Inversion Symmetry) is only valid for Momentum sectors k=0 or k=L/2."
];


Momentum::usage = "Momentum is an option for QuantumGameOfLifeHamiltonian. \
It specifies the translation symmetry sector (integer k).";

Parity::usage = "Parity is an option for QuantumGameOfLifeHamiltonian. \
It specifies the inversion symmetry sector (+1 or -1).";


(* ::Subsection:: *)
(*Spin Hamiltonians*)


TwoSpinHamiltonian::usage = FormatUsage[
	"TwoSpinHamiltonian[s_1, s_2, \[Epsilon]_1, \[Epsilon]_2, \[Theta], g] returns the Hamiltonian \
	matrix for two spins of magnitude ```s_1``` and ```s_2```.
	The Hamiltonian is given by:
	H = ```\[Epsilon]_1``` S_1^z + ```\[Epsilon]_2``` (Cos[```\[Theta]```]S_2^z + Sin[```\[Theta]```]S_2^x) + \
	(```g```/Sqrt(s_1 s_2)) S_1 \[CenterDot] S_2."
];


(* ::Section:: *)
(*Private definitions*)


(* ::Subsection:: *)
(*Spin chains*)


(* ::Subsubsection:: *)
(*Symmetries*)


SpinParityEigenvectors[L_]:=Module[{tuples,nonPalindromes,palindromes},
tuples=Tuples[{0,1},L];
nonPalindromes=Select[tuples,#!=Reverse[#]&];
palindromes=Complement[tuples,nonPalindromes];
nonPalindromes=DeleteDuplicatesBy[nonPalindromes,Sort[{#,Reverse[#]}]&];
Normal[
{
Join[SparseArray[FromDigits[#,2]+1->1.,2^L]&/@palindromes,Normalize[SparseArray[{FromDigits[#,2]+1->1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes],
Normalize[SparseArray[{FromDigits[#,2]+1->-1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes
}]
]


Translation[state_, L_] := BitShiftRight[state] + BitAnd[state, 1]*2^(L-1)
BinaryNecklaces[L_Integer] := Module[{tuples=Tuples[{0,1},L]},
	Union[Table[First[Sort[NestList[RotateLeft, t, L-1]]], {t,tuples}]]
]
\[Omega][L_,k_] := Exp[2Pi*k*I/L]

TranslationEigenvectorRepresentatives[L_Integer] := Module[
	{
	necklaces = FromDigits[#,2]& /@ BinaryNecklaces[L],
	orbits
	},
	
	orbits = DeleteDuplicates[
	Sort[NestWhileList[Translation[#, L]&, #, UnsameQ[##]&, All][[;;-2]]]& /@ necklaces
	 ];
	Catenate[
		Outer[
			If[Mod[Length[#1]*#2, L]==0, {First[#1], #2, Length[#1]}, Nothing]&, 
			orbits, 
			Range[0, L-1], 
			1
		]
	]
]


RepresentativesOddBasis[basis_]:=DeleteDuplicatesBy[Discard[basis,PalindromeQ],Sort[{#,Reverse[#]}]&]
RepresentativesEvenBasis[basis_]:=DeleteDuplicatesBy[basis,Sort[{#,Reverse[#]}]&]


Options[BlockDiagonalize] = {
  Symmetry -> "Translation" (*coming soon: \"Parity\"*)
};

BlockDiagonalize[A_, opts:OptionsPattern[]]:= 
Switch[OptionValue[Symmetry],
	"Translation",
		Module[
			{
				L = Log[2, Length[A]],
				repsgathered,
				P
			},
				(*Gather translation eigenvectors by their pseudomomentum k*)
				repsgathered = GatherBy[TranslationEigenvectorRepresentatives[L], #[[2]]&];
				
				(*change of basis matrix*)
				P = SparseArray[
						Catenate[
							Apply[
								Normalize[SparseArray[
									Thread[
										FoldList[Translation[#, L]&, #1, Range[#3 - 1]] + 1 -> 
										Power[\[Omega][L, #2], Range[0., #3 - 1]]
									]
										, 2^L
								]]&,
								repsgathered,
								{2}
							]
						]
					];
					
				BlockDiagonalMatrix[Chop[Conjugate[P] . A . Transpose[P]]]
		],
	"Parity",
		Module[
		{
			L = Log[2, Length[A]],
			basis, reps, rules, Heven, Hodd
		},
			basis = Tuples[{0, 1}, L];
			rules=AssociationThread[basis->Range[Length[basis]]];
			reps=Comap[{RepresentativesEvenBasis, RepresentativesOddBasis},basis];
			Heven=1/2 (A[[#1,#1]] + A[[#1,#2]] + A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[1]]], {2}]);
			Hodd=1/2 (A[[#1,#1]] - A[[#1,#2]] - A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[2]]], {2}]);
			Heven = # . Heven . #&[DiagonalMatrix[ReplacePart[ConstantArray[1.,Length[reps[[1]]]],Thread[Catenate[Position[reps[[1]],_?PalindromeQ,1]]->1/Sqrt[2.]]],TargetStructure->"Sparse"]];
			
			{Heven, Hodd}
		],
	_,
		Message[BlockDiagonalize::badSymmetry, OptionValue[Symmetry]];
		Return[$Failed];
]

(*Mensaje de error si la opci\[OAcute]n es inv\[AAcute]lida*)
BlockDiagonalize::badSymmetry = 
  "Option badSymmetry -> `1` is not valid. Valid options: \"Translation\".";


(* ::Subsubsection:: *)
(*Hamiltonians*)


Options[IsingHamiltonian] = {
  BoundaryConditions -> "Open"
};

IsingHamiltonian[hx_, hz_, J_, L_, opts:OptionsPattern[]] := Module[
	{NNIndices},
	NNIndices = Switch[OptionValue[BoundaryConditions],
		"Open",
			Normal[SparseArray[Thread[{#,#+1}->3],{L}]&/@Range[L-1]],
		"Periodic",
			Normal[SparseArray[Thread[{#,Mod[#+1,L,1]}->3],{L}]&/@Range[L]],
		_,
			Message[
                IsingHamiltonian::badBoundaryCondition, 
                OptionValue[BoundaryConditions]
            ];
            Return[$Failed];
	];
	
	Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]
]

(*Mensaje de error si la opci\[OAcute]n es inv\[AAcute]lida*)
IsingHamiltonian::badBoundaryCondition = 
  "Option BoundaryConditions -> `1` not valid. Valid options: \"Open\" o \"Periodic\".";


IsingNNOpenHamiltonian[args___] := Message[
	IsingNNOpenHamiltonian::replaced, "IsingNNOpenHamiltonian", "IsingHamiltonian"];
	
(*IsingNNOpenHamiltonian[hx_,hz_,J_,L_] := Module[{NNIndices},
	NNIndices=Normal[SparseArray[Thread[{#,#+1}->3],{L}]&/@Range[L-1]];
	N[Normal[Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]]]]*)


IsingNNClosedHamiltonian[args___] := Message[
	IsingNNClosedHamiltonian::replaced, "IsingNNClosedHamiltonian", "IsingHamiltonian"];

(*IsingNNClosedHamiltonian[hx_,hz_,J_,L_] := 
Module[{NNIndices},
	NNIndices=Normal[SparseArray[Thread[{#,Mod[#+1,L,1]}->3],{L}]&/@Range[L]];
	Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]
]*)


ClosedXXZHamiltonian[L_,\[CapitalDelta]_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L] ];
		N[Normal[-1/2*Total[Join[Pauli/@NNindices,Pauli/@(2NNindices),\[CapitalDelta] (Pauli[#]-IdentityMatrix[2^L])&/@(3NNindices)]]]]
	]


OpenXXZHamiltonian[L_,\[CapitalDelta]_,h1_,h2_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L-1] ];
		N[Normal[-1/2*Total[Join[Pauli/@NNindices,Pauli/@(2NNindices),\[CapitalDelta] (Pauli[#]-IdentityMatrix[2^L])&/@(3NNindices)]]  
		- 1/2*(h1 Pauli[Join[{1},ConstantArray[0,L-1]]] + h2*Pauli[Join[ConstantArray[0,L-1],{1}]])+ 1/2*(h1 + h2)IdentityMatrix[2^L]]]
	]


HamiltonianNN[Jxy_,Jz_,L_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L-1] ];
		N[Normal[(1/4)*Total[Join[Jxy*(Pauli/@NNindices),Jxy*(Pauli/@(2NNindices)),Jz*(Pauli[#]&/@(3NNindices))]]]]
	]

HamiltonianZ[\[Omega]_,\[Epsilon]d_,L_,d_]:=N[(1/2)*(\[Omega]*Total[Pauli/@(3*IdentityMatrix[L])]+\[Epsilon]d*Pauli[Normal[SparseArray[d->3,L]]])]

LeaSpinChainHamiltonian[Jxy_,Jz_,\[Omega]_,\[Epsilon]d_,L_,d_]:=HamiltonianNN[Jxy,Jz,L]+HamiltonianZ[\[Omega],\[Epsilon]d,L,d]
XXZOpenHamiltonian[Jxy_,Jz_,\[Omega]_,\[Epsilon]d_,L_,d_]:=HamiltonianNN[Jxy,Jz,L]+HamiltonianZ[\[Omega],\[Epsilon]d,L,d]


HeisenbergXXXwNoise[h_List,L_]:=
Module[{NNIndices,firstSum,secondSum},
(* \sum_{k=1}^{L-1} S_k^xS_{k+1}^x + S_k^zS_{k+1}^z + S_k^zS_{k+1}^z *)
NNIndices=Normal[SparseArray[Thread[{#,#+1}->1],{L}]&/@Range[L-1]];
firstSum=1/4*Total[Table[Pauli[i*#]&/@NNIndices,{i,3}],2];

(* \sum_{k=1}^{L} h_k^z S_k^z *)
secondSum=1/2*h . (Pauli/@DiagonalMatrix[ConstantArray[3,L]]);

firstSum+secondSum
]


(* Quantum Game of Life *)

Options[QuantumGameOfLifeHamiltonian] = {
    Momentum -> "All",
    Parity -> "None" (* "None", 1, or -1 *)
};

QuantumGameOfLifeHamiltonian[L_Integer, OptionsPattern[]] := Module[
    {
        k = OptionValue[Momentum],
        p = OptionValue[Parity],
        dim = 2^L,
        Hfull,
        rules,
        basisVectors,
        P (* Projection Matrix *)
    },

    (* 1. VALIDACIONES *)
    If[k =!= "All" && !IntegerQ[k], Return[Message[QuantumGameOfLifeHamiltonian::invMom, k]; $Failed]];
    If[p =!= "None" && !MemberQ[{1, -1}, p], Return[Message[QuantumGameOfLifeHamiltonian::invPar, p]; $Failed]];
    
    If[IntegerQ[k] && p =!= "None",
        If[k != 0 && k != L/2,
            Message[QuantumGameOfLifeHamiltonian::parityWarning, k];
        ]
    ];

    (* 2. CONSTRUCCI\[CapitalOAcute]N DEL HAMILTONIANO COMPLETO (Bitwise) *)
    rules = Flatten @ Table[
        Module[{n = state, flips = {}, neighborsSum},
            Do[
                neighborsSum = 
                    BitGet[n, Mod[i - 2, L]] + 
                    BitGet[n, Mod[i - 1, L]] + 
                    BitGet[n, Mod[i + 1, L]] + 
                    BitGet[n, Mod[i + 2, L]];
                
                If[neighborsSum == 2 || neighborsSum == 3,
                    AppendTo[flips, {n + 1, BitXor[n, 2^i] + 1} -> 1.0]
                ],
                {i, 0, L - 1}
            ];
            flips
        ],
        {state, 0, dim - 1}
    ];
    
    Hfull = SparseArray[rules, {dim, dim}];

    (* 3. RETORNO SI NO HAY SIMETR\[CapitalIAcute]AS *)
    If[k === "All", Return[Hfull]];

    (* 4. CONSTRUCCI\[CapitalOAcute]N DE LA BASE DE MOMENTO k *)
    basisVectors = Module[{necklaces, kBasis},
        (* Agrupar por collares *)
        necklaces = GroupBy[Range[0, dim - 1], 
            Function[x, Min[NestList[RotateLeftBits[#, L]&, x, L-1]]]
        ];
        
        kBasis = Reap[
            Do[
                (* AQU\[CapitalIAcute] ESTABA EL ERROR: Elimin\[EAcute] 'T_k' y limpi\[EAcute] variables *)
                Module[{orbit = members, R, vec},
                    R = Length[orbit];
                    
                    (* Condici\[OAcute]n de compatibilidad de momento *)
                    If[Divisible[k * R, L],
                        vec = SparseArray[{}, dim];
                        Do[
                            (* Construcci\[OAcute]n del vector de Bloch *)
                            vec = vec + SparseArray[{orbit[[j+1]] + 1} -> Exp[-I * 2. * Pi * k * j / L], dim];
                        , {j, 0, R-1}];
                        
                        Sow[Normalize[vec]];
                    ]
                ],
                {members, Values[necklaces]}
            ]
        ][[2]];
        
        If[kBasis === {}, Return[SparseArray[{}, {0, 0}]]];
        Flatten[kBasis, 1]
    ];

    (* 5. PROYECCI\[CapitalOAcute]N DE PARIDAD *)
    If[p =!= "None",
        basisVectors = Reap[
            Scan[Function[v,
                Module[{vInv, vProj},
                    vInv = SparseArray[
                        Thread[(ReverseBits[# - 1, L] + 1 & /@ v["NonzeroPositions"][[All, 1]]) -> v["NonzeroValues"]], 
                        dim
                    ];
                    
                    vProj = v + p * vInv;
                    
                    If[Norm[vProj] > 10^-10, Sow[Normalize[vProj]]];
                ]
            ], basisVectors]
        ][[2]];
        
        If[basisVectors =!= {},
             (* Eliminaci\[OAcute]n de duplicados num\[EAcute]ricos *)
             basisVectors = DeleteDuplicates[Flatten[basisVectors, 1], (Norm[#1 - #2] < 10^-8 || Norm[#1 + #2] < 10^-8) &];
        ,
             Return[SparseArray[{}, {0, 0}]]
        ];
    ];

    (* 6. CONSTRUCCI\[CapitalOAcute]N FINAL *)
    P = Transpose[SparseArray[basisVectors]];
    Chop[ConjugateTranspose[P] . Hfull . P]
];

(* Mensajes y Funciones Auxiliares *)
QuantumGameOfLifeHamiltonian::invMom = "Momentum `1` must be \"All\" or an integer 0 <= k < L.";
QuantumGameOfLifeHamiltonian::invPar = "Parity `1` must be \"None\", 1, or -1.";
QuantumGameOfLifeHamiltonian::parityWarning = "Warning: Inversion Parity is generally not a good quantum number for momentum k=`1` (unless k=0 or k=L/2).";

RotateLeftBits[n_, L_] := BitOr[BitShiftLeft[BitAnd[n, 2^(L-1)-1], 1], BitShiftRight[n, L-1]];
ReverseBits[n_, L_] := FromDigits[Reverse[IntegerDigits[n, 2, L]], 2];


(* ::Subsection::Closed:: *)
(*Spin Hamiltonians*)


(* Helper function to generate sparse spin operators for arbitrary S *)
GenerateSpinOperators[s_] := Module[
    {d, range, diagonal, upperDiag, Sz, Sp, Sm, Sx, Sy, Id},
    
    (* Dimension of the Hilbert space for spin s *)
    d = Round[2 s + 1];
    
    (* Validate s is integer or half-integer *)
    If[Abs[d - (2 s + 1)] > 10^-10, Return[$Failed]];

    (* Basis range from S to -S *)
    range = Range[s, -s, -1];

    (* S_z: Diagonal matrix with m values *)
    Sz = SparseArray[Band[{1, 1}] -> range, {d, d}];

    (* S_plus: Elements Sqrt[s(s+1) - m(m+1)] on superdiagonal *)
    (* Note: range contains 'm'. The element <m+1|S+|m> is at position corresponding to transition m -> m+1 *)
    (* In matrix indices 1..d, index i is m_i. i-1 is m_i + 1. So it acts on column i to row i-1 *)
    upperDiag = Table[
        Sqrt[(s - m) (s + m + 1)], 
        {m, range[[2 ;;]]} (* Exclude the first m=s as it cannot be raised from *)
    ];
    Sp = SparseArray[Band[{1, 2}] -> upperDiag, {d, d}];
    
    (* S_minus: Transpose of S_plus *)
    Sm = Transpose[Sp];

    (* S_x and S_y *)
    Sx = (Sp + Sm) / 2;
    (* Sy = (Sp - Sm) / (2 I); Unused here but standard definition *)
    
    (* Identity *)
    Id = IdentityMatrix[d, SparseArray];

    <| "z" -> Sz, "x" -> Sx, "+" -> Sp, "-" -> Sm, "Id" -> Id |>
];

(* Main Hamiltonian Function *)
TwoSpinHamiltonian[s1_, s2_, \[Epsilon]1_, \[Epsilon]2_, \[Theta]_, g_] := Module[
    {
        ops1, ops2,
        term1, term2, interaction,
        H
    },
    
    (* 1. Generate local operators *)
    ops1 = GenerateSpinOperators[s1];
    ops2 = GenerateSpinOperators[s2];
    
    If[FailureQ[ops1] || FailureQ[ops2], 
        Return[Message[TwoSpinHamiltonian::invalidSpin, "{s1, s2}"]]
    ];

    (* 2. Term 1: \[Epsilon]1 * S1_z (x) Id_2 *)
    term1 = \[Epsilon]1 * KroneckerProduct[ops1["z"], ops2["Id"]];

    (* 3. Term 2: \[Epsilon]2 * Id_1 (x) (Cos[\[Theta]] S2_z + Sin[\[Theta]] S2_x) *)
    term2 = \[Epsilon]2 * KroneckerProduct[
        ops1["Id"], 
        Cos[\[Theta]] * ops2["z"] + Sin[\[Theta]] * ops2["x"]
    ];

    (* 4. Interaction Term: (g / Sqrt[S1 S2]) * (S1 . S2) *)
    (* Decomposed as SzSz + 1/2 (S+S- + S-S+) *)
    interaction = (g / Sqrt[s1 * s2]) * (
        KroneckerProduct[ops1["z"], ops2["z"]] + 
        0.5 * (
            KroneckerProduct[ops1["+"], ops2["-"]] + 
            KroneckerProduct[ops1["-"], ops2["+"]]
        )
    );

    (* 5. Total Sum *)
    H = term1 + term2 + interaction;

    (* Return SparseArray, optionally chopped if purely numerical *)
    If[AllTrue[{s1, s2, \[Epsilon]1, \[Epsilon]2, \[Theta], g}, NumericQ],
        Chop[H],
        H
    ]
];

TwoSpinHamiltonian::invalidSpin = "Spin values `1` must be integers or half-integers.";
