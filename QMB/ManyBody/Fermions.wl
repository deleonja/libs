(* ::Package:: *)

BeginPackage["QMB`ManyBody`", {"QMB`"}];


Get["ForScience`"];


(* ::Section:: *)
(*Public definitions*)


Quiet[

SYKHamiltonian::usage = FormatUsage[
"SYKHamiltonian[\[ScriptCapitalN],J] yields the Hamiltonian of the q=4 Sachdev-Ye-Kitaev \
	(SYK) model for ```\[ScriptCapitalN]``` Majorana fermions and coupling variance J. 
SYKHamiltonian '''Options''': ParitySector."
];

MajoranaFermionOperator::usage = FormatUsage[
"MajoranaFermionOperator[Index,NumSpins] yields the matrix representation \
	of the Majorana fermion operator via the Jordan-Wigner transformation."
];

, {FrontEndObject::notavail, First::normal}];


(* ::Section:: *)
(*Private definitions*)


Begin["`Fermions`Private`"];


(* Memoized function to compute and cache symmetry sector indices *)
GetParitySectorIndices[NumSpins_Integer, Parity_Integer] := 
    GetParitySectorIndices[NumSpins, Parity] = Module[
        {Dim, AllIndices, BitCounts, SectorMask},
        
        Dim = 2^NumSpins;
        AllIndices = Range[Dim];
        
        (* Vectorized computation of Hamming weights *)
        BitCounts = Total[IntegerDigits[AllIndices - 1, 2, NumSpins], {2}];
        
        (* Mask creation: Parity == 1 implies Even, otherwise Odd *)
        SectorMask = If[Parity == 1,
            EvenQ[BitCounts],
            Not /@ EvenQ[BitCounts]
        ];
        
        (* Return the filtered indices *)
        Pick[AllIndices, SectorMask]
    ];


(* Map a Majorana index to its corresponding Pauli string via Jordan-Wigner *)
MajoranaFermionOperator[Index_Integer, NumSpins_Integer] := Module[
    {SiteIndex, IsEven, PauliList},
    
    SiteIndex = Ceiling[Index / 2];
    IsEven = EvenQ[Index];
    
    (* Pauli lists use QMB GeneralQM conventions: 0->I, 1->X, 2->Y, 3->Z *)
    PauliList = Join[
        ConstantArray[3, SiteIndex - 1],
        {If[IsEven, 2, 1]},
        ConstantArray[0, NumSpins - SiteIndex]
    ];
    
    (* Normalize to {chi_i, chi_j} = delta_ij *)
    (1.0 / Sqrt[2.0]) * Pauli[PauliList]
];



(*Options[SYKHamiltonian] = {
  ParitySector -> "None"
};

SYKHamiltonian[NumMajoranas_Integer, CouplingJ_?NumericQ, OptionsPattern[]] := Module[
    {Parity, NumSpins, Majoranas, Subsets4, Couplings, 
     HamiltonianTerm, FullHamiltonian, Dim, AllIndices, BitCounts, SectorMask, SectorIndices},
    
    Parity = OptionValue[ParitySector];
    NumSpins = NumMajoranas / 2;
    
    If[!IntegerQ[NumSpins],
        Print["Error: NumMajoranas must be an even integer."];
        Return[$Failed];
    ];
    
    (* 1. Pre-compute Majorana operators *)
    Majoranas = AssociationMap[
        MajoranaFermionOperator[#, NumSpins] &, 
        Range[NumMajoranas]
    ];
    
    (* 2. Generate 4-body interactions and sample couplings *)
    Subsets4 = Subsets[Range[NumMajoranas], {4}];
    Couplings = RandomVariate[
        NormalDistribution[0, CouplingJ * Sqrt[6.0 / (NumMajoranas^3)]], 
        Length[Subsets4]
    ];
    
    HamiltonianTerm = Function[{SubsetIndices, Jijkl},
        Jijkl * (Majoranas[SubsetIndices[[1]]] . Majoranas[SubsetIndices[[2]]] . 
                 Majoranas[SubsetIndices[[3]]] . Majoranas[SubsetIndices[[4]]])
    ];
    
    (* 3. Compute full sparse Hamiltonian *)
    FullHamiltonian = Chop @ Total @ MapThread[HamiltonianTerm, {Subsets4, Couplings}];
    
    (* 4. Return full matrix if no specific sector is requested *)
    If[Parity === "None",
        Return[FullHamiltonian]
    ];
    
    (* 5. Retrieve memoized sector indices and extract submatrix in O(1) *)
    FullHamiltonian[[#, #]] &[GetParitySectorIndices[NumSpins, Parity]]
]*)


MajoranaOperatorTable[NumSpins_Integer] :=
    MajoranaOperatorTable[NumSpins] = Module[
        {
            PauliX, PauliY, PauliZ, Id2, IdFull,
            ZPrefix, LocalX, LocalY, TrailingId,
            Prefactor, Result, k
        },

        (* Fixed 2x2 matrices; avoid repeated symbolic lookup. *)
        PauliX   = Pauli[1];
        PauliY   = Pauli[2];
        PauliZ   = Pauli[3];
        Id2      = Pauli[0];
        IdFull   = Pauli[ConstantArray[0, NumSpins]];
        Prefactor = 1.0 / Sqrt[2.0];

        Result = Association[];

        (* ZPrefix[[k]] = Z^{ox k} as a 2^k x 2^k matrix.               *)
        (* We grow it site by site with a single KroneckerProduct.        *)
        ZPrefix = {{{1.}}}; (* ZPrefix[[1]] is the 1x1 identity, i.e.    *)
                            (* no Z-string before site 1.                 *)

        Do[
            (* Z-string for sites 1..k-1, then local X or Y at site k,   *)
            (* then identity on remaining NumSpins-k sites.               *)
            LocalX  = KroneckerProduct[ZPrefix[[k]], PauliX,
                          IdentityMatrix[2^(NumSpins - k)]];
            LocalY  = KroneckerProduct[ZPrefix[[k]], PauliY,
                          IdentityMatrix[2^(NumSpins - k)]];

            Result[2k - 1] = Prefactor * LocalX;
            Result[2k]     = Prefactor * LocalY;

            (* Extend the Z-string by one factor for the next iteration.  *)
            AppendTo[ZPrefix, KroneckerProduct[ZPrefix[[k]], PauliZ]];
            ,
            {k, 1, NumSpins}
        ];

        Result
    ];


(* \[HorizontalLine]\[HorizontalLine] Step 2: Cached list of all 4-element index subsets \[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine] *)
(* Subsets[Range[N],{4}] is purely combinatorial; independent of J.          *)
FourBodySubsets[NumMajoranas_Integer] :=
    FourBodySubsets[NumMajoranas] = Subsets[Range[NumMajoranas], {4}];


(* \[HorizontalLine]\[HorizontalLine] Step 3: Cached operator basis \[LongDash] one matrix per 4-body term \[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine]\[HorizontalLine] *)
(* Each entry is gamma_i.gamma_j.gamma_k.gamma_l as a concrete sparse        *)
(* matrix. This is the expensive step (Binomial[N,4] matrix products),       *)
(* but it is fully deterministic and cached after the first call.            *)
(* On subsequent calls with the same N, the Hamiltonian is just a linear     *)
(* combination of these stored matrices, costing O(Binomial[N,4] * D^2).     *)
FourBodyBasisTerms[NumMajoranas_Integer] :=
    FourBodyBasisTerms[NumMajoranas] = Module[
        {NumSpins, Ops, Subsets4},

        NumSpins = NumMajoranas / 2;
        Ops      = MajoranaOperatorTable[NumSpins];
        Subsets4 = FourBodySubsets[NumMajoranas];

        Map[
            Function[S,
                Ops[S[[1]]] . Ops[S[[2]]] . Ops[S[[3]]] . Ops[S[[4]]]
            ],
            Subsets4
        ]
    ];


Options[SYKHamiltonian] = {
  ParitySector -> "None"
};

ClearAll[CompiledSYKSum];

CompiledSYKSum = Compile[
  {
    {coeffs, _Real, 1},
    {mats,   _Complex, 3}
  },
  Module[{res, i, n},
    n   = Length[coeffs];
    res = ConstantArray[0. + 0. I, Dimensions[mats[[1]]]];
    
    For[i = 1, i <= n, i++,
      res += coeffs[[i]] * mats[[i]];
    ];
    
    res
  ],
  CompilationTarget -> "C",
  RuntimeOptions -> "Speed"
];

(*SYKHamiltonian[NumMajoranas_Integer, CouplingJ_?NumericQ, OptionsPattern[]] :=
    Module[
        {Parity, NumSpins, BasisTerms, Sigma, Couplings, FullHamiltonian},

        Parity   = OptionValue[ParitySector];
        NumSpins = NumMajoranas / 2;

        If[!IntegerQ[NumSpins],
            Message[SYKHamiltonian::evenN];
            Return[$Failed]
        ];
		
		t1=AbsoluteTime[];
        (* Steps 1-3: O(1) after the first call for a given NumMajoranas.    *)
        BasisTerms = FourBodyBasisTerms[NumMajoranas];
        t2=AbsoluteTime[];

        (* Sample fresh couplings \[LongDash] the only source of randomness.           *)
        Sigma      = CouplingJ * Sqrt[6.0 / NumMajoranas^3];
        Couplings  = RandomVariate[
            NormalDistribution[0, Sigma],
            Length[BasisTerms]
        ];
        t3=AbsoluteTime[];

        (* Assemble H = sum_ijkl J_ijkl * (gamma_i gamma_j gamma_k gamma_l) *)
        FullHamiltonian = Chop @ Total[Couplings * BasisTerms];
        t4=AbsoluteTime[];

        (* Return the full matrix or the requested parity sector.            *)
        If[Parity === "None",
            Return[FullHamiltonian]
        ];

        FullHamiltonian[[#, #]] &[GetParitySectorIndices[NumSpins, Parity]]
    ];*)
SYKHamiltonian[NumMajoranas_Integer, CouplingJ_?NumericQ, OptionsPattern[]] :=
    Module[
        {Parity, NumSpins, BasisTerms, Sigma, Couplings, FullHamiltonian},

        Parity   = OptionValue[ParitySector];
        NumSpins = NumMajoranas / 2;

        If[!IntegerQ[NumSpins],
            Message[SYKHamiltonian::evenN];
            Return[$Failed]
        ];
		
		t1=AbsoluteTime[];
        (* Steps 1-3: O(1) after the first call for a given NumMajoranas.    *)
        (*BasisTerms = Developer`ToPackedArray @ N @ Normal @ FourBodyBasisTerms[NumMajoranas];*)
        BasisTerms = FourBodyBasisTerms[NumMajoranas];
        t2=AbsoluteTime[];

        (* Sample fresh couplings \[LongDash] the only source of randomness.           *)
        Sigma      = CouplingJ * Sqrt[6.0 / NumMajoranas^3];
        Couplings  = RandomVariate[
            NormalDistribution[0, Sigma],
            Length[BasisTerms]
        ];
        t3=AbsoluteTime[];

        (* Assemble H = sum_ijkl J_ijkl * (gamma_i gamma_j gamma_k gamma_l) *)
        (*FullHamiltonian = Chop @ CompiledSYKSum[Couplings, BasisTerms];*)
        (*FullHamiltonian = SparseArray[{}, Dimensions[BasisTerms[[1]]]];Print[Length[BasisTerms]]

		Do[
			FullHamiltonian += Couplings[[i]] * BasisTerms[[i]],
			{i, Length[Couplings]}
		];

		FullHamiltonian = Chop[FullHamiltonian];*)
		(**)
		FullHamiltonian = Total[
							MapThread[Times, {Couplings, BasisTerms}]
							];

		FullHamiltonian = Chop[FullHamiltonian];
        (**)

        (* Return the full matrix or the requested parity sector.            *)
        If[Parity === "None",
            Return[FullHamiltonian]
        ];

        FullHamiltonian[[#, #]] &[GetParitySectorIndices[NumSpins, Parity]]
    ];

SYKHamiltonian::evenN = "NumMajoranas must be an even integer.";


End[];

EndPackage[];
