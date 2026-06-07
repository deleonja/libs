(* ::Package:: *)

BeginPackage["QMB`ManyBody`", {"QMB`"}];


Get["ForScience`"];


(* 80-column reference: *)
(*0123456789012345678901234567890123456789012345678901234567890123456789012345*)


(* ::Section:: *)
(*Public definitions*)


(* ::Subsection::Closed:: *)
(*Periodic boundary conditions *)


Quiet[
PeriodicBoseHubbard::usage = FormatUsage[
"PeriodicBoseHubbard[N,L,J,U] returns the Hamiltonian of the 1D Bose Hubbard \
	chain with N sites and L bosons wtih periodic boundary conditions. J \
	and U are the hopping and interactions parameters, respectively."
	]
];

Quiet[
MomentumSectorPeriodicBoseHubbard::usage = FormatUsage[
"MomentumSectorPeriodicBoseHubbard[N,L,J,U,sector:All or interger] returns \
	the Hamiltonian of a 1D Bose Hubbard chain with N sites and L bosons \
	wtih periodic boundary conditions in block diagonal form due to \
	Momentum symmetry. Using the sector=All returns the full sectorized \
	hamiltonian whereas sector= i (i <= N) returns the i-th momentum sector."
	]
];


(* ::Subsection::Closed:: *)
(*Open boundary conditions*)


BoseHubbardHamiltonian::usage = FormatUsage[
"BoseHubbardHamiltonian[n,L,J,U,opts] returns the Bose-Hubbard \
	Hamiltonian with open boundary conditions of ```n``` bosons and \
	```L``` sites with hopping parameter ```J``` and interaction \
	parameter ```U```. Options: SymmetricSubspace."
];


SymmetricSubspace::usage = FormatUsage[
"SymmetricSubspace is an option for '''BoseHubbardHamiltonian'''. \
	Valid values are '''\"All\"''', '''\"EvenParity\"''', and \
	'''\"OddParity\"'''."
];


KineticTermBoseHubbardHamiltonian::usage = FormatUsage[
"KineticTermBoseHubbardHamiltonian[basis] returns the kinetic term \
	of the BH Hamiltonian with ```basis``` a list with Fock basis \
	elements.\n" <>
"KineticTermBoseHubbardHamiltonian[basis,SymmetricSubspace] returns \
	the kinetic term of the BH Hamiltonian in a symmetric subspace \
	with ```basis``` a list with Fock basis elements. Option \
	SymmetricSubspace takes the values \"All\" | \"EvenParity\" | \
	\"OddParity\"."
];


PotentialTermBoseHubbardHamiltonian::usage = FormatUsage[
"PotentialTermBoseHubbardHamiltonian[n, L, SymmetricSubspace] \
	returns the potential term of the Bose-Hubbard Hamiltonian for \
````n``` bosons and ```L``` sites.
PotentialTermBoseHubbardHamiltonian[basis] returns the potential \
	term of the Bose-Hubbard Hamiltonian given the ```basis``` of \
	the Fock space of the bosonic system.\n\n"
(*"**Notes**\n" <>
"- If you want the potential term in a parity symmetry sector, \
```basis``` should be a list containing only the representative \
	Fock states of that symmetric subspace.\n\n" <>
"**Examples of usage**\n" <>
"- '''PotentialTermBoseHubbardHamiltonian[{{3,0,0},{2,1,0},{2,0,1},\
	{1,2,0}}]''' returns the matrix in the odd parity sector of a \
	system with 3 bosons and 3 sites."
*)];


HilbertSpaceDim::replaced = "Function `1` has been replaced by `2`.";


BoseHubbardHilbertSpaceDimension::usage = FormatUsage[
"BoseHubbardHilbertSpaceDimension[n,L] returns the dimension of \
	Hilbert space of a Bose Hubbard system of N bosons and L sites."
];


FockBasis::usage =
"FockBasis[N,L] returns the lexicographical-sorted Fock basis of \
	```N``` bosons and ```L``` sites.";


SortFockBasis::usage = FormatUsage[
"SortFockBasis[fockBasis] returns fockBasis in ascending-order \
	according to the tag of Fock states."
];


Tag::usage = FormatUsage[
"Tag[kList] returns the tag `````\[Sum]_i \[Sqrt](100i+3) ki``` of the Fock state \
	|```kList```\[RightAngleBracket]."
];


FockBasisStateAsColumnVector::usage = FormatUsage[
"FockBasisStateAsColumnVector[fockState,N,L] returns the matrix \
	representation of ```fockState``` for ```N``` bosons and ```L``` \
	sites."
];


FockBasisIndex::usage = FormatUsage[
"FockBasisIndex[fockState, sortedTagsFockBasis] returns the position \
	of ```fockState``` in the tag-sorted Fock basis with tags \
	````sortedTagsFockBasis```."
];


BosonicPartialTrace::usage = FormatUsage[
"BosonicPartialTrace[\[Rho]] calculates the partial trace of \
	```\[Rho]```. Requires initialization via \
	'''InitializationBosonicPartialTrace'''."
];


InitializationBosonicPartialTrace::usage = FormatUsage[
"InitializationBosonicPartialTrace[{```i_1```,...,```i_k```},N,L] \
	initializes variables for '''BosonicPartialTrace''' to calculate \
	the reduced density matrix of sites {```i_1```,...,```i_k```} \
	for ```N``` bosons and ```L``` sites."
];


(* ::Section::Closed:: *)
(*Private definitions*)


Begin["`BoseHubbard`Private`"];


(* ::Subsection::Closed:: *)
(*Periodic boundary conditions*)


(*Secondary*)

listaBase[v_,n_]:=Table[If[i==1,v,0],{i,n}]; 
Repetir[f_,ini_,n_]:=NestList[f,ini,n-1];
Indice[lista_]:=Module[{n=Length[lista]},SelectFirst[Range[n-1],AllTrue[lista[[#+1;;n-1]],(#==0)&]&,Missing["NotFound"]]];
BaseN[lista_,n_]:=Module[{k=Indice[lista],copia=lista},
copia[[k]]-=1;
If[k+1<=Length[copia],copia[[k+1]]=n-Total[copia[[1;;k]]];];
Do[copia[[i]]=0,{i,k+2,Length[copia]}];
copia];
Basis[N_,M_]:=Module[{Base1=listaBase[N,M],n1=(M+N-1)!/(N!*(M-1)!)},Repetir[BaseN[#,N]&,Base1,n1]//Sort];

PolarRep[num_] := Module[
  {z, T},
  z = N[num]; (* fuerza forma num\[EAcute]rica *)
  T = ToPolarCoordinates[{Re[z], Im[z]}];
  Chop[T[[1]] Exp[I T[[2]]]]
];

ShiftOperator[n_,m_]:=Module[{Fock=Basis[n,m],IntFock,S},
IntFock=AssociationThread[Fock,Range[Fock//Length]];
S=RotateRight[#]&/@Fock;
{IntFock[S[[#]]],IntFock[Fock[[#]]]}->1&/@Range[Fock//Length]//SparseArray];

OperatorAij[i_,j_,lista_]:=Module[{copy=lista,valorI,valorII},If[i==j,Return[Null]];
valorI=copy[[i]];
valorII=copy[[j]];
copy[[i]]=valorI+1;
copy[[j]]=valorII-1;
If[Min[copy]<0,Return[Null]];
{copy,Sqrt[(valorI+1.)*valorII]}];

UnitaryAtoB[EigenvecA_,EigenvecB_]:=Module[{eigenvecB=EigenvecB,eigenvecA=EigenvecA,FockToIntH,FockToIntS},FockToIntH=AssociationThread[Range[eigenvecA//Length],eigenvecA];FockToIntS=AssociationThread[Range[eigenvecB//Length],eigenvecB];
Table[{#,i}->Chop[Dot[FockToIntH[#],FockToIntS[i]]]&/@Range[FockToIntH//Length],{i,FockToIntS//Length}]//Flatten//SparseArray];

DegenerateSpectrum[S_]:=Module[{Polars,Vecs,Eigens,EigenInfoS1,InfoS,D1},
EigenInfoS1=S//Eigensystem;
Eigens=EigenInfoS1[[1]];
Vecs=EigenInfoS1[[2]];
Polars=PolarRep[#]&/@Eigens;
InfoS={Polars[[#]],Vecs[[#]]}&/@Range[Length[Vecs]];
D1=GatherBy[InfoS,First];
Table[Join[{D1[[i,1,1]]//Chop},{D1[[i,#,2]]&/@Range[D1[[i]]//Length]}],{i,Length[D1]}]];

(*Main*)
MomentumSectorPeriodicBoseHubbard[n_,m_,J_,U_,sector_:All]:=Module[{S,base,DataS,vecs,U2,H},H=PeriodicBoseHubbard[n,m,J,U];
S=ShiftOperator[n,m];
base=IdentityMatrix[Length[S]];
DataS=DegenerateSpectrum[S];
vecs=Normalize/@Which[sector===All,Flatten[Last/@DataS,1],IntegerQ[sector]&&1<=sector<=Length[DataS],DataS[[sector,2]],True,(Message[SectorPeriodicBH::badsector,sector];Return[$Failed])];
U2=UnitaryAtoB[base,vecs];
Chop[ConjugateTranspose[U2] . H . U2]];

PeriodicBoseHubbard[n_,m_,J2_,U_]:=Module[{D1,Fock=Basis[n,m]//Sort,IntFock,Pares,Data,J1,len,diag,H2,Kin},
len=Fock//Length;
IntFock=AssociationThread[Fock,Range[Fock//Length]];
Pares=Table[{Mod[i+1,Fock[[1]]//Length,1],Mod[i,Fock[[1]]//Length,1]},{i,1,Fock[[1]]//Length}];
Data=Table[D1=DeleteCases[OperatorAij[#[[1]],#[[2]],Fock[[i]]]&/@Pares,Null];{IntFock[D1[[#,1]]],IntFock[Fock[[i]]]}->D1[[#,2]]&/@Range[Length[DeleteCases[OperatorAij[#[[1]],#[[2]],Fock[[i]]]&/@Pares,Null]]],{i,Fock//Length}]//Flatten;
J1=Data//SparseArray;
diag=(U/2)*(Total[#*(#-1)]&/@Fock);
H2=SparseArray[Band[{1,1}]->diag,{len,len}];
Kin=-J2(J1+ConjugateTranspose[J1]);
Kin+H2];


(* ::Subsection::Closed:: *)
(*Open boundary conditions*)


(* Mensajes para tipos incorrectos *)
BoseHubbardHamiltonian::int = 
    "n (`1`) and L (`2`) are expected to be integers.";
BoseHubbardHamiltonian::real = 
    "Hopping J (`1`) and interaction parameters U (`2`) are expected " <>
    "to be real numbers.";
BoseHubbardHamiltonian::badSymmetricSubspace = 
    "Opci\[OAcute]n SymmetricSubspace `1` inv\[AAcute]lida. " <>
    "Opciones v\[AAcute]lidas: \"All\", \"EvenParity\" o \"OddParity\".";

Options[BoseHubbardHamiltonian] = {
    SymmetricSubspace -> "All" (* "All"|"EvenParity"|"OddParity" *),
    Version -> "New"
};

OptionValuePatterns[BoseHubbardHamiltonian] = {
  SymmetricSubspace -> Alternatives["All", "EvenParity", "OddParity"],
  Version -> _
};

BoseHubbardHamiltonian[n_Integer, L_Integer, J_Real, U_Real, 
    OptionsPattern[]] := 
Module[
    {tags, basis, basiseven, rbasiseven, rbasisodd, basisodd, H, T, V, map},
    
    basis = N[FockBasis[n, L]];
    H = -J*KineticTermBoseHubbardHamiltonian[basis] + 
        U/2*PotentialTermBoseHubbardHamiltonian[basis];
    
    (* Para subespacios de simetria, obtenerlos a partir de H *)
    Switch[OptionValue[SymmetricSubspace],
        "All",
            Nothing,
            
        "EvenParity",
            basiseven = DeleteDuplicatesBy[basis, Sort[{#, Reverse[#]}]&];
            rbasiseven = Reverse /@ basiseven;
            map = AssociationThread[
                basis -> Range[BoseHubbardHilbertSpaceDimension[n, L]]];
            H = 1/2 # . (H[[#1,#1]] + H[[#1,#2]] + H[[#2,#1]] + H[[#2,#2]] & @@ 
                Map[map, {basiseven, rbasiseven}, {2}]) . # &[
                DiagonalMatrix[
                    ReplacePart[
                        ConstantArray[1., Length[basiseven]],
                        Thread[
                            Flatten[Position[
                                basiseven, 
                                _?(PalindromeQ[#] &), 
                                {1}
                            ]] -> 1/Sqrt[2.]
                        ]
                    ],
                    TargetStructure -> "Sparse"
                ]
            ],
            
        "OddParity",
            basisodd = DeleteDuplicatesBy[
                Discard[basis, PalindromeQ], 
                Sort[{#, Reverse[#]}]&];
            rbasisodd = Reverse /@ basisodd;
            map = AssociationThread[
                basis -> Range[BoseHubbardHilbertSpaceDimension[n, L]]];
            H = 1/2 (H[[#1,#1]] - H[[#1,#2]] - H[[#2,#1]] + H[[#2,#2]] & @@ 
                Map[map, {basisodd, rbasisodd}, {2}]),
                
        _,
            Message[
                BoseHubbardHamiltonian::badSymmetricSubspace, 
                OptionValue[SymmetricSubspace]
            ];
            Return[$Failed];
    ];
    
    H
]

(* Handle cases where arguments don't match the expected types *)
BoseHubbardHamiltonian[N_, L_, J_, U_, OptionsPattern[]] := Module[{},
  If[!IntegerQ[N] || !IntegerQ[L],
    Message[BoseHubbardHamiltonian::int, N, L];
    Return[$Failed];
  ];
  If[Head[J] =!= Real || Head[U] =!= Real,
    Message[BoseHubbardHamiltonian::real, J, U];
    Return[$Failed];
  ];
];

SyntaxInformation[BoseHubbardHamiltonian] = <|
  "ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}
|>;


ClearAll[PotentialTermBoseHubbardHamiltonian];

(* Mensajes de error *)
PotentialTermBoseHubbardHamiltonian::nint = 
    "El n\[UAcute]mero de part\[IAcute]culas n (`1`) debe ser un entero positivo.";
PotentialTermBoseHubbardHamiltonian::lint = 
    "El n\[UAcute]mero de sitios L (`1`) debe ser un entero positivo.";
PotentialTermBoseHubbardHamiltonian::int = 
    "Los par\[AAcute]metros n (`1`) y L (`2`) deben ser enteros positivos.";
PotentialTermBoseHubbardHamiltonian::empty = 
    "La base proporcionada est\[AAcute] vac\[IAcute]a.";
PotentialTermBoseHubbardHamiltonian::dim = 
    "Todos los estados en la base deben tener la misma longitud.";
    
PotentialTermBoseHubbardHamiltonian[n_Integer?Positive, L_Integer?Positive] := 
DiagonalMatrix[
    Total /@ ((#^2 - #) &[FockBasis[n, L]]),
    TargetStructure -> "Sparse"
]

(* Versi\[OAcute]n con validaci\[OAcute]n de basis *)
PotentialTermBoseHubbardHamiltonian[basis_List] := 
Module[
    {lengths},
    If[basis === {},
        Message[PotentialTermBoseHubbardHamiltonian::empty];
        Return[$Failed]
    ];
    
    lengths = Length /@ basis;
    If[!AllTrue[lengths, EqualTo[First[lengths]]],
        Message[PotentialTermBoseHubbardHamiltonian::dim];
        Return[$Failed]
    ];
    
    DiagonalMatrix[
        Total /@ ((#^2 - #) &[basis]),
        TargetStructure -> "Sparse"
    ]
]

(* Versi\[OAcute]n con validaci\[OAcute]n de par\[AAcute]metros *)
PotentialTermBoseHubbardHamiltonian[n_Integer, L_Integer] := 
If[n <= 0 || L <= 0,
    Message[PotentialTermBoseHubbardHamiltonian::int, n, L];
    Return[$Failed]
];

PotentialTermBoseHubbardHamiltonian[n_, L_] /; !IntegerQ[n] := 
(
    Message[PotentialTermBoseHubbardHamiltonian::nint, n];
    $Failed
)

PotentialTermBoseHubbardHamiltonian[n_, L_] /; !IntegerQ[L] := 
(
    Message[PotentialTermBoseHubbardHamiltonian::lint, L];
    $Failed
)

PotentialTermBoseHubbardHamiltonian[___] := 
(
    Message[PotentialTermBoseHubbardHamiltonian::usage];
    $Failed
)


KineticTermBoseHubbardHamiltonian::badSymmetricSubspace = 
    "Opci\[OAcute]n SymmetricSubspace `1` inv\[AAcute]lida. " <>
    "Opciones v\[AAcute]lidas: \"All\", \"EvenParity\" o \"OddParity\".";

Options[KineticTermBoseHubbardHamiltonian] = {
    SymmetricSubspace -> "All" (* "All"|"EvenParity"|"OddParity" *)
};

KineticTermBoseHubbardHamiltonian[basis_, OptionsPattern[]] := 
Module[
    {len = Length[basis], basisNumRange, T, basiseven, rbasiseven, 
     basisodd, rbasisodd, map},
    
    basisNumRange = Range[len];
    T = # + ConjugateTranspose[#] & [
        SparseArray[
            AssignationRulesKinetic[basis, 
             AssociationThread[basis -> basisNumRange], basisNumRange], 
            {len, len}
        ]
    ];
    
    Switch[OptionValue[SymmetricSubspace],
        "All",
            Nothing,
        "EvenParity",
            basiseven = DeleteDuplicatesBy[basis, Sort[{#, Reverse[#]}]&];
            rbasiseven = Reverse /@ basiseven;
            map = AssociationThread[basis -> Range[len]];
            T = 1/2 # . (T[[#1,#1]] + T[[#1,#2]] + T[[#2,#1]] + T[[#2,#2]] & @@ 
                Map[map, {basiseven, rbasiseven}, {2}]) . # &[
                DiagonalMatrix[
                    ReplacePart[
                        ConstantArray[1., Length[basiseven]],
                        Thread[
                            Flatten[Position[
                                basiseven, 
                                _?(PalindromeQ[#] &), 
                                {1}
                            ]] -> 1/Sqrt[2.]
                        ]
                    ],
                    TargetStructure -> "Sparse"
                ]
            ],
        "OddParity",
            basisodd = DeleteDuplicatesBy[
                Discard[basis, PalindromeQ], Sort[{#, Reverse[#]}]&];
            rbasisodd = Reverse /@ basisodd;
            map = AssociationThread[basis -> Range[len]];
            T = 1/2 (T[[#1,#1]] - T[[#1,#2]] - T[[#2,#1]] + T[[#2,#2]] & @@ 
                Map[map, {basisodd, rbasisodd}, {2}]),
            
        _,
            Message[
                KineticTermBoseHubbardHamiltonian::badSymmetricSubspace,
                OptionValue[SymmetricSubspace]
            ];
            Return[$Failed];
    ];
    
    T
]


AssignationRulesKinetic[basis_, positionMap_, basisNumRange_] := 
Catenate[
    MapThread[
        AssignationRulesKineticMapFunc,
        {
            Apply[
                {positionMap[#1], #2} &,
                DeleteCases[
                    Transpose[
                        {
                            StateAfterADaggerA[basis],
                            ValueAfterADaggerA[basis]
                        },
                        {3, 1, 2}
                    ],
                    {_, 0.},
                    {2}
                ],
                {2}
            ],
            basisNumRange
        }
    ]
]


StateAfterADaggerA[basis_] := 
Module[{len = Length[First[basis]]},
    Outer[
        Plus,
        basis,
        #,
        1
    ] & [
        Catenate[
            NestList[
                RotateRight,
                PadRight[#, len],
                len - 2
            ] & /@ {{1, -1}}
        ]
    ]
]


FockBasisStateAsColumnVector[FockBasisState_, N_, L_] :=
    Normal[SparseArray[
        Position[
            SortFockBasis[Normal[FockBasis[N, L]]][[2]], FockBasisState
        ] -> 1,
        Binomial[N + L - 1, L]
    ]]

(* ------------------------------------------------------
   FockBasis[N,M] returns the lexicographical-sorted Fock basis for
   N bosons in M sites.
   New implementation as of 15/Jun/2025 uses more efficient algorithm
   based on IntegerPartitions.
   ------------------------------------------------------ *)
FockBasis[N_, M_] :=
    ReverseSort[
        Catenate[Permutations[PadRight[#, M]] & /@ IntegerPartitions[N, M]]
    ]

(* Old implementation preserved for reference *)
(*
FockBasis[N_, M_] := Module[{k, fockState},
    k = 1;
    Normal[Join[
        {fockState = SparseArray[{1 -> N}, {M}]},
        Table[
            fockState = SparseArray[
                Join[Table[i -> fockState[[i]], {i, k - 1}],
                    {k -> fockState[[k]] - 1}],
                {M}
            ];
            fockState[[k + 1]] = N - Total[fockState[[1 ;; k]]];
            k = Assignationk[M, N, fockState];
            fockState,
        HilbertSpaceDim[N, M] - 1]
    ]]
]
*)

SortFockBasis[fockBasis_] :=
    Transpose[Sort[{Tag[#], #} & /@ fockBasis]]

ValueAfterADaggerA[basis_] :=
    MapApply[
        Sqrt[(#1 + 1.) * #2] &,
        Partition[#, 2, 1]
    ] & /@ basis

AssignationRulesKineticMapFunc[stateValuePairs_, index_] :=
    ({index, #1} -> #2) & @@@ stateValuePairs

Tag[FockBasisElement_] :=
    N[Round[
        Sum[
            Sqrt[100 i + 3] #[[i + 1]], 
            {i, 0, Length[#] - 1}
        ] & [FockBasisElement],
        10^-8
    ]]

Tag[Nothing] := Nothing

InitializationBosonicPartialTrace[SitesOfSubsystem_, n_, L_] :=
    Module[{
        SubsystemSize,
        SubsystemComplementSize,
        SystemFockBasisTags,
        SystemFockBasis,
        SubsystemFockBasisTags,
        SubsystemFockBasis,
        SubsystemComplementFockBasis,
        RulesForOrderedSystemBasis,
        RulesForOrderedSubsystemBasis,
        FockIndicesInRho,
        FockIndicesInReducedRho
        },

        SubsystemSize = Length[SitesOfSubsystem];
        SubsystemComplementSize = L - SubsystemSize;

        (* System's Fock basis *)
        {SystemFockBasisTags, SystemFockBasis} =
            SortFockBasis[FockBasis[n, L]];

        (* Subsystem's Fock basis *)
        {SubsystemFockBasisTags, SubsystemFockBasis} =
            SortFockBasis[
                Flatten[Table[FockBasis[k, SubsystemSize], {k, 0, n}], 1]
            ];

        (* Complement subsystem's Fock basis *)
        SubsystemComplementFockBasis =
            Map[
                ReplacePart[
                    ConstantArray[_, L],
                    Thread[Complement[Range[L], SitesOfSubsystem] -> #]
                ] &,
                Flatten[Table[
                    SortFockBasis[FockBasis[k, SubsystemComplementSize]][[2]],
                    {k, 0, n}
                ], 1]
            ]; (*<<<*)

        RulesForOrderedSystemBasis =
            Thread[Rule[
                SystemFockBasisTags,
                Range[BoseHubbardHilbertSpaceDimension[n, L]]
            ]];

        SubsystemHilbertSpaceDim = Length[SubsystemFockBasis];

        RulesForOrderedSubsystemBasis =
            Thread[
                Rule[SubsystemFockBasisTags, 
                Range[SubsystemHilbertSpaceDim]]
            ];

        FockIndicesInRho =
            Map[
                Tuples[{#, #}] & [
                    Extract[SystemFockBasis, 
                    Position[SystemFockBasis, #]]
                ] &, SubsystemComplementFockBasis
            ]; (*<<<*)

        FockIndicesInReducedRho =
            Map[#[[All, All, SitesOfSubsystem]] &, FockIndicesInRho]; (*<<<*)

        ComputationalIndicesInRho =
            ReplaceAll[
                Map[Tag, FockIndicesInRho, {3}],
                RulesForOrderedSystemBasis
            ]; (*<<<*)

        ComputationalIndicesInReducedRho =
            ReplaceAll[
                Map[Tag, FockIndicesInReducedRho, {3}],
                RulesForOrderedSubsystemBasis
            ]; (*<<<*)
    ]

BosonicPartialTrace[Rho_] :=
    Module[{MatrixElementsOfRho, rules},
        MatrixElementsOfRho = Extract[Rho, #] & /@ ComputationalIndicesInRho;
        rules = MapThread[
            Thread[Rule[#1, #2]] &,
            {ComputationalIndicesInReducedRho, MatrixElementsOfRho}
        ];
        Total[Map[
            SparseArray[
                #, 
                {SubsystemHilbertSpaceDim, SubsystemHilbertSpaceDim}
            ] &, rules
        ]]
    ]


End[];


EndPackage[];
