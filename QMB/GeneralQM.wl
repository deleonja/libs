(* ::Package:: *)

BeginPackage["QMB`"];


Get["ForScience`"];


(* ::Section::Closed:: *)
(*Public definitions*)


Quiet[
MatrixPartialTraceFromPureState::usage = FormatUsage[
"MatrixPartialTraceFromPureState[stateVector,dimA,dimB,sysToTrace] \
	computes the partial trace over ```sysToTrace``` of the density \
	matrix of a pure bipartite state ```stateVector```."
];

Purity::usage = FormatUsage[
"Purity[\[Rho]] calculates the purity of ```\[Rho]```."
];

HaarRandomState::usage = FormatUsage[
"HaarRandomState[dim] returns a Haar-random quantum state \
	vector of dimension ```dim```."
];

PartialTranspose::usage = FormatUsage[
"PartialTranspose[matrix,sysToTranspose,{dimA,dimB}] computes the \
	partial transpose of a bipartite ```matrix``` over the specified \
	```sysToTranspose``` (1 or 2).
PartialTranspose[matrix,sysToTranspose] assumes a bipartite system \
	of equal dimensions (e.g., two qubits)."
];

PartialTranspose::invSub = "Subsystem parameter must be 1 or 2.";
PartialTranspose::invDim = "Matrix dimensions do not match the \
	provided subsystem dimensions.";

\[Omega]::usage = FormatUsage[
"\[Omega][d] returns the first dth root of unity.
\[Omega][d,k] returns \[Omega][d]^k."
];

WeylMatrix::usage = FormatUsage[
"WeylMatrix[m,n,d] returns the Weyl matrix \
	U(```m```,```n```) of dimension ```d```.
WeylMatrix[mList,nList,dList] returns the multiparticle Weyl matrix \
	U(```mList```,```nList```) of dimensions ```dList```."
];

BlochVector::usage = FormatUsage[
"BlochVector[\[Rho]] returns the Bloch vector of a single-qubit density matrix \[Rho].
BlochVector[\[Psi]] returns the Bloch vector of a single-qubit pure state \[Psi]."
];

QubitStateFromBlochVector::usage = FormatUsage[
"QubitStateFromBlochVector[\[Theta],\[Phi]] returns the state cos(```\[Theta]```/2)|0\[RightAngleBracket] + \[ExponentialE]^```\[Phi]``` sin(```\[Theta]```/2)|1\[RightAngleBracket]"
];

SU2RotationParameters::usage = FormatUsage[
"SU2RotationParameters[Matrix] computes the spherical coordinates \[CurlyTheta] \
and \[Phi] of the rotation axis, and the rotation angle \[Theta] from an SU(2) \
matrix ```Matrix```."
];

RandomSU2Parameters::usage = FormatUsage[
"RandomSU2Parameters[] gives a random list {\[CurlyTheta],\[Phi],\[Theta]} of the spherical coordinates \
\[CurlyTheta] and \[Phi] of the rotation axis, and the rotation angle \[Theta] to generate a Haar \
random SU(2) operator in \[DoubleStruckCapitalC]^2 . 
RandomSU2Parameters[n] gives n random lists."
];

RandomSU2Rotation::usage = FormatUsage[
"RandomSU2Rotation[] gives a Haar random SU(2) operator in \[DoubleStruckCapitalC]^2 .
RandomSU2Rotation[n] gives ```n``` Haar random SU(2) operators. "
];

QuantumTexture::usage = FormatUsage[
"QuantumTexture[PkList,evalsList,t] computes the quantum-channel texture grand sum \
for a unitary evolution with Pk=|<Ek|f1>|^2."
];

, {FrontEndObject::notavail, First::normal}];


(* ::Subsection::Closed:: *)
(*Deprecated functions*)


Qubit::usage = "Qubit is deprecated and has been removed. Use \
QubitStateFromBlochVector instead.";

Qubit::deprecated = "Qubit is deprecated and has been removed. Use \
QubitStateFromBlochVector instead.";


(* ::Section::Closed:: *)
(*Private definitions*)


Begin["`GeneralQM`Private`"];


(* ================================================================= *)
(* 1. MOTORES DE C\[CapitalAAcute]LCULO PRIVADOS (Compilados en C)                  *)
(* ================================================================= *)

(* Motor A: Traza el primer subsistema (A), conserva B *)
coreTraceA = Compile[
    {{stateVector, _Complex, 1}, {dimA, _Integer}, {dimB, _Integer}},
    Module[{reducedRho, offset, valI, valJ},
        reducedRho = Table[0.0 + 0.0 I, {dimB}, {dimB}];
        Do[
            offset = (k - 1) * dimB;
            Do[
                valI = stateVector[[offset + i]];
                Do[
                    valJ = stateVector[[offset + j]];
                    reducedRho[[i, j]] += valI * Conjugate[valJ];
                , {j, 1, dimB}]
            , {i, 1, dimB}]
        , {k, 1, dimA}];
        reducedRho
    ],
    CompilationTarget -> "C",
    RuntimeAttributes -> {Listable},
    Parallelization -> True,
    RuntimeOptions -> "Speed"
];

(* Motor B: Traza el segundo subsistema (B), conserva A *)
coreTraceB = Compile[
    {{stateVector, _Complex, 1}, {dimA, _Integer}, {dimB, _Integer}},
    Module[{reducedRho, offsetI, offsetJ, valI, valJ},
        reducedRho = Table[0.0 + 0.0 I, {dimA}, {dimA}];
        Do[
            offsetI = (i - 1) * dimB;
            offsetJ = (j - 1) * dimB;
            Do[
                valI = stateVector[[offsetI + k]];
                valJ = stateVector[[offsetJ + k]];
                reducedRho[[i, j]] += valI * Conjugate[valJ];
            , {k, 1, dimB}]
        , {i, 1, dimA}, {j, 1, dimA}];
        reducedRho
    ],
    CompilationTarget -> "C",
    RuntimeAttributes -> {Listable},
    Parallelization -> True,
    RuntimeOptions -> "Speed"
];

(* ================================================================= *)
(* 2. INTERFAZ P\[CapitalUAcute]BLICA (Dispatcher / Overloading)                    *)
(* ================================================================= *)

(* Caso 1: El usuario quiere trazar el subsistema 1 *)
MatrixPartialTraceFromPureState[state_, dimA_Integer, dimB_Integer, 1] := 
    coreTraceA[state, dimA, dimB];

(* Caso 2: El usuario quiere trazar el subsistema 2 *)
MatrixPartialTraceFromPureState[state_, dimA_Integer, dimB_Integer, 2] := 
    coreTraceB[state, dimA, dimB];

(* Caso de Error: Cualquier otro valor para sysToTrace (Ej. 3, o letras) *)
MatrixPartialTraceFromPureState::invalidSys = "The argument `1` must be either 1 or 2.";
MatrixPartialTraceFromPureState[state_, dimA_Integer, dimB_Integer, sysToTrace_] := 
    (
        Message[MatrixPartialTraceFromPureState::invalidSys,sysToTrace];
        $Failed
    );


Purity[\[Rho]_] := Chop[Tr[\[Rho] . \[Rho]]]


HaarRandomState[Dim_Integer?Positive] := Module[
    {ComplexVector},
    
    (* Sample real and imaginary parts from a standard normal distribution *)
    ComplexVector = RandomVariate[NormalDistribution[], Dim] + 
                    I * RandomVariate[NormalDistribution[], Dim];
    
    (* Normalizing the isotropic complex Gaussian vector yields a Haar-distributed state *)
    Normalize[ComplexVector]
];


(* Main polymorphic definition handling generic bipartite dimensions *)
PartialTranspose[matrix_?MatrixQ, subsystem_Integer, dims_List] := Module[
    {dimA, dimB, tensorForm, permutedTensor},
    
    {dimA, dimB} = dims;

    (* Validation *)
    If[subsystem =!= 1 && subsystem =!= 2,
        Message[PartialTranspose::invSub];
        Return[$Failed]
    ];

    If[Length[matrix] != dimA * dimB || Length[matrix[[1]]] != dimA * dimB,
        Message[PartialTranspose::invDim];
        Return[$Failed]
    ];

    (* Step 1: Reshape into a 4-index tensor {outA, outB, inA, inB} *)
    tensorForm = ArrayReshape[matrix, {dimA, dimB, dimA, dimB}];

    (* Step 2: Swap the indices corresponding to the chosen subsystem *)
    permutedTensor = Transpose[
        tensorForm, 
        If[subsystem == 1, {3, 2, 1, 4}, {1, 4, 3, 2}]
    ];

    (* Step 3: Flatten back to a matrix *)
    ArrayReshape[permutedTensor, {dimA * dimB, dimA * dimB}]
];

(* Overload for symmetric bipartite systems (e.g., two identical spins/qubits) *)
PartialTranspose[matrix_?MatrixQ, subsystem_Integer] := Module[
    {subDim},
    
    (* Infer the dimension assuming dimA == dimB *)
    subDim = Round[Sqrt[Length[matrix]]];
    
    PartialTranspose[matrix, subsystem, {subDim, subDim}]
];


\[Omega][d_] := Exp[ 2*Pi*I / d ]
\[Omega][d_,k_]:=Exp[ 2k*Pi*I / d ]


WeylMatrix[m_Integer,n_Integer,d_Integer] := 
	Sum[SparseArray[Mod[{k,k+n},d]+1->\[Omega][d,k*m],{d,d}],{k,0,d-1}]
	
WeylMatrix[m_List,n_List,d_List] := 
	KroneckerProduct@@(WeylMatrix[#1,#2,#3]&@@@Transpose[Join[{m,n},{d}]])


BlochVector[\[Rho]_?MatrixQ] := Chop[Tr[Pauli[#] . \[Rho]] & /@ Range[3]]

BlochVector[\[Psi]_?VectorQ] := Chop[Braket[\[Psi], Pauli[#] . \[Psi]]] & /@ Range[3]


QubitStateFromBlochVector[\[Theta]_, \[Phi]_] := {Cos[\[Theta]/2], Exp[I*\[Phi]]*Sin[\[Theta]/2]}

QubitStateFromBlochVector[{r_, \[Theta]_, \[Phi]_}] /; Chop[r - 1] == 0 := 
    QubitStateFromBlochVector[\[Theta], \[Phi]]

QubitStateFromBlochVector[{r_, \[Theta]_, \[Phi]_}] /; r < 1 := 
    Module[{CartesianVector, IdentityPart, InteractionPart},
        (* Convert Spherical to Cartesian Bloch Vector *)
        CartesianVector = r * {
            Sin[\[Theta]] * Cos[\[Phi]], 
            Sin[\[Theta]] * Sin[\[Phi]], 
            Cos[\[Theta]]
        };
        
        (* Resulting Density Matrix: 1/2 * (I + b . sigma) *)
        Chop[0.5 * (Pauli[0] + Total[CartesianVector * (Pauli /@ Range[3])])]
    ]


SU2RotationParameters[uMatrix_?MatrixQ] := Module[
    {
        HalfAngleCos, 
        VectorX, VectorY, VectorZ, 
        NormV, 
        ThetaR, VarTheta, Phi
    },
    
    (* Extract the half-angle cosine from the real part of the first element *)
    HalfAngleCos = Re[uMatrix[[1, 1]]];
    ThetaR = 2.0 * ArcCos[HalfAngleCos];
    
    (* Extract the unnormalized Cartesian components of the rotation axis *)
    VectorX = -Im[uMatrix[[1, 2]]];
    VectorY = -Re[uMatrix[[1, 2]]];
    VectorZ = -Im[uMatrix[[1, 1]]];
    
    NormV = Norm[{VectorX, VectorY, VectorZ}];
    
    (* Handle the identity matrix case (zero rotation) *)
    If[Chop[NormV] == 0,
        Return[{0.0, 0.0, 0.0}]
    ];
    
    (* Calculate spherical coordinates from the normalized vector *)
    VarTheta = ArcCos[VectorZ / NormV];
    Phi = Mod[ArcTan[VectorX, VectorY], 2Pi];
    
    {VarTheta, Phi, ThetaR}
]


RandomSU2Parameters[] := Module[
    {PhiAngle, CosTheta, PolarAngle, RotationAngle, UniformSample, ArgVar},
    
    (* Step 1: Sample azimuthal angle uniformly in [0, 2Pi] *)
    PhiAngle = RandomReal[{0, 2 Pi}];
    
    (* Step 2: Sample polar angle uniformly on the sphere *)
    CosTheta = RandomReal[{-1, 1}];
    PolarAngle = ArcCos[CosTheta];
    
    (* Step 3: Sample rotation angle theta with PDF p(theta) = (1-Cos[theta])/(2Pi) *)
    (* We solve the CDF equation: (theta - Sin[theta])/(2Pi) == u *)
    UniformSample = RandomReal[];
    RotationAngle = ArgVar /. FindRoot[
        (ArgVar - Sin[ArgVar]) / (2 Pi) == UniformSample, 
        {ArgVar, Pi}
    ];
    
    (* Returns {PolarAngle, AzimuthalAngle, RotationAngle} *)
    {PolarAngle, PhiAngle, RotationAngle}
]

RandomSU2Parameters[n_] := Table[RandomSU2Parameters[], n]


RandomSU2Rotation[] := Module[
    {PolarAngle, PhiAngle, RotationAngle},
    
    {PolarAngle, PhiAngle, RotationAngle} = RandomSU2Parameters[];
    
    (* Reutilizing the existing SU2Rotation function from QMB *)
    SU2Rotation[{PolarAngle, PhiAngle}, RotationAngle]
]

RandomSU2Rotation[n_] := Table[RandomSU2Rotation[], n]


QuantumTexture[Pk_, evals_, t_] := Module[
	{d = Length[evals]}, 
	d * Abs[Total[Pk * Exp[-I * evals * t]]]^2
]


(* ::Subsection::Closed:: *)
(*Deprecated functions*)


Qubit[args___] := 
  Module[{},
    Message[Qubit::deprecated];
    Return[$Failed]
  ];


End[];


EndPackage[];
