(* ::Package:: *)

BeginPackage["QMB`"];


Get["ForScience`"];


(* ::Section:: *)
(*Public definitions*)


MatrixPartialTraceFromPureState::usage = FormatUsage[
"MatrixPartialTraceFromPureState[stateVector,dimA,dimB,sysToTrace] \
	computes the partial trace over ```sysToTrace``` of the density \
	matrix of a pure bipartite state ```stateVector```."
];


Purity::usage = FormatUsage[
"Purity[\[Rho]] calculates the purity of ```\[Rho]```."
];


HaarRandomState::usage = FormatUsage[
"RandomHaarState[dim] returns a Haar-random quantum state \
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


(* ::Section:: *)
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


End[];


EndPackage[];
