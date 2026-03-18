(* ::Package:: *)

BeginPackage["QMB`"];


Get["ForScience`"];


MatrixPartialTraceFromPureState::usage = FormatUsage[
"MatrixPartialTraceFromPureState[stateVector,dimA,dimB,sysToTrace] \
	computes the partial trace over ```sysToTrace``` of the density \
	matrix of a pure bipartite state ```stateVector```."
];


(* ::Section:: *)
(*Public definitions*)


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


End[];


EndPackage[];
