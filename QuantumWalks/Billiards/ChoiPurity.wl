(* ::Package:: *)

BeginPackage["QuantumWalks`Billiards`", {"QuantumWalks`"}];

ChoiPurity::usage =
    "ChoiPurity[U, positionState, t] gives the purity of the normalized " <>
    "Choi state of the reduced four-state coin channel after t steps. " <>
    "The Hilbert space is ordered as position then coin, and " <>
    "positionState must be normalized.";

ChoiPurityTrajectory::usage =
    "ChoiPurityTrajectory[U, positionState, times] returns {time, purity} " <>
    "pairs in the order given by times. " <>
    "ChoiPurityTrajectory[U, positionState, tMax] evaluates times " <>
    "0 through tMax.";

Begin["`ChoiPurity`Private`"];

(* Each position contributes one d x d Kraus block to the evolved states. *)
choiPurityFromStates[states_, nPos_Integer, d_Integer : 4] := Module[
    {krausVectors, choi},

    krausVectors = Transpose[ArrayReshape[states, {nPos, d^2}]];
    choi = krausVectors . ConjugateTranspose[krausVectors] / d;

    Chop[Tr[choi . choi]]
];

(* The columns are the evolved states of the four coin basis vectors. *)
initialCoinStates[positionState_, d_Integer : 4] := Transpose[
    Table[
        Flatten[KroneckerProduct[positionState, UnitVector[d, j]]],
        {j, d}
    ]
];

ChoiPurity[
    U_,
    positionState_?VectorQ,
    t_Integer?NonNegative
] := Module[{d = 4, nPos, states},
    nPos = Length[positionState];
    states = initialCoinStates[positionState, d];
    states = Nest[U . # &, states, t];

    choiPurityFromStates[states, nPos, d]
];

ChoiPurityTrajectory[
    U_,
    positionState_?VectorQ,
    times_List
] := Module[
    {d = 4, nPos, states, sortedTimes, currentTime = 0, purities = <||>},

    If[! AllTrue[times, IntegerQ[#] && NonNegative[#] &],
        Return[$Failed]
    ];

    nPos = Length[positionState];
    sortedTimes = Sort[DeleteDuplicates[times]];
    states = initialCoinStates[positionState, d];

    Do[
        states = Nest[U . # &, states, time - currentTime];
        AssociateTo[
            purities,
            time -> choiPurityFromStates[states, nPos, d]
        ];
        currentTime = time,
        {time, sortedTimes}
    ];

    Transpose[{times, Lookup[purities, times]}]
];

ChoiPurityTrajectory[
    U_,
    positionState_?VectorQ,
    tMax_Integer?NonNegative
] := ChoiPurityTrajectory[U, positionState, Range[0, tMax]];

End[];
EndPackage[];
