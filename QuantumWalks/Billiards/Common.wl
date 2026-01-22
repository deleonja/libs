(* ::Package:: *)

BeginPackage["QuantumWalks`Billiards`Common`"];

(* S\[IAcute]mbolos exportados gen\[EAcute]ricos *)
BuildGridShiftOperators::usage = "BuildGridShiftOperators[gridData] construye los operadores de desplazamiento (Wm, Wn) para cualquier geometr\[IAcute]a de rejilla 2D definida en gridData.";
GenericSU2Coin::usage = "GenericSU2Coin[\[Alpha], dim] devuelve el operador moneda SU(2) extendido al espacio de Hilbert de dimensi\[OAcute]n dim.";

Begin["`Private`"];

(* --- Construcci\[OAcute]n de Operadores (Gen\[EAcute]rico para cualquier Rejilla) --- *)
BuildGridShiftOperators[gridData_Association] := Module[
    {
        coords, map, dim,
        indicesUp, indicesDown,
        destUp, destDown,
        rulesWm, rulesWn  (* CORREGIDO: Eliminados guiones bajos *)
    },
    
    coords = gridData["Coords"];
    map = gridData["Mapping"];
    dim = gridData["Dimension"];

    (* Estructura del Espacio de Hilbert: \[CapitalIAcute]ndice = 2*(map[{m,n}] - 1) + spin + 1 *)
    indicesUp = 2 * (Range[dim] - 1) + 1;   (* |m,n, U> *)
    indicesDown = indicesUp + 1;            (* |m,n, D> *)

    (* --- Wm (Horizontal) --- *)
    destUp = map /@ (coords + ConstantArray[{1, 0}, dim]);
    destDown = map /@ (coords - ConstantArray[{1, 0}, dim]);
    
    rulesWm = Flatten[{
        (* Spin UP: Si choca (Missing), refleja. Si no, avanza. *)
        MapThread[If[MissingQ[#2], {#1, #1 + 1} -> 1, {#1, 2*(#2 - 1) + 1} -> 1] &, {indicesUp, destUp}],
        (* Spin DOWN *)
        MapThread[If[MissingQ[#2], {#1, #1 - 1} -> 1, {#1, 2*(#2 - 1) + 2} -> 1] &, {indicesDown, destDown}]
    }];

    (* --- Wn (Vertical) --- *)
    destUp = map /@ (coords + ConstantArray[{0, 1}, dim]);
    destDown = map /@ (coords - ConstantArray[{0, 1}, dim]);

    rulesWn = Flatten[{
        (* Spin UP *)
        MapThread[If[MissingQ[#2], {#1, #1 + 1} -> 1, {#1, 2*(#2 - 1) + 1} -> 1] &, {indicesUp, destUp}],
        (* Spin DOWN *)
        MapThread[If[MissingQ[#2], {#1, #1 - 1} -> 1, {#1, 2*(#2 - 1) + 2} -> 1] &, {indicesDown, destDown}]
    }];

    {
        SparseArray[rulesWm, {2 dim, 2 dim}],
        SparseArray[rulesWn, {2 dim, 2 dim}]
    }
]

(* --- Moneda Gen\[EAcute]rica --- *)
GenericSU2Coin[\[Alpha]_, gridDim_Integer] := 
    KroneckerProduct[
        IdentityMatrix[gridDim, SparseArray], 
        SparseArray[{{Cos[\[Alpha]], Sin[\[Alpha]]}, {-Sin[\[Alpha]], Cos[\[Alpha]]}}]
    ]

End[];
EndPackage[];
