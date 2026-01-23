(* ::Package:: *)

BeginPackage["QuantumWalks`"];

BuildShiftOperators::usage = "BuildShiftOperators[gridData] construye los operadores de \
desplazamiento (Wn, Wm) usando la l\[OAcute]gica com\[UAcute]n de rejilla.";

Begin["`Billiards`Common`Private`"];

(* --- Construcci\[OAcute]n de Operadores (Gen\[EAcute]rico para cualquier geometr\[IAcute]a) --- *)
BuildShiftOperators[gridData_Association] := Module[
    {
        coords, map, dim,
        indicesUp, indicesDown,
        destUp, destDown,
        rulesSx, rulesSy  
    },
    
    coords = gridData["Coords"];
    map = gridData["Mapping"];
    dim = gridData["Dimension"];

    (* Estructura del Espacio de Hilbert: \[CapitalIAcute]ndice = 2*(map[{x, y}] - 1) + spin + 1 *)
    indicesUp = 2 * (Range[dim] - 1) + 1;   (* |x,y,0> *)
    indicesDown = indicesUp + 1;            (* |x,y,1> *)

    (* --- Sx (Horizontal conditional shift) --- *)
    destUp = map /@ (coords + ConstantArray[{1, 0}, dim]);
    destDown = map /@ (coords - ConstantArray[{1, 0}, dim]);

    rulesSx = Flatten[{
        (* Spin UP: Si se sali\[OAcute] del grid (Missing), refleja. Si no, avanza. *)
        MapThread[If[MissingQ[#2], {#1 + 1, #1} -> 1, {2*(#2 - 1) + 1, #1} -> 1] &, {indicesUp, destUp}],
        (* Spin DOWN *)
        MapThread[If[MissingQ[#2], {#1 - 1, #1} -> 1, {2*(#2 - 1) + 2, #1} -> 1] &, {indicesDown, destDown}]
    }];


    (* --- Sy (Vertical) --- *)
    destUp = map /@ (coords + ConstantArray[{0, 1}, dim]);
    destDown = map /@ (coords - ConstantArray[{0, 1}, dim]);

    rulesSy = Flatten[{
        (* Spin UP *)
        MapThread[If[MissingQ[#2], {#1 + 1, #1} -> 1, {2*(#2 - 1) + 1, #1} -> 1] &, {indicesUp, destUp}],
        (* Spin DOWN *)
        MapThread[If[MissingQ[#2], {#1 - 1, #1} -> 1, {2*(#2 - 1) + 2, #1} -> 1] &, {indicesDown, destDown}]
    }];

    {
        SparseArray[rulesSx, {2 dim, 2 dim}],
        SparseArray[rulesSy, {2 dim, 2 dim}]
    }

]

End[];
EndPackage[];
