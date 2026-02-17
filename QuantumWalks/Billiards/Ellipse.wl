(* ::Package:: *)

BeginPackage["QuantumWalks`Billiards`", {"QuantumWalks`"}];

GenerateEllipseBasis::usage = "GenerateEllipseBasis[A, B] genera la base de coordenadas para un billar \
el\[IAcute]ptico. 
A: Semieje mayor (Horizontal). 
B: Semieje menor (Vertical). 
Si A == B, genera un billar circular.";

Begin["`Ellipse`Private`"];

(* --- Billar El\[IAcute]ptico --- *)
(* Ecuaci\[OAcute]n: (x/A)^2 + (y/B)^2 <= 1 *)
(* L\[OAcute]gica de rejilla: Para cada x, |y| <= B * Sqrt[1 - (x/A)^2] *)

GenerateEllipseBasis[A_Integer, B_Integer] := Module[
    {
        coords, 
        stateToIndex,
        validRegionQ
    },

    (* Predicado optimizado con Round para evitar puntas aisladas en los v\[EAcute]rtices *)
    validRegionQ = Function[{x, y},
        Abs[y] <= Round[B * Sqrt[Max[0, 1.0 - (x/A)^2]]]
    ];

    (* Bounding Box: [-A, A] x [-B, B] *)
    coords = Select[
        Flatten[
            Table[
                {x, y}, 
                {x, -A, A}, 
                {y, -B, B}
            ], 
            1
        ],
        Apply[validRegionQ]
    ];

    stateToIndex = AssociationThread[coords -> Range[Length[coords]]];

    <|
        "Coords" -> coords,
        "Mapping" -> stateToIndex,
        "Dimension" -> Length[coords],
        "Params" -> <|"A" -> A, "B" -> B|>,
        "Type" -> If[A == B, "Circle", "Ellipse"]
    |>
]

End[];
EndPackage[];
