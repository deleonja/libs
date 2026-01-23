(* ::Package:: *)

BeginPackage["QuantumWalks`"];

GenerateRectangleBasis::usage = "GenerateRectangleBasis[Lx, Ly] genera la base de coordenadas \
para un billar rectangular (sistema integrable). \
Lx: Ancho (o semiancho si se centra en 0, aqu\[IAcute] definimos longitud total Nx = Lx). \
Ly: Alto.";

Begin["`Billiards`Rectangle`Private`"];

(* --- Billar Rectangular (Integrable) --- *)
(* Definici\[OAcute]n: 0 <= x <= Lx, 0 <= y <= Ly *)
(* Nota: Se puede ajustar para estar centrado en 0 si se prefiere simetr\[IAcute]a *)

GenerateRectangleBasis[Lx_Integer, Ly_Integer] := Module[
    {
        coords, 
        stateToIndex
    },

    (* Generaci\[OAcute]n directa O(1) sin necesidad de Select *)
    coords = Flatten[
        Table[{x, y}, {x, 0, Lx}, {y, 0, Ly}],
        1
    ];

    stateToIndex = AssociationThread[coords -> Range[Length[coords]]];

    <|
        "Coords" -> coords,
        "Mapping" -> stateToIndex,
        "Dimension" -> Length[coords],
        "Params" -> <|"Lx" -> Lx, "Ly" -> Ly|>,
        "Type" -> "Rectangle"
    |>
]

End[];
EndPackage[];
