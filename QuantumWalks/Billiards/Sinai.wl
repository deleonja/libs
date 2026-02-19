(* ::Package:: *)

BeginPackage["QuantumWalks`Billiards`", {"QuantumWalks`"}];


<<ForScience`;


GenerateSinaiBasis::usage = FormatUsage[
"GenerateSinaiBasis[L, R] generates the coordinate basis and the mapping \
association for the fundamental domain (1/8) of the Sinai billiard \
with semi-side ```L``` of the square (total side is 2L), and \
radius ```R```." ];

GenerateFullSinaiBasis::usage = FormatUsage[
"GenerateFullSinaiBasis[L, R] generates the coordinate basis for the full \
Sinai billiard with semi-side ```L``` of the square (total side is 2L), and \
radius ```R```." ];


GenerateSinaiBasis::usage = FormatUsage[
"GenerateSinaiBasis[L,R] returns an ///Association/// with the \
grid data for the fundamental domain (1/8) of the Sinai billiard \
with semi-side ```L``` of the square (total side is 2L), and \
radius ```R```. The keys are: \"Coords\", \"Mapping\" , \"Dimension\" \
, \"Params\", \"Type\"."];

GenerateFullSinaiBasis::usage = FormatUsage[
"GenerateFullSinaiBasis[L,R] returns an ///Association/// with the \
grid data for the full Sinai billiard with semi-side ```L``` of \
the square (total side is 2L), and radius ```R```. The keys are: \
\"Coords\", \"Mapping\" , \"Dimension\" , \"Params\", \"Type\"."];


Begin["`Sinai`Private`"];

(* --- 1. Geometr\[IAcute]a del Dominio Fundamental (1/8 Sinai) --- *)

(* El dominio fundamental est\[AAcute] definido por:
   1. Dentro del cuadrado: 0 <= x <= L
   2. Simetr\[IAcute]a diagonal (y <= x): 0 <= y <= x
   3. Exclusi\[OAcute]n del disco: x^2 + y^2 >= R^2
*)

GenerateSinaiBasis[L_Integer, R_Integer] := Module[
    {
        coords, 
        stateToIndex,
        validRegionQ
    },
    
    (* Funci\[OAcute]n compilada para chequeo r\[AAcute]pido de la regi\[OAcute]n *)
    validRegionQ = Function[{x, y}, (x^2 + y^2 >= R^2)];

    (* Generaci\[OAcute]n optimizada de coordenadas *)
    coords = Flatten[
        Table[
            If[validRegionQ[x, y], {x, y}, Nothing],
            {x, 0, L}, 
            {y, 0, x}
        ],
        1
    ];

    (* Creaci\[OAcute]n del Mapeo O(1) *)
    stateToIndex = AssociationThread[coords -> Range[Length[coords]]];

    <|
        "Coords" -> coords,
        "Mapping" -> stateToIndex,
        "Dimension" -> Length[coords],
        "Params" -> <|"L" -> L, "R" -> R, "Symmetry" -> "1/8"|>,
        "Type" -> "Sinai"
    |>
];

(* --- 2. Geometr\[IAcute]a del Sinai Completo --- *)
(* Cuadrado de lado 2L centrado en 0, con un disco de radio R removido del centro *)

GenerateFullSinaiBasis[L_Integer, R_Integer] := Module[
    {
        coords, 
        stateToIndex,
        validRegionQ
    },

    (* L\[OAcute]gica: Dentro del cuadrado Y fuera del c\[IAcute]rculo *)
    validRegionQ = Function[{x, y},
        And[
            Abs[x] <= L,
            Abs[y] <= L,
            x^2 + y^2 >= R^2 (* >= incluye el borde del obst\[AAcute]culo *)
        ]
    ];

    (* Bounding Box: El cuadrado completo [-L, L] *)
    coords = Select[
        Flatten[
            Table[
                {x, y}, 
                {x, -L, L}, 
                {y, -L, L}
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
        "Params" -> <|"L" -> L, "R" -> R, "Symmetry" -> "None"|>,
        "Type" -> "SinaiFull"
    |>
]

End[];
EndPackage[];
