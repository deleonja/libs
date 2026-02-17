(* ::Package:: *)

BeginPackage["QuantumWalks`Billiards`", {"QuantumWalks`"}];

GenerateCardioidBasis::usage = "GenerateCardioidBasis[R] genera la base de coordenadas \
para un billar con forma de cardioide definido por r = 2R(1 - cos(\[Theta])). \
R: Par\[AAcute]metro de escala (Radio base).";

Begin["`Cardioid`Private`"];

(* --- Geometr\[IAcute]a del Cardioide --- *)
(* Definici\[OAcute]n polar: r = 2R(1 - Cos[\[Theta]]) *)
(* Orientaci\[OAcute]n: C\[UAcute]spide en el origen (0,0), cuerpo principal hacia x negativos. *)
(* Nota f\[IAcute]sica: El cardioide se extiende en X desde -4R hasta +0.5R. *)
(* Ecuaci\[OAcute]n Algebraica: (x^2 + y^2 + 2Rx)^2 <= 4R^2(x^2 + y^2) *)

GenerateCardioidBasis[R_Integer] := Module[
    {
        coords, 
        stateToIndex,
        validRegionAlgebraicQ,
        boundingBoxX, boundingBoxY
    },

    (* Predicado Algebraico *)
    validRegionAlgebraicQ = Function[{x, y},
        Module[{r2 = x^2 + y^2},
            (* La desigualdad define el interior del cardioide *)
            (r2 + 2 * R * x)^2 <= 4 * R^2 * r2
        ]
    ];

    (* Bounding Box Corregido *)
    (* Error previo: Cortaba en 0. El cardioide tiene l\[OAcute]bulos que llegan hasta x = R/2 *)
    (* X: [-4R, R] cubre todo el rango necesario *)
    (* Y: [-3R, 3R] cubre el ancho vertical (~2.6R) *)
    boundingBoxX = {-4 * R, R}; 
    boundingBoxY = {-3 * R, 3 * R};

    coords = Select[
        Flatten[
            Table[
                {x, y}, 
                {x, boundingBoxX[[1]], boundingBoxX[[2]]}, 
                {y, boundingBoxY[[1]], boundingBoxY[[2]]}
            ], 
            1
        ],
        Apply[validRegionAlgebraicQ]
    ];

    stateToIndex = AssociationThread[coords -> Range[Length[coords]]];

    <|
        "Coords" -> coords,
        "Mapping" -> stateToIndex,
        "Dimension" -> Length[coords],
        "Params" -> <|"R" -> R|>,
        "Type" -> "Cardioid"
    |>
]

End[];
EndPackage[];
