(* ::Package:: *)

BeginPackage["QuantumWalks`Billiards`", {"QuantumWalks`"}];

GenerateBunimovichBasis::usage = "GenerateBunimovichBasis[L, R] generates de position basis\
and association map for a stadium of lenght L and radius R.";

GenerateFullStadiumBasis::usage = "GenerateFullStadiumBasis[L, R] genera la base de coordenadas \
y la asociaci\[OAcute]n de mapeo para el billar de estadio completo (Full Stadium). \
L: Semilongitud de la parte recta. \
R: Radio de los semic\[IAcute]rculos.";

Begin["`Bunimovich`Private`"];

(* --- 1. Geometr\[IAcute]a del Estadio de Bunimovich desimetrizado (1/4) --- *)

(* Eqs A7 y A8 del paper de Alonso-Lobo (2025) *)
BoundaryY[m_, mc_, nu_] := If[m <= mc, nu, Round[Sqrt[nu^2 - (m - mc)^2]]];
BoundaryX[n_, mc_, nu_] := mc + Round[Sqrt[nu^2 - n^2]];

GenerateBunimovichBasis[xc_Integer, nu_Integer] := Module[
    {coords, stateToIndex},
    
    (* Generaci\[OAcute]n optimizada de coordenadas *)
    coords = Select[
        Flatten[Table[{m, n}, {m, 0, 2 xc + nu}, {n, 0, nu}], 1],
        Function[pos, 
            And[
                0 <= pos[[1]] <= BoundaryX[pos[[2]], xc, nu],
                0 <= pos[[2]] <= BoundaryY[pos[[1]], xc, nu]
            ]
        ]
    ];

    (* Mapeo O(1) *)
    stateToIndex = AssociationThread[coords -> Range[Length[coords]]];

    <|
        "Coords" -> coords,
        "Mapping" -> stateToIndex,
        "Dimension" -> Length[coords],
        "Params" -> {xc, nu},
        "Type" -> "Bunimovich"
    |>
]

(* --- 2. Geometr\[IAcute]a del Estadio Completo (Full Stadium) --- *)

GenerateFullStadiumBasis[L_Integer, R_Integer] := Module[
    {
        coords, 
        stateToIndex,
        validRegionQ
    },

    (* L\[OAcute]gica basada en Round para consistencia visual y topol\[OAcute]gica *)
    (* Esto asegura que el ancho del estadio en 'y' coincida con la definici\[OAcute]n de Bunimovich.wl original *)
    validRegionQ = Function[{x, y},
        Module[{dx = Abs[x], dy = Abs[y]},
            If[dx <= L,
                dy <= R, (* Regi\[OAcute]n Rectangular *)
                (* Regi\[OAcute]n Semicircular: Usamos Round para determinar la extensi\[OAcute]n en X dada Y *)
                (* Esto evita que el punto (L+R, 0) quede aislado visualmente *)
                dx <= L + Round[Sqrt[Max[0, R^2 - dy^2]]] 
            ]
        ]
    ];

    (* Generaci\[OAcute]n de coordenadas: Iteramos sobre el Bounding Box y seleccionamos *)
    (* El rango de Y es estricto [-R, R] para evitar Ra\[IAcute]ces negativas *)
    coords = Select[
        Flatten[
            Table[
                {x, y}, 
                {x, -L - R, L + R}, 
                {y, -R, R}
            ], 
            1
        ],
        Apply[validRegionQ]
    ];

    (* Mapeo O(1) *)
    stateToIndex = AssociationThread[coords -> Range[Length[coords]]];

    <|
        "Coords" -> coords,
        "Mapping" -> stateToIndex,
        "Dimension" -> Length[coords],
        "Params" -> <|"L" -> L, "R" -> R, "Symmetry" -> "None"|>,
        "Type" -> "BunimovichFull"
    |>
]

End[];
EndPackage[];
