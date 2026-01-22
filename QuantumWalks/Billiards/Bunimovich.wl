(* ::Package:: *)

(* ::Package:: *)

(* Cargamos Common para tener acceso a BuildGridShiftOperators y GenericSU2Coin *)
BeginPackage["QuantumWalks`Billiards`Bunimovich`", {"QuantumWalks`Billiards`Common`"}];

(* S\[IAcute]mbolos espec\[IAcute]ficos del Estadio *)
GenerateStadiumBasis::usage = "GenerateStadiumBasis[xc, nu] genera la base de posiciones y la asociaci\[OAcute]n de b\[UAcute]squeda para el estadio.";

(* Wrappers *)
BuildShiftOperators::usage = "Alias de BuildGridShiftOperators para el contexto Bunimovich.";
StadiumCoin::usage = "Alias de GenericSU2Coin para el contexto Bunimovich.";

Begin["`Private`"];

(* --- 1. Geometr\[IAcute]a Espec\[IAcute]fica del Estadio --- *)

(* Eqs A7 y A8 del paper *)
BoundaryF[m_, mc_, nu_] := If[m < mc, nu, Round[Sqrt[nu^2 - (m - mc)^2]]];
BoundaryW[n_, mc_, nu_] := mc + Round[Sqrt[nu^2 - n^2]];

GenerateStadiumBasis[xc_Integer, nu_Integer] := Module[
    {coords, stateToIndex},
    
    (* Generaci\[OAcute]n optimizada de coordenadas *)
    coords = Select[
        Flatten[Table[{m, n}, {m, 0, 2 xc + nu}, {n, 0, nu}], 1],
        Function[pos, 
            And[
                0 <= pos[[1]] <= BoundaryW[pos[[2]], xc, nu],
                0 <= pos[[2]] <= BoundaryF[pos[[1]], xc, nu]
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

(* --- 2. Enlaces a Common (Wrappers) --- *)

(* Simplemente pasamos los datos a la funci\[OAcute]n gen\[EAcute]rica optimizada de Common *)
(* Nota: gridData_ es correcto aqu\[IAcute] porque es un patr\[OAcute]n de argumento *)
BuildShiftOperators[gridData_] := BuildGridShiftOperators[gridData];

(* Lo mismo para la moneda *)
StadiumCoin[\[Alpha]_, dim_] := GenericSU2Coin[\[Alpha], dim];

End[];
EndPackage[];
