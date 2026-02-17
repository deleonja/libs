(* ::Package:: *)

BeginPackage["QuantumWalks`Billiards`", {"QuantumWalks`"}];

BuildShiftOperators::usage = "BuildShiftOperators[gridData, opts] construye \
los operadores de desplazamiento para la rejilla dada.

Opciones:
  CoinDimension -> 2 (Por defecto): Retorna lista {Sx, Sy} para caminata dividida.
  CoinDimension -> 4: Retorna matriz \[UAcute]nica S para caminata completa con 4 estados \
internos (Arriba, Abajo, Derecha, Izquierda).";

CoinDimension::usage = "Opci\[OAcute]n para especificar la dimensi\[OAcute]n del espacio de moneda (2 o 4).";

Options[BuildShiftOperators] = {
    CoinDimension -> 2
};

GetSymmetrySectorBasis::usage = "GetSymmetrySectorBasis[gridData, sector] construye \
la matriz de isometr\[IAcute]a V.";

Begin["`Common`Private`"];

(*Dispatcher*)
BuildShiftOperators[gridData_Association, opts:OptionsPattern[]] := 
    Module[{dimCoin},
        dimCoin = OptionValue[CoinDimension];
        
        Switch[dimCoin,
            2, BuildShiftOperators2State[gridData],
            4, BuildShiftOperators4State[gridData],
            _, Message[BuildShiftOperators::invOpt, dimCoin]; $Failed
        ]
    ];

BuildShiftOperators::invOpt = "Valor de CoinDimension `1` no soportado. Use 2 o 4.";

(* --- Construcci\[OAcute]n de Operadores (Gen\[EAcute]rico para cualquier geometr\[IAcute]a) --- *)
BuildShiftOperators2State[gridData_Association] := Module[
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

(* --- Nueva L\[OAcute]gica (4 Estados - Grover Walk style) --- *)
(* Retorna un \[UAcute]nico operador S global *)
(* Codificaci\[OAcute]n: 1:Arriba(+y), 2:Abajo(-y), 3:Derecha(+x), 4:Izquierda(-x) *)
BuildShiftOperators4State[gridData_] := Module[
    {
        coords, map, dim,
        idxUp, idxDown, idxRight, idxLeft,
        dstUp, dstDown, dstRight, dstLeft,
        rules
    },
    
    coords = gridData["Coords"];
    map = gridData["Mapping"];
    dim = gridData["Dimension"];
    
    (* \[CapitalIAcute]ndices base para cada estado *)
    idxUp    = 4 * (Range[dim] - 1) + 1;
    idxDown  = idxUp + 1;
    idxRight = idxUp + 2;
    idxLeft  = idxUp + 3;

    (* Destinos vectorizados *)
    dstUp    = map /@ (coords + ConstantArray[{0, 1}, dim]); (* y+1 *)
    dstDown  = map /@ (coords - ConstantArray[{0, 1}, dim]); (* y-1 *)
    dstRight = map /@ (coords + ConstantArray[{1, 0}, dim]); (* x+1 *)
    dstLeft  = map /@ (coords - ConstantArray[{1, 0}, dim]); (* x-1 *)

    (* Construcci\[OAcute]n de reglas con rebote unitario *)
    (* Si choca (MissingQ), invierte direcci\[OAcute]n: Up<->Down, Right<->Left *)
    rules = Flatten[{
        (* Up (1) -> Si choca, va a Down (2) en mismo sitio. Si no, va a Up (1) en destino *)
        MapThread[If[MissingQ[#2], {#1 + 1, #1} -> 1, {4*(#2 - 1) + 1, #1} -> 1] &, {idxUp, dstUp}],
        
        (* Down (2) -> Si choca, va a Up (1). Si no, va a Down (2) *)
        MapThread[If[MissingQ[#2], {#1 - 1, #1} -> 1, {4*(#2 - 1) + 2, #1} -> 1] &, {idxDown, dstDown}],
        
        (* Right (3) -> Si choca, va a Left (4). Si no, va a Right (3) *)
        MapThread[If[MissingQ[#2], {#1 + 1, #1} -> 1, {4*(#2 - 1) + 3, #1} -> 1] &, {idxRight, dstRight}],
        
        (* Left (4) -> Si choca, va a Right (3). Si no, va a Left (4) *)
        MapThread[If[MissingQ[#2], {#1 - 1, #1} -> 1, {4*(#2 - 1) + 4, #1} -> 1] &, {idxLeft, dstLeft}]
    }];

    (* Retorna una \[UAcute]nica matriz dispersa S *)
    SparseArray[rules, {4 dim, 4 dim}]
];

(* Helpers de Permutaci\[OAcute]n *)
GetCoinPermutations[] := <|
  "Px" -> {1, 2, 4, 3}, (* x->-x: Right<->Left *)
  "Py" -> {2, 1, 3, 4}  (* y->-y: Up<->Down *)
|>;

GetSectorSigns[sector_String] := Switch[sector,
  "A1", {1, 1},   "A2", {-1, -1},
  "B1", {1, -1},  "B2", {-1, 1},
  _, Print["Sector inv\[AAcute]lido"]; {$Failed, $Failed}
];

GetSymmetrySectorBasis[gridData_Association, sector_String] := Module[
  {
    coords, map, dim,
    signX, signY,
    permCoinX, permCoinY,
    seeds, rules, 
    colCounter, V
  },

  coords = gridData["Coords"];
  map = gridData["Mapping"];
  dim = gridData["Dimension"];
  
  {signX, signY} = GetSectorSigns[sector];
  If[FailureQ[signX], Return[$Failed]];

  {permCoinX, permCoinY} = Values[GetCoinPermutations[]];

  (* 1. Seleccionar Semillas (Dominio Fundamental: 1er Cuadrante) *)
  seeds = Select[coords, #[[1]] >= 0 && #[[2]] >= 0 &];

  colCounter = 0;
  
  rules = Reap[
    Do[
      Module[{
        seedCoord = seedPos,
        seedIdx = map[seedPos],
        localCandidates = {}
      },
        
        (* 2. Generar vectores candidatos para los 4 estados de la moneda *)
        Do[
            Module[{rawVec = {}},
                (* T1: Identidad *)
                AppendTo[rawVec, {seedIdx, coin} -> 1];
                
                (* T2: Px *)
                With[{p = {-seedCoord[[1]], seedCoord[[2]]}, c = permCoinX[[coin]]},
                    If[!MissingQ[map[p]], AppendTo[rawVec, {map[p], c} -> signX]]
                ];
                
                (* T3: Py *)
                With[{p = {seedCoord[[1]], -seedCoord[[2]]}, c = permCoinY[[coin]]},
                    If[!MissingQ[map[p]], AppendTo[rawVec, {map[p], c} -> signY]]
                ];
                
                (* T4: PxPy *)
                With[{p = {-seedCoord[[1]], -seedCoord[[2]]}, c = permCoinY[[permCoinX[[coin]]]]},
                    If[!MissingQ[map[p]], AppendTo[rawVec, {map[p], c} -> signX * signY]]
                ];
                
                (* Consolidar amplitudes: {key, val} *)
                rawVec = GatherBy[rawVec, First];
                rawVec = {First[#][[1]], Total[Last /@ #]} & /@ rawVec;
                
                (* Filtro num\[EAcute]rico b\[AAcute]sico *)
                rawVec = Select[rawVec, Abs[Last[#]] > 10^-8 &];
                
                If[Length[rawVec] > 0, AppendTo[localCandidates, rawVec]];
            ],
            {coin, 1, 4}
        ];

        (* 3. Ortogonalizaci\[OAcute]n Local (Correcci\[OAcute]n de BUG) *)
        If[Length[localCandidates] > 0,
            Module[{allKeys, denseMat, orthoDense, restoredRules},
                
                (* CORRECCI\[CapitalOAcute]N CR\[CapitalIAcute]TICA: Flatten nivel 1 para preservar llaves {idx, c} *)
                allKeys = Union[Flatten[localCandidates[[All, All, 1]], 1]];
                
                (* Construir matriz densa peque\[NTilde]a (aprox 4x16) *)
                (* Rule @@@ vec convierte {{k,v}...} a {k->v...} *)
                denseMat = Table[
                    Lookup[Association[Rule @@@ vec], allKeys, 0.],
                    {vec, localCandidates}
                ];
                
                (* Ortogonalizar sin Method espec\[IAcute]fico (usa default robusto) *)
                orthoDense = Orthogonalize[denseMat, Tolerance -> 10^-6];
                (* Eliminar vectores nulos resultantes *)
                orthoDense = Select[orthoDense, Norm[#] > 10^-6 &];
                
                (* Restaurar a reglas dispersas y sembrar *)
                Do[
                    colCounter++;
                    restoredRules = Transpose[{allKeys, orthoDense[[k]]}];
                    restoredRules = Select[restoredRules, Abs[Last[#]] > 10^-8 &];
                    
                    Sow[
                        ({4*(#[[1, 1]] - 1) + #[[1, 2]], colCounter} -> #[[2]]) & /@ restoredRules
                    ];
                , {k, Length[orthoDense]}]
            ]
        ];
      ],
      {seedPos, seeds}
    ]
  ][[2]]; (* Parte 2 del Reap contiene las listas sembradas *)

  (* Flatten[rules] fusiona las listas de reglas sembradas en cada paso *)
  V = SparseArray[Flatten[rules], {4 * dim, colCounter}];
  V
];

End[];
EndPackage[];
