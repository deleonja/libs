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


BuildSymmetryIsometry::usage = 
"BuildSymmetryIsometry[gridData, \"sector\"] construye la isometria V \
para el subespacio \"sector\" (\"A1\", \"B2\", \"Even\", etc.).
BuildSymmetryIsometry[gridData, {sx, sy}] construye la isometria para \
simetria C2v (2 ejes), con paridades espaciales sx y sy (1 o -1).
BuildSymmetryIsometry[gridData, sy] construye la isometria para Cs \
(eje Y), con paridad sy (1 o -1).
BuildSymmetryIsometry[gridData, chars] construye la isometria para C4v \
(8 isometrias) usando una Association de 8 caracteres.";


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


BuildSymmetryIsometry::invSector = "Sector invalido. Usa etiquetas \
(A1, Even, C4v_B1), signos (1, -1) o tuplas ({1, -1}).";

(* ---------------------------------------------------------------------- *)
(* 1. INTERFAZ POLIMORFICA UNIFICADA (STRINGS)                            *)
(* ---------------------------------------------------------------------- *)
BuildSymmetryIsometry[gridData_Association, sector_String] := Module[
    {SignosTraducidos},
    
    SignosTraducidos = Lookup[<|
        "A1" -> {1, 1},   "++" -> {1, 1},
        "A2" -> {-1, -1}, "--" -> {-1, -1},
        "B1" -> {1, -1},  "+-" -> {1, -1},
        "B2" -> {-1, 1},  "-+" -> {-1, 1},
        "Even" -> 1, "Odd" -> -1, "+" -> 1, "-" -> -1,
        
        "C4v_A1" -> <|"E"->1, "Px"->1,  "Py"->1,  "Pxy"->1,  "Pd"->1,  "Pad"->1,  "C4"->1,  "C4inv"->1|>,
        "C4v_A2" -> <|"E"->1, "Px"->-1, "Py"->-1, "Pxy"->1,  "Pd"->-1, "Pad"->-1, "C4"->1,  "C4inv"->1|>,
        "C4v_B1" -> <|"E"->1, "Px"->1,  "Py"->1,  "Pxy"->1,  "Pd"->-1, "Pad"->-1, "C4"->-1, "C4inv"->-1|>,
        "C4v_B2" -> <|"E"->1, "Px"->-1, "Py"->-1, "Pxy"->1,  "Pd"->1,  "Pad"->1,  "C4"->-1, "C4inv"->-1|>,
        "C4v_Ex" -> {1, -1},
        "C4v_Ey" -> {-1, 1}
    |>, sector, $Failed];
    
    If[FailureQ[SignosTraducidos], 
        Message[BuildSymmetryIsometry::invSector]; Return[$Failed]
    ];
    
    BuildSymmetryIsometry[gridData, SignosTraducidos]
];

(* ------------------------------------------------------------------------- *)
(* 2. L\[CapitalOAcute]GICA MATEM\[CapitalAAcute]TICA: SIMETR\[CapitalIAcute]A 2 EJES (C2v)                               *)
(* ------------------------------------------------------------------------- *)
BuildSymmetryIsometry[gridData_Association, {SignoX_Integer, SignoY_Integer}] := Module[
    {
        Coordenadas = gridData["Coords"], 
        Mapeo = gridData["Mapping"], 
        Dim = gridData["Dimension"],
        PermMonedaX = {1, 2, 4, 3}, PermMonedaY = {2, 1, 3, 4},
        Semillas, ReglasV, ContadorColumnas = 0, ReglasLimpias
    },
    
    Semillas = Select[Coordenadas, #[[1]] >= 0 && #[[2]] >= 0 &];
    
    ReglasV = Reap[
        Do[
            Module[{IndiceSemilla = Mapeo[PosicionSemilla], CandidatosLocales = {}},
                Do[
                    Module[{VectorCrudo = {}, Px, Py, Cx, Cy, Cxy},
                        Px = {-PosicionSemilla[[1]], PosicionSemilla[[2]]};
                        Py = {PosicionSemilla[[1]], -PosicionSemilla[[2]]};
                        Cx = PermMonedaX[[EstadoMoneda]];
                        Cy = PermMonedaY[[EstadoMoneda]];
                        Cxy = PermMonedaY[[Cx]];
                        
                        AppendTo[VectorCrudo, {IndiceSemilla, EstadoMoneda} -> 1.];
                        If[!MissingQ[Mapeo[Px]], AppendTo[VectorCrudo, {Mapeo[Px], Cx} -> SignoX]];
                        If[!MissingQ[Mapeo[Py]], AppendTo[VectorCrudo, {Mapeo[Py], Cy} -> SignoY]];
                        If[!MissingQ[Mapeo[-PosicionSemilla]], AppendTo[VectorCrudo, {Mapeo[-PosicionSemilla], Cxy} -> SignoX * SignoY]];
                        
                        VectorCrudo = {First[#][[1]], Total[Last /@ #]} & /@ GatherBy[VectorCrudo, First];
                        VectorCrudo = Select[VectorCrudo, Abs[Last[#]] > 10^-8 &];
                        If[Length[VectorCrudo] > 0, AppendTo[CandidatosLocales, VectorCrudo]];
                    ], {EstadoMoneda, 1, 4}
                ];
                ContadorColumnas = OrtogonalizarYSembrar[CandidatosLocales, ContadorColumnas];
            ], {PosicionSemilla, Semillas}
        ]
    ][[2]];
    
    ReglasLimpias = Cases[Flatten[ReglasV], _Rule];
    If[Length[ReglasLimpias] == 0, Return[SparseArray[{}, {4 * Dim, 0}]]];
    SparseArray[ReglasLimpias, {4 * Dim, ContadorColumnas}]
];

(* ------------------------------------------------------------------------- *)
(* 3. L\[CapitalOAcute]GICA MATEM\[CapitalAAcute]TICA: SIMETR\[CapitalIAcute]A 1 EJE (Bilateral)                          *)
(* ------------------------------------------------------------------------- *)
BuildSymmetryIsometry[gridData_Association, SignoY_Integer] := Module[
    {
        Coordenadas = gridData["Coords"], 
        Mapeo = gridData["Mapping"], 
        Dim = gridData["Dimension"],
        PermMonedaY = {2, 1, 3, 4},
        Semillas, ReglasV, ContadorColumnas = 0, ReglasLimpias
    },
    
    Semillas = Select[Coordenadas, #[[2]] >= 0 &];
    
    ReglasV = Reap[
        Do[
            Module[{IndiceSemilla = Mapeo[PosicionSemilla], CandidatosLocales = {}},
                Do[
                    Module[{VectorCrudo = {}, Py, Cy},
                        Py = {PosicionSemilla[[1]], -PosicionSemilla[[2]]};
                        Cy = PermMonedaY[[EstadoMoneda]];
                        
                        AppendTo[VectorCrudo, {IndiceSemilla, EstadoMoneda} -> 1.];
                        If[!MissingQ[Mapeo[Py]], AppendTo[VectorCrudo, {Mapeo[Py], Cy} -> SignoY]];
                        
                        VectorCrudo = {First[#][[1]], Total[Last /@ #]} & /@ GatherBy[VectorCrudo, First];
                        VectorCrudo = Select[VectorCrudo, Abs[Last[#]] > 10^-8 &];
                        If[Length[VectorCrudo] > 0, AppendTo[CandidatosLocales, VectorCrudo]];
                    ], {EstadoMoneda, 1, 4}
                ];
                ContadorColumnas = OrtogonalizarYSembrar[CandidatosLocales, ContadorColumnas];
            ], {PosicionSemilla, Semillas}
        ]
    ][[2]];
    
    ReglasLimpias = Cases[Flatten[ReglasV], _Rule];
    If[Length[ReglasLimpias] == 0, Return[SparseArray[{}, {4 * Dim, 0}]]];
    SparseArray[ReglasLimpias, {4 * Dim, ContadorColumnas}]
];

(* ========================================================================= *)
(* 4. L\[CapitalOAcute]GICA MATEM\[CapitalAAcute]TICA: SIMETR\[CapitalIAcute]A C4v (8 Isometr\[IAcute]as como el Billar de Sina\[IAcute]) *)
(* ========================================================================= *)
BuildSymmetryIsometry[gridData_Association, chars_Association] /; Length[chars] == 8 := Module[
    {
        Coordenadas = gridData["Coords"], 
        Mapeo = gridData["Mapping"], 
        Dim = gridData["Dimension"],
        PermMoneda, Semillas, ReglasV, ContadorColumnas = 0, ReglasLimpias
    },
    
    PermMoneda = <|
        "E"     -> {1, 2, 3, 4}, "Px"  -> {1, 2, 4, 3}, 
        "Py"    -> {2, 1, 3, 4}, "Pxy" -> {2, 1, 4, 3},
        "Pd"    -> {3, 4, 1, 2}, "Pad" -> {4, 3, 2, 1}, 
        "C4"    -> {4, 3, 1, 2}, "C4inv"->{3, 4, 2, 1}
    |>;
    
    Semillas = Select[
        Coordenadas, 
        #[[1]] >= 0 && #[[2]] >= 0 && #[[2]] <= #[[1]] &
    ];
    
    ReglasV = Reap[
        Map[
            Function[{PosicionSemilla},
                Module[{CandidatosLocales},
                    CandidatosLocales = DeleteMissing[
                        Map[
                            Function[{EstadoMoneda},
                                GenerarVectorC4v[PosicionSemilla, EstadoMoneda, Mapeo, PermMoneda, chars]
                            ],
                            {1, 2, 3, 4}
                        ]
                    ];
                    ContadorColumnas = OrtogonalizarYSembrar[CandidatosLocales, ContadorColumnas];
                ]
            ],
            Semillas
        ]
    ][[2]];
    
    (* El escudo definitivo: Cases garantiza que SparseArray NUNCA reciba artefactos *)
    ReglasLimpias = Cases[Flatten[ReglasV], _Rule];
    If[Length[ReglasLimpias] == 0, Return[SparseArray[{}, {4 * Dim, 0}]]];
    SparseArray[ReglasLimpias, {4 * Dim, ContadorColumnas}]
];

(* --- Funci\[OAcute]n pura delegada para generar vectores sin AppendTo --- *)
GenerarVectorC4v[Posicion_, Moneda_, Mapeo_, PermMoneda_, chars_] := 
  Module[
    {VectorCrudo, TransEspaciales},
    
    TransEspaciales = <|
        "E"     -> {Posicion[[1]], Posicion[[2]]},
        "Px"    -> {-Posicion[[1]], Posicion[[2]]},
        "Py"    -> {Posicion[[1]], -Posicion[[2]]},
        "Pxy"   -> {-Posicion[[1]], -Posicion[[2]]},
        "Pd"    -> {Posicion[[2]], Posicion[[1]]},
        "Pad"   -> {-Posicion[[2]], -Posicion[[1]]},
        "C4"    -> {-Posicion[[2]], Posicion[[1]]},
        "C4inv" -> {Posicion[[2]], -Posicion[[1]]}
    |>;
    
    VectorCrudo = DeleteMissing[
        Map[
            Function[{Operacion},
                Module[{IndiceTrans, MonedaTrans, Caracter},
                    IndiceTrans = Mapeo[TransEspaciales[Operacion]];
                    MonedaTrans = PermMoneda[Operacion][[Moneda]];
                    Caracter = chars[Operacion];
                    
                    If[MissingQ[IndiceTrans], Missing[], {IndiceTrans, MonedaTrans} -> Caracter]
                ]
            ],
            Keys[chars]
        ]
    ];
    
    VectorCrudo = {First[#][[1]], Total[Last /@ #]} & /@ GatherBy[VectorCrudo, First];
    VectorCrudo = Select[VectorCrudo, Abs[Last[#]] > 10^-8 &];
    
    If[Length[VectorCrudo] > 0, VectorCrudo, Missing[]]
];

(* --- Funci\[OAcute]n Auxiliar Blindada --- *)
OrtogonalizarYSembrar[CandidatosLocales_List, ContadorActual_Integer] := Module[
    {LlavesGlobales, MatrizDensa, MatrizOrto, NuevoContador = ContadorActual},
    
    If[Length[CandidatosLocales] == 0, Return[NuevoContador]];
    
    LlavesGlobales = Union[Flatten[CandidatosLocales[[All, All, 1]], 1]];
    MatrizDensa = Table[Lookup[Association[Rule @@@ vec], LlavesGlobales, 0.], {vec, CandidatosLocales}];
    
    (* Evita procesar si la matriz es completamente nula (ruido de m\[AAcute]quina) *)
    If[Total[Abs[Flatten[MatrizDensa]]] < 10^-8, Return[NuevoContador]];
    
    MatrizOrto = Select[Orthogonalize[MatrizDensa, Tolerance -> 10^-6], Norm[#] > 10^-6 &];
    
    Do[
        NuevoContador++;
        Sow[
            (* Chop remueve partes imaginarias o ruidos min\[UAcute]sculos garantizando n\[UAcute]meros puros *)
            ({4*(#[[1, 1]] - 1) + #[[1, 2]], NuevoContador} -> Chop[#[[2]]]) & /@ 
            Select[Transpose[{LlavesGlobales, MatrizOrto[[k]]}], Abs[Last[#]] > 10^-8 &]
        ];
    , {k, Length[MatrizOrto]}];
    
    NuevoContador
];


(*BuildSymmetryIsometry::invSector = "Sector invalido. Usa etiquetas \
(A1, Even, C4v_B1), signos (1, -1) o tuplas ({1, -1}).";

(* ---------------------------------------------------------------------- *)
(* 1. INTERFAZ POLIMORFICA UNIFICADA (STRINGS)                            *)
(* Traduce todas las etiquetas a representaciones matematicas             *)
(* ---------------------------------------------------------------------- *)
BuildSymmetryIsometry[gridData_Association, sector_String] := Module[
    {SignosTraducidos},
    
    SignosTraducidos = Lookup[<|
        (* Sectores de 2 ejes (C2v) *)
        "A1" -> {1, 1},   "++" -> {1, 1},
        "A2" -> {-1, -1}, "--" -> {-1, -1},
        "B1" -> {1, -1},  "+-" -> {1, -1},
        "B2" -> {-1, 1},  "-+" -> {-1, 1},
        
        (* Sectores de 1 eje (Bilateral en Y) *)
        "Even" -> 1, "Odd" -> -1, "+" -> 1, "-" -> -1,
        
        (* Sectores de 4 ejes (C4v - 8 isometrias) *)
        "C4v_A1" -> <|"E"->1, "Px"->1,  "Py"->1,  "Pxy"->1, 
                      "Pd"->1,  "Pad"->1,  "C4"->1,  "C4inv"->1|>,
        "C4v_A2" -> <|"E"->1, "Px"->-1, "Py"->-1, "Pxy"->1, 
                      "Pd"->-1, "Pad"->-1, "C4"->1,  "C4inv"->1|>,
        "C4v_B1" -> <|"E"->1, "Px"->1,  "Py"->1,  "Pxy"->1, 
                      "Pd"->-1, "Pad"->-1, "C4"->-1, "C4inv"->-1|>,
        "C4v_B2" -> <|"E"->1, "Px"->-1, "Py"->-1, "Pxy"->1, 
                      "Pd"->1,  "Pad"->1,  "C4"->-1, "C4inv"->-1|>,
        "C4v_E"  -> <|"E"->2, "Px"->0,  "Py"->0,  "Pxy"->-2, 
                      "Pd"->0,  "Pad"->0,  "C4"->0,  "C4inv"->0|>
    |>, sector, $Failed];
    
    If[FailureQ[SignosTraducidos], 
        Message[BuildSymmetryIsometry::invSector]; Return[$Failed]
    ];
    
    (* Recursion polimorfica: llama a la logica matematica exacta *)
    BuildSymmetryIsometry[gridData, SignosTraducidos]
];

(* ------------------------------------------------------------------------- *)
(* 2. L\[CapitalOAcute]GICA MATEM\[CapitalAAcute]TICA: SIMETR\[CapitalIAcute]A 2 EJES (C2v)                               *)
(* Se activa autom\[AAcute]ticamente si el sector es una lista de 2 enteros          *)
(* ------------------------------------------------------------------------- *)
BuildSymmetryIsometry[gridData_Association, {SignoX_Integer, SignoY_Integer}] := Module[
    {
        Coordenadas = gridData["Coords"], 
        Mapeo = gridData["Mapping"], 
        Dim = gridData["Dimension"],
        PermMonedaX = {1, 2, 4, 3}, PermMonedaY = {2, 1, 3, 4},
        Semillas, ReglasV, ContadorColumnas = 0
    },
    
    (* Dominio Fundamental: Primer Cuadrante *)
    Semillas = Select[Coordenadas, #[[1]] >= 0 && #[[2]] >= 0 &];
    
    ReglasV = Reap[
        Do[
            Module[{IndiceSemilla = Mapeo[PosicionSemilla], CandidatosLocales = {}},
                Do[
                    Module[{VectorCrudo = {}, Px, Py, Cx, Cy, Cxy},
                        (* Transformaciones C2v *)
                        Px = {-PosicionSemilla[[1]], PosicionSemilla[[2]]};
                        Py = {PosicionSemilla[[1]], -PosicionSemilla[[2]]};
                        Cx = PermMonedaX[[EstadoMoneda]];
                        Cy = PermMonedaY[[EstadoMoneda]];
                        Cxy = PermMonedaY[[Cx]];
                        
                        AppendTo[VectorCrudo, {IndiceSemilla, EstadoMoneda} -> 1.];
                        If[!MissingQ[Mapeo[Px]], AppendTo[VectorCrudo, {Mapeo[Px], Cx} -> SignoX]];
                        If[!MissingQ[Mapeo[Py]], AppendTo[VectorCrudo, {Mapeo[Py], Cy} -> SignoY]];
                        If[!MissingQ[Mapeo[-PosicionSemilla]], AppendTo[VectorCrudo, {Mapeo[-PosicionSemilla], Cxy} -> SignoX * SignoY]];
                        
                        VectorCrudo = {First[#][[1]], Total[Last /@ #]} & /@ GatherBy[VectorCrudo, First];
                        VectorCrudo = Select[VectorCrudo, Abs[Last[#]] > 10^-8 &];
                        If[Length[VectorCrudo] > 0, AppendTo[CandidatosLocales, VectorCrudo]];
                    ], {EstadoMoneda, 1, 4}
                ];
                
                (* Ortogonalizaci\[OAcute]n local delegada a rutina auxiliar (para no repetir c\[OAcute]digo) *)
                ContadorColumnas = OrtogonalizarYSembrar[CandidatosLocales, ContadorColumnas];
            ], {PosicionSemilla, Semillas}
        ]
    ][[2]];
    
    SparseArray[Flatten[ReglasV], {4 * Dim, ContadorColumnas}]
];

(* ------------------------------------------------------------------------- *)
(* 3. L\[CapitalOAcute]GICA MATEM\[CapitalAAcute]TICA: SIMETR\[CapitalIAcute]A 1 EJE (Bilateral)                          *)
(* Se activa autom\[AAcute]ticamente si el sector es un solo entero                  *)
(* ------------------------------------------------------------------------- *)
BuildSymmetryIsometry[gridData_Association, SignoY_Integer] := Module[
    {
        Coordenadas = gridData["Coords"], 
        Mapeo = gridData["Mapping"], 
        Dim = gridData["Dimension"],
        PermMonedaY = {2, 1, 3, 4},
        Semillas, ReglasV, ContadorColumnas = 0
    },
    
    (* Dominio Fundamental: Semiplano Superior *)
    Semillas = Select[Coordenadas, #[[2]] >= 0 &];
    
    ReglasV = Reap[
        Do[
            Module[{IndiceSemilla = Mapeo[PosicionSemilla], CandidatosLocales = {}},
                Do[
                    Module[{VectorCrudo = {}, Py, Cy},
                        (* Transformaciones Cs *)
                        Py = {PosicionSemilla[[1]], -PosicionSemilla[[2]]};
                        Cy = PermMonedaY[[EstadoMoneda]];
                        
                        AppendTo[VectorCrudo, {IndiceSemilla, EstadoMoneda} -> 1.];
                        If[!MissingQ[Mapeo[Py]], AppendTo[VectorCrudo, {Mapeo[Py], Cy} -> SignoY]];
                        
                        VectorCrudo = {First[#][[1]], Total[Last /@ #]} & /@ GatherBy[VectorCrudo, First];
                        VectorCrudo = Select[VectorCrudo, Abs[Last[#]] > 10^-8 &];
                        If[Length[VectorCrudo] > 0, AppendTo[CandidatosLocales, VectorCrudo]];
                    ], {EstadoMoneda, 1, 4}
                ];
                
                (* Ortogonalizaci\[OAcute]n local delegada *)
                ContadorColumnas = OrtogonalizarYSembrar[CandidatosLocales, ContadorColumnas];
            ], {PosicionSemilla, Semillas}
        ]
    ][[2]];
    
    If[ReglasV === {}, Return[SparseArray[{}, {4 * Dim, 0}]]];
    SparseArray[Flatten[ReglasV], {4 * Dim, ContadorColumnas}]
];

(* ========================================================================= *)
(* 4. L\[CapitalOAcute]GICA MATEM\[CapitalAAcute]TICA: SIMETR\[CapitalIAcute]A C4v (8 Isometr\[IAcute]as como el Billar de Sina\[IAcute]) *)
(* Se activa si recibe una Association con exactamente 8 caracteres          *)
(* ========================================================================= *)

BuildSymmetryIsometry[gridData_Association, chars_Association] /; 
  Length[chars] == 8 := Module[
    {
        Coordenadas = gridData["Coords"], 
        Mapeo = gridData["Mapping"], 
        Dim = gridData["Dimension"],
        PermMoneda, Semillas, ReglasV, ContadorColumnas = 0
    },
    
    PermMoneda = <|
        "E"     -> {1, 2, 3, 4}, "Px"  -> {1, 2, 4, 3}, 
        "Py"    -> {2, 1, 3, 4}, "Pxy" -> {2, 1, 4, 3},
        "Pd"    -> {3, 4, 1, 2}, "Pad" -> {4, 3, 2, 1}, 
        "C4"    -> {4, 3, 1, 2}, "C4inv"->{3, 4, 2, 1}
    |>;
    
    Semillas = Select[
        Coordenadas, 
        #[[1]] >= 0 && #[[2]] >= 0 && #[[2]] <= #[[1]] &
    ];
    
    ReglasV = Reap[
        Map[
            Function[{PosicionSemilla},
                Module[
                    {CandidatosLocales},
                    
                    (* Flujo puro: Mapeamos los 4 estados y eliminamos nulos *)
                    CandidatosLocales = DeleteMissing[
                        Map[
                            Function[{EstadoMoneda},
                                GenerarVectorC4v[
                                    PosicionSemilla, 
                                    EstadoMoneda, 
                                    Mapeo, 
                                    PermMoneda, 
                                    chars
                                ]
                            ],
                            {1, 2, 3, 4}
                        ]
                    ];
                    
                    ContadorColumnas = OrtogonalizarYSembrar[
                        CandidatosLocales, 
                        ContadorColumnas
                    ];
                ]
            ],
            Semillas
        ]
    ][[2]];
    
    SparseArray[Flatten[ReglasV], {4 * Dim, ContadorColumnas}]
];

(* --- Funci\[OAcute]n pura delegada para generar vectores sin AppendTo --- *)
GenerarVectorC4v[Posicion_, Moneda_, Mapeo_, PermMoneda_, chars_] := 
  Module[
    {VectorCrudo, TransEspaciales},
    
    TransEspaciales = <|
        "E"     -> {Posicion[[1]], Posicion[[2]]},
        "Px"    -> {-Posicion[[1]], Posicion[[2]]},
        "Py"    -> {Posicion[[1]], -Posicion[[2]]},
        "Pxy"   -> {-Posicion[[1]], -Posicion[[2]]},
        "Pd"    -> {Posicion[[2]], Posicion[[1]]},
        "Pad"   -> {-Posicion[[2]], -Posicion[[1]]},
        "C4"    -> {-Posicion[[2]], Posicion[[1]]},
        "C4inv" -> {Posicion[[2]], -Posicion[[1]]}
    |>;
    
    VectorCrudo = DeleteMissing[
        Map[
            Function[{Operacion},
                Module[
                    {IndiceTrans, MonedaTrans, Caracter},
                    IndiceTrans = Mapeo[TransEspaciales[Operacion]];
                    MonedaTrans = PermMoneda[Operacion][[Moneda]];
                    Caracter = chars[Operacion];
                    
                    If[MissingQ[IndiceTrans],
                        Missing[],
                        {IndiceTrans, MonedaTrans} -> Caracter
                    ]
                ]
            ],
            Keys[chars]
        ]
    ];
    
    (* Agrupamos amplitudes de sitios que coinciden *)
    VectorCrudo = {First[#][[1]], Total[Last /@ #]} & /@ 
        GatherBy[VectorCrudo, First];
        
    VectorCrudo = Select[VectorCrudo, Abs[Last[#]] > 10^-8 &];
    
    If[Length[VectorCrudo] > 0, VectorCrudo, Missing[]]
];

(* --- Funci\[OAcute]n Auxiliar para la Ortogonalizaci\[OAcute]n y Siembra (DRY) --- *)
OrtogonalizarYSembrar[CandidatosLocales_List, ContadorActual_Integer] := Module[
    {LlavesGlobales, MatrizDensa, MatrizOrto, NuevoContador = ContadorActual},
    
    If[Length[CandidatosLocales] == 0, Return[NuevoContador]];
    
    LlavesGlobales = Union[Flatten[CandidatosLocales[[All, All, 1]], 1]];
    MatrizDensa = Table[Lookup[Association[Rule @@@ vec], LlavesGlobales, 0.], {vec, CandidatosLocales}];
    MatrizOrto = Select[Orthogonalize[MatrizDensa, Tolerance -> 10^-6], Norm[#] > 10^-6 &];
    
    Do[
        NuevoContador++;
        Sow[
            ({4*(#[[1, 1]] - 1) + #[[1, 2]], NuevoContador} -> #[[2]]) & /@ 
            Select[Transpose[{LlavesGlobales, MatrizOrto[[k]]}], Abs[Last[#]] > 10^-8 &]
        ];
    , {k, Length[MatrizOrto]}];
    
    NuevoContador
];*)


End[];
EndPackage[];
