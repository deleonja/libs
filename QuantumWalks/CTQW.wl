(* ::Package:: *)

BeginPackage["QuantumWalks`CTQW`", {"QuantumWalks`"}];


BuildAdjacencyMatrix::usage = "BuildAdjacencyMatrix[gridData, VectoresConexion] construye \
la matriz de adyacencia usando b\[UAcute]squeda O(1). \
VectoresConexion por defecto asume vecinos pr\[OAcute]ximos en una red cuadrada 2D.";

CTQWHamiltonian::usage = "CTQWHamiltonian[gridData, \[Gamma], VectoresConexion] construye \
el Hamiltoniano H = -\[Gamma] L, donde L = D - A es el Laplaciano del grafo. \
\[Gamma] es el factor de salto (por defecto 1.0).";


Begin["`Private`"];

(* Si no se proveen vectores, usamos vecinos pr\[OAcute]ximos por defecto *)
OpcionesConexionPorDefecto = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

BuildAdjacencyMatrix[GridData_Association, VectoresConexion_List : OpcionesConexionPorDefecto] := 
    Module[{Coordenadas, Mapeo, DimensionGrid, ReglasAdyacencia},
        
        Coordenadas = GridData["Coords"];
        Mapeo = GridData["Mapping"];
        DimensionGrid = GridData["Dimension"];
        
        ReglasAdyacencia = Flatten @ MapIndexed[
            Module[{IndiceActual = #2[[1]], PosicionActual = #1, PosicionesVecinos, IndicesVecinos},
                
                (* 1. Calcular coordenadas te\[OAcute]ricas de los vecinos *)
                PosicionesVecinos = (# + PosicionActual) & /@ VectoresConexion;
                
                (* 2. Buscar en el diccionario O(1) y eliminar los faltantes *)
                IndicesVecinos = DeleteMissing[Mapeo /@ PosicionesVecinos];
                
                (* 3. Generar las reglas correctamente usando Map en lugar de Thread *)
                ({IndiceActual, #} -> 1.) & /@ IndicesVecinos
            ] &,
            Coordenadas
        ];
        
        (* Retornar la matriz *)
        SparseArray[ReglasAdyacencia, {DimensionGrid, DimensionGrid}]
    ];
    
CTQWHamiltonian[GridData_, Gamma_ : 1.0, VectoresConexion_List : OpcionesConexionPorDefecto] := 
    Module[{MatrizA, MatrizD, Laplaciano},
        
        (* 1. Obtener la matriz de adyacencia usando la funci\[OAcute]n previa *)
        MatrizA = BuildAdjacencyMatrix[GridData, VectoresConexion];
        
        (* 2. Calcular los grados de los v\[EAcute]rtices sumando las filas de A (dimensi\[OAcute]n 2) *)
        (* DiagonalMatrix reconoce autom\[AAcute]ticamente vectores dispersos y crea una matriz dispersa *)
        MatrizD = DiagonalMatrix[SparseArray[Total[MatrizA, {2}]]];
        
        (* 3. L = D - A. El Hamiltoniano es H = -\[Gamma] L = \[Gamma] (A - D) *)
        N[Gamma * (MatrizA - MatrizD)]
    ]; 

End[];
EndPackage[];
