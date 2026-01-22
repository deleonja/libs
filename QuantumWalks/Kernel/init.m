(* ::Package:: *)

(* QuantumWalks Initialization File 
   Se ejecuta autom\[AAcute]ticamente al cargar << QuantumWalks` 
*)

Block[{$ContextPath},
    Module[{rootPath},
        (* 1. Obtener la ruta ra\[IAcute]z del paquete *)
        (* Kernel/init.m -> Subir dos niveles -> Ra\[IAcute]z del paquete *)
        rootPath = DirectoryName[DirectoryName[$InputFileName]];
        
        (* 2. Cargar el m\[OAcute]dulo base (1D / General) *)
        Get[FileNameJoin[{rootPath, "DQWL.wl"}]];

        (* 3. Cargar librer\[IAcute]as de Billares (Respetando CamelCase) *)
        
        (* A. Cargar Common.wl primero (l\[OAcute]gica base) *)
        (* Nota: Carpeta "Billiards" y archivo "Common.wl" *)
        Get[FileNameJoin[{rootPath, "Billiards", "Common.wl"}]];
        
        (* B. Cargar Bunimovich.wl (depende de Common) *)
        (* Nota: Carpeta "Billiards" y archivo "Bunimovich.wl" *)
        Get[FileNameJoin[{rootPath, "Billiards", "Bunimovich.wl"}]];
        
        (* Futuros billares: *)
        (* Get[FileNameJoin[{rootPath, "Billiards", "Rectangle.wl"}]]; *)
    ]
];
