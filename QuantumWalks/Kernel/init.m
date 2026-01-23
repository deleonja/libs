(* ::Package:: *)

(* QuantumWalks Initialization File 
   Se ejecuta autom\[AAcute]ticamente al cargar << QuantumWalks` 
*)

(* QuantumWalks/Kernel/init.m *)
Module[{rootPath},
    rootPath = DirectoryName[DirectoryName[$InputFileName]];
    
    (* 1. Cargar Core *)
    Get[FileNameJoin[{rootPath, "DQWL.wl"}]];
    
    (* 2. Cargar Utilidades *)
    Get[FileNameJoin[{rootPath, "Billiards", "Common.wl"}]];
    
    (* 3. Cargar Geometr\[IAcute]as *)
    Get[FileNameJoin[{rootPath, "Billiards", "Bunimovich.wl"}]];
    Get[FileNameJoin[{rootPath, "Billiards", "Sinai.wl"}]];
    Get[FileNameJoin[{rootPath, "Billiards", "Cardioid.wl"}]];
    Get[FileNameJoin[{rootPath, "Billiards", "Rectangle.wl"}]];
    Get[FileNameJoin[{rootPath, "Billiards", "Ellipse.wl"}]];
];
