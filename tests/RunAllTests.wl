(* ::Package:: *)

(* ============================================================
   RunAllTests.wl
   Discovers and runs all *Tests.wl files in the repo.
   Usage: wolframscript -file tests/RunAllTests.wl
   ============================================================ *)

Get[FileNameJoin[{DirectoryName[$InputFileName], 
    "..", "QMB", "Tests", "TestUtils.wl"}]]

Print["============================================================"];
Print["libs Test Suite"];
Print["============================================================"];

(* Busca todos los *Tests.wl en el repo, excluyendo TestUtils.wl *)
testFiles = Select[
    FileNames[
        "*Tests.wl",
        DirectoryName[DirectoryName[$InputFileName]],
        Infinity
    ],
    !StringContainsQ[#, "TestUtils"] && 
    !StringContainsQ[#, "RunAllTests"] &
];

If[Length[testFiles] === 0,
    Print["No test files found."];
    Quit[]
];

Print["Found ", Length[testFiles], " test file(s):"];
Scan[Print["  ", FileNameTake[#]] &, testFiles];

Scan[
    Function[file,
        Print["\n>>> Running: ", FileNameTake[file]];
        Get[file]
    ],
    testFiles
];

PrintTestSummary[]
