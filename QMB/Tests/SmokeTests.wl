(* ::Package:: *)

(* ============================================================
   SmokeTests.wl
   Auto-discovers all public QMB functions and checks they
   exist and have ::usage defined.
   ============================================================ *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "TestUtils.wl"}]]
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "Kernel", "init.m"}]]

publicFunctions = Names["QMB`*"];

TestSection["Smoke \[LongDash] package exports functions"]
VerifyProperty[
    Length[publicFunctions] > 0,
    "QMB exports at least one function"
]

TestSection["Smoke \[LongDash] all functions have ::usage"]
Scan[
    Function[fn,
        VerifyProperty[
            StringQ[MessageName[Symbol[fn], "usage"]],
            fn <> " has ::usage defined"
        ]
    ],
    publicFunctions
]

PrintTestSummary[]
