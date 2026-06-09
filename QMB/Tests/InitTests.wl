(* ::Package:: *)

(* ============================================================
   InitTests.wl
   Tests for QMB/Kernel/init.m
   Run with: wolframscript -file QMB/Tests/InitTests.wl
   ============================================================ *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "TestUtils.wl"}]]
RootPath = DirectoryName[DirectoryName[$InputFileName]];
PacletPath = FileNameJoin[{RootPath, "PacletInfo.wl"}];
(* ------------------------------------------------------------ *)
TestSection["Step 1 \[LongDash] Local version reading"]
Print["  Checking PacletInfo.wl exists..."];
VerifyProperty[
    FileExistsQ[PacletPath],
    "PacletInfo.wl exists"
]
Print["  Reading local version..."];
VerifyProperty[
    StringMatchQ[
        Module[{data}, Version /. List @@ Import[PacletPath]],
        DigitCharacter.. ~~ "." ~~ DigitCharacter.. ~~ "." ~~ DigitCharacter..
    ],
    "Local version has valid semver format X.Y.Z"
]
(* ------------------------------------------------------------ *)
TestSection["Step 2 \[LongDash] Submodules load correctly"]
submodules = {
    "OldQMB.wl",
    "GeneralQM.wl",
    "RMT.wl",
    "QKT.wl",
    FileNameJoin[{"ManyBody", "SpinChains.wl"}],
    FileNameJoin[{"ManyBody", "BoseHubbard.wl"}],
    FileNameJoin[{"ManyBody", "Fermions.wl"}]
};
Print["  Checking submodules exist on disk..."];
Scan[
    Function[sub,
        Print["    -> ", sub, "..."];
        VerifyProperty[
            FileExistsQ[FileNameJoin[{RootPath, sub}]],
            sub <> " exists on disk"
        ]
    ],
    submodules
]
Print["  Loading init.m..."];
VerifyProperty[
    Quiet[Check[
        Get[FileNameJoin[{RootPath, "Kernel", "init.m"}]],
        $Failed
    ]] =!= $Failed,
    "init.m loads without errors"
]
(* ------------------------------------------------------------ *)
PrintTestSummary[]
