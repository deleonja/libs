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
TestSection["Step 1 \[LongDash] ForScience paclet"]

Print["  Checking ForScience installation..."];
VerifyProperty[
    Length[PacletFind["ForScience"]] > 0,
    "ForScience is installed after init"
]

(* ------------------------------------------------------------ *)
TestSection["Step 2 \[LongDash] Local version reading"]

Print["  Checking PacletInfo.wl exists..."];
VerifyProperty[
    FileExistsQ[PacletPath],
    "PacletInfo.wl exists"
]

Print["  Reading local version..."];
VerifyTest[
    Module[{data},
        data = Import[PacletPath];
        Version /. List @@ data
    ],
    "2.4.1",
    "Local version reads correctly from PacletInfo.wl"
]

Print["  Checking local version type..."];
VerifyProperty[
    StringQ[Version /. List @@ Import[PacletPath]],
    "Local version is a String"
]

(* ------------------------------------------------------------ *)
TestSection["Step 3 \[LongDash] Remote version check"]

Print["  Fetching remote version.txt (timeout: 3s)..."];
VerifyProperty[
    StringQ[Quiet[TimeConstrained[
        URLRead[
            "https://raw.githubusercontent.com/deleonja/libs/main/QMB/version.txt",
            "Body"
        ],
        3, $Failed
    ]]],
    "Remote version.txt is reachable within 3s"
]

Print["  Validating remote version format..."];
VerifyProperty[
    StringMatchQ[
        StringTrim @ Quiet[TimeConstrained[
            URLRead[
                "https://raw.githubusercontent.com/deleonja/libs/main/QMB/version.txt",
                "Body"
            ],
            3, $Failed
        ]],
        DigitCharacter.. ~~ "." ~~ DigitCharacter.. ~~ "." ~~ DigitCharacter..
    ],
    "Remote version has valid semver format X.Y.Z"
]

(* ------------------------------------------------------------ *)
TestSection["Step 3 \[LongDash] Timeout fallback"]

Print["  Simulating network timeout (may take up to 3s)..."];
VerifyTest[
    Quiet[TimeConstrained[
        URLRead["http://192.0.2.1", "Body"],
        3, $Failed
    ]],
    $Failed,
    "Returns $Failed when timeout is exceeded"
]

(* ------------------------------------------------------------ *)
TestSection["Step 4 \[LongDash] Version comparison"]

Print["  Checking =!= detects mismatch..."];
VerifyProperty[
    ("2.4.1" =!= "2.4.0"),
    "=!= correctly detects version mismatch"
]

Print["  Checking =!= detects match..."];
VerifyProperty[
    !("2.4.1" =!= "2.4.1"),
    "=!= correctly detects matching versions"
]

(* ------------------------------------------------------------ *)
TestSection["Step 5 \[LongDash] Submodules load correctly"]

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
