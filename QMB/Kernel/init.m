(* ::Package:: *)

(* QMB Initialization File *)
Module[{RootPath, VersionLocal, VersionRemote, GitHubUrl, CheckResult},
    RootPath = DirectoryName[DirectoryName[$InputFileName]];
    
    (* URL for the raw version text file on GitHub *)
    GitHubUrl = "https://raw.githubusercontent.com/deleonja/libs/main/" <> 
        "QMB/version.txt";

    (* 1. Handling ForScience paclet installation *)
    If[Length[PacletFind["ForScience"]] == 0, 
        PacletInstall[FileNameJoin[{RootPath, "ForScience-0.88.45.paclet"}]]
    ];

    (* 2. Modern Paclet version extraction (Fix for PacletInformation warning) *)
    VersionLocal = PacletObject["QMB"]["Version"];

    (* 3. Check for updates with 0.5s timeout *)
    CheckResult = TimeConstrained[
        URLRead[GitHubUrl, "Body"], 
        0.5, 
        $Failed
    ];

    (* 4. Notify user of discrepancies *)
    If[StringQ[CheckResult],
        VersionRemote = StringTrim[CheckResult];
        If[VersionLocal != VersionRemote,
            Print["[QMB Update] A newer version is available (Local: ", 
                VersionLocal, ", Remote: ", VersionRemote, ")."];
            Print["[QMB Update] Please pull from GitHub to stay up to date."]
        ]
    ];

    (* 5. Load OldQMB core and ManyBody submodules *)
    Get[FileNameJoin[{RootPath, "OldQMB.wl"}]]; 
    Get[FileNameJoin[{RootPath, "ManyBody", "SpinChains.wl"}]];
];
