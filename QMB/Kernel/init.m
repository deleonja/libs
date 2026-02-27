(* ::Package:: *)

(* QMB Initialization File *)
Module[{RootPath, PacletPath, PacletData, VersionLocal, 
    VersionRemote, GitHubUrl, CheckResult},
    
    RootPath = DirectoryName[DirectoryName[$InputFileName]];
    PacletPath = FileNameJoin[{RootPath, "PacletInfo.wl"}];
    GitHubUrl = "https://raw.githubusercontent.com/deleonja/libs/main/" <> 
        "QMB/version.txt?t=" <> ToString[Round[AbsoluteTime[]]];

    (* 1. Handling ForScience paclet installation *)
    If[Length[PacletFind["ForScience"]] == 0, 
        PacletInstall[FileNameJoin[{RootPath, 
            "ForScience-0.88.45.paclet"}]]
    ];

    (* 2. Read version directly from disk file to avoid cache issues *)
    VersionLocal = If[FileExistsQ[PacletPath],
        PacletData = Import[PacletPath];
        (* List @@ converts Paclet[a->b, c->d] to {a->b, c->d} *)
        Version /. List @@ PacletData,
        "0.0.0"
    ]; Print[VersionLocal];

    (* 3. Check for updates with 0.5s timeout *)
    CheckResult = TimeConstrained[URLRead[GitHubUrl, "Body"], 0.5, $Failed];

    (* 4. Notify user of discrepancies *)
    If[StringQ[CheckResult],
        VersionRemote = StringTrim[CheckResult];
        If[VersionLocal != VersionRemote,
            Print["[QMB Update] A newer version is available (Local: ", 
                VersionLocal, ", Remote: ", VersionRemote, ")."];
            Print["[QMB Update] You may want to pull from GitHub to update."]
        ]
    ];

    (* 5. Load submodules *)
    Get[FileNameJoin[{RootPath, "OldQMB.wl"}]]; 
    Get[FileNameJoin[{RootPath, "ManyBody", "SpinChains.wl"}]];
];
