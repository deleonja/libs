(* ::Package:: *)

(* QMB Initialization File  *)
Module[{RootPath, PacletPath, PacletData, VersionLocal, 
    VersionRemote, GitHubUrl, CheckResult},
    
    RootPath = DirectoryName[DirectoryName[$InputFileName]];
    PacletPath = FileNameJoin[{RootPath, "PacletInfo.wl"}];
    GitHubUrl = "https://raw.githubusercontent.com/deleonja/libs/main/" <> 
        "QMB/version.txt?timestamp=" <> ToString[IntegerPart[AbsoluteTime[]]];

    (* 1. Handling ForScience paclet installation *)
    If[Length[PacletFind["ForScience"]] == 0, 
        PacletInstall[FileNameJoin[{RootPath, "Kernel", 
            "ForScience-0.88.45.paclet"}]]
    ];

    (* 2. Read version directly from disk file to avoid cache issues *)
    VersionLocal = If[FileExistsQ[PacletPath],
        PacletData = Import[PacletPath];
        (* List @@ converts Paclet[a->b, c->d] to {a->b, c->d} *)
        Version /. List @@ PacletData,
        "0.0.0"
    ]; 

    (* 3. Check for updates with 0.5s timeout *)
    CheckResult = Quiet[
		TimeConstrained[URLRead[GitHubUrl, "Body"], 3, $Failed]
	];

    (* 4. Notify user *)
    If[StringQ[CheckResult],
        VersionRemote = StringTrim[CheckResult];
        If[VersionLocal =!= VersionRemote,
            Print["[QMB Update] A newer version is available (Local: ", 
                VersionLocal, ", Remote: ", VersionRemote, ")."];
            Print["[QMB Update] Pull updates from GitHub for newest version."]
        ]
    ];

	(* 5. Load submodules with error handling *)
    Scan[
        Function[file,
            If[FileExistsQ[file],
                Quiet[Get[file]],
                Print["[QMB Warning] Submodule not found: ", file]
            ]
        ],
        {
            FileNameJoin[{RootPath, "OldQMB.wl"}],
            FileNameJoin[{RootPath, "GeneralQM.wl"}],
            FileNameJoin[{RootPath, "RMT.wl"}],
            FileNameJoin[{RootPath, "QKT.wl"}],
            FileNameJoin[{RootPath, "ManyBody", "SpinChains.wl"}],
            FileNameJoin[{RootPath, "ManyBody", "BoseHubbard.wl"}],
            FileNameJoin[{RootPath, "ManyBody", "Fermions.wl"}]
        }
    ]
];
