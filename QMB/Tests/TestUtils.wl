(* ::Package:: *)

(* ============================================================
   TestUtils.wl
   Utility functions for the QMB test suite.
   ============================================================ *)

BeginPackage["TestUtils`"]

VerifyTest::usage = "VerifyTest[expr, expected, testName] checks \
    that expr === expected and reports the result.";
VerifyProperty::usage = "VerifyProperty[expr, testName] checks that \
    expr is True and reports the result.";
TestSection::usage = "TestSection[name] prints a section header.";
PrintTestSummary::usage = "PrintTestSummary[] prints the total \
    number of passed and failed tests.";

Begin["`Private`"]

$PassCount = 0;
$FailCount = 0;
$FailedTests = {};

TestSection[name_String] :=
    Print["\n--- ", name, " ---"]

VerifyTest[expr_, expected_, testName_String] :=
    Module[{result},
        result = (expr === expected);
        If[result,
            $PassCount++;
            Print["  PASS  ", testName],
            $FailCount++;
            AppendTo[$FailedTests, testName];
            Print["  FAIL  ", testName];
            Print["    Expected: ", expected];
            Print["    Got:      ", expr]
        ]
    ]

VerifyProperty[expr_, testName_String] :=
    VerifyTest[expr, True, testName]

PrintTestSummary[] :=
    Module[{total},
        total = $PassCount + $FailCount;
        Print["\n============================================================"];
        Print["Results: ", $PassCount, "/", total, " passed"];
        If[$FailCount > 0,
            Print["Failed tests:"];
            Scan[Print["  - ", #] &, $FailedTests];
            Print["============================================================"];
            Exit[1],
            Print["All tests passed."];
            Print["============================================================"];
            Exit[0]
        ]
    ]

End[]
EndPackage[]
