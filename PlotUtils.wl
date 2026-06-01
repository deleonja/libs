(* ::Package:: *)

BeginPackage["PlotUtils`", {"MaTeX`"}]


ResourceFunction["AddMatplotlibColors"][];


ExportToPDFAndCrop::usage = "ExportToPDFAndCrop[plot, fileName, opts] exporta a PDF y recorta bordes blancos.";

APSPlotStyle::usage = "APSPlotStyle es una lista de opciones predefinidas para figuras de Physical Review.";

ColumnWidth::usage = "Ancho de una columna en puntos (244.8 pt).";

FullWidth::usage = "Ancho de doble columna en puntos (504 pt).";

DefaultFontSize::usage = "Tama\[NTilde]o adecuado para art\[IAcute]culos.";

$PreviewScale::usage = "$PreviewScale es un factor multiplicador global que ajusta 
simult\[AAcute]neamente el ImageSize y el FontSize de las figuras en el Notebook. 
Permite trabajar con gr\[AAcute]ficos m\[AAcute]s grandes en pantalla sin perder la proporci\[OAcute]n visual 
que tendr\[AAcute]n en la publicaci\[OAcute]n final (donde $PreviewScale = 1).";

NiceFigure::usage = "NiceFigure[plot, opts] aplica el estilo APSPlotStyle al gr\[AAcute]fico. 
Convierte autom\[AAcute]ticamente los strings en FrameLabel usando MaTeX.";

ExportPublicationReady::usage = "";


Begin["`Private`"]

(* --- Constantes de Tama\[NTilde]o --- *)
ColumnWidth = 72 * 3.4; 
FullWidth = 72 * 7.0;
DefaultFontSize = 10.25;
$PreviewScale = 2.;

(* --- Configuraci\[OAcute]n Global de MaTeX --- *)
SetOptions[MaTeX, 
  "Preamble" -> {
    "\\usepackage{amsmath}", 
    "\\usepackage{physics}", 
    "\\usepackage{amsfonts}"
  },
  FontSize -> DefaultFontSize * $PreviewScale
];

(* --- Estilo Maestro --- *)
(* Definimos un estilo base que cumple con las reglas de la APS *)
APSPlotStyle := {
  BaseStyle -> {FontFamily -> "Latin Modern Math", FontSize -> DefaultFontSize * $PreviewScale},
  Frame -> True,
  FrameStyle -> Directive[Black, AbsoluteThickness[0.4]],
  GridLines -> None,
  ImageSize -> ColumnWidth * $PreviewScale, (* Por defecto una columna *)
  LabelStyle -> {Black}
};


(* --- Funci\[OAcute]n de Exportaci\[OAcute]n --- *)
Options[ExportToPDFAndCrop] = Options[Export];

ExportToPDFAndCrop[plot_, fileName_String, opts : OptionsPattern[]] := 
 Module[{fullPath, result},
  
  fullPath = If[StringEndsQ[fileName, ".pdf", IgnoreCase -> True], 
    fileName, 
    fileName <> ".pdf"
  ];
  
  (* Exportamos usando filtros para las opciones *)
  Export[fullPath, plot, FilterRules[{opts}, Options[Export]]];
  
  (* Ejecuci\[OAcute]n de pdfcrop *)
  result = RunProcess[{"pdfcrop", fullPath, fullPath}];
  
  If[result["ExitCode"] != 0,
    Message[ExportToPDFAndCrop::error, result["StandardError"]];
    $Failed,
    fullPath
  ]
];

ExportToPDFAndCrop::error = "Error en pdfcrop: `1`";


NiceFigure[plot_, opts : OptionsPattern[]] := Module[{processedOpts},
  
  (* 
    1. Tomamos todas las opciones dadas por el usuario: {opts}
    2. Buscamos la regla FrameLabel -> valores
    3. Si la encontramos, buscamos cualquier String dentro de sus valores
       y le aplicamos MaTeX, pas\[AAcute]ndole tambi\[EAcute]n el tama\[NTilde]o de fuente escalado.
  *)
  processedOpts = {opts} /. (FrameLabel -> val_) :> (
    FrameLabel -> (val /. s_String :> MaTeX[s, FontSize -> DefaultFontSize * $PreviewScale])
  );
  
  (* Combinamos el plot original, el estilo maestro y las opciones procesadas *)
  Show[plot, processedOpts, APSPlotStyle]
]


(* Le decimos a Mathematica que no eval\[UAcute]e el c\[OAcute]digo del gr\[AAcute]fico inmediatamente *)
SetAttributes[ExportPublicationReady, HoldFirst];

Options[ExportPublicationReady] = Options[Export];

ExportPublicationReady[plotCode_, fileName_String, opts : OptionsPattern[]] := 
  Block[{$PreviewScale = 1.},
    (* Dentro de este Block, toda funci\[OAcute]n que use $PreviewScale usar\[AAcute] 1 *)
    ExportToPDFAndCrop[plotCode, fileName, opts]
  ]


End[]
EndPackage[]
