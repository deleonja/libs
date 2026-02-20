(* ::Package:: *)

(*  If ForScience paclet not installed, install it. 
	See https://github.com/MMA-ForScience/ForScience *)
If[Length[PacletFind["ForScience"]] == 0, 
	PacletInstall[FileNameJoin[{DirectoryName[$InputFileName], "ForScience-0.88.45.paclet"}]]
];


(* ::Section:: *)
(*Begin package*)


BeginPackage["QMB`"];


(* For nice formatting of usage messages, see https://github.com/MMA-ForScience/ForScience *)
<<ForScience`;


(* ::Section:: *)
(*Notas*)


(* ::Text:: *)
(*Hay cosas en IAA_model.nb que 1) hay que migrar para ac\[AAcute] y 2) que hay que revisar si deber\[IAcute]a de poner ac\[AAcute]*)


(* ::Text:: *)
(*Hay cosas de Heisenberg meets fuzzy que tambi\[EAcute]n tengo que pasar para ac\[AAcute]*)


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsection::Closed:: *)
(*General quantum mechanics*)


(* All usage messages are evaluated quietly as FormatUsage[] requires FrontEnd. Therefore, if 
   QMB.wl is loaded in a .wls no error about FrontEnd pops up. *)


(*01234567890123456789012345678901234567890123456789012345678901234567890123456789*)
Quiet[
DensityMatrix::usage = FormatUsage[
"DensityMatrix[\[Psi]] returns the density matrix of state vector ```\[Psi]```."
];, {FrontEndObject::notavail, First::normal}];


Quiet[
Pauli::usage = FormatUsage[
	"Pauli[i] returns the ```i```-th Pauli matrix.\n"<>
	"Pauli[{i_1,...,i_N}] returns the ```N```-qubit Pauli string \
	'''Pauli[```i_1```]''' \[CircleTimes] '''...''' \[CircleTimes] '''Pauli[```i_N```]'''."];
, {FrontEndObject::notavail, First::normal}];


MatrixPartialTrace::usage = "MatrixPartialTrace[mat, n, d] calculates the partial trace of mat over the nth subspace, where all subspaces have dimension d.
MatrixPartialTrace[mat, n, {\!\(\*SubscriptBox[\(d\), \(1\)]\),\!\(\*SubscriptBox[\(d\), \(1\)]\),\[Ellipsis]}] calculates the partial trace of matrix mat over the nth subspace, where mat is assumed to lie in a space constructed as a tensor product of subspaces with dimensions {d1,d2,\[Ellipsis]}.";


VectorFromKetInComputationalBasis::usage = "VectorFromKetInComputationalBasis[ket] returns the matrix representation of ket.";


KetInComputationalBasisFromVector::usage = "KetInComputationalBasisFromVector[vector] returns the ket representation in computational basis of vector.";


Quiet[
RandomQubitState::usage = FormatUsage["RandomQubitState[] returns a Haar random qubit state."]
, {FrontEndObject::notavail, First::normal}];


Quiet[
RandomChainProductState::usage = FormatUsage[
"RandomChainProductState[L] returns a random ```L```-qubit product state."
], {FrontEndObject::notavail, First::normal}];


Quiet[
Dyad::usage = FormatUsage[
"Dyad[\[Psi]] returns ```|\[Psi]\[RightAngleBracket]\[LeftAngleBracket]\[Psi]|```.
Dyad[\[Psi],\[Phi]] returns ```|\[Psi]\[RightAngleBracket]\[LeftAngleBracket]\[Phi]|```."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
MatrixCommutator::usage = FormatUsage["MatrixCommutator[A,B] returns AB - BA."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
CommutationQ::usage = FormatUsage["CommutationQ[A,B] yields True if ```A``` and ```B``` commute, and False otherwise."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
MutuallyCommutingSetQ::usage= FormatUsage[
"MutuallyCommutingSetQ[{A,B,...}] yields True if all matrices ```{A,B,...}``` mutually commute, and False otherwise."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
Braket::usage = FormatUsage[
"Braket[\[Psi],\[Phi]] returns \[LeftAngleBracket]```\[Psi]```|```\[Phi]```\[RightAngleBracket]."
];
, {FrontEndObject::notavail, First::normal}];


FixCkForStateEvoultion::usage = "FixCkForStateEvoultion[\!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), { \!\(\*TemplateBox[{SubscriptBox[\"E\", \"k\"]},\n\"Ket\"]\) }] fixes \!\(\*SubscriptBox[\(c\), \(k\)]\) = \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"E\", \"k\"], \" \"}], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\n\"BraKet\"]\) for StateEvolution[]";


StateEvolution::usage = "StateEvolution[t, \!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), {E_i}, {\!\(\*TemplateBox[{\"E_i\"},\n\"Ket\"]\)} ] returns \!\(\*TemplateBox[{RowBox[{\"\[Psi]\", RowBox[{\"(\", \"t\", \")\"}]}]},\n\"Ket\"]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(\(\\ \)\(i\)\)]\) \!\(\*SuperscriptBox[\(\[ExponentialE]\), \(\(-\[ImaginaryI]\)\\  \*SubscriptBox[\(E\), \(i\)]\\  t\)]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\n\"BraKet\"]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"]},\n\"Ket\"]\).
StateEvolution[t, {\!\(\*SubscriptBox[\(E\), \(k\)]\)}] calculates \!\(\*TemplateBox[{RowBox[{\"\[Psi]\", RowBox[{\"(\", \"t\", \")\"}]}]},\"Ket\"]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(\(\\\\\)\(i\)\)]\) \!\(\*SuperscriptBox[\(\[ExponentialE]\), \(\(-\[ImaginaryI]\)\\\\\*SubscriptBox[\(E\), \(i\)]\\\\t\)]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\"BraKet\"]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"]},\"Ket\"]\) having fixed the \!\(\*SubscriptBox[\(c\), \(k\)]\)'s with FixCkForStateEvoultion[\!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), { \!\(\*TemplateBox[{SubscriptBox[\"E\", \"k\"]},\n\"Ket\"]\) }].";


Quiet[
BlochVector::usage = FormatUsage["BlochVector[\[Rho]] returns the Bloch vector of a single-qubit density matrix \[Rho]."];
, {FrontEndObject::notavail, First::normal}];


KroneckerVectorProduct::usage = "KroneckerVectorProduct[a,b] calculates \!\(\*TemplateBox[{\"a\"},\n\"Ket\"]\)\[CircleTimes]\!\(\*TemplateBox[{\"b\"},\n\"Ket\"]\).";


Quiet[
Purity::usage = FormatUsage[
"Purity[\[Rho]] calculates the purity of ```\[Rho]```."
];
, {FrontEndObject::notavail, First::normal}];


Quiet[
Concurrence::usage = FormatUsage[
"Concurrence[\[Rho]] returns the two-qubit concurrence of density matrix ```\[Rho]```."
];
, {FrontEndObject::notavail, First::normal}];


Quiet[
Qubit::usage = FormatUsage[
"Qubit[\[Theta],\[Phi]] returns the state cos(```\[Theta]```/2)|0\[RightAngleBracket] + \[ExponentialE]^\[Phi] sin(```\[Theta]```/2)|1\[RightAngleBracket]"
];
, {FrontEndObject::notavail, First::normal}];


Quiet[
SU2Rotation::usage = FormatUsage[
	"SU2Rotation[{\[Theta]_a,\[Phi]_a},\[Theta]_R] the SU(2) matrix rotation with axis ```(\[Theta]_a,\[Phi]_a)``` \
	and angle rotation ```\[Theta]_R```.
	'''SU2Rotation[{```x,y,z```},```\[Theta]_R```]''' the SU(2) matrix rotation with axis ```(x,y,z)``` \
	and angle rotation ```\[Theta]_R```."
];
, {FrontEndObject::notavail, First::normal}];


(* ::Text:: *)
(*Agregadas por Miguel*)


coherentstate::usage = "coherentstate[state,L] Generates a spin coherent state of L spins given a general single qubit state";


Quiet[
SU2Rotation::usage = FormatUsage[
	"SU2Rotation[{\[Theta]_a,\[Phi]_a},\[Theta]_R] devuelve la matriz de rotaci\[OAcute]n SU(2) \
	correspondiente a una rotaci\[OAcute]n de \[AAcute]ngulo ```\[Theta]_R``` alrededor del eje definido por las coordenadas \
	esf\[EAcute]ricas ```(\[Theta]_a,\[Phi]_a)```.\n" <>
	"SU2Rotation[{nx,ny,nz},\[Theta]_R] usa un vector cartesiano como eje."
];
, {FrontEndObject::notavail, First::normal}];


(* ::Subsection::Closed:: *)
(*Quantum chaos*)


MeanLevelSpacingRatio::usage = FormatUsage[
	"MeanLevelSpacingRatio[eigenvalues] returns \[LeftAngleBracket]r_n\[RightAngleBracket]."
];


Quiet[
	IPR::usage = 
		FormatUsage["IPR[\[Psi]] computes the Inverse Participation Ratio of  ```\[Psi]``` in computational basis."];
, {FrontEndObject::notavail, First::normal}];


kthOrderSpacings::usage = FormatUsage[
	"kthOrderSpacings[spectrum,k] returns the ```k```-th order level spacing of ```spectrum```."
];


SpacingRatios::usage = FormatUsage[
	"SpacingRatios[spectrum,k] returns the ```k```-th order level spacing ratios of ```spectrum```."
];


Unfold::usage = FormatUsage[
  "Unfold[spectrum] returns an Association containing the unfolded ```spectrum``` and \
smoothed density functions calculated via Kernel Density Estimation (KDE).
Unfold[spectrum, opts] allows specifying options for the kernel distribution.

**Output Association Keys**
* \"UnfoldedLevels\": List of unfolded eigenvalues (mean spacing = 1).
* \"SmoothCDF\": Pure function representing the cumulative mean level number N(E).
* \"SmoothPDF\": Pure function representing the mean level density \[Rho](E).
* \"Bandwidth\": The bandwidth parameter used in the KDE.
* \"OriginalLevels\": The sorted input spectrum.

**Options**
* \"Bandwidth\": (Default: Automatic) Controls the smoothing scale. Set a Real number to manually tune the separation between secular variation and fluctuations.
* \"Kernel\": (Default: \"Gaussian\") Specifies the kernel function type (e.g., \"Epanechnikov\", \"Rectangular\")."
];


ComplexSpacingRatios::usage =
"ComplexSpacingRatios[eigs_List] computes the complex spacing ratios (CSR) \
for a list of complex eigenvalues. \
The CSR is defined for each eigenvalue \[Lambda]_k as: \
z_k = (\[Lambda]_NN - \[Lambda]_k) / (\[Lambda]_NNN - \[Lambda]_k), \
where NN and NNN are the nearest and next-to-nearest neighbors in the \
complex plane."


SpectralFormFactor::usage = FormatUsage[
"SpectralFormFactor[spectrum, t] computes the Spectral Form Factor K(t) \
for a given energy ```spectrum``` at time ```t```.\nSpectralFormFactor[spectrum, tList] \
computes K(t) vectorized over a list of times."];


TimeAveragedSFF::usage = FormatUsage[
"TimeAveragedSFF[tList, sffValues, windowSize] returns \
the time-averaged SFF as a list of coordinate pairs {```tAveraged```, \
```sffAveraged```} using a moving window to smooth out rapid fluctuations."];


NumberVariance::usage = FormatUsage[
"NumberVariance[unfoldedSpectrum, L, numSamples] computes the spectral number variance \
\[CapitalSigma]^2(L) of an unfolded spectrum by sampling ```numSamples``` windows of length ```L```."
];

AnalyticalNumberVariancePoisson::usage = FormatUsage[
"AnalyticalNumberVariancePoisson[L] returns the analytical number \
variance for a Poissonian spectrum: L."
];

AnalyticalNumberVarianceGOE::usage = FormatUsage[
"AnalyticalNumberVarianceGOE[L] returns the asymptotic analytical \
number variance for the Gaussian Orthogonal Ensemble."
];

AnalyticalNumberVarianceGUE::usage = FormatUsage[
"AnalyticalNumberVarianceGUE[L] returns the asymptotic analytical \
number variance for the Gaussian Unitary Ensemble."
];


(* ::Subsection::Closed:: *)
(*RMT*)


RatiosDistribution::usage = FormatUsage[
	"RatiosDistribution[r,\[Beta]] represents the probability distribution of level spacing \
	ratios P_\[Beta](r)."
];


RatiosDistributionPoisson::usage = FormatUsage[
	"RatiosDistributionPoisson[r,k] represents the probability distribution of level spacing \
	ratios P(r) of a Poissonian spectrum."
];


NormalizedSpacingRatios::usage = FormatUsage[
"NormalizedSpacingRatios[spectrum, k] returns the list of spacing ratios \
mapped to the interval (0, 1] via min(r, 1/r)."
];


Quiet[
KLDivergence::usage = FormatUsage[
"KLDivergence[P, Q] computes the Kullback-Leibler divergence D_{*KL*}(P||Q) = \[CapitalSigma]_i P(i) log(P(i)/Q(i)).
KLDivergence[P, Q, b] computes the KL divergence using logarithms to base ```b```."];
, {FrontEndObject::notavail, First::normal}];

Quiet[
EmpiricalKLDivergence::usage = FormatUsage[
"EmpiricalKLDivergence[data, pdf, bins] computes the KL divergence between the histogram of ```data``` and a theoretical ```pdf``` function using specified ```bins```.
EmpiricalKLDivergence[data, pdf] automatically determines bins (Freedman-Diaconis rule).

**Example for RMT:**
EmpiricalKLDivergence[ratios, RatiosDistribution[#, 1]&, {0, 5, 0.1}] measures distance to GOE."];
, {FrontEndObject::notavail, First::normal}];


(* ::Subsubsection::Closed:: *)
(*Ginibre matrices*)


GenerateGinibreMatrix::usage =
"GenerateGinibreMatrix[n, OptionsPattern[]] generates a non-Hermitian \
random matrix of dimension n from a Ginibre ensemble.

Options:
  Ensemble -> \"Unitary\" (default), \"Orthogonal\", or \"Symplectic\".

- \"Unitary\" (GinUE, \[Beta]=2): Returns an n x n complex matrix.
- \"Orthogonal\" (GinOE, \[Beta]=1): Returns an n x n real matrix.
- \"Symplectic\" (GinSE, \[Beta]=4): Returns a 2n x 2n complex matrix \
(representing an n x n quaternion matrix)."


GenerateGinibreMatrix::invalidEnsemble = "Ensemble type `1` is not recognized. \
Use \"Unitary\", \"Orthogonal\", or \"Symplectic\".";


(* ::Subsection::Closed:: *)
(*Quantum channels*)


Reshuffle::usage = "Reshuffle[m] applies the reshuffle transformation to the matrix m with dimension \!\(\*SuperscriptBox[\(d\), \(2\)]\)\[Times]\!\(\*SuperscriptBox[\(d\), \(2\)]\).
Reshuffle[A,m,n] reshuffles matrix A, where dim(A) = mn.";


Quiet[
	SuperoperatorFromU::usage = FormatUsage["DensityMatrix[\[Psi]] returns the density matrix of state vector ```\[Psi]```."];
, {FrontEndObject::notavail, First::normal}];


(* ::Subsection::Closed:: *)
(*Bose-Hubbard*)


BoseHubbardHamiltonian::usage = FormatUsage[
	"BoseHubbardHamiltonian[n,L,J,U,opts] returns the Bose-Hubbard Hamiltonian for \
	```n``` bosons and ```L``` sites with hopping parameter ```J``` and interaction \
	parameter ```U```. Options: SymmetricSubspace."
];


SymmetricSubspace::usage = FormatUsage[
	"SymmetricSubspace is an option for '''BoseHubbardHamiltonian'''. Valid values are \
	'''\"All\"''', '''\"EvenParity\"''', and '''\"OddParity\"'''."
];


KineticTermBoseHubbardHamiltonian::usage = FormatUsage[
	"KineticTermBoseHubbardHamiltonian[basis] returns the kinetic term of the BH Hamiltonian \
	with ```basis``` a list with Fock basis elements.\n"<> 
	"KineticTermBoseHubbardHamiltonian[basis,SymmetricSubspace] returns the kinetic term of the BH \
	Hamiltonian in a symmetric subspace with ```basis``` a list with Fock basis elements. \
	Option SymmetricSubspace takes the values \"All\" | \"EvenParity\" | \"OddParity\"."
];


PotentialTermBoseHubbardHamiltonian::usage = FormatUsage[
  "PotentialTermBoseHubbardHamiltonian[n, L, SymmetricSubspace] returns the \
potential term of the Bose-Hubbard Hamiltonian for ```n``` bosons and ```L``` \
sites.
PotentialTermBoseHubbardHamiltonian[basis] returns the potential term of the \
Bose-Hubbard Hamiltonian given the ```basis``` of the Fock space of the bosonic \
system.\n\n"
  (*"**Notes**\n" <>
  "- If you want the potential term in a parity symmetry sector, ```basis``` \
should be a list containing only the representative Fock states of that \
symmetric subspace.\n\n" <>
  "**Examples of usage**\n" <>
  "- '''PotentialTermBoseHubbardHamiltonian[{{3,0,0},{2,1,0},{2,0,1},{1,2,0}}]''' \
returns the matrix in the odd parity sector of a system with 3 bosons and 3 \
sites."
*)];



BosonEscapeKrausOperators::usage = "BosonEscapeKrausOperators[N, L]: bosons escape to nearest neighbouring sites. N: bosons, L: site.";


BosonEscapeKrausOperators2::usage = "sdfa";


HilbertSpaceDim::replaced = "Function `1` has been replaced by `2`.";


BoseHubbardHilbertSpaceDimension::usage = FormatUsage["BoseHubbardHilbertSpaceDimension[n,L] returns the dimension"<> 
"of Hilbert space of a Bose Hubbard system of N bosons and L sites."];


FockBasis::usage = "FockBasis[N, M] returns the lexicographical-sorted Fock basis of N bosons in M sites.";


SortFockBasis::usage = "SortFockBasis[fockBasis] returns fockBasis in ascending-order according to the tag of Fock states.";


Tag::usage = "Tag[ { \!\(\*SubscriptBox[\(k\), \(1\)]\),\!\(\*SubscriptBox[\(k\), \(2\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(k\), \(L\)]\) } ] returns the tag \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\)\!\(\*SqrtBox[\(100  i + 3\)]\)\!\(\*SubscriptBox[\(k\), \(i\)]\) of Fock state \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"1\"], \",\", SubscriptBox[\"k\", \"2\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"L\"]}]},\n\"Ket\"]\).";


FockBasisStateAsColumnVector::usage = "FockBasisStateAsColumnVector[FockState, N, L] returns matrix representation of fockState={\!\(\*SubscriptBox[\(i\), \(1\)]\),\!\(\*SubscriptBox[\(i\), \(2\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(L\)]\)}. N: bosons, L: sites.";


FockBasisIndex::usage = "FockBasisIndex[fockState, sortedTagsFockBasis] returns the position of fockState in the tag-sorted Fock basis with tags sortedTagsFockBasis.";


RenyiEntropy::usage = "RenyiEntropy[\[Alpha], \[Rho]] computes the \[Alpha]-th order Renyi entropy of density matrix \[Rho].";


BosonicPartialTrace::usage = "BosonicPartialTrace[\[Rho]] calculates the partial trace of \[Rho]. Requires initialization.";


InitializationBosonicPartialTrace::usage = "InitializationBosonicPartialTrace[{\!\(\*SubscriptBox[\(i\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(k\)]\)}, N, L] initializes variables for BosonicPartialTrace[] to calculate the reduced density matrix of sites {\!\(\*SubscriptBox[\(i\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(k\)]\)}.";


(* ::Subsubsection::Closed:: *)
(*Fuzzy measurements in bosonic systems*)


InitializeVariables::usage = "InitializeVariables[n, L, boundaries, FMmodel] sets up the necessary variables for correct running of FuzzyMeasurement[\[Psi], \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\)]; boundaries: 'open' or 'closed'; FMmodel: '#NN'.";


FuzzyMeasurement::usage = "FuzzyMeasurement[\[Psi], \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\)] gives \[ScriptCapitalF](\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\)) = (1 - \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\))\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\) + \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\) \!\(\*UnderscriptBox[\(\[Sum]\), \(i\)]\) \!\(\*SubscriptBox[\(S\), \(i\)]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\)\!\(\*SubsuperscriptBox[\(S\), \(i\), \(\[Dagger]\)]\), where \!\(\*SubscriptBox[\(S\), \(i\)]\) must be initizalized runnning InitializeVariables[n, L, boundaries, FMmodel].";


(* ::Subsection::Closed:: *)
(*Spin chains*)


(* ::Subsubsection::Closed:: *)
(*Symmetries*)


SpinParityEigenvectors::usage = "SpinParityEigenvectors[L] gives a list of {even, odd} eigenvectors of the L-spin system parity operator P; P\!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"1\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"L\"]}]},\n\"Ket\"]\) = \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"L\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"1\"]}]},\n\"Ket\"]\), \!\(\*SubscriptBox[\(k\), \(i\)]\)=0,1.";


TranslationEigenvectorRepresentatives::usage = FormatUsage[
  "TranslationEigenvectorRepresentatives[L] returns a list of sublists, each containing:\
  the decimal representation of a bit-string representative eigenvector,\
  its pseudomomentum k, and the length of its translation orbit,\
  for a system of ```L``` qubits."
];


BlockDiagonalize::usage = FormatUsage[
	"BlockDiagonalize[matrix,opts] returns ```matrix``` in block-diagonal form. Only \
	option is '''Symmetry'''."
];


Symmetry::usage = FormatUsage[
	"Symmetry is an option for '''BlockDiagonalize''' to specify the symmetry according to which \
	```matrix``` in '''BlockDiagonalize'''[```matrix```] is block diagonalized. It takes the values: \
	\"Translation\". Default option is \"Translation\"."
];


(* ::Subsubsection::Closed:: *)
(*Hamiltonians*)


IsingHamiltonian::usage = FormatUsage[
	"IsingHamiltonian[h_x,h_z,J,L,opts] returns the Hamiltonian \ 
	H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - ```J``` \[Sum]_{*i=1*}^{*L-1*} \[Sigma]^z_i \[Sigma]^z_{*i+1*} \
	with boundary conditions specified by option BoundaryConditions (default is \"Open\")."
];


BoundaryConditions::usage = FormatUsage[
	"BoundaryConditions is an option for '''IsingHamiltonian''' to specify the boundary conditions. It \
	takes the values \"Open\" or \"Periodic\". Default option is \"Open\"."
];


IsingNNOpenHamiltonian::replaced = "Function `1` has been replaced by `2`.";


(*Quiet[
IsingNNOpenHamiltonian::usage = FormatUsage["IsingNNOpenHamiltonian[h_x,h_z,J,L] returns the Hamiltonian H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - ```J```\[Sum]_{*i=1*}^{*L-1*} \[Sigma]^z_i \[Sigma]^z_{*i+1*}.
IsingNNOpenHamiltonian[h_x,h_z,{J_1,...,J_L},L] returns the Hamiltonian H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - \[Sum]_{*i=1*}^{*L-1*} ```J_i``` \[Sigma]^z_i \[Sigma]^z_{*i+1*}."];
, {FrontEndObject::notavail, First::normal}];*)


IsingNNClosedHamiltonian::replaced = "Function `1` has been replaced by `2`.";


(*IsingNNClosedHamiltonian::usage = "IsingNNClosedHamiltonian[\!\(\*
StyleBox[SubscriptBox[\"h\", \"x\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"h\", \"z\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"J\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)] returns the Hamiltonian \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\)(\!\(\*SubscriptBox[\(h\), \(x\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(x\)]\) + \!\(\*SubscriptBox[\(h\), \(z\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\)) + \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), L]\) \!\(\*SubscriptBox[\(J\), \(i\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\)\!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(z\)]\) with \!\(\*SubscriptBox[\(\[Sigma]\), \(L + 1\)]\) = \!\(\*SubscriptBox[\(\[Sigma]\), \(1\)]\).";*)


ClosedXXZHamiltonian::usage = "ClosedXXZHamiltonian[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[CapitalDelta]\",\nFontSlant->\"Italic\"]\)] returns the closed XXZ 1/2-spin chain as in appendix A.1 of Quantum 8, 1510 (2024).";


OpenXXZHamiltonian::usage= "OpenXXZHamiltonian[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[CapitalDelta]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"h1\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"h2\",\nFontSlant->\"Italic\"]\)] returns the open XXZ 1/2-spin chain as in appendix A.2 of Quantum 8, 1510 (2024).";


Quiet[
LeaSpinChainHamiltonian::usage = FormatUsage["LeaSpinChainHamiltonian[J_{*xy*},J_z,\[Omega],\[Epsilon]_d,L,d] returns the spin-1/2 chain H = \[Sum]_{*i=1*}^{*L-1*} ```J_{*xy*}```(S^x_i S^x_{*i+1*} + S^y_i S^y_{*i+1*}) + ```J_z```S^z_i S^z_{*i+1*} + \[Sum]_{*i=1*}^{*L*} ```\[Omega]``` S^z_i + \[Epsilon]_d S^z_d. [Eq. (1) in Am. J. Phys. 80, 246\[Dash]251 (2012)]."];
, {FrontEndObject::notavail, First::normal}];


Quiet[
XXZOpenHamiltonian::usage = FormatUsage["XXZOpenHamiltonian[J_{*xy*},J_z,\[Omega],\[Epsilon]_d,L,d] returns the spin-1/2 chain \n H = \[Sum]_{*i=1*}^{*L-1*} ```J_{*xy*}```(S^x_i S^x_{*i+1*} + S^y_i S^y_{*i+1*}) + ```J_z```S^z_i S^z_{*i+1*} + \[Sum]_{*i=1*}^{*L*} ```\[Omega]``` S^z_i + \[Epsilon]_d S^z_d. \n [Eq. (1) in Am. J. Phys. 80, 246\[Dash]251 (2012)]."];
, {FrontEndObject::notavail, First::normal}];


HeisenbergXXXwNoise::usage="HeisenbergXXXwNoise[hz,L] returns the Heisenberg XXX spin 1/2 chain with noise: \!\(\*FormBox[\(H\\\  = \\\ \*FractionBox[\(1\), \(4\)]\\\ \(\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(L - 1\)]\\\ \((\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(x\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(x\)]\\\  + \\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(y\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(y\)]\\\  + \\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(z\)])\)\)\\\  + \\\ \*FractionBox[\(1\), \(2\)]\\\ \(\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\*SubsuperscriptBox[\(h\), \(i\), \(z\)]\\\ \*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\\\ \(\((open\\\ boundaries)\)\(.\)\)\)\),
TraditionalForm]\)";


QuantumGameOfLifeHamiltonian::usage = FormatUsage[
  "QuantumGameOfLifeHamiltonian[L, opts] returns the Hamiltonian for the Quantum Game of Life (QGL) \
  model for a spin-1/2 chain of length ```L``` (arXiv:2510.16570):
  H = \[Sum] \!\(\[Sigma]\_i\^x\) (\!\(\[ScriptN]\_i\^2\) + \!\(\[ScriptN]\_i\^3\)).
  
  Options:
  Momentum -> \"All\" (default) | Integer k (0 <= k < L).
  Parity -> \"None\" (default) | 1 | -1.
  
  Note: Parity projection (Inversion Symmetry) is only valid for Momentum sectors k=0 or k=L/2."
];


Momentum::usage = "Momentum is an option for QuantumGameOfLifeHamiltonian. \
It specifies the translation symmetry sector (integer k).";

Parity::usage = "Parity is an option for QuantumGameOfLifeHamiltonian. \
It specifies the inversion symmetry sector (+1 or -1).";


(* ::Subsection::Closed:: *)
(*Spin Hamiltonians*)


TwoSpinHamiltonian::usage = FormatUsage[
	"TwoSpinHamiltonian[s_1, s_2, \[Epsilon]_1, \[Epsilon]_2, \[Theta], g] returns the Hamiltonian \
	matrix for two spins of magnitude ```s_1``` and ```s_2```.
	The Hamiltonian is given by:
	H = ```\[Epsilon]_1``` S_1^z + ```\[Epsilon]_2``` (Cos[```\[Theta]```]S_2^z + Sin[```\[Theta]```]S_2^x) + \
	(```g```/Sqrt(s_1 s_2)) S_1 \[CenterDot] S_2."
];


(* ::Subsection::Closed:: *)
(*Light-matter systems*)


ARQMHamiltonian::usage = 
FormatUsage[
	"ARQMHamiltonian[\[Omega], \[CapitalDelta], g, \[Epsilon], \[Xi], nMax] calculates \
	the Hamiltonian matrix for the Anisotropically quantum Rabi \
	model (AQRM) in a truncated Fock space of dimension nMax + 1.

	The Hamiltonian is:
	H = \[Omega] a\[Dagger]a + \[CapitalDelta]*Z + g/2*[(1+\[Xi]) (a+a\[Dagger])*X \
	+ (1-\[Xi]) (a-a\[Dagger]) \[ImaginaryI]*Y] + \[Epsilon] X	

	Parameters are:
	\[Omega]: bosonic mode frequency (photons).
	\[CapitalDelta]: qubit energy splitting.
	g: qubit-boson coupling strength.
	\[Epsilon]: asymmetry term strength (qubit driving).
	\[Xi]: anisotropy parameter (ranging from -1 to 1).
	nMax: maximum photon number (n=0 to nMax)."
];


(* ::Subsection::Closed:: *)
(*Fuzzy measurement channels*)


SwapMatrix::usage = "SwapMatrix[targetSite, wrongSite, N, L] returns the swap matrix that exchanges site targetSite with wrongSite for a system of N bosons and L sites.";


FuzzyMeasurementChannel::usage = 
"FuzzyMeasurementChannel[\[Rho], p, PermMatrices] returns (1-p) \[Rho] + \!\(\*FractionBox[\(p\), \(N - 1\)]\) \!\(\*SubscriptBox[\(\[Sum]\), \(i\)]\)\!\(\*SubscriptBox[\(S\), \(i, i + 1\)]\)\!\(\*SubscriptBox[\(\[Rho]S\), \(i, i + 1\)]\).
FuzzyMeasurementChannel[\[Rho], {\!\(\*SubscriptBox[\(p\), \(totalError\)]\), \!\(\*SubscriptBox[\(p\), \(NN\)]\), \!\(\*SubscriptBox[\(p\), \(SNN\)]\)}, {{\!\(\*SubscriptBox[\(S\), \(i, i + 1\)]\)}, {\!\(\*SubscriptBox[\(S\), \(i, i + 2\)]\)}}] returns \[ScriptCapitalE](\[Rho])=(1-\!\(\*SubscriptBox[\(p\), \(totalError\)]\))\[Rho] + \!\(\*SubscriptBox[\(p\), \(totalError\)]\)(\!\(\*FractionBox[SubscriptBox[\(p\), \(NN\)], \(L - 1\)]\)) \!\(\*SuperscriptBox[SubscriptBox[\(\[Sum]\), \(i = 1\)], \(L - 1\)]\) \!\(\*SubscriptBox[\(S\), \(i, i + 1\)]\) \[Rho] \!\(\*SubscriptBox[\(S\), \(i, i + 1\)]\) + \!\(\*SubscriptBox[\(p\), \(totalError\)]\)(\!\(\*FractionBox[SubscriptBox[\(p\), \(SNN\)], \(L - 2\)]\)) \!\(\*SuperscriptBox[SubscriptBox[\(\[Sum]\), \(i = 1\)], \(L - 2\)]\) \!\(\*SubscriptBox[\(S\), \(i, i + 2\)]\) \[Rho] \!\(\*SubscriptBox[\(S\), \(i, i + 2\)]\).";


(* ::Section:: *)
(*Beginning of Package*)


Begin["`Private`"];


(* ::Section:: *)
(*Routine definitions*)


(*no poner los nombres de funciones p\[UAcute]blicas porque se joden la definici\[OAcute]n de uso*)
ClearAll[SigmaPlusSigmaMinus, SigmaMinusSigmaPlus];


(* ::Subsection::Closed:: *)
(*General quantum mechanics*)


DensityMatrix[\[Psi]_] := Outer[Times, \[Psi], Conjugate[\[Psi]]]


Pauli[0]=Pauli[{0}]=SparseArray[{{1,0}, {0,1}}]; 
Pauli[1]=Pauli[{1}]=SparseArray[{{0,1}, {1,0}}]; 
Pauli[2]=Pauli[{2}]=SparseArray[{{0,-I},{I,0}}]; 
Pauli[3]=Pauli[{3}]=SparseArray[{{1,0}, {0,-1}}];
Pauli[Indices_List] := KroneckerProduct @@ (Pauli /@ Indices)


MatrixPartialTrace=ResourceFunction["MatrixPartialTrace"];


VectorFromKetInComputationalBasis[ket_]:=Normal[SparseArray[FromDigits[ket,2]+1->1,Power[2,Length[ket]]]]


KetInComputationalBasisFromVector[vector_]:=IntegerDigits[Position[vector,1][[1,1]]-1,2,Log[2,Length[vector]]]


RandomQubitState[] := 
Module[{x,y,z,\[Theta],\[Phi]},
	{x,y,z} = RandomPoint[Sphere[]];
	{\[Theta],\[Phi]}={ArcCos[z],Sign[y]ArcCos[x/Sqrt[x^2+y^2]]};
	{Cos[\[Theta]/2],Exp[I \[Phi]]Sin[\[Theta]/2]}
]


RandomChainProductState[0] := {1}
RandomChainProductState[1] := RandomQubitState[]
RandomChainProductState[L_] := Flatten[KroneckerProduct@@Table[RandomQubitState[],L]]


Dyad[a_] := Outer[Times,a,Conjugate[a]]
Dyad[a_,b_] := Outer[Times,a,Conjugate[b]]


MatrixCommutator[A_,B_]:=A . B-B . A


ZeroMatrix[d_]:=ConstantArray[0,{d,d}]


CommutationQ[A_,B_]:=MatrixCommutator[A,B]==ZeroMatrix[Length[A]]


MutuallyCommutingSetQ[ListOfMatrices_]:=Module[{SetLength=Length[ListOfMatrices]},
AllTrue[Table[CommutationQ@@ListOfMatrices[[{i,j}]],{i,SetLength-1},{j,i+1,SetLength}],TrueQ,2]
]


Braket[a_,b_] := Conjugate[a] . b


StateEvolution[t_,psi0_List,eigenvals_List,eigenvecs_List]:=
(*|\[Psi](t)\[RightAngleBracket] = Underscript[\[Sum], k] Subscript[c, k]\[ExponentialE]^(-Subscript[\[ImaginaryI]E, k]t)|Subscript[E, k]\[RightAngleBracket], Subscript[c, k]=\[LeftAngleBracket]Subscript[E, k]\[VerticalSeparator] Subscript[\[Psi], 0]\[RightAngleBracket]*)
	Module[{ck},
		ck = N[Chop[Conjugate[eigenvecs] . psi0]];
		N[Chop[Total[ ck * Exp[-I*eigenvals*N[t]] * eigenvecs]]]
	]


FixCkForStateEvoultion[\[Psi]0_, eigenvecs_] :=
	Module[{},
		ck = N[ Chop[ Conjugate[eigenvecs] . \[Psi]0 ] ];
		Heigenvecs = eigenvecs;
	]


StateEvolution[t_,eigenvals_List]:=
(*|\[Psi](t)\[RightAngleBracket] = Underscript[\[Sum], k] Subscript[c, k]\[ExponentialE]^(-Subscript[\[ImaginaryI]E, k]t)|Subscript[E, k]\[RightAngleBracket], Subscript[c, k]=\[LeftAngleBracket]Subscript[E, k]\[VerticalSeparator] Subscript[\[Psi], 0]\[RightAngleBracket]*)
	N[Chop[Total[ ck * Exp[-I*eigenvals*N[t]] * Heigenvecs]]]


BlochVector[\[Rho]_]:=Chop[Tr[Pauli[#] . \[Rho]]&/@Range[3]]


KroneckerVectorProduct[a_,b_]:=Flatten[KroneckerProduct[a,b]]


Purity[\[Rho]_]:=Tr[\[Rho] . \[Rho]] 


Concurrence[\[Rho]_] :=
Module[{\[Rho]tilde, R, \[Lambda]},
	\[Rho]tilde=# . Conjugate[\[Rho]] . #&[Pauli[{2,2}]];
	R=MatrixPower[# . \[Rho]tilde . #&[MatrixPower[\[Rho],1/2]],1/2];
	\[Lambda]=ReverseSort[Chop[Eigenvalues[R]]];
	Max[0,\[Lambda][[1]]-Total[\[Lambda][[2;;]]]]
]


Qubit[\[Theta]_,\[Phi]_] := {Cos[\[Theta]/2],Exp[I \[Phi]]Sin[\[Theta]/2]}


coherentstate[state_,L_]:=Flatten[KroneckerProduct@@Table[state,L]]


SU2Rotation[sphCoord_List?VectorQ,\[Theta]_]/; Length[sphCoord]==2 :=
	Module[{n={Sin[#1]Cos[#2], Sin[#1]Sin[#2], Cos[#1]} & @@ sphCoord},
		MatrixExp[-I * \[Theta]/2 * n . (Pauli /@ Range[3])]
	]

SU2Rotation[n_List?VectorQ,\[Theta]_]/; Length[n]==3 && Chop[Norm[n]-1]==0 :=
	MatrixExp[-I * \[Theta]/2 * n . (Pauli /@ Range[3])]


SU2Rotation[sphCoord_List, \[Theta]R_] /; Length[sphCoord] == 2 := 
    Module[{nx, ny, nz, c, s, phase},
        (* Conversi\[OAcute]n Esf\[EAcute]ricas -> Cartesianas (Radio = 1) *)
        {nx, ny, nz} = {Sin[sphCoord[[1]]] Cos[sphCoord[[2]]], Sin[sphCoord[[1]]] Sin[sphCoord[[2]]], Cos[sphCoord[[1]]]};
        
        (* Pre-c\[AAcute]lculo de t\[EAcute]rminos trigonom\[EAcute]tricos *)
        c = Cos[\[Theta]R / 2];
        s = Sin[\[Theta]R / 2];
        
        (* Construcci\[OAcute]n directa de la matriz (Evita MatrixExp para velocidad) *)
        (* n . sigma = {{nz, nx - i ny}, {nx + i ny, -nz}} *)
        {{c - I * s * nz, -I * s * (nx - I * ny)}, 
         {-I * s * (nx + I * ny), c + I * s * nz}}
    ]

(* Soporte para vector Cartesiano normalizado *)
SU2Rotation[n_List, \[Theta]R_] /; Length[n] == 3 := 
    Module[{c, s, nx, ny, nz},
        {nx, ny, nz} = Normalize[n]; (* Asegurar normalizaci\[OAcute]n *)
        c = Cos[\[Theta]R / 2];
        s = Sin[\[Theta]R / 2];
        
        {{c - I * s * nz, -I * s * (nx - I * ny)}, 
         {-I * s * (nx + I * ny), c + I * s * nz}}
    ]


(* ::Subsection:: *)
(*Quantum chaos*)


MeanLevelSpacingRatio[eigenvalues_]:=Mean[Min/@Transpose[{#,1/#}]&[Ratios[Differences[Sort[eigenvalues]]]]]


IPR[\[Psi]_] := Total[\[Psi]^4]


kthOrderSpacings[spectrum_, k_] := RotateLeft[#, k] - # &[Sort[spectrum]][[;; -(k+1)]]


SpacingRatios[spectrum_, k_]:=RotateLeft[#, k]/# &[kthOrderSpacings[spectrum, k]][[;; -(k+1)]]


(*Unfold*)

(* Options for Unfold allow tuning the kernel smoothness *)
Options[Unfold] = {
    "Bandwidth" -> Automatic, (* Can be set to a real number for manual control *)
    "Kernel" -> "Gaussian"
};

Unfold[spectrum_List, opts : OptionsPattern[]] := Module[
    {
        sortedSpectrum, 
        nLevels, 
        skd, 
        unfoldedLevels, 
        smoothCDFFunc, 
        smoothPDFFunc,
        bwParam,
        kernelType
    },

    (* 1. Validation and Preparation *)
    (* RMT requires sorted, numeric eigenvalues *)
    sortedSpectrum = Sort[Select[spectrum, NumericQ]];
    nLevels = Length[sortedSpectrum];

    If[nLevels < 3, 
        Message[Unfold::notEnoughLevels, nLevels]; 
        Return[$Failed]
    ];

    (* 2. Algorithmic Logic: Kernel Density Estimation (KDE) *)
    (* We model the secular (smooth) part of the density of states.
       The cumulative mapping \[CurlyEpsilon]_i = N_smooth(E_i) gives the unfolded spectrum. *)
    
    bwParam = OptionValue["Bandwidth"];
    kernelType = OptionValue["Kernel"];

    (* Construct the distribution object *)
    skd = SmoothKernelDistribution[
        sortedSpectrum, 
        bwParam, 
        kernelType
    ];

    (* 3. Compute Unfolded Levels *)
    (* Map physics eigenvalues to mean level spacing of 1 *)
    unfoldedLevels = nLevels * CDF[skd, sortedSpectrum];

    (* 4. Functional Construction with Scope Injection *)
    With[{distObject = skd, n = nLevels},
        smoothCDFFunc = Function[e, n * CDF[distObject, e]];
        smoothPDFFunc = Function[e, n * PDF[distObject, e]];
    ];

    (* 5. Return Structured Data *)
    <|
        "UnfoldedLevels" -> unfoldedLevels,
        "SmoothCDF" -> smoothCDFFunc,      (* Cumulative Mean Level Number *)
        "SmoothPDF" -> smoothPDFFunc,      (* Mean Level Density *)
        "Bandwidth" -> skd["Bandwidth"],
        "OriginalLevels" -> sortedSpectrum
    |>
];

(* Error Message definition *)
Unfold::notEnoughLevels = "Spectrum length `1` is too short for unfolding. At least 3 levels are required.";


ComplexSpacingRatios[eigs_List?VectorQ] := Module[
  {
    nnFunction,    (* The optimized neighbor-search function *)
    neighborsList  (* List to store the {k, NN, NNN} for each point *)
  },

  (* 1. Validate input *)
  If[!AllTrue[eigs, NumericQ] || Length[eigs] < 3,
    Print["ComplexSpacingRatios::Error: Input must be a list of at least 3 numeric values."];
    Return[$Failed];
  ];

  (* 2. Create the NearestFunction. This is the key optimization, built ONCE. *)
  nnFunction = Nearest[eigs];

  (* 3. Get the 3-neighbor list for EACH eigenvalue. *)
  neighborsList = nnFunction[#, 3] & /@ eigs;

  (* 4. Calculate ratios: (NN - k) / (NNN - k)
     Quiet suppresses warnings from degeneracies (e.g., 0/0 -> Indeterminate)
  *)
  Quiet[
    ( (#[[2]] - #[[1]]) / (#[[3]] - #[[1]]) ) & /@ neighborsList,
    {Power::infy, Infinity::indet}
  ]
]


(* Base calculation for a single time step *)
(* O(N) complexity using packed array vectorization *)
SpectralFormFactor[spectrumList_List, t_?NumericQ] := 
    Abs[Total[Exp[-I * spectrumList * t]]]^2 / Length[spectrumList];

(* Vectorized and memory-efficient calculation for a list of times *)
(* Map is preferred over Outer here to prevent massive RAM spikes for large spectra *)
SpectralFormFactor[spectrumList_List, tList_List] := 
    Map[SpectralFormFactor[spectrumList, #] &, tList];


(* Time Averaging using purely functional, built-in C-level routines *)
TimeAveragedSFF[tList_List, sffValues_List, windowSize_Integer?Positive] := Module[
    {tAvg, sffAvg},
    
    (* MovingAverage perfectly handles the sliding window without For loops *)
    tAvg = MovingAverage[tList, windowSize];
    sffAvg = MovingAverage[sffValues, windowSize];
    
    (* Return as coordinate pairs ready for ListPlot *)
    Transpose[{tAvg, sffAvg}]
];


NumberVariance[unfoldedSpectrum_List, l_?NumericQ, numSamples_Integer : 1000] := Module[
    {
        sortedSpectrum, minEnergy, maxEnergy, validMax, 
        empiricalDist, starts, counts, totalLevels
    },
    
    sortedSpectrum = Sort[unfoldedSpectrum];
    totalLevels = Length[sortedSpectrum];
    minEnergy = First[sortedSpectrum];
    maxEnergy = Last[sortedSpectrum];
    validMax = maxEnergy - l;
    
    (* Validaci\[OAcute]n para asegurar que el intervalo L cabe en el espectro *)
    If[validMax <= minEnergy, 
        Return[Missing["NotEnoughLevels"]]
    ];
    
    (* Optimizaci\[OAcute]n O(1) usando la CDF emp\[IAcute]rica vectorizada pura *)
    empiricalDist = EmpiricalDistribution[sortedSpectrum];
    
    (* Muestreo uniforme de los puntos de inicio *)
    starts = Subdivide[minEnergy, validMax, numSamples];
    
    (* C\[AAcute]lculo vectorizado del n\[UAcute]mero de niveles en cada intervalo [x, x + L] *)
    counts = totalLevels * (CDF[empiricalDist, starts + l] - CDF[empiricalDist, starts]);
    
    Variance[counts]
];

AnalyticalNumberVariancePoisson[l_?NumericQ] := l;

AnalyticalNumberVarianceGOE[l_?NumericQ] := 
    (2.0 / Pi^2) * (Log[2.0 * Pi * l] + EulerGamma + 1.0 - Pi^2 / 8.0);

AnalyticalNumberVarianceGUE[l_?NumericQ] := 
    (1.0 / Pi^2) * (Log[2.0 * Pi * l] + EulerGamma + 1.0);


(* ::Subsection::Closed:: *)
(*RMT*)


RatiosDistribution[r_, \[Beta]_] := 
    (r + r^2)^\[Beta] / ((1 + r + r^2)^(1 + 3 \[Beta]/2) * RatiosNormalizationZ[\[Beta]])

(* Helper function for normalization constant using a dummy variable 'x' *)
(* Uses Memoization: Z[beta] is calculated only once and stored *)
RatiosNormalizationZ[\[Beta]_] := RatiosNormalizationZ[\[Beta]] = 
    NIntegrate[(x + x^2)^\[Beta]/(1 + x + x^2)^(1 + 3*\[Beta]/2), {x, 0, Infinity}]

RatiosDistributionPoisson[r_, k_:1] := (2k -1)!*Power[r, k - 1]/(((k - 1)!)^2 * (1 + r)^(2*k))


NormalizedSpacingRatios[spectrum_, k_:1] := 
    Module[{r},
        r = SpacingRatios[spectrum, k];
        (* Vectorized Min[r, 1/r] calculation *)
        Min[#, 1/#] & /@ r
    ]


(* Core KL Divergence calculation *)
KLDivergence[p_?VectorQ, q_?VectorQ, base_: E] := 
    Module[{pNorm, qNorm, validIndices, pClean, qClean},
        (* 1. Normalize vectors to ensure they form probability distributions *)
        pNorm = N[p / Total[p]];
        qNorm = N[q / Total[q]];
        
        (* 2. Safety Check: Dimensions *)
        If[Length[pNorm] != Length[qNorm], 
            Message[KLDivergence::shlen]; Return[$Failed]
        ];
        
        (* 3. Filter: Only indices where P > 0 matter (limit x->0 of x log x is 0) *)
        (* Using Pick is faster than Select/Cases for packed arrays *)
        validIndices = Pick[Range[Length[pNorm]], UnitStep[pNorm - $MachineEpsilon], 1];
        pClean = pNorm[[validIndices]];
        qClean = qNorm[[validIndices]];
        
        (* 4. Safety Check: Absolute continuity. If P > 0 but Q = 0, D_KL is infinity *)
        If[AnyTrue[qClean, # <= $MachineEpsilon &], 
            Return[Infinity]
        ];
        
        (* 5. Vectorized Calculation *)
        pClean . Log[base, pClean / qClean]
    ];

KLDivergence::shlen = "Vectors P and Q must have the same length.";


(* Helper for RMT Data Analysis *)
EmpiricalKLDivergence[data_List, targetPDF_, binSpec_: "FreedmanDiaconis"] := 
    Module[{binData, binEdges, binCenters, pEmpirical, qTheoretical},
        
        (* 1. Generate Histogram Data from experimental/numerical data *)
        {binEdges, binData} = HistogramList[data, binSpec];
        
        (* 2. Calculate P (Empirical Distribution) *)
        (* Note: binData contains counts. KLDivergence handles normalization internally. *)
        pEmpirical = N[binData]; 
        
        (* 3. Calculate Q (Theoretical Distribution) *)
        (* Evaluate PDF at bin centers to approximate the discrete probability mass *)
        binCenters = MovingAverage[binEdges, 2];
        qTheoretical = Map[targetPDF, binCenters];
        
        (* 4. Compute KL Divergence *)
        KLDivergence[pEmpirical, qTheoretical]
    ];


(* ::Subsubsection::Closed:: *)
(*Ginibre matrices*)


(* Default option definition *)
Options[GenerateGinibreMatrix] = {
  Ensemble -> "Unitary"
};

GenerateGinibreMatrix[n_Integer?Positive, opts : OptionsPattern[]] := Module[
  {
    type = OptionValue[Ensemble],
    mat, A, B
  },

  Switch[type,

    "Unitary",
    (* GinUE (beta=2): n x n complex.
       Entries are N(0, 1/2) + i*N(0, 1/2) *)
    mat = (1/Sqrt[2]) * (
      RandomVariate[NormalDistribution[], {n, n}] +
      I*RandomVariate[NormalDistribution[], {n, n}]
    ),

    "Orthogonal",
    (* GinOE (beta=1): n x n real. Entries are N(0, 1) *)
    mat = RandomVariate[NormalDistribution[], {n, n}],

    "Symplectic",
    (* GinSE (beta=4): 2n x 2n complex representation of n x n quaternion.
       Built from two independent n x n GinUE matrices A and B. *)
    A = (1/Sqrt[2]) * (
      RandomVariate[NormalDistribution[], {n, n}] +
      I*RandomVariate[NormalDistribution[], {n, n}]
    );
    B = (1/Sqrt[2]) * (
      RandomVariate[NormalDistribution[], {n, n}] +
      I*RandomVariate[NormalDistribution[], {n, n}]
    );
    (* Construct the 2n x 2n matrix: [[A, B], [-B^\[Dagger], A^\[Dagger]]] *)
    mat = ArrayFlatten[{
      {A, B},
      {-ConjugateTranspose[B], ConjugateTranspose[A]}
    }],

    _,
    (* Error handling for unknown type *)
    Message[GenerateGinibreMatrix::invalidEnsemble, type];
    mat = $Failed
  ];

  Return[mat]
]


(* ::Subsection::Closed:: *)
(*Quantum channels*)


(* ::Subsubsection::Closed:: *)
(*Reshuffle*)


Reshuffle[m_] := ArrayFlatten[ArrayFlatten/@Partition[Partition[ArrayReshape[#,{Sqrt[Dimensions[m][[1]]],Sqrt[Dimensions[m][[1]]]}]&/@m,Sqrt[Dimensions[m][[1]]]],Sqrt[Dimensions[m][[1]]]],1];


Reshuffle[A_,m_,n_] := ArrayFlatten[ArrayReshape[A, {m, n, m, n}]]


(* ::Subsection::Closed:: *)
(*Bosons*)


FockBasisStateAsColumnVector[FockBasisState_,N_,L_]:=Normal[SparseArray[Position[SortFockBasis[Normal[FockBasis[N,L]]][[2]],FockBasisState]->1,Binomial[N+L-1,L]]]


(* ------------------------------------------------------
FockBasis[N,M] returns the lexicographical-sorted Fock basis for N bosons in M sites.
New implementation as of 15/Jun/2025 uses more efficient algorithm based on IntegerPartitions.
------------------------------------------------------ *)
FockBasis[N_, M_] := ReverseSort[Catenate[Permutations[PadRight[#, M]] & /@ IntegerPartitions[N, M]]]

(* Old implementation preserved for reference *)
(*
FockBasis[N_, M_] := Module[{k, fockState},
    k = 1;
    Normal[Join[
        {fockState = SparseArray[{1 -> N}, {M}]},
        Table[
            fockState = SparseArray[Join[Table[i -> fockState[[i]], {i, k - 1}], {k -> fockState[[k]] - 1}], {M}];
            fockState[[k + 1]] = N - Total[fockState[[1 ;; k]]];
            k = Assignationk[M, N, fockState];
            fockState,
        HilbertSpaceDim[N, M] - 1]
    ]]
]
*)


SortFockBasis[fockBasis_]:=Transpose[Sort[{Tag[#],#}&/@fockBasis]]


(* ::Subsection::Closed:: *)
(*Bose-Hubbard*)


(* Mensajes para tipos incorrectos *)
BoseHubbardHamiltonian::int = 
    "n (`1`) and L (`2`) are expected to be integers.";
BoseHubbardHamiltonian::real = 
    "Hopping J (`1`) and interaction parameters U (`2`) are expected " <>
    "to be real numbers.";
BoseHubbardHamiltonian::badSymmetricSubspace = 
    "Opci\[OAcute]n SymmetricSubspace `1` inv\[AAcute]lida. " <>
    "Opciones v\[AAcute]lidas: \"All\", \"EvenParity\" o \"OddParity\".";

Options[BoseHubbardHamiltonian] = {
    SymmetricSubspace -> "All" (* "All"|"EvenParity"|"OddParity" *),
    Version -> "New"
};

OptionValuePatterns[BoseHubbardHamiltonian] = {
  SymmetricSubspace -> Alternatives["All", "EvenParity", "OddParity"],
  Version -> _
};

BoseHubbardHamiltonian[n_Integer, L_Integer, J_Real, U_Real, 
    OptionsPattern[]] := 
Module[
    {tags, basis, basiseven, rbasiseven, rbasisodd, basisodd, H, T, V, map},
    
    basis = N[FockBasis[n, L]];
    H = -J*KineticTermBoseHubbardHamiltonian[basis] + 
        U/2*PotentialTermBoseHubbardHamiltonian[basis];
    
    (* Para subespacios de simetria, obtenerlos a partir de H *)
    Switch[OptionValue[SymmetricSubspace],
        "All",
            Nothing,
            
        "EvenParity",
            basiseven = DeleteDuplicatesBy[basis, Sort[{#, Reverse[#]}]&];
            rbasiseven = Reverse /@ basiseven;
            map = AssociationThread[
                basis -> Range[BoseHubbardHilbertSpaceDimension[n, L]]];
            H = 1/2 # . (H[[#1,#1]] + H[[#1,#2]] + H[[#2,#1]] + H[[#2,#2]] & @@ 
                Map[map, {basiseven, rbasiseven}, {2}]) . # &[
                DiagonalMatrix[
                    ReplacePart[
                        ConstantArray[1., Length[basiseven]],
                        Thread[
                            Flatten[Position[
                                basiseven, 
                                _?(PalindromeQ[#] &), 
                                {1}
                            ]] -> 1/Sqrt[2.]
                        ]
                    ],
                    TargetStructure -> "Sparse"
                ]
            ],
            
        "OddParity",
            basisodd = DeleteDuplicatesBy[
                Discard[basis, PalindromeQ], 
                Sort[{#, Reverse[#]}]&];
            rbasisodd = Reverse /@ basisodd;
            map = AssociationThread[
                basis -> Range[BoseHubbardHilbertSpaceDimension[n, L]]];
            H = 1/2 (H[[#1,#1]] - H[[#1,#2]] - H[[#2,#1]] + H[[#2,#2]] & @@ 
                Map[map, {basisodd, rbasisodd}, {2}]),
                
        _,
            Message[
                BoseHubbardHamiltonian::badSymmetricSubspace, 
                OptionValue[SymmetricSubspace]
            ];
            Return[$Failed];
    ];
    
    H
]

(* Handle cases where arguments don't match the expected types *)
BoseHubbardHamiltonian[N_, L_, J_, U_, OptionsPattern[]] := Module[{},
  If[!IntegerQ[N] || !IntegerQ[L],
    Message[BoseHubbardHamiltonian::int, N, L];
    Return[$Failed];
  ];
  If[Head[J] =!= Real || Head[U] =!= Real,
    Message[BoseHubbardHamiltonian::real, J, U];
    Return[$Failed];
  ];
];

SyntaxInformation[BoseHubbardHamiltonian] = <|
  "ArgumentsPattern" -> {_, _, _, _, OptionsPattern[]}
|>;


ClearAll[PotentialTermBoseHubbardHamiltonian];

(* Mensajes de error *)
PotentialTermBoseHubbardHamiltonian::nint = 
    "El n\[UAcute]mero de part\[IAcute]culas n (`1`) debe ser un entero positivo.";
PotentialTermBoseHubbardHamiltonian::lint = 
    "El n\[UAcute]mero de sitios L (`1`) debe ser un entero positivo.";
PotentialTermBoseHubbardHamiltonian::int = 
    "Los par\[AAcute]metros n (`1`) y L (`2`) deben ser enteros positivos.";
PotentialTermBoseHubbardHamiltonian::empty = 
    "La base proporcionada est\[AAcute] vac\[IAcute]a.";
PotentialTermBoseHubbardHamiltonian::dim = 
    "Todos los estados en la base deben tener la misma longitud.";
    
PotentialTermBoseHubbardHamiltonian[n_Integer?Positive, L_Integer?Positive] := 
DiagonalMatrix[
    Total /@ ((#^2 - #) &[FockBasis[n, L]]),
    TargetStructure -> "Sparse"
]

(* Versi\[OAcute]n con validaci\[OAcute]n de basis *)
PotentialTermBoseHubbardHamiltonian[basis_List] := 
Module[
    {lengths},
    If[basis === {},
        Message[PotentialTermBoseHubbardHamiltonian::empty];
        Return[$Failed]
    ];
    
    lengths = Length /@ basis;
    If[!AllTrue[lengths, EqualTo[First[lengths]]],
        Message[PotentialTermBoseHubbardHamiltonian::dim];
        Return[$Failed]
    ];
    
    DiagonalMatrix[
        Total /@ ((#^2 - #) &[basis]),
        TargetStructure -> "Sparse"
    ]
]

(* Versi\[OAcute]n con validaci\[OAcute]n de par\[AAcute]metros *)
PotentialTermBoseHubbardHamiltonian[n_Integer, L_Integer] := 
If[n <= 0 || L <= 0,
    Message[PotentialTermBoseHubbardHamiltonian::int, n, L];
    Return[$Failed]
];

PotentialTermBoseHubbardHamiltonian[n_, L_] /; !IntegerQ[n] := 
(
    Message[PotentialTermBoseHubbardHamiltonian::nint, n];
    $Failed
)

PotentialTermBoseHubbardHamiltonian[n_, L_] /; !IntegerQ[L] := 
(
    Message[PotentialTermBoseHubbardHamiltonian::lint, L];
    $Failed
)

PotentialTermBoseHubbardHamiltonian[___] := 
(
    Message[PotentialTermBoseHubbardHamiltonian::usage];
    $Failed
)


KineticTermBoseHubbardHamiltonian::badSymmetricSubspace = 
    "Opci\[OAcute]n SymmetricSubspace `1` inv\[AAcute]lida. " <>
    "Opciones v\[AAcute]lidas: \"All\", \"EvenParity\" o \"OddParity\".";

Options[KineticTermBoseHubbardHamiltonian] = {
    SymmetricSubspace -> "All" (* "All"|"EvenParity"|"OddParity" *)
};

KineticTermBoseHubbardHamiltonian[basis_, OptionsPattern[]] := 
Module[
    {len = Length[basis], basisNumRange, T, basiseven, rbasiseven, 
     basisodd, rbasisodd, map},
    
    basisNumRange = Range[len];
    T = # + ConjugateTranspose[#] & [
        SparseArray[
            AssignationRulesKinetic[basis, 
             AssociationThread[basis -> basisNumRange], basisNumRange], 
            {len, len}
        ]
    ];
    
    Switch[OptionValue[SymmetricSubspace],
        "All",
            Nothing,
        "EvenParity",
            basiseven = DeleteDuplicatesBy[basis, Sort[{#, Reverse[#]}]&];
            rbasiseven = Reverse /@ basiseven;
            map = AssociationThread[basis -> Range[len]];
            T = 1/2 # . (T[[#1,#1]] + T[[#1,#2]] + T[[#2,#1]] + T[[#2,#2]] & @@ 
                Map[map, {basiseven, rbasiseven}, {2}]) . # &[
                DiagonalMatrix[
                    ReplacePart[
                        ConstantArray[1., Length[basiseven]],
                        Thread[
                            Flatten[Position[
                                basiseven, 
                                _?(PalindromeQ[#] &), 
                                {1}
                            ]] -> 1/Sqrt[2.]
                        ]
                    ],
                    TargetStructure -> "Sparse"
                ]
            ],
        "OddParity",
            basisodd = DeleteDuplicatesBy[
                Discard[basis, PalindromeQ], Sort[{#, Reverse[#]}]&];
            rbasisodd = Reverse /@ basisodd;
            map = AssociationThread[basis -> Range[len]];
            T = 1/2 (T[[#1,#1]] - T[[#1,#2]] - T[[#2,#1]] + T[[#2,#2]] & @@ 
                Map[map, {basisodd, rbasisodd}, {2}]),
            
        _,
            Message[
                KineticTermBoseHubbardHamiltonian::badSymmetricSubspace,
                OptionValue[SymmetricSubspace]
            ];
            Return[$Failed];
    ];
    
    T
]


AssignationRulesKinetic[basis_, positionMap_, basisNumRange_] := 
Catenate[
    MapThread[
        AssignationRulesKineticMapFunc,
        {
            Apply[
                {positionMap[#1], #2} &,
                DeleteCases[
                    Transpose[
                        {
                            StateAfterADaggerA[basis],
                            ValueAfterADaggerA[basis]
                        },
                        {3, 1, 2}
                    ],
                    {_, 0.},
                    {2}
                ],
                {2}
            ],
            basisNumRange
        }
    ]
]


StateAfterADaggerA[basis_] := 
Module[{len = Length[First[basis]]},
    Outer[
        Plus,
        basis,
        #,
        1
    ] & [
        Catenate[
            NestList[
                RotateRight,
                PadRight[#, len],
                len - 2
            ] & /@ {{1, -1}}
        ]
    ]
]


ValueAfterADaggerA[basis_] := 
MapApply[
    Sqrt[(#1 + 1.) * #2] &,
    Partition[#, 2, 1]
] & /@ basis


AssignationRulesKineticMapFunc[stateValuePairs_,index_] := 
({index, #1} -> #2) & @@@ stateValuePairs


Tag[FockBasisElement_]:=N[Round[Sum[Sqrt[100 i+3]#[[i+1]],{i,0,Length[#]-1}]&[FockBasisElement],10^-8]]
Tag[Nothing]:=Nothing


InitializationBosonicPartialTrace[SitesOfSubsystem_, n_, L_] :=
Module[{SubsystemSize,SubsystemComplementSize,SystemFockBasisTags,SystemFockBasis,SubsystemFockBasisTags,SubsystemFockBasis,SubsystemComplementFockBasis,RulesForOrderedSystemBasis,RulesForOrderedSubsystemBasis,FockIndicesInRho,FockIndicesInReducedRho},
SubsystemSize=Length[SitesOfSubsystem];
SubsystemComplementSize=L-SubsystemSize;
(*System's Fock basis*)
{SystemFockBasisTags,SystemFockBasis}=SortFockBasis[FockBasis[n,L]];
(*Subsystem's Fock basis*)
{SubsystemFockBasisTags,SubsystemFockBasis}=SortFockBasis[Flatten[Table[FockBasis[k,SubsystemSize],{k,0,n}],1]];
(*Complement subsystem's Fock basis*)
SubsystemComplementFockBasis=Map[ReplacePart[ConstantArray[_,L],Thread[Complement[Range[L],SitesOfSubsystem]->#]]&,Flatten[Table[SortFockBasis[FockBasis[k,SubsystemComplementSize]][[2]],{k,0,n}],1]];(*<<<*)
RulesForOrderedSystemBasis=Thread[Rule[SystemFockBasisTags,Range[BoseHubbardHilbertSpaceDimension[n,L]]]];
SubsystemHilbertSpaceDim=Length[SubsystemFockBasis];
RulesForOrderedSubsystemBasis=Thread[Rule[SubsystemFockBasisTags,Range[SubsystemHilbertSpaceDim]]];
FockIndicesInRho=Map[Tuples[{#,#}]&[Extract[SystemFockBasis,Position[SystemFockBasis,#]]]&,SubsystemComplementFockBasis];(*<<<*)
FockIndicesInReducedRho=Map[#[[All,All,SitesOfSubsystem]]&,FockIndicesInRho];(*<<<*)
ComputationalIndicesInRho=ReplaceAll[Map[Tag,FockIndicesInRho,{3}],RulesForOrderedSystemBasis];(*<<<*)
ComputationalIndicesInReducedRho=ReplaceAll[Map[Tag,FockIndicesInReducedRho,{3}],RulesForOrderedSubsystemBasis];(*<<<*)];


BosonicPartialTrace[Rho_] := 
	Module[{MatrixElementsOfRho,rules},
	
		MatrixElementsOfRho=Extract[Rho,#]&/@ComputationalIndicesInRho;
		rules=MapThread[Thread[Rule[#1,#2]]&,{ComputationalIndicesInReducedRho,MatrixElementsOfRho}];
		Total[Map[SparseArray[#,{SubsystemHilbertSpaceDim,SubsystemHilbertSpaceDim}]&,rules]]
	]


(* ::Subsubsection::Closed:: *)
(*Fuzzy measurements in bosonic systems*)


(* Initialization function *)
InitializeVariables[n_, L_, boundaries_, FMmodel_] := 
 Module[{basis, SwapAppliedToBasis},
  
  indices = Which[
    boundaries == "open",
    Table[ReplacePart[Range[L], {i -> Mod[i + 1, L, 1], Mod[i + 1, L, 1] -> i}], {i, L - ToExpression[StringTake[FMmodel, 1]]}],
    boundaries == "closed" || boundaries == "close",
    Print["Yet to come"]
  ];
  
  permutedBasisIndices = Table[
    basis = SortFockBasis[FockBasis[n, L]][[2]];
    SwapAppliedToBasis = #[[i]] & /@ basis;
    Flatten[Position[basis, #] & /@ SwapAppliedToBasis],
    {i, indices}
  ];
  
  numberOfPermutations = Length[permutedBasisIndices];
];


FuzzyMeasurement[\[Psi]_,pFuzzy_] := 
(1 - pFuzzy) Dyad[\[Psi]] + (pFuzzy/numberOfPermutations) * Total[ Table[ Dyad[\[Psi][[i]]], {i,permutedBasisIndices}]]


(* ::Subsection::Closed:: *)
(*BosonEscapeKrausOperators[n, L]*)


BosonEscapeKrausOperators[n_,L_] := Module[{dimH=BoseHubbardHilbertSpaceDimension[n,L],fockBasisTags,fockBasisStates,KrausOperators},
{fockBasisTags,fockBasisStates}=SortFockBasis[FockBasis[n,L]];

KrausOperators=1/Sqrt[2(n-1)]*Join[
Table[Normal[SparseArray[Thread[
Select[Table[Module[{fockState=fockBasisStates[[k]]},
{FockBasisIndex[SigmaPlusSigmaMinus[fockState,i,L],fockBasisTags],k}],
{k,dimH}],AllTrue[#,Positive]&]->1.],{dimH,dimH}]],
{i,n-1}],
Table[Normal[SparseArray[Thread[
Select[Table[Module[{fockState=fockBasisStates[[k]]},
{FockBasisIndex[SigmaMinusSigmaPlus[fockState,i,L],fockBasisTags],k}],
{k,dimH}],AllTrue[#,Positive]&]->1.],{dimH,dimH}]],
{i,n-1}]];
Append[KrausOperators,Sqrt[IdentityMatrix[dimH]-Sum[ConjugateTranspose[K] . K,{K,KrausOperators}]]]
]


SigmaPlusSigmaMinus[fockState_,i_,L_] := 
If[fockState[[i]]==L||fockState[[i+1]]==0,
Nothing,
ReplaceAt[ReplaceAt[fockState,x_:>x+1,i],x_:>x-1,i+1]
]


SigmaMinusSigmaPlus[fockState_,i_,L_] := 
If[fockState[[i]]==0||fockState[[i+1]]==L,
Nothing,
ReplaceAt[ReplaceAt[fockState,x_:>x-1,i],x_:>x+1,i+1]
]


(*Secondary routines*)


HilbertSpaceDim[args___] := Message[HilbertSpaceDim::replaced, "HilbertSpaceDim", "BoseHubbardHilbertSpaceDimension"];


BoseHubbardHilbertSpaceDimension[n_,L_]:=(n+L-1)!/(n!(L-1)!)


FockBasisIndex[fockState_,sortedTagsFockBasis_]:=
FromDigits[Flatten[Position[sortedTagsFockBasis,Tag[fockState]],1]]


(*Check definitions*)


Assignationk[M_,N_,n_]:=If[n[[1;;M-1]]==ConstantArray[0,M-1],M-1,FromDigits[Last[Position[Normal[n[[1;;M-1]]],x_ /;x!=0]]]]


RenyiEntropy[\[Alpha]_,\[Rho]_]:=1/(1-\[Alpha]) Log[Tr[MatrixPower[\[Rho],\[Alpha]]]]


(* ::Subsection:: *)
(*Spin chains*)


(* ::Subsubsection::Closed:: *)
(*Symmetries*)


SpinParityEigenvectors[L_]:=Module[{tuples,nonPalindromes,palindromes},
tuples=Tuples[{0,1},L];
nonPalindromes=Select[tuples,#!=Reverse[#]&];
palindromes=Complement[tuples,nonPalindromes];
nonPalindromes=DeleteDuplicatesBy[nonPalindromes,Sort[{#,Reverse[#]}]&];
Normal[
{
Join[SparseArray[FromDigits[#,2]+1->1.,2^L]&/@palindromes,Normalize[SparseArray[{FromDigits[#,2]+1->1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes],
Normalize[SparseArray[{FromDigits[#,2]+1->-1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes
}]
]


Translation[state_, L_] := BitShiftRight[state] + BitAnd[state, 1]*2^(L-1)
BinaryNecklaces[L_Integer] := Module[{tuples=Tuples[{0,1},L]},
	Union[Table[First[Sort[NestList[RotateLeft, t, L-1]]], {t,tuples}]]
]
\[Omega][L_,k_] := Exp[2Pi*k*I/L]

TranslationEigenvectorRepresentatives[L_Integer] := Module[
	{
	necklaces = FromDigits[#,2]& /@ BinaryNecklaces[L],
	orbits
	},
	
	orbits = DeleteDuplicates[
	Sort[NestWhileList[Translation[#, L]&, #, UnsameQ[##]&, All][[;;-2]]]& /@ necklaces
	 ];
	Catenate[
		Outer[
			If[Mod[Length[#1]*#2, L]==0, {First[#1], #2, Length[#1]}, Nothing]&, 
			orbits, 
			Range[0, L-1], 
			1
		]
	]
]


RepresentativesOddBasis[basis_]:=DeleteDuplicatesBy[Discard[basis,PalindromeQ],Sort[{#,Reverse[#]}]&]
RepresentativesEvenBasis[basis_]:=DeleteDuplicatesBy[basis,Sort[{#,Reverse[#]}]&]


Options[BlockDiagonalize] = {
  Symmetry -> "Translation" (*coming soon: \"Parity\"*)
};

BlockDiagonalize[A_, opts:OptionsPattern[]]:= 
Switch[OptionValue[Symmetry],
	"Translation",
		Module[
			{
				L = Log[2, Length[A]],
				repsgathered,
				P
			},
				(*Gather translation eigenvectors by their pseudomomentum k*)
				repsgathered = GatherBy[TranslationEigenvectorRepresentatives[L], #[[2]]&];
				
				(*change of basis matrix*)
				P = SparseArray[
						Catenate[
							Apply[
								Normalize[SparseArray[
									Thread[
										FoldList[Translation[#, L]&, #1, Range[#3 - 1]] + 1 -> 
										Power[\[Omega][L, #2], Range[0., #3 - 1]]
									]
										, 2^L
								]]&,
								repsgathered,
								{2}
							]
						]
					];
					
				BlockDiagonalMatrix[Chop[Conjugate[P] . A . Transpose[P]]]
		],
	"Parity",
		Module[
		{
			L = Log[2, Length[A]],
			basis, reps, rules, Heven, Hodd
		},
			basis = Tuples[{0, 1}, L];
			rules=AssociationThread[basis->Range[Length[basis]]];
			reps=Comap[{RepresentativesEvenBasis, RepresentativesOddBasis},basis];
			Heven=1/2 (A[[#1,#1]] + A[[#1,#2]] + A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[1]]], {2}]);
			Hodd=1/2 (A[[#1,#1]] - A[[#1,#2]] - A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[2]]], {2}]);
			Heven = # . Heven . #&[DiagonalMatrix[ReplacePart[ConstantArray[1.,Length[reps[[1]]]],Thread[Catenate[Position[reps[[1]],_?PalindromeQ,1]]->1/Sqrt[2.]]],TargetStructure->"Sparse"]];
			
			{Heven, Hodd}
		],
	_,
		Message[BlockDiagonalize::badSymmetry, OptionValue[Symmetry]];
		Return[$Failed];
]

(*Mensaje de error si la opci\[OAcute]n es inv\[AAcute]lida*)
BlockDiagonalize::badSymmetry = 
  "Option badSymmetry -> `1` is not valid. Valid options: \"Translation\".";


(* ::Subsubsection::Closed:: *)
(*Hamiltonians*)


Options[IsingHamiltonian] = {
  BoundaryConditions -> "Open"
};

IsingHamiltonian[hx_, hz_, J_, L_, opts:OptionsPattern[]] := Module[
	{NNIndices},
	NNIndices = Switch[OptionValue[BoundaryConditions],
		"Open",
			Normal[SparseArray[Thread[{#,#+1}->3],{L}]&/@Range[L-1]],
		"Periodic",
			Normal[SparseArray[Thread[{#,Mod[#+1,L,1]}->3],{L}]&/@Range[L]],
		_,
			Message[
                IsingHamiltonian::badBoundaryCondition, 
                OptionValue[BoundaryConditions]
            ];
            Return[$Failed];
	];
	
	Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]
]

(*Mensaje de error si la opci\[OAcute]n es inv\[AAcute]lida*)
IsingHamiltonian::badBoundaryCondition = 
  "Option BoundaryConditions -> `1` not valid. Valid options: \"Open\" o \"Periodic\".";


IsingNNOpenHamiltonian[args___] := Message[
	IsingNNOpenHamiltonian::replaced, "IsingNNOpenHamiltonian", "IsingHamiltonian"];
	
(*IsingNNOpenHamiltonian[hx_,hz_,J_,L_] := Module[{NNIndices},
	NNIndices=Normal[SparseArray[Thread[{#,#+1}->3],{L}]&/@Range[L-1]];
	N[Normal[Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]]]]*)


IsingNNClosedHamiltonian[args___] := Message[
	IsingNNClosedHamiltonian::replaced, "IsingNNClosedHamiltonian", "IsingHamiltonian"];

(*IsingNNClosedHamiltonian[hx_,hz_,J_,L_] := 
Module[{NNIndices},
	NNIndices=Normal[SparseArray[Thread[{#,Mod[#+1,L,1]}->3],{L}]&/@Range[L]];
	Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]
]*)


ClosedXXZHamiltonian[L_,\[CapitalDelta]_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L] ];
		N[Normal[-1/2*Total[Join[Pauli/@NNindices,Pauli/@(2NNindices),\[CapitalDelta] (Pauli[#]-IdentityMatrix[2^L])&/@(3NNindices)]]]]
	]


OpenXXZHamiltonian[L_,\[CapitalDelta]_,h1_,h2_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L-1] ];
		N[Normal[-1/2*Total[Join[Pauli/@NNindices,Pauli/@(2NNindices),\[CapitalDelta] (Pauli[#]-IdentityMatrix[2^L])&/@(3NNindices)]]  
		- 1/2*(h1 Pauli[Join[{1},ConstantArray[0,L-1]]] + h2*Pauli[Join[ConstantArray[0,L-1],{1}]])+ 1/2*(h1 + h2)IdentityMatrix[2^L]]]
	]


HamiltonianNN[Jxy_,Jz_,L_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L-1] ];
		N[Normal[(1/4)*Total[Join[Jxy*(Pauli/@NNindices),Jxy*(Pauli/@(2NNindices)),Jz*(Pauli[#]&/@(3NNindices))]]]]
	]

HamiltonianZ[\[Omega]_,\[Epsilon]d_,L_,d_]:=N[(1/2)*(\[Omega]*Total[Pauli/@(3*IdentityMatrix[L])]+\[Epsilon]d*Pauli[Normal[SparseArray[d->3,L]]])]

LeaSpinChainHamiltonian[Jxy_,Jz_,\[Omega]_,\[Epsilon]d_,L_,d_]:=HamiltonianNN[Jxy,Jz,L]+HamiltonianZ[\[Omega],\[Epsilon]d,L,d]
XXZOpenHamiltonian[Jxy_,Jz_,\[Omega]_,\[Epsilon]d_,L_,d_]:=HamiltonianNN[Jxy,Jz,L]+HamiltonianZ[\[Omega],\[Epsilon]d,L,d]


HeisenbergXXXwNoise[h_List,L_]:=
Module[{NNIndices,firstSum,secondSum},
(* \sum_{k=1}^{L-1} S_k^xS_{k+1}^x + S_k^zS_{k+1}^z + S_k^zS_{k+1}^z *)
NNIndices=Normal[SparseArray[Thread[{#,#+1}->1],{L}]&/@Range[L-1]];
firstSum=1/4*Total[Table[Pauli[i*#]&/@NNIndices,{i,3}],2];

(* \sum_{k=1}^{L} h_k^z S_k^z *)
secondSum=1/2*h . (Pauli/@DiagonalMatrix[ConstantArray[3,L]]);

firstSum+secondSum
]


(* Quantum Game of Life *)

Options[QuantumGameOfLifeHamiltonian] = {
    Momentum -> "All",
    Parity -> "None" (* "None", 1, or -1 *)
};

QuantumGameOfLifeHamiltonian[L_Integer, OptionsPattern[]] := Module[
    {
        k = OptionValue[Momentum],
        p = OptionValue[Parity],
        dim = 2^L,
        Hfull,
        rules,
        basisVectors,
        P (* Projection Matrix *)
    },

    (* 1. VALIDACIONES *)
    If[k =!= "All" && !IntegerQ[k], Return[Message[QuantumGameOfLifeHamiltonian::invMom, k]; $Failed]];
    If[p =!= "None" && !MemberQ[{1, -1}, p], Return[Message[QuantumGameOfLifeHamiltonian::invPar, p]; $Failed]];
    
    If[IntegerQ[k] && p =!= "None",
        If[k != 0 && k != L/2,
            Message[QuantumGameOfLifeHamiltonian::parityWarning, k];
        ]
    ];

    (* 2. CONSTRUCCI\[CapitalOAcute]N DEL HAMILTONIANO COMPLETO (Bitwise) *)
    rules = Flatten @ Table[
        Module[{n = state, flips = {}, neighborsSum},
            Do[
                neighborsSum = 
                    BitGet[n, Mod[i - 2, L]] + 
                    BitGet[n, Mod[i - 1, L]] + 
                    BitGet[n, Mod[i + 1, L]] + 
                    BitGet[n, Mod[i + 2, L]];
                
                If[neighborsSum == 2 || neighborsSum == 3,
                    AppendTo[flips, {n + 1, BitXor[n, 2^i] + 1} -> 1.0]
                ],
                {i, 0, L - 1}
            ];
            flips
        ],
        {state, 0, dim - 1}
    ];
    
    Hfull = SparseArray[rules, {dim, dim}];

    (* 3. RETORNO SI NO HAY SIMETR\[CapitalIAcute]AS *)
    If[k === "All", Return[Hfull]];

    (* 4. CONSTRUCCI\[CapitalOAcute]N DE LA BASE DE MOMENTO k *)
    basisVectors = Module[{necklaces, kBasis},
        (* Agrupar por collares *)
        necklaces = GroupBy[Range[0, dim - 1], 
            Function[x, Min[NestList[RotateLeftBits[#, L]&, x, L-1]]]
        ];
        
        kBasis = Reap[
            Do[
                (* AQU\[CapitalIAcute] ESTABA EL ERROR: Elimin\[EAcute] 'T_k' y limpi\[EAcute] variables *)
                Module[{orbit = members, R, vec},
                    R = Length[orbit];
                    
                    (* Condici\[OAcute]n de compatibilidad de momento *)
                    If[Divisible[k * R, L],
                        vec = SparseArray[{}, dim];
                        Do[
                            (* Construcci\[OAcute]n del vector de Bloch *)
                            vec = vec + SparseArray[{orbit[[j+1]] + 1} -> Exp[-I * 2. * Pi * k * j / L], dim];
                        , {j, 0, R-1}];
                        
                        Sow[Normalize[vec]];
                    ]
                ],
                {members, Values[necklaces]}
            ]
        ][[2]];
        
        If[kBasis === {}, Return[SparseArray[{}, {0, 0}]]];
        Flatten[kBasis, 1]
    ];

    (* 5. PROYECCI\[CapitalOAcute]N DE PARIDAD *)
    If[p =!= "None",
        basisVectors = Reap[
            Scan[Function[v,
                Module[{vInv, vProj},
                    vInv = SparseArray[
                        Thread[(ReverseBits[# - 1, L] + 1 & /@ v["NonzeroPositions"][[All, 1]]) -> v["NonzeroValues"]], 
                        dim
                    ];
                    
                    vProj = v + p * vInv;
                    
                    If[Norm[vProj] > 10^-10, Sow[Normalize[vProj]]];
                ]
            ], basisVectors]
        ][[2]];
        
        If[basisVectors =!= {},
             (* Eliminaci\[OAcute]n de duplicados num\[EAcute]ricos *)
             basisVectors = DeleteDuplicates[Flatten[basisVectors, 1], (Norm[#1 - #2] < 10^-8 || Norm[#1 + #2] < 10^-8) &];
        ,
             Return[SparseArray[{}, {0, 0}]]
        ];
    ];

    (* 6. CONSTRUCCI\[CapitalOAcute]N FINAL *)
    P = Transpose[SparseArray[basisVectors]];
    Chop[ConjugateTranspose[P] . Hfull . P]
];

(* Mensajes y Funciones Auxiliares *)
QuantumGameOfLifeHamiltonian::invMom = "Momentum `1` must be \"All\" or an integer 0 <= k < L.";
QuantumGameOfLifeHamiltonian::invPar = "Parity `1` must be \"None\", 1, or -1.";
QuantumGameOfLifeHamiltonian::parityWarning = "Warning: Inversion Parity is generally not a good quantum number for momentum k=`1` (unless k=0 or k=L/2).";

RotateLeftBits[n_, L_] := BitOr[BitShiftLeft[BitAnd[n, 2^(L-1)-1], 1], BitShiftRight[n, L-1]];
ReverseBits[n_, L_] := FromDigits[Reverse[IntegerDigits[n, 2, L]], 2];


(* ::Subsection::Closed:: *)
(*Spin Hamiltonians*)


(* Helper function to generate sparse spin operators for arbitrary S *)
GenerateSpinOperators[s_] := Module[
    {d, range, diagonal, upperDiag, Sz, Sp, Sm, Sx, Sy, Id},
    
    (* Dimension of the Hilbert space for spin s *)
    d = Round[2 s + 1];
    
    (* Validate s is integer or half-integer *)
    If[Abs[d - (2 s + 1)] > 10^-10, Return[$Failed]];

    (* Basis range from S to -S *)
    range = Range[s, -s, -1];

    (* S_z: Diagonal matrix with m values *)
    Sz = SparseArray[Band[{1, 1}] -> range, {d, d}];

    (* S_plus: Elements Sqrt[s(s+1) - m(m+1)] on superdiagonal *)
    (* Note: range contains 'm'. The element <m+1|S+|m> is at position corresponding to transition m -> m+1 *)
    (* In matrix indices 1..d, index i is m_i. i-1 is m_i + 1. So it acts on column i to row i-1 *)
    upperDiag = Table[
        Sqrt[(s - m) (s + m + 1)], 
        {m, range[[2 ;;]]} (* Exclude the first m=s as it cannot be raised from *)
    ];
    Sp = SparseArray[Band[{1, 2}] -> upperDiag, {d, d}];
    
    (* S_minus: Transpose of S_plus *)
    Sm = Transpose[Sp];

    (* S_x and S_y *)
    Sx = (Sp + Sm) / 2;
    (* Sy = (Sp - Sm) / (2 I); Unused here but standard definition *)
    
    (* Identity *)
    Id = IdentityMatrix[d, SparseArray];

    <| "z" -> Sz, "x" -> Sx, "+" -> Sp, "-" -> Sm, "Id" -> Id |>
];

(* Main Hamiltonian Function *)
TwoSpinHamiltonian[s1_, s2_, \[Epsilon]1_, \[Epsilon]2_, \[Theta]_, g_] := Module[
    {
        ops1, ops2,
        term1, term2, interaction,
        H
    },
    
    (* 1. Generate local operators *)
    ops1 = GenerateSpinOperators[s1];
    ops2 = GenerateSpinOperators[s2];
    
    If[FailureQ[ops1] || FailureQ[ops2], 
        Return[Message[TwoSpinHamiltonian::invalidSpin, "{s1, s2}"]]
    ];

    (* 2. Term 1: \[Epsilon]1 * S1_z (x) Id_2 *)
    term1 = \[Epsilon]1 * KroneckerProduct[ops1["z"], ops2["Id"]];

    (* 3. Term 2: \[Epsilon]2 * Id_1 (x) (Cos[\[Theta]] S2_z + Sin[\[Theta]] S2_x) *)
    term2 = \[Epsilon]2 * KroneckerProduct[
        ops1["Id"], 
        Cos[\[Theta]] * ops2["z"] + Sin[\[Theta]] * ops2["x"]
    ];

    (* 4. Interaction Term: (g / Sqrt[S1 S2]) * (S1 . S2) *)
    (* Decomposed as SzSz + 1/2 (S+S- + S-S+) *)
    interaction = (g / Sqrt[s1 * s2]) * (
        KroneckerProduct[ops1["z"], ops2["z"]] + 
        0.5 * (
            KroneckerProduct[ops1["+"], ops2["-"]] + 
            KroneckerProduct[ops1["-"], ops2["+"]]
        )
    );

    (* 5. Total Sum *)
    H = term1 + term2 + interaction;

    (* Return SparseArray, optionally chopped if purely numerical *)
    If[AllTrue[{s1, s2, \[Epsilon]1, \[Epsilon]2, \[Theta], g}, NumericQ],
        Chop[H],
        H
    ]
];

TwoSpinHamiltonian::invalidSpin = "Spin values `1` must be integers or half-integers.";


(* ::Subsection::Closed:: *)
(*Light-matter systems*)


(* The following definition calculates the Hamiltonian for the Asymmetric quantum 
Rabi model (AQRM) in a truncated Fock space. Definition from arXiv:2508.00765v2*)
ARQMHamiltonian[\[Omega]_, \[CapitalDelta]_, g_, \[Epsilon]_, \[Xi]_, nMax_] :=
    Module[{
        idQ = IdentityMatrix[2],         (* Identity in the qubit space *)
        idB = IdentityMatrix[nMax + 1],  (* Identity in the bosonic space *)
        a, adag, X, Y, Z, H              (* Operators *)
    },

    (* === PAULI OPERATORS DEFINITION === *)
    {X, Y, Z} = Pauli /@ Range[3];

    (* === BOSONIC OPERATORS DEFINITION (In the truncated Fock space) === *)
    adag = 
    SparseArray[{Band[{2, 1}] -> Table[Sqrt[n], {n, 1, nMax}]}, {nMax + 1, nMax + 1}];

    a = Transpose[adag];

    (* === HAMILTONIAN CONSTRUCTION TERM BY TERM === *)

    (* Term 1: Photon energy (\[Omega] a\[Dagger]a \otimes idQ) *)
    H1 = \[Omega] * KroneckerProduct[idQ, adag . a];

    (* Term 2: Qubit energy (\[CapitalDelta] Z \otimes idB) *)
    H2 = \[CapitalDelta] * KroneckerProduct[Z, idB];

    (* Term 3: Coupling term 1 (g/2) (1+\[Xi]) (X \otimes (a + a\[Dagger])) *)
    H3 = g/2 * (1 + \[Xi]) * KroneckerProduct[X, (a + adag)];

    (* Term 4: Coupling term 2 (g/2) (1-\[Xi]) (iY \otimes (a - a\[Dagger])) *)
    H4 = g/2 * (1 - \[Xi]) * KroneckerProduct[I * Y, (a - adag)];

    (* Term 5: Asymmetry term (\[Epsilon] X \otimes idB) *)
    H5 = \[Epsilon] * KroneckerProduct[X, idB];

    (* === TOTAL HAMILTONIAN === *)
    H = H1 + H2 + H3 + H4 + H5;

    (* Return the Hamiltonian as a dense, numerically evaluated matrix,
       with small imaginary parts rounded to zero (Chop). *)
    Return[Chop[N[Normal[H]]]];
];


(* ::Subsection::Closed:: *)
(*Fuzzy measurement channels*)


SwapMatrix[targetSite_, wrongSite_, N_, L_] := Module[
    {
        tagsOfSwappedFockBasis, (* List to hold tags of swapped Fock basis *)
        sortedTags,             (* List to hold sorted tags *)
        sortedFockBasis         (* List to hold sorted Fock basis *)
    },
    
    (* Step 1: Generate and sort the Fock basis (according to its tags) *)
    {sortedTags, sortedFockBasis} = Normal[SortFockBasis[FockBasis[N, L]]];
    
    (* Step 2: Generate tags for the swapped Fock basis *)
    tagsOfSwappedFockBasis = Tag /@ Map[SwappedFockState[#, targetSite, wrongSite] &, sortedFockBasis];
    
    (* Step 3: Create the sparse matrix corresponding to the swap matrix *)
    SparseArray[
        Table[{i, Position[sortedTags, tagsOfSwappedFockBasis[[i]]][[1, 1]]},
        {i, Length[tagsOfSwappedFockBasis]}] -> 1]
]


FuzzyMeasurementChannel[\[Rho]_, pTotalError_, SwapNNmatrices_]:=(1 - pTotalError) \[Rho] +
	pTotalError 1/Length[SwapNNmatrices] Total[# . \[Rho] . # & /@ SwapNNmatrices]


FuzzyMeasurementChannel[\[Rho]_, pErrors_List, SwapMatrices_] := Module[
	{
		{pTotalError, pNN, pSNN} = pErrors,
		{SwapNN, SwapSNN} = SwapMatrices
	},
	
	(1 - pTotalError) \[Rho] + pNN ( pTotalError/Length[SwapNN] ) Total[# . \[Rho] . # & /@ SwapNN]
	+ pNN ( pTotalError/Length[SwapSNN] ) Total[# . \[Rho] . # & /@ SwapSNN]
]


(* ::Subsubsection:: *)
(*Private functions*)


SwappedFockState[fockState_, targetSite_, wrongSite_] := 
  ReplacePart[fockState, 
    Thread[# -> Part[fockState, Reverse[#]] ] & [{targetSite, wrongSite}]
  ]


(* ::Section:: *)
(*End of Package*)


End[];


EndPackage[];
