(* ::Package:: *)

BeginPackage["QMB`"];


Get["ForScience`"]; (* nice formatting of usage definitions *)


(* ::Section::Closed:: *)
(*TODO*)


(* ::Text:: *)
(*Hay cosas en IAA_model.nb que 1) hay que migrar para ac\[AAcute] y 2) que hay que revisar si deber\[IAcute]a de poner ac\[AAcute]*)


(* ::Text:: *)
(*Hay cosas de Heisenberg meets fuzzy que tambi\[EAcute]n tengo que pasar para ac\[AAcute]*)


(* ::Section:: *)
(*Public definitions*)


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


RenyiEntropy::usage = "RenyiEntropy[\[Alpha], \[Rho]] computes the \[Alpha]-th order Renyi entropy of density matrix \[Rho].";


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


LevelSpacingDistribution::usage = FormatUsage[
"LevelSpacingDistribution[s, \[Beta]] returns the Wigner surmise probability \
density function P(s) for the normalized nearest-neighbor level spacings. \
The Dyson index ```\[Beta]``` determines the ensemble: 0 (Poisson), 1 (GOE), \
2 (GUE), or 4 (GSE)."];


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



HilbertSpaceDim::replaced = "Function `1` has been replaced by `2`.";


BoseHubbardHilbertSpaceDimension::usage = FormatUsage["BoseHubbardHilbertSpaceDimension[n,L] returns the dimension"<> 
"of Hilbert space of a Bose Hubbard system of N bosons and L sites."];


FockBasis::usage = "FockBasis[N, M] returns the lexicographical-sorted Fock basis of N bosons in M sites.";


SortFockBasis::usage = "SortFockBasis[fockBasis] returns fockBasis in ascending-order according to the tag of Fock states.";


Tag::usage = "Tag[ { \!\(\*SubscriptBox[\(k\), \(1\)]\),\!\(\*SubscriptBox[\(k\), \(2\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(k\), \(L\)]\) } ] returns the tag \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\)\!\(\*SqrtBox[\(100  i + 3\)]\)\!\(\*SubscriptBox[\(k\), \(i\)]\) of Fock state \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"1\"], \",\", SubscriptBox[\"k\", \"2\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"L\"]}]},\n\"Ket\"]\).";


FockBasisStateAsColumnVector::usage = "FockBasisStateAsColumnVector[FockState, N, L] returns matrix representation of fockState={\!\(\*SubscriptBox[\(i\), \(1\)]\),\!\(\*SubscriptBox[\(i\), \(2\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(L\)]\)}. N: bosons, L: sites.";


FockBasisIndex::usage = "FockBasisIndex[fockState, sortedTagsFockBasis] returns the position of fockState in the tag-sorted Fock basis with tags sortedTagsFockBasis.";


BosonicPartialTrace::usage = "BosonicPartialTrace[\[Rho]] calculates the partial trace of \[Rho]. Requires initialization.";


InitializationBosonicPartialTrace::usage = "InitializationBosonicPartialTrace[{\!\(\*SubscriptBox[\(i\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(k\)]\)}, N, L] initializes variables for BosonicPartialTrace[] to calculate the reduced density matrix of sites {\!\(\*SubscriptBox[\(i\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(k\)]\)}.";


(* ::Subsection::Closed:: *)
(*Fuzzy measurements in bosonic systems*)


BosonEscapeKrausOperators::usage = "BosonEscapeKrausOperators[N, L]: bosons escape to nearest neighbouring sites. N: bosons, L: site.";


InitializeVariables::usage = "InitializeVariables[n, L, boundaries, FMmodel] sets up the necessary variables for correct running of FuzzyMeasurement[\[Psi], \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\)]; boundaries: 'open' or 'closed'; FMmodel: '#NN'.";


FuzzyMeasurement::usage = "FuzzyMeasurement[\[Psi], \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\)] gives \[ScriptCapitalF](\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\)) = (1 - \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\))\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\) + \!\(\*SubscriptBox[\(p\), \(fuzzy\)]\) \!\(\*UnderscriptBox[\(\[Sum]\), \(i\)]\) \!\(\*SubscriptBox[\(S\), \(i\)]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\)\!\(\*SubsuperscriptBox[\(S\), \(i\), \(\[Dagger]\)]\), where \!\(\*SubscriptBox[\(S\), \(i\)]\) must be initizalized runnning InitializeVariables[n, L, boundaries, FMmodel].";


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
(*Private definitions*)


Begin["`Private`"]


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


RenyiEntropy[\[Alpha]_,\[Rho]_]:=1/(1-\[Alpha]) Log[Tr[MatrixPower[\[Rho],\[Alpha]]]]


(* ::Subsection::Closed:: *)
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


LevelSpacingDistribution[Spacing_, DysonIndex_] := 
    Switch[DysonIndex,
        0, 
            Exp[-Spacing],
        1, 
            (Pi / 2) * Spacing * Exp[-(Pi / 4) * Spacing^2],
        2, 
            (32 / Pi^2) * Spacing^2 * Exp[-(4 / Pi) * Spacing^2],
        4, 
            (262144 / (729 * Pi^3)) * Spacing^4 * Exp[-(64 / (9 * Pi)) * Spacing^2],
        _, 
            Message[LevelSpacingDistribution::invBeta, DysonIndex]; $Failed
    ];

LevelSpacingDistribution::invBeta = "Dyson index \[Beta]=`1` is not supported. \
Valid values are 0 (Poisson), 1 (GOE), 2 (GUE), and 4 (GSE).";


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
    Switch[OptionValue[SymmetricSubspace],BosonEscapeKrausOperators
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


(* ::Subsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*Private functions*)


SwappedFockState[fockState_, targetSite_, wrongSite_] := 
  ReplacePart[fockState, 
    Thread[# -> Part[fockState, Reverse[#]] ] & [{targetSite, wrongSite}]
  ]


End[]


EndPackage[]
