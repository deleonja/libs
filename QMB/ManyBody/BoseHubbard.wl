(* ::Package:: *)

BeginPackage["QMB`ManyBody`", {"QMB`"}];


Get["ForScience`"];


(* ::Section:: *)
(*Public definitions*)


(* ::Subsection:: *)
(*Bose Hubbard with periodic boundary conditions *)


Quiet[
PeriodicBoseHubbard::usage = FormatUsage[
"PeriodicBoseHubbard[N,L,J,U] returns the Hamiltonian of a 1D Bose Hubbard chain with \
	N sites and L bosons wtih periodic boundary conditions. J and U are the hopping \
	and interactions parameters, respectively. "]
];

Quiet[
MomentumSectorPeriodicBoseHubbard::usage = FormatUsage[
"SectorPeriodicBoseHubbard[N,L,J,U,sector:All or interger] returns the Hamiltonian of \
	a 1D Bose Hubbard chain with N sites and L bosons wtih periodic boundary conditions \
	in block diagonal form due to Momentum symmetry. Using the sector=All returns the \
	full sectorized hamiltonian whereas sector= i (i <= N) returns the i-th momentum sector."
	]
];


(* ::Section:: *)
(*Private definitions*)


Begin["`BoseHubbard`Private`"];


(* ::Subsection:: *)
(*Bose Hubbard with periodic boundary conditions*)


(*Secondary*)

listaBase[v_,n_]:=Table[If[i==1,v,0],{i,n}]; 
Repetir[f_,ini_,n_]:=NestList[f,ini,n-1];
Indice[lista_]:=Module[{n=Length[lista]},SelectFirst[Range[n-1],AllTrue[lista[[#+1;;n-1]],(#==0)&]&,Missing["NotFound"]]];
BaseN[lista_,n_]:=Module[{k=Indice[lista],copia=lista},
copia[[k]]-=1;
If[k+1<=Length[copia],copia[[k+1]]=n-Total[copia[[1;;k]]];];
Do[copia[[i]]=0,{i,k+2,Length[copia]}];
copia];
Basis[N_,M_]:=Module[{Base1=listaBase[N,M],n1=(M+N-1)!/(N!*(M-1)!)},Repetir[BaseN[#,N]&,Base1,n1]//Sort];

PolarRep[num_] := Module[
  {z, T},
  z = N[num]; (* fuerza forma num\[EAcute]rica *)
  T = ToPolarCoordinates[{Re[z], Im[z]}];
  Chop[T[[1]] Exp[I T[[2]]]]
];

ShiftOperator[n_,m_]:=Module[{Fock=Basis[n,m],IntFock,S},
IntFock=AssociationThread[Fock,Range[Fock//Length]];
S=RotateRight[#]&/@Fock;
{IntFock[S[[#]]],IntFock[Fock[[#]]]}->1&/@Range[Fock//Length]//SparseArray];

OperatorAij[i_,j_,lista_]:=Module[{copy=lista,valorI,valorII},If[i==j,Return[Null]];
valorI=copy[[i]];
valorII=copy[[j]];
copy[[i]]=valorI+1;
copy[[j]]=valorII-1;
If[Min[copy]<0,Return[Null]];
{copy,Sqrt[(valorI+1.)*valorII]}];

UnitaryAtoB[EigenvecA_,EigenvecB_]:=Module[{eigenvecB=EigenvecB,eigenvecA=EigenvecA,FockToIntH,FockToIntS},FockToIntH=AssociationThread[Range[eigenvecA//Length],eigenvecA];FockToIntS=AssociationThread[Range[eigenvecB//Length],eigenvecB];
Table[{#,i}->Chop[Dot[FockToIntH[#],FockToIntS[i]]]&/@Range[FockToIntH//Length],{i,FockToIntS//Length}]//Flatten//SparseArray];

DegenerateSpectrum[S_]:=Module[{Polars,Vecs,Eigens,EigenInfoS1,InfoS,D1},
EigenInfoS1=S//Eigensystem;
Eigens=EigenInfoS1[[1]];
Vecs=EigenInfoS1[[2]];
Polars=PolarRep[#]&/@Eigens;
InfoS={Polars[[#]],Vecs[[#]]}&/@Range[Length[Vecs]];
D1=GatherBy[InfoS,First];
Table[Join[{D1[[i,1,1]]//Chop},{D1[[i,#,2]]&/@Range[D1[[i]]//Length]}],{i,Length[D1]}]];

(*Main*)
MomentumSectorPeriodicBoseHubbard[n_,m_,J_,U_,sector_:All]:=Module[{S,base,DataS,vecs,U2,H},H=PeriodicBH[n,m,J,U];
S=ShiftOperator[n,m];
base=IdentityMatrix[Length[S]];
DataS=DegenerateSpectrum[S];
vecs=Normalize/@Which[sector===All,Flatten[Last/@DataS,1],IntegerQ[sector]&&1<=sector<=Length[DataS],DataS[[sector,2]],True,(Message[SectorPeriodicBH::badsector,sector];Return[$Failed])];
U2=UnitaryAtoB[base,vecs];
Chop[ConjugateTranspose[U2] . H . U2]];

PeriodicBoseHubbard[n_,m_,J2_,U_]:=Module[{D1,Fock=Basis[n,n]//Sort,IntFock,Pares,Data,J1,len,diag,H2,Kin},
len=Fock//Length;
IntFock=AssociationThread[Fock,Range[Fock//Length]];
Pares=Table[{Mod[i+1,Fock[[1]]//Length,1],Mod[i,Fock[[1]]//Length,1]},{i,1,Fock[[1]]//Length}];
Data=Table[D1=DeleteCases[OperatorAij[#[[1]],#[[2]],Fock[[i]]]&/@Pares,Null];{IntFock[D1[[#,1]]],IntFock[Fock[[i]]]}->D1[[#,2]]&/@Range[Length[DeleteCases[OperatorAij[#[[1]],#[[2]],Fock[[i]]]&/@Pares,Null]]],{i,Fock//Length}]//Flatten;
J1=Data//SparseArray;
diag=(U/2)*(Total[#*(#-1)]&/@Fock);
H2=SparseArray[Band[{1,1}]->diag,{len,len}];
Kin=-J2(J1+ConjugateTranspose[J1]);
Kin+H2];


End[];


EndPackage[];
