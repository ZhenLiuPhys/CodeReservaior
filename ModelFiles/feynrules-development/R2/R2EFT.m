(* ::Package:: *)

(* ::Section:: *)
(*Some printout*)


NLO$Version="0.2.1";

Print[" - R2 - "];
Print["Version: "<>NLO$Version];
Print["Authors: C. Degrande"];


(* ::Subtitle::Closed:: *)
(*Simplification rules*)


red4v::usage="replacement rules for four 4-vector with loop momentum with lorentz indices not contracted together to the squared of the scalar product 
of the loop momentum times the appropriate lorentz structure and coefficient in M$dim dimensions"
red2v::usage="replacement rules for four 2-vector with loop momentum with lorentz indices not contracted together to the scalar product 
of the loop momentum times the appropriate lorentz structure and coefficient in M$dim dimensions"


red4v[lm_]:={FV[lm,a_]FV[lm,b_]FV[lm,c_]FV[lm,d_]->1/(M$dim+2)/M$dim SP[lm,lm]^2(ME[a,b]ME[c,d]+ME[a,c]ME[b,d]+ME[a,d]ME[b,c]),
             LeviCivita[lm,a_,c_,d_]FV[lm,b_]FV[lm,e_]FV[lm,f_]->1/M$dim/(M$dim+2)* SP[lm,lm]^2(LeviCivita[b,a,c,d]ME[e,f]+
             LeviCivita[e,a,c,d]ME[b,f]+LeviCivita[f,a,c,d]ME[e,b])};
red2v[lm_]:={FV[lm,a_]FV[lm,b_]->1/M$dim* SP[lm,lm]ME[a,b],SP[lm,a_]FV[lm,b_]:>1/M$dim* FV[a,b]SP[lm,lm]/;Not[a===lm],
             FCh[a__,NonCommutative[DiracSlash[lm]],b__]FV[lm,c_]->1/M$dim* FCh[a,NonCommutative[DiracMatrix[c]],b]SP[lm,lm],
             LeviCivita[lm,a_,c_,d_]FV[lm,b_]->1/M$dim* SP[lm,lm]LeviCivita[b,a,c,d]};


(* ::Subsubtitle::Closed:: *)
(*LeviCivita*)


LeviCivita[d___,FourMomentum[Incoming,a_],FourMomentum[Internal,b_],c___]:=-LeviCivita[d,FourMomentum[Internal,b],FourMomentum[Incoming,a],c];

LeviCivita[d___,FourMomentum[Internal,a_],FourMomentum[Internal,b_],c___]:=-LeviCivita[d,FourMomentum[Internal,b],FourMomentum[Incoming,a],c]/;a>b;
LeviCivita[d___,FourMomentum[Incoming,a_],FourMomentum[Incoming,b_],c___]:=-LeviCivita[d,FourMomentum[Incoming,b],FourMomentum[Incoming,a],c]/;a>b;

LeviCivita[c___,Index[Lorentz,d_],a_FourMomentum,b___]:=-LeviCivita[c,a,Index[Lorentz,d],b];

LeviCivita[c___,Index[Lorentz,d_],Index[Lorentz,a_],b___]:=-LeviCivita[c,Index[Lorentz,a],Index[Lorentz,d],b]/;d>a;


(* ::Subsubtitle::Closed:: *)
(*Metric*)


ME::usage="Metric tensor in M$dim"


Attributes[ME]=Orderless;
ME/:ME[Index[Lorentz,a_],Index[Lorentz,b_]]*ME[Index[Lorentz,a_],Index[Lorentz,c_]]:=ME[Index[Lorentz,c],Index[Lorentz,b]];
ME/:ME[Index[Lorentz,a_],Index[Lorentz,b_]]^2:=M$dim;
ME/:ME[Index[Lorentz,a_],Index[Lorentz,b_]]*f_[c___,Index[Lorentz,a_],d___]:=f[c,Index[Lorentz,b],d]/;Not[f===PolarizationVector];
ME/:ME[Index[Lorentz,a_],Index[Lorentz,b_]]*FCh[c___,NonCommutative[DiracMatrix[Index[Lorentz,a_]]],d___]:=FCh[c,NonCommutative[DiracMatrix[Index[Lorentz,b]]],d];
ME/:ME[Index[Lorentz,a_],Index[Lorentz,b_]]*ME[Index[Lorentz,a_],Index[Lorentz,b_]]:=M$dim;
ME[Index[Lorentz,a_],Index[Lorentz,a_]]:=M$dim;


(* ::Subsubtitle::Closed:: *)
(*Four vector*)


FV::usage="four vector"


FV[-a_,b_]:=-FV[a,b];
FV[a_+b_,c_]:=FV[a,c]+FV[b,c];
FV[a_?(FreeQ[#,FourMomentum]&)*b_,c_]:=a*FV[b,c];
FV/:FV[a_,b_]*FV[c_,b_]:=SP[a,c];
FV/:FV[a_,b_]^2:=SP[a,a];
FV/:FV[FourMomentum[a__],Index[Lorentz,b_]]*FCh[c___,NonCommutative[DiracMatrix[Index[Lorentz,b_]]],d___]:=FCh[c,NonCommutative[DiracSlash[FourMomentum[a]]],d];


(* ::Subsubtitle::Closed:: *)
(*Dirac Trace*)


DTr::usage="Trace over gamma matices"


(*mass*)
DTr[a___,NonCommutative[DiracSlash[b__]+Mass[yy__]],c___]:=DTr[a,NonCommutative[DiracSlash[b]],c]+Mass[yy]DTr[a,c];
DTr[a___,NonCommutative[DiracSlash[b__]]+Mass[yy__],c___]:=DTr[a,NonCommutative[DiracSlash[b]],c]+Mass[yy]DTr[a,c];

DTr[a___,Mass[yy__],c___]:=Mass[yy]DTr[a,c];

(*sum and product of momentum in dirac slash*)
DTr[a___,NonCommutative[DiracSlash[d_?(FreeQ[#,FourMomentum]&)*b__]],c___]:=d*DTr[a,NonCommutative[DiracSlash[b]],c];
DTr[a___,NonCommutative[DiracSlash[d_+b_]],c___]:=DTr[a,NonCommutative[DiracSlash[b]],c]+DTr[a,NonCommutative[DiracSlash[d]],c];
DTr[a___,NonCommutative[DiracSlash[0]],c___]:=0;

(*adjacent dirac slash with same mometum*)
DTr[a___,NonCommutative[DiracSlash[b_]],NonCommutative[DiracSlash[b_]],c___]:=DTr[a,c]SP[b,b];
DTr[NonCommutative[DiracSlash[b_]],a___,NonCommutative[DiracSlash[b_]],NonCommutative[ChiralityProjector[pm_]]]:=DTr[a,NonCommutative[ChiralityProjector[-pm]]]SP[b,b];

(*sum and product*)
DTr[a___,b_+c_,d___]:=DTr[a,b,d]+DTr[a,c,d];
DTr[a___,b_*c_?(FreeQ[#,NonCommutative]&),d___]:=c*DTr[a,b,d];
DTr[a___,b_*G[k_][l_][m__][n__],d___]:=G[k][l][m][n]*DTr[a,b,d];
(*Depth of one for the feynman parameters*)
DTr[a___,c_?((FreeQ[#,NonCommutative]&&Depth[#]<2)&),d___]:=c*DTr[a,d];
DTr[a___,NonCommutative[b_,c__],d___]:=DTr[a,NonCommutative[b],NonCommutative[c],d];

(*Projector rules*)
DTr[a___,NonCommutative[ChiralityProjector[pm_]],NonCommutative[b:_DiracSlash|_DiracMatrix],d___]:=DTr[a,NonCommutative[b],NonCommutative[ChiralityProjector[-pm]],d];
DTr[a___,NonCommutative[ChiralityProjector[1]],NonCommutative[ChiralityProjector[-1]],d___]:=0;
DTr[a___,NonCommutative[ChiralityProjector[-1]],NonCommutative[ChiralityProjector[1]],d___]:=0;
DTr[a___,NonCommutative[ChiralityProjector[pm_]],NonCommutative[ChiralityProjector[pm_]],d___]:=DTr[a,NonCommutative[ChiralityProjector[pm]],d];

(*1 and 3 gamma*)
DTr[NonCommutative[a:_DiracMatrix|_DiracSlash],NonCommutative[b:_DiracMatrix|_DiracSlash],NonCommutative[c:_DiracMatrix|_DiracSlash],NonCommutative[ChiralityProjector[pm_]]]:=0;
DTr[NonCommutative[a:_DiracMatrix|_DiracSlash],NonCommutative[ChiralityProjector[pm_]]]:=0;

DTr[NonCommutative[ChiralityProjector[pm_]]]:=2;

(*gathering gamma with the same index or momentum*)
DTr[a___,NonCommutative[DiracSlash[lm_]],NonCommutative[DiracMatrix[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2FV[lm,b]*DTr[a,NonCommutative[DiracSlash[lm]],c]-
SP[lm,lm]*DTr[a,NonCommutative[DiracMatrix[b]],c];
DTr[NonCommutative[DiracSlash[lm_]],a___,NonCommutative[DiracSlash[lm_]],NonCommutative[DiracMatrix[b_]],NonCommutative[ChiralityProjector[pm_]]]:=
2FV[lm,b]*DTr[a,NonCommutative[DiracSlash[lm]],NonCommutative[ChiralityProjector[-pm]]]-SP[lm,lm]*DTr[a,NonCommutative[DiracMatrix[b]],NonCommutative[ChiralityProjector[-pm]]];
(*generic*)
DTr[a___,NonCommutative[DiracMatrix[lm_]],d___,NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[lm_]],c___]:=2ME[lm,b]*DTr[a,NonCommutative[DiracMatrix[lm]],d,c]-
DTr[a,NonCommutative[DiracMatrix[lm]],d,NonCommutative[DiracMatrix[lm]],NonCommutative[DiracMatrix[b]],c];
DTr[a___,NonCommutative[DiracSlash[lm_]],d___,NonCommutative[DiracMatrix[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2FV[lm,b]*DTr[a,NonCommutative[DiracSlash[lm]],d,c]-
DTr[a,NonCommutative[DiracSlash[lm]],d,NonCommutative[DiracSlash[lm]],NonCommutative[DiracMatrix[b]],c];
DTr[a___,NonCommutative[DiracMatrix[lm_]],d___,NonCommutative[DiracSlash[b_]],NonCommutative[DiracMatrix[lm_]],c___]:=2FV[b,lm]*DTr[a,NonCommutative[DiracMatrix[lm]],d,c]-
DTr[a,NonCommutative[DiracMatrix[lm]],d,NonCommutative[DiracMatrix[lm]],NonCommutative[DiracSlash[b]],c];
DTr[a___,NonCommutative[DiracSlash[lm_]],d___,NonCommutative[DiracSlash[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2SP[lm,b]*DTr[a,NonCommutative[DiracSlash[lm]],d,c]-
DTr[a,NonCommutative[DiracSlash[lm]],d,NonCommutative[DiracSlash[lm]],NonCommutative[DiracSlash[b]],c];


(*Only valid in 4 dimension (to do : just for external momenta and indices) *)

DTr[NonCommutative[DiracMatrix[a_]],NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(ME[a,b]ME[c,d]-ME[a,c]ME[b,d]+ME[a,d]ME[b,c]+ I pm*LeviCivita[a,b,c,d]);
DTr[NonCommutative[DiracSlash[a__]],NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(FV[a,b]ME[c,d]-FV[a,c]ME[b,d]+FV[a,d]ME[b,c]+ I pm*LeviCivita[a,b,c,d]);
DTr[NonCommutative[DiracMatrix[a_]],NonCommutative[DiracSlash[b__]],NonCommutative[DiracMatrix[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(FV[b,a]ME[c,d]-ME[a,c]FV[b,d]+ME[a,d]FV[b,c]+ I pm*LeviCivita[a,b,c,d]);
DTr[NonCommutative[DiracMatrix[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[DiracSlash[a__]],NonCommutative[DiracMatrix[b_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(FV[a,b]ME[c,d]-FV[a,c]ME[b,d]+FV[a,d]ME[b,c]+ I pm*LeviCivita[a,b,c,d]);
DTr[NonCommutative[DiracMatrix[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[DiracMatrix[a_]],NonCommutative[DiracSlash[b__]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(FV[b,a]ME[c,d]-ME[a,c]FV[b,d]+ME[a,d]FV[b,c]+ I pm*LeviCivita[a,b,c,d]);

DTr[NonCommutative[DiracSlash[a__]],NonCommutative[DiracSlash[b_]],NonCommutative[DiracMatrix[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(SP[a,b]ME[c,d]-FV[a,c]FV[b,d]+FV[a,d]FV[b,c]+ pm*I LeviCivita[a,b,c,d])/;FreeQ[a,Internal]||FreeQ[b,Internal];
DTr[NonCommutative[DiracMatrix[a_]],NonCommutative[DiracSlash[b__]],NonCommutative[DiracSlash[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(FV[b,a]FV[c,d]-FV[c,a]FV[b,d]+ME[a,d]SP[b,c]+ pm*I LeviCivita[a,b,c,d])/;FreeQ[c,Internal]||FreeQ[b,Internal];
DTr[NonCommutative[DiracMatrix[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[DiracSlash[a__]],NonCommutative[DiracSlash[b_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(SP[a,b]ME[c,d]-FV[a,c]FV[b,d]+FV[a,d]FV[b,c]+ pm*I LeviCivita[a,b,c,d])/;FreeQ[a,Internal]||FreeQ[b,Internal];
DTr[NonCommutative[DiracSlash[c_]],NonCommutative[DiracMatrix[d_]],NonCommutative[DiracMatrix[a_]],NonCommutative[DiracSlash[b__]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(FV[b,a]FV[c,d]-FV[c,a]FV[b,d]+ME[a,d]SP[b,c]+ pm*I LeviCivita[a,b,c,d])/;FreeQ[c,Internal]||FreeQ[b,Internal];

DTr[NonCommutative[DiracSlash[a__]],NonCommutative[DiracMatrix[c_]],NonCommutative[DiracSlash[b_]],NonCommutative[DiracMatrix[d_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(-SP[a,b]ME[c,d]+FV[a,c]FV[b,d]+FV[a,d]FV[b,c]- pm*I LeviCivita[a,b,c,d])/;FreeQ[a,Internal]||FreeQ[b,Internal];
DTr[NonCommutative[DiracMatrix[a_]],NonCommutative[DiracSlash[b__]],NonCommutative[DiracMatrix[d_]],NonCommutative[DiracSlash[c_]],NonCommutative[ChiralityProjector[pm_]]]:=
2*(FV[b,a]FV[c,d]+FV[c,a]FV[b,d]-ME[a,d]SP[b,c]- pm*I LeviCivita[a,b,c,d])/;FreeQ[c,Internal]||FreeQ[b,Internal];

DTr[NonCommutative[DiracMatrix[a_]],NonCommutative[DiracMatrix[b_]],NonCommutative[ChiralityProjector[pm_]]]:=2*(ME[b,a]);
DTr[NonCommutative[DiracMatrix[a_]],NonCommutative[DiracSlash[b_]],NonCommutative[ChiralityProjector[pm_]]]:=2*(FV[b,a]);
DTr[NonCommutative[DiracSlash[b_]],NonCommutative[DiracMatrix[a_]],NonCommutative[ChiralityProjector[pm_]]]:=2*(FV[b,a]);
DTr[NonCommutative[DiracSlash[a_]],NonCommutative[DiracSlash[b_]],NonCommutative[ChiralityProjector[pm_]]]:=2*SP[a,b];/;FreeQ[a,Internal]||FreeQ[b,Internal];


(* ::Subsubtitle::Closed:: *)
(*Fermion Chain*)


FCh::usage="Chain with gamma matices between two spinors"


FCh[a___,NonCommutative[DiracSlash[b__]+Mass[yy__]],c___]:=FCh[a,NonCommutative[DiracSlash[b]],c]+Mass[yy]FCh[a,c];
FCh[a___,NonCommutative[DiracSlash[b__]]+Mass[yy__],c___]:=FCh[a,NonCommutative[DiracSlash[b]],c]+Mass[yy]FCh[a,c];

FCh[a___,NonCommutative[DiracSlash[b_]],NonCommutative[DiracSlash[b_]],c___]:=FCh[a,c]SP[b,b];
FCh[a___,NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[b_]],c___]:=FCh[a,c]M$dim;

FCh[a___,NonCommutative[DiracSlash[d_?(FreeQ[#,FourMomentum]&)*b__]],c___]:=d*FCh[a,NonCommutative[DiracSlash[b]],c];
FCh[a___,NonCommutative[DiracSlash[d_+b_]],c___]:=FCh[a,NonCommutative[DiracSlash[b]],c]+FCh[a,NonCommutative[DiracSlash[d]],c];

FCh[a___,NonCommutative[DiracSlash[-b_]],c___]:=-FCh[a,NonCommutative[DiracSlash[b]],c];
FCh[a___,NonCommutative[DiracSlash[0]],c___]:=0;

FCh[a___,b_+c_,d___]:=FCh[a,b,d]+FCh[a,c,d];
FCh[a___,b_*c_?(FreeQ[#,NonCommutative]&),d___]:=c*FCh[a,b,d];
FCh[a___,b_*G[k_][l_][m__][n__],d___]:=G[k][l][m][n]*FCh[a,b,d];
FCh[a___,c_?((FreeQ[#,NonCommutative]&&Depth[#]<2)&),d___]:=c*FCh[a,d];

FCh[a___,NonCommutative[b_,c__],d___]:=FCh[a,NonCommutative[b],NonCommutative[c],d];

FCh[a___,NonCommutative[ChiralityProjector[pm_]],NonCommutative[b:_DiracSlash|_DiracMatrix],d___]:=FCh[a,NonCommutative[b],NonCommutative[ChiralityProjector[-pm]],d];
FCh[a___,NonCommutative[ChiralityProjector[1]],NonCommutative[ChiralityProjector[-1]],d___]:=0;
FCh[a___,NonCommutative[ChiralityProjector[-1]],NonCommutative[ChiralityProjector[1]],d___]:=0;
FCh[a___,NonCommutative[ChiralityProjector[pm_]],NonCommutative[ChiralityProjector[pm_]],d___]:=FCh[a,NonCommutative[ChiralityProjector[pm]],d];

FCh[a___,NonCommutative[DiracSlash[lm_]],NonCommutative[DiracMatrix[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2FV[lm,b]*FCh[a,NonCommutative[DiracSlash[lm]],c]-
SP[lm,lm]*FCh[a,NonCommutative[DiracMatrix[b]],c];

FCh[a___,NonCommutative[DiracSlash[lm_]],NonCommutative[DiracSlash[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2SP[lm,b]*FCh[a,NonCommutative[DiracSlash[lm]],c]-
SP[lm,lm]*FCh[a,NonCommutative[DiracSlash[b]],c];

FCh[a___,NonCommutative[DiracMatrix[l_]],NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[l_]],c___]:=(2-M$dim)*FCh[a,NonCommutative[DiracMatrix[b]],c];

FCh[a___,NonCommutative[DiracMatrix[l_]],NonCommutative[DiracSlash[b_]],NonCommutative[DiracMatrix[l_]],c___]:=(2-M$dim)*FCh[a,NonCommutative[DiracSlash[b]],c];

(*generic case*)

FCh[a___,NonCommutative[DiracSlash[lm_]],d__,NonCommutative[DiracMatrix[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2FV[lm,b]*FCh[a,NonCommutative[DiracSlash[lm]],d,c]-
FCh[a,NonCommutative[DiracSlash[lm]],d,NonCommutative[DiracSlash[lm]],NonCommutative[DiracMatrix[b]],c];

FCh[a___,NonCommutative[DiracSlash[lm_]],d__,NonCommutative[DiracSlash[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2SP[lm,b]*FCh[a,NonCommutative[DiracSlash[lm]],d,c]-
FCh[a,NonCommutative[DiracSlash[lm]],d,NonCommutative[DiracSlash[lm]],NonCommutative[DiracSlash[b]],c];

FCh[a___,NonCommutative[DiracMatrix[l_]],d__,NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[l_]],c___]:=2ME[b,l]*FCh[a,NonCommutative[DiracMatrix[l]],d,c]-
FCh[a,NonCommutative[DiracMatrix[l]],d,NonCommutative[DiracMatrix[l]],NonCommutative[DiracMatrix[b]],c];

FCh[a___,NonCommutative[DiracMatrix[l_]],d__,NonCommutative[DiracSlash[b_]],NonCommutative[DiracMatrix[l_]],c___]:=2FV[b,l]*FCh[a,NonCommutative[DiracMatrix[l]],d,c]-
FCh[a,NonCommutative[DiracMatrix[l]],d,NonCommutative[DiracMatrix[l]],NonCommutative[DiracSlash[b]],c];


(* ::Subsubtitle::Closed:: *)
(*Distributive Fermion Chain*)


FChDist::usage="Chain with gamma matices between two spinors, only distribution of the sum, no gamma algebra"


(*FCh[a___,NonCommutative[DiracSlash[b__]+Mass[yy__]],c___]:=FCh[a,NonCommutative[DiracSlash[b]],c]+Mass[yy]FCh[a,c];
FCh[a___,NonCommutative[DiracSlash[b__]]+Mass[yy__],c___]:=FCh[a,NonCommutative[DiracSlash[b]],c]+Mass[yy]FCh[a,c];

FCh[a___,NonCommutative[DiracSlash[b_]],NonCommutative[DiracSlash[b_]],c___]:=FCh[a,c]SP[b,b];
FCh[a___,NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[b_]],c___]:=FCh[a,c]M$dim;

FCh[a___,NonCommutative[DiracSlash[d_?(FreeQ[#,FourMomentum]&)*b__]],c___]:=d*FCh[a,NonCommutative[DiracSlash[b]],c];
FCh[a___,NonCommutative[DiracSlash[d_+b_]],c___]:=FCh[a,NonCommutative[DiracSlash[b]],c]+FCh[a,NonCommutative[DiracSlash[d]],c];

FCh[a___,NonCommutative[DiracSlash[-b_]],c___]:=-FCh[a,NonCommutative[DiracSlash[b]],c];
FCh[a___,NonCommutative[DiracSlash[0]],c___]:=0;*)

FChDist[a___,b_+c_,d___]:=FChDist[a,b,d]+FChDist[a,c,d];
FChDist[a___,b_*c_?(FreeQ[#,NonCommutative]&),d___]:=c*FChDist[a,b,d];
FChDist[a___,b_*G[k_][l_][m__][n__],d___]:=G[k][l][m][n]*FChDist[a,b,d];
(*FCh[a___,c_?((FreeQ[#,NonCommutative]&&Depth[#]<2)&),d___]:=c*FCh[a,d];*)

(*FCh[a___,NonCommutative[b_,c__],d___]:=FCh[a,NonCommutative[b],NonCommutative[c],d];

FCh[a___,NonCommutative[ChiralityProjector[pm_]],NonCommutative[b:_DiracSlash|_DiracMatrix],d___]:=FCh[a,NonCommutative[b],NonCommutative[ChiralityProjector[-pm]],d];
FCh[a___,NonCommutative[ChiralityProjector[1]],NonCommutative[ChiralityProjector[-1]],d___]:=0;
FCh[a___,NonCommutative[ChiralityProjector[-1]],NonCommutative[ChiralityProjector[1]],d___]:=0;
FCh[a___,NonCommutative[ChiralityProjector[pm_]],NonCommutative[ChiralityProjector[pm_]],d___]:=FCh[a,NonCommutative[ChiralityProjector[pm]],d];

FCh[a___,NonCommutative[DiracSlash[lm_]],NonCommutative[DiracMatrix[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2FV[lm,b]*FCh[a,NonCommutative[DiracSlash[lm]],c]-
SP[lm,lm]*FCh[a,NonCommutative[DiracMatrix[b]],c];

FCh[a___,NonCommutative[DiracSlash[lm_]],NonCommutative[DiracSlash[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2SP[lm,b]*FCh[a,NonCommutative[DiracSlash[lm]],c]-
SP[lm,lm]*FCh[a,NonCommutative[DiracSlash[b]],c];

FCh[a___,NonCommutative[DiracMatrix[l_]],NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[l_]],c___]:=(2-M$dim)*FCh[a,NonCommutative[DiracMatrix[b]],c];

FCh[a___,NonCommutative[DiracMatrix[l_]],NonCommutative[DiracSlash[b_]],NonCommutative[DiracMatrix[l_]],c___]:=(2-M$dim)*FCh[a,NonCommutative[DiracSlash[b]],c];

(*generic case*)

FCh[a___,NonCommutative[DiracSlash[lm_]],d__,NonCommutative[DiracMatrix[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2FV[lm,b]*FCh[a,NonCommutative[DiracSlash[lm]],d,c]-
FCh[a,NonCommutative[DiracSlash[lm]],d,NonCommutative[DiracSlash[lm]],NonCommutative[DiracMatrix[b]],c];

FCh[a___,NonCommutative[DiracSlash[lm_]],d__,NonCommutative[DiracSlash[b_]],NonCommutative[DiracSlash[lm_]],c___]:=2SP[lm,b]*FCh[a,NonCommutative[DiracSlash[lm]],c]-
FCh[a,NonCommutative[DiracSlash[lm]],d,NonCommutative[DiracSlash[lm]],NonCommutative[DiracSlash[b]],c];

FCh[a___,NonCommutative[DiracMatrix[l_]],d__,NonCommutative[DiracMatrix[b_]],NonCommutative[DiracMatrix[l_]],c___]:=2ME[b,l]*FCh[a,NonCommutative[DiracMatrix[l]],d,c]-
FCh[a,NonCommutative[DiracMatrix[l]],d,NonCommutative[DiracMatrix[l]],NonCommutative[DiracMatrix[b]],c];

FCh[a___,NonCommutative[DiracMatrix[l_]],d__,NonCommutative[DiracSlash[b_]],NonCommutative[DiracMatrix[l_]],c___]:=2FV[b,l]*FCh[a,NonCommutative[DiracMatrix[l]],d,c]-
FCh[a,NonCommutative[DiracMatrix[l]],d,NonCommutative[DiracMatrix[l]],NonCommutative[DiracSlash[b]],c];*)


(* ::Subsubtitle::Closed:: *)
(*ScalarProduct*)


SP::usage="scalar product"


Attributes[SP]=Orderless;
SP[-a_,b_]:=-SP[a,b];
SP[-a_,a_]:=-SP[a,a];
SP[-a_,-a_]:=SP[a,a];
SP[a_,b_+c_]:=SP[a,b]+SP[a,c];
SP[a_,b_*FourMomentum[c__]]:=b*SP[a,FourMomentum[c]];
SP[a_,b_?(FreeQ[#,FourMomentum]&)*c_]:=b*SP[a,c];


(* ::Subtitle:: *)
(*UV Tools*)


(* ::Subsubtitle:: *)
(*AddWf*)


AddWf[vert_]:=Block[{partlist,incr,a,wflist,pos},
  partlist=vert[[1,All,1]]/.{V[a_,__]->V[a],F[a_,__]->F[a],S[a_,__]->S[a]}/.{-V[a_]->V[a],-S[a_]->S[a],-F[a_]->F[a]}; 
  wflist=UV$Wftlist[[All,1]]/.{V[a_,__]->V[a],F[a_,__]->F[a],S[a_,__]->S[a]}/.{-V[a_]->V[a],-S[a_]->S[a],-F[a_]->F[a]};
  pos=(If[FreeQ[wflist,#],0,Position[wflist,#][[1,1]]]&)/@partlist;
  pos = DeleteCases[pos,0];
  (*remove part if interaction comes from the kinetic term like VVV, VVVV, FFV, SSV, SSVV since all vector are gauge boson*)
  {vert[[1]],vert[[2]],-1/2*vert[[3]]*(Total[UV$Wftlist[[pos,2]]])}
];


(* ::Subsubtitle::Closed:: *)
(*MergeVertList*)


MergeVertList[vl1_,vl2_]:=Block[{vsh,vlg,pos,verk},
  vsh=If[Length[vl1]>Length[vl2],vl2,vl1];
  vlg=If[Length[vl1]>Length[vl2],vl1,vl2];
  For[verk=1,verk<=Length[vsh],verk++,
    If[Not[FreeQ[vlg,vsh[[verk,1]]]],
      pos=Position[vlg[[All,1]],vsh[[verk,1]]][[1]];
      vlg[[pos,2]]=vlg[[pos,2]]+vsh[[verk,2]];
      vlg[[pos,3]]=vlg[[pos,3]]+vsh[[verk,3]];,
      vlg=Append[vlg,vsh[[verk]]];
    ];
  ];
  vlg
];


(* ::Subsubtitle:: *)
(*XIntegrate*)


XIntegrate[num_, del_,x_]:=Block[{m1,m2,p2,c0,c1,c2,nolog,mext,withlog},

m1 = (del/.{x->1})[[1]]/FR$MU;
m2 = (del/.{x->0})[[1]]/FR$MU;
p2 = Coefficient[del,x^2]/FR$MU^2;
nolog = Integrate[Total[Cases[Expand[num],_?(FreeQ[#,Log]&)]],{x,0,1}];
withlog = Total[Cases[Expand[num],_?(Not[FreeQ[#,Log]]&)]]/.Log[__]->1;

nolog+Sum[IntxnLog[m1,m2,p2,nx]Coefficient[withlog,x,nx],{nx,0,2}]
];


IntxnLog::usage="IntxnLog[m1,m2,p2,n] is the integral over x of x^n Log[m2^2+x(m1^2-m2^2)+p2 x(x-1)], m1>=0,m2>=0,p2>=0,n=0,1 or 2"


IntxnLog[0,0,0,n_Integer]:=-2/FR$Eps*FR$IR/(n+1);

IntxnLog[0,0,p2_,0]:=-2 + If[p2\[Element]Reals,0,1]*Re[Log[-p2]] + If[p2\[Element]Reals,1,0]Re[(-I*Pi + Log[p2])]/;FreeQ[p2,FourMomentum];
IntxnLog[0,0,p2_,1]:=(-2 +  (If[p2\[Element]Reals,0,1]*Re[Log[-p2]] + If[p2\[Element]Reals,1,0]Re[(-I*Pi + Log[p2])]))/2/;FreeQ[p2,FourMomentum];
IntxnLog[0,0,p2_,2]:=(-13 + 6*(If[p2\[Element]Reals,0,1]*Re[Log[-p2]] + If[p2\[Element]Reals,1,0]Re[(-I*Pi + Log[p2])]))/18/;FreeQ[p2,FourMomentum];

IntxnLog[0,m2_,0,0]:=-1 + 2*Re[Log[m2]]/;FreeQ[m2,Mass];
IntxnLog[0,m2_,0,1]:=-3/4 + Re[Log[m2]]/;FreeQ[m2,Mass];
IntxnLog[0,m2_,0,2]:=-11/18 + (2*Re[Log[m2]])/3/;FreeQ[m2,Mass];

IntxnLog[m1_,0,0,0]:=-1 + 2*Re[Log[m1]]/;FreeQ[m1,Mass];
IntxnLog[m1_,0,0,1]:=-1/4 + Re[Log[m1]]/;FreeQ[m1,Mass];
IntxnLog[m1_,0,0,2]:=(-1 + 6*Re[Log[m1]])/9/;FreeQ[m1,Mass];

IntxnLog[m_,m_,0,0]:=Re[Log[m^2]]/;FreeQ[m,Mass];
IntxnLog[m_,m_,0,1]:=Re[Log[m]]/;FreeQ[m,Mass];
IntxnLog[m_,m_,0,2]:=(2*Re[Log[m]])/3/;FreeQ[m,Mass];

IntxnLog[m_/FR$MU,0,m_^2/FR$MU^2,0]:=2*(-1 + Re[Log[m/FR$MU]])/;FreeQ[m,Mass];
IntxnLog[m_/FR$MU,0,m_^2/FR$MU^2,1]:=-1/2 + Re[Log[m/FR$MU]]/;FreeQ[m,Mass];
IntxnLog[m_/FR$MU,0,m_^2/FR$MU^2,2]:=(2*(-1 + 3*Re[Log[m/FR$MU]]))/9/;FreeQ[m,Mass];

IntxnLog[0,m_/FR$MU,m_^2/FR$MU^2,0]:=2*(-1 + Re[Log[m/FR$MU]])/;FreeQ[m,Mass];
IntxnLog[0,m_/FR$MU,m_^2/FR$MU^2,1]:=-3/2 + Re[Log[m/FR$MU]]/;FreeQ[m,Mass];
IntxnLog[0,m_/FR$MU,m_^2/FR$MU^2,2]:=(-11 + 6*Re[Log[m/FR$MU]])/9/;FreeQ[m,Mass];

IntxnLog[m1_,m2_,0,0]:=Re[(-m1^2 + m2^2 + 2*m1^2*Log[m1] - 2*m2^2*Log[m2])/(m1^2 - m2^2)]/;Not[m1===m2]&&FreeQ[m1,Mass]&&FreeQ[m2,Mass];
IntxnLog[m1_,m2_,0,1]:=Re[(-m1^4 + 4*m1^2*m2^2 - 3*m2^4 + 4*(m1^4 - 2*m1^2*m2^2)*Log[m1] + 4*m2^4*Log[m2])/(4*(m1^2 - m2^2)^2)]/;Not[m1===m2]&&FreeQ[m1,Mass]&&FreeQ[m2,Mass];
IntxnLog[m1_,m2_,0,2]:=Re[(-2*m1^6 + 9*m1^4*m2^2 - 18*m1^2*m2^4 + 11*m2^6 + 12*(m1^6 - 3*m1^4*m2^2 + 3*m1^2*m2^4)*Log[m1] - 12*m2^6*Log[m2])/(18*(m1^2 - m2^2)^3)]/;Not[m1===m2]&&FreeQ[m1,Mass]&&FreeQ[m2,Mass];

IntxnLog[m1_,0,p2_,0]:=Re[(-2*p2 + 2*m1^2*Log[m1] + (-m1^2 + p2)*(If[m1^2<p2&&p2\[Element]Reals,1,0]*(Log[-m1^2 + p2]+I Pi)+If[m1^2<p2&&p2\[Element]Reals,0,1]*Log[m1^2 - p2]))/p2]/;FreeQ[p2,FourMomentum];
IntxnLog[m1_,0,p2_,1]:=Re[1/2/p2(m1^2(2Log[m1]-1)-(IntxnLog[m1,0,p2,0]/.Re->Identity)(m1^2-p2))]/;FreeQ[m1,Mass]&&FreeQ[p2,FourMomentum];
IntxnLog[m1_,0,p2_,2]:=Re[((m1^2-p2/3)/2-2(m1^2-p2)(IntxnLog[m1,0,p2,1]/.Re->Identity)-m1^2(1-2Log[m1]))/3/p2]/;FreeQ[m1,Mass]&&FreeQ[p2,FourMomentum];

IntxnLog[0,m2_,p2_,0]:=Re[(-2*p2 + 2*m2^2*Log[m2] + (-m2^2 + p2)*(If[m2^2<p2&&p2\[Element]Reals,1,0]*(Log[-m2^2 + p2]+I Pi)+If[m2^2<p2&&p2\[Element]Reals,0,1]*Log[m2^2 - p2]))/p2]/;FreeQ[p2,FourMomentum];
IntxnLog[0,m2_,p2_,1]:=Re[1/2/p2(-m2^2(2Log[m2]-1)-(IntxnLog[0,m2,p2,0]/.Re->Identity)(-m2^2-p2))]/;FreeQ[m2,Mass]&&FreeQ[p2,FourMomentum];
IntxnLog[0,m2_,p2_,2]:=Re[((m2^2-p2/3)/2-m2^2 (IntxnLog[0,m2,p2,0]/.Re->Identity)-2(-m2^2-p2)(IntxnLog[0,m2,p2,1]/.Re->Identity))/3/p2]/;FreeQ[m2,Mass]&&FreeQ[p2,FourMomentum];

IntxnLog[m_/FR$MU,m_/FR$MU,m_^2/FR$MU^2,0]:=-2 + Pi/Sqrt[3] + 2*Re[Log[m/FR$MU]]/;FreeQ[m,Mass];
IntxnLog[m_/FR$MU,m_/FR$MU,m_^2/FR$MU^2,1]:=-1 + Pi/(2*Sqrt[3]) + Re[Log[m/FR$MU]]/;FreeQ[m,Mass];
IntxnLog[m_/FR$MU,m_/FR$MU,m_^2/FR$MU^2,2]:=-1/18 + (2*Re[Log[m/FR$MU]])/3/;FreeQ[m,Mass];

IntxnLog[m_,m_,p2_,0]:=Re[-2 + Log[m^2] - (Sqrt[p2*(-4*m^2 + p2)]*(If[m\[Element]Reals&&p2\[Element]Reals&&m>0&&p2>4m^2,0,1]*Log[(2*m^2 - p2 + Sqrt[p2*(-4*m^2 + p2)])/(2*m^2)]+ 
  If[m\[Element]Reals&&p2\[Element]Reals&&m>0&&p2>4m^2,1,0](Log[-(2*m^2 - p2 + Sqrt[p2*(-4*m^2 + p2)])/(2*m^2)]+I Pi)))/p2]/;FreeQ[p2,FourMomentum];

IntxnLog[m1_/FR$MU,m2_/FR$MU,m1_^2/FR$MU^2,0]:= Re[-2 + (-1 + m2^2/m1^2)*Log[m2/m1] + Log[(m1*m2)/FR$MU^2] - (Sqrt[-4*m1^2*m2^2 + m2^4]*
  Log[(m2^2 + Sqrt[-4*m1^2*m2^2 + m2^4])/(2*m1*m2)])/m1^2]/;FreeQ[m1,Mass]&&FreeQ[m2,Mass];

IntxnLog[m1_/FR$MU,m2_/FR$MU,m2_^2/FR$MU^2,0]:=IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,0];
IntxnLog[m1_/FR$MU,m2_/FR$MU,m2_^2/FR$MU^2,1]:=IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,0]-IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,1];
IntxnLog[m1_/FR$MU,m2_/FR$MU,m2_^2/FR$MU^2,2]:=IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,0]-2 IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,1]+IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,2];

IntxnLog[m1_,m2_,p2_,0]:=Re[-2 - ((m1^2 - m2^2)*Log[m2/m1])/p2 + Log[(m1*m2)/FR$MU^2] + 
 (m1*m2*(-(m1^2 + m2^2 + Sqrt[-4*m1^2*m2^2 + (m1^2 + m2^2 - p2)^2] - p2)/(2*m1*m2) - (-m1^2 - m2^2 + Sqrt[-4*m1^2*m2^2 + (m1^2 + m2^2 - p2)^2] + p2)/(2*m1*m2))*
   (If[p2>=(m1+m2)^2&&p2\[Element]Reals&&m1\[Element]Reals&&m2\[Element]Reals,0,1]Log[(m1^2 + m2^2 + Sqrt[-4*m1^2*m2^2 + (m1^2 + m2^2 - p2)^2] - p2)/(2*m1*m2)]+
    If[p2>=(m1+m2)^2&&p2\[Element]Reals&&m1\[Element]Reals&&m2\[Element]Reals,1,0](Log[(m1^2 + m2^2 + Sqrt[-4*m1^2*m2^2 + (m1^2 + m2^2 - p2)^2] - p2)/(2*m1*m2)]+I Pi)))/p2]/;FreeQ[m1,Mass]&&FreeQ[m2,Mass]&&FreeQ[p2,FourMomentum];

IntxnLog[m1_,m2_,p2_,1]:=Re[1/2/p2(-m2^2(2Log[m2]-1)+m1^2(2Log[m1]-1)-(IntxnLog[m1,m2,p2,0]/.Re->Identity)(-m2^2+m1^2-p2))]/;FreeQ[m1,Mass]&&FreeQ[m2,Mass]&&
FreeQ[p2,FourMomentum];

IntxnLog[m1_,m2_,p2_,2]:=Re[((m1^2+m2^2-p2/3)/2-m2^2 (IntxnLog[m1,m2,p2,0]/.Re->Identity)-2(-m2^2+m1^2-p2)(IntxnLog[m1,m2,p2,1]/.Re->Identity)-m1^2(1-2Log[m1]))/3/p2]/;
FreeQ[m1,Mass]&&FreeQ[m2,Mass]&&FreeQ[p2,FourMomentum];


UVLog[0]:=0;(*from the tadpole from x Log[x], x->0*)
UVLog[x_?(Not[Head[#]===Mass]&)^2/FR$MU^2]:=Re[Log[x^2/FR$MU^2]];


(* ::Subsubtitle:: *)
(*Dp2IntxnLog(to finish)*)


Dp2IntxnLog::usage="Dp2IntxnLog[m1,m2,p2,n] is the derivative with respect to p2 of the real part of the integral over x of x^n Log[m2^2+x(m1^2-m2^2)+p2x(x-1)], 
m1>=0,m2>=0,p2>=0,n=0,1 or 2"


(*warning d=4+eps almost until the end*)

Dp2IntxnLog[0,0,0,n_]:=0;

Dp2IntxnLog[0,0,p2_,0]:=Re[1/p2]/;FreeQ[p2,FourMomentum];
Dp2IntxnLog[0,0,p2_,1]:=Re[1/2/p2]/;FreeQ[p2,FourMomentum];
Dp2IntxnLog[0,0,p2_,2]:=Re[1/3/p2]/;FreeQ[p2,FourMomentum];

Dp2IntxnLog[m1_/FR$MU,0,m1_^2/FR$MU^2,0]:=-1(1/FR$Eps*FR$IR/Re[m1^2]+Re[(-1+Log[m1/FR$MU])/m1^2])FR$MU^2/;FreeQ[m1,FourMomentum];(*Infared divergent*)
Dp2IntxnLog[m1_/FR$MU,0,m1_^2/FR$MU^2,1]:=-1/2/Re[m1^2]*FR$MU^2/;FreeQ[m1,FourMomentum];
Dp2IntxnLog[m1_/FR$MU,0,m1_^2/FR$MU^2,2]:=-1/Re[m1^2]/6*FR$MU^2/;FreeQ[m1,FourMomentum];

Dp2IntxnLog[0,m2_/FR$MU,m2_^2/FR$MU^2,0]:=-(1/FR$Eps*FR$IR/Re[m2^2]+Re[(-1+Log[m2/FR$MU])/m2^2])*FR$MU^2/;FreeQ[m2,FourMomentum];
Dp2IntxnLog[0,m2_/FR$MU,m2_^2/FR$MU^2,1]:=-(1/FR$Eps*FR$IR/Re[m2^2]+Re[(-3/2+Log[m2/FR$MU])/m2^2])*FR$MU^2/;FreeQ[m2,FourMomentum];
Dp2IntxnLog[0,m2_/FR$MU,m2_^2/FR$MU^2,2]:=-(1/FR$Eps*FR$IR/Re[m2^2]+Re[(-11/6+Log[m2/FR$MU])/m2^2])*FR$MU^2/;FreeQ[m2,FourMomentum];

Dp2IntxnLog[m1_,0,p2_,0]:=Re[(p2 - 2*m1^2*Log[m1] + m1^2*(If[m1^2<p2&&p2\[Element]Reals,1,0]*(Log[-m1^2 + p2]+I Pi)+If[m1^2<p2&&p2\[Element]Reals,0,1]*Log[m1^2 - p2]))/p2^2]/;FreeQ[p2,FourMomentum];
Dp2IntxnLog[m1_,0,p2_,1]:=Re[-(IntxnLog[m1,0,p2,1]/.Re->Identity)/p2+(IntxnLog[m1,0,p2,0]-(m1^2-p2)Dp2IntxnLog[m1,0,p2,0]/.Re->Identity)/2/p2]/;FreeQ[m1,Mass]&&
FreeQ[p2,FourMomentum];
Dp2IntxnLog[m1_,0,p2_,2]:=Re[-(IntxnLog[m1,0,p2,2]/Re->Identity)/p2+(-1/6+2IntxnLog[m1,0,p2,1]-2(m1^2-p2)Dp2IntxnLog[m1,0,p2,1]/.Re->Identity)/3/p2]/;
FreeQ[m1,Mass]&&FreeQ[p2,FourMomentum];

Dp2IntxnLog[0,m2_,p2_,0]:=Dp2IntxnLog[m2,0,p2,0]/;FreeQ[p2,FourMomentum];
Dp2IntxnLog[0,m2_,p2_,1]:=Re[- IntxnLog[0,m2,p2,1]/p2+(IntxnLog[0,m2,p2,0]-(-m2^2-p2)Dp2IntxnLog[0,m2,p2,0])/2/p2/.Re->Identity]/;FreeQ[m2,Mass]&&
FreeQ[p2,FourMomentum];
Dp2IntxnLog[0,m2_,p2_,2]:=Re[-IntxnLog[0,m2,p2,2]/p2+(-1/6-m2^2 Dp2IntxnLog[0,m2,p2,0]+2IntxnLog[0,m2,p2,1]-2(-m2^2-p2)Dp2IntxnLog[0,m2,p2,1])/3/p2/.Re->Identity]/;
FreeQ[m2,Mass]&&FreeQ[p2,FourMomentum];

Dp2IntxnLog[m_/FR$MU,m_/FR$MU,0,n_Integer]:=-1/m^2*FR$MU^2/(n+2)/(n+3)/;FreeQ[m,Mass];

Dp2IntxnLog[m_/FR$MU,m_/FR$MU,m_^2/FR$MU^2,0]:=(9 - 2*Sqrt[3]*Pi)/(9*Re[m^2])*FR$MU^2/;FreeQ[m,Mass];
Dp2IntxnLog[m_/FR$MU,m_/FR$MU,m_^2/FR$MU^2,1]:=(9 - 2*Sqrt[3]*Pi)/(18*Re[m^2])*FR$MU^2/;FreeQ[m,Mass];
Dp2IntxnLog[m_/FR$MU,m_/FR$MU,m_^2/FR$MU^2,2]:=(-6 + Sqrt[3]*Pi)/(9*Re[m^2])*FR$MU^2/;FreeQ[m,Mass];

Dp2IntxnLog[m1_,m2_,0,0]:=Re[(-m1^4 + m2^4 + 4*m1^2*m2^2*Log[m1/m2])/(2*(m1^2 - m2^2)^3)]/;Not[m1===m2]&&FreeQ[m1,Mass]&&FreeQ[m2,Mass];
Dp2IntxnLog[m1_,m2_,0,1]:=Re[-(m1^6 - 6*m1^4*m2^2 + 3*m1^2*m2^4 + 2*m2^6 - 12*m1^2*m2^4*Log[m2/m1])/(6*(m1^2 - m2^2)^4)]/;Not[m1===m2]&&FreeQ[m1,Mass]&&FreeQ[m2,Mass];
Dp2IntxnLog[m1_,m2_,0,2]:=Re[(-m1^8 + 6*m1^6*m2^2 - 18*m1^4*m2^4 + 10*m1^2*m2^6 + 3*m2^8 + 24*m1^2*m2^6*Log[m1/m2])/(12*(m1^2 - m2^2)^5)]/;Not[m1===m2]&&FreeQ[m1,Mass]&&FreeQ[m2,Mass];

Dp2IntxnLog[m_,m_,p2_,0]:=Re[(Sqrt[p2*(-4*m^2 + p2)] - 2*m^2*(If[m\[Element]Reals&&p2\[Element]Reals&&m>0&&p2>4m^2,0,1]*Log[(2*m^2 - p2 + Sqrt[p2*(-4*m^2 + p2)])/(2*m^2)]+ 
  If[m\[Element]Reals&&p2\[Element]Reals&&m>0&&p2>4m^2,1,0](Log[-(2*m^2 - p2 + Sqrt[p2*(-4*m^2 + p2)])/(2*m^2)]+I Pi)))/(p2*Sqrt[p2*(-4*m^2 + p2)])]/;FreeQ[p2,FourMomentum];

Dp2IntxnLog[m1_/FR$MU,m2_/FR$MU,m1_^2/FR$MU^2,0]:= FR$MU^2*Re[((m1^2*Sqrt[-4*m1^2*m2^2 + m2^4] + (m1^2 - m2^2)*Sqrt[-4*m1^2*m2^2 + m2^4]*Log[m2/m1] + 
(-3*m1^2*m2^2 + m2^4)*Log[(m2^2 + Sqrt[-4*m1^2*m2^2 + m2^4])/(2*m1*m2)])/(m1^4*Sqrt[-4*m1^2*m2^2 + m2^4]))]/;FreeQ[m1,Mass]&&FreeQ[m2,Mass];

Dp2IntxnLog[m1_/FR$MU,m2_/FR$MU,m2_^2/FR$MU^2,0]:=Dp2IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,0];
Dp2IntxnLog[m1_/FR$MU,m2_/FR$MU,m2_^2/FR$MU^2,1]:=Dp2IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,0]-Dp2IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,1];
Dp2IntxnLog[m1_/FR$MU,m2_/FR$MU,m2_^2/FR$MU^2,2]:=Dp2IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,0]-2 Dp2IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,1]+Dp2IntxnLog[m2/FR$MU,m1/FR$MU,m2^2/FR$MU^2,2];

Dp2IntxnLog[m1_,m2_,p2_,0]:=Re[((p2*Sqrt[m1^4 + (-m2^2 + p2)^2 - 2*m1^2*(m2^2 + p2)] + (m1^2 - m2^2)*Sqrt[m1^4 + (-m2^2 + p2)^2 - 2*m1^2*(m2^2 + p2)]*Log[m2/m1] + 
  (m1^4 + m2^2*(m2^2 - p2) - m1^2*(2*m2^2 + p2))*
   (If[p2>=(m1+m2)^2&&p2\[Element]Reals&&m1\[Element]Reals&&m2\[Element]Reals,0,1]Log[(m1^2 + m2^2 + Sqrt[-4*m1^2*m2^2 + (m1^2 + m2^2 - p2)^2] - p2)/(2*m1*m2)]+
    If[p2>=(m1+m2)^2&&p2\[Element]Reals&&m1\[Element]Reals&&m2\[Element]Reals,1,0](Log[(m1^2 + m2^2 + Sqrt[-4*m1^2*m2^2 + (m1^2 + m2^2 - p2)^2] - p2)/(2*m1*m2)]+I Pi)))/
 (Sqrt[-4*m1^2*m2^2 + (m1^2 + m2^2 - p2)^2]*p2^2))]/;Not[m1^2===p2]&&Not[m2^2===p2]&&FreeQ[m1,Mass]&&FreeQ[m2,Mass]&&FreeQ[p2,FourMomentum];

Dp2IntxnLog[m1_,m2_,p2_,1]:=Re[(- IntxnLog[m1,m2,p2,1]/p2+(IntxnLog[m1,m2,p2,0]-(m1^2-m2^2-p2)Dp2IntxnLog[m1,m2,p2,0])/2/p2)/.Re->Identity]/;FreeQ[m1,Mass]&&FreeQ[m2,Mass]&&
FreeQ[p2,FourMomentum];

Dp2IntxnLog[m1_,m2_,p2_,2]:=Re[-IntxnLog[m1,m2,p2,2]/p2+(-1/6-m2^2 Dp2IntxnLog[m1,m2,p2,0]+2IntxnLog[m1,m2,p2,1]-2(m1^2-m2^2-p2)Dp2IntxnLog[m1,m2,p2,1])/3/p2/.Re->Identity]/;
FreeQ[m1,Mass]&&FreeQ[m2,Mass]&&FreeQ[p2,FourMomentum];


(* ::Subsubtitle:: *)
(*SolveDelta*)


SolveDelta[vertCT_,vertUV_,verttype_]:=Module[{fullUV,fl,fr,fsr,fsl,kuv,eqlist,deq,varlist,realvar,assumlist(*should be an option?*),res,tRule,Cond,extra,ReInt,
  Ordered2VertQ},


assumlist={};(*{2MTA<MH, 2MB<MH, 2MT>MH,MTA<MW, MT>(MW+MB),2MTA<MZ, 2MB<MZ, 2MT>MZ,MH<2MW,MH<2MZ,MW<MZ,MZ<2MW,MTA>0,MB>0,MT>0,MH>0,MZ>0,MW>0};*)
(*remove p1 using momentum conservation*)
fullUV=MergeVertList[vertCT,vertUV][[All,{1,3}]];
fullUV=fullUV/.{TensDot[a___,SlashedP[FourMomentum[1]], b___][Index[Spin, Ext[1]], Index[Spin, Ext[2]]]->
              -TensDot[a,SlashedP[FourMomentum[2]], b][Index[Spin, Ext[1]], Index[Spin, Ext[2]]]};
fullUV=fullUV/.{FourMomentum[1]->-FourMomentum[2]};

Ordered2VertQ[v1_,v2_]:=Block[{vert1,vert2},
  vert1=Sort[(v1[[1]]/.-x_->x)[[All,1,1]]];
  vert2=Sort[(v2[[1]]/.-x_->x)[[All,1,1]]];
  If[vert1[[1]]<vert2[[1]],True,If[vert1[[1]]==vert2[[1]],vert1[[2]]<vert2[[2]],False]]];

fullUV=Sort[fullUV,Ordered2VertQ];
res={};

For[kuv=1,kuv<=Length[fullUV]-0,kuv++,Print[ToString[kuv]<>" is "<>ToString[fullUV[[kuv,1]]]];
   If[kuv==1||Not[Sort[(fullUV[[kuv,1]]/.-x_->x)[[All,1,1]]]===Sort[(fullUV[[kuv-1,1]]/.-x_->x)[[All,1,1]]]],eqlist={};];
   Switch[verttype,
     {F,F},
(*Print[InputForm[Union[Cases[fullUV[[kuv,2]],_IntxnLog,\[Infinity]]]]];*)
     (*fermions*)
     fl=Coefficient[fullUV[[kuv,2]],TensDot[SlashedP[FourMomentum[2]], ProjM][Index[Spin, Ext[1]], Index[Spin, Ext[2]]]];
     fr=Coefficient[fullUV[[kuv,2]],TensDot[SlashedP[FourMomentum[2]], ProjP][Index[Spin, Ext[1]], Index[Spin, Ext[2]]]];
     fsl=Coefficient[fullUV[[kuv,2]],ProjM[Index[Spin, Ext[1]], Index[Spin, Ext[2]]]];
     fsr=Coefficient[fullUV[[kuv,2]],ProjP[Index[Spin, Ext[1]], Index[Spin, Ext[2]]]];

     (*check a piori on the form of the vertex*)
     If[Not[Simplify[fullUV[[kuv,2]]-fl*TensDot[SlashedP[FourMomentum[2]], ProjM][Index[Spin, Ext[1]], Index[Spin, Ext[2]]]-
                    fr*TensDot[SlashedP[FourMomentum[2]], ProjP][Index[Spin, Ext[1]], Index[Spin, Ext[2]]]-
                    fsl ProjM[Index[Spin, Ext[1]], Index[Spin, Ext[2]]]- fsr ProjP[Index[Spin, Ext[1]], Index[Spin, Ext[2]]]]===0],
       Print["error: FF vertex does not match expected form"]];

     eqlist=Append[eqlist,If[TheMass[fullUV[[kuv,1,1,1]]]===0&&TheMass[fullUV[[kuv,1,2,1]]]===0,
         Normal[Series[Refine[(Simplify[fl==0,Assumptions->{IndexDelta[Index[Colour, Ext[1]], Index[Colour, Ext[2]]]!=0}]/.
           SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,2,1]]]^2),Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]],
         Normal[Series[Refine[(Simplify[fl*TheMass[fullUV[[kuv,1,1,1]]]+fsr==0,Assumptions->{IndexDelta[Index[Colour, Ext[1]], Index[Colour, Ext[2]]]!=0}]/.
           SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,2,1]]]^2),Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]]]];
     eqlist=Append[eqlist,
         Normal[Series[Refine[(Simplify[fr*TheMass[fullUV[[kuv,1,1,1]]]+fsl==0,Assumptions->{IndexDelta[Index[Colour, Ext[1]], Index[Colour, Ext[2]]]!=0}]/.
           SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,1,1]]]^2),Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]]];
     (*If[Not[FreeQ[fullUV[[kuv,1,1]],fullUV[[kuv,1,2,1]]]], be careful for non diagonal piece*)
      deq=Simplify[D[(TheMass[fullUV[[kuv,1,1,1]]](fr+fl)+fsl+fsr),SP[FourMomentum[2], FourMomentum[2]]]];
      deq=(deq/.Derivative[0, 0, 1, 0][IntxnLog][aa_, bb_, cc_, dd_]:>Dp2IntxnLog[aa,bb,cc,dd]);
      deq=Simplify[2*TheMass[fullUV[[kuv,1,1,1]]]*deq+(fr+fl)==0,Assumptions->{IndexDelta[Index[Colour, Ext[1]], Index[Colour, Ext[2]]]!=0}]/.SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,1,1]]]^2;
      deq=Normal[Series[Refine[deq,Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]];

      If[DeleteCases[fullUV[[kuv,1,1,1]],Index[Gluon|Colour,_],\[Infinity]]===DeleteCases[fullUV[[kuv,1,2,1]],Index[Gluon|Colour,_],\[Infinity]]||
        DeleteCases[fullUV[[kuv,1,1,1]],Index[Gluon|Colour,_],\[Infinity]]===-DeleteCases[fullUV[[kuv,1,2,1]],Index[Gluon|Colour,_],\[Infinity]](*Same particle or particle and antiparticle*),
        eqlist=Append[eqlist,deq];
      ];
     ,
	 {S,S},
     If[DeleteCases[fullUV[[kuv,1,1,1]],Index[Gluon|Colour,_],\[Infinity]]===DeleteCases[fullUV[[kuv,1,2,1]],Index[Gluon|Colour,_],\[Infinity]]||
          DeleteCases[fullUV[[kuv,1,1,1]],Index[Gluon|Colour,_],\[Infinity]]===-DeleteCases[fullUV[[kuv,1,2,1]],Index[Gluon|Colour,_],\[Infinity]](*Same particle or particle and antiparticle*),
       If[FreeQ[FR$GoldstoneList,fullUV[[kuv,1,1,1]]]&&FreeQ[FR$GoldstoneList,-fullUV[[kuv,1,1,1]]](*Not a goldstone boson*),
         eqlist=Append[eqlist,Normal[Series[Refine[(fullUV[[kuv,2]]/.SP[FourMomentum[2],FourMomentum[2]]->
                TheMass[fullUV[[kuv,1,1,1]]]^2),Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]]==0];
       ];
       deq=D[(fullUV[[kuv,2]]),SP[FourMomentum[2], FourMomentum[2]]];
       deq=(deq/.Derivative[0, 0, 1, 0][IntxnLog][aa_, bb_, cc_, dd_]:>Dp2IntxnLog[aa,bb,cc,dd])/.SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,1,1]]]^2;
       deq=Normal[Series[Refine[deq,Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]];
       eqlist=Append[eqlist,deq==0];, 
       (*different scalars*)
       eqlist=Append[eqlist,Normal[Series[Refine[(fullUV[[kuv,2]]/.SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,1,1]]]^2),Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]]==0];
       eqlist=Append[eqlist,Normal[Series[Refine[(fullUV[[kuv,2]]/.SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,2,1]]]^2),Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]]==0];
     ];,
     {V,V},
     If[FreeQ[DeleteCases[fullUV[[kuv,1,1,1]],Index[Gluon|Colour,xx_],\[Infinity]],DeleteCases[fullUV[[kuv,1,2,1]],Index[Gluon|Colour,xx_],\[Infinity]]],
     (*different particles*)
     deq = Coefficient[fullUV[[kuv,2]],ME[Index[Lorentz,Ext[1]],Index[Lorentz,Ext[2]]]];
     eqlist=Append[eqlist,(Normal[Series[Refine[deq/.SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,1,1]]]^2,
              Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]])==0];
     eqlist=Append[eqlist,(Normal[Series[Refine[deq/.SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,2,1]]]^2,
              Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]])==0];,
     (*Same particle*)
     deq = Coefficient[fullUV[[kuv,2]],ME[Index[Lorentz,Ext[1]],Index[Lorentz,Ext[2]]]];
     If[Not[TheMass[fullUV[[kuv,1,1,1]]]===0],
       (*Massive*)
       eqlist=Append[eqlist,Normal[Series[Refine[(deq/.SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,1,1]]]^2),
         Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]]==0]];
     deq=D[deq,SP[FourMomentum[2], FourMomentum[2]]];
     deq=(deq/.Derivative[0, 0, 1, 0][IntxnLog][aa_, bb_, cc_, dd_]:>Dp2IntxnLog[aa,bb,cc,dd]);
     deq=deq/.SP[FourMomentum[2],FourMomentum[2]]->TheMass[fullUV[[kuv,1,1,1]]]^2;
     deq=Normal[Series[Refine[deq,Assumptions->Join[{FR$MU>0,FR$IR>0,FR$Eps>0},assumlist]],{FR$Eps,0,0}]];
     eqlist=Append[eqlist,deq==0];
    ];
  ];(*end switch*)
  varlist=Union[If[TheMass[fullUV[[kuv,1,1,1]]]===TheMass[fullUV[[kuv,1,2,1]]],Cases[eqlist,FR$delta[{TheMass[fullUV[[kuv,1,1,1]]]},{}],\[Infinity]],{}],Cases[eqlist,_FR$deltaZ,\[Infinity]]];
  realvar=Cases[varlist,FR$deltaZ[{x_,x_},b__]];
  extra=Complement[Union[Cases[fullUV[[kuv,2]],_FR$deltaZ,\[Infinity]]],varlist];
  eqlist=(Simplify[#/.((Rule[Conjugate[#],#]&)/@realvar)(*/.FR$IR->1*),Assumptions->Append[assumlist,FR$MU>0],TimeConstraint->0.02]&)/@eqlist;
Print[InputForm[ToSolve[eqlist/.{If->Cond,Re->ReInt},varlist]]];
(*Print["before join"];Print[InputForm[If[Length[extra]>0,(#->0&)/@extra,{}]]];*)
  If[kuv==Length[fullUV]||Not[Sort[(fullUV[[kuv,1]]/.-x_->x)[[All,1,1]]]===Sort[(fullUV[[kuv+1,1]]/.-x_->x)[[All,1,1]]]],
    res=Join[res,Solve[eqlist/.{If->Cond,Re->ReInt},varlist][[1]],If[Length[extra]>0,(#->0&)/@extra,{}]];  (*set to zero no usefull delta(Z)*)
  ];Print["solve done"];
];(*end for*)
Print["end of the for"];

res=Flatten[res/.Rule->tRule/.{tRule[FR$deltaZ[{a_,b_},xx__],-Conjugate[FR$deltaZ[{b_,a_},xx__]]]->{tRule[FR$deltaZ[{a,b},xx],0],tRule[FR$deltaZ[{b,a},xx],0]}}]/.
  tRule->Rule;
Print["before last simplification"];
(Simplify[#,Assumptions->Append[assumlist,FR$MU>0],TimeConstraint->0.02]&)/@(res/.{ReInt->Re,Cond->If})

];


(* ::Subtitle:: *)
(*R2*)


(*Performs q integration*)


(* ::Subsubtitle::Closed:: *)
(*Tadpoles*)


R2Tadpoles::usage="compute the R2 for a Tadpole amplitude after the non contributing terms have been removed"


R2Tadpoles[num_,del_,next_] := Block[{tmp, coef},
   coef = Times @@ Cases[num, _?((FreeQ[#, MatrixTrace] && FreeQ[#, MetricTensor] && FreeQ[#, FourVector]) &)];
   tmp = num/coef;
   tmp = Expand[tmp /. {MatrixTrace -> DTr, MetricTensor -> ME, FourVector -> FV}];
   Simplify[-2 I \[Pi]^2*coef*del*Coefficient[Normal[Series[tmp /. {M$dim -> 4 + M$eps}, {M$eps, 0, 1}]], M$eps]]*FR$R2-2 I \[Pi]^2*coef*del*Simplify[(FR$UV(1/FR$Eps*
   Normal[Series[(tmp/.M$dim->4+FR$Eps),{FR$Eps,0,If[next===2,1,0]}]]+ If[next===2,1,0]*UVLog[del/FR$MU^2]/2*tmp/.{M$dim->4}))] /. {DTr -> MatrixTrace, ME -> MetricTensor, 
   FV -> FourVector}
];


(* ::Subsubtitle::Closed:: *)
(*Bubbles*)


R2BubblesF::usage="compute the R2 for a bubble amplitude with no occurence of the loop momentum on the numerator after the non contributing terms have 
been removed"
R2BubblesQ2::usage="compute the R2 for a bubble amplitude with two occurence of the loop momentum on the numerator after the non contributing terms have 
been removed"


R2BubblesF[num_,del_,next_]:=Block[{tmp,coef},
  coef=Times@@Cases[num,_?((FreeQ[#,MatrixTrace]&&FreeQ[#,MetricTensor]&&FreeQ[#,FourVector]&&FreeQ[#,FermionChain])&)];
  tmp=num/coef;
  tmp=Expand[tmp/.{MatrixTrace->DTr,MetricTensor->ME,FourVector->FV,FermionChain->FCh}];
  tmp=-2I \[Pi]^2*(Simplify[Coefficient[Normal[Series[tmp/.{M$dim->4+M$eps},{M$eps,0,1}]],M$eps]]FR$R2 + 
               FR$UV*Simplify[1/FR$Eps*Normal[Series[(tmp/.M$dim->4+FR$Eps),{FR$Eps,0,If[next===2,1,0]}]]+ If[next===2,1,0]*(Normal[Series[Log[del/FR$MU^2]/2*tmp/.{M$dim->4+FR$Eps},{FR$Eps,0,1}]])]/.
               {DTr->MatrixTrace,ME->MetricTensor,FV->FourVector,FCh->FermionChain});
  coef*tmp
];

R2BubblesQ2[num_,lm_,del_,next_]:=Block[{tmp,coef},
  coef=Times@@Cases[num,_?((FreeQ[#,MatrixTrace]&&FreeQ[#,MetricTensor]&&FreeQ[#,FourVector]&&FreeQ[#,FermionChain])&)];
  tmp=num/coef;
  tmp=Expand[tmp/.{MatrixTrace->DTr,MetricTensor->ME,FourVector->FV,FermionChain->FCh}];
  tmp=tmp/.red2v[lm];
  tmp=Simplify[-4I \[Pi]^2*Coefficient[Normal[Series[tmp/.{M$dim->4+M$eps},{M$eps,0,1}]],M$eps]-I \[Pi]^2tmp/.{M$dim->4}]*del*FR$R2 - 
                I \[Pi]^2*del FR$UV*Simplify[(1/FR$Eps*Normal[Series[((4-FR$Eps)tmp/.M$dim->4+FR$Eps),{FR$Eps,0,If[next===2,1,0]}]]+ If[next===2,1,0]*Normal[Series[(4-FR$Eps)Log[del/FR$MU^2]/2*tmp/.{M$dim->4+FR$Eps},{FR$Eps,0,1}]])]/.
               {SP[lm,lm]->1}/.{DTr->MatrixTrace,ME->MetricTensor,FV->FourVector,FCh->FermionChain};
  coef*tmp
];


(* ::Subsubtitle::Closed:: *)
(*Triangles*)


R2FTriangles::usage="compute the R2 for a Triangle amplitude with only fermions in the loop after the non contributing terms have been removed"
R2FBTriangles::usage="compute the R2 for a Triangle amplitude with both fermions and bosons in the loop after the non contributing terms have been removed"
R2BTriangles::usage="compute the R2 for a Triangle amplitude with only bosons in the loop after the non contributing terms have been removed"


R2FTriangles[num_,lm_]:=Block[{tmp,coef},
  tmp=Cases[num,_MatrixTrace][[1]];
  coef=Times@@Cases[num,Except[_MatrixTrace]];
  tmp=tmp/.{MatrixTrace->DTr};
  tmp=Expand[Expand[tmp]/.red2v[lm]];
  I \[Pi]^2*coef (Simplify[-2Coefficient[Normal[Series[tmp/.{M$dim->4+M$eps},{M$eps,0,1}]],M$eps]-1/2(tmp/.M$dim->4)]FR$R2- 
               2 FR$UV*Simplify[(tmp/.{M$dim->4})]/FR$Eps)/.SP[lm,lm]->1/.{ME->MetricTensor,DTr->MatrixTrace,FV->FourVector}
];

R2FBTriangles[num_,lm_]:=Block[{tmp,coef},
  tmp=Times@@Cases[num,_?(Not[FreeQ[#,lm]]&)];
  coef=num/tmp;
  tmp=Times@@Cases[coef,_?(Not[FreeQ[#,MetricTensor]]&)]*tmp;
  coef=num/tmp;
  tmp=Expand[tmp/.{FermionChain->FCh,MetricTensor->ME,FourVector->FV}];
  tmp=Expand[Expand[tmp]/.red2v[lm]];
  I \[Pi]^2*coef (Simplify[-2Coefficient[Normal[Series[tmp/.{M$dim->4+M$eps},{M$eps,0,1}]],M$eps]-1/2(tmp/.M$dim->4)]FR$R2- 
               2 FR$UV*Simplify[(tmp/.{M$dim->4})]/FR$Eps)/.SP[lm,lm]->1/.{ME->MetricTensor,DTr->MatrixTrace,FV->FourVector,FCh->FermionChain}
];

R2BTriangles[num_,lm_]:=Block[{coef, tmp},
  tmp=Times@@Cases[num,_?(Not[FreeQ[#,lm]]&)];
  coef=num/tmp;
  tmp=Times@@Cases[coef,_?(Not[FreeQ[#,MetricTensor]]&)]*tmp;
  coef=num/tmp;
  tmp=Expand[tmp]/.{MetricTensor->ME}/.{FourVector->FV};
  tmp=Expand[tmp]/.red2v[lm];
  tmp=Expand[tmp];
  I \[Pi]^2*coef (Simplify[-2Coefficient[Normal[Series[tmp/.{M$dim->4+M$eps},{M$eps,0,1}]],M$eps]-1/2(tmp/.M$dim->4)]FR$R2- 
               2 FR$UV*Simplify[(tmp/.{M$dim->4})]/FR$Eps)/.SP[lm,lm]->1/.{ME->MetricTensor,DTr->MatrixTrace,FV->FourVector}
];


(* ::Subsubtitle::Closed:: *)
(*Boxes*)


R2FBoxes::usage="compute the R2 for a boxe amplitude with only fermions in the loop after the non contributing terms have been removed"
R2BBoxes::usage="compute the R2 for a boxe amplitude with only bosons in the loop after the non contributing terms have been removed"


R2FBoxes[num_,lm_]:=Block[{tmp,coef},
  tmp=Cases[num,_MatrixTrace][[1]];
  coef=Times@@Cases[num,Except[_MatrixTrace]];
  tmp=tmp/.{MatrixTrace->DTr};
  tmp=Simplify[Expand[Expand[tmp]/.red4v[lm]/.red2v[lm]]];
  I \[Pi]^2*coef (Simplify[-2Coefficient[Normal[Series[tmp/.{M$dim->4+M$eps},{M$eps,0,1}]],M$eps]-5/6(tmp/.M$dim->4)]FR$R2- 
               2 FR$UV*Simplify[(tmp/.{M$dim->4})]/FR$Eps)/.SP[lm,lm]->1/.{ME->MetricTensor,FV->FourVector}
];

R2BBoxes[num_,lm_]:=Block[{coef, tmp},
  tmp=Times@@Cases[num,_?(Not[FreeQ[#,lm]]&)];
  coef=num/tmp;
  tmp=Times@@Cases[coef,_?(Not[FreeQ[#,MetricTensor]]&)]*tmp;
  coef=num/tmp;
  tmp=Expand[tmp]/.{MetricTensor->ME}/.{FourVector->FV}/.red4v[lm]/.red2v[lm];
  tmp=Expand[tmp];
  I \[Pi]^2*coef (Simplify[-2Coefficient[Normal[Series[tmp/.{M$dim->4+M$eps},{M$eps,0,1}]],M$eps]-5/6(tmp/.M$dim->4)]FR$R2- 
               2 FR$UV*Simplify[(tmp/.{M$dim->4})]/FR$Eps)/.SP[lm,lm]->1/.{ME->MetricTensor,FV->FourVector}
];


(* ::Subtitle:: *)
(*GetR2*)


GetR2::usage="compute the R2 for a FeynArts amplitude at the generic level"


GetR2[amp_,next_]:=Block[{temp,loopmom,nden,den,num,kk,extmom,extmom2,intmom,lmcoef,DisFC,DisMT,DisDS,DisNC,nlmom,x,y,z,rl,lr,delta,noextpos,feynp,PW,pairtmp,
                    perm,mep,mep2,fvp,nq,LoopInt,Epsi,sumind,feynpar,RemoveHDim,FChain,evenOnly=True,maxDim=6(*should be option*),tempUV,FullLoopInt},

(*pole of the integral of (k^2)^nm/(k^2-dlta)^nm*Gamma[2-d/2]*)
  LoopInt[nm_,nd_,dlta_]:=(-1)^(nm+nd)*(\[Pi])^2 I Product[M$dim+2kk,{kk,0,nm-1}]/2^(nm)dlta^(2-nd+nm)/Gamma[nd]*
                          If[nd-nm-2>0,0,1/Product[nd-2-nm+ll,{ll,0,1-nd+nm}]];
(*the integral of (k^2)^nm/(k^2-dlta)^nm*Gamma[2-d/2]*)
  FullLoopInt[nm_,nd_,dlta_]:=(-1)^(nm+nd)*(\[Pi])^2 I Product[M$dim+2kk,{kk,0,nm-1}]/2^(nm)dlta^(2-nd+nm)/Gamma[nd]*
                          If[nd-nm-2>0,Gamma[nd-2-nm],1/(2-M$dim/2)(1-(2-M$dim/2)Log[dlta/FR$MU^2])/Product[nd-M$dim/2-nm+ll,{ll,0,1-nd+nm}]];

  RemoveHDim[am_]:= If[If[evenOnly,
    Total[Ceiling[(Cases[am,G[_][_][x__][y__]:>Count[{x},F[__],\[Infinity]]*3/2+Count[{x},V[__],\[Infinity]]+Count[{x},S[__],\[Infinity]]+Count[{y},Mom[_],\[Infinity]]]-4)/2]*2],
    Total[(Cases[am,G[_][_][x__][y__]:>Count[{x},F[__],\[Infinity]]*3/2+Count[{x},V[__],\[Infinity]]+Count[{x},S[__],\[Infinity]]+Count[{y},Mom[_],\[Infinity]]]-4)]]>maxDim-4,0,am];

  loopmom=amp[[2,1]];
  den=Cases[amp[[3]],_PropagatorDenominator,\[Infinity]];
  nden=Length[den];
  num=amp[[3]]/.{FeynAmpDenominator[__]->1};
(*compute momentum shift and delta*)
  Switch[nden,
    1,
    delta = den[[1,2]]^2;
    ,
    _,
    feynpar=Table[ToExpression["x"<>ToString[fcount]],{fcount,1,nden-1}]; (*list of feynman parameters*)
    extmom2=(Coefficient[#,loopmom]&/@Cases[den,\!\(TraditionalForm\`PropagatorDenominator\)[a_,b_]->a])*Cases[den,\!\(TraditionalForm\`PropagatorDenominator\)[a_,b_]->a]-loopmom;(*vector of external momenta*)
    noextpos=Position[extmom2,0][[1,1]];
    feynp=Insert[feynpar,(1-Total[feynpar]),noextpos];
    extmom=extmom2.feynp;(*momentum shift*)
    delta = -Simplify[Cases[den,\!\(TraditionalForm\`PropagatorDenominator\)[a_,b_]->SP[a,a]-b^2].feynp-SP[loopmom+extmom,loopmom+extmom]];
    num=num/.loopmom ->loopmom-extmom;
  ];(*end of Switch*)
  Attributes[pairtmp]=Orderless;
  num = Expand[num/.{MetricTensor->ME,FourVector->FV,MatrixTrace->DTr,FermionChain->FChDist,ScalarProduct->SP,MatrixTrace->DTr}];
(*remove terms with sum(dimension of the vertices-4)>maxDim-4*)
  If[Head[num]===Plus,num=RemoveHDim/@num;,num=RemoveHDim[num];];

  num = Expand[num/.FChDist->FCh];Print["after total expand"];Print[Length[num]];

  num = (num/.FCh[aa__,NonCommutative[DiracSlash[loopmom]],bb__]:>
             FCh[aa,NonCommutative[DiracMatrix[Index[Lorentz,sumind]]],bb]*pairtmp[loopmom,Index[Lorentz, sumind]]);
  num=(num/.{FV[loopmom,aa_]->pairtmp[loopmom,aa],
             SP[loopmom,aa_?(FreeQ[#,Internal]&)]->pairtmp[loopmom,aa]});
(*expand power for the counting of loop momenta*)
  num = num/.{Power[pairtmp[loopmom,aa_?(FreeQ[#,loopmom]&)],n_Integer]:>
             Product[pairtmp[loopmom,PW[kk]*aa],{kk,1,n}],Power[SP[loopmom,loopmom],n_Integer]:>
             Product[SP[loopmom,PW[loopmom,kk]],{kk,1,n}]};
(*remove odd power of loop momenta and finite terms except for two point functions*)
  num=DeleteCases[num,_?(OddQ[Count[#,loopmom,\[Infinity]]]&)];
  If[next>2,num=DeleteCases[num,_?(Count[#,loopmom,\[Infinity]]<2nden-4&)];];
  num=num/.PW[loopmom,_]->loopmom;

(*If[nden>3,Print[Style["removing odd",Blue]];Print[num];];*)

(*replace k^mu1 ... k^mu2n by (k^2)^n*metric combination*)
  For[nq=Max[(Count[#,loopmom,\[Infinity]]&)/@List@@num],nq>0,nq=nq-2,
    perm=Permutations[Table[ToExpression["ind"<>ToString[i]],{i,nq}]];
    mep=Simplify[Plus@@Table[Times@@Table[ME[Index[Lorentz,ToExpression["ind"<>ToString[2*i-1]]],Index[Lorentz,ToExpression["ind"<>ToString[2*i]]]],{i,nq/2}]/.Table[perm[[1,kk]]->perm[[nn,kk]],{kk,nq}],{nn,Length[perm]}]/Factorial[nq/2]/2^(nq/2)];
    mep2=Simplify[mep/Expand[mep*If[Head[mep]===Plus,mep[[1]],mep]]]*SP[loopmom,loopmom]^(nq/2);
    fvp=Times@@(pairtmp[loopmom,Pattern[#,_]]&/@perm[[1]]);
    num = num/.LeviCivita[a___,loopmom,b___]->Module[{lcind},LeviCivita[a,Index[Lorentz,lcind],b]pairtmp[loopmom,Index[Lorentz,lcind]]];
    num=num/.{fvp->(mep2/.{Index[Lorentz,aa_]->aa})};
  ];

  num=num/.{PW[_]->1};
  num=num/.{ME[a_?(Not[FreeQ[#,FourMomentum]]&),Index[Lorentz,b_]]->FV[a,Index[Lorentz,b]],
            ME[a_?(Not[FreeQ[#,FourMomentum]]&),a_?(Not[FreeQ[#,FourMomentum]]&)]->SP[a,b]};

  temp=0;
  tempUV=0;
  For[kk=Max[Join[Cases[num,SP[loopmom, loopmom]^nexp_Integer->nexp,\[Infinity]],Cases[num,SP[loopmom, loopmom]->1,\[Infinity]],{0}]],kk>=0,kk--,
  (*-2 factor is from -2/Epsilon*)
    temp=temp-2 *Coefficient[Normal[Series[Coefficient[num,SP[loopmom,loopmom],kk]*LoopInt[kk,nden,delta]/.M$dim->4+Epsi,{Epsi,0,1}]],Epsi];
    tempUV=tempUV+If[next>2,-2/FR$Eps Normal[Series[Coefficient[num,SP[loopmom,loopmom],kk]*LoopInt[kk,nden,delta]/.M$dim->4+Epsi,{Epsi,0,0}]],
                    Normal[Series[Coefficient[num,SP[loopmom,loopmom],kk]*FullLoopInt[kk,nden,delta]/.M$dim->4+FR$Eps,{FR$Eps,0,1}]]/.Log[Mass[x__]^2/FR$MU^2]->UVLog[Mass[x]^2/FR$MU^2]];
  ];

  Switch[nden,
    1,
    num=temp*FR$R2+tempUV*FR$UV;,
    _,
    num = If[next>2,
            Factorial[nden-1]*Integrate[FR$R2*temp+tempUV*FR$UV,Sequence@@Table[{feynpar[[kk]],0,1-Total[feynpar[[1;;kk-1]]]},{kk,nden-1}]],
		    XIntegrate[FR$R2*temp+tempUV*FR$UV,delta,x1]
          ];
  ];(*end second switch*)
  num/.{FCh->FermionChain,FV->FourVector,ME->MetricTensor,SP->ScalarProduct}
];(*end of getR2*)


(* ::Subtitle:: *)
(*Get Several R2*)


(* from generic to class in FR format*)


(* ::Subsubtitle::Closed:: *)
(*R2atClass*)


R2atClass::usage="compute the R2 for a FeynArts amplitude #1 and a topology #2 at the class level. All the contributions are summed but the wave functions keep track of the external particles.
#3 is the list of the generic external particles (F,S,V). If #4 is True, each contribution is multiplyed by IPL[list of the particles in the loop]. The last argument is the list of indices to
be kept in IPL like generation."


R2atClass[Ampl_,INTopo_,verttype_,lab_,kept_,qcd_]:=Block[{res,tmp,tmp2,top,gen,rf,rl,fc,scalar,intern,nExtC,ll}, 
  res=0;
  For[kk=1,kk<=Length[Ampl],kk++,
    top=Ampl[[kk,1,1,2]];
    gen=Ampl[[kk,1,2,2]];
    tmp=GetR2[Ampl[[kk]],Length[verttype]](*ToExpression["amp"<>ToString[kk]]*);(*Print[Ampl[[kk]]];Print[InputForm[Simplify[tmp/.FR$R2->0]]];*)
    If[Not[tmp===0],
      For[ll=1,ll<=Length[Ampl[[kk,4,2]]],ll++,
        scalar=1;
        fc=INTopo[[top,2,gen,2,ll]];
        rf={};
        If[OrderedQ[fc[[1;;Length[verttype],2]],(PartOrder[#1,#2]&)],
          For[nn=1,nn<=Length[Ampl[[0,1,2,1]]],nn++,
            rf=Append[rf,Rule[Ampl[[0,1,2,1,nn,1]][Index[Generic,fc[[nn,1,1]]]],fc[[nn,2]]]];
            If[Ampl[[0,1,2,1,nn,1]]===S,
              scalar=scalar*SWF[fc[[nn,2]],fc[[nn,1,1]]];
            ];
            If[Ampl[[0,1,2,1,nn,1]]===F,
              rf=Append[rf,Rule[SpinorType[Index[Generic,fc[[nn,1,1]]]],SpinorType[fc[[nn,2]]]]];
            ];
          ];
          rl=Table[Rule[Ampl[[kk,4,1,nn]],Ampl[[kk,4,2,ll,nn]]],{nn,1,Length[Ampl[[kk,4,1]]]}];
          rl=rl/.{FourVector->FV}/.{FV->FourVector};
          intern=If[lab,IPL[Sort[List@@Union[DeleteCases[fc[[Length[verttype]+1;;,2]]/.{-1->1},Index[Except[Alternatives@@kept],__],\[Infinity]]],PartOrder]],1];
          tmp2=scalar*intern*tmp/.{FourVector->FV}/.{FV->FourVector}/.rl/.rf/.If[Length[M$FACouplings]>0,M$FACouplings,{}];
          If[qcd,
            nExtC=Length[Union[Cases[tmp2,Index[Colour,a_?(#<=Length[verttype]&)],\[Infinity]],Cases[tmp2,Index[Gluon,a_?(#<=Length[verttype]&)],\[Infinity]]]];
            tmp2= Total[Table[Coefficient[tmp2,GS,nn]*GS^nn,{nn,Max[{2,nExtC}],4}]];
          ];
          res=res+tmp2;
        ];
      ];
    ];
  ];
  res
]


(* ::Subsubtitle::Closed:: *)
(*CTatClass*)


CTatClass::usage="compute the CT for a FeynArts amplitude #1 and a topology #2 at the class level. All the contributions are summed but the wave functions keep track of the external particles.
#3 is the list of the generic external particles (F,S,V). "


CTatClass[Ampl_,INTopo_,verttype_]:=Block[{res,tmp,tmp2,top,gen,rf,rl,fc,scalar,ll}, 
  res=0;
  For[kk=1,kk<=Length[Ampl],kk++,
    top=Ampl[[kk,1,1,2]];
    gen=Ampl[[kk,1,2,2]];
    tmp=Simplify[Expand[Ampl[[kk,3]]/.{FermionChain->FCh}]]/.{FCh->FermionChain};
    
    If[Not[tmp===0],
      For[ll=1,ll<=Length[Ampl[[kk,4,2]]],ll++,
        scalar=1;
        fc=INTopo[[top,2,gen,2,ll]];
        rf={};
        If[OrderedQ[fc[[1;;Length[verttype],2]],(PartOrder[#1,#2]&)],
          For[nn=1,nn<=Length[Ampl[[0,1,2,1]]],nn++,
            rf=Append[rf,Rule[Ampl[[0,1,2,1,nn,1]][Index[Generic,fc[[nn,1,1]]]],fc[[nn,2]]]];
            If[Ampl[[0,1,2,1,nn,1]]===S,
              scalar=scalar*SWF[fc[[nn,2]],fc[[nn,1,1]]];
            ];
            If[Ampl[[0,1,2,1,nn,1]]===F,
              rf=Append[rf,Rule[SpinorType[Index[Generic,fc[[nn,1,1]]]],SpinorType[fc[[nn,2]]]]];
            ];
          ];
          rl=Table[Rule[Ampl[[kk,4,1,nn]],Ampl[[kk,4,2,ll,nn]]],{nn,1,Length[Ampl[[kk,4,1]]]}];
          rl=rl/.{FourVector->FV}/.{FV->FourVector};
          tmp2=scalar*tmp/.{FourVector->FV}/.{FV->FourVector}/.rl/.rf/.If[Length[M$FACouplings]>0,M$FACouplings,{}];
          
          res=res+tmp2;
        ];
      ];
    ];
  ];
  res*FR$UV
]


(* ::Subsubtitle::Closed:: *)
(*PartOrder*)


PartOrder::usage="Ordering function for the particles, fermion goes before vector which goes before scalar, antiparticles goes before particles of the same type, remaining ordering is done 
according to the particle numbers"


PartOrder[V[x_,ex1___],-V[y_,ex2___]]:=False;
PartOrder[V[x_,ex1___],V[y_,ex2___]]:=False/;y<x;
PartOrder[-V[x_,ex1___],-V[y_,ex2___]]:=False/;y<x;
PartOrder[-S[x_,ex1___]|S[x_,ex1___],V[y_,ex2___]|-V[y_,ex2___]]:=False;
PartOrder[S[x_,ex1___],-S[y_,ex2___]]:=False;
PartOrder[S[x_,ex1___],S[y_,ex2___]]:=False/;y<x;
PartOrder[-S[x_,ex1___],-S[y_,ex2___]]:=False/;y<x;
PartOrder[-F[x_,ex1___],-F[y_,ex2___]]:=False/;y<x;
PartOrder[F[x_,ex1___],F[y_,ex2___]]:=False/;y<x;
PartOrder[S[y_,ex1___]|V[y_,ex1___],F[x_,ex2___]]:=False;
PartOrder[-S[y_,ex1___]|-V[y_,ex1___],F[x_,ex2___]]:=False;
PartOrder[-Except[F,X_][y_,ex1___],-F[x_,ex2___]]:=False;
PartOrder[Except[Times,X_][y_,ex1___],-F[x_,ex2___]]:=False;
PartOrder[F[y_,ex1___],-F[x_,ex2___]]:=False;

PartOrder[V[x_,ex1___],-V[x_,ex2___]]:=False;
PartOrder[-S[x_,ex1___]|S[x_,ex1___],V[x_,ex2___]|-V[x_,ex2___]]:=False;
PartOrder[S[x_,ex1___],-S[x_,ex2___]]:=False;
PartOrder[S[x_,ex1___]|V[x_,ex1___],F[x_,ex2___]]:=False;
PartOrder[-S[x_,ex1___]|-V[x_,ex1___],F[x_,ex2___]]:=False;
PartOrder[X_[x_,ex1___],-F[x_,ex2___]]:=False;
PartOrder[F[x_,ex1___],-F[x_,ex2___]]:=False;

PartOrder[x_,x_]:=True;
PartOrder[x_,y_]:=True;


(* ::Subsubtitle::Closed:: *)
(*R2vertlist*)


R2vertlist::usage="rewrite the result from Expand[R2atClass] into a list of vertices similar to those of FeynRules"


R2vertlist[vertsum_]:=Block[{vertList,lab,ver,pos,ss,v,ind,indrep,f,indlist,ind2,xx,ll,kk,nn,ff,mm,UVmass,deltarule,indl,aa,bb,s1,s2,extsum,wft,inveps,res},
  vertList={};
  If[Head[vertsum]===Times,res={vertsum};,res=vertsum];
  For[kk=1,kk<=Length[res],kk++,
    lab=Cases[res[[kk]],PolarizationVector[v_,FourMomentum[Incoming,ind_],Index[Lorentz,ind2_]]->{v,ind,ind2}];
    lab=Join[lab,Cases[res[[kk]],SWF[ss_,ind_]->{ss,ind}]];
    lab=Join[lab,Cases[res[[kk]],FermionChain[NonCommutative[SpinorType[f_][-FourMomentum[Incoming,ind_]|FourMomentum[Incoming,ind_],__]],xx__]->{f,ind,1},\[Infinity]]];
    lab=Join[lab,Cases[res[[kk]],FermionChain[xx__,NonCommutative[SpinorType[f_][-FourMomentum[Incoming,ind_]|FourMomentum[Incoming,ind_],__]]]->{f,ind,2},\[Infinity]]];
    lab=Sort[lab,PartOrder[#1[[1]],#2[[1]]]&];
    (*Majorana reordering*)
    If[lab[[1,1,0]]===F&&lab[[2,1,0]]===F&&lab[[1,1,1]]===lab[[2,1,1]]&&lab[[1,2]]>lab[[2,2]],lab[[1;;2]]=lab[[{2,1}]];];
    indrep=Table[FourMomentum[Incoming,lab[[ll,2]]]->FourMomentum[Incoming,ll],{ll,1,Length[lab]}];
    For[nn=1,nn<=Length[lab],nn++,
      If[MatchQ[lab[[nn,1]],V[__]|-V[__]],indrep=Append[indrep,  Index[Lorentz,lab[[nn,3]]]->Index[Lorentz,Ext[nn]]  ]; ];
      If[MatchQ[lab[[nn,1]],F[__]|-F[__]],indrep=Append[indrep,  Index[Spin,Ext[lab[[nn,3]]]]->Index[Spin,Ext[nn]]  ]; ];
      indlist=Cases[lab[[nn,1]],Index[xx_,ind2_],\[Infinity]];
      For[ll=1,ll<=Length[indlist],ll++,indrep=Append[indrep,indlist[[ll]]->Index[indlist[[ll,1]],Ext[nn]]];];
    ];
    ver=DeleteCases[res[[kk]],PolarizationVector[v_,FourMomentum[Incoming,ind_],Index[Lorentz,ind2_]]];
    ver=DeleteCases[ver,SWF[ss_,ind_]];
    ver=Replace[ver,{FermionChain[NonCommutative[SpinorType[-F[f__]|F[f__]][-FourMomentum[Incoming,ind_]|FourMomentum[Incoming,ind_],__]],xx___,
                     NonCommutative[SpinorType[F[f2__]|-F[f2__]][-FourMomentum[Incoming,ind2_]|FourMomentum[Incoming,ind2_],__]]]->
                      -TensDot[xx][Index[Spin,Ext[1]],Index[Spin,Ext[2]]]},\[Infinity]];

    ver=ver/.indrep;
(*Print[ver];*)
    lab[[All,2]]=Table[ll,{ll,1,Length[lab]}];
    lab[[All,1]]=lab[[All,1]]/.indrep;
    lab=Transpose[{lab[[All,1]],lab[[All,2]]}];

    If[FreeQ[vertList[[All,1]],lab],
      vertList=Append[vertList,{lab,(ver/.{FR$UV->0,FR$R2->1}),(ver/.{FR$R2->0,FR$UV->1})}];
      If[Length[lab]===2,UV$Wftlist=Append[UV$Wftlist,{lab[[1,1]],wft}]];
      ,
      pos=Position[vertList[[All,1]],lab][[1]];
      vertList[[pos,2]]=vertList[[pos,2]]+(ver/.{FR$UV->0,FR$R2->1});
      vertList[[pos,3]]=vertList[[pos,3]]+(ver/.{FR$UV->1,FR$R2->0});
      If[Length[lab]===2,
        pos=Position[UV$Wftlist[[All,1]],lab[[1,1]]][[1]];
        UV$Wftlist[[pos,2]]=UV$Wftlist[[pos,2]]+wft;
      ];
    ];
  ];
  vertList
]


(* ::Subsubtitle:: *)
(*ModelR2*)


R2vertlist::usage="generate the list of R2 vertices of a model and a generic model. If #3 is True, each contribution is multiplyed by IPL[list of the particles in the loop]. The last argument 
is the list of indices to be kept in IPL like generation."


ModelR2[mod_,gen_,lab_,kept_,qcd_]:=Block[{vertlist,genver,kkpl,topo,INTopo,Ampl,tmp,tmpCT,x,ind,indl,rl1,rlf,tdrl,a,b,temp,ckt,topoCT,INTopoCT,AmplCT,CTSol,
vrl,vrl2,momk,vertCTlist},
temp=SessionTime[];
  genver={{S},{F,F},{S,S},{V,V},{F,F,S},{F,F,V},{V,V,V},{V,V,S}(*,{V,S,S},{S,S,S},{V,V,V,V},{V,V,S,S},{S,S,S,S}*),{F,F,V,S}};
  vertlist={};
  vertCTlist={};
  CTSol={};
  tdrl={TensDot[ProjP][s1_,s2_]->ProjP[s1,s2],TensDot[ProjM][s1_,s2_]->ProjM[s1,s2],TensDot[Ga[mu_]][s1_,s2_]->Ga[mu,s1,s2]};
  rlf={TensDot[gm1:_Ga|_SlashedP,gm2:_Ga|_SlashedP,ProjP][Index[Spin,Ext[2]],Index[Spin,Ext[1]]]->-TensDot[gm2,gm1,ProjP][Index[Spin,Ext[1]],Index[Spin,Ext[2]]],
       TensDot[gm1:_Ga|_SlashedP,gm2:_Ga|_SlashedP,ProjM][Index[Spin,Ext[2]],Index[Spin,Ext[1]]]->-TensDot[gm2,gm1,ProjM][Index[Spin,Ext[1]],Index[Spin,Ext[2]]],
       TensDot[gm:_Ga|_SlashedP,ProjP][Index[Spin,Ext[2]],Index[Spin,Ext[1]]]->TensDot[gm,ProjM][Index[Spin,Ext[1]],Index[Spin,Ext[2]]],
       TensDot[gm:_Ga|_SlashedP,ProjM][Index[Spin,Ext[2]],Index[Spin,Ext[1]]]->TensDot[gm,ProjP][Index[Spin,Ext[1]],Index[Spin,Ext[2]]],
       ProjP[Index[Spin,Ext[2]],Index[Spin,Ext[1]]]->-ProjP[Index[Spin,Ext[1]],Index[Spin,Ext[2]]],
       ProjM[Index[Spin,Ext[2]],Index[Spin,Ext[1]]]->-ProjM[Index[Spin,Ext[1]],Index[Spin,Ext[2]]]};
  
  vrl={MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[3]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[3]],Index[Gluon,Ext[2]],Index[Gluon,Ext[4]]]->
       MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[3]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[2]],Index[Gluon,Ext[3]],Index[Gluon,Ext[4]]]+
       MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[3]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[4]],Index[Gluon,Ext[2]],Index[Gluon,Ext[3]]],
       MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[2]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[2]],Index[Gluon,Ext[3]],Index[Gluon,Ext[4]]]->
       MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[2]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[3]],Index[Gluon,Ext[2]],Index[Gluon,Ext[4]]]-
       MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[2]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[4]],Index[Gluon,Ext[2]],Index[Gluon,Ext[3]]],
       MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[4]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[4]],Index[Gluon,Ext[2]],Index[Gluon,Ext[3]]]->
       -MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[4]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[2]],Index[Gluon,Ext[3]],Index[Gluon,Ext[4]]]+
       MetricTensor[Index[Lorentz, Ext[1]], Index[Lorentz, Ext[4]]]SUF[Index[Gluon,Ext[1]],Index[Gluon,Ext[3]],Index[Gluon,Ext[2]],Index[Gluon,Ext[4]]]};

  For[kkpl=1,kkpl<=Length[genver]-0,kkpl++,Print[kkpl];
    rl1={MetricTensor->ME,ScalarProduct[FourMomentum[Incoming,a_],FourMomentum[Incoming,b_]]->SP[FourMomentum[a],FourMomentum[b]],ScalarProduct[FourMomentum[Incoming,a_],FourMomentum[Incoming,a_]]->
        SP[FourMomentum[a],FourMomentum[a]],SP[FourMomentum[Incoming,a_],FourMomentum[Incoming,b_]]->SP[FourMomentum[a],FourMomentum[b]],
        SP[FourMomentum[Incoming,a_],FourMomentum[Incoming,a_]]->SP[FourMomentum[a],FourMomentum[a]],LeviCivita[x__]->Eps[x],
	    FourVector[FourMomentum[Incoming,ind_],Index[Lorentz,Ext[indl_]]]:>FV[FourMomentum[ind],Index[Lorentz,Ext[indl]]]/;indl<=Length[genver[[kkpl]]],
	    FourVector[-FourMomentum[Incoming,ind_],Index[Lorentz,Ext[indl_]]]:>-FV[FourMomentum[ind],Index[Lorentz,Ext[indl]]]/;indl<=Length[genver[[kkpl]]],
        NonCommutative[ChiralityProjector[1]]->ProjP,NonCommutative[ChiralityProjector[-1]]->ProjM,NonCommutative[DiracSlash[FourMomentum[Incoming,x_]]]->
        SlashedP[FourMomentum[x]],NonCommutative[DiracMatrix[Index[Lorentz,x_]]]->Ga[Index[Lorentz,x]],SumOver[x__,External]->1,EL->ee,GS->gs};

    topo = CreateTopologies[1,Length[genver[[kkpl]]]->0,Adjacencies->{3,4},ExcludeTopologies->{Internal}];
    INTopo = InsertFields[topo,genver[[kkpl]]->{},Model->mod,GenericModel->gen,ExcludeParticles->Join[If[kkpl===13,{U},{}],{}],
                          ExcludeFieldPoints->If[kkpl>10,{FieldPoint[_][S,S,S],FieldPoint[_][S,V,V]},{}]];
    If[qcd,INTopo=DeleteCases[DeleteCases[DeleteCases[INTopo,FeynmanGraph[_,Classes==_][__?(FreeQ[#,Gluon]&)],\[Infinity]],
                              _->Insertions[Classes][],{3}],_->Insertions[Generic][],{1}];];
    Ampl= CreateFeynAmp[INTopo];
    For[ckt=1,ckt<=Length[Ampl],ckt++,Ampl[[ckt,4,2]]=Ampl[[ckt,4,2]]/.{aa_[bb__?(FreeQ[#,Colour]&),Index[Colour,cc_]]->aa[bb]};(*remove colour index in masses*)];
    Print["FA finished after "<>ToString[SessionTime[]-temp]];
    tmp=R2atClass[Ampl,INTopo,genver[[kkpl]],lab,kept,qcd];
    Print["R2atClass finished after "<>ToString[SessionTime[]-temp]];
    tmp=R2vertlist[Expand[Expand[Expand[tmp,_SWF],_PolarizationVector],FermionChain]];

	(*Counterterm amplitude*)
    topoCT=CreateCTTopologies[1,Length[genver[[kkpl]]]->0,Adjacencies->{3,4},ExcludeTopologies->{Internal}];
    INTopoCT = InsertFields[topoCT,genver[[kkpl]]->{},Model->mod,GenericModel->gen];
    AmplCT= CreateFeynAmp[INTopoCT];
    tmpCT=CTatClass[AmplCT,INTopoCT,genver[[kkpl]]];
    tmpCT=R2vertlist[Expand[Expand[Expand[tmpCT,_SWF],_PolarizationVector],FermionChain]];

    Print["R2vertlist finished after "<>ToString[SessionTime[]-temp]];
    tmp=Expand[Expand[tmp/.{SUNT->SUT,SUNF->SUF},_Colour],_Gluon];
    tmpCT=Expand[Expand[tmpCT/.{SUNT->SUT,SUNF->SUF},_Colour],_Gluon];
     Print["Expand finished after "<>ToString[SessionTime[]-temp]];
    (* Esthetic replacements*)
    vrl2={FourVector[FourMomentum[Incoming, a_Integer], Index[Lorentz, Ext[a_Integer]]]->-(Total[Table[FourVector[FourMomentum[Incoming,momk ], 
          Index[Lorentz, Ext[a]]],{momk,1,Length[genver[[kkpl]]]}]]-FourVector[FourMomentum[Incoming, a], Index[Lorentz, Ext[a]]])};
    tmp=tmp/.vrl/.vrl2/.SUT[Index[Colour, Ext[a_]], Index[Colour, Ext[b_]]]->IndexDelta[Index[Colour, Ext[a]], Index[Colour, Ext[b]]];
    tmpCT=tmpCT/.vrl/.vrl2/.SUT[Index[Colour, Ext[a_]], Index[Colour, Ext[b_]]]->IndexDelta[Index[Colour, Ext[a]], Index[Colour, Ext[b]]];
    Print["SU3 finished after "<>ToString[SessionTime[]-temp]];


    If[Length[tmp]>0,
      tmp[[All,2]]=I*(Simplify[ReplaceRepeated[#/.rl1/.tdrl/.rlf,SumOver[ind_?((FreeQ[#,Gluon]&&FreeQ[#,Colour])&),x_Integer]*a__:>Total[Table[(Times[a])/.{ind->kkpl},{kkpl,1,x}]]],TimeConstraint->10]&)/@tmp[[All,2]];

      tmp[[All,3]]=(ReplaceRepeated[#/.rl1/.tdrl/.rlf,SumOver[ind_?((FreeQ[#,Gluon]&&FreeQ[#,Colour])&),x_Integer]*a__:>Total[Table[(Times[a])/.{ind->kkpl},{kkpl,1,x}]]]&)/@tmp[[All,3]];
      tmpCT[[All,3]]=(ReplaceRepeated[#/.rl1/.tdrl/.rlf,SumOver[ind_?((FreeQ[#,Gluon]&&FreeQ[#,Colour])&),x_Integer]*a__:>Total[Table[(Times[a])/.{ind->kkpl},{kkpl,1,x}]]]&)/@tmpCT[[All,3]];

      If[Length[genver[[kkpl]]]===2,Print["before CTSol"];
        (*Only the finite piece in the renormalization constants*)
        CTSol=Join[CTSol,SolveDelta[tmpCT,({#1,#2,Simplify[#3-0*Coefficient[#3,inveps,1]*inveps/.inveps->1/FR$Eps]}&)@@@(tmp/.FR$Eps->1/inveps),genver[[kkpl]]]];

      ];
      (*Only the 1/Eps in the vertlist if not tadpole, minus sign because counter terms have opposite sign than amplitudes*)
      If[Length[genver[[kkpl]]]===1,
        tmp[[All,3]]=-I*(Simplify[#,TimeConstraint->10]&)/@tmp[[All,3]];,
        tmp[[All,3]]=If[Length[genver[[kkpl]]]===2,0*tmp[[All,3]],-I*(Simplify[Coefficient[#/.FR$Eps->1/inveps,inveps,1]/FR$Eps,TimeConstraint->10]&)/@tmp[[All,3]]];
      ];
      tmpCT[[All,3]]=I*(Simplify[#/.{Conjugate[FR$deltaZ[{x_,x_},b__]]->FR$deltaZ[{x,x},b]},TimeConstraint->10]&)/@tmpCT[[All,3]];Print["end lenght tmp not 0"];,(*Length[tmp]=0*)
Print["length tmp = 0"];
      If[Length[tmpCT]>0&&Length[genver[[kkpl]]]===2,
        tmpCT[[All,3]]=(ReplaceRepeated[#/.rl1/.tdrl/.rlf,SumOver[ind_?((FreeQ[#,Gluon]&&FreeQ[#,Colour])&),x_Integer]*a__:>Total[Table[(Times[a])/.{ind->kkpl},{kkpl,1,x}]]]&)/@tmpCT[[All,3]];
        CTSol=Join[CTSol,SolveDelta[tmpCT,{},genver[[kkpl]]]];
      ];
    ];(*end if Length[tmp]>0*)
    Print["Simplify finished after "<>ToString[SessionTime[]-temp]];
    
    vertlist=Join[vertlist,tmp];
    vertCTlist=Join[vertCTlist,tmpCT];
    Print[ToString[genver[[kkpl]]]<>" finished after "<>ToString[SessionTime[]-temp]];

  ];

Print[InputForm[CTSol]];

Print["before merge"];Print[Length[vertlist]];Print[Length[vertCTlist]];
(*Print[InputForm[vertCTlist]];*)
vertlist=MergeVertList[vertlist,({#[[1]],#[[2]],#[[3]]-If[Length[#[[1]]]===2,0,Coefficient[Coefficient[Refine[FR$Eps*#[[3]],Assumptions->{FR$Eps>=0,FR$IR>0}],FR$IR,0],FR$Eps,0]/FR$Eps]}&)/@(vertCTlist//.CTSol)]/.{FourMomentum[xx_]->xx};
Print["merging done"];
(*Print[InputForm[Cases[vertlist,{a_?(Not[FreeQ[#,9]]&&FreeQ[#,S[2]]&&FreeQ[#,S[3]]&),b__}]]];*)
vertlist/.{SUT->SUNT,SUF->SUNF}
]


(* ::Subsubtitle::Closed:: *)
(*WriteR2*)


WriteR2::usage="Write the list of the R2 vertices of the model and generic model in the file named #3.fr2. Two options are LabelInternal and KeptIndices"
LabelInternal::usage="Option of WriteR2. If True, IPL[list of internal particles] multiply each contribution. Default value is True"
KeptIndices::usage="Option of WriteR2. Contain the list of indices that have to be kept in IPL. Default value is {}"
QCDonly::usage="Option of WriteR2. Keep only the QCD corrections if set to True. Default value is False"


Options[WriteR2] = {LabelInternal -> True, KeptIndices->{},QCDonly->False};


WriteR2[mod_,gen_,out_,options___]:=Block[{totlist,outfile,LabInt,KeptInd,qcd,r2list,uvlist,InvEps,ipl1m,wsp,info},

  LabInt=LabelInternal/.{options}/.Options[WriteR2];
  KeptInd=KeptIndices/.{options}/.Options[WriteR2];
  qcd=QCDonly/.{options}/.Options[WriteR2];

  UV$Wftlist = {};

  totlist=ModelR2[mod,gen,LabInt,KeptInd,qcd];

  r2list = Transpose[Delete[Transpose[totlist],{3}]];
  uvlist = Transpose[Delete[Transpose[totlist],{2}]];
  r2list=DeleteCases[r2list,{a_,0},1];
  uvlist=DeleteCases[uvlist,{a_,0},1]/.FR$Eps->-2FR$Eps;(*change from d=4+eps to d=4-2eps*)

  uvlist=uvlist/.(Rule[-Append[#2,xx___],anti[#1]]&)@@@FR$ClassesTranslation/.(Rule[Append[#2,xx___],#1]&)@@@FR$ClassesTranslation;
  r2list=r2list/.(Rule[-Append[#2,xx___],anti[#1]]&)@@@FR$ClassesTranslation/.(Rule[Append[#2,xx___],#1]&)@@@FR$ClassesTranslation;

  FR$InteractionOrderPerturbativeExpansion=If[qcd,FR$InteractionOrderPerturbativeExpansion/.{{QCD,0}->{QCD,1}},FR$InteractionOrderPerturbativeExpansion/.{0->1}];
  
  wsp="                                                                                                                             ";
  info = (ToString[#1]<>" : "<>If[Head[#2]===List,StringDrop[StringJoin[(#<>", "&)/@#2],-2],#2]&)@@@FR$ModelInformation;
  outfile=out<>".nlo";
  Print["Writing R2 from "<>mod<>" in "<>outfile];
  OpenWrite[outfile];
  WriteString[outfile,"(****************************************************************************************************************)\n"];
  WriteString[outfile,"(* Model automalically generated by NLO-"<>StringTake[NLO$Version<>wsp,72]<>"*)\n"];
  WriteString[outfile,"(* date and time of generation : "<>StringTake[DateString[]<>wsp,79]<>"*)\n"];
  WriteString[outfile,"(* FeynRules model information : "<>StringTake[wsp,79]<>"*)\n"];
  (WriteString[outfile,"(*   "<>StringTake[#<>wsp,106]<>" *)\n"]&)/@info;
  WriteString[outfile,"(****************************************************************************************************************)\n"];
  WriteString[outfile,"\n"];
  WriteString[outfile,"(* FeynArts model files used to generated the output *)\n"];
  WriteString[outfile,"\n"];
  WriteString[outfile,"R2$Model = "<>mod<>";\n"];
  WriteString[outfile,"R2$GenericModel = "<>gen<>";\n"];
  WriteString[outfile,"\n"];
  WriteString[outfile,"(* Perturbed orders *)\n"];
  WriteString[outfile,"\n"];
  WriteString[outfile,"FR$InteractionOrderPerturbativeExpansion = "<>ToString[InputForm[FR$InteractionOrderPerturbativeExpansion]]<>";\n"];
  WriteString[outfile,"\n"];
  WriteString[outfile,"(*R2 vertices*)\n"];
  WriteString[outfile,"\n"];
  WriteString[outfile,"R2$vertlist = "];
  WriteString[outfile,ToString[InputForm[r2list]]<>";\n"];
  WriteString[outfile,"\n"];
  WriteString[outfile,"(*UV vertices*)\n"];
  WriteString[outfile,"\n"];
  WriteString[outfile,"UV$vertlist = "];
  WriteString[outfile,ToString[InputForm[uvlist]]<>";\n"];
  Close[outfile];
  Print["done"];
]


(* ::Subtitle::Closed:: *)
(*SUN*)


(*special unitary group simplification rules*)


(*IndexLessQ*)
ILQ[Ext[a_Integer],b_Integer]:=True;
ILQ[b_Integer,Ext[a_Integer]]:=False;
ILQ[a_Integer,b_Integer]:=a<b;
ILQ[Ext[a_Integer],Ext[b_Integer]]:=a<b;


(* ::Subsubtitle::Closed:: *)
(*SUT Trace*)


SUT[Index[Gluon,b_],Index[Gluon,a_],Index[Gluon,c_]]:=SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,b]]/;ILQ[a,c]&&ILQ[a,b]&&ILQ[c,b];
SUT[Index[Gluon,c_],Index[Gluon,b_],Index[Gluon,a_]]:=SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,b]]/;ILQ[a,c]&&ILQ[a,b]&&ILQ[c,b];
SUT[Index[Gluon,b_],Index[Gluon,a_],Index[Gluon,c_]]:=SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]-
I/2 SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]/;ILQ[a,c]&&ILQ[a,b]&&ILQ[b,c];
SUT[Index[Gluon,c_],Index[Gluon,b_],Index[Gluon,a_]]:=SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]-
I/2 SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]/;ILQ[a,c]&&ILQ[a,b]&&ILQ[b,c];
SUT[Index[Gluon,a_],Index[Gluon,c_],Index[Gluon,b_]]:=SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]-
I/2 SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]/;ILQ[a,c]&&ILQ[a,b]&&ILQ[b,c];

SUT[Index[Gluon,b_],Index[Gluon,a_],Index[Gluon,c_],Index[Gluon,d_]]:=SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,d],Index[Gluon,b]]/;ILQ[a,c]&&ILQ[a,b]&&ILQ[a,d];
SUT[Index[Gluon,d_],Index[Gluon,b_],Index[Gluon,a_],Index[Gluon,c_]]:=SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,d],Index[Gluon,b]]/;ILQ[a,c]&&ILQ[a,b]&&ILQ[a,d];
SUT[Index[Gluon,c_],Index[Gluon,d_],Index[Gluon,b_],Index[Gluon,a_]]:=SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,d],Index[Gluon,b]]/;ILQ[a,c]&&ILQ[a,b]&&ILQ[a,d];

SUT[Index[Gluon,a_],Index[Gluon,d_],Index[Gluon,c_],Index[Gluon,b_]]:=-1/2 SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c],Index[Gluon,d]]-
SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c],Index[Gluon,d]]+SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,d],Index[Gluon,c]]+
SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,d],Index[Gluon,b]]/;ILQ[a,b]&&ILQ[b,c]&&ILQ[c,d];

SUT[Index[Gluon,a_],Index[Gluon,b_],Index[Gluon,d_],Index[Gluon,c_]]:=1/2 SUF[Index[Gluon,a],Index[Gluon,c],Index[Gluon,b],Index[Gluon,d]]-
SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,d],Index[Gluon,b]]+SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,b],Index[Gluon,d]]+
SUT[Index[Gluon,a],Index[Gluon,d],Index[Gluon,b],Index[Gluon,c]]/;ILQ[a,b]&&ILQ[b,c]&&ILQ[c,d];

(*SUT[Index[Gluon,a_],Index[Gluon,c_],Index[Gluon,b_],Index[Gluon,d_]]:=-1/2 SUF[Index[Gluon,a],Index[Gluon,d],Index[Gluon,b],Index[Gluon,c]]-
SUT[Index[Gluon,a],Index[Gluon,d],Index[Gluon,b],Index[Gluon,c]]+SUT[Index[Gluon,a],Index[Gluon,d],Index[Gluon,c],Index[Gluon,b]]+
SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c],Index[Gluon,d]]/;ILQ[a,b]&&ILQ[b,c]&&ILQ[c,d];*)


(* ::Subsubtitle::Closed:: *)
(*SUT*)


SUT/:SUT[a__,Index[Colour,i_],Index[Colour,j_]]*SUT[b__,Index[Colour,j_],Index[Colour,k_]]*SumOver[Index[Colour,j_],n_]:=SUT[a,b,Index[Colour,i],Index[Colour,k]];
SUT/:SUT[a__,Index[Colour,i_],Index[Colour,j_]]*SUT[b__,Index[Colour,j_],Index[Colour,i_]]*SumOver[Index[Colour,j_],n_]:=SUT[a,b]/SumOver[Index[Colour,i],n];

SUT/:SUT[a__,Index[Colour,i_],Index[Colour,i_]]*SumOver[Index[Colour,i_],n_]:=SUT[a];

SUT/:SUT[a___,Index[Gluon,c_],Index[Gluon,c_],b___]*SumOver[Index[Gluon,c_],n2m1_Integer]:=n2m1/2/Sqrt[n2m1+1]*SUT[a,b];

SUT/:SumOver[Index[Gluon,a_],n2m1_Integer]*SUT[Index[Gluon,a_],Index[Gluon,a_],Index[Colour,i_],Index[Colour,l_]]:=(n2m1/2/Sqrt[n2m1+1])IndexDelta[Index[Colour,i],Index[Colour,l]];

SUT/:SumOver[Index[Gluon,a_],n2m1_Integer]*SUT[Index[Gluon,a_],Index[Gluon,b_],Index[Gluon,a_],Index[Colour,i_],Index[Colour,l_]]:= 
-1/2/(Sqrt[n2m1+1]) SUT[Index[Gluon,b],Index[Colour,i],Index[Colour,l]];

SUT/:SumOver[Index[Colour,j_],n_Integer]SUT[Index[Gluon,a_],Index[Gluon,a_],Index[Colour,i_],Index[Colour,l_]]:=
1/2(n-1/n) IndexDelta[Index[Colour,i],Index[Colour,l]]/SumOver[Index[Gluon,a],n^2-1];

SUT[Index[Gluon,a_],Index[Colour,i_],Index[Colour,i_]]:=0;

SUT[Index[Gluon,a_]]:=0;

SUT/:SUT[Index[Gluon,c_],Index[Gluon,d_]]:= IndexDelta[Index[Gluon,c],Index[Gluon,d]]/2;

SUT/:IndexDelta[Index[Colour, a_], Index[Colour, b_]]SumOver[Index[Colour, b_], n_]*SUT[h__,Index[Colour, b_],k___]:=SUT[h,Index[Colour, a],k];

SUT/:SUF[Index[Gluon,a_],Index[Gluon,b_],Index[Gluon,c_]]*SUT[Index[Gluon,b_],Index[Gluon,c_],Index[Colour,j_],Index[Colour,i_]]*SumOver[Index[Gluon,c_],n2m1_]:=
I Sqrt[n2m1+1]*SUT[Index[Gluon,a],Index[Colour,j],Index[Colour,i]]/2/(SumOver[Index[Gluon,b],n2m1]);

SUT/:SUF[Index[Gluon,a_],Index[Gluon,b_],Index[Gluon,c_]]*SUT[Index[Gluon,c_],Index[Gluon,b_],Index[Colour,i_],Index[Colour,j_]]*SumOver[Index[Gluon,c_],n2m1_]:=
-I Sqrt[n2m1+1]*SUT[Index[Gluon,a],Index[Colour,i],Index[Colour,j]]/2/(SumOver[Index[Gluon,b],n2m1]);

SUT/: SUT[Index[Gluon, m_],Index[Gluon,l_],Index[Gluon,k_],Index[Gluon,m_], Index[Colour, e_], Index[Colour, d_]]*SumOver[Index[Gluon, m_],n2m1_]:= 
  1/2(IndexDelta[Index[Colour, e],Index[Colour, d]]IndexDelta[Index[Gluon, l],Index[Gluon, k]]/2-
          1/Sqrt[n2m1+1]  SUT[Index[Gluon,l],Index[Gluon,k], Index[Colour, e], Index[Colour, d]]);


(* ::Subsubtitle::Closed:: *)
(*SUF*)


SUF[Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,a_]]:=Signature[{e,c,a}]*Signature[Sort[{e,c,a},ILQ]]*
SUF[Sequence@@Sort[{Index[Gluon,e],Index[Gluon,c],Index[Gluon,a]},(ILQ[#1[[2]],#2[[2]]]&)]]/;Not[OrderedQ[{e,c,a},ILQ]];


SUF[Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,b_]]:=- SUF[Index[Gluon,e],Index[Gluon,c],Index[Gluon,b],Index[Gluon,a]]/;ILQ[b,a]&&ILQ[e,c]&&ILQ[e,b];
SUF[Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,b_]]:=- SUF[Index[Gluon,c],Index[Gluon,e],Index[Gluon,a],Index[Gluon,b]]/;ILQ[c,e]&&ILQ[a,b]&&ILQ[c,a];
SUF[Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,b_]]:= SUF[Index[Gluon,c],Index[Gluon,e],Index[Gluon,b],Index[Gluon,a]]/;ILQ[c,e]&&ILQ[b,a]&&ILQ[c,b];
SUF[Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,b_]]:= SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,e],Index[Gluon,c]]/;ILQ[a,e]&&ILQ[e,c]&&ILQ[a,b];
SUF[Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,b_]]:=- SUF[Index[Gluon,b],Index[Gluon,a],Index[Gluon,e],Index[Gluon,c]]/;ILQ[b,e]&&ILQ[e,c]&&ILQ[b,a];
SUF[Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,b_]]:= SUF[Index[Gluon,b],Index[Gluon,a],Index[Gluon,c],Index[Gluon,e]]/;ILQ[b,c]&&ILQ[c,e]&&ILQ[b,a];
SUF[Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,b_]]:=- SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c],Index[Gluon,e]]/;ILQ[a,c]&&ILQ[c,e]&&ILQ[a,b];


SUF/:SumOver[Index[Gluon,a_],n2m1_Integer]SUF[Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,b_]]*
SUF[Index[Gluon,e_],Index[Gluon,a_],Index[Gluon,b_]]:=(Sqrt[n2m1+1])*IndexDelta[Index[Gluon,c],Index[Gluon,e]]/SumOver[Index[Gluon,b],n2m1];

SUF/:SumOver[Index[Gluon,a_],n2m1_Integer]SUF[Index[Gluon,c_],Index[Gluon,a_],Index[Gluon,e_],Index[Gluon,a_]]:=
(Sqrt[n2m1+1])*IndexDelta[Index[Gluon,c],Index[Gluon,e]];

SUF/:SumOver[Index[Gluon,a_],n2m1_Integer]*SUF[Index[Gluon,d_],Index[Gluon,a_],Index[Gluon,b_]]*SUF[Index[Gluon,e_],Index[Gluon,a_],Index[Gluon,c_]]*
SUF[Index[Gluon,f_],Index[Gluon,b_],Index[Gluon,c_]]:=-3/2 SUF[Index[Gluon,d],Index[Gluon,f],Index[Gluon,e]]/SumOver[Index[Gluon,b],n2m1]/SumOver[Index[Gluon,c],n2m1];

SUF/:SUF[Index[Gluon, c_], Index[Gluon, d_], Index[Gluon, e_]]*SUF[Index[Gluon, a_], Index[Gluon, d_], Index[Gluon, b_], Index[Gluon, e_]]*
  SumOver[Index[Gluon, d_], n2m1_Integer]:=3/2 SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]/SumOver[Index[Gluon, e], n2m1];

SUF/:SUF[Index[Gluon, c_], Index[Gluon, e_], Index[Gluon, d_]]*SUF[Index[Gluon, a_], Index[Gluon, d_], Index[Gluon, b_], Index[Gluon, e_]]*
  SumOver[Index[Gluon, d_], n2m1_Integer]:=-3/2 SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]/SumOver[Index[Gluon, e], n2m1];

SUF/:SUF[Index[Gluon, c_], Index[Gluon, d_], Index[Gluon, e_]]*SUF[Index[Gluon, a_], Index[Gluon, b_], Index[Gluon, d_], Index[Gluon, e_]]*
  SumOver[Index[Gluon, d_], n2m1_Integer]:= Sqrt[n2m1+1]SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c]]/SumOver[Index[Gluon, e], n2m1];

SUF/:SUF[Index[Gluon,a_],Index[Gluon,d_],Index[Gluon,e_],Index[Gluon,f_]]*SUF[Index[Gluon,b_],Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,f_]]*
SumOver[Index[Gluon,e_],n2m1_Integer]:= 3/2/SumOver[Index[Gluon,f],n2m1]SUF[Index[Gluon,a],Index[Gluon,d],Index[Gluon,b],Index[Gluon,c]];

SUF/:SUF[Index[Gluon,a_],Index[Gluon,d_],Index[Gluon,f_],Index[Gluon,e_]]*SUF[Index[Gluon,b_],Index[Gluon,e_],Index[Gluon,c_],Index[Gluon,f_]]*
SumOver[Index[Gluon,e_],n2m1_Integer]:=-3/2/SumOver[Index[Gluon,f],n2m1]SUF[Index[Gluon,a],Index[Gluon,d],Index[Gluon,b],Index[Gluon,c]];

SUF/:SUF[Index[Gluon,a_],Index[Gluon,d_],Index[Gluon,e_],Index[Gluon,f_]]*SUF[Index[Gluon,b_],Index[Gluon,c_],Index[Gluon,e_],Index[Gluon,f_]]*
SumOver[Index[Gluon,e_],n2m1_Integer]:= Sqrt[n2m1+1]/SumOver[Index[Gluon,f],n2m1]SUF[Index[Gluon,a],Index[Gluon,d],Index[Gluon,b],Index[Gluon,c]];

SUF/:SUF[Index[Gluon,c_],Index[Gluon,e_],Index[Gluon,f_]]*SUF[Index[Gluon,d_],Index[Gluon,e_],Index[Gluon,h_]]*
SUF[Index[Gluon,a_],Index[Gluon,b_],Index[Gluon,f_],Index[Gluon,h_]]*SumOver[Index[Gluon,e_],n2m1_Integer]:=-3/2*
SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,d],Index[Gluon,c]]/SumOver[Index[Gluon,f],n2m1]/SumOver[Index[Gluon,h],n2m1];

SUF/:SUF[Index[Gluon,a_],Index[Gluon,e_],Index[Gluon,f_]]*SUF[Index[Gluon,d_],Index[Gluon,e_],Index[Gluon,h_]]*
SUF[Index[Gluon,b_],Index[Gluon,f_],Index[Gluon,c_],Index[Gluon,h_]]*SumOver[Index[Gluon,e_],n2m1_Integer]:=(IndexDelta[Index[Gluon,a],Index[Gluon,c]]*
IndexDelta[Index[Gluon,d],Index[Gluon,b]]/2+IndexDelta[Index[Gluon,a],Index[Gluon,d]]IndexDelta[Index[Gluon,c],Index[Gluon,b]]/2+
IndexDelta[Index[Gluon,a],Index[Gluon,b]]*IndexDelta[Index[Gluon,c],Index[Gluon,d]]/2+(Sqrt[n2m1+1])*(SUT[Index[Gluon,a],Index[Gluon,d],Index[Gluon,c],Index[Gluon,b]]+
SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c],Index[Gluon,d]]))/(SumOver[Index[Gluon,f],n2m1]*SumOver[Index[Gluon,h],n2m1]);

SUF/:SUF[Index[Gluon, a_], Index[Gluon, e_], Index[Gluon, b_], Index[Gluon, f_]]*
SUF[Index[Gluon, c_], Index[Gluon, e_], Index[Gluon, d_], Index[Gluon, f_]]*SumOver[Index[Gluon, e_], n2m1_Integer]:=
(IndexDelta[Index[Gluon,a],Index[Gluon,c]]IndexDelta[Index[Gluon,d],Index[Gluon,b]]/2+
IndexDelta[Index[Gluon,a],Index[Gluon,d]]IndexDelta[Index[Gluon,c],Index[Gluon,b]]/2+IndexDelta[Index[Gluon,a],Index[Gluon,b]]IndexDelta[Index[Gluon,c],Index[Gluon,d]]/2+
(Sqrt[n2m1+1])*(SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,d],Index[Gluon,b]]+SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,d],Index[Gluon,c]]))/
SumOver[Index[Gluon, f], n2m1];

SUF/:SUF[Index[Gluon,a_],Index[Gluon,f_],Index[Gluon,b_],Index[Gluon,e_]]*SUF[Index[Gluon,c_],Index[Gluon,e_],Index[Gluon,d_],Index[Gluon,f_]]*
SumOver[Index[Gluon,e_], n2m1_Integer]:=(IndexDelta[Index[Gluon,a],Index[Gluon,c]]IndexDelta[Index[Gluon,d],Index[Gluon,b]]/2+IndexDelta[Index[Gluon,a],Index[Gluon,d]]*
IndexDelta[Index[Gluon,c],Index[Gluon,b]]/2+IndexDelta[Index[Gluon,a],Index[Gluon,b]]IndexDelta[Index[Gluon,c],Index[Gluon,d]]/2+(Sqrt[n2m1+1])*
(SUT[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c],Index[Gluon,d]]+SUT[Index[Gluon,a],Index[Gluon,d],Index[Gluon,c],Index[Gluon,b]]))/SumOver[Index[Gluon,f], n2m1];

SUF/:SUF[Index[Gluon,a_],Index[Gluon,e_],Index[Gluon,f_]]*SUF[Index[Gluon,b_],Index[Gluon,g_],Index[Gluon,h_]]*SUF[Index[Gluon,c_],Index[Gluon,e_],Index[Gluon,g_]]*
SUF[Index[Gluon,d_],Index[Gluon,f_],Index[Gluon,h_]]*SumOver[Index[Gluon,e_],n2m1_Integer]:=(IndexDelta[Index[Gluon,a],Index[Gluon,c]]*
IndexDelta[Index[Gluon,d],Index[Gluon,b]]/2+IndexDelta[Index[Gluon,a],Index[Gluon,d]]IndexDelta[Index[Gluon,c],Index[Gluon,b]]/2+
IndexDelta[Index[Gluon,a],Index[Gluon,b]]IndexDelta[Index[Gluon,c],Index[Gluon,d]]/2+(Sqrt[n2m1+1])*(SUT[Index[Gluon,a],Index[Gluon,d],Index[Gluon,b],Index[Gluon,c]]+
SUT[Index[Gluon,a],Index[Gluon,c],Index[Gluon,b],Index[Gluon,d]]))/(SumOver[Index[Gluon,f],n2m1]*SumOver[Index[Gluon,g],n2m1]*SumOver[Index[Gluon,h],n2m1]);

SUF/:SUF[Index[Gluon, a_], Index[Gluon, c_], Index[Gluon,b_], Index[Gluon, d_]]*SUT[Index[Gluon, c_],Index[Gluon, d_], Index[Colour, i_], Index[Colour, j_]]*
SumOver[Index[Gluon, d_], n2m1_]:=(Sqrt[n2m1+1]/2 SUT[Index[Gluon, a],Index[Gluon, b], Index[Colour, i], Index[Colour, j]]+
IndexDelta[Index[Gluon, a],Index[Gluon, b]]IndexDelta[Index[Colour, i], Index[Colour, j]]/4)/SumOver[Index[Gluon, c], n2m1];

SUF/:SUF[Index[Gluon, a_], Index[Gluon, b_], Index[Gluon, c_], Index[Gluon, d_]]*SUT[Index[Gluon, d_], Index[Gluon, b_],Index[Colour, i_], Index[Colour, j_]]*
SumOver[Index[Gluon, d_], n2m1_]:=(Sqrt[n2m1+1]/2 SUT[Index[Gluon, c], Index[Gluon, a],Index[Colour, i], Index[Colour, j]]+
IndexDelta[Index[Gluon, c], Index[Gluon, a]]IndexDelta[Index[Colour, i], Index[Colour, j]]/4)/SumOver[Index[Gluon, b], n2m1];

SUF/:SUF[Index[Gluon,a_],Index[Gluon,e_],Index[Gluon,b_]]*SUF[Index[Gluon,c_],Index[Gluon,e_],Index[Gluon,d_]]*SumOver[Index[Gluon,e_],n2m1_]:=
 SUF[Index[Gluon,a],Index[Gluon,b],Index[Gluon,c],Index[Gluon,d]];

SUF/:SumOver[Index[Gluon, a_], n2m1_Integer]*SUF[Index[Gluon, c_], Index[Gluon, a_], Index[Gluon, b_]]* 
SUT[Index[Gluon, a_], Index[Gluon, b_], Index[Gluon, d_], Index[Colour,i_],Index[Colour,j_]]:= 
 Sqrt[n2m1+1] I SUT[Index[Gluon, c],Index[Gluon, d],Index[Colour,i],Index[Colour,j]]/2/SumOver[Index[Gluon, b], n2m1];

SUF/:SumOver[Index[Gluon, a_], n2m1_Integer]*SUF[Index[Gluon, c_], Index[Gluon, a_], Index[Gluon, b_]]* 
SUT[Index[Gluon, d_], Index[Gluon, b_],Index[Gluon, a_], Index[Colour,i_], Index[Colour,j_]]:= 
- Sqrt[n2m1+1] I SUT[Index[Gluon, c],Index[Gluon, d],Index[Colour,i],Index[Colour,j]]/2/SumOver[Index[Gluon, b], n2m1];

SUF/:SumOver[Index[Gluon, b_],n2m1_Integer]SUF[Index[Gluon, a_],Index[Gluon, b_],Index[Gluon, c_]]*
SUT[Index[Gluon, b_],Index[Gluon, d_],Index[Gluon, c_],Index[Colour,i_],Index[Colour,j_]]:= -I  IndexDelta[Index[Gluon, a],Index[Gluon, d]]*
IndexDelta[Index[Colour,i],Index[Colour,j]]/4/SumOver[Index[Gluon, c],n2m1];

SUF/:SumOver[Index[Gluon, b_],n2m1_Integer]SUF[Index[Gluon, a_],Index[Gluon, b_],Index[Gluon, c_]]*
SUT[Index[Gluon, c_],Index[Gluon, d_],Index[Gluon, b_],Index[Colour,i_],Index[Colour,j_]]:= I IndexDelta[Index[Gluon, a],Index[Gluon, d]]*
IndexDelta[Index[Colour,i],Index[Colour,j]]/4/SumOver[Index[Gluon, c],n2m1];


(* ::Subtitle:: *)
(*Protection*)


Protect[ME,SUF,SUT,SUNTr,GetR2,R2Tadpoles,R2BubblesF,R2BubblesQ2,R2FTriangles,R2BTriangles,R2FBTriangles,R2BBoxes,R2FBoxes,ILQ,WriteR2,ModelR2,PartOrder,R2atClass,
SP,FCh,DTr,FV,LeviCivita,M$dim,red2v,red4v,FR$UV,FR$R2,FR$MU,FR$Eps,MergeVertList,UVwfatClass,AddWf];
