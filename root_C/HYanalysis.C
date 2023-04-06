#define HYanalysis_cxx
#include "HYanalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

inline Double_t deltaphi(Double_t phi1, Double_t phi2){
  Double_t result = TMath::Abs(phi1-phi2);
  if(result>TMath::Pi()) result=2*TMath::Pi()-result;
  return result;
}

inline Double_t deltaR(TLorentzVector v1, TLorentzVector v2){
  Double_t result;
  result=TMath::Sqrt(TMath::Power(deltaphi(v1.Phi(),v2.Phi()),2) + TMath::Power(v1.Eta()-v2.Eta(),2));
  return result;
}

void HYanalysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L eezhDim6.C
//      Root > eezhDim6 t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;


   TFile *output = new TFile("output.root","recreate");
   TTree *mytree = new TTree("mytree","mytree");
   Double_t hpt,b1pt,b2pt,b3pt,b4pt,a1pt,a2pt,db12,db34,da12,j1pt,j2pt,wpt,lbpt,ljpt,lapt;

   mytree->Branch("hpt",&hpt,"hpt/D");
   mytree->Branch("b1pt",&b1pt,"b1pt/D");
   mytree->Branch("b2pt",&b2pt,"b2pt/D");
   mytree->Branch("b3pt",&b3pt,"b3pt/D");
   mytree->Branch("b4pt",&b4pt,"b4pt/D");
   mytree->Branch("a1pt",&a1pt,"a1pt/D");
   mytree->Branch("a2pt",&a2pt,"a2pt/D");
   mytree->Branch("db12",&db12,"db12/D");
   mytree->Branch("db34",&db34,"db34/D");
   mytree->Branch("da12",&da12,"da12/D");
   mytree->Branch("j1pt",&j1pt,"j1pt/D");
   mytree->Branch("j2pt",&j2pt,"j2pt/D");
   mytree->Branch("wpt",&wpt,"wpt/D");
   mytree->Branch("lbpt",&lbpt,"lbpt/D");
   mytree->Branch("ljpt",&ljpt,"ljpt/D");
   mytree->Branch("lapt",&lapt,"lapt/D");

Int_t na, nbj, b1mother, b2mother, b3mother, b4mother, nj;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

na=0;
nbj=0;
nj=0;
b1mother=0;
b2mother=1;
b3mother=2;
b4mother=4;
hpt=0;
b1pt=0;
b2pt=0;
b3pt=0;
b4pt=0;
a1pt=0;
a2pt=0;
db12=0;
db34=0;
da12=0;
j1pt=0;
j2pt=0;
wpt=0;
lbpt=0;
ljpt=0;
lapt=0;
      TLorentzVector a1, a2, b1, b2, b3, b4;

      for(Int_t i = 0 ; i<Particle_; i++){
	
	if(Particle_Status[i]==3&&Particle_PID[i]==9000025){
	  na++;
	  if(na==1) {
		a1.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
		a1pt=Particle_PT[i];
		if(Particle_PT[i]>=lapt) lapt=Particle_PT[i];}
	  else if(na==2){
		a2.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
		a2pt=Particle_PT[i];
		if(Particle_PT[i]>=lapt) lapt=Particle_PT[i];}
	  else{
		cout << "something not right about A\n";}
	}

	if(Particle_Status[i]==3&&TMath::Abs(Particle_PID[i])==5){
	  nbj++;
	  if(nbj==1) {
		b1.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
		b1mother=Particle_M1[i];
		b1pt=Particle_PT[i];
		if(Particle_PT[i]>=lbpt) lbpt=Particle_PT[i];}
	  else if(nbj==2){
		b2.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
		b2mother=Particle_M1[i];
		b2pt=Particle_PT[i];
		if(Particle_PT[i]>=lbpt) lbpt=Particle_PT[i];}
	  else if(nbj==3){
		b3.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
		b3mother=Particle_M1[i];
		b3pt=Particle_PT[i];
		if(Particle_PT[i]>=lbpt) lbpt=Particle_PT[i];}
	  else if(nbj==4){
		b4.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
		b4mother=Particle_M1[i];
		b4pt=Particle_PT[i];
		if(Particle_PT[i]>=lbpt) lbpt=Particle_PT[i];}
	  else{
		cout << "something not right about bjets\n";}
	}
	else if(Particle_Status[i]==3&&Particle_PID[i]==25){
		hpt=Particle_PT[i];
	}
	else if(Particle_Status[i]==3&&TMath::Abs(Particle_PID[i])==24){
		wpt=Particle_PT[i];
	}
	else if(Particle_Status[i]==3&&(TMath::Abs(Particle_PID[i])<=4||Particle_PID[i]==21)){
	  nj++;
	  if(nj==1) {j1pt=Particle_PT[i];
		if(j1pt>=ljpt) ljpt=j1pt;}
	  else if(nj==2) {j2pt=Particle_PT[i];
		if(j2pt>=ljpt) ljpt=j2pt;}
	  else {cout << "something not right about jets\n";}
	}
	//else {cout << "additional particles not so expected\n";}
      }

	if(na==2&&nbj==4){
		da12=a1.DeltaR(a2);
		//if(a2.DeltaR(a1)<da12) da12=a2.DeltaR(a1);
		if(b1mother==b2mother&&b3mother==b4mother){
			if((b1+b2).Pt()>(b3+b4).Pt()){
			db12=b1.DeltaR(b2);
			db34=b3.DeltaR(b4);}
			else {
			db34=b1.DeltaR(b2);
			db12=b3.DeltaR(b4);}
		}
		else if(b1mother==b3mother&&b2mother==b4mother){
			if((b1+b3).Pt()>(b2+b4).Pt()){
			db12=b1.DeltaR(b3);
			db34=b2.DeltaR(b4);}
			else {
			db34=b1.DeltaR(b3);
			db12=b2.DeltaR(b4);}
		}
		else if(b1mother==b4mother&&b2mother==b3mother){
			if((b1+b4).Pt()>(b3+b2).Pt()){
			db12=b1.DeltaR(b4);
			db34=b3.DeltaR(b2);}
			else {
			db34=b1.DeltaR(b4);
			db12=b3.DeltaR(b2);}
		}
		else {
		cout << "something not right about pairing bjets\n";}
	}

      mytree->Fill();

   }

   mytree->Write();
   output->Close();
}
