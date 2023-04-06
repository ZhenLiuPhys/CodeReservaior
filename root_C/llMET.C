#define llMET_cxx
// The class definition in llMET.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("llMET.C")
// Root > T->Process("llMET.C","some options")
// Root > T->Process("llMET.C+")
//

#include "llMET.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <TH2.h>
#include <TStyle.h>
#include "TLorentzVector.h"


int r1 = 37851287;
int r2 = 890628713;
int r3 = 78257863;
int r4 = 197536741;
const double MZ = 91.1987;

//void sort_by_pt(vector<TLorentzVector> & vec);
Double_t get_MT(TLorentzVector visible, TLorentzVector invisible);
//TLorentzVector get_missing(TLorentzVector dil, TLorentzVector inv);

TFile *f1;
TTree *mytree;

int total = 0;

   Int_t n_jets;//number of jets
   Int_t n_leptons;//number of leptons
   Int_t n_b;//number of tagged b quark in the pythia event file
   Int_t n_a;//number of photons
   Float_t LEPPT;
   Double_t DRTLEP;
   Double_t MTLEP;
   Double_t MTALL;
   Double_t MLTAU;
   Double_t MET;
   Float_t JETPT[20];
   Float_t BJETPT[20];

void llMET::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
   
   f1 = new TFile("output.root", "recreate");
   mytree = new TTree("mytree", "ino variables");

   mytree->Branch("n_jets",&n_jets,"n_jets/I");
   mytree->Branch("n_leptons",&n_leptons,"n_leptons/I");
   mytree->Branch("n_b",&n_b,"n_b/I");
   mytree->Branch("n_a",&n_a,"n_a/I");

   mytree->Branch("LEPPT",&LEPPT,"LEPPT/F");
   mytree->Branch("DRTLEP",&DRTLEP,"DRTLEP/D");
   mytree->Branch("MTLEP",&MTLEP,"MTLEP/D");
   mytree->Branch("MTALL",&MTALL,"MTALL/D");
   mytree->Branch("MLTAU",&MLTAU,"MLTAU/D");
   mytree->Branch("MET",&MET,"MET/D");

}

void llMET::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t llMET::Process(Long64_t entry)
{
  fChain->GetTree()->GetEntry(entry);
  total++;

	
  vector<TLorentzVector> leptons;
  vector<TLorentzVector> jets;
  vector<TLorentzVector> bjets;  
  vector<TLorentzVector> taus;
  vector<TLorentzVector> photons;

  TLorentzVector temp(0,0,0,0);
  TLorentzVector dilep(0,0,0,0);
  TLorentzVector invisible(0,0,0,0);
  cout<< Electron_size << " " << Muon_size << " " << Photon_size << " " << Jet_size << "@" << entry << "\n";

  for (int i = 0; i < Electron_size && Electron_size >=2 ; i++){
    temp.SetPtEtaPhiM(Electron_PT[i],Electron_Eta[i], Electron_Phi[i], 0.0);
    if (temp.Pt() > 20.0 && fabs(temp.Eta() ) < 2.5){
      leptons.push_back(temp);
    }
  }

  for (int i = 0; i < Muon_size && Muon_size >=2 ; i++){
    temp.SetPtEtaPhiM(Muon_PT[i],Muon_Eta[i], Muon_Phi[i], 0.1);
    if (temp.Pt() > 20.0 && fabs(temp.Eta() ) < 2.5){
      leptons.push_back(temp);
    }
  }

  for (int i = 0; i < Photon_size; i++){
    temp.SetPtEtaPhiM(Photon_PT[i],Photon_Eta[i], Photon_Phi[i], 0.0);
    if (temp.Pt() > 20.0 && fabs(temp.Eta() ) < 2.5){
      photons.push_back(temp);
    }
  }

  for (int i = 0; i < Jet_size && i < 20; i++){
    temp.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i], Jet_Phi[i],Jet_Mass[i]);
    if (Jet_TauTag[i] > 0.0){
      if (temp.Pt() > 20.0 && fabs(temp.Eta() ) < 2.5 ){
	taus.push_back(temp);
      } 
    } else if (Jet_BTag[i] > 0.0){
      if (temp.Pt() > 20.0 && fabs(temp.Eta() ) < 2.5){
	bjets.push_back(temp);
      }
    } else {
	if (temp.Pt() > 30.0 && fabs(temp.Eta() ) < 2.8 ){
	jets.push_back(temp);
	}
    }
  }
      
  invisible.SetPtEtaPhiM(MissingET_MET[0], MissingET_Eta[0], MissingET_Phi[0],0.0);


//  sort_by_pt(taus);
//  sort_by_pt(leptons);
//  sort_by_pt(jets);
//  sort_by_pt(bjets);


  if ( (Electron_size > 1 || Muon_size > 1) && leptons.size() > 1){

    for (unsigned int i = 0; i < jets.size(); i++){
      JETPT[i] = jets[i].Pt();
    }

    for (unsigned int i = 0; i < bjets.size(); i++){
      BJETPT[i] = bjets[i].Pt();
    }
    
    n_jets=Jet_size;
    n_b=bjets.size();
    n_a=0;
    n_leptons=Electron_size + Muon_size;
    
    LEPPT = leptons.at(0).Pt();
    DRTLEP = leptons[0].DeltaR(leptons[1]);
    
    MTLEP = get_MT(leptons[0], invisible);

    dilep = leptons[0] + leptons[1];

    MTALL = get_MT(invisible, dilep);
    
    MLTAU = dilep.M();
    MET = invisible.Pt();
    mytree->Fill();
  }
 
  return kTRUE;
}

void llMET::SlaveTerminate()
{

}

void llMET::Terminate()
{
  cout << total << endl;
  f1->cd();
  mytree->Write();
  f1->Close(); 
}

double get_MT(TLorentzVector visible, TLorentzVector invisible){

  // NB assumes that the invisible system has no mass

  double ETvis = sqrt(visible.M()*visible.M() + visible.Pt()*visible.Pt());
  double MTSQ = visible.M()*visible.M() + 2.0*(invisible.Pt()*ETvis - visible.Px()*invisible.Px() - visible.Py()*invisible.Py() );
  return TMath::Sqrt(MTSQ);

}
