#define dzp2b_cxx
#include "dzp2b.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TEnv.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include "TObject.h"
#include "TMath.h"
#include "TRandom.h"
#include "TThread.h"
#include "TChain.h"
#include "TProof.h"
#include "TProofLog.h"
#include "TProofMgr.h"
#include "TClonesArray.h"
#include "TH1D.h"

#include "TLorentzVector.h"
#include "TRandom.h"

void dzp2b::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L dzp2m.C
//      Root > dzp2m t
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
//   gSystem->Load("/scratch/novar/trunk/libDelphes");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Int_t n_events, i, bflag, first_b, second_b;
   n_events = 0;

   Double_t y_ll, cosstar_lzp, m_eff, mpt, y_ll_c, cosstar_lzp_b, bpt, bpt2;	

   Long64_t nbytes = 0, nb = 0;

   TFile *output = new TFile("output.root","recreate");
   TTree *mytree = new TTree("mytree","mytree");

   mytree->Branch("n_events",&n_events,"n_events/I");
   mytree->Branch("y_ll",&y_ll, "y_ll/D");
   mytree->Branch("y_ll_c",&y_ll_c, "y_ll_c/D");
   mytree->Branch("cosstar_lzp",&cosstar_lzp, "cosstar_lzp/D");
   mytree->Branch("cosstar_lzp_b",&cosstar_lzp_b, "cosstar_lzp_b/D");
   mytree->Branch("m_eff",&m_eff,"m_eff/D");
   mytree->Branch("bflag",&bflag,"bflag/I");
   mytree->Branch("bpt",&bpt,"bpt/D");
   mytree->Branch("bpt2",&bpt2,"bpt2/D");

//   mytree->Branch("y_ll_c",&y_ll_c, "y_ll_c/D");
//   mytree->Branch("cosstar_lzp_c",&cosstar_lzp_c, "cosstar_lzp_c/D");

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

	y_ll = 0;
	y_ll_c = 0;
	cosstar_lzp = 0;
	cosstar_lzp_b = 0;
	m_eff = 0;
	bflag = 0;
	Int_t n_jet = 0;
	bpt = 0;
	bpt2 = 0;
	first_b = 0;
	second_b = 0;

	TLorentzVector jet1, jet1_cm;
	TLorentzVector jet_all;
	TLorentzVector jet[10];
	TLorentzVector bjet, bjet_cm;

	if(Jet_size>=2){
		for(i=0; i<Jet_size; i++){
			if(Jet_BTag[i]==1&&Jet_PT[i]>=100) {
				if(bflag==0) first_b = i;
				if(bflag==1) second_b = i;
				bflag++;
			}
			if(Jet_PT[i]>=100) n_jet++;
		}
		if(bflag>=1){
			for(i=0; i<n_jet; i++){
	     			jet[i].SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
				jet_all = jet_all + jet[i];
				if(i == first_b) bjet=jet[i];
				if(i == second_b && second_b != 0) bpt2=jet[i].Pt(); 
			}
			n_events++;
//			jet_all = jet[0]+jet[1]+jet[2];
//			lepton_all = lepton1 + lepton2;

			y_ll = TMath::Abs(TMath::Log((jet_all.E()-jet_all.Pz())/(jet_all.E()+jet_all.Pz())))/2.;
     			y_ll_c = -TMath::Log((jet_all.E()-jet_all.Pz())/(jet_all.E()+jet_all.Pz()))/2.;
      
			jet1 = jet[0];
			jet1_cm = jet1;
			bjet_cm = bjet;

      			jet1_cm.Boost( - jet_all.BoostVector());
			bjet_cm.Boost( - jet_all.BoostVector());

      			cosstar_lzp = jet1_cm.CosTheta()*TMath::Abs(jet_all.Pz())/jet_all.Pz();
		      	cosstar_lzp_b = bjet_cm.CosTheta()*TMath::Abs(jet_all.Pz())/jet_all.Pz();

       			m_eff = jet_all.Mag();
			bpt = bjet.Pt();

        		mytree->Fill();
		}
	}

      // if (Cut(ientry) < 0) continue;
   }
 mytree->Write();
 output->Close();
}
