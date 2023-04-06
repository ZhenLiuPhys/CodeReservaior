#define s8_cxx
#include "s8.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void s8::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L s8.C
//      Root > s8 t
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
   Double_t deta;

   mytree->Branch("deta",&deta,"deta/D");

   for (Long64_t jentry=0; jentry<nentries;jentry++) {



      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      deta=0;
      Int_t ng=0;

      for(Int_t i = 0 ; i<Particle_; i++){
	Double_t eta1, eta2;
	
	if(Particle_Status[i]==1){
	  ng++;
	  if(ng==1) eta1=Particle_Eta[i];
	  else if(ng==2)eta2=Particle_Eta[i];
	}
      }

	if(ng==2)deta=eta1-eta2;
	else deta=-100;

      mytree->Fill();

   }

   mytree->Write();
   output->Close();
}
