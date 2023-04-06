#define mmh_cxx
#include "mmh.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <vector>

inline TLorentzVector psmearing(TLorentzVector v, TRandom2 *randgen, Double_t a, Double_t b){
     TLorentzVector v_s;
     v_s.SetPx(randgen->Gaus(1,a*v.Px()+b)*v.Px());
     v_s.SetPy(randgen->Gaus(1,a*v.Py()+b)*v.Py());
     v_s.SetPz(randgen->Gaus(1,a*v.Pz()+b)*v.Pz());
     v_s.SetE(TMath::Sqrt(v_s.Px()**2+v_s.Py()**2+v_s.Pz()**2));
     return v_s;
}

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

void mmh::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L mmh.C
//      Root > mmh t
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
   TRandom2 *randgen = new TRandom2();

   Double_t m_dl;
   Double_t pt_dl;
   Double_t m_recoil;   

   Double_t m_dl_s;
   Double_t pt_dl_s;
   Double_t m_recoil_s;  

   mytree->Branch("m_dl",&m_dl,"m_dl/D");
   mytree->Branch("pt_dl",&pt_dl,"pt_dl/D");
   mytree->Branch("m_recoil",&m_recoil,"m_recoil/D");

   mytree->Branch("m_dl_s",&m_dl_s,"m_dl_s/D");
   mytree->Branch("pt_dl_s",&pt_dl_s,"pt_dl_s/D");
   mytree->Branch("m_recoil_s",&m_recoil_s,"m_recoil_s/D");

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      m_dl = 0;
      pt_dl = 0;
      m_recoil = 0;

      m_dl_s = 0;
      pt_dl_s = 0;
      m_recoil_s = 0;
      

      TLorentzVector ele0,pos0,mup,mum,higgs;

     for(Int_t i = 0; i<Particle_;i++){
          if(Particle_Status[i]==-1){
               if(Particle_PID[i]==-11){
                    pos0.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
               }
               else if(Particle_PID[i]==11){
                    ele0.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
               }          
          }
    	     else if(Particle_Status[i]==1){
	          if(Particle_PID[i]==13){
	               mup.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	               }
	          else if(Particle_PID[i]==-13){
	               mum.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	               }
	          else if(Particle_PID[i]==25)
                    higgs.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	     }
      }

      m_dl = (mup+mum).Mag();
      pt_dl = (mup+mum).Pt();
      
      TLorentzVector recoil;
      recoil.SetPxPyPzE(0.-mup.Px()-mum.Px(),0.-mup.Py()-mum.Py(),0.-mup.Pz()-mum.Pz(),250.0.-mup.E()-mum.E());

      m_recoil = recoil.Mag();

//Begin Smearing according to Fig6.3, Pg.72 of Li's Thesis. Eq.(3.1) Pg28

      TLorentzVector mup_s,mum_s,higgs_s;

      mup_s=psmearing(mup,randgen,2.5E-5,8E-4);
      mum_s=psmearing(mum,randgen,2.5E-5,8E-4);
//      higgs_s=psmearing(higgs,randgen,2.5E-5,8E-4);

      m_dl_s = (mup_s + mum_s).Mag();
      pt_dl_s = (mup_s + mum_s).Pt();

      TLorentzVector recoil_s;
      recoil_s.SetPxPyPzE(0.-mup_s.Px()-mum_s.Px(),0.-mup_s.Py()-mum_s.Py(),0.-mup_s.Pz()-mum_s.Pz(),250.0.-mup_s.E()-mum_s.E());

      m_recoil_s = recoil_s.Mag();      

      // if (Cut(ientry) < 0) continue;

      mytree->Fill();
   }

 mytree->Write();
 output->Close();
}
