#define stau_analysis_cxx
#include "stau_analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double get_MT(TLorentzVector visible, TLorentzVector invisible){

  // NB assumes that the invisible system has no mass
  // ZL probably different from V. Barger's definition. Be carefull with this function

  double ETvis = sqrt(visible.M()*visible.M() + visible.Pt()*visible.Pt());
  double MTSQ = visible.M()*visible.M() + 2.0*(invisible.Pt()*ETvis - visible.Px()*invisible.Px() - visible.Py()*invisible.Py() );
  return TMath::Sqrt(MTSQ);

}

void sort_by_pt(vector<TLorentzVector> & vec){

  TLorentzVector swap;
  double ptmax;

  if (vec.size() > 0){

  for (unsigned i = 0; i < vec.size()-1; i++){

    ptmax = vec[i].Pt();

    for (unsigned j = i+1; j < vec.size(); j++){

      if (vec[j].Pt() > ptmax){

	ptmax = vec[j].Pt();
	swap = vec[i];
	vec[i] = vec[j];
	vec[j] = swap;
      }
    }
  }
  }

  return;
}

void stau_analysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L stau_analysis.C
//      root> stau_analysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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

   Int_t n_jets;//number of jets
   Int_t n_leptons;//number of leptons
   Int_t n_b;//number of tagged b quark in the pythia event file
   Int_t n_a;//number of photons
   Int_t n_ltau;
   Int_t n_2tau;
   Int_t n_bbltau;
   Int_t n_bb2tau;
   Float_t LEPPT;
   Double_t DRTLEP;
   Double_t MTLEP;
   Double_t MTALL;
   Double_t MLTAU;
   Double_t MET;
   Double_t MCT;
   Double_t DLPT;
   Double_t MT, MTZZ, PTmissjet, dphiMETjet, dphiMETll, etadilep; 
   Float_t JETPT[20];
   Float_t BJETPT[20];

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *output = new TFile("output.root","recreate");
   TTree *mytree = new TTree("mytree","mytree");

   mytree->Branch("n_jets",&n_jets,"n_jets/I");
   mytree->Branch("n_leptons",&n_leptons,"n_leptons/I");
   mytree->Branch("n_b",&n_b,"n_b/I");
   mytree->Branch("n_a",&n_a,"n_a/I");

   mytree->Branch("LEPPT",&LEPPT,"LEPPT/F");
   mytree->Branch("DRTLEP",&DRTLEP,"DRTLEP/D");
   mytree->Branch("MTLEP",&MTLEP,"MTLEP/D");
   mytree->Branch("MTALL",&MTALL,"MTALL/D");
   mytree->Branch("MLTAU",&MLTAU,"MLTAU/D");
   mytree->Branch("MCT",&MCT,"MCT/D");
   mytree->Branch("MET",&MET,"MET/D");
   mytree->Branch("DLPT",&DLPT,"DLPT/D");
   mytree->Branch("MT",&MT,"MT/D");
   mytree->Branch("MTZZ",&MTZZ,"MTZZ/D");
   mytree->Branch("PTmissjet",&PTmissjet,"PTmissjet/D");
   mytree->Branch("dphiMETjet",&dphiMETjet,"dphiMETjet/D");
   mytree->Branch("dphiMETll",&dphiMETll,"dphiMETll/D");
   mytree->Branch("etadilep",&etadilep,"etadilep/D");

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

        vector<TLorentzVector> leptons;
        vector<TLorentzVector> jets;
        vector<TLorentzVector> bjets;  
        vector<TLorentzVector> taus;
        vector<TLorentzVector> photons;
        vector<int> taucharges;

        TLorentzVector temp(0,0,0,0);
        TLorentzVector dilep(0,0,0,0);
        TLorentzVector invisible(0,0,0,0);
        TLorentzVector jetall(0,0,0,0);
        if(jentry%1000==0) cout<< Electron_size << " " << Muon_size << " " << Photon_size << " " << Jet_size << "@" << jentry << "\n";

        for (int i = 0; i < Electron_size; i++){
            temp.SetPtEtaPhiM(Electron_PT[i],Electron_Eta[i], Electron_Phi[i], 0.0);
            if (temp.Pt() > 20.0 && fabs(temp.Eta() ) < 2.5){
            leptons.push_back(temp);
            }
        }

        for (int i = 0; i < Muon_size; i++){
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
            taucharges.push_back(Jet_Charge[i]);
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

        sort_by_pt(leptons);


        if ( taus.size()>=1 && leptons.size() >= 1){

            if(leptons.at(0).Pt()< 30) continue;
            if(TMath::Abs(leptons.at(0).Eta())>2.4) continue; //missed small eta cut between 1.37 and 1.52
            if(taus.at(0).Pt()<25) continue;
            if(Electron_size > 1 && Electron_PT[1]>15) continue;
            if(Muon_size >1 && Muon_PT[1]>7 && TMath::Abs(Muon_Eta[1])<2.5) continue;
            int leptoncharge=0;
            if(Muon_size==0) leptoncharge=Electron_Charge[0];
            else if(Electron_size==0) leptoncharge=Muon_Charge[0];
            else if(Electron_PT[0]>Muon_PT[0]) leptoncharge=Electron_Charge[0];
            else if(Electron_PT[0]<Muon_PT[0]) leptoncharge=Muon_Charge[0];

            if(leptoncharge + taucharges.at(0)!=0) continue;
            if(TMath::Abs(taus.at(0).DeltaPhi(leptons.at(0)))<2.4) continue;
            if(get_MT(leptons.at(0), invisible)>40) continue;
            if((leptons.at(0)+taus.at(0)).M()<110&&(leptons.at(0)+taus.at(0)).M()>80) continue;
            
            if(bjets.size()>0) n_bbltau++;
            else n_ltau++;

            // for (unsigned int i = 0; i < jets.size(); i++){
            // JETPT[i] = jets[i].Pt();
            // jetall = jetall + jets[i];
            // }

            // for (unsigned int i = 0; i < bjets.size(); i++){
            // BJETPT[i] = bjets[i].Pt();
            // }
            
            // n_jets=Jet_size;
            // n_b=bjets.size();
            // n_a=0;
            // n_leptons=Electron_size + Muon_size;
            
            // LEPPT = leptons.at(0).Pt();
            // DRTLEP = leptons[0].DeltaR(leptons[1]);
            
            // MTLEP = get_MT(leptons[0], invisible);

            // dilep = leptons[0] + leptons[1];
            
            // DLPT = dilep.Pt();

            // etadilep = dilep.Rapidity();

            // MTALL = get_MT(invisible, dilep);
            
            // MLTAU = dilep.M();
            // MET = invisible.Pt();
            
            // MCT = TMath::Sqrt((TMath::Sqrt(dilep.Pt()*dilep.Pt()+dilep.M()*dilep.M())+invisible.Pt())*(TMath::Sqrt(dilep.Pt()*dilep.Pt()+dilep.M()*dilep.M())+invisible.Pt())-(dilep+invisible).Pt()*(dilep+invisible).Pt());
            // // Following the definitions in ATLAS-CONF-2016-056.pdf
            // MTZZ = TMath::Sqrt((TMath::Sqrt(dilep.Pt()*dilep.Pt()+dilep.M()*dilep.M())+TMath::Sqrt(MET*MET+91.187*91.187))*(TMath::Sqrt(dilep.Pt()*dilep.Pt()+dilep.M()*dilep.M())+TMath::Sqrt(MET*MET+91.187*91.187))-(dilep+invisible).Pt()*(dilep+invisible).Pt());
            // MT = TMath::Sqrt(2*dilep.Pt()*MET*(1-TMath::Cos(dilep.DeltaPhi(invisible))));
            // PTmissjet = (invisible + jetall).Pt();
            // dphiMETjet = invisible.DeltaPhi(jetall);
            // dphiMETll = invisible.DeltaPhi(dilep);
        }
        mytree->Fill();
   }
    mytree->Write();
    output->Close();
}


