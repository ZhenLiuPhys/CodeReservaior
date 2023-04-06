//Analysing events in root format
//For charged higgs pair production and decays via dibosons
//Created by Zhen Liu
//Dec. 2011 UW-Madison
//zliu57@wisc.edu

#define Sig_cxx
#include "Sig.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include "TLorentzVector.h"

const double Pi = 3.14159265359;


inline Int_t hardronicw(vector<TLorentzVector> jets){
  //inline function to see if at least one pair of jets/leptons reconstructed as w/z shows the behavior of decay products of vector bosons by requiring them to be close to their mother boson direction due to the $1+cos\theta$ amplitudes.
  Int_t njets = jets.size();
  Int_t flag = 0;
  TLorentzVector temp1, temp2, temp12;
  for(Int_t i = 0; i < njets-1  && flag != 1; i++ ){
    temp1 = jets.at(i);
    for(Int_t j = i+1; j < njets && flag != 1; j++){
      temp2 = jets.at(j);
      temp12 = temp1 + temp2;
      if(temp12.Mag()>=75 && temp12.Mag()<=85){
	temp12=temp1+temp2;
	TVector3 boost = temp12.BoostVector();
	temp1.Boost(-boost);
	TVector3 vect1 = temp1.Vect();
	Double_t angle = vect1.Angle(boost);
	if(4*angle<=TMath::Pi()||4/3*angle>=TMath::Pi()) flag = 1;
      }
    }
  }
  return flag;
}

inline Double_t deltaphi(Double_t phi1, Double_t phi2){
  Double_t result = TMath::Abs(phi1-phi2);
  if(result>TMath::Pi()) result=2*TMath::Pi()-result;
  return result;
}

inline Double_t deltaR(TLorentzVector v1, TLorentzVector v2){
  Double_t result;
  result=TMath::Sqrt(pow(deltaphi(v1.Phi(),v2.Phi()),2) + pow(v1.Eta()-v2.Eta(),2));
  return result;
}


inline Int_t separationcut(vector<TLorentzVector> leptons, Double_t r_cut){
  //Function returns where all the momentums of the vector<TLorentzVector> have speration smaller than the r_cut. Important for the spherical behavior check.
  Int_t flag = 1;
  for(Int_t i = 0 ; i< leptons.size() && flag!=0; i++){
    for(Int_t j = i + 1 ; j< leptons.size()-1 && flag!=0; j++){
      if(deltaR(leptons.at(i),leptons.at(j))> r_cut){
	flag = 0;
	break;
      }
    }
  }
  return flag;
}

inline Double_t separationllmax(vector<TLorentzVector> leptons){
  //  Int_t n = leptons.size()*(leptons.size()-1)/2-1;
  //  Double_t result=0;
  //  if(n<=1)return result;
  
  Int_t n=20;
  Double_t sep[20];
  for(Int_t i = 0 ; i<= n; i++)sep[i]=0;
  sep[0]=100;
  Int_t k=0;
  for(Int_t i = 0 ; i< leptons.size(); i++){
    for(Int_t j = i + 1 ; j< leptons.size()-1; j++){
      sep[k]=deltaR(leptons.at(i),leptons.at(j));
      k++;
    }
  }
  return TMath::MaxElement(n,sep);
}

inline Double_t separationljmin(vector<TLorentzVector> leptons, vector<TLorentzVector> jets){
  //  Int_t n = leptons.size()*(leptons.size()-1)/2-1;
  //  Double_t result=0;
  //  if(n<=1)return result;
  
  Int_t n=100;
  Double_t sep[100];
  for(Int_t i = 0 ; i<= n; i++)sep[i]=100;
  sep[0]=-1;
  Int_t k=0;
  for(Int_t i = 0 ; i< leptons.size(); i++){
    for(Int_t j = 0 ; j< jets.size(); j++){
      sep[k]=deltaR(leptons.at(i),jets.at(j));
      k++;
    }
  }
  return TMath::MinElement(n,sep);
}

double uniformRandom()
{
  return rand()*1./RAND_MAX;
}

double Normal()
{
  double U = uniformRandom();
  double V = uniformRandom();
  return sqrt(-2*log(U))*cos(2*Pi*V);  // Box-Muller method
}

double Smear(double E, double Ra, double Rb)// smearing, a/Sqrt(E)\oplus b
{
  double temp_a = Ra*TMath::Sqrt(E);
  double temp_b = Rb*E;
  return E + Normal()*TMath::Sqrt(temp_a*temp_a + temp_b*temp_b);
}


void Sig::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L nparticle.C
//      Root > nparticle t
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
   TFile *output = new TFile("Sigoutput.root","recreate");
   TTree *mytree = new TTree("mytree","mytree");
   TTree *atree = new TTree("atree","atree");
   
   Double_t beam_R=0.001/TMath::Sqrt(2.);
   Double_t a_Ra=0.17;
   Double_t a_Rb=0.01;
   
   Int_t n_jets;//number of jets
   Int_t n_leptons;//number of leptons
   Int_t n_vl;//number of neutrinos
   Int_t n_b;//number of tagged b quark in the pythia event file
   Int_t n_a;//number of photons
   Double_t m_recoil;//invirant mass of the different combinations of pair of jets
   Double_t E_a;
   Double_t PT_a;
   Double_t Pz_a;
   Double_t E_a_min;
   Double_t PT_a_min;

   mytree->Branch("n_jets",&n_jets,"n_jets/I");
   mytree->Branch("n_leptons",&n_leptons,"n_leptons/I");
   mytree->Branch("n_vl",&n_vl,"n_vl/I");
   mytree->Branch("n_b",&n_b,"n_b/I");
   mytree->Branch("n_a",&n_a,"n_a/I");

   mytree->Branch("E_a_min",&E_a_min,"E_a_min/D");
   mytree->Branch("PT_a_min",&PT_a_min,"PT_a_min/D");
   
   atree->Branch("m_recoil",&m_recoil,"m_recoil/D");
   atree->Branch("E_a",&E_a,"E_a/D");
   atree->Branch("PT_a",&PT_a,"PT_a/D");
   atree->Branch("Pz_a",&Pz_a,"Pz_a/D");


   n_events = 0;
   
   srand(11232698);


for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
    
      n_jets = 0;
      n_leptons = 0;
      n_vl = 0;
      n_b = 0;
      n_a = 0;
      
      E_a_min;
      PT_a_min;


      TLorentzVector temp_all;
      TLorentzVector jet_all;
      TLorentzVector b_all; //Already know only two b's, be careful when generalize
      TLorentzVector lepton_all;
      TLorentzVector photon_all;
      TLorentzVector vl;
      vector<TLorentzVector> jet;
      vector<TLorentzVector> lepton;
      vector<TLorentzVector> bjet;
      vector<TLorentzVector> photon;
      vector<Int_t> leptypes;


      for(Int_t i = 0; i<Particle_;i++){
	if(Particle_Status[i]==1){
	  if(Particle_PID[i]==21||TMath::Abs(Particle_PID[i])<=4){
	     n_jets++;
	     TLorentzVector jet_temp;
	     jet_temp.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	     jet.push_back(jet_temp);
	     jet_all.SetPxPyPzE(jet_all.Px()+Particle_Px[i],jet_all.Py()+Particle_Py[i],jet_all.Pz()+Particle_Pz[i],jet_all.E()+Particle_E[i]);
	   }
	  elseif(TMath::Abs(Particle_PID[i])==5){
	    n_b++;
	    TLorentzVector jet_temp;
	    jet_temp.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	    bjet.push_back(jet_temp);
	    b_all.SetPxPyPzE(b_all.Px()+Particle_Px[i],b_all.Py()+Particle_Py[i],b_all.Pz()+Particle_Pz[i],b_all.E()+Particle_E[i]);
	  }
	  else if(TMath::Abs(Particle_PID[i])==11||TMath::Abs(Particle_PID[i])==13||TMath::Abs(Particle_PID[i])==15){
	     n_leptons++;
	     leptypes.push_back(Particle_PID[i]);
	     TLorentzVector lepton_temp;
	     lepton_temp.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	     lepton.push_back(lepton_temp);
	     lepton_all.SetPxPyPzE(lepton_all.Px()+Particle_Px[i],lepton_all.Py()+Particle_Py[i],lepton_all.Pz()+Particle_Pz[i],lepton_all.E()+Particle_E[i]);
	   }
	  else if(TMath::Abs(Particle_PID[i])==12||TMath::Abs(Particle_PID[i])==14||TMath::Abs(Particle_PID[i])==16){
	     n_vl++;
	     vl.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	   }
	  else if(TMath::Abs(Particle_PID[i])==22){
	     TLorentzVector photon_temp;
	     n_a++;
	     photon_temp.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	     photon.push_back(photon_temp);
	   }
	}
      }


      for(int i = 0; i < photon.size(); i++){
        PT_a = Smear(photon.at(i).Pt(), a_Ra, a_Rb);
        //Pz_a = Smear(photon.at(i).Pz(), a_Ra, a_Rb);
        E_a = Smear(photon.at(i).E(), a_Ra, a_Rb);
        double shat = Smear(3000., 0., beam_R);
        m_recoil = TMath::Sqrt((shat - E_a)*(shat - E_a) - E_a*E_a);
        if(PT_a < PT_a_min || i == 0) PT_a_min = PT_a;
        if(E_a < E_a_min || i == 0) E_a_min = E_a;
        atree->Fill();
      }

      mytree->Fill();
      // if (Cut(ientry) < 0) continue;
   }
 atree->Write();
 mytree->Write();
 output->Close();
}
