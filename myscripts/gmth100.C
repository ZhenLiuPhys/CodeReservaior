//Analysing in root format
//For charged higgs pair production and decays via dibosons
//Created by Zhen Liu
//Dec. 2011 UW-Madison
//zliu57@wisc.edu
//Modified to execute the cuts descirbed on Table 1 in the paper 0911.3656
//WW longtitugnal mode

#define gmth100_cxx
#include "gmth100.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include "TLorentzVector.h"




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
	if(4*angle<=TMath::Pi()||4./3*angle>=TMath::Pi()) flag = 1;
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
  result=TMath::Sqrt(TMath::Power(deltaphi(v1.Phi(),v2.Phi()),2) + TMath::Power(v1.Eta()-v2.Eta(),2));
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



void gmth100::Loop()
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
   TFile *output = new TFile("output.root","recreate");
   TTree *mytree = new TTree("mytree","mytree");
   Int_t n_jets;//number of jets
   Int_t n_leptons;//number of leptons
   Int_t n_vl;//number of neutrinos
   Int_t n_b;//number of tagged b quark in the pythia event file
   Double_t inv_12;//invirant mass of the different combinations of pair of jets
   Double_t inv_34;
   Double_t inv_13;
   Double_t inv_24;
   Double_t inv_14;
   Double_t inv_23;
   Double_t m_eff;//invariant mass of all the visible particles
   Double_t m_eff_jets;//invariant mass of all the jets
   Double_t m_eff_leptons;//invariant mass of all the leptons
   Double_t h_t;//h_t of all visible particles
   Double_t h_t_jets;//h_t of all the jets
   Double_t h_t_leptons;//h_t of all the leptons
   Double_t p_t;//p_t of all visible particles
   Double_t p_t_jets;//p_t of all the jets
   Double_t p_t_leptons;//p_t of all the leptons
   Double_t vl_pt;//p_t of the neutrino in the event
   Int_t n_fake;//different combinations that fall into the range of w/z mass. 2-12/34, 4-13/24, 6-14/23, 11-doesn't fall into range but still with 3 leptons and 4 jets in final state.
   Int_t n_pass;//number of events with particles other than jets/leptons/neutrinos.Should be 0.
   Int_t n_w;//number of events with at least one reconstructed hadronic w/z with vector behavior described in hadronicw function.
   Double_t m_tr;//negative of the transverse mass, where the missing p_t is the opposite of that of all visible particles in the final state.
   Double_t m_trr;//transverse mass, where the missing p_t is the neutrino particle in the events file.
   Double_t m_tr_clustered;
   Double_t m_minimize;
   Double_t n_sep;//This tags whether all the leptons have separations smaller than pi/2. If it equals 1, basically all the leptons comes from a same mother charged higgs.
/*   Int_t n_passll;//Added for the study of Tao Han, David Krohn, Lian-Tao Wang, Wenhan Zhu's paper
   Int_t n_passjj;
   Int_t n_passll1;
   Int_t n_passll2;
   Int_t n_passll3;
   Int_t n_passll4;
   Int_t n_passjj1;
   Int_t n_passjj2;
*/
//   mytree->Branch("n_jets",&n_jets,"n_jets/I");
//   mytree->Branch("n_leptons",&n_leptons,"n_leptons/I");
//   mytree->Branch("n_vl",&n_vl,"n_vl/I");
//   mytree->Branch("n_b",&n_b,"n_b/I");
//   mytree->Branch("n_fake",&n_fake,"n_fake/I");
   mytree->Branch("m_eff",&m_eff,"m_eff/D");
   mytree->Branch("m_eff_jets",&m_eff_jets,"m_eff_jets/D");
   mytree->Branch("m_eff_leptons",&m_eff_leptons,"m_eff_leptons/D");
   mytree->Branch("h_t",&h_t,"h_t/D");
   mytree->Branch("h_t_jets",&h_t_jets,"h_t_jets/D");
   mytree->Branch("h_t_leptons",&h_t_leptons,"h_t_leptons/D");
   mytree->Branch("p_t",&p_t,"p_t/D");
   mytree->Branch("p_t_jets",&p_t_jets,"p_t_jets/D");
   mytree->Branch("p_t_leptons",&p_t_leptons,"p_t_leptons/D");

/*   mytree->Branch("inv_12",&inv_12,"inv_12/D");
   mytree->Branch("inv_34",&inv_34,"inv_34/D");
   mytree->Branch("inv_13",&inv_13,"inv_13/D");
   mytree->Branch("inv_24",&inv_24,"inv_24/D");
   mytree->Branch("inv_14",&inv_14,"inv_14/D");
   mytree->Branch("inv_23",&inv_23,"inv_23/D");
*/

//   mytree->Branch("n_pass",&n_pass,"n_pass/I");
//   mytree->Branch("n_w",&n_w,"n_w/I");

//   mytree->Branch("vl_pt",&vl_pt,"vl_pt/D");
   mytree->Branch("m_tr",&m_tr,"m_tr/D");
   mytree->Branch("m_trr",&m_trr,"m_trr/D");
   mytree->Branch("m_tr_clustered",&m_tr_clustered,"m_tr_clustered/D");
   mytree->Branch("m_minimize",&m_minimize,"m_minimize/D");
   mytree->Branch("n_sep",&n_sep,"n_sep/D");

/*
   mytree->Branch("n_passll",&n_passll, "n_passll/I");
   mytree->Branch("n_passjj",&n_passjj, "n_passjj/I");
   mytree->Branch("n_passll1",&n_passll1, "n_passll1/I");
   mytree->Branch("n_passll2",&n_passll2, "n_passll2/I");
   mytree->Branch("n_passll3",&n_passll3, "n_passll3/I");
   mytree->Branch("n_passll4",&n_passll4, "n_passll4/I");
   mytree->Branch("n_passjj1",&n_passjj1, "n_passjj1/I");
   mytree->Branch("n_passjj2",&n_passjj2, "n_passjj2/I");
*/

for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
    
      n_jets = 0;
      n_leptons = 0;
      n_vl = 0;
      n_b = 0;
      n_fake = 0;
/*
      inv_12 = 0;
      inv_34 = 0;
      inv_13 = 0;
      inv_24 = 0;
      inv_14 = 0;
      inv_23 = 0;
*/

      h_t = 0;
      h_t_jets = 0;
      h_t_leptons = 0;
      m_eff = 0;
      m_eff_jets = 0;
      m_eff_leptons = 0;
      p_t = 0;
      p_t_jets = 0;
      p_t_leptons = 0;
      n_pass=0;
      vl_pt = 0;
      m_tr = 0;
      m_trr = 0;
      m_minimize = 0;
      m_tr_clustered = 0;
      n_sep = 0;
/*
      n_passll = 0;
      n_passjj = 0;
      n_passll1 = 0;
      n_passll2 = 0;
      n_passll3 = 0;
      n_passll4 = 0;
      n_passjj1 = 0;
      n_passjj2 = 0;
*/

      TLorentzVector temp_all;
      TLorentzVector jet_all;
      TLorentzVector lepton_all;
      TLorentzVector vl;
      vector<TLorentzVector> jet;
      vector<TLorentzVector> lepton;
      vector<Int_t> leptypes;


      for(Int_t i = 0; i<Particle_;i++){
	if(Particle_Status[i]==1){
	  if(Particle_PID[i]==21||TMath::Abs(Particle_PID[i])<=5){
	     if(TMath::Abs(Particle_PID[i])==5) n_b++;
	     n_jets++;
	     h_t_jets = h_t_jets + Particle_PT[i];
	     TLorentzVector jet_temp;
	     jet_temp.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	     jet.push_back(jet_temp);
	     jet_all.SetPxPyPzE(jet_all.Px()+Particle_Px[i],jet_all.Py()+Particle_Py[i],jet_all.Pz()+Particle_Pz[i],jet_all.E()+Particle_E[i]);
	   }
	   else if(TMath::Abs(Particle_PID[i])==11||TMath::Abs(Particle_PID[i])==13||TMath::Abs(Particle_PID[i])==15){
	     n_leptons++;
	     leptypes.push_back(Particle_PID[i]);
	     h_t_leptons = h_t_leptons + Particle_PT[i];
	     TLorentzVector lepton_temp;
	     lepton_temp.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	     lepton.push_back(lepton_temp);
	     lepton_all.SetPxPyPzE(lepton_all.Px()+Particle_Px[i],lepton_all.Py()+Particle_Py[i],lepton_all.Pz()+Particle_Pz[i],lepton_all.E()+Particle_E[i]);
	   }
	   else if(TMath::Abs(Particle_PID[i])==12||TMath::Abs(Particle_PID[i])==14||TMath::Abs(Particle_PID[i])==16){
	     n_vl++;
	     vl_pt=Particle_PT[i];
	     vl.SetPxPyPzE(Particle_Px[i],Particle_Py[i],Particle_Pz[i],Particle_E[i]);
	   }
	   else n_pass++;
	}
      }

      h_t = h_t_leptons+h_t_jets;

      temp_all.SetPxPyPzE(jet_all.Px()+lepton_all.Px(),jet_all.Py()+lepton_all.Py(),jet_all.Pz()+lepton_all.Pz(),jet_all.E()+lepton_all.E());

      m_eff = temp_all.Mag();
      p_t = temp_all.Pt();
      m_eff_jets = jet_all.Mag();
      p_t_jets = jet_all.Pt();
      m_eff_leptons = lepton_all.Mag();
      p_t_leptons = lepton_all.Pt();

      TLorentzVector met, lept;

      met.SetPx(-temp_all.Px());
      met.SetPy(-temp_all.Py());
      met.SetE(met.Pt());

      lept.SetPx(lepton_all.Px());
      lept.SetPy(lepton_all.Py());
      lept.SetE(lept.Pt());

      m_tr=(met+lept).Mag();
      m_trr=(lepton_all+vl).Mag();
      lept.SetE(TMath::Sqrt(TMath::Power(lept.Pt(),2) + lepton_all.Mag2()));
      m_tr_clustered=(lept+met).Mag();

      if(hardronicw(jet)==1) n_w = 1;
      n_sep = separationcut(lepton,TMath::Pi()/2);

/*
      if(n_leptons == 2){
	Int_t n = lepton.size();
	TLorentzVector temp1, temp2;
	temp1 = lepton.at(n-2);
	temp2 = lepton.at(n-1);
	if(temp1.Pt() > 100 && temp2.Pt() > 100) n_passll = 1;
	else n_passll = 0;
	if((temp1-temp2).Pt()>440) n_passll1 = 1;
	else n_passll1 = 0;
	if(cos(deltaphi(temp1.Phi(),temp2.Phi()))<-0.8) n_passll2 = 1;
	else n_passll2 = 0;
	if(TMath::Abs(temp1.Eta())<2 && TMath::Abs(temp2.Eta())<2) n_passll3 = 1;
	else n_passll3 = 0;
	if((temp1+temp2).Mag()>250) n_passll4 = 1;
	else n_passll4 = 0;
      }
      else n_passll = 0;

      if(n_jets == 2){
	Int_t n = jet.size();
	TLorentzVector temp1, temp2, temptag, tempveto;
	temp1 = jet.at(n-2);
	temp2 = jet.at(n-1);
	if(temp1.E()>800 && TMath::Abs(temp1.Eta())>3 && TMath::Abs(temp1.Eta()) <5 && temp1.Pt()>40){
	  n_passjj1 = 1;
	  temptag = temp1;
	  tempveto = temp2;
	}
	else if(temp2.E()>800 && TMath::Abs(temp2.Eta())>3 && TMath::Abs(temp2.Eta()) <5 && temp2.Pt()>40){
	    n_passjj1 = 1;
	    temptag = temp2;
	    tempveto = temp1;
	  }
	  else n_passjj1 = 0;

     	if((tempveto.Pt() < 30 || TMath::Abs(tempveto.Eta()) > 3)) n_passjj2 = 1;
          else n_passjj2 = 0;
        }
	else n_passjj2 = 0;
*/
      
      /*        TLorentzVector w1a, w1b , z1, w2, z2, wlep;

	if(n_jets==4&&n_leptons==3){
	  TLorentzVector temp, temp1, temp2, temp3;
	  TLorentzVector temp_12, temp_34, temp_13, temp_24, temp_14, temp_23;
	  Int_t n = jet.size();
	  Int_t lepflag = 0;
	  temp_12 = jet.at(n-1) + jet.at(n-2);
	  inv_12 = temp_12.Mag();
	  temp_34 = jet.at(n-3) + jet.at(n-4);	  
	  inv_34 = temp_34.Mag();
	  temp_13 = jet.at(n-1) + jet.at(n-3);
	  inv_13 = temp_13.Mag();
	  temp_24 = jet.at(n-2) + jet.at(n-4);
	  inv_24 = temp_24.Mag();
	  temp_14 = jet.at(n-1) + jet.at(n-4);
	  inv_14 = temp_14.Mag();
	  temp_23 = jet.at(n-2) + jet.at(n-3);
	  inv_23 = temp_23.Mag();

	    if(inv_12<=85&&inv_12>=75&&inv_34<=95&&inv_34>=85) {
	      n_fake = 2;
	      w2 = temp_12;
	      z2 = temp_34;
	    }
	    else if(inv_34<=85&&inv_34>=75&&inv_12<=95&&inv_12>=85) {
	      n_fake = 3;
	      w2 = temp_34;
	      z2 = temp_12;
	    }
	    else if(inv_13<=85&&inv_13>=75&&inv_24<=95&&inv_24>=85) {
	      n_fake = 4;
	      w2 = temp_13;
	      z2 = temp_24;
	    }
	    else if(inv_24<=85&&inv_24>=75&&inv_13<=95&&inv_13>=85) {
	      n_fake = 5;
	      w2 = temp_24;
	      z2 = temp_13;
	    }
	    else if(inv_14<=85&&inv_14>=75&&inv_23<=95&&inv_23>=85) {
	      n_fake = 6;
	      w2 = temp_14;
	      z2 = temp_13;
	    }
	    else if(inv_23<=85&&inv_23>=75&&inv_14<=95&&inv_14>=85) {
	      n_fake = 7;
	      w2 = temp_23;
	      z2 = temp_14;
	    }
	    else n_fake=11;

	    Int_t nlep = lepton.size();

	    if(TMath::Abs(leptypes.at(nlep-3))==TMath::Abs(leptypes.at(nlep-2))&&TMath::Abs(leptypes.at(nlep-3))==TMath::Abs(leptypes.at(nlep-1))){
	      temp1 = lepton.at(nlep-1) + lepton.at(nlep-2);
	      temp2 = lepton.at(nlep-1) + lepton.at(nlep-3);
	      temp3 = lepton.at(nlep-2) + lepton.at(nlep-3);
	      if(temp1.Mag()<=95&&temp1.Mag()>=85){
		z1 = temp1;
		wlep = lepton.at(nlep-3);
	      }
	      else if(temp2.Mag()<=95&&temp2.Mag()>=85){
		z1 = temp2;
		wlep = lepton.at(nlep-2);
	      }
	      else if(temp3.Mag()<=95&&temp3.Mag()>=85){
		z1 = temp3;
		wlep = lepton.at(nlep-1);
	      }
	      else lepflag = 1;
	    }
	    else if(TMath::Abs(leptypes.at(nlep-3))==TMath::Abs(leptypes.at(nlep-2))){
	      z1 = temp3;
	      wlep = lepton.at(nlep-1);
	    }
	    else if(TMath::Abs(leptypes.at(nlep-3))==TMath::Abs(leptypes.at(nlep-1))){
	      z1 = temp2;
	      wlep = lepton.at(nlep-2);
	    }
	    else if(TMath::Abs(leptypes.at(nlep-2))==TMath::Abs(leptypes.at(nlep-1))){
	      z1 = temp1;
	      wlep = lepton.at(nlep-3);
	    }
	    else lepflag = 1;
	}

	//start my decision tree
	if(n_fake<=7&&n_fake>=2&&lepflag==0){
	  //reconstruct the neutrino momentum; there are two solutions.
	  TLorentzVector vl1, vl2;
	  vl1.SetPx(-temp_all.Px());
	  vl1.SetPy(-temp_all.Py());
	  vl1.SetPz(sqrt(pow(80.398,4) - pow(wlep.Pt()*vl1.Pt(),2) + pow(vl1.Dot(wlep),2) - pow(wlep.Pz()*vl1.Pt(),2)));
	  vl1.SetE(sqrt(pow(vl1.Px(),2) + pow(vl1.Py(),2) + pow(vl1.Pz(),2)));
	  vl2 = vl1;
	  vl2.SetPz(-vl1.E());

	  w1a = vl1 + wlep;
	  w1b = vl2 + wlep;

	  vector<Double_t> masses1;
	  vector<Double_t> masses2;
	  masses1.push_back((w1a+z1).Mag());
	  masses2.push_back((w2+z2).Mag());
	  masses1.push_back((w1a+z2).Mag());
	  masses2.push_back((w2+z1).Mag());
	  masses1.push_back((w1b+z1).Mag());
	  masses2.push_back((w2+z2).Mag());
	  masses1.push_back((w1b+z2).Mag());
	  masses2.push_back((w2+z1).Mag());

	  for(Int_t i = 0; i < 3; i++) {
	    if(masses1.at(i)-masses2.at(i) < masses1.at(i+1)-masses2.at(i+1))
	      m_minimize=(masses1.at(i)+masses2.at(i))/2;
	    else
	      m_minimize=(masses1.at(i+1)+masses2.at(i+1))/2;
	  }
	}
      */

      mytree->Fill();
      // if (Cut(ientry) < 0) continue;
   }
 mytree->Write();
 output->Close();
}
