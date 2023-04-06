#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <strstream>
#include "newstyle.C"

void comparer_200(){

  newstyle(0);
  
  TFile *fb = new TFile("higgsinos_ttbar.root");
  TFile *fs = new TFile("higgsinos_400_200.root");
  TFile *fs1 = new TFile("higgsinos_500_200.root");
  TFile *fs2 = new TFile("higgsinos_300_200.root");
 
  TTree *tb = (TTree*)fb->Get("t1");
  TTree *ts = (TTree*)fs->Get("t1");
  TTree *ts1 = (TTree*)fs1->Get("t1");
  TTree *ts2 = (TTree*)fs2->Get("t1");

  char* standcuts = "MT22 > 0";
  //  char* standcuts = "MT22 > 0 && BJMULT > 1";
  //  char* standcuts = "MT22 > 80";
  //char* standcuts = "MT22 > 0 && BJMULT > 1 && MET > 50.";
  //  char* standcuts = "MT22 > 0 && BJMULT > 1";
  //  char* standcuts = "BJETPT[0]> 100. && MET > 80.";

  vector<TH1F*> bhists;
  vector<TH1F*> shists;
  vector<TH1F*> shists1;
  vector<TH1F*> shists2;

  TH1F *b_met = new TH1F("b_met", ";MET (GeV);area norm", 30, 0, 300);
  TH1F *b_mt2 = new TH1F("b_mt2", ";mT2 (GeV);area norm", 20, 0, 200);
  TH1F *b_tapt = new TH1F("b_tapt", ";pT,#tau (GeV);area norm", 30, 0,300);
  TH1F *b_leppt = new TH1F("b_leppt", ";pT,lep (GeV);area norm", 30, 0,300);
  TH1F *b_nb = new TH1F("b_nb", ";# b jets;area norm", 5, 0, 5);
  TH1F *b_nj = new TH1F("b_nj", ";# light jets;area norm", 8, 0, 8);
  TH1F *b_bpt = new TH1F("b_bpt", ";pT,lead b (GeV);area norm", 30, 0,300);
  
  TH1F *s_met = new TH1F("s_met", ";MET (GeV);area norm", 30, 0, 300);
  TH1F *s_mt2 = new TH1F("s_mt2", ";mT2 (GeV);area norm", 20, 0, 200);
  TH1F *s_tapt = new TH1F("s_tapt", ";pT,#tau (GeV);area norm", 30, 0,300);
  TH1F *s_leppt = new TH1F("s_leppt", ";pT,lep (GeV);area norm", 30, 0,300);
  TH1F *s_nb = new TH1F("s_nb", ";# b jets;area norm", 5, 0, 5);
  TH1F *s_nj = new TH1F("s_nj", ";# light jets;area norm", 8, 0, 8);
  TH1F *s_bpt = new TH1F("s_bpt", ";pT,lead b (GeV);area norm", 30, 0,300);

  TH1F *s1_met = new TH1F("s1_met", ";MET (GeV);area norm", 30, 0, 300);
  TH1F *s1_mt2 = new TH1F("s1_mt2", ";mT2 (GeV);area norm", 20, 0, 200);
  TH1F *s1_tapt = new TH1F("s1_tapt", ";pT,#tau (GeV);area norm", 30, 0,300);
  TH1F *s1_leppt = new TH1F("s1_leppt", ";pT,lep (GeV);area norm", 30, 0,300);
  TH1F *s1_nb = new TH1F("s1_nb", ";# b jets;area norm", 5, 0, 5);
  TH1F *s1_nj = new TH1F("s1_nj", ";# light jets;area norm", 8, 0, 8);
  TH1F *s1_bpt = new TH1F("s1_bpt", ";pT,lead b (GeV);area norm", 30, 0,300);

  TH1F *s2_met = new TH1F("s2_met", ";MET (GeV);area norm", 30, 0, 300);
  TH1F *s2_mt2 = new TH1F("s2_mt2", ";mT2 (GeV);area norm", 20, 0, 200);
  TH1F *s2_tapt = new TH1F("s2_tapt", ";pT,#tau (GeV);area norm", 30, 0,300);
  TH1F *s2_leppt = new TH1F("s2_leppt", ";pT,lep (GeV);area norm", 30, 0,300);
  TH1F *s2_nb = new TH1F("s2_nb", ";# b jets;area norm", 5, 0, 5);
  TH1F *s2_nj = new TH1F("s2_nj", ";# light jets;area norm", 8, 0, 8);
  TH1F *s2_bpt = new TH1F("s2_bpt", ";pT,lead b (GeV);area norm", 30, 0,300);
  
  tb->Draw("MET >> b_met", standcuts); bhists.push_back(b_met);
  tb->Draw("MT22 >> b_mt2", standcuts); bhists.push_back(b_mt2);
  tb->Draw("TAUPT >> b_tapt", standcuts); bhists.push_back(b_tapt);
  tb->Draw("LEPPT >> b_leppt", standcuts); bhists.push_back(b_leppt);
  tb->Draw("BJMULT >> b_nb", standcuts); bhists.push_back(b_nb);
  tb->Draw("JMULT >> b_nj", standcuts); bhists.push_back(b_nj);
  tb->Draw("BJETPT[0] >> b_bpt", standcuts); bhists.push_back(b_bpt);
  
  ts->Draw("MET >> s_met", standcuts); shists.push_back(s_met);
  ts->Draw("MT22 >> s_mt2", standcuts); shists.push_back(s_mt2);
  ts->Draw("TAUPT >> s_tapt", standcuts); shists.push_back(s_tapt);
  ts->Draw("LEPPT >> s_leppt", standcuts); shists.push_back(s_leppt);
  ts->Draw("BJMULT >> s_nb", standcuts); shists.push_back(s_nb);
  ts->Draw("JMULT >> s_nj", standcuts); shists.push_back(s_nj);
  ts->Draw("BJETPT[0] >> s_bpt", standcuts); shists.push_back(s_bpt);
  
  ts1->Draw("MET >> s1_met", standcuts); shists1.push_back(s1_met);
  ts1->Draw("MT22 >> s1_mt2", standcuts); shists1.push_back(s1_mt2);
  ts1->Draw("TAUPT >> s1_tapt", standcuts); shists1.push_back(s1_tapt);
  ts1->Draw("LEPPT >> s1_leppt", standcuts); shists1.push_back(s1_leppt);
  ts1->Draw("BJMULT >> s1_nb", standcuts); shists1.push_back(s1_nb);
  ts1->Draw("JMULT >> s1_nj", standcuts); shists1.push_back(s1_nj);
  ts1->Draw("BJETPT[0] >> s1_bpt", standcuts); shists1.push_back(s1_bpt);

  ts2->Draw("MET >> s2_met", standcuts); shists2.push_back(s2_met);
  ts2->Draw("MT22 >> s2_mt2", standcuts); shists2.push_back(s2_mt2);
  ts2->Draw("TAUPT >> s2_tapt", standcuts); shists2.push_back(s2_tapt);
  ts2->Draw("LEPPT >> s2_leppt", standcuts); shists2.push_back(s2_leppt);
  ts2->Draw("BJMULT >> s2_nb", standcuts); shists2.push_back(s2_nb);
  ts2->Draw("JMULT >> s2_nj", standcuts); shists2.push_back(s2_nj);
  ts2->Draw("BJETPT[0] >> s2_bpt", standcuts); shists2.push_back(s2_bpt);
  

  for (int i = 0; i < shists.size(); i++){
    shists[i]->Scale(1.0/shists[i]->Integral() ); shists[i]->SetLineColor(kBlue);
    shists1[i]->Scale(1.0/shists1[i]->Integral() ); shists1[i]->SetLineColor(kBlue-9);
    shists2[i]->Scale(1.0/shists2[i]->Integral() ); shists2[i]->SetLineColor(kBlue+2);
    bhists[i]->Scale(1.0/bhists[i]->Integral() ); bhists[i]->SetLineColor(kRed);
  }

  Double_t theBR = 1000*151.70789*b_met->GetEntries()/2800000;
  Double_t theS = 1000*1.03385*s_met->GetEntries()/60000;
  Double_t theS1 =  1000*0.287268*s_met->GetEntries()/60000;
  Double_t theS2 =  1000*4.91704*s_met->GetEntries()/60000;

  Double_t S1 = theS1/sqrt(theBR);
  Double_t S0 = theS/sqrt(theBR);
  Double_t S2 = theS2/sqrt(theBR);

  cout << endl;
  cout << "sigmaBR (fb):\t" << theBR << endl;
  cout << "**************" << endl;
  cout << "sigma(mst300):\t" <<  theS1 << "\t" << theS1/theBR << "\t" << theS1/sqrt(theBR) << "\t" << 25/(S1*S1) << endl;
  cout << "sigma(mst400):\t" << theS << "\t" << theS/theBR << "\t" << theS/sqrt(theBR) << "\t" << 25/(S0*S0) <<endl;
  cout << "sigma(mst500):\t" << theS2 << "\t" << theS2/theBR << "\t" << theS2/sqrt(theBR) << "\t" << 25/(S2*S2) << endl;
  cout << endl;
  
  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->cd();
  c3->SetLogy();
  
  for (int j = 0; j < shists.size(); j++){
    bhists[j]->SetTitle("13 TeV, ttbar (red), mstop = 500, 400, 300 GeV, m_{#tilde{#tau}}= 200 GeV");
    bhists[j]->Draw();
    shists[j]->Draw("same");
    shists1[j]->Draw("same");
    shists2[j]->Draw("same");
  if (j == 0){
      c3->SaveAs("higgsino_shapes_200.pdf(", "pdf");
    } else if ( j == shists.size() -1){
      c3->SaveAs("higgsino_shapes_200.pdf)", "pdf");
    } else {
      c3->SaveAs("higgsino_shapes_200.pdf","pdf");
    }

    
  }
}
