 void plot(){
  THStack hs("ZPrime","");

 TFile *_file0 = TFile::Open("ZPChi.root");
 TH1F *mll_bkgd = new TH1F("mll_bkgd","mll_bkgd",25,-1.,1.);
 mytree->Draw("cosstar_lzp>>mll_bkgd");
 mll_bkgd->SetLineColor(kBlack);
 mll_bkgd->Scale(1./mytree->GetEntries()*25./2.);
 mll_bkgd->GetYaxis()->SetRange(0.,1.05);
 hs.Add(mll_bkgd);

 TFile *_file1 = TFile::Open("ZPEta.root");
 TH1F *mll_signal = new TH1F("mll_signal","mll_signal",25,-1.,1.);
 mytree->Draw("cosstar_lzp>>mll_signal");
 mll_signal->SetLineColor(kRed);
 mll_signal->Scale(1./mytree->GetEntries()*25./2.);
 mll_signal->GetYaxis()->SetRange(0.,1.05);
 hs.Add(mll_signal); 

 TFile *_file2 = TFile::Open("ZPPsi.root");
 TH1F *mll_ZPPsi = new TH1F("mll_ZPPsi","mll_ZPPsi",25,-1.,1.);
 mytree->Draw("cosstar_lzp>>mll_ZPPsi");
 mll_ZPPsi->SetLineColor(kGreen);
 mll_ZPPsi->Scale(1./mytree->GetEntries()*25./2.);
 mll_ZPPsi->GetYaxis()->SetRange(0.,1.05);
 hs.Add(mll_ZPPsi); 

 TFile *_file3 = TFile::Open("ZPLR.root");
 TH1F *mll_ZPLR = new TH1F("mll_ZPLR","mll_ZPLR",25,-1.,1.);
 mytree->Draw("cosstar_lzp>>mll_ZPLR");
 mll_ZPLR->SetLineColor(kBlue);
 mll_ZPLR->Scale(1./mytree->GetEntries()*25./2.);
 mll_ZPLR->GetYaxis()->SetRange(0.,1.05);
 hs.Add(mll_ZPLR); 

 
 TLegend *leg = new TLegend(0.3,0.7,0.4,0.9);
 leg->AddEntry(mll_signal,"#eta");
 leg->AddEntry(mll_bkgd,"#Chi");
 leg->AddEntry(mll_ZPPsi,"#Psi");
 leg->AddEntry(mll_ZPLR,"LR");
 leg->SetTextSize(0.05);

 hs.SetMinimum(0.15); 
 TCanvas *myC = new TCanvas();
// myC->Divide(3,1);
// myC->cd(3);
 hs.Draw("nostack");
 hs.GetXaxis()->SetTitle("cos(#theta*)");
 hs.GetYaxis()->SetTitle("1/#sigma d#sigma/dcos(#theta*)"); 
 leg->Draw();
 hs.GetXaxis()->SetTitleSize(0.05);
 hs.GetXaxis()->SetLabelSize(0.05);
 hs.GetYaxis()->SetTitleSize(0.05);
 hs.GetYaxis()->SetLabelSize(0.05);
 myC->SaveAs("mll.eps");
 }

