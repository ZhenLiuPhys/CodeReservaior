#include <iostream>
#include <fstream>

const double cs[24]={15.7131,17.3917,18.6378,20.424,23.015,26.952,32.903,42.607,59.789,74.121,95.944,123.873,143.849,170.566,208.46,233.34,265.31,305.88,360.52,435.62,551,740.56,948.3,1222.4};
const double csbkg=2.6*1000*1000;

const double mass[24]={500,1000,1200,1400,1600,1800,2000,2200,2400,2500,2600,2680,2720,2760,2800,2820,2840,2860,2880,2900,2920,2940,2960,2980};

double signif(int countsig, int countbkg, int i){
if(countbkg!=0){
return (double)countsig*0.8*0.96/10000.*6.4*cs[i]/TMath::Sqrt((double)countbkg/100000.*6.4*csbkg);}
else{
return 0.1;}
}


void analysis2v3(){

ofstream myfile;
myfile.open ("significances2v3.dat");

TFile *f2 = TFile::Open("BKG2v3.root");
TTree *tree2 = (TTree*) f2->Get("atree");

int upper, upper_final;
int lower, lower_final;
double center;
double significance;

int countsig;
int countbkg;

for(int j = 0; j < 72; j ++){

	i=j/3;

	stringstream input;
	stringstream in_filename;

	in_filename << "Sigoutput" << j << ".root";

	TFile *f1 = TFile::Open(&in_filename.str()[0]);
	TTree *tree1 = (TTree*) f1->Get("atree");
	in_filename.str("");

	center = (3000.*3000.-mass[i]*mass[i])/6000.;
	lower = (int) center;
	upper = (int) center+1;
	upper_final=upper;
	lower_final=lower;
	int fails = 0;

	input << "E_a<" << upper << "&&E_a>" << lower << "&&PT_a>10";
	countsig = tree1->GetEntries(&input.str()[0]);
	countbkg = tree2->GetEntries(&input.str()[0]);
	printf("%s, sig=%d, bkg=%d\n",&input.str()[0],countsig, countbkg);
	significance=signif(countsig,countbkg,i);
	lower--;
	input.str("");

	while(fails<=6&&lower>10){
	input << "E_a<" << upper << "&&E_a>" << lower << "&&PT_a>10";
	countsig = tree1->GetEntries(&input.str()[0]);
	countbkg = tree2->GetEntries(&input.str()[0]);
	printf("%s, sig=%d, bkg=%d, sig=%f\n",&input.str()[0],countsig, countbkg, significance);
	if(signif(countsig,countbkg,i) < significance) {
		fails++;
		if(fails%2==0&&lower>11) {lower--;}
		else {upper++;}
	}
	else{
		significance=signif(countsig,countbkg,i);
		upper_final=upper;
		lower_final=lower;	
		if(fails%2==0&&lower>11) {lower--;}
		else {upper++;}
	}
	input.str("");
	}

	myfile << mass[i] << "\t" << (upper_final+lower_final)/2 << "\t" << lower_final << "\t" << upper_final << "\t" << significance <<endl ;

}
myfile.close();
}
