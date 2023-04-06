//
//   Example of a program to fit non-equidistant data points
//   =======================================================
//
//   The fitting function fcn is a simple chisquare function
//   The data consists of 5 data points (arrays x,y,z) + the errors in errorsz
//   More details on the various functions or parameters for these functions
//   can be obtained in an interactive ROOT session with:
//    Root > TMinuit *minuit = new TMinuit(10);
//    Root > minuit->mnhelp("*")  to see the list of possible keywords
//    Root > minuit->mnhelp("SET") explains most parameters
//Author: Rene Brun

#include "TMinuit.h"

//Float_t zh,zhbb1,zhcc1,zhgg1,zhww1,zhtata1,zhzz1,zhgaga1;
Float_t erzh,erzhbb1,erzhcc1,erzhgg1,erzhww1,erzhtata1,erzhzz1,erzhgaga1;

//Float_t zhbb2,zhcc2,zhgg2,zhww2,zhtata2,zhzz2,zhgaga2;
Float_t erzh2, erzhbb2,erzhcc2,erzhgg2,erzhww2,erzhtata2,erzhzz2,erzhgaga2;

//Float_t vhbb,vhcc,vhgg,vhww,vhtata,vhzz,vhgaga;
Float_t ervhbb,ervhcc,ervhgg,ervhww,ervhtata,ervhzz,ervhgaga;
Float_t ervhbb1t,ervhcc1t,ervhgg1t,ervhww1t,ervhtata1t,ervhzz1t,ervhgaga1t;

//LHC
Float_t erPghzz, erPghww, erPghgaga, erPwwzz, erPgagazz;
Float_t erPCvbfbb, erPCghww, erPCvbftata, erPCghzz, erPCghgaga;
Float_t erPC2vbfbb, erPC2ghww, erPC2vbftata, erPC2ghzz, erPC2ghgaga;

Float_t brhbb, brhcc, brhgg, brhww, brhtata, brhzz, brhgaga, brhmm;

Float_t ervhbb1, erzhmm, erPghmm, erPvbfww, erPvbftata, erPvbfgaga;

//______________________________________________________________________________
//Double_t func(float x,float y,Double_t *par)
//{
// Double_t value=( (par[0]));
// return value;
//}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
//calculate chisquare
   Double_t chisq = 0;
   Double_t delta0,delta1,delta2,delta3,delta4,delta5,delta6,delta7,delta8,delta9,delta10;
   Double_t delta11,delta12,delta13,delta14,delta15,delta16,delta17;
   Double_t delta21,delta22,delta23,delta24,delta25,delta26,delta27,delta20;
   Double_t delta31,delta32,delta33,delta34,delta35,delta36,delta37;
   Double_t deltaP4,deltaP6,deltaP7,deltaP9,deltaPa,deltaPb;
   Double_t deltaPC1,deltaPC4,deltaPC5,deltaPC6,deltaPC7;
   Double_t deltaP4v, deltaP5v, deltaP7v;
   Double_t total;

   total = ((brhbb*par[1]*par[1]+brhcc*par[2]*par[2]+brhgg*par[3]*par[3]+brhww*par[4]*par[4]+brhtata*par[5]*par[5]+brhzz*par[6]*par[6]+brhgaga*par[7]*par[7]+brhmm*par[8]*par[8])*(1-0.01*par[0])+0.01*par[0]);

   delta0  = (1-par[6]*par[6]/total)/erzh;
   delta1  = (1-par[6]*par[6]*par[1]*par[1]/total)/erzhbb1;
   delta2  = (1-par[6]*par[6]*par[2]*par[2]/total)/erzhcc1;
   delta3  = (1-par[6]*par[6]*par[3]*par[3]/total)/erzhgg1;
   delta4  = (1-par[6]*par[6]*par[4]*par[4]/total)/erzhww1;
   delta5  = (1-par[6]*par[6]*par[5]*par[5]/total)/erzhtata1;
   delta6  = (1-par[6]*par[6]*par[6]*par[6]/total)/erzhzz1;
   delta7  = (1-par[6]*par[6]*par[7]*par[7]/total)/erzhgaga1;
   delta8  = (1-par[4]*par[4]*par[1]*par[1]/total)/ervhbb1;
   delta9  = (1-par[6]*par[6]*par[8]*par[8]/total)/erzhmm;
   delta10 = (par[6]*par[6]*par[0])/0.52;
   chisq += delta10*delta10;
//   chisq += (par[0]-total)/(total);

   chisq += delta0*delta0;
   chisq += delta1*delta1;
   chisq += delta2*delta2;
   chisq += delta3*delta3;
   chisq += delta4*delta4;
   chisq += delta5*delta5;
   chisq += delta6*delta6;
   chisq += delta7*delta7;
   chisq += delta9*delta9;
   chisq += delta8*delta8;
/*

//   deltaP4  = (1-par[3]*par[3]*par[4]*par[4]/total)/erPghww;
   deltaP6  = (1-par[3]*par[3]*par[6]*par[6]/total)/erPghzz;
//   deltaP7  = (1-par[3]*par[3]*par[7]*par[7]/total)/erPghgaga;
   deltaP9  = (1-par[3]*par[3]*par[8]*par[8]/total)/erPghmm;

//   chisq += deltaP4*deltaP4;
   chisq += deltaP6*deltaP6;
//   chisq += deltaP7*deltaP7;
   chisq += deltaP9*deltaP9;
*/
   deltaPC1  = (1-par[6]*par[6]*par[1]*par[1]/total)/erPC2vbfbb;
   deltaPC4  = (1-par[3]*par[3]*par[4]*par[4]/total)/erPC2ghww;
   deltaPC5  = (1-(0.5*par[4]*par[4]+0.5*par[3]*par[3])*par[5]*par[5]/total)/erPC2vbftata;
   deltaPC6  = (1-par[3]*par[3]*par[6]*par[6]/total)/erPC2ghzz;
   deltaPC7  = (1-par[3]*par[3]*par[7]*par[7]/total)/erPC2ghgaga;
//   deltaP4v  = (1-par[]*par[3]*par[4]*par[4]/par[0])/erPCghww;
//   deltaP5v  = (1-par[3]*par[3]*par[5]*par[5]/par[0])/erPCghtata;
//   deltaP7v  = (1-par[3]*par[3]*par[7]*par[7]/par[0])/erPCghgaga;



   chisq += deltaPC1*deltaPC1;
   chisq += deltaPC4*deltaPC4;
   chisq += deltaPC5*deltaPC5;
   chisq += deltaPC6*deltaPC6;
   chisq += deltaPC7*deltaPC7;

/*
   deltaPa = (1-par[4]*par[4]/par[6]/par[6])/erPwwzz;
   deltaPb = (1-par[7]*par[7]/par[6]/par[6])/erPgagazz;
   chisq += deltaPa*deltaPa;
   chisq += deltaPb*deltaPb;
*/


   delta11  = (1-par[6]*par[6]*par[1]*par[1]/total)/erzhbb2;
   delta12  = (1-par[6]*par[6]*par[2]*par[2]/total)/erzhcc2;
   delta13  = (1-par[6]*par[6]*par[3]*par[3]/total)/erzhgg2;
   delta14  = (1-par[6]*par[6]*par[4]*par[4]/total)/erzhww2;
   delta15  = (1-par[6]*par[6]*par[5]*par[5]/total)/erzhtata2;
   delta16  = (1-par[6]*par[6]*par[6]*par[6]/total)/erzhzz2;
   delta17  = (1-par[6]*par[6]*par[7]*par[7]/total)/erzhgaga2;

   chisq += delta11*delta11;
   chisq += delta12*delta12;
   chisq += delta13*delta13;
   chisq += delta14*delta14;
   chisq += delta15*delta15;
   chisq += delta16*delta16;
   chisq += delta17*delta17;

   delta21  = (1-par[4]*par[4]*par[1]*par[1]/total)/ervhbb;
   delta22  = (1-par[4]*par[4]*par[2]*par[2]/total)/ervhcc;
   delta23  = (1-par[4]*par[4]*par[3]*par[3]/total)/ervhgg;
   delta24  = (1-par[4]*par[4]*par[4]*par[4]/total)/ervhww;
   delta25  = (1-par[4]*par[4]*par[5]*par[5]/total)/ervhtata;
   delta26  = (1-par[4]*par[4]*par[6]*par[6]/total)/ervhzz;
   delta27  = (1-par[4]*par[4]*par[7]*par[7]/total)/ervhgaga;

   chisq += delta21*delta21;
   chisq += delta22*delta22;
   chisq += delta23*delta23;
   chisq += delta24*delta24;
   chisq += delta25*delta25;
   chisq += delta26*delta26;
   chisq += delta27*delta27;

   delta20  = (1-par[6]*par[6])/erzh2;
   chisq += delta20*delta20;

   delta31  = (1-par[4]*par[4]*par[1]*par[1]/total)/ervhbb1t;
   delta32  = (1-par[4]*par[4]*par[2]*par[2]/total)/ervhcc1t;
   delta33  = (1-par[4]*par[4]*par[3]*par[3]/total)/ervhgg1t;
   delta34  = (1-par[4]*par[4]*par[4]*par[4]/total)/ervhww1t;
   delta35  = (1-par[4]*par[4]*par[5]*par[5]/total)/ervhtata1t;
   delta36  = (1-par[4]*par[4]*par[6]*par[6]/total)/ervhzz1t;
   delta37  = (1-par[4]*par[4]*par[7]*par[7]/total)/ervhgaga1t;

   chisq += delta31*delta31;
   chisq += delta32*delta32;
   chisq += delta33*delta33;
   chisq += delta34*delta34;
   chisq += delta35*delta35;
   chisq += delta36*delta36;
   chisq += delta37*delta37;

   f = chisq;
}

//______________________________________________________________________________
void IfitSM_inv()
{
//SM Br for mh=120 GeV Higgs from Table 2.4.2 ILC DBD Physics Volumn

   brhbb=0.56;
   brhcc=0.028;
   brhgg=0.085;
   brhww=0.23;
   brhtata=0.062;
   brhzz=0.029;
   brhgaga=0.0023;

//LHC uncertainties on ggh with PDF/QCD scale uncertainties @300fb^-1, from ATL-PHYS-PUB-2012-004
   erPghzz=0.156;
   erPghww=0.289;
   erPghgaga=0.145;
   erPwwzz=0.254;
   erPgagazz=0.110;
//ZL combined 300 fbi
   erPCvbfbb=0.14;
   erPCghww=0.1039;
   erPCvbftata=0.1354;
   erPCghzz=0.0963;
   erPCghgaga=0.1153;
//3000 fbi
   erPC2vbfbb=0.07;
   erPC2ghww=0.0654;
   erPC2vbftata=0.0777;
   erPC2ghzz=0.0640;
   erPC2ghgaga=0.0748;

//DBD Table 2.5.3, errors on individual cross sections
//DBD Table 2.5.3, errors on individual cross sections

   erzh=0.026;
   erPvbftata=0.227;
   erzh2=0.0302;

//bb
   erzhbb1=0.012;
   erzhbb2=0.018;
   ervhbb=0.0066;
   ervhbb1=0.11;
   ervhbb1t=0.005;
//cc
   erzhcc1=0.083;
   erzhcc2=0.13;
   ervhcc=0.062;
   ervhcc1t=0.031;
//gg
   erzhgg1=0.07;
   erzhgg2=0.11;
   ervhgg=0.041;
   ervhgg1t=0.016;
//ww
   erzhww1=0.064;
   erzhww2=0.060;
   ervhww=0.024;
   ervhww1t=0.031;
//tata
   erzhtata1=0.042;
   erzhtata2=0.054;
   ervhtata=0.09;
   ervhtata1t=0.023;
//zz
   erzhzz1=0.19;
   erzhzz2=0.25;
   ervhzz=0.082;
   ervhzz1t=0.041;
//gaga
   erzhgaga1=0.38;
   erzhgaga2=0.38;
   ervhgaga=0.26;
   ervhgaga1t=0.085;

//mumu
   erzhmm=1.00;
   erPghmm=0.525;
   

   TMinuit *gMinuit = new TMinuit(9);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

// Set starting values and step sizes for parameters
   static Double_t vstart[9] = {1,1,1,1,1,1,1,1,1};
//   static Double_t vstart[9] = {0.9,1.1,0.9,1.1,0.9,1.1,0.9,1.1,0.9};
   static Double_t step[10] = {0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00005,0.00005};
   gMinuit->mnparm(0, "ghinv", vstart[0], step[0], 0,0,ierflg);
   gMinuit->mnparm(1, "ghbb", vstart[1], step[1], 0,0,ierflg);
   gMinuit->mnparm(2, "ghcc", vstart[2], step[2], 0,0,ierflg);
   gMinuit->mnparm(3, "ghgg", vstart[3], step[3], 0,0,ierflg);
   gMinuit->mnparm(4, "ghww", vstart[4], step[4], 0,0,ierflg);
   gMinuit->mnparm(5, "ghtata", vstart[5], step[5], 0,0,ierflg);
   gMinuit->mnparm(6, "ghzz", vstart[6], step[6], 0,0,ierflg);
   gMinuit->mnparm(7, "ghgaga", vstart[7], step[7], 0,0,ierflg);
   gMinuit->mnparm(8, "ghmm", vstart[8], step[8], 0,0,ierflg);
//   gMinuit->mnparm(9, "Br_inv", vstart[9], step[9], 0,0,ierflg);

// Now ready for minimization step
   arglist[0] = 50000000;
   arglist[1] = 1.;
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

// Print results
   Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   //gMinuit->mnprin(3,amin);

}
