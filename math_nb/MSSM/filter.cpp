#include<stdio.h>
#include<cstdlib>
#include<complex>
#include<time.h>
#include "CFeynHiggs.h"
#include "FHCouplings.h"
//#include "CSLHA.h"
//#include "FHRecord.h"
//#include "PDG.h"
//#include "SLHA.h"
//#include "SLHADefs.h"
//#include"../sources/micromegas.h"
//#include"../sources/micromegas_aux.h"
//#include"lib/pmodel.h"

extern "C" { 
  void initialize_higgsbounds_(int&, int&, char[5]);
  void finish_higgsbounds_();
  void higgsbounds_neutral_input_part_(double[3], double[3], int[3], double[3], double[3], double[3], double[3][3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3], double[3][3]);
  void  higgsbounds_charged_input_(double[1], double[1], double[1], double&, double[1], double[1], double[1], double[1]);
  void run_higgsbounds_(int&, int&, double&, int&);
}


double rnd(const double &min, const double &max){
  return min+(max-min)*(double)rand()/(double)RAND_MAX;
}

class FHPars {
public:
  int err, fail;
  /*                  FeynHiggs Parameters                 */
  //Input pars
  double TB;
  double MA0;
  double M3SQ;
  double M3SU;
  double MSUSY;
  double MUE;
  double At;
  double m1;
  double m2;
  double m3;
  //FHHiggsCorr
  int iHiggs;
  double MHiggs[4];
  double_complex SAeff, UHiggs[3][3], ZHiggs[3][3];
  //FHSelectUZ
  int uzint, uzext;
  //FHCouplings
  double_complex couplings[ncouplings], couplingsms[ncouplingsms];
  double gammas[ngammas], gammasms[ngammasms];
  //FHFlavour
  double bsgammaMSSM, bsgammaSM;
  double deltaMsMSSM, deltaMsSM;
  double  bsmumuMSSM, bsmumuSM;
  //FHHiggsProd
  double prodxs[nprodxs];
  double gghAA7,gghAA7SM, gghWW7,gghWW7SM, gghZZ7,gghZZ7SM;
  double ggH0AA7,ggH0AA7SM, ggH0WW7,ggH0WW7SM, ggH0ZZ7,ggH0ZZ7SM;
  double ggH1AA7,ggH1AA7SM, ggH1WW7,ggH1WW7SM, ggH1ZZ7,ggH1ZZ7SM;
  double ggh7[3], ggh8[3], ggh14[3];
  double gghSM7[3], gghSM8[3], gghSM14[3];
  double bbh7[3], bbh8[3], bbh14[3];
  double bbhSM7[3], bbhSM8[3], bbhSM14[3];
  double qqh14[3], Wh14[3], Zh14[3], tth14[3];
  double qqhSM14[3], WhSM14[3], ZhSM14[3], tthSM14[3];
  double StSth14[3], tHm14;
  //FHGetPara
  int nmfv;
  double MASf[4][6], MCha[2], MNeu[4];
  double_complex UASf[4][6][6], UCha[2][2], VCha[2][2], ZNeu[4][4];
  double_complex Deltab;
  double MGl;
  double MHtree[4], SAtree;
  //FHRetrieveSMPara
  double invAlfa, AlfasMZ, GF;
  double ME, MU, MD, MM, MC, MS, ML, MB;
  double MW, MZ;
  double CKMlambda, CKMA, CKMrhobar, CKMetabar;
  /*                      HiggsBound Parameters                  */
  //Neutral input part
  double Mh[3], GammaTotal_hj[3];
  int CP_value[3];
  double CS_lep_hjZ_ratio[3],
    CS_lep_bbhj_ratio[3], CS_lep_tautauhj_ratio[3],
    CS_lep_hjhi_ratio_nHbynH[3][3],
    CS_gg_hj_ratio[3], CS_bb_hj_ratio[3],
    CS_bg_hjb_ratio[3],
    CS_ud_hjWp_ratio[3], CS_cs_hjWp_ratio[3],
    CS_ud_hjWm_ratio[3], CS_cs_hjWm_ratio[3],
    CS_gg_hjZ_ratio[3],
    CS_dd_hjZ_ratio[3], CS_uu_hjZ_ratio[3],
    CS_ss_hjZ_ratio[3], CS_cc_hjZ_ratio[3],
    CS_bb_hjZ_ratio[3],
    CS_tev_vbf_ratio[3], CS_tev_tthj_ratio[3],
    CS_lhc7_vbf_ratio[3], CS_lhc7_tthj_ratio[3],
    BR_hjss[3], BR_hjcc[3],
    BR_hjbb[3], BR_hjmumu[3], BR_hjtautau[3],
    BR_hjWW[3], BR_hjZZ[3], BR_hjZga[3],
    BR_hjgaga[3], BR_hjgg[3],
    BR_hjinvisible[3], BR_hjhihi_nHbynH[3][3];
  //Charged input part
  double Mhplus[1], GammaTotal_Hpj[1],
    CS_lep_HpjHmj_ratio[1],
    BR_tWpb, BR_tHpjb[1],
    BR_Hpjcs[1], BR_Hpjcb[1], BR_Hpjtaunu[1]; 
  //run_HiggsBounds
  int HBresultL, chanL, ncombinedL;
  double obsratioL;
  int HBresultH, chanH, ncombinedH;
  double obsratioH;
  //misc:
  double norm, CW2, Pi;
  double g2hjbb[3], g2hjWW[3], g2hjZZ[3],
    g2hjgg[3], g2hjhiZ_nHbynH[3][3];
  double g2hjbb_s[3], g2hjbb_p[3];
  double g2hjtautau_s[3], g2hjtautau_p[3];
  int sneutrino_lspcandidate_number;
  bool invisible_lsp;
  double lspcandidate_mass;
  //For exporting
  char filename[100];
  FILE *fp1, *fp2;

  FHPars(char *flnm):
    err(0),MSUSY(3000),Pi(3.1415926535897932384626433832795029){
    sprintf(filename,"%s",flnm);
  };
  FHPars(FILE *fp1tmp, FILE *fp2tmp):
    err(0),MSUSY(3000),Pi(3.1415926535897932384626433832795029){
    fp1=fp1tmp;
    fp2=fp2tmp;
    sprintf(filename,"");
  };
  
  void checkHiggsMass(double Mmin, double Mmax){
    iHiggs=1;
    //fprintf(stderr,"Mh=%G\tMH=%G\tMA=%G\tMHp=%G\n",MHiggs[0],MHiggs[1],MHiggs[2],MHiggs[3]);
    if(Mmin<MHiggs[0] && MHiggs[0]<Mmax && Mmin<MHiggs[1] && MHiggs[1]<Mmax) iHiggs=3;
    else if(Mmin<MHiggs[0] && MHiggs[0]<Mmax) iHiggs=1;
    else if(Mmin<MHiggs[1] && MHiggs[1]<Mmax) iHiggs=2;
    else if(fail==0) fail=3;
  }
  
  void getOtherPars(){
    FHGetPara(&err,&nmfv,MASf,UASf,MCha,UCha,VCha,MNeu,ZNeu,&Deltab,&MGl,MHtree,&SAtree);
  }
  void getSMPars(){
    FHRetrieveSMPara(&err,&invAlfa,&AlfasMZ,&GF,&ME,&MU,&MD,&MM,&MC,&MS,&ML,&MB,&MW,&MZ,&CKMlambda,&CKMA,&CKMrhobar,&CKMetabar);
  }
  
  void calcCouplings(){
    uzint=2;
    uzext=2;
    FHSelectUZ(&err,uzint,uzext);
    FHCouplings(&err,couplings,couplingsms,gammas,gammasms,1);
    if(err==0) FHFlavour(&err, &bsgammaMSSM, &bsgammaSM, &deltaMsMSSM, &deltaMsSM, &bsmumuMSSM, &bsmumuSM);
    if(err==0) bsgammaMSSM-=7.6e-5;
  }

  void calcProd(double sqrts){
    FHHiggsProd(&err,sqrts,prodxs);
    if(sqrts==7)
      for(int i=0;i<3;i++){
	ggh7[i]=ggh(i+1);
	gghSM7[i]=gghSM(i+1);
	bbh7[i]=bbh(i+1);
	bbhSM7[i]=bbhSM(i+1);
      }
    else if(sqrts==8)
      for(int i=0;i<3;i++){
	ggh8[i]=ggh(i+1);
	gghSM8[i]=gghSM(i+1);
	bbh8[i]=bbh(i+1);
	bbhSM8[i]=bbhSM(i+1);
      }
    else if(sqrts==14){
      for(int i=0;i<3;i++){
	ggh14[i]=ggh(i+1);
	gghSM14[i]=gghSM(i+1);
	bbh14[i]=bbh(i+1);
	bbhSM14[i]=bbhSM(i+1);
	qqh14[i]=qqh(i+1);
	qqhSM14[i]=qqhSM(i+1);
	Wh14[i]=Wh(i+1);
	WhSM14[i]=WhSM(i+1);
	Zh14[i]=Zh(i+1);
	ZhSM14[i]=ZhSM(i+1);
	tth14[i]=tth(i+1);
	tthSM14[i]=tthSM(i+1);
	StSth14[i]=StSth(i+1);
      }
      tHm14=tHm;
    }
  }

  void FH2HB(){
    /*             Taken from HBwithFH.F
     *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *      
     c Set variables needed by HiggsBounds (using results from FeynHiggs).
     c See HiggsBounds documentation for definition of variables used
     c as arguments to HiggsBounds_neutral_input_part and run_HiggsBounds
     c and FeynHiggs documentation for all other variables.
     
     c Note: It is slightly more accurate to use the subroutine HiggsBounds_neutral_input_part
     c rather than the subroutine HiggsBounds_neutral_input_effC because the SM branching ratios
     c used internally in HiggsBounds (from HDecay) are not identical to the SM branching
     c ratios used in FeynHiggs
    */

    
    for(int i=1;i<=3;i++){

	Mh[i-1]=MHiggs[i-1];
	GammaTotal_hj[i-1] = GammaTot(i);

	BR_hjss[i-1]       = BR(H0FF(i,4,2,2));
	BR_hjcc[i-1]       = BR(H0FF(i,3,2,2));
	BR_hjbb[i-1]       = BR(H0FF(i,4,3,3));
	BR_hjmumu[i-1]     = BR(H0FF(i,2,2,2));
	BR_hjtautau[i-1]   = BR(H0FF(i,2,3,3)); 
	
	BR_hjWW[i-1]     = BR(H0VV(i,4));
	BR_hjgaga[i-1]   = BR(H0VV(i,1));
	BR_hjZga[i-1]    = BR(H0VV(i,2));
	BR_hjZZ[i-1]     = BR(H0VV(i,3));
	BR_hjgg[i-1]     = BR(H0VV(i,5));
	
	if(GammaSM(H0FF(i,4,3,3))<0)
          g2hjbb[i-1]=0;
	else
          g2hjbb[i-1]=Gamma(H0FF(i,4,3,3))/GammaSM(H0FF(i,4,3,3));

	/*
	  c Note that this is currently equivalent to
	  c         g2hjbb(i)= bbh(i)/bbhSM(i)
	  c         g2hjbb(i)= btagbh(i)/btagbhSM(i)
	  c as long as MH>80 GeV
	*/
	CS_bg_hjb_ratio[i-1] = g2hjbb[i-1];
	CS_bb_hj_ratio[i-1]  = g2hjbb[i-1];
	
	g2hjbb_s[i-1]=pow(abs(RCoupling(H0FF(i,4,3,3))/RCouplingSM(H0FF(i,4,3,3)) + LCoupling(H0FF(i,4,3,3))/LCouplingSM(H0FF(i,4,3,3)))/2.0,2.0);
	g2hjbb_p[i-1]=pow(abs(RCoupling(H0FF(i,4,3,3))/RCouplingSM(H0FF(i,4,3,3)) - LCoupling(H0FF(i,4,3,3))/LCouplingSM(H0FF(i,4,3,3)))/2.0,2.0);
	g2hjtautau_s[i-1]=pow(abs(RCoupling(H0FF(i,2,3,3))/RCouplingSM(H0FF(i,2,3,3)) + LCoupling(H0FF(i,2,3,3))/LCouplingSM(H0FF(i,2,3,3)))/2.0,2.0);
	g2hjtautau_p[i-1]=pow(abs(RCoupling(H0FF(i,2,3,3))/RCouplingSM(H0FF(i,2,3,3)) - LCoupling(H0FF(i,2,3,3))/LCouplingSM(H0FF(i,2,3,3)))/2.0,2.0);

	if(        g2hjbb_p[i-1]<1.0E-10)
	  CP_value[i-1] = 1;
	else if(   g2hjbb_s[i-1]<1.0E-10)
	  CP_value[i-1] = -1;
	else
	  CP_value[i-1] = 0;


	CS_lep_bbhj_ratio[i-1]     = g2hjbb_s[i-1]+g2hjbb_p[i-1];
	CS_lep_tautauhj_ratio[i-1] = g2hjtautau_s[i-1]+g2hjtautau_p[i-1];
	
	g2hjWW[i-1]= pow(real( Coupling(H0VV(i,4)) / CouplingSM(H0VV(i,4)) ),2.0)
	           + pow(imag( Coupling(H0VV(i,4)) / CouplingSM(H0VV(i,4)) ),2.0);
	/*c Note that this is currently equivalent to
	  c         g2hjWW(i)= WhTev(i)/WhTevSM(i
	  c	  g2hjWW(i)= qqhTev(i)/qqhTevSM(i)
	  c as long as MH>80 GeV and uzint=uzext
	*/
	g2hjZZ[i-1]= pow(real(  Coupling(H0VV(i,3)) / CouplingSM(H0VV(i,3)) ),2.0)
	         + pow(imag( Coupling(H0VV(i,3)) / CouplingSM(H0VV(i,3)) ),2.0);
	/*c Note that this is currently equivalent to
	  c         g2hjZZ(i)= ZhTev(i)/ZhTevSM(i)
	  c as long as MH>80 GeV and uzint=uzext
	  c It is also equivalent to g2hjWW(i)
	*/
	CS_lep_hjZ_ratio[i-1]        = g2hjZZ[i-1];
	
	CS_gg_hjZ_ratio[i-1]     = 0.0;
	CS_dd_hjZ_ratio[i-1]     = g2hjZZ[i-1];
	CS_uu_hjZ_ratio[i-1]     = g2hjZZ[i-1];
	CS_ss_hjZ_ratio[i-1]     = g2hjZZ[i-1];
	CS_cc_hjZ_ratio[i-1]     = g2hjZZ[i-1];
	CS_bb_hjZ_ratio[i-1]     = g2hjZZ[i-1];
	
	CS_ud_hjWp_ratio[i-1]    = g2hjZZ[i-1];
	CS_cs_hjWp_ratio[i-1]    = g2hjZZ[i-1];
	CS_ud_hjWm_ratio[i-1]    = g2hjZZ[i-1];
	CS_cs_hjWm_ratio[i-1]    = g2hjZZ[i-1];

	CS_tev_vbf_ratio[i-1]    = g2hjZZ[i-1];
	CS_lhc7_vbf_ratio[i-1]   = g2hjZZ[i-1];


	if(tthSM(i)>0)
	  CS_tev_tthj_ratio[i-1]    = tth(i)/tthSM(i);  
	else
	  CS_tev_tthj_ratio[i-1]    = 0.0;

	CS_lhc7_tthj_ratio[i-1] = CS_tev_tthj_ratio[i-1];

	//c tevatron gluon fusion XS is not calculated in FH is MH<90 geV
	if(Mh[i-1]>90){
	  if(gghSM(i)>0)
            CS_gg_hj_ratio[i-1] = ggh(i)/gghSM(i);
	  else
            CS_gg_hj_ratio[i-1] = 0.0;
	}
	else{
	  if(GammaSM(H0VV(i,5))<0)
            CS_gg_hj_ratio[i-1]=0.0;
	  else
            CS_gg_hj_ratio[i-1]= Gamma(H0VV(i,5))/GammaSM(H0VV(i,5));
	}

    }
    norm=GF*sqrt(2.0)*MZ*MZ;
    
    for(int j=1;j<=3;j++)
      for(int i=1;i<=3;i++){
	g2hjhiZ_nHbynH[i-1][j-1]= ( pow(real(  Coupling(H0HV(j,i)) ),2.0) + pow(imag( Coupling(H0HV(j,i)) ),2.0) )/norm;
	CS_lep_hjhi_ratio_nHbynH[i-1][j-1] = g2hjhiZ_nHbynH[i-1][j-1];
	BR_hjhihi_nHbynH[i-1][j-1]=BR(H0HH(j,i,i));
      }

    /*c higgs->neutralino1 neutralino1 contributes the invisible Higgs decay width
      c when neutralino1 or sneutrino is the LSP
    */
    for(int i=1;i<=3;i++){
      sneutrino_lspcandidate_number=0;
      invisible_lsp=true;
      
      /*c first determine whether lightest sneutrino is lighter than the lightest neutralino
	c
	c sneutrino_lspcandidate_number=0 indicates that lightest neutralino is 
	c lighter than all the sneutrinos
      */
      lspcandidate_mass=MNeu[1-1];
      for(int as=1;as<=3;as++)
	if( MASf[1-1][as-1] < lspcandidate_mass ){
	  lspcandidate_mass=MASf[1-1][as-1];
	  sneutrino_lspcandidate_number=as;
	}
      
      
      if( MCha[1-1] < lspcandidate_mass )
	invisible_lsp=false;
      else if( MGl < lspcandidate_mass )
	invisible_lsp=false;
      else
	for(int as=1;as<=6;as++)
	  for(int t=2;t<=4;t++)
	    if( MASf[t-1][as-1] < lspcandidate_mass )
	      invisible_lsp=false;
      
      
      
      if(invisible_lsp){
	if(sneutrino_lspcandidate_number==0)
	  BR_hjinvisible[i-1] = BR(H0NeuNeu(i,1,1));
	else
	  BR_hjinvisible[i-1] = BR(H0SfSf(i,1,1,1,6));
      }
      else
	BR_hjinvisible[i-1] = 0.0;
    
      
    }

    //*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

    Mhplus[1-1]              = MHiggs[4-1]; 
    GammaTotal_Hpj[1-1]      = GammaTot(4);  
    CS_lep_HpjHmj_ratio[1-1] = 1.0;
    BR_tWpb                  = BR( tBF(1) );
    BR_tHpjb[1-1]            = BR( tBF(2) );
    BR_Hpjcs[1-1]            = BR( HpFF(2,2,2) );
    BR_Hpjcb[1-1]            = BR( HpFF(2,2,3) );
    BR_Hpjtaunu[1-1]         = BR( HpFF(1,3,3) );

    //*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

    /*c calls to HiggsBounds_neutral_input_part,HiggsBounds_charged_input,
      c which give input to HiggsBounds
      
      print*,'calling HiggsBounds_neutral_input_part'
    */
    higgsbounds_neutral_input_part_(Mh,GammaTotal_hj,CP_value, 
				   CS_lep_hjZ_ratio,
				   CS_lep_bbhj_ratio,CS_lep_tautauhj_ratio,
				   CS_lep_hjhi_ratio_nHbynH,
				   CS_gg_hj_ratio,CS_bb_hj_ratio,
				   CS_bg_hjb_ratio,
				   CS_ud_hjWp_ratio,CS_cs_hjWp_ratio,
				   CS_ud_hjWm_ratio,CS_cs_hjWm_ratio,
				   CS_gg_hjZ_ratio,
				   CS_dd_hjZ_ratio,CS_uu_hjZ_ratio,
				   CS_ss_hjZ_ratio,CS_cc_hjZ_ratio,
				   CS_bb_hjZ_ratio,
				   CS_tev_vbf_ratio,CS_tev_tthj_ratio,
				   CS_lhc7_vbf_ratio,CS_lhc7_tthj_ratio,
				   BR_hjss,BR_hjcc,
				   BR_hjbb,BR_hjmumu,BR_hjtautau,
				   BR_hjWW,BR_hjZZ,BR_hjZga, BR_hjgaga,BR_hjgg,
				   BR_hjinvisible,BR_hjhihi_nHbynH              );
     
      //print*,'calling HiggsBounds_charged_input'
      higgsbounds_charged_input_(Mhplus,GammaTotal_Hpj,
				CS_lep_HpjHmj_ratio,
				BR_tWpb,BR_tHpjb,
				BR_Hpjcs,BR_Hpjcb,BR_Hpjtaunu);


      //*  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
    
  }

  void checkHB(){
    int nHiggsneut=3, nHiggsplus=1;
    char whichanalyses[10];
    
    sprintf(whichanalyses,"%s ","onlyL");
    initialize_higgsbounds_(nHiggsneut, nHiggsplus, whichanalyses);
    FH2HB();
    run_higgsbounds_(HBresultL ,chanL, obsratioL, ncombinedL);
    finish_higgsbounds_();    
    if(HBresultL==-1) err=1;
    if(HBresultL==0&&fail==0) fail=1;

    sprintf(whichanalyses,"%s ","onlyH");
    initialize_higgsbounds_(nHiggsneut, nHiggsplus, whichanalyses);
    FH2HB();
    run_higgsbounds_(HBresultH ,chanH, obsratioH, ncombinedH);
    finish_higgsbounds_();    
    if(HBresultH==-1) err=1;
    if(HBresultH==0&&fail==0) fail=2;


  }
    
  void checkLHC(){
    calcProd(7);
    if(err==0){
      //    gg   ->   h   ->   AA
      double gghAAmin=0.8, gghAAmax=30;
      ggH0AA7=ggh(1)*BR(H0VV(1,1));
      ggH0AA7SM=gghSM(1)*BRSM(H0VV(1,1));
      ggH1AA7=ggh(2)*BR(H0VV(2,1));
      ggH1AA7SM=gghSM(2)*BRSM(H0VV(2,1));
      if(iHiggs==3){
	gghAA7=ggH0AA7+ggH1AA7;
	gghAA7SM=(ggH0AA7SM+ggH1AA7SM)/2;
      }
      else{
	gghAA7=ggh(iHiggs)*BR(H0VV(iHiggs,1));
	gghAA7SM=gghSM(iHiggs)*BRSM(H0VV(iHiggs,1));
      }
      if(gghAAmin>gghAA7/gghAA7SM || gghAA7/gghAA7SM>gghAAmax) if(fail==0) fail=4;
      
      //    gg   ->   h   ->   ZZ
      double gghZZmin=0, gghZZmax=30;
      ggH0ZZ7=ggh(1)*BR(H0VV(1,3));
      ggH0ZZ7SM=gghSM(1)*BRSM(H0VV(1,3));
      ggH1ZZ7=ggh(2)*BR(H0VV(2,3));
      ggH1ZZ7SM=gghSM(2)*BRSM(H0VV(2,3));
      if(iHiggs==3){
	gghZZ7=ggH0ZZ7+ggH1ZZ7;
	gghZZ7SM=(ggH0ZZ7SM+ggH1ZZ7SM)/2;
      }
      else{
	gghZZ7=ggh(iHiggs)*BR(H0VV(iHiggs,3));
	gghZZ7SM=gghSM(iHiggs)*BRSM(H0VV(iHiggs,3));
      }
      if(gghZZmin>gghZZ7/gghZZ7SM || gghZZ7/gghZZ7SM>gghZZmax) if(fail==0) fail=4;
      
      //    gg   ->   h   ->   WW
      double gghWWmin=0, gghWWmax=30;
      ggH0WW7=ggh(1)*BR(H0VV(1,4));
      ggH0WW7SM=gghSM(1)*BRSM(H0VV(1,4));
      ggH1WW7=ggh(2)*BR(H0VV(2,4));
      ggH1WW7SM=gghSM(2)*BRSM(H0VV(2,4));
      if(iHiggs==3){
	gghWW7=ggH0WW7+ggH1WW7;
	gghWW7SM=(ggH0WW7SM+ggH1WW7SM)/2;
      }
      else{
	gghWW7=ggh(iHiggs)*BR(H0VV(iHiggs,4));
	gghWW7SM=gghSM(iHiggs)*BRSM(H0VV(iHiggs,4));
      }
      if(gghWWmin>gghWW7/gghWW7SM || gghWW7/gghWW7SM>gghWWmax) if(fail==0) fail=4; 
    }
  }

  void checkBSGamma(){
    if(bsgammaMSSM<2.79e-4 || bsgammaMSSM>4.31e-4) if(fail==0) fail=5;
  }

  void exprt(){
    FILE *fp;
    if(fail==0){
      fprintf(stdout,"fail=%d\n",fail);
      fprintf(stdout,"TB=%G\tMA0=%G\tM3SQ=%G\tM3SU=%G\tMSUSY=%G\tMUE=%G\tAt=%G\n",TB,MA0,M3SQ,M3SU,MSUSY,MUE,At);
      fprintf(stdout,"\tMh=%G\tMH=%G\tMA=%G\tMHp=%G\n",MHiggs[0],MHiggs[1],MHiggs[2],MHiggs[3]);
      fprintf(stdout,"\tgghAA7/gghAA7SM=%G\tgghZZ7/gghZZ7SM=%G\tgghWW7/gghWW7SM=%G\n",
	      gghAA7/gghAA7SM, gghZZ7/gghZZ7SM, gghWW7/gghWW7SM);
      fprintf(stdout,"\tbsgamma/bsgammaSM=%G\n",bsgammaMSSM/bsgammaSM);
    }

    if(strcmp(filename,"")==0) fp=fp1;
    else fp = fopen(filename,"a");
    fprintf(fp,"{");
    fprintf(fp,"TB->%G, MA0->%G, M3SQ->%G, M3SU->%G, MSUSY->%G, MUE->%G, At->%G, ",TB,MA0,M3SQ,M3SU,MSUSY,MUE,At);
    fprintf(fp,"fail->%d, HBresultL->%d, chanL->%d, obsratioL->%G, ncombinedL->%d, HBresultH->%d, chanH->%d, obsratioH->%G, ncombinedH->%d, ",fail,HBresultL,chanL, obsratioL, ncombinedL,HBresultH,chanH, obsratioH, ncombinedH);
    fprintf(fp,"Mh->%G, MH->%G, MA->%G, MHp->%G,",MHiggs[0],MHiggs[1],MHiggs[2],MHiggs[3]);
    fprintf(fp,"MNeu1->%G, MNeu2->%G, MNeu3->%G, MNeu4->%G,",MNeu[0],MNeu[1],MNeu[2],MNeu[3]);
    fprintf(fp,"MCha1->%G, MCha2->%G, ",MCha[0],MCha[1]);
    fprintf(fp,"Mst1->%G, Mst2->%G, ",MASf[3][3],MASf[3][6]);
    fprintf(fp,"gghAAR->%G, gghZZR->%G, gghWWR->%G, ",gghAA7/gghAA7SM,gghZZ7/gghZZ7SM,gghWW7/gghWW7SM);
    fprintf(fp,"ggH0AAR->%G, ggH0ZZR->%G, ggH0WWR->%G, ",ggH0AA7/ggH0AA7SM,ggH0ZZ7/ggH0ZZ7SM,ggH0WW7/ggH0WW7SM);
    fprintf(fp,"ggH1AAR->%G, ggH1ZZR->%G, ggH1WWR->%G, ",ggH1AA7/ggH1AA7SM,ggH1ZZ7/ggH1ZZ7SM,ggH1WW7/ggH1WW7SM);
    fprintf(fp,"bsgammaMSSM->%G, bsgammaSM->%G, ",bsgammaMSSM,bsgammaSM);
    fprintf(fp,"ggh07->%G, ggH07->%G, ggA7->%G, bbh07->%G, bbH07->%G, bbA7->%G, ",ggh7[0],ggh7[1],ggh7[2],bbh7[0],bbh7[1],bbh7[2]);
    fprintf(fp,"ggh08->%G, ggH08->%G, ggA8->%G, bbh08->%G, bbH08->%G, bbA8->%G, ",ggh8[0],ggh8[1],ggh8[2],bbh8[0],bbh8[1],bbh8[2]);
    fprintf(fp,"ggh014->%G, ggH014->%G, ggA14->%G, bbh014->%G, bbH014->%G, bbA14->%G, ",ggh14[0],ggh14[1],ggh14[2],bbh14[0],bbh14[1],bbh14[2]);
    fprintf(fp,"qqh014->%G, qqH014->%G, qqA14->%G, Wh014->%G, WH014->%G, WA14->%G, Zh014->%G, ZH014->%G, ZA14->%G, ",qqh14[0],qqh14[1],qqh14[2],Wh14[0],Wh14[1],Wh14[2],Zh14[0],Zh14[1],Zh14[2]);
    fprintf(fp,"tth014->%G, ttH014->%G, ttA14->%G, StSth014->%G, StStH014->%G, StStA14->%G, tHm14->%G, ",tth14[0],tth14[1],tth14[2],StSth14[0],StSth14[1],StSth14[2],tHm14);
    fprintf(fp,"BRh0tt->%G, BRH0tt->%G, BRAtt->%G, BRh0bb->%G, BRH0bb->%G, BRAbb->%G, BRh0tautau->%G, BRH0tautau->%G, BRAtautau->%G, ",BR(H0FF(1,3,3,3)),BR(H0FF(2,3,3,3)),BR(H0FF(3,3,3,3)),BR(H0FF(1,4,3,3)),BR(H0FF(2,4,3,3)),BR(H0FF(3,4,3,3)),BR(H0FF(1,2,3,3)),BR(H0FF(2,2,3,3)),BR(H0FF(3,2,3,3)));
    fprintf(fp,"BRh0gaga->%G, BRH0gaga->%G, BRAgaga->%G, BRh0gaZ->%G, BRH0gaZ->%G, BRAgaZ->%G, BRh0ZZ->%G, BRH0ZZ->%G, BRAZZ->%G, BRh0WW->%G, BRH0WW->%G, BRAWW->%G, BRh0gg->%G, BRH0gg->%G, BRAgg->%G, ",BR(H0VV(1,1)),BR(H0VV(2,1)),BR(H0VV(3,1)),BR(H0VV(1,2)),BR(H0VV(2,2)),BR(H0VV(3,2)),BR(H0VV(1,3)),BR(H0VV(2,3)),BR(H0VV(3,3)),BR(H0VV(1,4)),BR(H0VV(2,4)),BR(H0VV(3,4)),BR(H0VV(1,5)),BR(H0VV(2,5)),BR(H0VV(3,5)));
    fprintf(fp,"BRh0st1st1->%G, BRH0st1st1->%G, BRAst1st1->%G, BRh0st1st2->%G, BRH0st1st2->%G, BRAst1st2->%G, BRh0st2st2->%G, BRH0st2st2->%G, BRAst2st2->%G, ",BR(H0SfSf(1,1,1,3,3)),BR(H0SfSf(2,1,1,3,3)),BR(H0SfSf(3,1,1,3,3)),BR(H0SfSf(1,1,2,3,3)),BR(H0SfSf(2,1,2,3,3)),BR(H0SfSf(3,1,2,3,3)),BR(H0SfSf(1,2,2,3,3)),BR(H0SfSf(2,2,2,3,3)),BR(H0SfSf(3,2,2,3,3)));
    fprintf(fp,"BRhptb->%G, BRhptaunu->%G, ",BR(HpFF(2,3,3)),BR(HpFF(1,3,3)));
    fprintf(fp,"gZAh0->%G, gZAH0->%G, gWAHp->%G",Coupling(H0HV(3,1)),Coupling(H0HV(3,2)),Coupling(HpHV(3)));
    fprintf(fp,"},\n");
    if(strcmp(filename,"")!=0) fclose(fp);
    

    if(strcmp(filename,"")==0) fp=fp2;
    else fp = fopen(filename,"a");
    fprintf(fp,"%.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E", MA0, TB, MSUSY, M3SQ, M3SU, MSUSY, MSUSY, MSUSY, At, 0., 0., MUE, m1, m2, m3);//15
    fprintf(fp,"  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E",  MASf[3-1][3-1], MASf[3-1][6-1], real(UASf[3-1][3-1][3-1]), real(UASf[3-1][6-1][3-1]), MASf[4-1][3-1], MASf[4-1][6-1], real(UASf[4-1][3-1][3-1]), real(UASf[4-1][6-1][3-1]), MASf[2-1][3-1], MASf[2-1][6-1], real(UASf[2-1][3-1][3-1]), real(UASf[2-1][6-1][3-1]), MASf[1-1][3-1], MCha[1-1], MCha[2-1], MNeu[1-1], MNeu[2-1], MNeu[3-1], MNeu[4-1], MGl);//20+15=35
    fprintf(fp,"  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E", MHiggs[0], MHiggs[1], MHiggs[2], MHiggs[3], real(SAeff), real(Deltab), bsgammaMSSM, bsgammaSM, bsmumuMSSM, bsmumuSM);//10+35=45
    fprintf(fp,"  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E", ggh(1),gghSM(1), ggh(2), gghSM(2), ggh(3), gghSM(3), bbh(1), bbhSM(1), bbh(2), bbhSM(2), bbh(3), bbhSM(3));//12+45=57
    fprintf(fp,"  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E", CS_gg_hj_ratio[1-1],CS_gg_hj_ratio[2-1], CS_gg_hj_ratio[3-1], CS_bb_hj_ratio[1-1], CS_bb_hj_ratio[2-1], CS_bb_hj_ratio[3-1], CS_lhc7_vbf_ratio[1-1], CS_lhc7_vbf_ratio[2-1], CS_lhc7_vbf_ratio[3-1], CS_dd_hjZ_ratio[1-1], CS_dd_hjZ_ratio[2-1], CS_dd_hjZ_ratio[3-1], CS_ud_hjWp_ratio[1-1], CS_ud_hjWp_ratio[2-1], CS_ud_hjWp_ratio[3-1], CS_lhc7_tthj_ratio[1-1], CS_lhc7_tthj_ratio[2-1], CS_lhc7_tthj_ratio[3-1]);//57+18=75
    fprintf(fp,"  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E",GammaTot(1), GammaSMTot(1), GammaTot(2), GammaSMTot(2),GammaTot(3), GammaSMTot(3), GammaTot(4));//7+75=82
    fprintf(fp,
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E",
	    BR(H0FF(1,3,3,3)), BR(H0FF(2,3,3,3)), BR(H0FF(3,3,3,3)), 
	    BR(H0FF(1,4,3,3)), BR(H0FF(2,4,3,3)), BR(H0FF(3,4,3,3)),
	    BR(H0FF(1,2,3,3)), BR(H0FF(2,2,3,3)), BR(H0FF(3,2,3,3)),
	    BR(H0VV(1,1)), BR(H0VV(2,1)), BR(H0VV(3,1)),
	    BR(H0VV(1,2)), BR(H0VV(2,2)), BR(H0VV(3,2)),
	    BR(H0VV(1,3)), BR(H0VV(2,3)), BR(H0VV(3,3)),
	    BR(H0VV(1,4)), BR(H0VV(2,4)), BR(H0VV(3,4)),
	    BR(H0VV(1,5)), BR(H0VV(2,5)), BR(H0VV(3,5)),
	    BR(H0SfSf(1,1,1,3,3)), BR(H0SfSf(2,1,1,3,3)), BR(H0SfSf(3,1,1,3,3)),
	    BR(H0SfSf(1,1,2,3,3)), BR(H0SfSf(2,1,2,3,3)), BR(H0SfSf(3,1,1,3,3)),
	    BR(HpFF(2,3,3)),BR(HpFF(1,3,3)),BR(HpFF(2,2,2)),
	    BRSM(H0FF(1,4,3,3)), BRSM(H0FF(2,4,3,3)), BRSM(H0FF(3,4,3,3)),
	    BRSM(H0FF(1,2,3,3)), BRSM(H0FF(2,2,3,3)), BRSM(H0FF(3,2,3,3)),
	    BRSM(H0VV(1,1)), BRSM(H0VV(2,1)), BRSM(H0VV(3,1)), 
	    BRSM(H0VV(1,3)), BRSM(H0VV(2,3)), BRSM(H0VV(3,3)),
	    BRSM(H0VV(1,4)), BRSM(H0VV(2,4)), BRSM(H0VV(3,4))
	    );//48+82=130
    fprintf(fp,"  %d  %d  %.4E  %d  %d  %d  %.4E  %d", HBresultL, chanL, obsratioL, ncombinedL, HBresultH, chanH, obsratioH, ncombinedH);//8+130=138
    fprintf(fp,
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E",
	    ggh7[0],gghSM7[0],ggh7[1],gghSM7[1],ggh7[2],gghSM7[2],
	    bbh7[0],bbhSM7[0],bbh7[1],bbhSM7[1],bbh7[2],bbhSM7[2],
	    ggh8[0],gghSM8[0],ggh8[1],gghSM8[1],ggh8[2],gghSM8[2],
	    bbh8[0],bbhSM8[0],bbh8[1],bbhSM8[1],bbh8[2],bbhSM8[2],
	    ggh14[0],gghSM14[0],ggh14[1],gghSM14[1],ggh14[2],gghSM14[2],
	    bbh14[0],bbhSM14[0],bbh14[1],bbhSM14[1],bbh14[2],bbhSM14[2],
	    qqh14[0],qqhSM14[0],qqh14[1],qqhSM14[1],qqh14[2],qqhSM14[2],
	    Wh14[0],WhSM14[0],Wh14[1],WhSM14[1],Wh14[2],WhSM14[2],
	    Zh14[0],ZhSM14[0],Zh14[1],ZhSM14[1],Zh14[2],ZhSM14[2],
	    tth14[0],tthSM14[0],tth14[1],tthSM14[1],tth14[2],tthSM14[2],
	    StSth14[0],StSth14[1],StSth14[2],tHm14);//64+138=202
    fprintf(fp,
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E"
	    "  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E  %.4E",
	    real(Coupling(H0HV(1,1))), imag(Coupling(H0HV(1,1))),
	    real(Coupling(H0HV(1,2))), imag(Coupling(H0HV(1,2))),
	    real(Coupling(H0HV(1,3))), imag(Coupling(H0HV(1,3))),
	    real(Coupling(H0HV(2,1))), imag(Coupling(H0HV(2,1))), 
	    real(Coupling(H0HV(2,2))),  imag(Coupling(H0HV(2,2))),
	    real(Coupling(H0HV(2,3))), imag(Coupling(H0HV(2,3))),
	    real(Coupling(H0HV(3,1))), imag(Coupling(H0HV(3,1))), 
	    real(Coupling(H0HV(3,2))), imag(Coupling(H0HV(3,2))), 
	    real(Coupling(H0HV(3,3))), imag(Coupling(H0HV(3,3))),
	    real(Coupling(HpHV(1))), imag(Coupling(HpHV(1))), 
	    real(Coupling(HpHV(2))), imag(Coupling(HpHV(2))), 
	    real(Coupling(HpHV(3))), imag(Coupling(HpHV(3)))
	    );
    fprintf(fp,"\n");
    if(strcmp(filename,"")!=0) fclose(fp);
  }

  /*
    fail codes:
    0 : pass all
    1 : fail HB LEP
    2 : fail HB Had
    3 : fail Mass
    4 : fail cs
   */


  int next(double m1, double m2, double m3, double mu, double tb, double ma, double m3q, double m3u, double at, int lvl){
  FHSetPara(&err,1,173.1,tb,ma,-1,MSUSY,MSUSY,m3q,m3u,MSUSY,MSUSY,MSUSY,MSUSY,MSUSY,MSUSY,MSUSY,MSUSY,MSUSY,MSUSY,MSUSY,mu,MSUSY,at,MSUSY,0,0,0,0,0,0,m1,m2,m3,0,0,0);
    TB=tb;
    MA0=ma;
    M3SU=m3u;
    M3SQ=m3q;
    MUE=mu;
    At=at;
    if(err==0) FHHiggsCorr(&err,MHiggs,&SAeff,UHiggs,ZHiggs);
    if(err==0) getOtherPars();
    if(err==0) calcCouplings();
    if(err==0) getSMPars();
    if(err==0) calcProd(14); 
    if(err==0) checkHB();
//    if(err==0) checkHiggsMass(Mhmin, Mhmax);
    if(err==0) checkLHC(); 
    if(err==0) calcProd(8);
    if(err==0&&(fail==0||fail>=lvl)) exprt();
    if(err==0&&fail==0) return 0;
    return 1;
  }
  
   

};


