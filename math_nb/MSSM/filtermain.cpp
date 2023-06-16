#include "filter.cpp"

int main(int argc, char* argv[]){
  double Mhmin, Mhmax, MAmin, MAmax;
  double datapoint[28];
  double m1, m2, mu, tb, at, ma, m3q, m3u, m3;
  int Atsgn=1, lvl=3, N;
  if(argc!=4){
    fprintf(stderr,"*****************************************************************************\n");
    fprintf(stderr,"Usage: ./filter inputfile output_level filename\n");
    fprintf(stderr,"argc=%i\n %s, %s\n",argc, argv[0], argv[1]);
    fprintf(stderr,"output_level =\n");
    fprintf(stderr,"  1 : all points printed.\n");
    fprintf(stderr,"  2 : points that pass HB LEP are printed.\n");
    fprintf(stderr,"  3 : points that pass HB LEP and HB Hadron colliders are printed.\n");
    fprintf(stderr,"  4 : points that pass HB and Mhmin<M<Mhmax are printed.\n");
    fprintf(stderr,"  5 : points that pass HB, Mhmin<M<Mhmax and cs/csSM(gg->h->ga,ga)>=0.8.\n");
    fprintf(stderr,"*****************************************************************************\n");
    exit(0);
  }
  else{   
    lvl=atoi(argv[2]);
    fprintf(stderr,"  output level = %d\n",lvl);
  }

  srand(time(NULL));
  
  //Initialize
  int err=0;
  FHSetFlags(&err, 4, 0, 0, 2, 0, 2, 1, 1, 0);
  if(err!=0)exit(1);
  /*int nHiggsneut=3, nHiggsplus=1;
  char whichanalyses[10];
  sprintf(whichanalyses,"%s ","LandH");
  initialize_higgsbounds_(nHiggsneut, nHiggsplus, whichanalyses);
  */
  
  //FHPars pars=FHPars(argv[1]);
  FILE *fp1, *fp2;
  char filename[100];
  sprintf(filename,"%s.math",argv[3]);
  fp1=fopen(filename,"a");
  sprintf(filename,"%s.dat",argv[3]);
  fp2=fopen(filename,"a");
  FHPars pars=FHPars(fp1,fp2);
  int iPass=0, fail;
//  ifstream indata;
//  indata.open(argv[1]);
  sprintf(filename,"%s.dat",argv[1]);
  FILE *fp3=fopen(filename,"r");
//  if(!indata){
  if(fp3==NULL){
     fprintf(stderr,"Fail to open input file");
     exit(0);
  }
  int i=0;
  char c;
  c=fgetc(fp3);
  while(c!=EOF){
//    for(int j=0;j<28;j++) indata >> datapoint[j];
     for(int j=0;j<27;j++) fscanf(fp3,"%le",&datapoint[j]);
     c=fgetc(fp3);
     c=fgetc(fp3);
     err=0;
     fail=0;
     i++;     
     m1=datapoint[4];
     m2=datapoint[5];
     mu=datapoint[6];
     tb=datapoint[7];
     at=datapoint[8];
     ma=datapoint[9];
     m3q=datapoint[10];
     m3u=datapoint[11];
     m3=datapoint[12];
     fprintf(stderr,"m1=%f, m2=%f, m3=%f, mu=%f, tb=%f, at=%f, m3q=%f, m3u=%f, ma=%f\n", m1, m2, m3, mu, tb, at, m3q, m3u, ma);
    fail=pars.next(m1,m2,m3,mu,tb,ma,m3q,m3u,at,lvl);
    if(i%50==0)
      fprintf(stderr,"%d/%d done, %d passed.\n",i,N,iPass);
    //fprintf(stderr,"i=%d\tfail=%d\n",i,fail);
    if(fail==0)iPass++;
  }
  fclose(fp1);
  fclose(fp2);
  //finish_higgsbounds_();
 

  fprintf(stderr,"\nFinished: %d passed!\n",iPass);
  return 0;
}
