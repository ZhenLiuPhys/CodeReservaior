C
C   SAMPLE-MONTE-CARLO FOR pp(bar) --> W --> l neutrino 
C   (1). IP=1,2 for PP & PPBAR;
C   (2). Include W channels ;
C   DATE :   3/3/94
c	include 'hanlib.inc'
       implicit none


	real xlam,ptmax,ymax,Qinc!,Pwknown(4),pwprec(4)
	real R(1),DG(2,40,1),wgt,tchi2,tmean,tsigma
       real*8 meq,rtss,s,tau0,x1,x2,q,ss,sighat,as
       real*8 u1,u2,d1,d2,c1,c2,s1,s2,ub1,ub2,db1,db2,sb1,sb2,cb1
       real*8 cb2,g1,g2,pi,conv,ctq6pdf,fc,b1,b2,bb1,bb2
       integer iseed,iread,iwrite,many,itn,ir,idat,kn,i
	logical bcorr,neutcorr 
	
c	EXTERNAL PLLCOMsol
c
	DATA CONV/0.389379323e9/, PI/3.14159265/

c

	common/iseed/iseed

	common /sctrl/Tmean,tsigma,iread,iwrite



	ITN   = 1
c	PRINT *, 'Enter IReactn(1=pp;2=ppbar) & Idat(no=0):'
c	READ  *,  Ir, Idat 
        Ir = 1
        Idat = 0
        MANY = 10





c       INITIALIZE SAMPLE VARIABLES (SOME TO MIMIC PHENO)
c        NDIM = 4              ! DIMENSION OF HYPERCUBE
c        ITN  = 5			! NUMBER OF BINNED ITERATIONS
        KN   = 0                ! ITERATION COUNTER

c        VOL = 0.0
        TMEAN = 0.0
        TSIGMA = 0.0
c        CHI2 = 0D0
c        XNG = 0d0
c        MEAN = 0.0
c        SIGMA = 0.0
c       call SetCtq6(4)


	wgt = 0.0

  100   CALL SAMPLE(R,WGT,1,MANY,ITN,KN,5,DG)

       x1 = R(1)

       write(*,*) x1

       wgt = wgt*x1





 200   	IF (KN.LT.MANY) GOTO 100
       write(*,*) 'tmean=',tmean

       close(20)




	STOP
	END

