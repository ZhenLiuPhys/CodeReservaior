C
C   SAMPLE-MONTE-CARLO FOR pp(bar) --> W --> l neutrino 
C   (1). IP=1,2 for PP & PPBAR;
C   (2). Include W channels ;
C   DATE :   3/3/94
c	include 'hanlib.inc'
       implicit none


	real xlam,ptmax,ymax,Qinc!,Pwknown(4),pwprec(4)
	real R(1),DG(2,40,1),wgt,tchi2,tmean,tsigma
       real*8 meq,rtss,s,tau0,x1,x2,q,ss,sighat,as, y
       real*8 u1,u2,d1,d2,c1,c2,s1,s2,ub1,ub2,db1,db2,sb1,sb2,cb1
       real*8 cb2,g1,g2,pi,conv,ctq6pdf,fc,b1,b2,bb1,bb2
       integer iseed,iread,iwrite,many,itn,ir,idat,kn,i
	logical bcorr,neutcorr
	real*8 fac1, fac2 
	
c	EXTERNAL PLLCOMsol
c
	DATA CONV/0.389379323e9/, PI/3.14159265/

c

	common/iseed/iseed

	common /sctrl/Tmean,tsigma,iread,iwrite



	ITN   = 10
c	PRINT *, 'Enter IReactn(1=pp;2=ppbar) & Idat(no=0):'
c	READ  *,  Ir, Idat 
        Ir = 1
        Idat = 0
        MANY = 1000





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
       call SetCtq6(4)


      open(20,file="y_dd.dat")
      write(20,*) '#zp mass'//
     .        '#unweighted false perentage'
c      write(20,*) '#mass    x-sect'
	meq = 3000
      do i=0,146
      y = 0+0.01*i

	wgt = 1.
c	meq = 3000.
	rtss = 13000.
	tau0 = meq**2./rtss**2.

c  100   CALL SAMPLE(R,WGT,1,MANY,ITN,KN,40,DG)

c       x1 = R(1)

c       write(*,*) x1

       x1 = dsqrt(exp(2.*y)*tau0)

	X2=Tau0/X1

	S=meq**2.

       Q = meq

	    U1 = Ctq6Pdf(1,x1,meq)
		D1 = Ctq6Pdf(2,x1,meq)
		S1 = Ctq6Pdf(3,x1,meq)
		C1 = Ctq6Pdf(4,x1,meq)
          b1 = ctq6pdf(5,x1,meq)
		UB1 = Ctq6Pdf(-1,x1,meq)
		DB1 = Ctq6Pdf(-2,x1,meq)
		SB1 = Ctq6Pdf(-3,x1,meq)
		CB1 = Ctq6Pdf(-4,x1,meq)
          bb1 = ctq6pdf(-5,x1,meq)
          g1 = ctq6pdf(0,x1,meq)

	    U2 = Ctq6Pdf(1,x2,meq)
		D2 = Ctq6Pdf(2,x2,meq)
		S2 = Ctq6Pdf(3,x2,meq)
		C2 = Ctq6Pdf(4,x2,meq)
          b2 = ctq6pdf(5,x2,meq)
		UB2 = Ctq6Pdf(-1,x2,meq)
		DB2 = Ctq6Pdf(-2,x2,meq)
		SB2 = Ctq6Pdf(-3,x2,meq)
		CB2 = Ctq6Pdf(-4,x2,meq)
          bb2 = ctq6pdf(-5,x2,meq)
          g2 = ctq6pdf(0,x2,meq)

      if(x1 .GT. x2) then
c		fac1 = u1*ub2
c	.	+c1*cb2
c		fac2 = u2*ub1
c	.	+c2*cb1
		fac1 = d1*db2
c	.	+s1*sb2
		fac2 = d2*db1
c	.	+s2*sb1
c		fac1 = s1*sb2+c1*cb2
c		fac2 = s2*sb1+c2*cb1
	else
c		fac1 = u2*ub1
c	.	+c2*cb1
c		fac2 = u1*ub2
c	.	+c1*cb2
		fac1 = d2*db1
c	.	+s2*sb1
		fac2 = d1*db2
c	.	+s1*sb2
c		fac2 = s1*sb2+c1*cb2
c		fac1 = s2*sb1+c2*cb1
	endif

	wgt = wgt * (fac2)/(fac1+fac2)


c 200   	IF (KN.LT.MANY) GOTO 100
c       write(*,*) 'tmean=',tmean

       write(20,70) y,wgt
 70    format(2(F20.5,2X))
       enddo
       close(20)



	STOP
	END

