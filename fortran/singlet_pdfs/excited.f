C   Sample-Monte-Carlo for pp-->New Color Resonance-->dijet
C	Zhen Liu UW-Madison zliu57@wisc.edu
C	Starting date : 9/1/2011
C	Modified base upon the following program:
C   SAMPLE-MONTE-CARLO FOR pp(bar) --> W --> l neutrino 
C   (1). IP=1,2 for PP & PPBAR;
C   (2). Include W channels ;
C   DATE :   3/3/94
c	include 'hanlib.inc'
       implicit none


	real xlam,ptmax,ymax,Qinc
	real R(4),DG(650),wgt,tchi2,tmean,tsigma
      real*8 meq,rtss,s,tau0,x1,x2,q,ss,sighat,as,tau, flux, aem
      real*8 u1,u2,d1,d2,c1,c2,s1,s2,ub1,ub2,db1,db2,sb1,sb2,cb1
      real*8 cb2,g1,g2,pi,conv,ctq6pdf,fc,b1,b2,bb1,bb2
	real*8 p1(4), p2(4), p01, p02, cost, sint, phi, pt, pcm(4)
	real*8 ptmin,rpmax,drmin,drmax,sparton,wid, twid
	real*8 rap1, rap2, drap, dr, lorentzcolor, jacobian
      integer iseed,iread,iwrite,many,itn,ir,idat,kn,i
	integer model,eventgen, initialid, spin, color
	logical bcorr,neutcorr
	character modelname 
	
c	EXTERNAL PLLCOMsol
c
	DATA CONV/0.389379323e9/, PI/3.14159265/

c

	common/iseed/iseed

	common /sctrl/Tmean,tsigma,iread,iwrite


	iseed=1000
	eventgen=0
	initialid=0

c	PRINT *, 'Enter many, Rtss'!, Excited quark mass:'
c	READ  *,  MANY, RTSS!
	MANY=100000
	RTSS=7000
       
	Meq = 3000
	ptmin = 200
	rpmax = 5
	drmin = 0.4
	drmax = 100

	write(*,*)"Generate events? 1. Yes; Other. No"
	read(*,*)eventgen
c
	open(20,file="out/dst332.dat")
c      write(20,*) 'Event File for Colored Resonance'//
c     .        'quark production'

	write(*,432) 'Events=',MANY
	write(*,433) 'Hadronic Center of Mass Energy=',RTSS
      write(*,433) 'Mass=', Meq
	write(*,433) 'Pt_min=',ptmin
	write(*,433) 'mu_F=mu_R=',Meq
	write(*,433) 'jet rapidity_min=', rpmax
	write(*,433) 'jets distances_min=', drmin
	write(*,433) 'jets distances_max=', drmax

432	format (A,I10)
433	format (A,F5.0)
434	format (I7,I3,I2,I2,I2,I2,E20.11,E20.11,E20.11,E20.11,
     !E20.11,F3.0,F4.0)
435	format (I2,	I3,	E17.7,	E17.7,	E17.7,	E17.7)

	ITN   = 10
      Ir = 1
      Idat = 0
      KN   = 0                ! ITERATION COUNTER

      TMEAN = 0.0
      TSIGMA = 0.0
      call SetCtq6(4)
	



	wgt = 0.0
      ss = rtss**2.
	tau0= 4*ptmin**2/ss


      call alphasrun(as,meq,4)
	aem=0.

	
	wid=twid/0.08716058*as/2000*meq


	write(*,*) 'itn=',itn

  100 CALL SAMPLE(R,WGT,4,MANY,ITN,KN,40,DG)

	tau=(1.0-tau0)*R(1)+tau0
	x2=(1-tau)*R(2)+tau
	x1=tau/x2

	p01=rtss/2.*x1
	p02=rtss/2.*tau/x1
	pcm(1)=0.
	pcm(2)=0.
	pcm(3)=p01-p02
	pcm(4)=p01+p02

	sparton=tau*ss

	cost=1.-2.*R(3)
	sint=dsqrt(1-cost**2)
	phi=2*pi*R(4)

	p1(4)=dsqrt(sparton)/2.
	p1(3)=p1(4)*cost
	pt=p1(4)*sint
	p1(2)=pt*dcos(phi)
	p1(1)=pt*dsin(phi)

	Do i=1, 3
	p2(i)=-p1(i)
	enddo

	p2(4)=p1(4)

	call	boostlab(pcm,p1)
	call	boostlab(pcm,p2)

	if(pt .lt. ptmin) then
	wgt=0.0
	goto 200
	endif

	call rapidity(p1, rap1)
	call rapidity(p2, rap2)

	if((rap1 .gt. rpmax) .or. (rap2 .gt. rpmax)) then
	wgt=0.0
	goto 200
	endif

	dr=dsqrt((rap1-rap2)**2+pi**2)

	if((dr .gt. drmax) .or. (dr .lt. drmin)) then
	wgt=0.0
	goto 200
	endif

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

	select case(initialid)
		case(1)
		FC = u1*u2+ub1*ub2
		case(2)
		FC = d1*d2+db1*db2
		case(3)
		FC = ub1*db2+ub2*db1+u1*d2+u2*d1
		case(4)
		FC = u1*g2+u2*g1+ub1*g2+ub2*g1
		case(5)
		FC = d1*g2+d2*g1+db1*g2+db2*g1
		case(6)
		FC = g1*g2
		case(7)
		FC = u1*ub2+u2*ub1
		case(8)
		FC = d1*db2+d2*db1
		case(9)
		FC = u1*db2+u2*db1+d2*ub1+d1*ub2
	end select

	sighat=sparton/((sparton-meq**2)**2+wid**2*meq**2)/2


	select case(model)
	case(1)
		!Diquark Scalar Triplet
		lorentzcolor=1./24.
	case(2)
		!Diquark Vector Sextet
		lorentzcolor=1./3.*((1-cost)/2.)**2
	case(3)
		!Excited Quark Triplet, Lambda=2M		
		lorentzcolor=1./9.*(1+cost)/2/8*sparton**2/meq**4
	case(4)
		!Excited Quark Triplet, spin 3/2, lambda=2M
		lorentzcolor=(sparton**2*(1+cost)**3
     !+meq**2*sparton*(1-cost)**3)/Meq**4/72.
	case(5)
	!Scalar Octet GG
		lorentzcolor=2.*5**2/3**2*sparton**2/Meq**4/2
	case(6)
	!Tensor Octet GG
		lorentzcolor=5.**2/3.**2/4./8.*Sparton**2*(
     !1+6*cost**2+cost**4+4.*(1)**4/9./meq**8*
     !(meq**2-sparton)**2*(2*meq**2+sparton)**2)/meq**4/2
	case(7)
	!Neutral Vector Octet QQb
		lorentzcolor=1/18.*(1+cost**2)
	end select

	Jacobian=1./x2*(1-tau)*(1-tau0)*conv*as**2*(4*Pi)**2/8/Pi
	
	wgt=sighat*FC*wgt*jacobian*lorentzcolor

c	Out put of Events. mimicing lhe format, but lots of null numbers like PDG#, spin, color
	If(eventgen .eq. 1) then
	write(20,*)'<event>'
	write(20,435)5,	1, wgt, meq, as, aem 
	write(20,434)9000001,	-1,	0,	0,	0,	0,	0.,	0.,	p01, p01,	0.	
     !,	0.,	-1.
	write(20,434)9000001,	-1,	0,	0,	0,	0,	0.,	0., -p02, p02,	0.
     !,	0.,	-1.
	write(20,434)9000006,	2,	0,	0,	0,	0,	0.,	0.,	
     !p01-p02, p01+p02,	meq,	0.,	-1.
	write(20,434)9000001,	1,	3,	0,	0,	0, p1(1), p1(2),	
     !p1(3), p1(4),	0.,	0.,	1.
	write(20,434)9000001,	1,	3,	0,	0,	0, p2(1), p2(2),	
     !p2(3), p2(4),	0.,	0.,	1.	
	write(20,*)'</event>'
	endif



 200   	IF (KN.LT.MANY) GOTO 100
	write(*,*) 'tmean=',tmean
       	write(20,70) meq,tmean
	write(20,*)'</LesHouchesEvents>'
 70    	format(2(F12.5,2X))
	
       	close(20)




	STOP
	END
ccccccccccccccccccccccc       
	subroutine alphasrun(alphas,mu,pdf)


	real*8 mu,lam,beta0,beta1,beta2,nf,pi
	real*8 dlog,lnmu,alphas,dabs,c1,c2,c3,c4
	integer trip, pdf

	Pi = 3.14159265d0
	nf = 5.d0
        if(pdf.eq.4) then
                lam = 0.165d0
        else
        	lam = 0.226d0
        endif


	lnmu = dlog(mu**2.d0/lam**2.d0)
	beta0 = 11.d0-2.d0/3.d0*nf
        if(pdf.eq.4) then
                beta1 = 0.d0
        else
                beta1 = 51.d0-19.d0/3.d0*nf
        endif
	beta2 = 2857.d0 -5033.d0/9.d0*nf+325.d0/27.d0*nf**2.d0
c	beta2 = 0.d0
	c1 = -beta0/(2.d0*Pi)
	c2 = -beta1/(8.d0*Pi**2.d0)

	c3 = -2.d0/c1
	c4 = -2.d0/c1**2.d0
	alphas = 4.d0*Pi/beta0/lnmu*(
     .	  1.d0-2.d0*beta1/beta0**2.d0*dlog(lnmu)/lnmu)


	return
	end

cccccccccccccccccccccccccc
c	This subroutine is to boost 4 moumentum p from p_ref's rest frame to p_ref's ref(lab) frame.
c	For simplicity, assume p_ref is going along z axis.
	subroutine boostlab(p_ref, p)
	real*8 p_ref(4), p(4), beta_ref, gamma_ref, dsqrt, p_tmp(4)
	beta_ref=-p_ref(3)/p_ref(4)
	gamma_ref=1./dsqrt(1-beta_ref**2)
	Do i=1, 4
	p_tmp(i)=p(i)
	enddo

	p(3)=gamma_ref*(p_tmp(3)-beta_ref*p_tmp(4))
	p(4)=gamma_ref*(p_tmp(4)-beta_ref*p_tmp(3))

	return
	end
	
cccccccccccccccccccccccccccc
	subroutine rapidity(p, rap)
	real*8 p(4), rap, dlog, dabs
	rap=1./2.*dabs(dlog((p(3)+p(4))/(p(4)-p(3))))

	return
	end
