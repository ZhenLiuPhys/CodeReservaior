c
c This is a complete collection of supporting codes, 
c as a library FORTRAN source file
c 
c It can be compiled/linked as hanlib.f,
c OR for convenience, compile/link the individual
c parts kinematics.f; spinors.f; integration.f; 
c Hbook.f or theory.f
c
c compiled for FNAL    on 10/20/1992 
c modified for UCDHTF  on 01/21/1994 
c modified for egret   on 04/08/1998 <--first time for unix
c split to subcodes    on 05/08/1998
c For OSX system       on 08/01/2006 <-- based on hanlib_linux.f on HEP
c     eliminate old PDF's and add in smearing code
c modified to include EM_invt etc/
C
c [ in metric (1,2,3,4) = (Px,y,z,E) ]
c
c It's a collection of the following parts:
c
c   (1). kinematics.f [Kinematics Related codes]
c       Angl;  Anglphi;  Boost;  CDot; CosA; DAngl; DBoost; DDot; 
c        Dec2; Dec3; DeltR; DLab; Dot; Dot2; Dot3; DDot3; 
c        Epsilon; CEpsilon; Lab;  Rapid; DRapid; XLam
c       EM_invt, E_tran, P_tran
c
c   (2). Spinors.f [Helicity amplitude subcodes]
c        a: Two-spinor Helicity Method (massless only)
c           Chi; 
c           Eps; Eps_c[these two are more general than spinor2]; 
c          Gama3; CGama3; CGama3_ns; Gama4; CGama4; 
c          Sig4; String9; CString9; CSleft9; CSright9
c        b: Four-spinor for massive (or massless) fermions 
c           Helpack8, 16 = epsn; u_spinor; v_spinor; slash; slash_c;
c           m0,1,2,3..9; matrix2,3; l_mult1,2; r_mult1,2; ab_gamma5;
c           current; current_g; generate_qcd; gamma3
c
c   (3). integration.f [Integration Subcodes]
c        Gaussx; ; HTuple; RN;
c
c   (4). hbook.f [HBook Related]
c        Hbook1; Hcurvesc; Hfill; hbook2
c
c   (5). theory.f [Theoretical collections, including math codes]
c        Gamma-function; alphas; alpha_s;
c
c   (6). Experiments related codes:
c        smearing in E and in pT
c
c -------------------------------------------------------------------
c   (1). kinematics.f [Kinematics Related codes] 
c -------------------------------------------------------------------
      FUNCTION ANGL(A,B,ID)
C       THIS FUNCTION RETURNS THE ANGLE BETWEEN TWO VECTORS
C       FOR ID = 1, ANGLE IN RADIAN ;
C       FOR ID NOT TO BE 1 , ANGLE IN DEGREE.
C      BY T. HAN, 11/15/1987
C
      REAL A(4), B(4)
      DATA PI/3.14159265/
      PM(PX,PY,PZ) = SQRT(PX**2 + PY**2 + PZ**2)
C
      DAB = DOT3(A,B)
      IF (DAB.EQ.0.0E0) DAB = 1.0E-30
      IF (DAB.EQ.1.0E0) DAB = 0.9999999
      AM  = AMAX1(1.0E-30,PM(A(1),A(2),A(3)))
      BM  = AMAX1(1.0E-30,PM(B(1),B(2),B(3)))
      COSA = DAB/(AM*BM)
      COSA = AMAX1(-1.0,COSA)
      COSA = AMIN1(1.0,COSA)
      ANGL= ACOS(COSA)
      IF (ID.EQ.1) RETURN
      ANGL= 180./PI*ANGL
      RETURN
      END
c -------------------------------------------------------------------
      FUNCTION ANGLPHI(A,B,ID)
C       THIS FUNCTION RETURNS THE AZIMUTHAL ANGLE BETWEEN TWO VECTORS
C      OR, THE ANGLE IN A TRANSVERSE PLAN
C ***   0 < ANGLE OF PHI < PI   ***
C       FOR ID = 1, ANGLE IN RADIAN ;
C       FOR ID NOT TO BE 1 , ANGLE IN DEGREE.
C      BY T. HAN, 1/4/1988
C
      REAL A(4), B(4)
      DATA PI/3.14159265/
      PM(PX,PY) = SQRT(PX**2 + PY**2)
C
      DAB  = DOT2(A,B)
C      IF (DAB.EQ.0.0E0) DAB = 1.0E-30
C      IF (DAB.EQ.1.0E0) DAB = 0.9999999
      AM   = AMAX1(1.0E-30,PM(A(1),A(2)))
      BM   = AMAX1(1.0E-30,PM(B(1),B(2)))
      COSA = DAB/(AM*BM)
      COSA = AMAX1(-1.0,COSA)
      COSA = AMIN1(1.0,COSA)
      ANGLPHI  = ACOS(COSA)
      IF (ID.EQ.1) RETURN
      ANGLPHI  = 180./PI*ANGLPHI
      RETURN
      END
c -------------------------------------------------------------------
      SUBROUTINE BOOST(P,MP,Q)
C ***
C      THIS SUBROUTINE CALCULATES BOOSTED 4-MOMENTA 
C      Q TO FRAME OF PARTICLE P
C      EQ=Q(4)
C      PQ=Q(I) ; I=1,3
C      MP = MASS OF PARTICLE P            -H.BAER-
C ***
      REAL P(4),Q(4),G,VQ,HQ,MP
      G=P(4)/MP
      VQ=(Q(1)*P(1)+Q(2)*P(2)+Q(3)*P(3))/P(4)
      HQ=VQ*G*G/(1.+G)/P(4)+Q(4)/MP
      Q(4)=G*(Q(4)+VQ)
      DO I=1,3
            Q(I)=Q(I)+HQ*P(I)
      END DO
      RETURN
      END
c -------------------------------------------------------------------
      FUNCTION CDOT(A,B)
C      THIS FUNCTION IS THE COMPLEX DOT PRODUCT OF 
C      TWO COMPLEX-4-VECTORS A AND B.
      COMPLEX A(4), B(4), CDOT
      CDOT = A(4)*B(4) - A(1)*B(1) - A(2)*B(2) - A(3)*B(3)
      RETURN
      END
c -------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DANGL(A,B,ID)
C       THIS FUNCTION RETURNS THE ANGLE BETWEEN TWO VECTORS
C       FOR ID = 1, ANGLE IN RADIAN ;
C       FOR ID NOT TO BE 1 , ANGLE IN DEGREE.
C      BY T. HAN, 11/15/1987
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(4), B(4)
      DATA PI/3.14159265/
      PM(PX,PY,PZ) = DSQRT(PX**2 + PY**2 + PZ**2)
C
      DAB = DDOT3(A,B)
      IF (DAB.EQ.0.0D0) DAB = 1.0D-40
      IF (DAB.EQ.1.0D0) DAB = 0.9999999
      AM  = DMAX1(1.0D-40,PM(A(1),A(2),A(3)))
      BM  = DMAX1(1.0D-40,PM(B(1),B(2),B(3)))
      COSA = DAB/(AM*BM)
      COSA = DMAX1(-1.0D0,COSA)
      COSA = DMIN1(1.0D0,COSA)
      DANGL= DACOS(COSA)
      IF (ID.EQ.1) RETURN
      DANGL= 180./PI*DANGL
      RETURN
      END
c -------------------------------------------------------------------
      SUBROUTINE DBOOST(P,MP,Q)
C ***
C      THIS SUBROUTINE CALCULATES BOOSTED 4-MOMENTA IN DOUBLE PRECISION
C      Q TO FRAME OF PARTICLE P
C      EQ=Q(4)
C      PQ=Q(I) ; I=1,3
C      MP = MASS OF PARTICLE P            -H.BAER-
C ***
      REAL*8 P(4),Q(4),G,VQ,HQ,MP
      G=P(4)/MP
      VQ=(Q(1)*P(1)+Q(2)*P(2)+Q(3)*P(3))/P(4)
      HQ=VQ*G*G/(1.+G)/P(4)+Q(4)/MP
      Q(4)=G*(Q(4)+VQ)
      DO I=1,3
            Q(I)=Q(I)+HQ*P(I)
      END DO
      RETURN
      END
c -------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DDOT(A,B)
C      THIS FUNCTION IS THE DOT PRODUCT OF 4-VECTORS A AND B.
      REAL*8 A(4), B(4)
      DDOT = A(4)*B(4) - A(1)*B(1) - A(2)*B(2) - A(3)*B(3)
      RETURN
      END
c -------------------------------------------------------------------
C ***
      SUBROUTINE DEC2(W,EM,EM1,EM2,P1,P2)
C
C      KINEMATIC GENERATOR OF TWO-BODY DECAY EVENTS 
C         WITH THE ABSOLUTE PHASE SPACE WEIGHT W
C
C      IMPLICIT REAL*8(A-H,O-Y)
      REAL P1(4), P2(4)
        REAL*8      RN
      DATA PI/3.14159265/, ISEED/9025/
c
       PMAX  = SQRT((EM**2-(EM1-EM2)**2)*(EM**2-(EM1+EM2)**2))/(2.*EM)
C      
C      RANDOMIZE THE OREITATION OF THE MOMENTA IN M-REST FRAME
C
      CT    = 1.-2.*RN(ISEED)
      ST    = SQRT(1.-CT*CT)
      PHI   = 2.*PI*RN(ISEED)
C   FOR THE DECAY PRODUCT 1 :
      P1(4) = (EM**2+EM1**2-EM2**2)/(2.*EM)
      P1(1) = PMAX*ST*COS(PHI)
      P1(2) = PMAX*ST*SIN(PHI)
      P1(3) = PMAX*CT
C   FOR THE DECAY PRODUCT 2 : 
      P2(4) = EM - P1(4)
      P2(1) = -P1(1)
      P2(2) = -P1(2)
      P2(3) = -P1(3)
C
C SET THE PHASE SPACE WEIGHT
      W = PI*PMAX/EM
C
      RETURN
      END
c -------------------------------------------------------------------
C***
      SUBROUTINE DEC3(W0,W3,EM0,EM1,EM2,EM3,P1,P2,P3)
C
C   KINEMATIC GENERATOR OF THREE-BODY DECAY EVENTS 
C      IN EM-REST FRAME WITH WEIGHT ---- 9/30/1988
C***
C   Special Remark: The normalized weight W0 is only for m1=m2=m3=0;
C      or m2=m2=0 with 0<m1<m0
C      All massive case is very complicated for the weight.
C      The absolute phase space weight is in W3
C***
C   INPUT  MASSES : EM0; EM1,EM2,EM3;
C   OUTPUT THREE  4-MOMENTA P1,P2,P3 IN EM-REST FRAME;
C         and NORMALIZED & ABSOLUTE WEIGHTS W0, W3
C
      REAL*4 P1(4), P2(4), P3(4), P23(4)
        REAL*8 RN
      DATA   PI/3.14159265/, ISEED/13451/
C   EXTERNAL FUNCTION :
      ALAM(X,Y,Z) = (X - Y - Z)**2 - 4.*Y*Z 
c
      EM02  = EM0**2
      EM12  = EM1**2
      EM22  = EM2**2
      EM32  = EM3**2
       PMAX1 = SQRT(ALAM(EM02,EM12,(EM2+EM3)**2))/(2.*EM0)
      PM1   = PMAX1*RN(ISEED)
C      
C   RANDOMIZE THE OREITATION OF THE P1 IN M-REST FRAME
C
      CT    = 1.-2.*RN(ISEED)
      ST    = SQRT(1.-CT*CT)
      PHI   = 2.*PI*RN(ISEED)
C   FOR THE DECAY PRODUCT 1 :
      P1(1) = PM1*ST*COS(PHI)
      P1(2) = PM1*ST*SIN(PHI)
      P1(3) = PM1*CT
      P1(4) = SQRT(PM1**2 + EM12)
C   BY THE WAY, THE SYSTEM IN EM-REST FRAME:
      DO I = 1, 3
         P23(I) = -P1(I)
      END DO
           P23(4) =  EM0 - P1(4)
      EM232 = EM02 + EM12 - 2.*EM0*P1(4)
      IF (EM232.LE.0.0) THEN
            PRINT *, ' M(2,3)**2 = ', EM232
            EM232 = 1.0E-16
      END IF
      EM23  = SQRT(EM232)
      IF (EM23.Eq.0.0) THEN
            PRINT *, ' M(2,3) = ', EM23
            EM23 = 1.0E-8
      END IF
C   GET THE P2, P3 IN 23-REST FRAME :
       PM2   = SQRT(ALAM(EM232,EM22,EM32))/(2.*EM23)
C      
C   RANDOMIZE THE OREITATION OF THE P2 IN 23-REST FRAME
C
      CT    = 1.-2.*RN(ISEED)
      ST    = SQRT(1.-CT*CT)
      PHI   = 2.*PI*RN(ISEED)
C   FOR THE DECAY PRODUCT 2 :
      P2(1) = PM2*ST*COS(PHI)
      P2(2) = PM2*ST*SIN(PHI)
      P2(3) = PM2*CT
      P2(4) = SQRT(PM2**2 + EM22)
C   FOR THE DECAY PRODUCT 3 :
      DO I = 1, 3
         P3(I) = -P2(I)
      END DO
           P3(4) =  EM23 - P2(4)
C   NOW BOOST P2, P3 BACK TO EM-REST FRAME :
      CALL BOOST(P23,EM23,P2)
      CALL BOOST(P23,EM23,P3)
C
C Normalized phase space weight:
C
      if (em1.eq.0.0.and.em2.eq.0.0.and.em3.eq.0.0) then
            W0 = 4.0*pm1/em0
      else if (em1.ne.0.0) then
            beta = 1. - em12/em02
            W0 = 4.0*pm1**2/em0/P1(4)*beta/
     1         (1. - (em12/em02)**2 - 4.*em12/em02*alog(em0/em1))
      else
            W0 = 1.0
      end if
C
C ABSOLUTE PHASE SPACE WEIGHT :
C
      W3 = 2.0*PI**2*PMAX1*PM1**2*PM2/P1(4)/EM23
C
      RETURN
      END
c -------------------------------------------------------------------
      REAL FUNCTION DELTR(A,B)
C
CALCULATE THE SEPERATION, DELTA_R, FOR THE TWO FOUR-VECTORS A & B
C
      REAL*4 A(4), B(4)
C
      DELTPHI = ANGLPHI(A,B,1)
      DELTY   = RAPID(A) - RAPID(B)
C
      DELTR   = SQRT(DELTY**2 + DELTPHI**2)
C
      RETURN
      END
c -------------------------------------------------------------------
c
      Function Epsilon(A,B,C,D)
c
c      This Function evaluates the Levi_Civita epsilon tensor
c      contracted with 4 four-vectors A,B,C,D:
c                eps^(mu,nu,rho,sig)*A_mu B_nu C_rho D_sig
c      Following VB & RP's book:
c              eps^(4,1,2,3) = -eps_(4,1,2,3) = +1
c              p_mu = (p4, -p1, -p2, -p3)
c      where 1,2,3,4 correspond to px,py,pz,E
cked by Han on 11/09/93
C
      Implicit None
c
      Real*4 A(4), B(4), C(4), D(4)
      Real*4 Epsilon, Term1, Term2, Term3, Term4
c
      Term4 =   A(4)*B(1)*C(2)*D(3) - A(4)*B(1)*C(3)*D(2)
     &       - A(4)*B(2)*C(1)*D(3) + A(4)*B(2)*C(3)*D(1)
     &       + A(4)*B(3)*C(1)*D(2) - A(4)*B(3)*C(2)*D(1)

      Term1 = - A(1)*B(2)*C(3)*D(4) + A(1)*B(2)*C(4)*D(3)
     &       + A(1)*B(3)*C(2)*D(4) - A(1)*B(3)*C(4)*D(2)
     &       - A(1)*B(4)*C(2)*D(3) + A(1)*B(4)*C(3)*D(2)

      Term2 =   A(2)*B(1)*C(3)*D(4) - A(2)*B(1)*C(4)*D(3)
     &       - A(2)*B(3)*C(1)*D(4) + A(2)*B(3)*C(4)*D(1)
     &       + A(2)*B(4)*C(1)*D(3) - A(2)*B(4)*C(3)*D(1)

      Term3 = - A(3)*B(1)*C(2)*D(4) + A(3)*B(1)*C(4)*D(2)
     &       + A(3)*B(2)*C(1)*D(4) - A(3)*B(2)*C(4)*D(1)
     &       - A(3)*B(4)*C(1)*D(2) + A(3)*B(4)*C(2)*D(1)
c
      Epsilon = Term1 + Term2 + Term3 + Term4
change the sign to give epsilon(SUPERscript mu,nu,rho,sig): 2/1/94
      Epsilon = - Epsilon 
c
      Return
      End
c -------------------------------------------------------------------
c
      Complex Function Cepsilon(A,B,C,D)
c
c      This Function evaluates the Levi_Civita epsilon tensor
c      contracted with 4 complex four-vectors A,B,C,D:
c                eps^(mu,nu,rho,sig)*A_mu B_nu C_rho D_sig
c      Following VB & RP's book:
c              eps^(4,1,2,3) = -eps_(4,1,2,3) = +1
c              p_mu = (p4, -p1, -p2, -p3)
c      where 1,2,3,4 correspond to px,py,pz,E
cked by Han on 11/09/93
C
      Implicit None
c
      Complex A(4),  B(4),  C(4),  D(4)
      Complex Term1, Term2, Term3, Term4
c
      Term4 =   A(4)*B(1)*C(2)*D(3) - A(4)*B(1)*C(3)*D(2)
     &       - A(4)*B(2)*C(1)*D(3) + A(4)*B(2)*C(3)*D(1)
     &       + A(4)*B(3)*C(1)*D(2) - A(4)*B(3)*C(2)*D(1)

      Term1 = - A(1)*B(2)*C(3)*D(4) + A(1)*B(2)*C(4)*D(3)
     &       + A(1)*B(3)*C(2)*D(4) - A(1)*B(3)*C(4)*D(2)
     &       - A(1)*B(4)*C(2)*D(3) + A(1)*B(4)*C(3)*D(2)

      Term2 =   A(2)*B(1)*C(3)*D(4) - A(2)*B(1)*C(4)*D(3)
     &       - A(2)*B(3)*C(1)*D(4) + A(2)*B(3)*C(4)*D(1)
     &       + A(2)*B(4)*C(1)*D(3) - A(2)*B(4)*C(3)*D(1)

      Term3 = - A(3)*B(1)*C(2)*D(4) + A(3)*B(1)*C(4)*D(2)
     &       + A(3)*B(2)*C(1)*D(4) - A(3)*B(2)*C(4)*D(1)
     &       - A(3)*B(4)*C(1)*D(2) + A(3)*B(4)*C(2)*D(1)
c
      Cepsilon = Term1 + Term2 + Term3 + Term4
change the sign to give epsilon(SUPERscript mu,nu,rho,sig): 2/1/94
      Cepsilon = - Cepsilon 
c
      Return
      End
c -------------------------------------------------------------------
      SUBROUTINE DLAB(P,X1,X2)
C      THIS SUBROUTINE DOES A LONGITUDINAL BOOST ON 4-VECTOR P 
C       IN DOUBLE PRECISION
C      FROM THE PARTON C.M. FRAME TO THE COLLIDER LAB FRAME.      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(4)
      BETA       = (X1 - X2)/(X1 + X2)
      GAMMA   = 1.D0/(DSQRT(1.D0 - BETA**2))
      E      = P(4)
      PZ      = P(3)
      P(4)      = GAMMA*(E  + BETA*PZ)
      P(3)      = GAMMA*(PZ + BETA*E )
      RETURN
      END
c -------------------------------------------------------------------
      FUNCTION DOT(A,B)
C      THIS FUNCTION IS THE DOT PRODUCT OF 4-VECTORS A AND B.
      REAL A(4), B(4)
      DOT = A(4)*B(4) - A(1)*B(1) - A(2)*B(2) - A(3)*B(3)
      RETURN
      END
c -------------------------------------------------------------------
      FUNCTION DOT2(A,B)
C      THIS FUNCTION RETURNS THE TRANSVERSE MONENTUM DOT PRODUCT.
      REAL A(4), B(4)
      DOT2      = A(1)*B(1) + A(2)*B(2) 
      RETURN
      END
c -------------------------------------------------------------------
      FUNCTION DOT3(A,B)
C      THIS FUNCTION RETURNS THE 3-VECTER DOT PRODUCT.
      REAL A(4), B(4)
      DOT3      = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      RETURN
      END
c -------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DDOT2(A,B)
C      THIS FUNCTION RETURNS THE transverse momentum DOT PRODUCT.
      REAL*8 A(4), B(4)
      DDOT2      = A(1)*B(1) + A(2)*B(2)
      RETURN
      END
c -------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DDOT3(A,B)
C      THIS FUNCTION RETURNS THE 3-VECTER DOT PRODUCT.
      REAL*8 A(4), B(4)
      DDOT3      = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      RETURN
      END
c -------------------------------------------------------------------
      SUBROUTINE LAB(P,X1,X2)
C      THIS SUBROUTINE DOES A LONGITUDINAL BOOST ON 4-VECTOR P 
C      FROM THE PARTON C.M. FRAME TO THE COLLIDER LAB FRAME.      
      REAL*8 BETA,GAMMA,E,PZ
      REAL*4 P(4),X1,X2
      BETA       = (X1 - X2)/(X1 + X2)
      GAMMA   = 1.D0/(DSQRT(1.D0 - BETA**2))
      E      = P(4)
      PZ      = P(3)
      P(4)      = GAMMA*(E  + BETA*PZ)
      P(3)      = GAMMA*(PZ + BETA*E )
      RETURN
      END
c -------------------------------------------------------------------
      REAL FUNCTION RAPID(A)
C
CALCULATE THE RAPIDITY OF REAL 4-MOMENTUM A(4)
C
      REAL*4 A(4)
      RAPID = 0.5*ALOG(AMAX1(1E-20,(A(4)+A(3)))/
     1                   AMAX1(1E-20,(A(4)-A(3))))
C
      RETURN
      END
c -------------------------------------------------------------------
      Double Precision FUNCTION DRAPID(A)
C
CALCULATE THE RAPIDITY OF REAL 4-MOMENTUM A(4)
C
      REAL*8 A(4)
      DRAPID = 0.5D0*DLOG(DMAX1(1D-20,(A(4)+A(3)))/
     1                      DMAX1(1D-20,(A(4)-A(3))))
C
      RETURN
      END
c -------------------------------------------------------------------
      FUNCTION XLAM(A,B,C)

C      THIS FUNCTION IS THE 2-BODY PHASE SPACE LAMDA FUNCTION.
c
      XLAM = A*A + B*B + C*C - 2.*(A*B + A*C + B*C)
c
      IF (XLAM .LT. 0.) XLAM=0.
c
      RETURN
      END

c -------------------------------------------------------------------
        REAL FUNCTION P_Tran(A)
C
CALCULATE THE transverse momentum OF a REAL 4-MOMENTUM A(4)
C
        REAL*4 A(4)
c
        SQ = DOT2(A,A)
        IF (SQ.GE.0.) THEN
           P_Tran = SQRT(SQ)
        ELSE
c         WRITE(*,*) 'PT SQ=' ,SQ
           P_Tran = 0.
        ENDIF
C
        RETURN
        END

c -------------------------------------------------------------------
        REAL FUNCTION E_Tran(A)
C
CALCULATE THE transverse energy OF a REAL 4-MOMENTUM A(4)
C
        REAL*4 A(4)
c
        ptsq = ( p_tran(A) )**2
        SQ   = DOT(A,A) + ptsq
        IF (SQ.GE.0.) THEN
           E_Tran = SQRT(SQ)
        ELSE
c         WRITE(*,*) 'PT SQ=' ,SQ
           E_Tran = 0.
        ENDIF
C
        RETURN
        END

c -------------------------------------------------------------------
        REAL FUNCTION EM_invt(A)
C
CALCULATE THE IVARIANT MASS OF REAL 4-MOMENTUM A(4)
C
        REAL*4 A(4)
c
        SQ  = DOT(A,A)
        IF (SQ.GE.0.) THEN
            EM_invt = SQRT(SQ)
        ELSE
c       WRITE(*,*) 'IVMASS SQ=' ,SQ
            EM_invt = 0.
        ENDIF
C
        RETURN
        END

c -------------------------------------------------------------------
c    (2). a: Massless Helicity Method [a la Hagiwara-Zeppenfeld]
c -------------------------------------------------------------------
C   EVALUATE THE TWO COMPONENT SPINORS : CHI = CHI(+,P) OR CHI(-,P)
C   AS FUNCTIONS OF A FOUR MOMENTUM P.  
C   SEE EQ.(3.22) AND EQ.(3.23) IN K.H. & D.Z.'S PAPER.
C
C      INPUT P & ALPH(+-1)
C       OUTPUT CHII AS A COMPLEX TWO DIMENSIONAL ARRAY
C   HELICITY INDEX ALPH > 0 FOR CHI+ ;  ELSE FOR CHI-
C
C   ALL FOUR VECTORS ARE ASSUMED TO HAVE TIME-COMPONENT AS FOURTH
C   ENTRY P(4).
C-----------------------------------------------------------------------
C
      SUBROUTINE CHI(P,ALPH,CHII)
      INTEGER   ALPH
      REAL      P(4)
      COMPLEX   CHII(2), CHIAU
      DATA EP/1E-30/
C
      PT2  = P(1)**2 + P(2)**2
      PABS = SQRT(PT2 + P(3)**2)
      PAPZ = PABS + ABS(P(3))
      IF (P(3).LT.0.) PAPZ = PT2/PAPZ
      IF(PAPZ.GT.EP*PABS) THEN
           DEN = 1./SQRT(2.*PABS*PAPZ)
           CHII(1) = DEN * PAPZ
           CHII(2) = DEN * CMPLX( P(1), P(2) )
      ELSE
           CHII(1) = 0.
           CHII(2) = 1.
      END IF
      IF (ALPH.GT.0)  RETURN
      CHIAU   =  CHII(1)
      CHII(1) = -CONJG(CHII(2))
      CHII(2) =  CHIAU
      RETURN
      END
c -------------------------------------------------------------------
C   MASSIVE BOSON POLARIZATION VECTORS IN RECTANGULAR BASIS
C   ARE CALCULATED HERE.
C
C   INPUT : P, EM & L --- FOUR-MOMENTUM, MASS & POLARIZATION INDEX
C   OUTPUT: FOUR-VECTOR EPS(MU) WITH THE POLARIZATION INDEX L(=1,2,3) 
C
      SUBROUTINE EPS(P,EM,L,EPSL)
      INTEGER  L
      REAL  P(4), EM, EPSL(4)
C
      PT2  = P(1)**2  + P(2)**2
      PT   = SQRT(PT2)
      PABS = SQRT(PT2 + P(3)**2)
      EPSL(4) = 0.0
      IF (L.EQ.1) THEN
        IF (PT.LE.1E-20*ABS(P(3))) THEN
          EPSL(1) = P(3)/ABS(P(3))
          EPSL(2) = 0.0
          EPSL(3) = 0.0
        ELSE
          FAC1 = 1.0/(PABS*PT)
          EPSL(1) = P(1)*P(3)*FAC1
          EPSL(2) = P(2)*P(3)*FAC1
          EPSL(3) = -PT2*FAC1
        END IF
      ELSE IF (L.EQ.2) THEN
        IF (PT.LE.1E-20*ABS(P(3))) THEN
          EPSL(1) = 0.0
          EPSL(2) = 1.0
          EPSL(3) = 0.0
        ELSE
          FAC2 = 1.0/PT
          EPSL(1) = -P(2)*FAC2
          EPSL(2) =  P(1)*FAC2
          EPSL(3) = 0.0
        END IF
      ELSE
      FAC = P(4)/(EM*PABS)
      DO I = 1, 3
        EPSL(I) = P(I)*FAC
      END DO
        EPSL(4) = PABS/EM
      END IF
      RETURN
      END
c -------------------------------------------------------------------
c
      SUBROUTINE EPS_C(K,KM,SIGMA,LAMBDA,CEPS)
c
C------------------------------------------------------------------------
C
C
C   INPUT:
C
C            K(4)        real*4 physical four momentum:
C                        [kx, ky, kz, E]
C            KM          real*4 mass of the boson: K.K = KM**2
C
C            sigma       integer sign factor for vector boson
C                        = +1 outgoing
C                        = -1 incoming
C                        [For outgoing bosons (sigma=+1) the complex
C                         conjugate polarisation vector is returned.]
C
C            lambda      integer: helicity index
C                        = -1,0,1 (with 0 = the longitudinal)
C
C   OUTPUT:
C
C            ceps          complex  array(4)
C
C   Tested Against HELVEC in VVSERVICE.for on 10/16/93 --- Han
C------------------------------------------------------------------------
C
C   arguments:
C
      Implicit none
C
      REAL  K(4), KM
      COMPLEX CEPS(4)
      INTEGER SIGMA, LAMBDA
C
C   local variables:
C
      REAL  A0, A1, A2, A3, B0, B1, B2, B3
      REAL  KT, KMAGNI, NORMAL
C
      CEPS(4) = 0.
      IF (ABS(LAMBDA) .EQ. 1) THEN
          KT = SQRT(2.D0*( K(1)**2 + K(2)**2 ))
          IF ( KT .GT. 0.0 ) THEN
             IF ( KM .EQ. 0.0 ) THEN
                  KMAGNI = K(4)
             ELSE
                  KMAGNI = SQRT( K(1)**2 + K(2)**2 + K(3)**2 )
             ENDIF
            NORMAL = -LAMBDA/(KMAGNI*KT)
            A1 =  K(1)*K(3)*NORMAL
            A2 =  K(2)*K(3)*NORMAL
            A3 = -KT**2*0.5D0*NORMAL
            B1 = -K(2)/KT
            B2 =  K(1)/KT
            CEPS(1) = CMPLX(A1, SIGMA*B1)
            CEPS(2) = CMPLX(A2, SIGMA*B2)
            CEPS(3) = A3
          ELSE
            NORMAL  = -LAMBDA/SQRT(2.)
            CEPS(1) = SIGN(1.0, K(3))*NORMAL
            CEPS(2) = CMPLX(0.0, -LAMBDA*NORMAL*sigma)
            CEPS(3) = 0.0
          ENDIF
      ELSE IF (LAMBDA .EQ. 0) THEN
            IF ( KM .EQ. 0) THEN
C
C                  Mass = 0 not allowed for longitudinal polarization
C
                  WRITE(*,*) 'Mass = 0 for Lambda = 0 in EPS_C'
                  CEPS(1) = 0.D0
                  CEPS(2) = 0.D0
                  CEPS(3) = 0.D0
                  RETURN
            END IF
            KMAGNI  = SQRT( K(1)**2+K(2)**2+K(3)**2 )
            NORMAL  = K(4)/( KM*KMAGNI )
          CEPS(4) = KMAGNI/KM      
            CEPS(1) = K(1)*NORMAL
            CEPS(2) = K(2)*NORMAL
            CEPS(3) = K(3)*NORMAL
      ELSE
C
C           Unrecognised value of Lambda
C
            WRITE (*,*) 'Invalid Lambda in EPS_C: Lambda = ',LAMBDA
                  CEPS(1) = 0.D0
                  CEPS(2) = 0.D0
                  CEPS(3) = 0.D0
            END IF
      RETURN
      END
c -------------------------------------------------------------------
CALCULATIONS FOR TRIPLE-VECTOR-BOSON COUPLING IN TERMS OF
C THE TWO MOMENTA AND THE POLARIZATIONS OF THE OUTGOING BOSONS
C
CONTRACTIONS WITH A PROPAGATOR IN UNITARY OR FEYNMANN GAUGES
C INPUT :
C      P1(4), P2(4) --- TWO MOMENTA OF OUTGOING BOSONS
C      E1(4), E2(4) --- TWO POLARIZATION VECTORS OF THEM
C      EMV2 --- SQUARED VECTOR BOSON MASS APPEARED IN THE UNITARY GAUGE
C      GI   --- GAUGE PARAMETER : 
C             GI=1: FEYNMANN GAUGE;
C             GI=0: UNITARY GAUGE.
C
C OUTPUT:
C      GAMA(DELT) --- GAMA(MU,NU,LAM)*E1(MU)*E2(NU)*
C                  (G(DELT,LAM)-PROP(DELT)*PROP(LAM)/EMV2)
C
C NOTE1: THE ORDER OF THE TWO INPUT MOMENTA OF OUTGOING BOSONS IS IMPORTANT!!
C      P1,P2 = P(W+),P(W-);   P1,P2 = P(Z),P(W+);   P1,P2 = P(W-), P(Z) .
C      AND ALL OF THE ABOVE MOMENTA ARE PHYSICAL.
C NOTE2: WE HAVE USED THE IDENTITY :
C            P1.E1 = P2.E2 = 0
C NOTATIONS FOLLOWING THE BOOK :
C            TRIPLE VERTEX = -i*g_v*GAMA(MU,NU,LAM)
C      UPDATE 4/3/88
      SUBROUTINE GAMA3(P1,P2,E1,E2,EMV2,GI,GAMA)
      REAL P1(4), P2(4), E1(4), E2(4), GAMA(4),
     1       P12(4), P21(4), PP(4)
C
      DO I= 1, 4
         P12(I) = 2.*P1(I) + P2(I)
         P21(I) = 2.*P2(I) + P1(I)
      END DO
      DEE   = DOT(E1,E2)
      DP12E2 = DOT(P12,E2)
      DP21E1 = DOT(P21,E1)
      DO I  = 1, 4
         GAMA(I) = (P2(I)-P1(I))*DEE+DP12E2*E1(I)-DP21E1*E2(I)
      END DO
      IF (GI.EQ.1.0) RETURN
      DO I = 1, 4
            PP(I)  = P1(I) + P2(I)
      END DO
      DPG  = DOT(PP,GAMA)
      DPP  = DOT(PP,PP)
      DO I = 1, 4
         GAMA(I) = GAMA(I) + (1.0-GI)*DPG/(GI*DPP-EMV2)*PP(I)
      END DO
      RETURN
      END
c -------------------------------------------------------------------
C CALCULATIONS FOR TRIPLE-VECTOR-BOSON COUPLING IN TERMS OF
C THE TWO MOMENTA AND THE POLARIZATIONS OF THE OUTGOING BOSONS
C
CONTRACTIONS WITH A PROPAGATOR IN UNITARY OR FEYNMANN GAUGES
C INPUT :
C      P1(4), P2(4) --- (real*4)TWO MOMENTA OF OUTGOING BOSONS
C      E1(4), E2(4) --- (complex)TWO POLARIZATION VECTORS OF THEM
C      EMV2 --- SQUARED VECTOR BOSON MASS APPEARED IN THE UNITARY GAUGE
C      GI   --- GAUGE PARAMETER : 
C             GI=1: FEYNMANN GAUGE;
C             GI=0: UNITARY GAUGE.
C
C OUTPUT:
C      (complex)GAMA(DELT) --- GAMA(MU,NU,LAM)*E1(MU)*E2(NU)*
C                  (G(DELT,LAM)-PROP(DELT)*PROP(LAM)/EMV2)
C
C NOTE1: THE ORDER OF THE TWO INPUT MOMENTA OF OUTGOING BOSONS IS IMPORTANT!!
C      P1,P2 = P(W+),P(W-);   P1,P2 = P(Z),P(W+);   P1,P2 = P(W-), P(Z) .
C      AND ALL OF THE ABOVE MOMENTA ARE PHYSICAL.
C NOTATIONS FOLLOWING THE BOOK :
C            TRIPLE VERTEX = -i*g_v*GAMA(MU,NU,LAM)
C      UPDATE 4/3/88
      SUBROUTINE CGAMA3(P1,P2,E1,E2,EMV2,GI,GAMA)
c
      implicit none
c
      integer i
      REAL*4  P1(4), P2(4), GI, EMV2
      complex cdot, dee, dp12e2, dp21e1, dpg, dpp
      complex E1(4),E2(4),CP1(4),CP2(4),P12(4),P21(4),PP(4),Gama(4)
C
      Do i=1,4
        CP1(I) = CMPLX(P1(I),0.)
        CP2(I) = CMPLX(P2(I),0.)
      END DO
      DO I= 1, 4
         P12(I) = 2.*CP1(I) + CP2(I)
         P21(I) = 2.*CP2(I) + CP1(I)
      END DO
      DEE    = CDOT(E1,E2)
      DP12E2 = CDOT(P12,E2)
      DP21E1 = CDOT(P21,E1)
      DO I   = 1, 4
         GAMA(I) = (P2(I)-P1(I))*DEE + DP12E2*E1(I) - DP21E1*E2(I)
      END DO
      IF (GI.EQ.1.0) RETURN      ! Feynman Gauge here
      DO I = 1, 4
            PP(I)  = CP1(I) + CP2(I)
      END DO
      DPG  = CDOT(PP,GAMA)
      DPP  = CDOT(PP,PP)
      DO I = 1, 4
         GAMA(I) = GAMA(I) + (1.0-GI)*DPG/(GI*DPP-EMV2)*PP(I)
      END DO
c
      RETURN
      END
c -------------------------------------------------------------------
C
      Subroutine cgama3_ns(k2, k3, eps2, eps3, mv2, kappa, lambda,
     1                 lambdap, gauge, gama)
c
      Implicit none
C Input variables:
      Real*4 k2(4), k3(4), mv2, gauge, kappa, lambda, lambdap
C
      Complex eps2(4), eps3(4), gama(4)
C 
C FUNCTIONAL DESCRIPTION:      
C 
C    CGAMA3_NS calculates the three gauge boson vertex including non-standard
C      model couplings with CP-conservation.  
C      It dots in the inital state boson propagator with an arbitary gauge
C      as well as the 2 fermion currents from final state bosons so as 
C       to return a complex 4-vector.  
C
C REFERENCES:
C
C      See Eq.2.1 of Hagiwara et al. in Nucl Phys B282 (1987) 253-307;
C      Or See Eq.2b of M.Hayashi et al. in Lett Nuov Cimemto 44 (1986) 455.
C 
C      No relation like eps.p=0 is used, so that this code is valid for
C         any real or virtual V's
C 
C DUMMY ARGUMENTS:
C 
C    k1, k2   : real 4-momenta of the outgoing vector bosons
C  eps1, eps2 : complex fermion currents from the outgoing boson decays
C     mv2     : mass of the incoming gauge boson.  Used in the propagator
C    kappa    : anomolous magnetic moment of the W. Is 1 in Standard Model
C   lambda    : anomolous coupling.  Is 0 in the S.M. 
C            This lambda should be divided by MW**2 BEFORE it is passed in.
C   lambdap   : anomolous coupling.  Is 0 in the S.M. 
C            This lambdap should be divided by MW**2 BEFORE it is passed in.
C   gauge     : gives the gauge we are in.  0 for Unitary, 1 for Feynman
C    gama     : returns a complex 4 vector.
C 
C IMPLICIT INPUTS:
C 
C      The order of the 2 input momenta is important.  It can only be
C      k2,k3 = k(W+)k(W-), k(Z)k(W+) or k(W-)k(Z)
C
C Local variables:
C
      Real*4  kp(4), km(4), k2k2, k3k3, k2k3, denom, dot
      Complex ckp(4),ckm(4),ck2(4),ck3(4),k2e2,k3e3,k2e3,k3e2, e2e3
      complex cdot, cepsilon
      Integer i
C
      k2k2 = dot(k2, k2)
      k3k3 = dot(k3, k3)
      k2k3 = dot(k2, k3)
C
      do i = 1, 4
         kp(i)  = k2(i) + k3(i)
         km(i)  = k3(i) - k2(i)
         ckp(i) = cmplx(kp(i))
         ckm(i) = cmplx(km(i))
         ck2(i) = cmplx(k2(i))
         ck3(i) = cmplx(k3(i))
      enddo
C
      k2e2 = cdot(ck2, eps2)
      k3e3 = cdot(ck3, eps3)
      k2e3 = cdot(ck2, eps3)
      k3e2 = cdot(ck3, eps2)
      e2e3 = cdot(eps2, eps3)
C
      do  i = 1, 4
          gama(i) =  ckm(i) * e2e3
     1          + ((kappa + 1.) *  k2e3 + kappa * k3e3) * eps2(i)
     2          - ((kappa + 1.) *  k3e2 + kappa * k2e2) * eps3(i)
c
     3           + lambda *(e2e3 * (k2k2 + k2k3) * ck3(i)
     4                 - e2e3 * (k3k3 + k2k3) * ck2(i)
     5                 +(k3k3 *  k2e3 - k2k3  * k3e3) * eps2(i)
     6                 -(k2k2 *  k3e2 - k2k3  * k2e2) * eps3(i)
     7                 + k3e2 * (k2e3 + k3e3) * ck2(i)
     8                 - k2e3 * (k2e2 + k3e2) * ck3(i) )
c
     9           + 0.5*ckm(i) * lambdap * cepsilon(eps3,eps2,ckp,ckm)
      end do
C
      if (gauge .eq. 1.0) return      ! Feynman gauge
C
      denom = (1. - gauge) / (gauge * dot(kp,kp) - mv2)
C
      do i = 1, 4
          gama(i) = gama(i) + denom*cdot(gama,ckp)*ckp(i)
      end do
C
      Return
      End
c -------------------------------------------------------------------
CALCULATIONS FOR FOUR-VECTOR-BOSON COUPLING IN TERMS OF
C THE THREE POLARIZATIONS OF THE OUTGOING BOSONS
C
CONTRACTIONS WITH A PROPAGATOR IN UNITARY OR FEYNMANN GAUGES
C INPUT :
C      E1(4), E2(4), E3(4) --- THREE POLARIZATION VECTORS OF THEM
C      P(4)   --- MOMENTUM OF THE VIRTUAL BOSONS
C      EMV2   --- SQUARED VECTOR BOSON MASS APPEARED IN PROPAGATOR
C      GI     --- FLAG OF GAUGE CHOICE: 
C             GI=1: FEYNMANN GAUGE;
C              GI=0: UNITARY GAUGE;
C
C OUTPUT:
C      GAMA(RHO) --- S(LAM,MU,NU,DELT)*E1(MU)*E2(NU)*E3(DELT)*
C                   (G(RHO,LAM)-PROP(RHO)*PROP(LAM)/EMV2)
C
C NOTE: THE ORDER OF THE FOUR POLARIZATIONS OF BOSONS IS IMPORTANT!!
C      E1,E2; E3,E4 = Z(OR PHOTON), Z(OR PHOTON); W+, W-;
C NOTATIONS FOLLOWING THE BOOK :
C            FOUR-FOLD-VERTEX = -i*g_v*S(MU,NU;DELT,LAM)
      SUBROUTINE GAMA4(E1,E2,E3,P,EMV2,GI,GAMA)
      REAL P(4), E1(4), E2(4), E3(4), GAMA(4)
C
      D12  = DOT(E1,E2)
      D13  = DOT(E1,E3)
      D23  = DOT(E2,E3)
      DO I = 1, 4
         GAMA(I) = 2.*E1(I)*D23 - E2(I)*D13 - E3(I)*D12
      END DO
      IF (GI.EQ.1.0) RETURN
      DPP  = DOT(P,P)
      DPG  = DOT(P,GAMA)
      DO I = 1, 4
         GAMA(I) = GAMA(I) + (1.0-GI)*DPG/(GI*DPP-EMV2)*P(I)
      END DO
      RETURN
      END
c -------------------------------------------------------------------
CALCULATIONS FOR FOUR-VECTOR-BOSON COUPLING IN TERMS OF
C THE OUTGOING BOSON-COUPLED CURRENTS
C
C IN FEYNMANN GAUGE (IG=1) ONLY
C INPUT :
C      E1(4), E2(4), E3(4) --- THREE CURRENTS COUPLING TO THE THREE
C       VECTOR BOSONS ( COMPLEX )
C OUTPUT:
C      GAMA(RHO) --- S(LAM,MU,NU,DELT)*E1(MU)*E2(NU)*E3(DELT)*
C                   (G(RHO,LAM)-PROP(RHO)*PROP(LAM)/EMV2)
C                      ( COMPLEX )
C NOTE: THE ORDER OF THE FOUR POLARIZATIONS OF BOSONS IS IMPORTANT!!
C      E1,E2; E3,E4 = Z(OR PHOTON), Z(OR PHOTON); W+, W-;
C NOTATIONS FOLLOWING THE BOOK :
C            FOUR-FOLD-VERTEX = -i*g_v*S(MU,NU;DELT,LAM)
      SUBROUTINE CGAMA4(E1,E2,E3,CGAMA)
      COMPLEX E1(4), E2(4), E3(4), CGAMA(4),
     1      D12, D13, D23, CDOT
C
      D12  = CDOT(E1,E2)
      D13  = CDOT(E1,E3)
      D23  = CDOT(E2,E3)
      DO I = 1, 4
         CGAMA(I) = 2.*E1(I)*D23 - E2(I)*D13 - E3(I)*D12
      END DO
      RETURN
      END
c -------------------------------------------------------------------
C ***
C SUBROUTINE TO FORM THE COMPLEX-FOUR-VECTOR OF
C            CHII(DAGGER)*(SIGMA(4)_pauli)*CHIF
C WHERE CHII IS A LEFT-MOST TWO-COMPONENT-HELICITY EIGENSTATE; 
C      CHIF IS A RIGHT-MOST ONE; AND SIGM(4)_pauli IS THE 2*2 PAULI
C      MATRICES.
C
      SUBROUTINE SIG4(CHII,CHIF,SIG,ALPHA)
      INTEGER  ALPHA
      COMPLEX CHII(2), CHIF(2), SIG(4), CHIH(2), ZI
      PARAMETER ( ZI=(0E0,1E0) )

      CHIH(1) = CONJG(CHII(1))
      CHIH(2) = CONJG(CHII(2))

      SIG(1) =     CHIH(2) * CHIF(1) + CHIH(1) * CHIF(2)
      SIG(2) = ZI*(CHIH(2) * CHIF(1) - CHIH(1) * CHIF(2))
      SIG(3) =     CHIH(1) * CHIF(1) - CHIH(2) * CHIF(2)
      SIG(4) =     CHIH(1) * CHIF(1) + CHIH(2) * CHIF(2)

      IF (ALPHA.GT.0) RETURN

      DO I = 1, 3
         SIG(I) = -SIG(I)
      END DO
      RETURN
      END
c -------------------------------------------------------------------
C-----------------------------------------------------------------------
C   THESE PROGRAMS ARE BEING PREPARED BY K. HAGIWARA AND D. ZEPPENFELD.
C                  FIRST VERSION  :  DEC. 19, 1985   BY K. HAGIWARA
C                REVISED VERSION  :  JAN. 15, 1986   BY D.ZEPPENFELD
C-----------------------------------------------------------------------
C   FIRST, ONE SHOULD EVALUATE THE TWO COMPONENT SPINORS
C   CHIP = CHI(+,P) AND CHIM = CHI(-,P) AS FUNCTIONS OF A FOUR MOMENTUM
C   P.  SEE EQ.(3.22) AND EQ.(3.23).  THIS IS DONE BY THE SUBROUTINE
C
C        SUBROUTINE CHI(P,CHIP,CHIM) ( ACTUALLY, CHI(P,ALPH,CHII)-HAN)
C
C   ONCE WE HAVE EVALUATED ALL THE EXTERNAL SPINORS BY USING THE SUBROUTINE
C   CHI, WE ARE ABLE TO CALCULATE ARBITRARY FERMION STRINGS BY THE FUNCTION SN:
C
C        FUNCTION SN(CHII,A1,A2,...,AN,CHIF,ALPHA)
C
C   WITH MULIPLE ENTRY POINTS FOR S0, S1, ..., S9.  THE VALUE OF THE
C   STRING AS DEFINED BY EQ.(3.28) IS CALCULATED BY DIRECT MATRIX
C   MULTIPLICATION METHOD.  THE CHI'S CAN SHOULD BE SPECIFIED TO ONE OF
C
C         CHIPI = CHI(+,PI),  CHIPF = CHI(+,PF)
C         CHIMI = CHI(-,PI),  CHIMF = CHI(-,PF)
C
C   IN ORDER TO ACHIEVE NUMERICAL EFFICIENCY, ALL THE FOUR VECTORS
C   A1 TO AN MUST BE REAL IN THIS PROGRAM.
C
C   ALL FOUR VECTORS ARE ASSUMED TO HAVE TIME-COMPONENT AS FOURTH
C   ENTRY P(4).
C-----------------------------------------------------------------------
      FUNCTION S9(CHII,A1,A2,A3,A4,A5,A6,A7,A8,A9,CHIF,ALPHA)
      INTEGER  ALPHA, ALP, N, NN
      COMPLEX  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX  S0, S1, S2, S3, S4, S5, S6, S7, S8, S9
      REAL  A1(4), A2(4), A3(4), A4(4), A5(4), A6(4), A7(4), A8(4),
     *      A9(4), AUX(4,9)
C
      N = 9
      Nnum = 9
      GOTO 10
C
      ENTRY S0(CHII,CHIF,ALPHA)
      N = 0
      GOTO 20
C
      ENTRY S1(CHII,A1,CHIF,ALPHA)
      N = 1
      Nnum = 1
      GOTO 10
C
      ENTRY S2(CHII,A1,A2,CHIF,ALPHA)
      N = 2
      Nnum = 2
      GOTO 10
C
      ENTRY S3(CHII,A1,A2,A3,CHIF,ALPHA)
      N = 3
      Nnum = 3
      GOTO 10
C
      ENTRY S4(CHII,A1,A2,A3,A4,CHIF,ALPHA)
      N = 4
      Nnum = 4
      GOTO 10
C
      ENTRY S5(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      N = 5
      Nnum = 5
      GOTO 10
C
      ENTRY S6(CHII,A1,A2,A3,A4,A5,A6,CHIF,ALPHA)
      N = 6
      Nnum = 6
      GOTO 10
C
      ENTRY S7(CHII,A1,A2,A3,A4,A5,A6,A7,CHIF,ALPHA)
      N = 7
      Nnum = 7
      GOTO 10
C
      ENTRY S8(CHII,A1,A2,A3,A4,A5,A6,A7,A8,CHIF,ALPHA)
      N = 8
      Nnum = 8
      GOTO 10
C
 10   CONTINUE
         I = 1
         IF (Nnum .EQ. 1) THEN
            GOTO 1
         ELSE IF (Nnum .EQ. 2) THEN
            GOTO 2
         ELSE IF (Nnum .EQ. 3) THEN
            GOTO 3
         ELSE IF (Nnum .EQ. 4) THEN
            GOTO 4
         ELSE IF (Nnum .EQ. 5) THEN
            GOTO 5
         ELSE IF (Nnum .EQ. 6) THEN
            GOTO 6
         ELSE IF (Nnum .EQ. 7) THEN
            GOTO 7
         ELSE IF (Nnum .EQ. 8) THEN
            GOTO 8
         ELSE IF (Nnum .EQ. 9) THEN
            GOTO 9
         ENDIF
  9      AUX(I,9) = A9(I)
  8      AUX(I,8) = A8(I)
  7      AUX(I,7) = A7(I)
  6      AUX(I,6) = A6(I)
  5      AUX(I,5) = A5(I)
  4      AUX(I,4) = A4(I)
  3      AUX(I,3) = A3(I)
  2      AUX(I,2) = A2(I)
  1      AUX(I,1) = A1(I)
         I = I + 1
         IF (I.LE.4) THEN
          IF (Nnum .EQ. 1) THEN
            GOTO 1
          ELSE IF (Nnum .EQ. 2) THEN
            GOTO 2
          ELSE IF (Nnum .EQ. 3) THEN
            GOTO 3
          ELSE IF (Nnum .EQ. 4) THEN
            GOTO 4
          ELSE IF (Nnum .EQ. 5) THEN
            GOTO 5
          ELSE IF (Nnum .EQ. 6) THEN
            GOTO 6
          ELSE IF (Nnum .EQ. 7) THEN
            GOTO 7
          ELSE IF (Nnum .EQ. 8) THEN
            GOTO 8
          ELSE IF (Nnum .EQ. 9) THEN
            GOTO 9
          ENDIF
         ENDIF
   20 CONTINUE
      CHIAUX(1) = CONJG(CHII(1))
      CHIAUX(2) = CONJG(CHII(2))
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(4,I) - AUX(3,I)
            ASLASH(1,2) = - CMPLX( AUX(1,I),-AUX(2,I) )
            ASLASH(2,1) = - CMPLX( AUX(1,I), AUX(2,I) )
            ASLASH(2,2) = AUX(4,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(4,I) + AUX(3,I)
            ASLASH(1,2) = CMPLX( AUX(1,I),-AUX(2,I) )
            ASLASH(2,1) = CMPLX( AUX(1,I), AUX(2,I) )
            ASLASH(2,2) = AUX(4,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      S9  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
c
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C   THESE PROGRAMS ARE BEING PREPARED BY K. HAGIWARA AND D. ZEPPENFELD.
C                  FIRST VERSION  :  DEC. 19, 1985   BY K. HAGIWARA
C                REVISED VERSION  :  JAN. 15, 1986   BY D.ZEPPENFELD
C   MODIFIED BY T. HAN TO DO THE STRINGS FOR COMPLEX 4-VECTORS
C-----------------------------------------------------------------------
C   FIRST, ONE SHOULD EVALUATE THE TWO COMPONENT SPINORS
C   CHIP = CHI(+,P) AND CHIM = CHI(-,P) AS FUNCTIONS OF A FOUR MOMENTUM
C   P.  SEE EQ.(3.22) AND EQ.(3.23).  THIS IS DONE BY THE SUBROUTINE
C
C        SUBROUTINE CHI(P,CHIP,CHIM) ( ACTUALLY, CHI(P,ALPH,CHII)-HAN)
C
C   WITH MULIPLE ENTRY POINTS FOR CS0, CS1, ..., CS9.  THE VALUE OF THE
C   STRING AS DEFINED BY EQ.(3.28) IS CALCULATED BY DIRECT MATRIX
C   MULTIPLICATION METHOD.  THE CHI'S CAN BE SPECIFIED TO ONE OF
C
C         CHIPI = CHI(+,PI),  CHIPF = CHI(+,PF)
C         CHIMI = CHI(-,PI),  CHIMF = CHI(-,PF)
C
C   ALL THE FOUR VECTORS A1 TO AN ARE COMPLEX IN THIS PROGRAM &
C   ALL FOUR VECTORS ARE ASSUMED TO HAVE TIME-COMPONENT AS FOURTH
C   ENTRY P(4).
C-----------------------------------------------------------------------
      FUNCTION CS9(CHII,A1,A2,A3,A4,A5,A6,A7,A8,A9,CHIF,ALPHA)
      INTEGER  ALPHA, ALP, N, NN, Nnum
      COMPLEX  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX  CS0, CS1, CS2, CS3, CS4, CS5, CS6, CS7, CS8, CS9
      COMPLEX  A1(4), A2(4), A3(4), A4(4), A5(4), A6(4), A7(4),
     1      A8(4), A9(4), AUX(4,9), ZI
      PARAMETER ( ZI=(0D0,1D0) )
C
      N = 9
      Nnum = 9
      GOTO 10
C
      ENTRY CS0(CHII,CHIF,ALPHA)
      N = 0
      GOTO 20
C
      ENTRY CS1(CHII,A1,CHIF,ALPHA)
      N = 1
      Nnum = 1
      GOTO 10
C
      ENTRY CS2(CHII,A1,A2,CHIF,ALPHA)
      N = 2
      Nnum = 2
      GOTO 10
C
      ENTRY CS3(CHII,A1,A2,A3,CHIF,ALPHA)
      N = 3
      Nnum = 3
      GOTO 10
C
      ENTRY CS4(CHII,A1,A2,A3,A4,CHIF,ALPHA)
      N = 4
      Nnum = 4
      GOTO 10
C
      ENTRY CS5(CHII,A1,A2,A3,A4,A5,CHIF,ALPHA)
      N = 5
      Nnum = 5
      GOTO 10
C
      ENTRY CS6(CHII,A1,A2,A3,A4,A5,A6,CHIF,ALPHA)
      N = 6
      Nnum = 6
      GOTO 10
C
      ENTRY CS7(CHII,A1,A2,A3,A4,A5,A6,A7,CHIF,ALPHA)
      N = 7
      Nnum = 7
      GOTO 10
C
      ENTRY CS8(CHII,A1,A2,A3,A4,A5,A6,A7,A8,CHIF,ALPHA)
      N = 8
      Nnum = 8
      GOTO 10
C
 10   CONTINUE
         I = 1
         IF (Nnum .EQ. 1) THEN
            GOTO 1
         ELSE IF (Nnum .EQ. 2) THEN
            GOTO 2
         ELSE IF (Nnum .EQ. 3) THEN
            GOTO 3
         ELSE IF (Nnum .EQ. 4) THEN
            GOTO 4
         ELSE IF (Nnum .EQ. 5) THEN
            GOTO 5
         ELSE IF (Nnum .EQ. 6) THEN
            GOTO 6
         ELSE IF (Nnum .EQ. 7) THEN
            GOTO 7
         ELSE IF (Nnum .EQ. 8) THEN
            GOTO 8
         ELSE IF (Nnum .EQ. 9) THEN
            GOTO 9
         ENDIF
  9      AUX(I,9) = A9(I)
  8      AUX(I,8) = A8(I)
  7      AUX(I,7) = A7(I)
  6      AUX(I,6) = A6(I)
  5      AUX(I,5) = A5(I)
  4      AUX(I,4) = A4(I)
  3      AUX(I,3) = A3(I)
  2      AUX(I,2) = A2(I)
  1      AUX(I,1) = A1(I)
         I = I + 1
         IF (I.LE.4) THEN
          IF (Nnum .EQ. 1) THEN
            GOTO 1
          ELSE IF (Nnum .EQ. 2) THEN
            GOTO 2
          ELSE IF (Nnum .EQ. 3) THEN
            GOTO 3
          ELSE IF (Nnum .EQ. 4) THEN
            GOTO 4
          ELSE IF (Nnum .EQ. 5) THEN
            GOTO 5
          ELSE IF (Nnum .EQ. 6) THEN
            GOTO 6
          ELSE IF (Nnum .EQ. 7) THEN
            GOTO 7
          ELSE IF (Nnum .EQ. 8) THEN
            GOTO 8
          ELSE IF (Nnum .EQ. 9) THEN
            GOTO 9
          ENDIF
         ENDIF
   20 CONTINUE
      CHIAUX(1) = CONJG(CHII(1))
      CHIAUX(2) = CONJG(CHII(2))
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(4,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(4,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(4,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(4,I) - AUX(3,I)
         ENDIF
         CHIDUM    = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP       = - ALP
 30   CONTINUE
C
      CS9  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END
C-----------------------------------------------------------------------
C   THESE PROGRAMS ARE PREPARED BY T. Han
C   First VERSION  :  Oct. 28, 1993
C-----------------------------------------------------------------------
C   FIRST, ONE SHOULD HAVE EVALUATED THE TWO COMPONENT SPINORS
C   CHIP = CHI(+,P) AND CHIM = CHI(-,P) AS FUNCTIONS OF A FOUR MOMENTUM
C   P.  SEE EQ.(3.22) AND EQ.(3.23).  THIS IS DONE BY THE SUBROUTINE
C
C        SUBROUTINE CHI(P,CHIP,CHIM) ( ACTUALLY, CHI(P,ALPH,CHII)-HAN)
C
C   WITH MULIPLE ENTRY POINTS FOR CSleft0, CSleft1, ..., CSleft9.
C   THE VALUE OF THE STRING IS CALCULATED BY DIRECT MATRIX MULTIPLICATION 
C   THE CHI'S CAN BE SPECIFIED TO ONE OF
C
C         CHIPI = CHI(+,PI),  CHIMI = CHI(-,PI)
C
C   OUTPUT :
C        CHIO = [ CHI^DAGGER*AN... ]^DAGGER
C         READY TO BE USED AS A NEW CHII
C
C   ALL THE FOUR VECTORS A1 TO AN ARE COMPLEX IN THIS PROGRAM &
C   ALL FOUR VECTORS ARE ASSUMED TO HAVE TIME-COMPONENT AS FOURTH
C   ENTRY P(4).
C-----------------------------------------------------------------------
      SUBROUTINE CSleft9(CHII,A1,A2,A3,A4,A5,A6,A7,A8,A9,CHIO,ALPHA)
C
      Implicit None
C
      INTEGER  ALPHA, ALP, N, NN, I, Nnum
      COMPLEX  CHII(2), CHIO(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX  A1(4), A2(4), A3(4), A4(4), A5(4), A6(4), A7(4),
     1      A8(4), A9(4), AUX(4,9), ZI
      PARAMETER ( ZI=(0D0,1D0) )
C
      N = 9
      Nnum = 9
      GOTO 10
C
      ENTRY CSLEFT1(CHII,A1,CHIO,ALPHA)
      N = 1
      Nnum = 1
      GOTO 10
C
      ENTRY CSLEFT2(CHII,A1,A2,CHIO,ALPHA)
      N = 2
      Nnum = 2
      GOTO 10
C
      ENTRY CSLEFT3(CHII,A1,A2,A3,CHIO,ALPHA)
      N = 3
      Nnum = 3
      GOTO 10
C
      ENTRY CSLEFT4(CHII,A1,A2,A3,A4,CHIO,ALPHA)
      N = 4
      Nnum = 4
      GOTO 10
C
      ENTRY CSLEFT5(CHII,A1,A2,A3,A4,A5,CHIO,ALPHA)
      N = 5
      Nnum = 5
      GOTO 10
C
      ENTRY CSLEFT6(CHII,A1,A2,A3,A4,A5,A6,CHIO,ALPHA)
      N = 6
      Nnum = 6
      GOTO 10
C
      ENTRY CSLEFT7(CHII,A1,A2,A3,A4,A5,A6,A7,CHIO,ALPHA)
      N = 7
      Nnum = 7
      GOTO 10
C
      ENTRY CSLEFT8(CHII,A1,A2,A3,A4,A5,A6,A7,A8,CHIO,ALPHA)
      N = 8
      Nnum = 8
      GOTO 10
C
 10   CONTINUE
         I = 1
         IF (Nnum .EQ. 1) THEN
            GOTO 1
         ELSE IF (Nnum .EQ. 2) THEN
            GOTO 2
         ELSE IF (Nnum .EQ. 3) THEN
            GOTO 3
         ELSE IF (Nnum .EQ. 4) THEN
            GOTO 4
         ELSE IF (Nnum .EQ. 5) THEN
            GOTO 5
         ELSE IF (Nnum .EQ. 6) THEN
            GOTO 6
         ELSE IF (Nnum .EQ. 7) THEN
            GOTO 7
         ELSE IF (Nnum .EQ. 8) THEN
            GOTO 8
         ELSE IF (Nnum .EQ. 9) THEN
            GOTO 9
         ENDIF
  9      AUX(I,9) = A9(I)
  8      AUX(I,8) = A8(I)
  7      AUX(I,7) = A7(I)
  6      AUX(I,6) = A6(I)
  5      AUX(I,5) = A5(I)
  4      AUX(I,4) = A4(I)
  3      AUX(I,3) = A3(I)
  2      AUX(I,2) = A2(I)
  1      AUX(I,1) = A1(I)
         I = I + 1
         IF (I.LE.4) THEN
          IF (Nnum .EQ. 1) THEN
            GOTO 1
          ELSE IF (Nnum .EQ. 2) THEN
            GOTO 2
          ELSE IF (Nnum .EQ. 3) THEN
            GOTO 3
          ELSE IF (Nnum .EQ. 4) THEN
            GOTO 4
          ELSE IF (Nnum .EQ. 5) THEN
            GOTO 5
          ELSE IF (Nnum .EQ. 6) THEN
            GOTO 6
          ELSE IF (Nnum .EQ. 7) THEN
            GOTO 7
          ELSE IF (Nnum .EQ. 8) THEN
            GOTO 8
          ELSE IF (Nnum .EQ. 9) THEN
            GOTO 9
          ENDIF
         ENDIF
   20 CONTINUE
      CHIAUX(1) = CONJG(CHII(1))
      CHIAUX(2) = CONJG(CHII(2))
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) =  AUX(4,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) =  AUX(4,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(4,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(4,I) - AUX(3,I)
         ENDIF
         CHIDUM    = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP       = - ALP
 30   CONTINUE
C
      CHIO(1) = CONJG( CHIAUX(1) )
      CHIO(2) = CONJG( CHIAUX(2) )
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C   THESE PROGRAMS ARE BEING PREPARED BY T. Han
C   First VERSION  :  Oct. 28, 1993
C-----------------------------------------------------------------------
C   FIRST, ONE SHOULD HAVE EVALUATED THE TWO COMPONENT SPINORS
C   CHIP = CHI(+,P) AND CHIM = CHI(-,P) AS FUNCTIONS OF A FOUR MOMENTUM
C   P.  SEE EQ.(3.22) AND EQ.(3.23).  THIS IS DONE BY THE SUBROUTINE
C
C        SUBROUTINE CHI(P,CHIP,CHIM) ( ACTUALLY, CHI(P,ALPH,CHII)-HAN)
C
C   WITH MULIPLE ENTRY POINTS FOR CSright0, CSright1, ..., CSright9.
C   THE VALUE OF THE STRING IS CALCULATED BY DIRECT MATRIX MULTIPLICATION 
C   THE CHI'S CAN BE SPECIFIED TO ONE OF
C
C         CHIPI = CHI(+,PI),  CHIMI = CHI(-,PI)
C
C   OUTPUT :
C        CHIO = AN...A1*CHI 
C         READY TO BE USED AS A NEW CHIF on the right-end of the string
C
C   ALL THE FOUR VECTORS A1 TO AN ARE COMPLEX IN THIS PROGRAM &
C   ALL FOUR VECTORS ARE ASSUMED TO HAVE TIME-COMPONENT AS FOURTH
C   ENTRY P(4).
C-----------------------------------------------------------------------
      SUBROUTINE CSright9(CHIO,A9,A8,A7,A6,A5,A4,A3,A2,A1,CHII,ALPHA)
C
      Implicit None
C
      INTEGER  ALPHA, ALP, N, NN, I, Nnum
      COMPLEX  CHII(2), CHIO(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX  A1(4), A2(4), A3(4), A4(4), A5(4), A6(4), A7(4),
     1      A8(4), A9(4), AUX(4,9), ZI
      PARAMETER ( ZI=(0D0,1D0) )
C
      N = 9
      Nnum = 9
      GOTO 10
C
      ENTRY CSright1(CHIO,A1,CHII,ALPHA)
      N = 1
      Nnum = 1
      GOTO 10
C
      ENTRY CSright2(CHIO,A2,A1,CHII,ALPHA)
      N = 2
      Nnum = 2
      GOTO 10
C
      ENTRY CSright3(CHIO,A3,A2,A1,CHII,ALPHA)
      N = 3
      Nnum = 3
      GOTO 10
C
      ENTRY CSright4(CHIO,A4,A3,A2,A1,CHII,ALPHA)
      N = 4
      Nnum = 4
      GOTO 10
C
      ENTRY CSright5(CHIO,A5,A4,A3,A2,A1,CHII,ALPHA)
      N = 5
      Nnum = 5
      GOTO 10
C
      ENTRY CSright6(CHIO,A6,A5,A4,A3,A2,A1,CHII,ALPHA)
      N = 6
      Nnum = 6
      GOTO 10
C
      ENTRY CSright7(CHIO,A7,A6,A5,A4,A3,A2,A1,CHII,ALPHA)
      N = 7
      Nnum = 7
      GOTO 10
C
      ENTRY CSright8(CHIO,A8,A7,A6,A5,A4,A3,A2,A1,CHII,ALPHA)
      N = 8
      Nnum = 8
      GOTO 10
C
 10   CONTINUE
         I = 1
         IF (Nnum .EQ. 1) THEN
            GOTO 1
         ELSE IF (Nnum .EQ. 2) THEN
            GOTO 2
         ELSE IF (Nnum .EQ. 3) THEN
            GOTO 3
         ELSE IF (Nnum .EQ. 4) THEN
            GOTO 4
         ELSE IF (Nnum .EQ. 5) THEN
            GOTO 5
         ELSE IF (Nnum .EQ. 6) THEN
            GOTO 6
         ELSE IF (Nnum .EQ. 7) THEN
            GOTO 7
         ELSE IF (Nnum .EQ. 8) THEN
            GOTO 8
         ENDIF
  9      AUX(I,9) = A9(I)
  8      AUX(I,8) = A8(I)
  7      AUX(I,7) = A7(I)
  6      AUX(I,6) = A6(I)
  5      AUX(I,5) = A5(I)
  4      AUX(I,4) = A4(I)
  3      AUX(I,3) = A3(I)
  2      AUX(I,2) = A2(I)
  1      AUX(I,1) = A1(I)
         I = I + 1
         IF (I.LE.4) THEN
          IF (Nnum .EQ. 1) THEN
            GOTO 1
          ELSE IF (Nnum .EQ. 2) THEN
            GOTO 2
          ELSE IF (Nnum .EQ. 3) THEN
            GOTO 3
          ELSE IF (Nnum .EQ. 4) THEN
            GOTO 4
          ELSE IF (Nnum .EQ. 5) THEN
            GOTO 5
          ELSE IF (Nnum .EQ. 6) THEN
            GOTO 6
          ELSE IF (Nnum .EQ. 7) THEN
            GOTO 7
          ELSE IF (Nnum .EQ. 8) THEN
            GOTO 8
          ENDIF
         ENDIF
   20 CONTINUE
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) =  AUX(4,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) =  AUX(4,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(4,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(4,I) - AUX(3,I)
         ENDIF
         CHIDUM    = ASLASH(1,1)*CHIAUX(1) + ASLASH(1,2)*CHIAUX(2)
         CHIAUX(2) = ASLASH(2,1)*CHIAUX(1) + ASLASH(2,2)*CHIAUX(2)
         CHIAUX(1) = CHIDUM
         ALP       = - ALP
 30   CONTINUE
C
      CHIO(1) = CHIAUX(1)
      CHIO(2) = CHIAUX(2)
C
      RETURN
      END
C-----------------------------------------------------------------------
c   (2). b: 4x4 formulism [a la Alan Stange's Hlpack8.for package]
C-----------------------------------------------------------------------
c
complex*8 ( which is equal to complex=[real*4,real*4] );
c     with real*8 throughout these codes;
c     but 
c      1). mass "m" in mi's is in real*4 
c          --- inconsistency! well, OK for complex*8
c      2). ab_gamma & current have real a, b, x, y
c
comments made on 06/18/93
c
      subroutine epsn(p,m,l,epsilon)
      implicit none
      integer l                 ! as input vector boson helicity 
                                ! in Cartesian basis 1,2,3(3=longitd)
      real*8 p(4), m, epsilon(4)
      real*8 pt2, pt, pabs, fac
                                !D. Zeppenfeld's eps() routine copied here.
                                !rewritten for readability and real*8
      pt2  = p(1)**2 + p(2)**2
      pt   = sqrt(pt2)
      pabs = sqrt(pt2 + p(3)**2)
      epsilon(4) = 0d0
      if (l .eq. 1) then
         if (pt .le. 1d-20*abs(p(3))) then
            epsilon(1) = p(3) / abs(p(3))
            epsilon(2) = 0d0
            epsilon(3) = 0d0
         else
            fac = 1d0 / (pabs * pt)
            epsilon(1) = p(1) * p(3) * fac
            epsilon(2) = p(2) * p(3) * fac
            epsilon(3) = - pt2 * fac
         endif
      else if (l .eq. 2) then
         if (pt .le. 1d-20*abs(p(3))) then
            epsilon(1) = 0d0
            epsilon(2) = 1d0
            epsilon(3) = 0d0
         else
            fac = 1d0 / pt
            epsilon(1) = - p(2) * fac
            epsilon(2) =   p(1) * fac
            epsilon(3) = 0d0
         endif
      else if (l .eq. 3) then
         fac = p(4) / (m * pabs)
         epsilon(1) = p(1) * fac
         epsilon(2) = p(2) * fac
         epsilon(3) = p(3) * fac
         epsilon(4) = pabs / m
      else
         write(*,*) 'eps(): Undefined value: l', l
      endif
      end


      subroutine u_spinor(p, spin, u)
calculate the Dirac Spinor u as a 4by1 complex matrix
c for a given helicity called spin=+-1 and momentum real*8 p(mu)
      implicit none
      real*8 p(4)
      integer spin        ! as input helicity of the fermion lambd=+-1
      complex u(4)
      
      real*8 q(4), p_mag, omega_U, omega_D, eps, pabz, coeff
      parameter (eps = 1d-20)
      
      if ((spin .ne. 1) .and. (spin .ne. -1)) then
         write(*,*) 'u_spinor: bad spin value', spin
         call exit(1)
      endif
      if (p(4) .lt. 0d0) then
         q(1) = - p(1)
         q(2) = - p(2)
         q(3) = - p(3)
         q(4) = - p(4)
         call v_spinor(q, -spin, u)
      else
         p_mag = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
         pabz = p_mag + p(3)
         omega_U = sqrt(max(0d0, p(4) - spin * p_mag))
         omega_D = sqrt(max(0d0, p(4) + spin * p_mag))
         if (pabz .gt. eps*p_mag) then
          coeff = 1d0 / sqrt(2 * p_mag * pabz)
          if (spin .eq. 1) then
               u(1) = omega_U * coeff * pabz
               u(2) = omega_U * coeff * cmplx(p(1), p(2))
               u(3) = omega_D * coeff * pabz
               u(4) = omega_D * coeff * cmplx(p(1), p(2))
          else
               u(1) = omega_U * coeff * cmplx(-p(1), p(2))
               u(2) = omega_U * coeff * pabz
               u(3) = omega_D * coeff * cmplx(-p(1), p(2))
               u(4) = omega_D * coeff * pabz
          endif
         else
          if (spin .eq. 1) then
               u(1) = 0d0
               u(2) = omega_U
               u(3) = 0d0
               u(4) = omega_D
          else
               u(1) = - omega_U
               u(2) = 0d0
               u(3) = - omega_D
               u(4) = 0d0
          endif
         endif
      endif
      end

      
      subroutine v_spinor(p, spin, v)
calculate the Dirac Spinor v as a 4by1 complex matrix
c for a given helicity called spin=+-1 and momentum real*8 p(mu)
      implicit none
      real*8 p(4)
      integer spin
      complex v(4)
      
      real*8 q(4), p_mag, omega_U, omega_D, eps, pabz, coeff
      parameter (eps = 1d-20)
      
      if ((spin .ne. 1) .and. (spin .ne. -1)) then
         write(*,*) 'v_spinor: bad spin value', spin
         call exit(1)
      endif
      if (p(4) .lt. 0d0) then
         q(1) = - p(1)
         q(2) = - p(2)
         q(3) = - p(3)
         q(4) = - p(4)
         call u_spinor(q, -spin, v)
      else
         p_mag = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
         pabz = p_mag + p(3)
         omega_U = sqrt(max(0d0, p(4) + spin * p_mag))
         omega_D = sqrt(max(0d0, p(4) - spin * p_mag))
         if (pabz .gt. eps*p_mag) then
          coeff = 1d0 / sqrt(2 * p_mag * pabz)
          if (spin .eq. 1) then
               v(1) = - omega_U * coeff * cmplx(-p(1), p(2))
               v(2) = - omega_U * coeff * pabz
               v(3) =   omega_D * coeff * cmplx(-p(1), p(2))
               v(4) =   omega_D * coeff * pabz
          else
               v(1) =   omega_U * coeff * pabz
               v(2) =   omega_U * coeff * cmplx(p(1), p(2))
               v(3) = - omega_D * coeff * pabz
               v(4) = - omega_D * coeff * cmplx(p(1), p(2))
          endif
         else
          if (spin .eq. 1) then
               v(1) =   omega_U
               v(2) =   0d0
               v(3) = - omega_D
               v(4) =   0d0
          else
               v(1) =   0d0
               v(2) =   omega_U
               v(3) =   0d0
               v(4) =   - omega_D
          endif
         endif
      endif
      end


      subroutine slash(p, ps, mass)
calculate:  p(mu)*gamma(mu) + mass as a complex 4by4 pslash 
c for a given input p(mu), a real*8 arbitrary 4-vector
      implicit none
      complex ps(4, 4)
      real*8 p(4), mass
      
      complex i, zero
      parameter (i = (0., 1.), zero = (0., 0.))
      
      ps(1,1) = mass
      ps(1,2) = zero
      ps(1,3) = p(4) - p(3)
      ps(1,4) = - (p(1) - i * p(2))
      
      ps(2, 1) = zero
      ps(2, 2) = mass
      ps(2, 3) = -(p(1) + i * p(2))
      ps(2, 4) = p(4) + p(3)
      
      ps(3, 1) = p(4) + p(3)
      ps(3, 2) = p(1) - i * p(2)
      ps(3, 3) = mass
      ps(3, 4) = zero
      
      ps(4, 1) = p(1) + i * p(2)
      ps(4, 2) = p(4) - p(3)
      ps(4, 3) = zero
      ps(4, 4) = mass
      
      end
      
      
      subroutine slash_c(p, ps, mass)
calculate:  p(mu)*gamma(mu) + mass as a complex 4by4 pslash 
c for a given input p(mu), an arbitrary complex 4-vector
      implicit none
      complex p(4), ps(4, 4)
      real*8 mass
      
      complex i, zero
      parameter (i = (0., 1.), zero = (0., 0.))
      
      ps(1,1) = mass
      ps(1,2) = zero
      ps(1,3) = p(4) - p(3)
      ps(1, 4) = - (p(1) - i * p(2))
      
      ps(2, 1) = zero
      ps(2, 2) = mass
      ps(2, 3) = -(p(1) + i * p(2))
      ps(2, 4) = p(4) + p(3)
      
      ps(3, 1) = p(4) + p(3)
      ps(3, 2) = p(1) - i * p(2)
      ps(3, 3) = mass
      ps(3, 4) = zero
      
      ps(4, 1) = p(1) + i * p(2)
      ps(4, 2) = p(4) - p(3)
      ps(4, 3) = zero
      ps(4, 4) = mass
      
      end


      complex function m0(psi1, psi2)
calculate:    psi1^bar * psi2 as a complex number:
c for given 4by1 complex spinors psi1, psi2, and could be from,
c     say, calling "l_mult..." or "r_mult..."
c note that psi^bar = psi^dagger * gamma0
      implicit none
      complex psi1(4), psi2(4)
      
      m0 = conjg(psi1(3)) * psi2(1) + conjg(psi1(4)) * psi2(2) +
     .     conjg(psi1(1)) * psi2(3) + conjg(psi1(2)) * psi2(4)
      end


      complex function m1(psi1, a, q)
c
calculate:    psi1^bar * (a^slash+m) * q as a complex number:
comments :    NOTE that 
c (1) in m1--->m3: inputs "a, b, c" are NOT arbitrary 4by4 matrices!
c     rather a contraction a(mu)*gamma(mu) + m, obtained from, 
c     say, calling "slash"
c (2) psi1 & q(4) are arbitrary complex vectors and could be from,
c     say, calling "l_mult..." & "r_mult..."
c
      implicit none
      complex psi1(4), q(4), a(4, 4)  
      real m

      m = a(1,1)
      m1 = conjg(psi1(3)) * (a(1,3)*q(3) + a(1,4)*q(4) + m * q(1)) +
     .     conjg(psi1(4)) * (a(2,3)*q(3) + a(2,4)*q(4) + m * q(2)) +
     .     conjg(psi1(1)) * (a(3,1)*q(1) + a(3,2)*q(2) + m * q(3)) +
     .     conjg(psi1(2)) * (a(4,1)*q(1) + a(4,2)*q(2) + m * q(4))
      end


      complex function mm1(psi1, a, q)
      implicit none
comments on "m1" still valid !!!
      complex psi1(4), q(4), a(4, 4)  
                                !same as m1() but with the mass terms
                                !removed for speed reasons.  mm1() can
                                !often be used since all polarization
                                !vectors and currents need no mass term
      mm1 = conjg(psi1(3)) * (a(1,3)*q(3) + a(1,4)*q(4)) +
     .      conjg(psi1(4)) * (a(2,3)*q(3) + a(2,4)*q(4)) +
     .      conjg(psi1(1)) * (a(3,1)*q(1) + a(3,2)*q(2)) +
     .      conjg(psi1(2)) * (a(4,1)*q(1) + a(4,2)*q(2))
      end


      complex function m2(p1, b, a, p2)
calculate:    p1^bar * (b^slash+m1) (a^slash+m2) * p2 as a complex number:
comments on "m1" still valid !!!
      implicit none
      complex p1(4), a(4,4), b(4,4), p2(4)

      real m
      complex p(4)

      m = a(1,1)
      p(1) = m*p2(1) + a(1,3)*p2(3) + a(1,4)*p2(4)
      p(2) = m*p2(2) + a(2,3)*p2(3) + a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + m*p2(3)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + m*p2(4)
      
      m = b(1,1)
      m2 = conjg(p1(3)) * (m*p(1) + b(1,3)*p(3) + b(1,4)*p(4)) +
     .     conjg(p1(4)) * (m*p(2) + b(2,3)*p(3) + b(2,4)*p(4)) +
     .     conjg(p1(1)) * (b(3,1)*p(1) + b(3,2)*p(2) + m*p(3)) +
     .     conjg(p1(2)) * (b(4,1)*p(1) + b(4,2)*p(2) + m*p(4))
      end


      complex function m3(p1, c, b, a, p2)
calculate p1^bar*(c^slash+m1)(b^slash+m2)(a^slash+m3)*p2 as a complex number:
comments on "m1" still valid !!!
      implicit none
      complex p1(4), c(4,4), b(4,4), a(4,4), p2(4)

      real m
      complex y(4), x(4)

      m = a(1,1)
      y(1) = m*p2(1) + a(1,3)*p2(3) + a(1,4)*p2(4)
      y(2) = m*p2(2) + a(2,3)*p2(3) + a(2,4)*p2(4)
      y(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + m*p2(3)
      y(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + m*p2(4)
      
      m = b(1,1)
      x(1) = m*y(1) + b(1,3)*y(3) + b(1,4)*y(4)
      x(2) = m*y(2) + b(2,3)*y(3) + b(2,4)*y(4)
      x(3) = b(3,1)*y(1) + b(3,2)*y(2) + m*y(3)
      x(4) = b(4,1)*y(1) + b(4,2)*y(2) + m*y(4)
      
      m = c(1,1)
      m3 = conjg(p1(3))*(m*x(1) + c(1,3)*x(3) + c(1,4)*x(4)) +
     .     conjg(p1(4))*(m*x(2) + c(2,3)*x(3) + c(2,4)*x(4)) +
     .     conjg(p1(1))*(c(3,1)*x(1) + c(3,2)*x(2) + m*x(3)) +
     .     conjg(p1(2))*(c(4,1)*x(1) + c(4,2)*x(2) + m*x(4))
      end


      complex function matrix1(p1, a, p2)
c More general version of "m1",
c here "a" is an arbitrary 4by4 matrix, no need to be a slashed:
calculate the psi1^bar * a * q as a complex number:
c
      implicit none
      complex p1(4), a(4,4), p2(4), p(4)

      p(1) = a(1,1)*p2(1) + a(1,2)*p2(2) + a(1,3)*p2(3) + a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1) + a(2,2)*p2(2) + a(2,3)*p2(3) + a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + a(3,3)*p2(3) + a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + a(4,3)*p2(3) + a(4,4)*p2(4)
      matrix1 = conjg(p1(3)) * p(1) + conjg(p1(4)) * p(2) +
     .          conjg(p1(1)) * p(3) + conjg(p1(2)) * p(4)
      end


      complex function matrix2(p1, b, a, p2)
c More general version of "m2"
      implicit none
      complex p1(4), a(4,4), b(4,4), p2(4), p(4), matrix1

      p(1) = a(1,1)*p2(1) + a(1,2)*p2(2) + a(1,3)*p2(3) + a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1) + a(2,2)*p2(2) + a(2,3)*p2(3) + a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + a(3,3)*p2(3) + a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + a(4,3)*p2(3) + a(4,4)*p2(4)
      
      matrix2 = matrix1(p1, b, p)
      end


      complex function matrix3(p1, c, b, a, p2)
c More general version of "m3"
      implicit none
      complex p1(4), a(4,4), b(4,4), c(4,4), p2(4), p(4), matrix2

      p(1) = a(1,1)*p2(1) + a(1,2)*p2(2) + a(1,3)*p2(3) + a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1) + a(2,2)*p2(2) + a(2,3)*p2(3) + a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + a(3,3)*p2(3) + a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + a(4,3)*p2(3) + a(4,4)*p2(4)
      
      matrix3 = matrix2(p1, c, b, p)
      end


      complex function m4(psi1, d, c, b, a, p2)
      implicit none
c From this m4 on till m9, there is NO assumption for a^slash form
c because it calls the arbitrary 4by4 MATRIX3 etc, so it
calculates   psi1^bar * d*c*b*a * p2
c
      complex psi1(4), d(4,4), c(4,4), b(4,4), a(4,4), p2(4), p(4)
      complex matrix3

      p(1) = a(1,1)*p2(1)+a(1,2)*p2(2) + a(1,3)*p2(3)+a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1)+a(2,2)*p2(2) + a(2,3)*p2(3)+a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1)+a(3,2)*p2(2) + a(3,3)*p2(3)+a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1)+a(4,2)*p2(2) + a(4,3)*p2(3)+a(4,4)*p2(4)
      
      m4 = matrix3(psi1, d, c, b, p)
      end


      complex function m5(psi1, e, d, c, b, a, p2)
      implicit none
      complex psi1(4), e(4,4), d(4,4), c(4,4), b(4,4), a(4,4), p2(4), 
     1        p(4), m4
      
      p(1) = a(1,1)*p2(1) + a(1,2)*p2(2) + a(1,3)*p2(3)+a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1) + a(2,2)*p2(2) + a(2,3)*p2(3)+a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + a(3,3)*p2(3)+a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + a(4,3)*p2(3)+a(4,4)*p2(4)
      
      m5 = m4(psi1, e, d, c, b, p)
      end


      complex function m6(psi1, f, e, d, c, b, a, p2)
      implicit none
      complex psi1(4), e(4,4), d(4,4), c(4,4), b(4,4), a(4,4), p2(4)
      complex f(4, 4), p(4), m5     

      p(1) = a(1,1)*p2(1) + a(1,2)*p2(2) + a(1,3)*p2(3)+a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1) + a(2,2)*p2(2) + a(2,3)*p2(3)+a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + a(3,3)*p2(3)+a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + a(4,3)*p2(3)+a(4,4)*p2(4)

      m6 = m5(psi1, f, e, d, c, b, p)
      end


      complex function m7(psi1, g, f, e, d, c, b, a, p2)
      implicit none
      complex psi1(4), e(4,4), d(4,4), c(4,4), b(4,4), a(4,4), p2(4)
      complex f(4, 4), g(4, 4), p(4), m6
      
      p(1) = a(1,1)*p2(1) + a(1,2)*p2(2) + a(1,3)*p2(3)+a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1) + a(2,2)*p2(2) + a(2,3)*p2(3)+a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + a(3,3)*p2(3)+a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + a(4,3)*p2(3)+a(4,4)*p2(4)
      
      m7 = m6(psi1, g, f, e, d, c, b, p)
      end


      complex function m8(psi1, h, g, f, e, d, c, b, a, p2)
      implicit none
      complex psi1(4), e(4,4), d(4,4), c(4,4), b(4,4), a(4,4), p2(4)
      complex h(4, 4), g(4, 4), f(4, 4), p(4), m7
      
      p(1) = a(1,1)*p2(1) + a(1,2)*p2(2) + a(1,3)*p2(3)+a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1) + a(2,2)*p2(2) + a(2,3)*p2(3)+a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + a(3,3)*p2(3)+a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + a(4,3)*p2(3)+a(4,4)*p2(4)
      
      m8 = m7(psi1, h, g, f, e, d, c, b, p)
      end


      complex function m9(psi1, i, h, g, f, e, d, c, b, a, p2)
      implicit none
      complex psi1(4), e(4,4), d(4,4), c(4,4), b(4,4), a(4,4), p2(4)
      complex i(4, 4), h(4, 4), g(4, 4), f(4, 4), p(4), m8
      
      p(1) = a(1,1)*p2(1) + a(1,2)*p2(2) + a(1,3)*p2(3)+a(1,4)*p2(4)
      p(2) = a(2,1)*p2(1) + a(2,2)*p2(2) + a(2,3)*p2(3)+a(2,4)*p2(4)
      p(3) = a(3,1)*p2(1) + a(3,2)*p2(2) + a(3,3)*p2(3)+a(3,4)*p2(4)
      p(4) = a(4,1)*p2(1) + a(4,2)*p2(2) + a(4,3)*p2(3)+a(4,4)*p2(4)
      
      m9 = m8(psi1, i, h, g, f, e, d, c, b, p)
      end



      subroutine l_mult1(p, a, out)
calculate the useful Left-part complex 4-vector as
c out = [p^bar*(a^slash+m)]^bar to be ready as an input new "psi1" to mi's
comments :    NOTE that in these "l_mult.." and "r_mult...": 
c     like in m1--->m3, inputs "a, b" are NOT arbitrary 4by4 matrices!
c     rather a contraction a(mu)*gamma(mu) + m, obtained from, 
c     say, calling "slash"
c    
      implicit none
      complex p(4), a(4,4), out(4), u1, u2, u3, u4
      real m

      u1 = conjg(p(3))
      u2 = conjg(p(4))
      u3 = conjg(p(1))
      u4 = conjg(p(2))
      
      m = a(1,1)
      out(1) = conjg(m * u3 + u1*a(1,3) + u2*a(2,3))
      out(2) = conjg(m * u4 + u1*a(1,4) + u2*a(2,4))
      out(3) = conjg(m * u1 + u3*a(3,1) + u4*a(4,1))
      out(4) = conjg(m * u2 + u3*a(3,2) + u4*a(4,2))
      end


      subroutine r_mult1(a, p, out)
calculate the useful right-part complex 4-vector as
c (a^slash+m)*p to be ready as an input new "psi2" to mi's
      implicit none
      complex p(4), a(4,4), out(4)
      real m

      m = a(1,1)
      out(1) = m * p(1) + a(1,3)*p(3) + a(1,4)*p(4)
      out(2) = m * p(2) + a(2,3)*p(3) + a(2,4)*p(4)
      out(3) = m * p(3) + a(3,1)*p(1) + a(3,2)*p(2)
      out(4) = m * p(4) + a(4,1)*p(1) + a(4,2)*p(2)
      end


      subroutine l_mult2(u, a, b, out)
comments in l_mult1 still valid !!!
      implicit none
      complex u(4), a(4,4), b(4,4), out(4), p(4), x1, x2, x3, x4
      real m

      p(1) = conjg(u(3))
      p(2) = conjg(u(4))
      p(3) = conjg(u(1))
      p(4) = conjg(u(2))

      m = a(1,1)
      x1 = m*p(1) + p(3)*a(3,1) + p(4)*a(4,1)
      x2 = m*p(2) + p(3)*a(3,2) + p(4)*a(4,2)
      x3 = m*p(3) + p(1)*a(1,3) + p(2)*a(2,3)
      x4 = m*p(4) + p(1)*a(1,4) + p(2)*a(2,4)
      
      m = b(1,1)
      out(1) = conjg(m*x3 + x1*b(1,3) + x2*b(2,3))
      out(2) = conjg(m*x4 + x1*b(1,4) + x2*b(2,4))
      out(3) = conjg(m*x1 + x3*b(3,1) + x4*b(4,1))
      out(4) = conjg(m*x2 + x3*b(3,2) + x4*b(4,2))
      end


      subroutine r_mult2(a, b, p, out)
comments in l_mult1 still valid !!!
      implicit none
      complex p(4), a(4,4), b(4,4), out(4), x1, x2, x3, x4
      real m

      m = b(1,1)
      x1 = m * p(3) + b(3,1) * p(1) + b(3,2) * p(2)
      x2 = m * p(4) + b(4,1) * p(1) + b(4,2) * p(2)
      x3 = m * p(1) + b(1,3) * p(3) + b(1,4) * p(4)
      x4 = m * p(2) + b(2,3) * p(3) + b(2,4) * p(4)
      
      m = a(1,1)
      out(1) = m * x3 + a(1,3) * x1 + a(1,4) * x2
      out(2) = m * x4 + a(2,3) * x1 + a(2,4) * x2
      out(3) = m * x1 + a(3,1) * x3 + a(3,2) * x4
      out(4) = m * x2 + a(4,1) * x3 + a(4,2) * x4
      end


      subroutine ab_gamma5(in, a, b, out)
calculate out = in * (a + b*gamma5)
c with arbitrary 4by4 "in" and couplings can be gv=a, ga=b
      implicit none
      complex in(4,4), out(4,4)
      real a, b, x, y
      integer i
      
      x = a - b
      y = a + b
      
      do i = 1, 4
         out(i, 1) = x * in(i, 1)
         out(i, 2) = x * in(i, 2)
         out(i, 3) = y * in(i, 3)
         out(i, 4) = y * in(i, 4)
      end do
      end


      subroutine current(psi1, a, b, v, out)
calculate  out(mu) = psi1^bar*(a - b*gamma5)*gamma(mu)*v
c arbitrary psi1 and v
c [ with the convention of neutral & charged couplings in the book, say:
c     for W: a = 1;  b = -1
c     for Z: a = gv, b = ga
c
      implicit none
      complex psi1(4), v(4), out(4)
      real a, b
      
      real x, y
      complex u1, u2, u3, u4, I
      parameter (I = (0., 1.))
      
      x = a + b
      y = a - b
      
      u1 = conjg(psi1(3)) * x
      u2 = conjg(psi1(4)) * x
      u3 = conjg(psi1(1)) * y
      u4 = conjg(psi1(2)) * y
      
      out(1) =   u1*v(4) + u2*v(3) - u3*v(2) - u4*v(1)
      out(2) = (-u1*v(4) + u2*v(3) + u3*v(2) - u4*v(1)) * I
      out(3) =   u1*v(3) - u2*v(4) - u3*v(1) + u4*v(2)
      out(4) =   u1*v(3) + u2*v(4) + u3*v(1) + u4*v(2)
      
      end

      subroutine current_g(psi1, v, out)
calculate  out(mu) = psi1^bar*gamma(mu)*v; with arbitrary psi1 and v; 
c this is a simplified code of "current" with a=1, b=0
c --- 4by4 form of my old "Sig4" ---
      implicit none
      complex psi1(4), v(4), out(4)
      
      complex I, u1, u2, u3, u4
      parameter (I = (0., 1.))
      
      u1 = conjg(psi1(3))
      u2 = conjg(psi1(4))
      u3 = conjg(psi1(1))
      u4 = conjg(psi1(2))
      
      out(1) =   u1*v(4) + u2*v(3) - u3*v(2) - u4*v(1)
      out(2) = (-u1*v(4) + u2*v(3) + u3*v(2) - u4*v(1)) * I
      out(3) =   u1*v(3) - u2*v(4) - u3*v(1) + u4*v(2)
      out(4) =   u1*v(3) + u2*v(4) + u3*v(1) + u4*v(2)
      end


      subroutine generate_qcd
calculate all QCD color factors T's and f_abc's numerically
      implicit none
      
      complex*16 ZI
      parameter( ZI = (0.d0, 1.d0))
      
      common /qcd_matrices/ T, f
      complex*16 T(8, 3, 3)
      real*8 f(8, 8, 8)
      
      integer i, j, k
      
      do i = 1, 8
         do j = 1, 3
            do k = 1, 3
               T(i, j, k) = dcmplx(0.)
            end do
         end do
         do j = 1, 8
            do k = 1, 8
               f(i, j, k) = 0.d0
            end do
         end do
      end do
      
      T(1, 1, 2) = 0.5d0
      T(1, 2, 1) = 0.5d0
      
      T(2, 1, 2) = - ZI / 2.d0
      T(2, 2, 1) =   ZI / 2.d0
      
      T(3, 1, 1) =  0.5d0
      T(3, 2, 2) = -0.5d0
      
      T(4, 1, 3) = 0.5d0
      T(4, 3, 1) = 0.5d0
      
      T(5, 1, 3) = - ZI / 2.d0
      T(5, 3, 1) =   ZI / 2.d0
      
      T(6, 2, 3) = 0.5d0
      T(6, 3, 2) = 0.5d0
      
      T(7, 2, 3) = - ZI / 2.d0
      T(7, 3, 2) =   ZI / 2.d0
      
      T(8, 1, 1) =   1.d0 / sqrt(3.d0) / 2.d0
      T(8, 2, 2) =   1.d0 / sqrt(3.d0) / 2.d0
      T(8, 3, 3) = - 2.d0 / sqrt(3.d0) / 2.d0
      
      f(1, 2, 3) =  1.d0
      f(1, 3, 2) = -1.d0
      f(1, 4, 7) =  0.5d0
      f(1, 7, 4) = -0.5d0
      f(1, 5, 6) = -0.5d0
      f(1, 6, 5) =  0.5d0
      
      f(2, 1, 3) = -1.d0
      f(2, 3, 1) =  1.d0
      f(2, 4, 6) =  0.5d0
      f(2, 6, 4) = -0.5d0
      f(2, 5, 7) =  0.5d0
      f(2, 7, 5) = -0.5d0
      
      f(3, 2, 1) = -1.d0
      f(3, 1, 2) =  1.d0
      f(3, 4, 5) =  0.5d0
      f(3, 5, 4) = -0.5d0
      f(3, 6, 7) = -0.5d0
      f(3, 7, 6) =  0.5d0
      
      f(4, 1, 7) = -0.5d0
      f(4, 7, 1) =  0.5d0
      f(4, 2, 6) = -0.5d0
      f(4, 6, 2) =  0.5d0
      f(4, 3, 5) = -0.5d0
      f(4, 5, 3) =  0.5d0
      f(4, 5, 8) =  sqrt(3.d0) / 2.d0
      f(4, 8, 5) = -sqrt(3.d0) / 2.d0
      
      f(5, 2, 7) = -0.5d0
      f(5, 7, 2) =  0.5d0
      f(5, 4, 3) = -0.5d0
      f(5, 3, 4) =  0.5d0
      f(5, 1, 6) =  0.5d0
      f(5, 6, 1) = -0.5d0
      f(5, 4, 8) = -sqrt(3.d0) / 2.d0
      f(5, 8, 4) =  sqrt(3.d0) / 2.d0
      
      f(6, 4, 2) = -0.5d0
      f(6, 2, 4) =  0.5d0
      f(6, 1, 5) = -0.5d0
      f(6, 5, 1) =  0.5d0
      f(6, 3, 7) =  0.5d0
      f(6, 7, 3) = -0.5d0
      f(6, 7, 8) =  sqrt(3.d0) / 2.d0 
      f(6, 8, 7) = -sqrt(3.d0) / 2.d0
      
      f(7, 4, 1) = -0.5d0
      f(7, 1, 4) =  0.5d0
      f(7, 5, 2) = -0.5d0
      f(7, 2, 5) =  0.5d0
      f(7, 3, 6) = -0.5d0
      f(7, 6, 3) =  0.5d0
      f(7, 6, 8) = -sqrt(3.d0) / 2.d0
      f(7, 8, 6) =  sqrt(3.d0) / 2.d0
      
      f(8, 5, 4) = -sqrt(3.d0) / 2.d0
      f(8, 4, 5) =  sqrt(3.d0) / 2.d0
      f(8, 7, 6) = -sqrt(3.d0) / 2.d0
      f(8, 6, 7) =  sqrt(3.d0) / 2.d0
      
      end
ccc
      subroutine gamma3(p1, p2, eps1, eps2, gamma)
c
c this does a triple coupling factor:
c      Gamma(mu) = (p2 - p1)*dot(eps1,eps2) 
c                 + 2*p1.eps2*eps1 
c                 - 2*p2.eps1*eps2
c      note that the condition p.eps = 0 has been used 
c       (valid for on-shell vector bosons or massless fermion couplings)
      implicit none
      real*8 p1(4), p2(4), eps1(4), eps2(4), gamma(4)

      real*8 Ddot
      integer i

      do i = 1, 4
           gamma(i) = 2.*Ddot(p1, eps2)*eps1(i) -
     1                2.*Ddot(p2, eps1)*eps2(i) +
     1                (p2(i) - p1(i)) * Ddot(eps1, eps2)
      end do
c
      return
      end
C-----------------------------------------------------------------------
c    (3). Integration.f [Gauss; Sample; random#'s Subcodes]
C-----------------------------------------------------------------------
c
C     THIS IS AN N-POINT GAUSS INTEGRATION SUB-PROGRAM 
C     UPTO-DATE 07/25/1991
c       N=8,10,14,20,28,40 OR 64
c
      subroutine gaussx(n,a,b,h,nm,ux,wx)
      implicit real*8(a-h,o-z)
      dimension ux(n),wx(n),u(1055),w(1055)    
C------------------------------------------------------------------------
C     FOR  8 GAUSSIAN POINTS:
      data (u(i),i=16,19)/.1834346,.5255324,.7966665,.9602899/
      data (w(i),i=16,19)/.3626838,.3137066,.2223810,.1012285/
C     FOR 10 GAUSSIAN POINTS:
       data (u(i),i=25,29)/.1488743,.4333954,.6794096,.8650634,.9739065/
       data (w(i),i=25,29)/.2955242,.2692667,.2190864,.1494514,.0666713/
C     FOR 14 GAUSSIAN POINTS:
      DATA (u(i),i=49,55)/.1080550D0,.3191124D0,.5152486D0,
     1           .6872929D0,.8272013D0,.9284349D0,.9862838D0/
      DATA (w(i),i=49,55)/.2152639D0,.2051985D0,.1855384D0,
     1          .1572032D0,.1215186D0,.8015809D-1,.3511946D-1/
C     FOR 20 GAUSSIAN POINTS:
      DATA (u(i),i=100,109)/.0765265D0,.2277859D0,.3737061D0,
     1                        .5108670D0,.6360537D0,.7463319D0,
     2             .8391170D0,.9122344D0,.9639719D0,.9931286D0/
      DATA (w(i),i=100,109)/.1527534D0,.1491730D0,.1420961D0,
     1                            .1316886D0,.1181945D0,.1019301D0,
     2          .8327674D-1,.6267205D-1,.4060143D-1,.1761401D-1/
C     FOR 28 GAUSSIAN POINTS:
      DATA (u(i),i=196,209)/.5507929D-1,.1645693D0,
     1       .2720616D0,.3762515D0,.4758742D0,.5697205D0,
     2       .6566511D0,.7356109D0,.8056414D0,.8658925D0,
     3       .9156330D0,.9542593D0,.9813032D0,.9964423D0/
      DATA (w(i),i=196,209)/.1100470D0,.1087112D0,
     1 .1060558D0,.1021130D0,.9693066D-1,.9057174D-1,
     2 .8311342D-1,.7464621D-1,.6527292D-1,.5510735D-1,
     3 .4427294D-1,.3290143D-1,.2113211D-1,.9124283D-2/
C     FOR 40 GAUSSIAN POINTS:
      DATA (u(i),i=400,419)/
     1  .3877242D-1,.1160841D0,.1926976D0,.2681522D0,.3419941D0,
     2  .4137792D0,.4830758D0,.5494671D0,.6125539D0,.6719567D0,
     3  .7273183D0,.7783057D0,.8246122D0,.8659595D0,.9020988D0,
     4  .9328128D0,.9579168D0,.9772599D0,.9907262D0,.9982377D0/
      DATA (w(i),i=400,419)/
     1  .7750595D-1,.7703982D-1,.7611036D-1,.7472317D-1,.7288658D-1,
     2  .7061165D-1,.6791205D-1,.6480401D-1,.6130624D-1,.5743977D-1,
     3  .5322785D-1,.4869581D-1,.4387091D-1,.3878217D-1,.3346020D-1,
     4  .2793701D-1,.2224585D-1,.1642106D-1,.1049828D-1,.4521277D-2/
C     FOR 64 GAUSSIAN POINTS:
      DATA (u(i),i=1024,1055)/.2435029D-1,.7299312D-1,
     1  .1214628D0,.1696444D0,.2174236D0,.2646872D0,.3113229D0,
     2  .3572202D0,.4022702D0,.4463660D0,.4894031D0,.5312795D0,
     3  .5718956D0,.6111554D0,.6489655D0,.6852363D0,.7198819D0,
     4  .7528199D0,.7839724D0,.8132653D0,.8406293D0,.8659994D0,
     5  .8893154D0,.9105221D0,.9295692D0,.9464114D0,.9610088D0,
     6  .9733268D0,.9833363D0,.9910134D0,.9963401D0,.9993050D0/
      DATA (w(i),i=1024,1055)/.4869096D-1,.4857547D-1,
     1  .4834476D-1,.4799939D-1,.4754017D-1,.4696818D-1,.4628480D-1,
     2  .4549163D-1,.4459056D-1,.4358372D-1,.4247352D-1,.4126256D-1,
     3  .3995374D-1,.3855015D-1,.3705513D-1,.3547221D-1,.3380516D-1,
     4  .3205793D-1,.3023466D-1,.2833967D-1,.2637747D-1,.2435270D-1,
     5  .2227017D-1,.2013482D-1,.1795172D-1,.1572603D-1,.1346305D-1,
     6  .1116814D-1,.8846760D-2,.6504458D-2,.4147033D-2,.1783281D-2/
C--------------------------------------------------------------------
      m=n/2
      jf=m*m
      jl=n-2*m
      l=jl*(jf-1)
      jl=jf+m-1
      m=abs(b-a)/h+.1
      IF(m .LE. 0) THEN
         GO TO 4   
      ELSE                             
         GO TO 2 
      END IF
 4      m=1
 2      hd=(b-a)/float(m)
      hd2=hd/2.
      x=a-hd2
      do 20 i=1,m
      k0=n*(i-1)+(n+1)/2
      x=x+hd
      do 20 j=jf,jl
      t=u(j)*hd2
      jj=j-jf+1
      k1=k0+jj
      k2=k0-jj
      if(l.eq.0) k2=k2+1
      ux(k1)=x+t
      ux(k2)=x-t
      wx(k1)=w(j)*hd2
      wx(k2)=w(j)*hd2
 20      continue
      nm=n*m
      return
      end
C-----------------------------------------------------------------------
c
      SUBROUTINE SAMPLE(X,WGT,NDIM,MANY,ITN,KN,NG,DG)
c-----------------------------------------------
c  locally modified by Alan Stange on Sept. 15 to return
c  the mean, sigma and chi2 in a common block for use 
c  in the user program
c--------------------------------------------------
C-----X IS THE SAMPLE PT. IN UNIT HYPERCUBE-----
C-----WGT(INPUT) IS PREVIOUS IMPORTANCE
C                   INFORMATION-----------------
C-----WGT(OUTPUT) IS   EQUAL    TO THE VOLUME
C               OCCUPIED BY EACH SAMPLE --------
C-----NDIM IS DIMENSION OF THE HYPERCUBE--------
C-----NG IS GRID NUMBER ON EACH AXIS------------
C-----DG(1,.,.) IS HISTOGRAM OF IMPORTANCE  IN-
C     FORMATION ON EACH AXIS OF HYPERCUBE-------
C-----D(2,.,.) IS GRID POSITION-----------------
      DIMENSION DG(2,NG,NDIM),X(NDIM)

      common /samplexx_common/ samplexx_mean, samplexx_sigma, 
     1        samplexx_chi2
      real*4 samplexx_mean, samplexx_sigma, samplexx_chi2

      REAL      DUM(50),MEAN
      INTEGER I
      character*4 GRID 
c
      COMMON/SCTRL/TMEAN,TSIGMA,IREAD,IWRITE
      DATA ID/0/,IREAD/1/,IWRITE/1/,MEAN/0./,SIGMA/0./
  50  FORMAT(A4)
  60  FORMAT('GRID')
  70  FORMAT(17HIITERATION NUMBER,I3,3X,6HMEAN--,G10.4,
     +          16H,  FLUCTUATION--,G10.4)
  80  FORMAT(/1X,79(1H-)/1X,23HACCUMULATED RESULTS:   ,
     ,10HINTEGRAL =,G10.4/24X,10HSTD DEV  =,G10.4/
     1      13X,21HCHI**2 PER DOF.     =,G10.4/1X,79(1H-))
 100  FORMAT(10F8.4)
C
C
C-----INITIATION OF GRIDS-----------------------
C
      KN=KN+1
c      IF(ID)  19,8,9
      IF(ID .LT. 0)  GOTO 19
      IF(ID .GT. 0)  GOTO 9
   8  IF(ITN .EQ. 0) THEN
        ID=-1
      ELSE
        ID=1
        KNT=(MANY-1)/IABS(ITN)
      ENDIF
   3  XNG=NG
      NGM=NG-1
      IT=1
      KN=1
      VOL=1.0/MANY
      TMEAN=0.0
      TSIGMA=0.0
      CHI2=0.0
C
c      IF (IREAD) 19,10,4
      IF (IREAD .EQ. 0) GOTO 10
      IF (IREAD .LT. 0) GOTO 19
4      open(unit=25 ,err=999, status='old')      
      READ(25,50)   GRID
      print *,'fort.25 as ', GRID
      IF (GRID .NE. 'GRID') GO TO 998
      READ(25,100) ((DG(2,I,J),I=1,NG),J=1,NDIM)
            GO TO 19
999     print *,' NO unit 25 file: Generate Uniform Grid'
            goto 10
998      print *,' Empty unit 25 file: Generate Uniform Grid'
10      CONTINUE
      DO 12 I=1,NG
      DO 11 J=1,NDIM
  11  DG(2,I,J)=I/XNG
  12  CONTINUE
      GO TO 19
C
C-----DUMP  PREVIOUS INPORTANCE INFORMATION
C       IN HISTOGRAM----------------------------
C
   9  CONTINUE
      MEAN=MEAN+WGT
      SIGMA=SIGMA+WGT*WGT
      DO 16 J=1,NDIM
      IP=IFIX(DUM(J))+1
  16  DG(1,IP,J)=DG(1,IP,J)+ABS(WGT)
C
      IF(KN.LT.KNT*IT) GO TO 19
      MEAN=MEAN/VOL/FLOAT(KNT)
C     SIGMA=(SIGMA/VOL/VOL-KNT*MEAN*MEAN)/FLOAT((KNT-1)*KNT)
      SIGMA=(SIGMA/VOL/VOL-KNT*MEAN*MEAN)/FLOAT(KNT-1)/FLOAT(KNT)
      TMEAN=TMEAN+MEAN*(MEAN*MEAN/SIGMA)
      TSIGMA=TSIGMA+(MEAN*MEAN/SIGMA)
      CHI2=CHI2+MEAN*MEAN*(MEAN*MEAN/SIGMA)
      SIGMA=SQRT(SIGMA)
      IF (IWRITE.LT.0)GO TO 17
      WRITE(6,70)  IT,MEAN,SIGMA
17    MEAN=0.0
      SIGMA=0.0
      IF(ITN.LT.0) GO TO 29
C-----REFINE GRIDS------------------------------
C
      DO 28 J=1,NDIM
      XO=DG(1,1,J)
      XN=DG(1,2,J)
      DG(1,1,J)=(XO+XN)/2.0
      X(J)=DG(1,1,J)
      DO 22 I=2,NGM
      DG(1,I,J)=XO+XN
      XO=XN
      IP=I+1
      XN=DG(1,IP,J)
      DG(1,I,J)=(DG(1,I,J)+XN)/3.
  22  X(J)=X(J)+DG(1,I,J)
      DG(1,NG,J)=(XN+XO)/2.0
      X(J)=X(J)+DG(1,NG,J)
C
      RC=0.0
  
      DO 24 I=1,NG
      XO=(1.0E-30)+DG(1,I,J)/X(J)
      DG(1,I,J)=((XO-1.0)/ALOG(XO))**1.5
  24  RC=RC+DG(1,I,J)
  
      RC=RC/XNG
      K=0
      XN=0.
      DR=0.0
      I=0

  25  K=K+1
      DR=DR+DG(1,K,J)
      XO=XN
      XN=DG(2,K,J)
  26  IF(RC.GT.DR) GO TO 25

      I=I+1
      DR=DR-RC
      DUM(I)=XN-(XN-XO)*DR/DG(1,K,J)
      IF(I.LT.NGM) GO TO 26
      DO 27 I=1,NGM
      DG(1,I,J)=0.0
  27  DG(2,I,J)=DUM(I)
      DG(1,NG,J)=0.0
  28  DG(2,NG,J)=1.0
  29  IT=IT+1
      IF(IT.LE.IABS(ITN)) GO TO 19
      TMEAN=TMEAN/TSIGMA
      CHI2=(CHI2/TMEAN/TMEAN-TSIGMA)/FLOAT((IABS(ITN)-1))**2
     ,    *(IABS(ITN)-1)
      TSIGMA=TMEAN/SQRT(TSIGMA)
c      IF (IWRITE) 32,31,30
      IF (IWRITE .GT. 0) THEN
  30    WRITE(26,60)
        WRITE(26,100) ((DG(2,I,J),I=1,NG),J=1,NDIM)
      ENDIF
      IF (IWRITE .GE. 0) THEN
  31    WRITE(6,80) TMEAN,TSIGMA,CHI2
c=====================================================
c   change made to return the results in a common block
      samplexx_mean = tmean
      samplexx_sigma = tsigma
      samplexx_chi2 = chi2
c=====================================================
      ENDIF
  32  ID=-1
C
C-----LOCATE SAMPLE PT. IN THE HYPERCUBE--------
  19  CONTINUE
      WGT=VOL
      DO 15 J=1,NDIM
c      CALL HTUPLE(DUM(J),0.0,XNG,J)
      CALL htuple(DUM(J),0.0,XNG,J)
      IM=DUM(J)
      IP=IM+1
      IF(IP.GT.1) GO TO 13
      XO=DG(2,IP,J)
      GO TO 14
  13  XO=DG(2,IP,J)-DG(2,IM,J)
  14  X(J)=DG(2,IP,J)-XO*(IP-DUM(J))
      WGT=WGT*XO*XNG
  15  CONTINUE
      IF(KN.GE.MANY) ID=0
      RETURN
      END
c --------------------------------------------------------------
CCCC  5/3/88 FOR REPLACING HTUPLE.FOR :
      subroutine htuple(x,a,b,i)

      integer ndim,maxj,mxbase
      parameter (ndim = 25, maxj = 18, mxbase = 97)
      integer base(ndim),mix(mxbase,ndim),k(ndim,maxj),l(ndim)
      real*4 accum(ndim,maxj)
      real*8 pbase(ndim,maxj)

      data pbase/2.d0,3.d0,5.d0,7.d0,11.d0,13.d0,17d0,19.d0,23.d0,29.d0,
     1       31.d0,37.d0,41.d0,43.d0,47.d0,53.d0,59.d0,61.d0,67d0,71.d0,
     2       73.d0,79.d0,83.d0,89.d0,97.d0,425*0.d0/
      data l/25*1/, base/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,
     1       53,59,61,67,71,73,79,83,89,97/

      data mix/1,96*0,
     * 1,2,95*0,
     * 3,1,4,2,93*0,
     * 4,2,6,1,5,3,91*0,
     * 5,8,2,10,3,6,1,9,7,4,87*0,
     * 6,10,2,8,4,12,1,9,5,11,3,7,85*0,
     * 8,13,3,11,5,16,1,10,7,14,4,12,2,15,6,9,81*0,
     * 9,14,3,17,6,11,1,15,7,12,4,18,8,2,16,10,5,13,79*0,
     * 11,17,4,20,7,13,2,22,9,15,5,18,1,14,10,21,6,16,3,19,8,12,75*0,
     * 15,7,24,11,20,2,27,9,18,4,22,13,26,5,16,10,23,1,19,28,6,14,17,3,
     +       25,12,8,21,69*0,
     * 15,23,5,27,9,18,2,29,12,20,7,25,11,17,3,30,14,22,1,21,8,26,10,16,
     +       28,4,19,6,24,13,67*0,
     * 18,28,6,23,11,34,3,25,14,31,8,20,36,1,16,27,10,22,13,32,4,29,17,
     +       7,35,19,2,26,12,30,9,24,15,33,5,21,61*0,
     * 20,31,7,26,12,38,3,23,34,14,17,29,5,40,10,24,1,35,18,28,9,33,15,
     +           21,4,37,13,30,8,39,22,2,27,16,32,11,25,6,36,19,57*0,
     * 21,32,7,38,13,25,3,35,17,28,10,41,5,23,30,15,37,1,19,33,11,26,
     +       42,8,18,29,4,39,14,22,34,6,24,12,40,2,31,20,27,9,36,16,
     1       55*0,
     * 24,12,39,6,33,20,44,3,29,16,36,10,42,22,8,31,26,14,46,1,35,18,28,
     +       5,40,19,37,11,25,43,4,30,15,34,9,45,21,2,32,17,41,13,27,
     +       7,38,23,51*0,
     * 26,40,9,33,16,49,4,36,21,45,12,29,6,51,23,38,14,43,1,30,19,47,10,
     2       34,24,      
     +       42,3,27,52,15,18,39,7,46,31,11,35,20,48,2,28,41,8,22,50,13,
     3       32,17,
     +       44,5,37,25,45*0,
     * 29,44,10,52,18,34,4,48,23,38,13,57,7,32,41,20,54,2,26,46,15,36,
     4       24,50,8,
     +       40,16,56,5,30,43,21,51,11,33,1,58,27,37,14,47,19,28,45,6,
     5       53,12,
     +       35,22,42,3,55,25,31,9,49,17,39,39*0,
     * 30,46,10,38,18,56,4,42,24,52,14,33,21,59,6,40,27,49,2,35,16,54,
     6       12,44,
     +       26,50,8,32,58,19,1,41,29,48,13,36,22,60,7,45,23,53,9,34, 
     7       17,55,
     +       3,39,28,47,15,37,20,57,5,43,25,51,11,31,37*0,
     * 33,50,11,59,20,39,5,54,26,44,15,64,23,36,2,57,30,47,9,62,18,41,
     8       13,52,
     +       28,37,4,66,24,46,8,55,31,17,60,34,1,48,21,43,63,12,38,25,
     9       53,7,
     +       49,16,58,29,6,42,65,19,35,10,51,27,56,3,40,32,61,14,45,  
     1       22,31*0,
     * 35,53,12,62,21,41,5,67,28,46,16,56,25,8,50,38,65,2,32,59,19,44,
     2       14,
     +       70,30,48,7,39,58,22,10,63,33,26,52,1,55,18,43,68,13,36,
     3       47,4,
     +       61,24,40,29,66,9,51,17,57,23,37,3,69,31,45,15,60,11,49,
     4       34,20,
     +       64,6,54,27,42,27*0,
     * 36,55,12,46,22,67,5,41,61,18,30,52,8,70,27,43,15,59,33,2,64,38,
     5       24,
     +       50,10,72,20,48,31,57,4,63,25,40,14,54,35,68,7,45,17,60,28,
     +       1,66,39,21,51,11,71,32,47,13,56,26,44,3,65,34,19,
     +       58,9,49,37,69,16,29,53,6,62,23,42,25*0,
     * 39,59,13,69,24,46,6,74,31,51,18,63,9,42,55,27,77,2,35,65,21,48, 
     6       15,
     +       71,33,53,4,61,29,43,17,75,37,10,67,49,22,57,7,72,26,40,56,
     +       1,64,30,45,14,78,20,52,34,11,68,41,60,5,36,73,23,50,16,62,
     +       28,3,76,44,25,58,12,66,38,19,54,32,70,8,47,19*0,
     * 41,62,14,73,25,48,6,67,32,54,19,80,10,44,58,29,76,2,37,64,
     +       22,51,16,71,35,56,8,82,27,46,12,69,39,60,4,50,24,78,31,65,       
     +       17,42,74,1,53,21,61,34,11,79,43,28,68,7,55,38,75,15,47,20,
     +       70,5,57,33,81,26,49,9,63,36,66,18,45,3,77,30,59,23,52,13,      
     +       72,40,15*0,
     * 44,67,15,56,27,82,6,50,74,22,36,63,10,86,33,53,18,77,40,2,
     +       70,47,29,80,12,60,38,65,20,88,4,51,31,72,24,58,8,78,42,46,
     +       16,84,34,62,1,69,26,55,19,76,41,11,83,49,30,66,7,59,37,87,
     +       14,54,25,73,21,68,43,3,79,35,57,13,81,45,28,64,5,75,32,52,
     +       17,85,39,9,61,71,23,48,9*0,
     * 48,73,16,61,29,89,7,55,81,34,22,69,41,94,3,52,77,19,38,85,
     +       12,64,44,26,91,58,9,71,32,79,14,50,66,24,96,1,46,83,36,59,
     +       18,75,30,87,5,54,42,68,21,92,10,62,39,80,27,56,6,86,47,72,
     +       15,35,93,43,65,2,76,25,53,84,17,37,67,11,90,49,31,74,20,60,
     +       95,4,45,63,28,82,13,57,40,78,8,88,33,51,23,70,0/
c
c
      j = 1
20      k(i,j)=k(i,j)+1
      if (k(i,j) .lt. base(i) ) then
            accum(i,j)=accum(i,j+1) + mix(k(i,j),i)/pbase(i,j)
            do  jj=2,j-1
                  accum(i,jj) = accum(i,j)
                  end do
            x = a + (b-a) * accum(i,j)
            return
      else
            k(i,j)=0
            j=j+1
            if (j.gt.l(i)) then
                  if (j.ge.maxj) then
                        do j=1,maxj
                              k(i,j)=0
                              accum(i,j)=0.
                              end do
                        l(i)=1
                        j=1
                  else
                        l(i)=j
                        pbase(i,j)=pbase(i,j-1)*pbase(i,1)
                        end if
                  end if
            go to 20
            end if
c--------
      entry htuple_reset
change & cked on 2/1/94, [do i=1,] --> [do ii=1,]: i has been used as an input
      do ii=1,ndim
            l(ii) = 1
            end do
      do ii=2,ndim*maxj
            k(ii,1)=0
            accum(ii,1)=0.
            end do
      return
      end
C-----------------------------------------------------------------------
c
c Uniformly Distributed Random Number Generator 
c Output:
c     Double Precision RN from 0 to 1:
c Input (changed by Han 2/13/92):
c     Iseed --- random seed for initialization
c     the range should be 0<= Iseed <=31328
c
      FUNCTION RN(ISEED)
      REAL*8 RN,RAN
      INTEGER ISEED
      SAVE INIT
      DATA INIT /1/
      IF (INIT.EQ.1) THEN
        INIT=0
        CALL RMARIN(ISEED,23456)
      END IF
*
  10  CALL RANMAR(RAN)
      IF (RAN.LT.1D-16) GOTO 10
      RN=RAN
*
      END
*
      SUBROUTINE RANMAR(RVEC)
*     -----------------
* Universal random number generator proposed by Marsaglia and Zaman
* in report FSU-SCRI-87-50
* In this version RVEC is a double precision variable.
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
      COMMON/ RASET2 / IRANMR,JRANMR
      SAVE /RASET1/,/RASET2/
      UNI = RANU(IRANMR) - RANU(JRANMR)
      IF(UNI .LT. 0D0) UNI = UNI + 1D0
      RANU(IRANMR) = UNI
      IRANMR = IRANMR - 1
      JRANMR = JRANMR - 1
      IF(IRANMR .EQ. 0) IRANMR = 97
      IF(JRANMR .EQ. 0) JRANMR = 97
      RANC = RANC - RANCD
      IF(RANC .LT. 0D0) RANC = RANC + RANCM
      UNI = UNI - RANC
      IF(UNI .LT. 0D0) UNI = UNI + 1D0
      RVEC = UNI
      END
 
      SUBROUTINE RMARIN(IJ,KL)
*     -----------------
* Initializing routine for RANMAR, must be called before generating
* any pseudorandom numbers with RANMAR. The input values should be in
* the ranges 0<=ij<=31328 ; 0<=kl<=30081
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
      COMMON/ RASET2 / IRANMR,JRANMR
      SAVE /RASET1/,/RASET2/
* This shows correspondence between the simplified input seeds IJ, KL
* and the original Marsaglia-Zaman seeds I,J,K,L.
* To get the standard values in the Marsaglia-Zaman paper (i=12,j=34
* k=56,l=78) put ij=1802, kl=9373
      I = MOD( IJ/177 , 177 ) + 2
      J = MOD( IJ     , 177 ) + 2
      K = MOD( KL/169 , 178 ) + 1
      L = MOD( KL     , 169 )
      DO 300 II = 1 , 97
        S =  0D0
        T = .5D0
        DO 200 JJ = 1 , 24
          M = MOD( MOD(I*J,179)*K , 179 )
          I = J
          J = K
          K = M
          L = MOD( 53*L+1 , 169 )
          IF(MOD(L*M,64) .GE. 32) S = S + T
          T = .5D0*T
  200   CONTINUE
        RANU(II) = S
  300 CONTINUE
      RANC  =   362436D0 / 16777216D0
      RANCD =  7654321D0 / 16777216D0
      RANCM = 16777213D0 / 16777216D0
      IRANMR = 97
      JRANMR = 33
      END
C-----------------------------------------------------------------------
c    (4). HBook Related:
C-----------------------------------------------------------------------
C
C      Routine to initialize a one-independent-variable histogram
C
      subroutine hbook1(id,inlabel,nx,xmin,xmax,zinitial)
C
C            id = integer used to identify histogram to HFILL
C            inlabel = label to be written on the output by the plotting
C                        program (character of len <=40)
C            nx = number of x bins (integer)
C            xmin = min x value (real)
C            xmax = max x value (real)
C            zinitial = initial value for each bin (real)
C
ccc      include 'hbook.inc'
C            Internal common blocks for pheno/hbook routines
C
C            LABELS(i) = label for histogram i
C            nhist        = number of histograms (starts as 0, max is 200)
C            idnumber(i) = code number to identify histograms in HFILL
C            pointer(i)  = index of beginning of data for histo # i
C            single dim(i) = .true. if single variable, .false. if double
C
      parameter (nhistmax=200,nhistmax1=nhistmax+1)
      real data(50000)
      integer pointer(nhistmax1),id number(nhistmax)
      logical single dim(nhistmax)
      character*40 label(nhistmax)

      common /hbooklabel/ label
c      common /hbookarrays/ nhist,id number, pointer, single dim, data
      common /hbookarrays/ nhist,id number, pointer, single dim
      common /hbookdata/ data
ccc
      character*(*) inlabel

      if (nhist .eq. nhistmax) then
            write(*,*)' Maximum number of histograms exceeded'
      else
            nhist = nhist+1
            label(nhist) = inlabel
            id number(nhist) = id
            single dim(nhist) = .true.
            k=pointer(nhist)
            pointer(nhist+1) = nx+3+k
            data(k)=nx
            data(k+1)=xmin
            data(k+2)=xmax
            do i=k+3,pointer(nhist+1)-1
                  data(i)=zinitial
                  end do
            end if
      return
      end
C-----------------------------------------------------------------------
C
C            Routine to dump histogram data to a file
C            hcurvesc is the same as hcurve except that
C            it scales the vertical axis by a factor scale.
      subroutine hcurvesc(id,filename,scale)      
C
C            Dumps all current histograms to file 'filename' and
C            clears all current histograms.
C            scale = vertical axis scale factor (input)
C
ccc      include 'hbook.inc'
C            Internal common blocks for pheno/hbook routines
C
C            LABELS(i) = label for histogram i
C            nhist        = number of histograms (starts as 0, max is 200)
C            idnumber(i) = code number to identify histograms in HFILL
C            pointer(i)  = index of beginning of data for histo # i
C            single dim(i) = .true. if single variable, .false. if double
C
      parameter (nhistmax=200,nhistmax1=nhistmax+1)
      real data(50000)
      integer pointer(nhistmax1),id number(nhistmax)
      logical single dim(nhistmax)
      character*40 label(nhistmax)

      common /hbooklabel/ label
c      common /hbookarrays/ nhist,id number, pointer, single dim, data
      common /hbookarrays/ nhist,id number, pointer, single dim
      common /hbookdata/ data
c
      character*(*) filename

      if (nhist .eq. 0) return
      open(unit=69,FILE=filename,status='new')
c,carriagecontrol='list')
c     1        carriagecontrol='list')
      do i = 1, nhist
            if (id .eq. idnumber(i)) go to 10
            end do
      return
10      continue
      k = pointer(i)
      nx = int(data(k)+.1)
      xmin = data(k+1)
      xmax = data(k+2)
      xbinsize = (xmax-xmin)/nx
      if (single dim(i)) then
            write (69,300) label(i)(1:labelleng(label(i)))
            write (69,400) (xmin+(m-.5)*xbinsize,data(k+2+m)*
     1                          scale,m=1,nx)
      else
            ny = int(data(k+3) + .1)
            ymin = data(k+4)
            ymax = data(k+5)
            ybinsize = (ymax-ymin)/ny
            write (69,300) label(i)(1:labelleng(label(i)))
            k = k + 5
            do n=1,ny
                  fixed y =  ymin + (n-.5)*ybinsize
                  write (69,500) 
     1    (xmin+(m-.5)*xbinsize,fixed y,data(k+m)*scale,m=1,nx)
                  k = k + nx
                  end do
            end if
      close (unit=69)
      return
300      format (' Histogram ',a)
400      format (1x,2g15.6)
500      format (1x,3g15.6)
      end
C
C
      function labelleng(string)
      character*(*) string

      do i=len(string),1,-1
            if (string(i:i) .ne. ' ') then
                  labelleng=i
                  return
                  end if
            end do
      labelleng=1
      return
      end
C-----------------------------------------------------------------------
C
C            Routine to add zincrement to a bin in a histogram
C
      subroutine hfill(id,x,y,zincrement)
C
C            id = integer associated with the histogram
C            x  = x value to locate bin (real)
C            y  = y value to locate bin (real) [ignored for 1-dim histo]
C            zincrement = value to be added to bin specified by (x,y)
C
ccc      include 'hbook.inc'
C            Internal common blocks for pheno/hbook routines
C
C            LABELS(i) = label for histogram i
C            nhist        = number of histograms (starts as 0, max is 200)
C            idnumber(i) = code number to identify histograms in HFILL
C            pointer(i)  = index of beginning of data for histo # i
C            single dim(i) = .true. if single variable, .false. if double
C
      parameter (nhistmax=200,nhistmax1=nhistmax+1)
      real data(50000)
      integer pointer(nhistmax1),id number(nhistmax)
      logical single dim(nhistmax)
      character*40 label(nhistmax)

      common /hbooklabel/ label
c      common /hbookarrays/ nhist,id number, pointer, single dim, data
      common /hbookarrays/ nhist,id number, pointer, single dim
      common /hbookdata/ data
c
      data nhist/0/,pointer(1)/1/

      do i=1,nhist
            if (id number(i) .eq. id) go to 10
            end do
      write(*,*) ' id number ',id,' not belong to a current histogram'
      return
10      continue
      k = pointer(i)
      nx=data(k)+.1
      ixoff = (x-data(k+1))/(data(k+2)-data(k+1))*data(k)+1
      if (ixoff .le. 0 .or. ixoff .gt. nx) return
      if (single dim(i)) then
            data(k+2+ixoff)=data(k+2+ixoff)+zincrement
      else
            ny=data(k+3)+.1
            iyoff = (y-data(k+4))/
     1                  (data(k+5)-data(k+4))*data(k+3)+1
            if (iyoff .le. 0 .or. iyoff .gt. ny) return
            ioff = nx*(iyoff-1)+ixoff
            data(k+5+ioff)=data(k+5+ioff)+zincrement
            end if
      return
      end
C-----------------------------------------------------------------------
C      Routine to initialize a two-independent-variable histogram
C
      subroutine hbook2(id,inlabel,nx,xmin,xmax,ny,ymin,ymax,zinitial)
C
C            id = integer used to identify histogram to HFILL
C            inlabel = label to be written on the output by the plotting
C                        program (character of len <=40)
C            nx = number of x bins (integer)
C            xmin = min x value (real)
C            xmax = max x value (real)
C            ny,ymin,ymax = same for y values
C            zinitial = initial value for each bin (real)
C--------------------------
c      include 'hbook.inc':
C            Internal common blocks for pheno/hbook routines
C
C            LABELS(i) = label for histogram i
C            nhist        = number of histograms (starts as 0, max is 200)
C            idnumber(i) = code number to identify histograms in HFILL
C            pointer(i)  = index of beginning of data for histo # i
C            single dim(i) = .true. if single variable, .false. if double
C
      parameter (nhistmax=200,nhistmax1=nhistmax+1)
      real data(50000)
      integer pointer(nhistmax1),id number(nhistmax)
      logical single dim(nhistmax)
      character*40 label(nhistmax)

      common /hbooklabel/ label
c      common /hbookarrays/ nhist,id number, pointer, single dim, data
      common /hbookarrays/ nhist,id number, pointer, single dim
      common /hbookdata/ data
c--------------------------
      character*(*) inlabel

      if (nhist .eq. nhistmax) then
            write(*,*) ' Maximum number of histograms exceeded'
      else
            nhist = nhist+1
            label(nhist) = inlabel
            id number(nhist) = id
            single dim(nhist) = .false.
            k=pointer(nhist)
            pointer(nhist+1) = nx*ny+6+k
            data(k)=nx
            data(k+1)=xmin
            data(k+2)=xmax
            data(k+3)=ny
            data(k+4)=ymin
            data(k+5)=ymax
            do i=k+6,pointer(nhist+1)-1
                  data(i)=zinitial
                  end do
            end if
      return
      end
C-----------------------------------------------------------------------
c    (5). Codes from Theories & math (now only QCD)
c         It includes: Gamma-function; alphas; alpha_s
C-----------------------------------------------------------------------
C
C        SUBROUTINE GAMMA
C
C        PURPOSE
C           COMPUTES THE GAMMA FUNCTION FOR A GIVEN ARGUMENT
C
C        USAGE
C           CALL GAMMA(XX,GX,IER)
C
C        DESCRIPTION OF PARAMETERS
C           XX -THE ARGUMENT FOR THE GAMMA FUNCTION
C           GX -THE RESULTANT GAMMA FUNCTION VALUE
C           IER-RESULTANT ERROR CODE WHERE
C               IER=0  NO ERROR
C               IER=1  XX IS WITHIN .000001 OF BEING A NEGATIVE INTEGER
C               IER=2  XX GT 34.5, OVERFLOW, GX SET TO 1.0E38
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           THE RECURSION RELATION AND POLYNOMIAL APPROXIMATION
C           BY C.HASTINGS,JR., 'APPROXIMATIONS FOR DIGITAL COMPUTERS',
C           PRINCETON UNIVERSITY PRESS, 1955
C
C     ..................................................................
C
      SUBROUTINE GAMMA(XX,GX,IER)
      IF(XX-34.5 .GT. 0.0) THEN
        IER=2
        GX=1.E38
        RETURN
      ENDIF
    6 X=XX
      ERR=1.0E-6
      IER=0
      GX=1.0

      IF(X-2.0 .LE. 0) THEN
         GO TO 50
      ELSE              
         GO TO 15
      END IF
   10 IF(X-2.0 .LE. 0) THEN
         GO TO 110  
      ELSE                             
         GO TO 15 
      END IF
   15 X=X-1.0
      GX=GX*X
      GO TO 10
   50 IF(X-1.0 .LT. 0) THEN
         GO TO 60  
      ELSE IF(X-1.0 .EQ. 0) THEN
         GO TO 120 
      ELSE                             
         GO TO 110 
      END IF
C
C        SEE IF X IS NEAR NEGATIVE INTEGER OR ZERO
C
   60 IF(X-ERR .LE. 0) THEN
         GO TO 62  
      ELSE                             
         GO TO 80  
      END IF
   62 Y=FLOAT(INT(X))-X
      IF((ABS(Y)-ERR .LE. 0.0) .OR. (1.0-Y-ERR .LE. 0.0)) GOTO 130
C
C        X NOT NEAR A NEGATIVE INTEGER OR ZERO
C
   70 IF(X-1.0 .LE. 0) THEN
         GO TO 80  
      ELSE                             
         GO TO 110 
      END IF
   80   GX=GX/X
        X=X+1.0
        GO TO 70
  110 Y=X-1.0
      GY=1.0+Y*(-0.5771017+Y*(+0.9858540+Y*(-0.8764218+Y*(+0.8328212+
     1Y*(-0.5684729+Y*(+0.2548205+Y*(-0.05149930)))))))
      GX=GX*GY
  120 RETURN
  130 IER=1
      RETURN
      END
C**************************************************************************
C       Date: 2/11/1991, modified by Han on 5/3/91
C**************************************************************************
c     Input: Q-square; ALAMqcd; NF--NO. FLAVORS
c
      function alphas(qsq,alam,nf)  !ALAM-LAMBDA MASS, NF--NO. FLAVORS
      beta0=11.-2.*float(nf)/3.
      alogg=alog(qsq/(alam*alam))
      alphas=4.*3.1415927/(beta0*alogg)
      return
      end
C---------------------------------------------------------------
      FUNCTION ALPHA_S(QSQ,ALAM4,TMASS,LOOP)
C---------------------------------------------------------------
C
C      THIS FUNCTION RETURNS THE 1, 2 OR 3-LOOP VALUE OF ALPHA_S
C
C      INPUT: 
C            QSQ      = Q**2 (REAL)
C            ALAM4        = LAMBDA FOR 4 ACTIVE QUARK FLAVORS (REAL)
C            TMASS   = TOP QUARK MASS TO DETERMINE LAMBDA-6 (REAL)
C            LOOP      = NUMBER OF LOOPS FOR ALPHA_S (INTEGER = 1, 2 OR 3)
C
C       PARAMETRIZATION OF THE STRONG COUPLING CONSTANT ACCORDING TO
C      LOOP = 1, 2 : FROM THE BOOK;
C       LOOP = 3:     W. J. MARCIANO, PHYS. REV. D 29 (1984) 580.
C      NOTE : THRESHOLD AT 2*Mq
ccc     9/11/91: lambda5,6 from the matching are always based on 2-loop 
ccc              alpha_s formula so that loop=1,3 will not be realy smooth
ccc     it seems that it's too hard to smooth the 3-loop case; 
ccc     but 1-loop case should be easily modified to an even simpler form
C-----------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*4 ALPHA_S, QSQ, ALAM4, TMASS      ! INPUTS
      INTEGER LOOP
      REAL*4 PI, BMASS                  ! DATA
      REAL*4 ANF, ALAM, ALAMSQ, ALAM5, T, TT, B0, B1, B2, X, ALPHAS
      DATA PI/3.1415927/, BMASS/5.0/
C
        IF (QSQ .LT. 4.0*BMASS**2) THEN
           ANF   = 4.0
           ALAM  = ALAM4
        ELSE IF (QSQ .LT. 4.0*TMASS**2) THEN
           ANF   = 5.0
           ALAM  = ALAM4*(ALAM4/(2.0*BMASS))**(2.0/23.0)
     1           *(ALOG(4.0*BMASS**2/ALAM4**2))**(-963.0/13225.0)
        ELSE
           ANF   = 6.0
           ALAM5 = ALAM4*(ALAM4/(2.0*BMASS))**(2.0/23.0)
     1           *(ALOG(4.0*BMASS**2/ALAM4**2))**(-963.0/13225.0)
           ALAM  = ALAM5*(ALAM5/(2.0*TMASS))**(2.0/21.0)
     1           *(ALOG(4.0*TMASS**2/ALAM5**2))**(-107.0/1127.0)
        END IF
        B0       = 11.0-2.0/3.0*ANF
      ALAMSQ   = ALAM**2
      T        = ALOG(QSQ/ALAMSQ)
      IF (T .LE. 1.0) T = ALOG(4.0/ALAMSQ)
      ALPHAS   = 4*PI/B0/T
      IF (LOOP .EQ. 1) THEN                  ! 1 LOOP ALPHA_S
         ALPHA_S = ALPHAS
      ELSE IF (LOOP .EQ. 2) THEN            ! 2 LOOP ALPHA_S
           B1 = 102.0-38.0/3.0*ANF
         X  = B1/(B0**2*T)
         TT = ALOG(T)
         ALPHA_S  = ALPHAS*(1.0-X*TT)
      ELSE IF (LOOP .EQ. 3) THEN            ! 3 LOOP ALPHA_S
           B1 = 102.0-38.0/3.0*ANF
           B2 = 0.5*(2857.0-5033.0/9.0*ANF+325.0/27.0*ANF**2)
           X  = B1/(B0**2*T)
           TT = ALOG(T)
        ALPHA_S = ALPHAS*(1.0-X*TT+X**2*((TT-0.5)**2+B2*B0/B1**2-1.25))
      ELSE
       PRINT *, ' WRONG LOOP NUMBER IN ALPHA-S EVALUATION !'
       STOP
      END IF
C
      RETURN
      END
C-----------------------------------------------------------------------
c    (6). Experiment related codes:
c       a. smearing: 
c        a1. CAL energy smearing:
c            SUBROUTINE SMEAR(P,A,B)
c        a2. Tracking smearing on pT:
c            SUBROUTINE SMEAR_PT(P,A,B,sint)
c
C ----------------------------------------------------------------------
c
ccc smearing code:
c
        SUBROUTINE SMEAR(P,A,B)
c
C***************************************************************************
C
C       THIS SUBROUTNE SMEARS THE MAGNITUDE OF THE 3-MOMENTUM COMPONENT
C       OF THE 4-VECTOR P(4)
C
C       INPUT:
C               P   : 4-VECTOR TO BE SMEARED (REAL*4)
C               A,B : SMEARING PARAMETERS    (REAL*4)
c               Formula in quadrature:
c               Del = sqrt[ (a*sqrt(E))**2 + (b*E)**2 ]
c
C               COMMON/ISEED/ISEED TO KEEP AN INPUT ISEED
C
C       OUTPUT:
C               P   : SMEARED 4-VECTOR
C
C****************************************************************************
 
        COMMON/ISEED/ISEED
        REAL*4 P(4), A, B, DEL, DP
c
change on 07/28/2006:
c        DEL = a*SQRT( P(4) ) + b*P(4)
c 
        DEL = sqrt( a**2*P(4) + (b*P(4))**2 )

        CALL RGAUSS(DEL,DP)
        DO I =  1, 3
                P(I) = P(I)*(1. + DP/P(4))
        END DO
                P(4) = P(4)*(1. + DP/P(4))
 
        IF (P(4) .LT. 0D0) THEN
                DO J=1, 4
                 P(J) = 0D0
                END DO
        END IF
 
        RETURN
        END
 
c---------------------------------------------------------------------------
        SUBROUTINE SMEAR_PT(P,A,B)
c
C***************************************************************************
C
C       THIS SUBROUTNE SMEARS THE MAGNITUDE OF THE transverse MOMENTUM,
C       assigned as pT = px = P(1)
C
C       INPUT:
C               P   : 4-VECTOR TO BE SMEARED (REAL*4)
C               A,B : SMEARING PARAMETERS    (REAL*4)
c                 A needs to be in GeV-1: a_ATLAS=0.36/1000
c               Formula in quadrature:
c               Del = sqrt[ (a*pT**2)**2 + (b/sqrt-sint * pT)**2 ]
c
C               COMMON/ISEED/ISEED TO KEEP AN INPUT ISEED
C
C       OUTPUT:
C               P   : SMEARED 4-VECTOR
C
C****************************************************************************
 
        COMMON/ISEED/ISEED
        REAL*4 P(4), A, B, DEL, DP
c
        pT  = sqrt( p(1)**2 + p(2)**2 )
        pm  = sqrt(   pT**2 + p(3)**2 )
        sint= max(1e-3, pT/pm)
        b2  = b**2/sint
c
        DEL = sqrt( a**2*pT**4 + b2*pT**2 )
c
        CALL RGAUSS(DEL,DP)
        DO I =  1, 4
                P(I) = P(I)*( 1. + DP/PT )
        END DO
  
        RETURN
        END
 
C***************************************************************************
        SUBROUTINE RGAUSS(DEL,DP)
C-------YIELDS GAUSSIAN DISTRIBUTION OF RANDOM NUMBERS
        Real*8 RN
        COMMON /ISEED/ ISEED
c
 5      DP    = 3*DEL*(2.*RN(Iseed) -1.)
        GAUSS = EXP(-DP**2/2./DEL**2)
        IF (RN(Iseed) - GAUSS .gt. 0d0) goto 5
 10           RETURN
        END
