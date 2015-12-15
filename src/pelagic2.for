C     Last change:  BE   15 Dec 2015    9:05 pm
C
C     PELAGIC2.FOR                          (21 June 1992)
C
C     This module contains most of the extraneous subroutines need to run
C     the implementation trials.
C
C     The subroutines in this module:
C       RESEED  - Subroutine to reseed the random number generators.
C       RAN1    - the uniform random number generator.
C       XNORM   - Normal random number generator.
C       SPOISS  - Generate a Scaled Poisson random variable
C       POISSON - Poisson random number generator.
C       GAMMA   - Generate a Gamma random variable
C       SORT2   - Sort an array of real numbers
c       ZBRAC   - Selecting a bracketing for root-finding
C       ZBRENT  - Solve an equation
C       MINIM   - Non-linear minimization routine
C       GENBET  - BETA random number generator
C       GENMUL  - Multinomial random number generator
C
C
C ------------------------------------------------------------------------------
C
      SUBROUTINE RESEED(ISEEDS,N)
C
C     Subroutine RESEED reseeds the random number generators.
C
C      USE ARYPAR
C      USE DECLARECOM
C
      INTEGER ISEEDS(400,6),N
C
C     Set new negative values for the random number generator seeds -
C     the generator is reset whenever a negative seed is used.
      ISEED1 = ISEEDS(N,1)
      ISEED2 = ISEEDS(N,2)
      ISEED3 = ISEEDS(N,3)
      ISEED4 = ISEEDS(N,4)
      ISEED5 = ISEEDS(N,5)
      ISEED6 = ISEEDS(N,6)
C
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      FUNCTION RAN1 (ISEED,MA,INEXT,INEXTP)

C     Copyright Numerical Recipes Software
C     Reference: Press, Flannery, Teukolsky & Vetterling: Numerical Recipes,
C     CUP, Cambridge 1986 (Page 199).
C
C     Function returns a uniform random deviate between 0.0 & 1.0.
C     Set ISEED to any negative value to reinitialize the sequence.
C     Note: ISEED MUST be -ve to initialise routine. It is reset to 1 on exit.
C
      DOUBLE PRECISION MA(55),MJ,MK,MBIG,MSEED,MZ,FAC,RAN1
      INTEGER ISEED,INEXT,INEXTP,I,K,II
      PARAMETER (MBIG=4000000.D0,MSEED=1618033.D0,MZ=0.D0,FAC=1.D0/MBIG)

      IF (ISEED.LT.0.D0) THEN
        MJ = MSEED-IABS(ISEED)
        MJ = MOD(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        DO 11 I = 1,54
          II = MOD(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          IF (MK.LT.MZ) MK = MK+MBIG
          MJ = MA(II)
   11   CONTINUE
        DO 12 K = 1,4
        DO 12 I = 1,55
          MA(I) = MA(I)-MA(1+MOD(I+30,55))
          IF (MA(I).LT.MZ) MA(I) = MA(I)+MBIG
   12   CONTINUE
        INEXT = 0
        INEXTP = 31
        ISEED = 1
      ENDIF
      INEXT = INEXT+1
      IF (INEXT.EQ.56) INEXT = 1
      INEXTP = INEXTP+1
      IF (INEXTP.EQ.56) INEXTP = 1
      MJ = MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ = MJ+MBIG
      MA(INEXT) = MJ
      RAN1 = MJ*FAC
      END


C ------------------------------------------------------------------------------
C ------------------------------------------------------------------------------

      FUNCTION XNORM (SIGG,MEAN,ISEED,MA,INEXT,INEXTP)

C     Uniform to Normal conversion routine

      DOUBLE PRECISION XNORM, Z1, Z2, MEAN, SIGG, RAN1, MA(55)
      INTEGER ISEED, INEXT,INEXTP
      EXTERNAL RAN1

      Z1 = RAN1(ISEED,MA,INEXT,INEXTP)
    1 Z2 = RAN1(ISEED,MA,INEXT,INEXTP)
      IF (Z2.LE.0.D0) GO TO 1

      XNORM = SIN(6.238319D0*Z1) * SQRT(-2.D0*LOG(Z2)) * SIGG+MEAN

      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      FUNCTION POISSN (XM,RANNO)

      DOUBLE PRECISION F,SS,R,XM,RANNO,POISSN
      INTEGER I

      R = RANNO
      F = EXP(-XM)
      SS = F
      I = 0
   99 IF (SS.GE.R) GOTO 100
      I = I + 1
      F = F*XM/DBLE(I)
      SS = SS + F
      GOTO 99
  100 POISSN = DBLE(I)

      END
C
C--------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION spoiss(x,s,ISEED)
C
C     Scaled Poisson. Parameters are the mean and s.d. as well as a random
C     number seed.
C
C     Global variables
      DOUBLE PRECISION X,S,RAN1,MA(55)
      INTEGER ISEED,INEXT,INEXTP
      EXTERNAL RAN1
C
C     Local variables
      INTEGER i
      DOUBLE PRECISION p
C
      IF (S.GT.0.d0) THEN
        p = x/s + LOG(RAN1(ISEED,MA,INEXT,INEXTP))
        i = 0
100     IF (p.LE.0.d0) THEN
          spoiss = s*i
          RETURN
        ENDIF
        p = p + LOG(RAN1(ISEED,MA,INEXT,INEXTP))
        i = i + 1
        GOTO 100
      ELSE
        spoiss = x
      ENDIF
C
      RETURN
      END
C
C--------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION GAMMA(PARX,PARS,ISEED)
C
C     Approximate random gamma generator: x = mean; s= sd; ISEED is a seed
C
      EXTERNAL RAN1
      DOUBLE PRECISION PARX,PARS,X,S,RAN1,P,F,FN,MA(55)
      INTEGER I,INEXT,INEXTP,ISEED
C
      X = PARX
      S = PARS
      S = S*S/X
C
      P = 0.d0
      IF (S.GT.0.d0) THEN
        X = X/S
        I = INT(X)
        FN = AINT(X)
        IF (I.GT.0) THEN
C
C         Count down from i to zero
100       IF (I.NE.0) THEN
            F = RAN1(ISEED,MA,INEXT,INEXTP)
            IF (F.LE.0.d0) THEN
              WRITE(*,*) 'WARNING F=0 in GAMMA',F
              F = 0.00001d0
            ENDIF
            P = P + LOG(F)
            I = I - 1
            GOTO 100
          ENDIF
        ENDIF
C
C       Usually this is the case (unless x is an integer)
        IF (X.GT.FN) THEN
          F = RAN1(ISEED,MA,INEXT,INEXTP)
          IF (F.LE.0.d0) THEN
            WRITE(*,*) 'WARNING F=0 in GAMMA',F
            F = 0.00001d0
          ENDIF
          GAMMA = -S*( (X-FN)*LOG(F) + P )
        ELSE
          GAMMA = -S*P
        ENDIF

      ELSE
C
        GAMMA = X
C
      ENDIF
C
      RETURN
      END
C
C ***********************************************************************
C
      SUBROUTINE SORT2 (ARRAY, N)
      DOUBLE PRECISION ARRAY(0:4000),TEMP
      INTEGER N,K,IR,I,J
C
      IF (N .LT. 2) RETURN
      K = N/2
      IR = N - 1
C
  10  IF (K .NE. 0) THEN
         K = K - 1
         TEMP = ARRAY(K)
      ELSE
         TEMP = ARRAY(IR)
         ARRAY(IR) = ARRAY(0)
         IR = IR - 1
         IF (IR .EQ. 0) THEN
            ARRAY(0) = TEMP
            RETURN
         ENDIF
      ENDIF
      I = K
      J = K + K + 1
   20 IF (J .LE. IR) THEN
         IF (J .LT. IR .AND. ARRAY(J) .LT. ARRAY(J + 1)) J = J + 1
         IF (TEMP .LT. ARRAY(J)) THEN
            ARRAY(I) = ARRAY(J)
            I = J
            J = J + I + 1
          ELSE
            J = IR + 1
          ENDIF
          GOTO 20
       ENDIF
       ARRAY(I) = TEMP
       GOTO 10
C
       END
C
C ******************************************************************************
C
      SUBROUTINE ZBRAC(FUNC,PARS,M,X1,X2,F1,F2,SUCCES)
C
C  This subroutine is based on the Numerical Recipes Software Library
C  (Refer, W.H.Press, B.P.Flannery, S.A.Teukolsky & W.T.Vetterling 1986.  
C  Numerical Recipes: The Art of Scientific Computing.  CUP, Cambridge. 818pp.
C
C  Given a function FUNC and an initial guessed range X1 to X2, the
C  routine expands the range geometrically until a root is bracketted
C  by the return values X1 and X2 (in which case SUCCES returns as
C  .TRUE.) or until the range becomes unacceptably large (in which
C  case SUCCES returns as .FALSE.).  The parameters to FUNC are passed
C  through the vector PARS.  X1 and X2 are passed to FUNC through
C  PARS(M).  The corresponding function values F1 and F2 are returned.
C  Success is guaranteed for a function which has opposite sign for
C  sufficiently large and small arguments.
C
      EXTERNAL FUNC
      DOUBLE PRECISION FUNC,PARS(4),FACTOR,X1,X2,F1,F2
      INTEGER NTRY,M,IERR,J
      LOGICAL SUCCES
      PARAMETER (FACTOR=1.6D0,NTRY=55)
C
      IF(X1.EQ.X2) STOP '** ERROR: You have to guess an initial range'
      PARS(M) = X1
      F1 = FUNC(PARS,IERR)
      IF (IERR .NE. 0) THEN
        SUCCES = .FALSE.
        RETURN
      ENDIF
      PARS(M) = X2
      F2 = FUNC(PARS,IERR)
      IF (IERR .NE. 0) THEN
        SUCCES = .FALSE.
        RETURN
      ENDIF
      SUCCES=.TRUE.
      DO 11 J=1,NTRY
        IF(F1*F2.LT.0.D0)RETURN
        IF(ABS(F1).LT.ABS(F2))THEN
          X1=X1+FACTOR*(X1-X2)
          PARS(M) = X1
          F1 = FUNC(PARS,IERR)
          IF (IERR .NE. 0) THEN
            SUCCES = .FALSE.
            RETURN
          ENDIF
        ELSE
          X2=X2+FACTOR*(X2-X1)
          PARS(M) = X2
          F2 = FUNC(PARS,IERR)
          IF (IERR .NE. 0) THEN
            SUCCES = .FALSE.
            RETURN
          ENDIF
        ENDIF
11    CONTINUE
      SUCCES=.FALSE.
      RETURN
      END
C
C  ***********************************************************************
C
      SUBROUTINE ZBRENT(FUNC,PARS,M,X1,X2,FX1,FX2,TOL,IERR)

C  This subroutine is based on the Numerical Recipes Software Library
C  (Refer, W.H.Press, B.P.Flannery, S.A.Teukolsky & W.T.Vetterling 1986.  
C  Numerical Recipes: The Art of Scientific Computing.  CUP, Cambridge. 818pp.
C
C  Using Brent's method, find the root of a function FUNC known to
C  lie between X1 and X2 (with function values FX1 and FX2).  The
C  parameters to FUNC are passed through the vector PARS.  The root
C  is returned through PARS(M)
C
      EXTERNAL FUNC
      DOUBLE PRECISION FUNC,PARS(4),X1,X2,FX1,FX2,TOL,EPS,A,FA,B,FB,C,
     +       D,E,TOL1,XM,S,P,Q,R,FC
      INTEGER M,ITMAX,ITER,IERR
      PARAMETER (ITMAX=100,EPS=3.E-5)
C
      A=X1
      FA = FX1
      B=X2
      FB = FX2
C
      IF(FB*FA>0.D0) STOP '** ERROR: Root must be bracketed for ZBRENT'
      FC=FB
      DO 11 ITER=1,ITMAX
C
        IF(FB*FC.GT.0.D0) THEN
C  Rename A, B, C and adjust bounding interval D
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
C
C  Check for convergence
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.D0*EPS*ABS(B)+0.5D0*TOL
        XM=.5D0*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.D0)THEN
          PARS(M)=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
C  Attempt inverse quadratic interpolation
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.D0*XM*S
            Q=1.D0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.D0*XM*Q*(Q-R)-(B-A)*(R-1.D0))
            Q=(Q-1.D0)*(R-1.D0)*(S-1.D0)
          ENDIF
C
C  Check whether in bounds
          IF(P.GT.0.D0) Q=-Q
          P=ABS(P)
          IF(2.D0*P .LT. MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
C  Accept interpolation
            E=D
            D=P/Q
          ELSE
C  Interpolation failed, use bisection
            D=XM
            E=D
          ENDIF
        ELSE
C  Bounds decreasing too slowly, use bisection
          D=XM
          E=D
        ENDIF
C  Save latest best guess
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
C  Evaluate new trial root
        PARS(M) = B
        FB = FUNC(PARS,IERR)
        IF (IERR .NE. 0) RETURN
11    CONTINUE
C
      STOP '** ERROR: ZBRENT exceeding maximum iterations.'
      PARS(M)=B
      RETURN
      END
C
C  ***********************************************************************
C  ***********************************************************************
C
      SUBROUTINE MINIM(F,STEP,NOP,FUNC,NMAX,IPRINT,STOPCR,NLOOP,IQUAD,
     *    SIMP,VAR,FUNCTN,IW,IFAULT)
C
C   A PROGRAM FOR FUNCTION MINIMISATION USING THE SIMPLEX METHOD.
C
C   AS THE PROGRAM IS CURRENTLY SET UP, IT WILL DEAL WITH UP TO 55
C   PARAMETERS
C
C   FOR DETAILS, SEE NELDER & MEAD, THE COMPUTER JOURNAL,JANUARY, 1965
C
C***  PROGRAMMED BY D.E.SHAW
C***  CSIRO, DIVISION OF MATHEMATICS AND STATISTICS,
C***  P.O.BOX 218, LINDFIELD.    2070
C
C***  AMENDED BY R.W.M.WEDDERBURN
C***  ROTHAMSTED EXPERIMENTAL STATION
C***  HARPENDEN, HERTFORDSHIRE, ENGLAND.
C
C   THE PARAMETER LIST IS AS FOLLOWS -
C   F         ON ENTRY, THE STARTING VALUES OF THE PARAMETERS
C             ON EXIT, THE PARAMETER VALUES SPECIFYING THE MINIMUM
C   STEP      THE INITIAL STEP SIZES
C   NOP       THE NUMBER OF PARAMETERS(INCLUDING THOSE HELD FIXED)
C   FUNC      ON EXIT, THE FUNCTION VALUE AT THE MINIMUM
C   NMAX      THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS TO BE ALLOWED
C   IPRINT    PARAMETER TO CONTROL OUTPUT FROM MINIM
C             LESS THAN 0   NO OUTPUT
C             EQUAL TO 0    REPORTS OF INITIAL EVIDENCE OF CONVERGENCE,
C                           FINAL CONVERGENCE ,WITH PARAMETER AND
C                           FUNCTION VALUES THERE), OVERRUNNING OF NMAX,
C                           AND THE FITTING OF THE QUADRATIC SURFACE IF
C                            CALLED FOR.
C             GREATER THAN 0 AS FOR IPRINT = 0, PLUS A PROGRESS REPORT
C                            ON THE MINIMISATION EVERY IPRINT FUNCTION
C                            EVALUATIONS.
C   STOPCR    STOPPING CRITERION
C   NLOOP     CONVERGENCE IS TESTED FOR EVERY NLOOP TIMES THE PROCESS
C             CHANGES THE SIMPLEX. AFTER INITIAL CONVERGENCE, NLOOP
C             FURTHER CHANGES ARE ALLOWED BEFORE TESTING FOR FINAL
C             CONVERGENCE
C   IQUAD     =1  IF FITTING OF QUADRATIC SURFACE REQUIRED
C             =0  IF NOT
C   SIMP      CRITERION FOR EXPANDING SIMPLEX TO OVERCOME ROUNDING
C             ERRORS BEFORE FITTING QUADRATIC SURFACE
C   VAR       ON EXIT, CONTAINS THE DIAGONAL ELEMENTS OF THE INVERSE OF
C             THE INFORMATION MATRIX
C   FUNCTN    A SUBROUTINE FUNCTN(F,FUNC) WHICH GIVEN AN ARRAY F OF
C             PARAMETER VALUES RETURNS THE FUNCTION VALUE FUNC
C   IFAULT    ON EXIT -
C                  IFAULT=1 IF NO. OF FUNCTION EVALUATIONS EXCEEDS NMAX,
C                  IFAULT=2 IF INFORMATION MATRIX IS NOT POSITIVE
C                           SEMI-DEFINITE
C                  IFAULT=3 IF NOP<1
C                  IFAULT=4 IF NLOOP<1
C                  IFAULT=0 OTHERWISE
C
C
C     F,STEP,NOP,IPRINT,STOPCR,NLOOP,IQUAD AND SIMP (SEE ** BELOW)
C     MUST BE SET BY THE CALLING PROGRAM
C
C
C **  IF IQUAD = 0 , SIMP IS NOT USED
C
C     F,STEP ( AND VAR IF IQUAD = 1 ) MUST HAVE DIMENSION AT LEAST NOP
C     IN THE CALLING PROGRAM
C
      IMPLICIT REAL(8) (A-H,O-Z)
      INTEGER NOP,NMAX,IPRINT,NLOOP,IQUAD,IW,IFAULT
      REAL(8) F(*),STEP(*)
      REAL(8) G(56,55),H(51),PBAR(55),PSTAR(55),PSTST(55)
      REAL(8)  AVAL(55),BMAT(10),PMIN(55),VC(10),VAR(*),TEMP(55)
      EQUIVALENCE(PMIN(1),PSTAR(1)),(AVAL(1),PBAR(1))

      REAL(8) A,B,C
      INTEGER NAP,LOOP,IFLAG,I,IROW,J,NP1,NEVAL,IMAX,IMIN,I1,I2,J1,K,L,
     +        NULLTY,IRANK,IJ,JK,IJK,II,IM1,JJ
C
C     A IS REFLECTION COEFFICIENT, B IS CONTRACTION COEFFICIENT, AND
C     C IS EXPANSION COEFFICIENT
      A=1.0D0
      B=0.5D0
      C=2.0D0
C
C     IF PROGRESS REPORTS DESIRED, PRINT HEADING FOR THEM
      IF(IPRINT)5,5,55
   55 WRITE(IW,100)IPRINT
  100 FORMAT(//' PROGRESS REPORT EVERY',I4,' FUNCTION EVALUATIONS'//
     * 'EVAL. NO.  FUNC. VALUE ',10X,'PARAMETERS')
C
    5 IFAULT = 0
      IF(NOP.LE.0) IFAULT=3
      IF(NLOOP.LE.0) IFAULT = 4
      IF (IFAULT.NE.0) RETURN
C
C     NAP IS THE NUMBER OF PARAMETERS TO BE VARIED,I.E. WITH STEP NOT 0.
      NAP=0
      LOOP=0
      IFLAG=0
      DO 1 I=1,NOP
      IF(STEP(I).NE.0.0D0) NAP=NAP+1
    1 CONTINUE
C
C     IF NAP=0 EVALUATE FUNCTION AT STARTING POINT AND RETURN
      IF(NAP)3,4,3
    4 CALL FUNCTN (F,FUNC)
      RETURN
C
C     SET UP INITIAL SIMPLEX
    3 DO 6 I=1,NOP
    6 G(1,I)=F(I)
      IROW=2
      DO 7 I=1,NOP
      IF(STEP(I).EQ.0.0D0) GO TO 7
      DO 9 J = 1,NOP
    9 G(IROW,J)=F(J)
      G(IROW,I)=G(IROW,I)+STEP(I)
      IROW=IROW+1
    7 CONTINUE
      NP1=NAP+1
      NEVAL=0
      DO 10 I=1,NP1
      DO 11 J=1,NOP
   11 F(J)=G(I,J)
      CALL FUNCTN(F,H(I))
      NEVAL=NEVAL+1
C
C     ALL POINTS OF THE INITIAL SIMPLEX ARE OUTPUT IF IPRINT> 0
      IF(IPRINT)10,10,12
   12 WRITE(IW,101)NEVAL,H(I),(F(J),J=1,NOP)
  101 FORMAT(3X,I4,3X,7G14.6/(24X,6G14.6))
   10 CONTINUE
C
C     NOW FOLLOWS THE BASIC LOOP, I.E. GIVEN A SIMPLEX, TO DETERMINE
C     NEW SIMPLEX AND TEST FOR CONVERGENCE AS REQUIRED (FOLLOWING THE
C     FLOW CHART GIVEN IN NELDER AND MEAD)
C
C     TO STATEMENT 13 , DETERMINE MAXIMUM AND MINIMUM POINTS OF
C     CURRENT SIMPLEX (FUNCTION VALUES ARE HMAX AND HMIN)
   45 LOOP=LOOP+1
      IMAX=1
      IMIN=1
      DO 13 I=2,NP1
      IF(H(I)-H(IMAX))15,15,14
   14 IMAX=I
   15 IF(H(I)-H(IMIN))16,13,13
   16 IMIN = I
   13 CONTINUE
      HMAX = H(IMAX)
      HMIN = H(IMIN)
C
C     TO STATEMENT 18 , FIND CENTROID OF ALL VERTICES, EXCLUDING THE
C     MAXIMUM
      DO 17 I=1,NOP
   17 PBAR(I)=0.0D0
      DO 18 I=1,NP1
      IF(I-IMAX)19,18,19
   19 DO 20 J=1,NOP
   20 PBAR(J) = PBAR(J)+G(I,J)
   18 CONTINUE
      DO 602 J=1,NOP
  602 PBAR(J) = PBAR(J)/NAP
C
C     REFLECT MAXIMUM THROUGH PBAR TO PSTAR, AND EVALUATE FUNCTION AT
C     PSTAR (GIVING HSTAR)
      DO 21 I=1,NOP
   21 PSTAR(I)=A*(PBAR(I)-G(IMAX,I))+PBAR(I)
      CALL FUNCTN (PSTAR,HSTAR)
C
C     NEXT 3 STATEMENTS TEST IF PROGRESS REPORT IS REQUIRED AT PRESENT,
C     AND IF SO, PROVIDE ONE
C     THIS PROCEDURE OCCURS FREQUENTLY IN THE PROGRAM
      NEVAL=NEVAL+1
      IF(IPRINT)57,57,56
   56 IF(MOD(NEVAL,IPRINT).EQ.0) WRITE(IW,101)NEVAL,HSTAR,(PSTAR(J),
     1    J=1,NOP)
   57 IF(HSTAR-HMIN)22,23,23
C
C     IF HSTAR LESS THAN HMIN REFLECT PBAR THROUGH PSTAR TO GIVE PSTST,
C     AND EVALUATE FUNCTION THERE (GIVING HSTST)
   22 DO 24 I=1,NOP
   24 PSTST(I)=C*(PSTAR(I)-PBAR(I))+PBAR(I)
      CALL FUNCTN (PSTST,HSTST)
      NEVAL=NEVAL+1
      IF(IPRINT)60,60,59
   59 IF(MOD(NEVAL,IPRINT).EQ.0) WRITE(IW,101)NEVAL,HSTST,(PSTST(J),
     1     J=1,NOP)
   60 IF(HSTST-HMIN)25,26,26
C
C     IF HSTST LESS THAN HMIN REPLACE MAXIMUM POINT OF CURRENT SIMPLEX
C     BY PSTST AND HMAX BY HSTAR, THEN TEST (STATEMENT 26 ONWARD)
   25 DO 27 I=1,NOP
      IF(STEP(I).NE.0.0D0) G(IMAX,I)=PSTST(I)
   27 CONTINUE
      H(IMAX)=HSTST
      GO TO 41
C
C     IF HSTAR NOT LESS THAN HMIN, TEST IF HSTAR GREATER THAN FUNCTION
C     VALUE AT ALL VERTICES OTHER THAN MAXIMUM ONE
   23 DO 28 I=1,NP1
      IF(I-IMAX)29,28,29
   29 IF(HSTAR-H(I))26,28,28
   28 CONTINUE
C
C     IF IT IS LESS THAN (AT LEAST) ONE OF THE VERTICES, REPLACE MAXIMUM
C     POINT OF CURRENT SIMPLEX BY PSTAR AND HMAX BY HSTAR, THEN TEST
C     (STATEMENT 26 ONWARD)
C
C     IF HSTAR GREATER THAN ALL FUNCTION VALUES EXCLUDING THE MAXIMUM,
C     TEST IF HSTAR GREATER THAN HMAX
C     IF NOT, REPLACE MAXIMUM POINT BY PSTAR AND HMAX BY HSTAR FOR
C     WHICHEVER SIMPLEX NOW IN STORE (I.E. DEPENDING ON WHETHER HSTAR
C     WAS GREATER OR LESS THAN HMAX), CALCULATE THE CONTRACTED POINT
C     PSTST AND THE FUNCTION VALUE THERE, HSTST
      IF(HSTAR-HMAX)30,30,31
   30 DO 32 I=1,NOP
      IF(STEP(I).NE.0.0D0) G(IMAX,I)=PSTAR(I)
   32 CONTINUE
      HMAX=HSTAR
      H(IMAX)=HSTAR
   31 DO 33 I=1,NOP
   33 PSTST(I)=B*G(IMAX,I)+(1.0D0-B)*PBAR(I)
      CALL FUNCTN(PSTST,HSTST)
      NEVAL=NEVAL+1
      IF(IPRINT) 63,63,62
   62 IF(MOD(NEVAL,IPRINT).EQ.0) WRITE(IW,101)NEVAL,HSTST,(PSTST(J),
     1    J=1,NOP)
   63 IF(HSTST-HMAX)35,35,34
C
C     IF HSTST LESS THAN HMAX, REPLACE MAXIMUM POINT BY PSTST AND HMAX
C     BY HSTST, THEN TEST (STATEMENT 41 ONWARD)
   35 DO 36 I=1,NOP
      IF(STEP(I).NE.0.0D0) G(IMAX,I)=PSTST(I)
   36 CONTINUE
      H(IMAX)=HSTST
      GO TO 41
C
C     IF HSTST NOT LESS THAN HMAX, REPLACE EACH POINT OF THE CURRENT
C     SIMPLEX BY A POINT MIDWAY BETWEEN ITS CURRENT POSITION AND THE
C     POSITION OF THE MINIMUM POINT OF THE CURRENT SIMPLEX. EVALUATE
C     FUNCTION AT EACH NEW VERTEX THEN TEST (STATEMENT 41 ONWARD)
   34 DO 38 I=1,NP1
      IF(I.EQ.IMIN) GO TO 38
      DO 39 J=1,NOP
      IF(STEP(J).NE.0.0D0) G(I,J)=(G(I,J)+G(IMIN,J))/2.0D0
   39 F(J) = G(I,J)
      CALL FUNCTN (F,H(I))
      NEVAL=NEVAL+1
      IF(IPRINT)38,38,65
   65 IF(MOD(NEVAL,IPRINT).EQ.0)WRITE(IW,101)NEVAL,H(I),(F(J),J=1,NOP)
   38 CONTINUE
      GO TO 41
   26 DO 40 I=1,NOP
      IF(STEP(I).NE.0.0D0) G(IMAX,I)=PSTAR(I)
   40 CONTINUE
      H(IMAX)=HSTAR
C
C     IF LOOP = NLOOP, BEGIN TESTS FOR CONVERGENCE
C     OTHERWISE, GO BACK TO BEGINNING OF BASIC LOOP
   41 IF(LOOP-NLOOP)45,46,45
C
C     TEST FOR CONVERGENCE
C     CALCULATE MEAN AND STANDARD DEVIATION OF FUNCTION VALUES OF
C     CURRENT SIMPLEX
   46 HSTD=0.0D0
      HMEAN=0.0D0
      DO 42 I=1,NP1
   42 HMEAN=HMEAN+H(I)
      HMEAN=HMEAN/NP1
      DO 601 I=1,NP1
  601 HSTD = HSTD+(H(I)-HMEAN)**2
      HSTD = SQRT(HSTD/NP1)
C
C     CALCULATE CENTROID OF CURRENT SIMPLEX, F, AND FUNCTION VALUE
C     THERE, FUNC
      DO 53 I=1,NOP
      IF(STEP(I).EQ.0.0D0) GO TO 53
      F(I)=0.0D0
      DO 54 J=1,NP1
   54 F(I)=F(I)+G(J,I)
      F(I)=F(I)/NP1
   53 CONTINUE
      CALL FUNCTN (F,FUNC)
      NEVAL=NEVAL+1
      IF(IPRINT.LE.0) GO TO 700
      IF(MOD(NEVAL,IPRINT).EQ.0)WRITE(IW,101)NEVAL,FUNC,(F(J),J=1,NOP)
C
C     IF NUMBER OF FUNCTION EVALUATIONS TO DATE HAS OVERRUN THE LIMIT
C     SET (NMAX), SET IFAULT = 1 AND RETURN
  700 IF(NEVAL-NMAX) 44,44,43
   43 IF(IPRINT)68,67,67
   67 WRITE(IW,102)NMAX
  102 FORMAT(40H NUMBER OF FUNCTION EVALUATIONS EXCEEDS ,I4)
      WRITE(IW,103)HSTD
  103 FORMAT(51H STANDARD ERROR OF FUNCTION VALUES OF LAST SIMPLEX ,
     1   G13.6)
      WRITE(IW,104)(F(I),I=1,NOP)
  104 FORMAT(27H CENTROID OF LAST SIMPLEX  ,6G14.6/(28X,6G14.6))
      WRITE(IW,105)FUNC
  105 FORMAT(31H  FUNCTION VALUE AT CENTROID   ,G13.6)
   68 IFAULT = 1
      RETURN
   44 IF(HSTD-STOPCR)72,48,48
C
C     IF THE STANDARD DEVIATION CALCULATED ABOVE IS NOT LESS THAN THE
C     CRITERION (STOPCR), SET IFLAG AND LOOP TO ZERO AND BEGIN BASIC
C     LOOP AGAIN
   48 IFLAG=0
      LOOP=0
      GO TO 45
   72 IF(IPRINT)47,70,70
   70 WRITE(IW,106)
  106 FORMAT(2H */33H  INITIAL EVIDENCE OF CONVERGENCE)
      WRITE(IW,104)(F(I),I=1,NOP)
      WRITE(IW,105)FUNC
C
C     IF THE STANDARD DEVIATION IS LESS THAN THE CRITERION TEST IFLAG
C     IFLAG=0 THERE WAS NO EVIDENCE OF CONVERGENCE ON LAST TEST
C     IFLAG=1 THERE WAS EVIDENCE OF CONVERGENCE ON LAST TEST
   47 IF(IFLAG)49,50,49
C
C     IF IFLAG=0, SET IFLAG=1 AND SAVE MEAN OF FUNCTION VALUES OF
C     CURRENT SIMPLEX. GO TO BEGINNING OF BASIC LOOP.
   50 IFLAG=1
   51 SAVEMN = HMEAN
      LOOP=0
      GO TO 45
C
C     IF IFLAG=1, TEST IF CHANGE IN MEAN IS LESS THAN THE CRITERION
C     (STOPCR). IF IT IS, PROCESS HAS CONVERGED
   49 IF(ABS(SAVEMN-HMEAN).GE.STOPCR)GO TO 51
      IF(IPRINT)74,73,73
   73 WRITE(IW,107)NEVAL
  107 FORMAT(5(/),36H PROCESS CONVERGES ON MINIMUM AFTER ,I4,21H FUNCTIO
     1N EVALUATIONS///)
      WRITE(IW,108)(F(I),I=1,NOP)
  108 FORMAT(14H MINIMUM AT   ,6G14.6/(14X,6G14.6))
      WRITE(IW,109)FUNC
  109 FORMAT(//26H MINIMUM FUNCTION VALUE   ,G13.6)
      WRITE(IW,110)
  110 FORMAT(///16H END  OF  SEARCH/1X,15(1H*))
   74 CONTINUE
C
C
C***********************************************************************
C
C
C     IF THE USER DOES NOT REQUIRE THE FITTING OF A QUADRATIC
C     SURFACE IN THE REGION OF THE MINIMUM ( FOR ESTIMATION
C     OF THE VARIANCE - COVARIANCE MATRIX), THE REMAINDER OF
C     THE PROGRAM DECK (INCLUDING AND AFTER THE CARD OF ***S)
C     MAY BE DETACHED. AN END CARD MUST BE INCLUDED IN THE
C     SHORTENED DECK .
C
      IF(IQUAD)200,75,200
   75 RETURN
  200 IF(IPRINT)233,232,232
  232 WRITE(IW,301)
  301 FORMAT(1H1,10(/),' FITTING OF QUADRATIC SURFACE IN REGION OF
     *MINIMUM'///)
C
C     THE FITTING OF THE QUADRATIC SURFACE FOLLOWS EXACTLY THE
C     PROCEDURE OUTLINED BY NELDER AND MEAD , AND WHERE
C     POSSIBLE THE NOTATION IN THE COMMENT CARDS WILL CORRESPOND
C     TO THAT OF NELDER AND MEAD
  233 NEVAL=0
C
C     FURTHER FUNCTION EVALUATIONS REQUIRED ARE COUNTED IN NEVAL
C
C     THE FINAL SIMPLEX IS EXPANDED TO OVERCOME ROUNDING ERRORS
      DO 201 I=1,NP1
  204 TEST=ABS(H(I)-FUNC)
      IF(TEST-SIMP)202,201,201
  202 DO 203 J=1,NOP
      IF(STEP(J).NE.0.0D0) G(I,J)=(G(I,J)-F(J))+G(I,J)
  203 PSTST(J)=G(I,J)
      CALL FUNCTN(PSTST,H(I))
      NEVAL=NEVAL+1
      GO TO 204
  201 CONTINUE
      A0 = H(1)
C
C     THE FUNCTION VALUES Y0(I) ARE CALCULATED AND STORED IN AVAL
      DO 205 I=1,NAP
      I1=I+1
      DO 206 J=1,NOP
  206 PSTAR(J)=(G(1,J)+G(I1,J))/2.0D0
      CALL FUNCTN(PSTAR,AVAL(I))
      NEVAL=NEVAL+1
  205 CONTINUE
C
C     THE MATRIX B(I,J) IS CALCULATED , AND THE LOWER DIAGONAL SECTION
C     STORED IN THE VECTOR BMAT
C     BMAT CONTAINS THE INFORMATION MATRIX
      DO 207 I=1,NAP
      I1=I-1
      I2=I+1
      IF(I1.LT.1)GOTO207
      DO 208 J=1,I1
      J1=J+1
      DO 209 K=1,NOP
  209 PSTST(K)=(G(I2,K)+G(J1,K))/2.0D0
      CALL FUNCTN(PSTST,HSTST)
      NEVAL=NEVAL+1
      L=I*(I-1)/2+J
      BMAT(L)=2.0D0*(HSTST+A0-AVAL(I)-AVAL(J))
  208 CONTINUE
  207 CONTINUE
      L=0
      DO 210 I=1,NAP
      I1=I+1
      L=L+I
      BMAT(L)=2.0D0*(H(I1)+A0-2.0D0*AVAL(I))
  210 CONTINUE
C
C     THE VECTOR A(I) IS CALCULATED , AND STORED IN AVAL
      DO 237 I=1,NAP
      I1=I+1
  237 AVAL(I)=2.0D0*AVAL(I)-(H(I1)+3.0D0*A0)/2.0D0
C
C     THE MATRIX Q IS CALCULATED , AND STORED IN THE MATRIX G
C     (IN FACT , CONSIDERING THE USUAL ORIENTATION OF ROWS AND COLUMNS,
C     TRANS(Q) IS STORED IN G)
      DO 219 I=1,NOP
  219 PMIN(I)=G(1,I)
      DO 211 I=1,NAP
      I1=I+1
      DO 211 J=1,NOP
      G(I1,J)=G(I1,J)-G(1,J)
  211 CONTINUE
      DO 212 I=1,NAP
      I1=I+1
      DO 212 J=1,NOP
      G(I,J)=G(I1,J)
  212 CONTINUE
C
      CALL SYMINV (BMAT,NAP,BMAT,TEMP,NULLTY,IFAULT,RMAX)
      IF(IFAULT.NE.0) GO TO 450
      IRANK = NAP-NULLTY
      GO TO 441
  450 IF(IPRINT.GE.0) WRITE(IW,302)
  302 FORMAT(49H MATRIX TO BE INVERTED NOT POSITIVE SEMI-DEFINITE)
      IFAULT = 2
      RETURN
C
C     (B**-1)*A*0.5 IS CALCULATED,AND STORED IN H
  441 DO 213 I=1,NAP
      H(I)=0.0D0
      DO 214 J=1,NAP
      IF(J-I)216,216,215
  216 IJ=I*(I-1)/2+J
      GO TO 217
  215 IJ=J*(J-1)/2+I
  217 H(I)=H(I)+BMAT(IJ)*AVAL(J)
  214 CONTINUE
  213 CONTINUE
C
C     YMIN (THE ESTIMATED MINIMUM VALUE) AND PMIN (ITS POSITION)
C     ARE CALCULATED
      YMIN=0
      DO 218 I=1,NAP
  218 YMIN=YMIN+H(I)*AVAL(I)
      YMIN=A0-YMIN
      DO 220 I=1,NOP
      PSTST(I)=0
      DO 220 J=1,NAP
  220 PSTST(I)=PSTST(I)+H(J)*G(J,I)
      DO 221 I=1,NOP
  221 PMIN(I)=PMIN(I)-PSTST(I)
      IF(IPRINT)223,222,222
  222 WRITE(IW,303)YMIN,(PMIN(I),I=1,NOP)
  303 FORMAT(42H MINIMUM OF FITTED QUADRATIC SURFACE IS   ,G13.6,4H  AT/
     *(6G14.6))
      WRITE(IW,304)FUNC,(F(I),I=1,NOP)
  304 FORMAT(/42H COMPARE WITH MINIMUM FOUND BY ITERATION  ,G13.6,4H  AT
     */(6G14.6))
      WRITE(IW,305)
  305 FORMAT(/58H IF DIFFERENCE IS LARGE , INFORMATION MATRIX IS INACCUR
     *ATE///)
C
C     Q*(B**-1)*TRANS(Q)*0.5 IS CALCULATED,AND ITS LOWER DIAGONAL
C     SECTION STORED IN THE VECTOR VC
  223 DO 224 I=1,NOP
      DO 225 J=1,NAP
      H(J)=0.0D0
      DO 226 K=1,NAP
      IF(K-J)227,227,228
  227 JK=J*(J-1)/2+K
      GO TO 229
  228 JK=K*(K-1)/2+J
  229 H(J) = H(J)+BMAT(JK)*G(K,I)*0.5D0
  226 CONTINUE
  225 CONTINUE
      DO 230 J=I,NOP
      IJ = J*(J-1)/2+I
      VC(IJ)=0.0D0
      DO 231 K=1,NAP
  231 VC(IJ)=VC(IJ)+H(K)*G(K,J)
  230 CONTINUE
  224 CONTINUE
C
C     THE DIAGONAL ELEMENTS OF VC ARE STORED IN VAR FOR RETURN TO THE
C     CALLING PROGRAMME
      J = 0
      DO 234 I=1,NOP
      J = J+I
  234 VAR(I)=VC(J)
      IF(IPRINT)235,236,236
  236 WRITE(IW,306)IRANK
  306 FORMAT(28H INFORMATION MATRIX HAS RANK,I3//
     1    42H GENERALIZED INVERSE OF INFORMATION MATRIX/)
      IJK=1
      GO TO 710
  717 CONTINUE
      WRITE(IW,308)
  308 FORMAT(///113H IF THE FUNCTION MINIMISED WAS -LOG(LIKELIHOOD) , TH
     *IS MATRIX IS THE VARIANCE-COVARIANCE MATRIX OF THE PARAMETERS/)
      WRITE(IW,309)
  309 FORMAT(120H IF THE FUNCTION MINIMISED WAS THE SUM OF SQUARES OF RE
     *SIDUALS , THIS MATRIX MUST BE MULTIPLIED BY TWICE THE ESTIMATE OF
     */79H THE RESIDUAL VARIANCE TO GIVE THE VARIANCE-COVARIANCE MATRIX
     *OF THE PARAMETERS///)
      CALL SYMINV (VC,NAP,BMAT,TEMP,NULLTY,IFAULT,RMAX)
C
C     BMAT NOW CONTAINS THE INFORMATION MATRIX
C
C     PRINTING OUT THE INFORMATION MATRIX
      WRITE(IW,320)
  320 FORMAT(19H INFORMATION MATRIX/)
      IJK=3
      GO TO 710
  711 IJK=2
      II = 0
      IJ = 0
      DO 701 I = 1,NOP
      II = II+I
      IF(VC(II).LE.0.0D0) GO TO 702
      VC(II) = 1.0D0/SQRT(VC(II))
      GO TO 703
  702 VC(II) = 0.0D0
  703 IF(I.EQ.1) GO TO 701
      IM1 = I-1
      JJ = 0
      DO 705 J=1,IM1
      JJ = JJ+J
      IJ = IJ+1
  705 VC(IJ) = VC(IJ)*VC(II)*VC(JJ)
  701 IJ = IJ+1
      WRITE(IW,312)
  312 FORMAT(///19H CORRELATION MATRIX/)
      II = 0
      DO 706 I=1,NOP
      II = II+I
      IF(VC(II).NE.0.0D0) VC(II) =1.0D0
  706 CONTINUE
      GO TO 710
  712 WRITE(IW,310)NEVAL
  310 FORMAT(11H A FURTHER ,I4,36H FUNCTION EVALUATIONS HAVE BEEN USED//
     */)
      WRITE(IW,311)
  311 FORMAT(37H END  OF  QUADRATIC  SURFACE  FITTING/1X,36(1H*))
  235 RETURN
  710 L=1
  716 IF(L.GT.NOP) GO TO (717,712,711),IJK
      II=L*(L-1)/2
      DO 713 I=L,NOP
      I1 = II+L
      II=II+I
      I2=MIN0(II,I1+5)
      IF(IJK.EQ.3) GO TO 718
      WRITE(IW,714)(VC(J),J=I1,I2)
      GO TO 713
  718 WRITE(IW,714)(BMAT(J),J=I1,I2)
  713 CONTINUE
  714 FORMAT(15X,6G14.6)
      WRITE(IW,715)
  715 FORMAT(1H0)
      L=L+6
      GO TO 716
      END
C
      SUBROUTINE SYMINV (A,N,C,W,NULLTY,IFAULT,RMAX)
C
C     ALGORITHM AS 7 J.R.STATIST.SOC. C, (1968) VOL.17, NO.2.
C
C     FORMS IN C( ) AS LOWER TRIANGLE, A GENERALISED INVERSE
C     OF THE POSITIVE SEMI-DEFINITE SYMMETRIC MATRIX A( )
C     ORDER N, STORED AS LOWER TRIANGLE.
C     C( ) MAY COINCIDE WITH A( ) . NULLTY IS RETURNED AS THE NULLITY
C     OF A( ). IFAULT IS RETURNED AS 1 IF N.LT.1, OTHERWISE ZERO
C     W( ) IS A WORK ARRAY OF LENGTH AT LEAST N THAT IS ALLOCATED BY
C     THE CALLING ROUTINE
C
C
      IMPLICIT REAL(8) (A-H,O-Z)
      DOUBLE PRECISION A(*),C(*),W(*)
      DOUBLE PRECISION RMAX
      INTEGER N,IFAULT,NULLTY,NROW,NN,IROW,NDIAG,L,ICOL,JCOL,I,MDIAG,K,J
C
      NROW=N
      IFAULT=1
      IF(NROW.LE.0) GO TO 100
      IFAULT=0
      CALL CHOLA(A,NROW,C,NULLTY,IFAULT,RMAX,W)
      IF(IFAULT.NE.0) GO TO 100
      NN=(NROW*(NROW+1))/2
      IROW=NROW
      NDIAG=NN
   16 IF(C(NDIAG).EQ.0.0D0) GO TO 11
      L=NDIAG
      DO 10 I=IROW,NROW
      W(I)=C(L)
      L=L+I
   10 CONTINUE
      ICOL=NROW
      JCOL=NN
      MDIAG=NN
   15 L=JCOL
      X=0.0D0
      IF(ICOL.EQ.IROW) X=1.0D0/W(IROW)
      K=NROW
   13 IF(K.EQ.IROW) GO TO 12
      X=X-W(K)*C(L)
      K=K-1
      L=L-1
      IF(L.GT.MDIAG) L=L-K+1
      GO TO 13
   12 C(L)=X/W(IROW)
      IF(ICOL.EQ.IROW) GO TO 14
      MDIAG=MDIAG-ICOL
      ICOL=ICOL-1
      JCOL=JCOL-1
      GO TO 15
   11 L=NDIAG
      DO 17 J=IROW,NROW
      C(L)=0.0D0
      L=L+J
   17 CONTINUE
   14 NDIAG=NDIAG-IROW
      IROW=IROW-1
      IF(IROW.NE.0) GO TO 16
  100 RETURN
      END
C
      SUBROUTINE CHOLA(A,N,U,NULLTY,IFAULT,RMAX,R)
C
C     ALGORITHM AS 6 J.R.STATIST.SOC. C, (1968) VOL.17 , NO.2
C
C     GIVEN A SYMMETRIC MATRIX ORDER N AS LOWER TRIANGLE IN A( ),
C     CALCULATES AN UPPER TRIANGLE, U( ) , SUCH THAT UPRIME*U = A.
C     U( ) MAY COINCIDE WITH A( ). A( ) MUST BE POSITIVE SEMI-
C     DEFINITE.
C     ETA IS SET TO MULTIPLYING FACTOR DETERMINING EFFECTIVE
C     ZERO FOR PIVOT.
C     NULLTY IS RETURNED AS NO. OF EFFECTIVE ZERO PIVOTS.
C     IFAULT IS RETURNED AS 1 IF N.LE.0, 2 IF A( ) IS NOT POSITIVE
C     SEMI-DEFINITE WITHIN THE TOLERANCE DEFINED BY ETA, OTHERWISE
C     ZERO.
C
C
C     CORRECTION ASR12 J.R.STATIST.SOC.C.(1974) ADDED
C
C     MODIFIED TEST FOR A SINGULARITY AND ACCURACY ESTIMATE BY
C     A.J.MILLER ADDED
C
C
      IMPLICIT REAL(8) (A-H,O-Z)
      DOUBLE PRECISION A(*),U(*),R(*)
      DOUBLE PRECISION ETA,RMAX,W,RSQ
      INTEGER N,IFAULT,NULLTY,J,K,L,ICOL,IROW,I,M
C
      ETA=1.0D-24
C
C     THE VALUE OF ETA WILL DEPEND ON THE WORD LENGTH OF THE
C     COMPUTER BEING USED
C
C
      IFAULT=1
      IF(N.LE.0) GO TO 100
      IFAULT=2
      NULLTY=0
      RMAX = ETA
      R(1) = ETA
      J=1
      K=0
      DO 10 ICOL=1,N
      L=0
      DO 11 IROW=1,ICOL
      K=K+1
      W=A(K)
      IF(IROW.EQ.ICOL) RSQ = (W*ETA)**2
      M=J
      DO 12 I=1,IROW
      L=L+1
      IF(I.EQ.IROW) GO TO 13
      W=W-U(L)*U(M)
      IF(IROW.EQ.ICOL) RSQ = RSQ + (U(L)**2*R(I))**2
      M=M+1
   12 CONTINUE
   13 IF(IROW.EQ.ICOL) GO TO 14
      IF(U(L).EQ.0.0D0) GO TO 21
      U(K)=W/U(L)
      GO TO 11
   21 U(K)=0.0D0
      IF(ABS(W).GT.ABS(RMAX*A(K))) GO TO 100
   11 CONTINUE
   14 RSQ = SQRT(RSQ)
      IF(ABS(W).LE.5.D0*RSQ) GO TO 20
      IF(W.LT.0.0D0) GO TO 100
      U(K)=SQRT(W)
      R(I) = RSQ/W
      IF(R(I).GT.RMAX) RMAX = R(I)
      GO TO 15
   20 U(K)=0.0D0
      NULLTY=NULLTY+1
   15 J=J+ICOL
   10 CONTINUE
      IFAULT=0
  100 RETURN
      END
C
C--------------------------------------------------------------------------
C
      FUNCTION GENBET(AA,BB,ISEED)
C     Always call with a negative seed to reset the RAN1 random generator

C                              LEGALITIES
C
C Code that appeared in an ACM publication is subject to their
C algorithms policy:
C
C     Submittal of an algorithm for publication in one of the ACM
C     Transactions implies that unrestricted use of the algorithm within a
C     computer is permissible . General permission to copy and distribute
C     the algorithm without fee is granted provided that the copies are not
C     made or distributed for direct commercial advantage.  The ACM
C     copyright notice and the title of the publication and its date appear,
C     and notice is given that copying is by permission of the Association
C     for Computing Machinery. To copy otherwise, or to republish, requires
C     a fee and/or specific permission.
C
C     Krogh, F. Algorithms Policy. ACM Tran. Math. Softw. 13(1987), 183-186.
C
C We place the Ranlib code that we have written in the public domain.  
C
C
C                                 NO WARRANTY
C     
C     WE PROVIDE ABSOLUTELY  NO WARRANTY  OF ANY  KIND  EITHER  EXPRESSED OR
C     IMPLIED,  INCLUDING BUT   NOT LIMITED TO,  THE  IMPLIED  WARRANTIES OF
C     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK
C     AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS  WITH YOU.  SHOULD
C     THIS PROGRAM PROVE  DEFECTIVE, YOU ASSUME  THE COST  OF  ALL NECESSARY
C     SERVICING, REPAIR OR CORRECTION.
C     
C     IN NO  EVENT  SHALL THE UNIVERSITY  OF TEXAS OR  ANY  OF ITS COMPONENT
C     INSTITUTIONS INCLUDING M. D.   ANDERSON HOSPITAL BE LIABLE  TO YOU FOR
C     DAMAGES, INCLUDING ANY  LOST PROFITS, LOST MONIES,   OR OTHER SPECIAL,
C     INCIDENTAL   OR  CONSEQUENTIAL DAMAGES   ARISING   OUT  OF  THE USE OR
C     INABILITY TO USE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA OR
C     ITS ANALYSIS BEING  RENDERED INACCURATE OR  LOSSES SUSTAINED  BY THIRD
C     PARTIES) THE PROGRAM.
C     
C**********************************************************************
C
C     REAL FUNCTION GENBET( A, B )
C               GeNerate BETa random deviate
C
C                              Function
C
C     Returns a single random deviate from the beta distribution with
C     parameters A and B.  The density of the beta is
C               x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
C
C                              Method
C     R. C. H. Cheng
C     Generating Beta Variatew with Nonintegral Shape Parameters
C     Communications of the ACM, 21:317-322  (1978)
C     (Algorithms BB and BC)
C
C**********************************************************************
C     .. Parameters ..
C     EXPMAX  Close to the largest number that can be exponentiated
C     INFNTY  Close to the largest representable single precision number
      DOUBLE PRECISION GENBET, EXPMAX, INFNTY
      PARAMETER (EXPMAX=89.0D0, INFNTY=1.0E38)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION AA,BB
      INTEGER ISEED
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,ALPHA,B,BETA,DELTA,GAMMA,K1,K2,OLDA,OLDB,R,S,T,
     +       U1,U2,V,W,Y,Z
      LOGICAL QSAME
C     ..
C     .. External & Intrinsic Functions ..
      DOUBLE PRECISION RAN1, MA(55)
      EXTERNAL RAN1
      INTEGER INEXT,INEXTP
      INTRINSIC EXP,LOG,MAX,MIN,SQRT
C     ..
C     .. Save statement ..
      SAVE OLDA,OLDB,ALPHA,BETA,GAMMA,K1,K2
C     ..
C     .. Data statements ..
      DATA OLDA,OLDB/-1,-1/
C     ..
C     .. Executable Statements ..
      QSAME = (OLDA.EQ.AA) .AND. (OLDB.EQ.BB)
      IF (QSAME) GO TO 20
      IF (.NOT. (AA.LE.0.0.OR.BB.LE.0.0)) GO TO 10
      WRITE (*,*) ' AA OR BB <= 0 IN GENBET - ABORT!'
      WRITE (*,*) ' AA: ',AA,' BB ',BB
      STOP ' AA OR BB <= 0 IN GENBET - ABORT!'

   10 OLDA = AA
      OLDB = BB
   20 IF (.NOT. (MIN(AA,BB).GT.1.0)) GO TO 100


C     Alborithm BB

C     Initialize
C
      IF (QSAME) GO TO 30
      A = MIN(AA,BB)
      B = MAX(AA,BB)
      ALPHA = A + B
      BETA = SQRT((ALPHA-2.0)/ (2.0*A*B-ALPHA))
      GAMMA = A + 1.0/BETA
   30 CONTINUE
   40 U1 = RAN1(ISEED,MA,INEXT,INEXTP)
C
C     Step 1
C
      U2 = RAN1(ISEED,MA,INEXT,INEXTP)
      V = BETA*LOG(U1/ (1.0-U1))
      IF (.NOT. (V.GT.EXPMAX)) GO TO 50
      W = INFNTY
      GO TO 60

   50 W = A*EXP(V)
   60 Z = U1**2*U2
      R = GAMMA*V - 1.3862944
      S = A + R - W
C
C     Step 2
C
      IF ((S+2.609438).GE. (5.0*Z)) GO TO 70
C
C     Step 3
C
      T = LOG(Z)
      IF (S.GT.T) GO TO 70
C
C     Step 4
C
      IF ((R+ALPHA*LOG(ALPHA/ (B+W))).LT.T) GO TO 40
C
C     Step 5
C
   70 IF (.NOT. (AA.EQ.A)) GO TO 80
      GENBET = W/ (B+W)
      GO TO 90

   80 GENBET = B/ (B+W)
   90 GO TO 230


C     Algorithm BC
C
C     Initialize
C
  100 IF (QSAME) GO TO 110
      A = MAX(AA,BB)
      B = MIN(AA,BB)
      ALPHA = A + B
      BETA = 1.0/B
      DELTA = 1.0 + A - B
      K1 = DELTA* (0.0138889+0.0416667*B)/ (A*BETA-0.777778)
      K2 = 0.25 + (0.5+0.25/DELTA)*B
  110 CONTINUE
  120 U1 = RAN1(ISEED,MA,INEXT,INEXTP)
C
C     Step 1
C
      U2 = RAN1(ISEED,MA,INEXT,INEXTP)
      IF (U1.GE.0.5) GO TO 130
C
C     Step 2
C
      Y = U1*U2
      Z = U1*Y
      IF ((0.25*U2+Z-Y).GE.K1) GO TO 120
      GO TO 170
C
C     Step 3
C
  130 Z = U1**2*U2
      IF (.NOT. (Z.LE.0.25)) GO TO 160
      V = BETA*LOG(U1/ (1.0-U1))
      IF (.NOT. (V.GT.EXPMAX)) GO TO 140
      W = INFNTY
      GO TO 150

  140 W = A*EXP(V)
  150 GO TO 200

  160 IF (Z.GE.K2) GO TO 120
C
C     Step 4
C
C
C     Step 5
C
  170 V = BETA*LOG(U1/ (1.0-U1))
      IF (.NOT. (V.GT.EXPMAX)) GO TO 180
      W = INFNTY
      GO TO 190

  180 W = A*EXP(V)
  190 IF ((ALPHA* (LOG(ALPHA/ (B+W))+V)-1.3862944).LT.LOG(Z)) GO TO 120
C
C     Step 6
C
  200 IF (.NOT. (A.EQ.AA)) GO TO 210
      GENBET = W/ (B+W)
      GO TO 220

  210 GENBET = B/ (B+W)
  220 CONTINUE
  230 RETURN

      END
      
C -----------------------------------------------------------------------------

      SUBROUTINE INVM(A,AI,npp,NN)

C     This subroutine computes the inverse of a matrix

C     Global variables
      INTEGER npp,NN,T,TT
      REAL(8) A(npp,npp),AI(npp,npp)

C     Local variables
      REAL(8) INDEX(1000),D

C     Set up an identity matrix
      DO 12 T = 1,NN
        DO 11 TT = 1,NN
          AI(T,TT) = 0.0
   11   CONTINUE
        AI(T,T) = 1.0
   12 CONTINUE

C     Perform a LU-decomposition
      CALL LUDCMP(A,NN,npp,INDEX,D)

C     Find the inverse
      DO 13 TT = 1,NN
        CALL LUBKSB(A,NN,npp,INDEX,AI(1,TT))
   13 CONTINUE

      RETURN
      END

C -----------------------------------------------------------------------------

      SUBROUTINE LUBKSB(A,NN,npp,INDX,B)

      integer npp,NN,ttt,T,ll,TT
      REAL(8) A(NPP,NPP),INDX(NN),B(NN),SUM

      ttt=0
      DO 12 T=1,NN
       LL=INDX(T)
       SUM=B(LL)
       B(LL)=B(T)
       IF (ttt.NE.0)THEN
         DO 11 TT=ttt,T-1
           SUM=SUM-A(T,TT)*B(TT)
   11    CONTINUE
       ELSE IF (SUM.NE.0.) THEN
         ttt=T
       ENDIF
       B(T)=SUM
   12 CONTINUE
      DO 14 T=NN,1,-1
       SUM=B(T)
       IF(T.LT.NN)THEN
         DO 13 TT=T+1,NN
           SUM=SUM-A(T,TT)*B(TT)
   13    CONTINUE
       ENDIF
       B(T)=SUM/A(T,T)
   14 CONTINUE

      RETURN
      END

C -----------------------------------------------------------------------------

      SUBROUTINE LUDCMP(A,NN,NPP,INDX,D)

      INTEGER NMAX,NN,NPP
      REAL(8) TINY
      parameter(nmax=200,tiny=1.0E-20)
      REAL(8) A(NPP,NPP),INDX(NN),VV(NMAX),SUM,AAMAX,DUM,D
      INTEGER t,imax,tt,kK

      D=1.
      DO 12 t=1,NN
        AAMAX=0.
        DO 11 tt=1,NN
         IF (ABS(A(t,tt)).GT.AAMAX) AAMAX=ABS(A(t,tt))
   11   CONTINUE
       IF (AAMAX.EQ.0.) WRITE(*,*) 'Singular matrix.'
       VV(t)=1./AAMAX
   12 CONTINUE

      DO 19 tt=1,NN
        IF (tt.GT.1) THEN
         DO 14 t=1,tt-1
           SUM=A(t,tt)
           IF (t.GT.1)THEN
             DO 13 KK=1,T-1
               SUM=SUM-A(t,KK)*A(KK,tt)
   13        CONTINUE
             A(t,tt)=SUM
           ENDIF
   14    CONTINUE
        ENDIF
       AAMAX=0.
       DO 16 t=tt,NN
         SUM=A(t,tt)
         IF (tt.GT.1)THEN
           DO 15 KK=1,tt-1
             SUM=SUM-A(t,KK)*A(KK,tt)
   15      CONTINUE
           A(t,tt)=SUM
         ENDIF
         DUM=VV(t)*ABS(SUM)
         IF (DUM.GE.AAMAX) THEN
           IMAX=t
           AAMAX=DUM
         ENDIF
   16  CONTINUE
       IF (tt.NE.IMAX)THEN
         DO 17 KK=1,NN
           DUM=A(IMAX,KK)
           A(IMAX,KK)=A(tt,KK)
           A(tt,KK)=DUM
   17    CONTINUE
         D=-D
         VV(IMAX)=VV(TT)
       ENDIF
       INDX(tt)=IMAX
       IF(tt.NE.NN)THEN
         IF(A(tt,tt).EQ.0.)A(tt,tt)=TINY
         DUM=1./A(tt,tt)
         DO 18 t=tt+1,NN
           A(t,tt)=A(t,tt)*DUM
   18    CONTINUE
       ENDIF
   19 CONTINUE
      IF(A(NN,NN).EQ.0.)A(NN,NN)=TINY

      RETURN
      END
      




      SUBROUTINE GENMUL(prop,TC,C,WithRep,ISEED,MA,INEXT,INEXTP)

C     This subroutine generates a multinomial random variable
C     Prop: expected proportion
C     TC: sample size
C     C: output
C     ISEED3: Seed (MUST BE NEGATIVE)
C     WithRep: FALSE for without replacement

      IMPLICIT NONE
C
C     Global variables
      REAL*8 prop(2000),C(2000)
      INTEGER ISEED,TC
      LOGICAL WithRep
C
C     Local variables
      REAL*8 tot,NTOT,MA(55)
      INTEGER A,A1,NP,DUM(100010),I,ICOUNT,INEXT,INEXTP
      REAL*8 RAN1
      EXTERNAL RAN1

C     Check if some catches were made!
      DO 9 a = 1,2000
       C(a) =  0.0
9     CONTINUE

      IF (TC.eq.0) RETURN

C     First total the probabilities and then normalise
      tot = 0.0
      DO 10 a = 1,2000
       tot = tot + prop(a)
10    CONTINUE
      IF (Tot.EQ.0) RETURN
      DO 11 a = 1,2000
       prop(a) = prop(a)/tot
11    CONTINUE
C
C     Adjsut to handle "without-replacement" sampling
      NTOT = 100000.0
      IF (WithRep.AND.tot.LT.14000) NTOT = tot
C
C     Do special set up for index matrix
      NP = 0
      DO 13 a = 1,2000
       ICOUNT = NINT(prop(a)*NTOT)
       IF (ICOUNT.GT.0) THEN
        DO 14 I = 1,ICOUNT
         DUM(NP+I) = a
14      CONTINUE
        NP = NP + ICOUNT
       ENDIF
13    CONTINUE

C     Do actual generation
      DO 15 I = 1,TC
C
C      Generate an age
C# FIX THIS
16     a1 = INT(RAN1(ISEED,ma,inext,inextp)*NP+1)
       a = DUM(a1)
       IF (a.LE.0.AND.WithRep) GOTO 16
       C(a) = C(a) + 1
       IF (WithRep) DUM(a1) = -1
C
15    CONTINUE

      RETURN
      END 


