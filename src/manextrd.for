C     Last change:  CA   12 Jun 2001    1:51 pm
C ------------------------------------------------------------------------------
C ------------------------------------------------------------------------------

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

      DOUBLE PRECISION FUNCTION XNORM (SIGG,MEAN,ISEED,MA,INEXT,INEXTP)

C     Uniform to Normal conversion routine

      DOUBLE PRECISION Z1, Z2, MEAN, SIGG, RAN1, MA(55)
      INTEGER ISEED, INEXT,INEXTP
      EXTERNAL RAN1

      Z1 = RAN1(ISEED,MA,INEXT,INEXTP)
    1 Z2 = RAN1(ISEED,MA,INEXT,INEXTP)
      IF (Z2.LE.0.D0) GO TO 1

      XNORM = SIN(6.238319D0*Z1) * SQRT(-2.D0*LOG(Z2)) * SIGG+MEAN

      RETURN
      END

C ------------------------------------------------------------------------------
C ------------------------------------------------------------------------------

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

C --------------------------------------------------------------------------
C --------------------------------------------------------------------------
C
      FUNCTION GAMMA (PARX,PARS,ISEED)
C
C     Approximate random gamma generator: PARX = mean; PARS= sd
C
      EXTERNAL RAN1
      DOUBLE PRECISION PARX,PARS,X,S,RAN1,MA(55),P,F,FN,GAMMA
      INTEGER ISEED,INEXT,INEXTP,I
C
      X = PARX
      S = PARS
C
C     Convert to gamma parameters
      s = s*s / x
C
      P = 0.D0
      IF (S.GT.0.D0) THEN
        X = X/S
        I = INT(X)
        FN = AINT(X)
        IF (I.GT.0) THEN
C
C         Count down from i to zero
100       IF (I.NE.0) THEN
            F = RAN1(ISEED,MA,INEXT,INEXTP)
            IF (F.LE.0.D0) THEN
              WRITE(*,*) 'WARNING F=0 in GAMMA',F
              F = 0.00001D0
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
          IF (F.LE.0.D0) THEN
            WRITE(*,*) 'WARNING F=0 in GAMMA',F
            F = 0.00001D0
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

C  ***********************************************************************

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
      PARAMETER (FACTOR=1.6D0,NTRY=50)
C
      IF(X1.EQ.X2)PAUSE 'You have to guess an initial range'
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
      IF(FB*FA.GT.0.D0) PAUSE 'Root must be bracketed for ZBRENT.'
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
      PAUSE 'ZBRENT exceeding maximum iterations.'
      PARS(M)=B
      RETURN
      END

