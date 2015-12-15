C     Last change:  CA   21 Feb 2011    3:53 pm

C    Version for use with combination trials to overcome 'PROB
C     INTEGRATESTO ZERO' problem when running trials with very low CVs.

C    Subroutines RDPARS & PHOUT moved Feb 2011
C    Changed to double precision 5/10/97, 15/1/98 + 9/2/98
C    Additional Prior parameters deleted 14/3/92
C    'D' threshold parameter added 14/3/92, deleted 5/5/93
C    Subroutine GENDAT removed: data now input directly  1/6/92
C    Individual CVs used in place of the average value over all years June92
C    Phaseout rule changed to match specification 1/11/92
C    Statements forbidding assessments in successive years deleted 19/11/92
C    Option to run assessments at fixed intervals or when data arrives.  15/1/93
C    Zero abundance estimate allowed for.   23/2/93
C    Rearranged for real applications (and dlm front end deleted): 5/5/93
C    Bias steps changed to mid points
C
C
C ***********************************************************************
C ***********************************************************************
C
C    CLC version 7
C
C ***********************************************************************
C
C     31 January 1994
C
C
C ***********************************************************************
C ***********************************************************************

      FUNCTION CLIMIT (CATCH,SIGHT,FMATRX,ISYR,IZYR,ZMULT,POP,G)

C     Run the CLA to obtain the nominal catch limit
C
      COMMON /MANPAR/ PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
      DOUBLE PRECISION PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
C
      COMMON /MANDAT/ ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER ISTART,IYEAR,NS,NZ
      DOUBLE PRECISION RKLO,RKHI
C
      INTEGER MAXSIM,MAXSTP,ISYR(0:*),IZYR(0:*),J,N,N2,I
      PARAMETER (MAXSIM = 10000000, MAXSTP=800)
      DOUBLE PRECISION CLIMIT,PRES(0:MAXSIM),QRES(0:MAXSIM),SS0,SS1,SS2,
     +       SS3,CATCH(0:*),SIGHT(0:*),FMATRX(0:*),ZMULT(0:*),POP(0:*),
     +       G(0:*),STKSIM,CONTRL,BINC,YINC,PTOT,R,DK,D,DD,P,RK,Q
C
C     Local variables:
      DOUBLE PRECISION SF,Y(0:MAXSTP),B(0:MAXSTP),RLGB(0:MAXSTP)
      DOUBLE PRECISION RSCAL, MSCAL,TERM
      INTEGER NB,NR, II
C     NB     Number of bias steps
C     NR     Number of steps for Y (the productivity parameter)
C     SF     Deviance scale factor (SF = .5 / PSCALE**2)
C
      CLIMIT = 0.D0
      IF (NS .LE. 0) RETURN
C
C     Set deviance scale factor  S = 1/PSCALE**2
      SF = 0.5D0 / (PSCALE * PSCALE)
C
C     Check the sizes of the Y and B arrays are large enough
      IF (PNBSTP.GT.MAXSTP .OR. PNYSTP.GT.MAXSTP) STOP 
     +                          'Y &/or B array sizes too small'
C
C     Set sightings bias step sizes & their log values. BINC = Bias increment
      NB = NINT(PNBSTP)
      BINC = (PBMAX - PBMIN) / DBLE(NB)
      DO 10 I = 0, NB - 1
        B(I) = PBMIN + (I+0.5D0)*BINC
        RLGB(I) = -LOG(B(I))
   10 CONTINUE
C
C     Set productivity parameter step sizes (midpoints)  [Documentation Eqn 2]
      NR = NINT(PNYSTP)
      YINC = PYMAX / DBLE(NR)
      DO 20 I = 0, NR - 1
        Y(I) = (I + 0.5D0)*YINC
   20 CONTINUE
C     
C
C     Reset the scalar and its maximum
      RSCAL = 0.D0
      MSCAL = -1.D0*1.0D20
C
C     Do two versions of the integral
      DO 55 II = 1,2

      PTOT = 0.D0
      N = 0
      DO 50 I = 0, NR - 1

C       Set R from the productivity parameter, Y   [Documentation Eqn 5]
        R = 1.4184D0*Y(I)

C       Step size for K
        DK = PKSTEP
        D = 1.D0
        RK = RKHI

C       Use function STKSIM to set up the Nth population trajectory
C       i.e. set up the POP array
   30   IF (RK .LE. RKLO .OR. STKSIM(RK, R,POP,CATCH).LE. 0.D0) GOTO 40
          IF (N .GE. MAXSIM) STOP 'ERROR: TOO MANY SIMULATIONS'

C         How much depletion covered?  [Documentation Eqns 6 and 12]
          DD = D - POP(IYEAR)/RK
          D  = POP(IYEAR)/RK
          P  = 0.D0

          IF (DD.GT. 0.D0) THEN

C           Compute the internal catch limit corresponding to D and Y(I)
            QRES(N) = CONTRL(D, Y(I),PLEVEL,PSLOPE)*POP(IYEAR)

C           Calculate deviance [Documentation Eqn 9]
            CALL DEVIAN (SS0,SS1,SS2,SS3, SIGHT,FMATRX,ISYR,IZYR,ZMULT,
     +                   POP,G)

C           Scale likelihood and integrate over values for the bias  
C           [Documentation Eqns 9, 10 & 11]
            DO 35 J = 0, NB - 1
              TERM = -SF * (SS0 + RLGB(J)*(SS1 + RLGB(J)*SS2)+SS3*B(J))
              IF (TERM.GT.MSCAL) MSCAL = TERM
              P = P + EXP(TERM-RSCAL)
   35       CONTINUE

C           Calculate the weight for this point (& total area under likelihood)
C           [Documentation Eqn 11]
            PRES(N) = P*DD
            PTOT = PTOT + PRES(N)

C           Update counter
            N = N + 1

C           Find the next K  [Documentation Eqn 13]
            DK = DK*PDSTEP/DD
            IF (DK .GT. PKSTEP) DK = PKSTEP

          ELSE
C           If DD=0 change the step size only [Documentation Eqn 13]
            DK = PKSTEP
          ENDIF

C         Set the new value of K [Documentation Eqn 14]
          RK = RK/(1.D0 + DK)

        GOTO 30
   40   CONTINUE
   50 CONTINUE

       IF (PTOT .LE. 0.D0 .AND. II .EQ. 2)
     +     STOP 'ERROR: PROB INTEGRATES TO ZERO **'
       RSCAL = MSCAL
C
   55 CONTINUE

C     Sort the QRES and PRES arrays in ascending order of QRES.
C     (Note: N2=N is beyond the end of the array but SORT
C            sorts the values from 0 to N-1 inclusive)
      N2 = N
      CALL SORTC (QRES, PRES, N2)

C     Normalise the relative likelihoods   [Documentation Eqn 15]
      DO 60 I = 0, N-1
        PRES(I) = PRES(I)/PTOT
   60 CONTINUE

C     Extract the desired probability level: the nominal catch limit (NCL) 
C     is the lower PPROB% of the distribution.  
C     First calculate PRES(I), the probability that the NCL is between 
C        QRES(I) & QRES(I+1).  [Documentation Eqn 16]
      P = 0.D0
      DO 70 I = 0, N-1
        P = P + PRES(I)
        IF (P .GT. PPROB) GOTO 80
   70 CONTINUE

C     Interpolate to set the nominal catch limit  [Documentation Eqn 17]
   80 IF (I .GE. N - 1) THEN
        Q = QRES(N-1)
      ELSE
        Q = (QRES(I+1)*(PPROB - P + PRES(I)) + QRES(I)*(P - PPROB))
     1             /PRES(I)
      ENDIF

      CLIMIT = Q

      RETURN
      END

C ***********************************************************************

      FUNCTION CONTRL (D, Y, PLEVEL, PSLOPE)
      DOUBLE PRECISION CONTRL,D,Y,PLEVEL,PSLOPE

C     Catch control law   [Documentation Eqn 7]

      IF (D .LT. PLEVEL) THEN
        CONTRL = 0.D0
      ELSEIF (D .LT. 1.D0) THEN
        CONTRL = PSLOPE*Y*(D - PLEVEL)
      ELSE
        CONTRL = PSLOPE*Y*(1.D0-PLEVEL)
      ENDIF

      END

C ***********************************************************************

      SUBROUTINE DEVIAN (SS0,SS1,SS2,SS3,SIGHT,FMATRX,ISYR,IZYR,ZMULT,
     +                   POP,G)

C     Calculate deviance (-2 log likelihood) in terms of coefficients for
C     the bias and log bias parameters  [Documentation Eqn 9]

      COMMON /MANDAT/ ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER ISTART,IYEAR,NS,NZ
      DOUBLE PRECISION RKLO,RKHI

      DOUBLE PRECISION SS0,SS1,SS2,SS3,SIGHT(0:*),FMATRX(0:*),
     +       ZMULT(0:*),POP(0:*),G(0:*)
      INTEGER ISYR(0:*),IZYR(0:*),N,K,J

      SS0 = 0.D0
      SS1 = 0.D0
      SS2 = 0.D0
      SS3 = 0.D0

      DO 100 N = 0, NS-1

C       [Documentation Eqn 8]
        G(N) = LOG( SIGHT(N)/POP(ISYR(N)) )
        K = N*(N+1)/2
        DO 10 J = 0,N-1
C         1st add non diagonal contributions (which are doubled up)
C         [Documentation Eqn 9]
          SS0 = SS0 + 2.D0 *  G(J) * G(N) * FMATRX(K)
          SS1 = SS1 + 2.D0 * (G(J) + G(N)) * FMATRX(K)
          SS2 = SS2 + FMATRX(K) + FMATRX(K)
          K = K + 1
  10    CONTINUE
C       Now add  diagonal contributions
        SS0 = SS0 + G(N) * G(N) * FMATRX(K)
        SS1 = SS1 + 2.D0 * G(N) * FMATRX(K)
        SS2 = SS2 + FMATRX(K)
 100  CONTINUE

C     Now do the zero estimates   [Documentation Eqn 9a]
      DO 200 N = 0, NZ-1
        SS3 = SS3 + 2.D0*POP(IZYR(N))/ZMULT(N)
 200  CONTINUE
         
      RETURN
      END

C ***********************************************************************

      FUNCTION STKSIM (RK, R,POP,CATCH)

C     Calculate the stock trajectory with parameters RK and R: return 
C     the current stock size     [Documentation Eqns 3 & 4]
C     RK = notional carrying capacity; R = productivity parameter * 1.4184

      COMMON /MANDAT/ ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER ISTART,IYEAR,NS,NZ
      DOUBLE PRECISION RKLO,RKHI

      DOUBLE PRECISION STKSIM,RK,R,CATCH(0:*),POP(0:*),D 
      INTEGER I

      POP(ISTART) = RK
      DO 10 I = ISTART + 1, IYEAR
        D = POP(I-1)/RK
        POP(I) = POP(I-1)*(1.D0 + R*(1.D0 - D*D)) - CATCH(I-1)
        IF (POP(I) .LE. 0.D0) GOTO 20
   10 CONTINUE
      STKSIM = POP(IYEAR)
      RETURN

  20  STKSIM = 0.D0

      END

C ***********************************************************************

      SUBROUTINE SORTC (ARRAY, ARRAY2, N)

C     SORT sorts a pair of arrays in ascending order of the first array
C     (using the heapsort algorithm)

      DOUBLE PRECISION ARRAY(0:*), ARRAY2(0:*), TEMP,TEMP2
      INTEGER N,K,IR,I,J
      IF (N .LT. 2) RETURN
      K = N/2
      IR = N - 1
  10  IF (K .NE. 0) THEN
         K = K - 1
         TEMP = ARRAY(K)
         TEMP2 = ARRAY2(K)
      ELSE
         TEMP = ARRAY(IR)
         TEMP2 = ARRAY2(IR)
         ARRAY(IR) = ARRAY(0)
         ARRAY2(IR) = ARRAY2(0)
         IR = IR - 1
         IF (IR .EQ. 0) THEN
            ARRAY(0) = TEMP
            ARRAY2(0) = TEMP2
            RETURN
         ENDIF
      ENDIF
      I = K
      J = K + K + 1
   20 IF (J .LE. IR) THEN
         IF (J .LT. IR .AND. ARRAY(J) .LT. ARRAY(J + 1)) J = J + 1
         IF (TEMP .LT. ARRAY(J)) THEN
            ARRAY(I) = ARRAY(J)
            ARRAY2(I) = ARRAY2(J)
            I = J
            J = J + I + 1
          ELSE
            J = IR + 1
          ENDIF
          GOTO 20
       ENDIF
       ARRAY(I) = TEMP
       ARRAY2(I) = TEMP2
       GOTO 10
       END

