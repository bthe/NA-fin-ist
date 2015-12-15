C     Last change:  CA   25 Sep 2015    2:51 pm
C
C     PELAGIC1.FOR                          (17 March 1993)
C
C     This module contains most of the subroutines which could be used for
C     other implementation trials.
C
C     The subroutines in this module:
C        SETO     - Set up a logistic ogive
C        SURV     - Set up survivorship profile
C        TRFORM   - Convert to a transitional form
C        MSYPAR   - Calculate A and z from MSYR and MSYL
C        ZCALC    - Used when calculating A and z from MSYR and MSYL
C        SMSYR    - Used to calculate MSYR on mature from MSY/N1+
C
C     -------------------------------------------------------------------
C
      SUBROUTINE SETO (V,SIG,MEAN,MAXAGE)
C
C     SETO computes the ogive corresponding to the given parameters
C
      INTEGER MAXAGE,L
      REAL(8) V(0:MAXAGE),SIG,MEAN
C
C     Knife-edged function
      IF (SIG.LE.0) THEN
       DO 9 L = 0,MAXAGE
        IF (L.LT.MEAN) THEN
         V(L) = 0
        ELSE
         V(L) = 1
        ENDIF
9      CONTINUE
       RETURN
      ENDIF          
C
C     Loop over all ages
      DO 10 L = 1, MAXAGE-2
        IF ((DBLE(L)-MEAN)/SIG.GT.10.d0) THEN
          V(L) = 1.d0
        ELSE IF ((DBLE(L)-MEAN)/SIG.LT.-10.d0) THEN
          V(L) = 0.d0
        ELSE
          V(L) = 1.d0/(1.d0 + EXP(-(DBLE(L)-MEAN)/SIG))
        ENDIF
   10 CONTINUE
      V(0) = 0.d0
      V(MAXAGE-1) = 1.d0
      V(MAXAGE) = 1.d0
C
      RETURN
      END
C
C     -------------------------------------------------------------------
C
      REAL(8) FUNCTION SURV(M,A)
C
C     This function computes survival as function of age
C
      REAL(8) M(3),BETA,ALPHA
      INTEGER A
C
      IF (M(3).NE.-1.d0) THEN
        SURV = EXP(-M(3))
        RETURN
      ENDIF

      BETA = (M(2)-M(1))/16.d0
      ALPHA = M(1) - 4.d0*BETA
      IF (A.LE.4) THEN
        SURV = EXP(-M(1))
      ELSE
        SURV = EXP(-(ALPHA+BETA*A))
      ENDIF
C
      END
C
C     -------------------------------------------------------------------
C
      SUBROUTINE TRFORM(V,MAXAGE)
C
C     Adjust an ogive to transition form, that is so that V(L) =
C     the proportion of animals in a given class at age A-1 which make the
C     transition to a different class age A
C
      INTEGER MAXAGE,L
      REAL(8) V(0:MAXAGE),RM,D
C
      RM   = V(0)
      DO 90 L = 1,MAXAGE
        IF (RM .LT. 1.d0) THEN
          D = RM
          RM = V(L)
          V(L) = (RM - D)/(1.d0 - D)
        ELSE
          RM = 1.d0
          V(L) = 1.d0
        ENDIF
   90 CONTINUE
C
      RETURN
      END
C
c     ---------------------------------------------------------------------
C     -------------------------------------------------------------------
C
      SUBROUTINE MSYPAR(IPNT)

C    Subroutine calculates the resilience A & density dependent exponent Z
C    to give the required MSY at the specified MSYL.  It equates the
C    slope in the balancing per capita birthrates at two levels of
C    fishing mortality just above and below MSY, with the derivative of
C    per capita birth rate at MSYL.  A root finding routine using Brent's
C    method is needed to solve for the density dependent exponent.

      USE ARYPAR
      USE DECLARECOM

      EXTERNAL ZCALC,SURV
      REAL(8) P(0:MXAGE2),UF(2),PARS(4),ZCALC,
     +     FECMSY,SURV,FMAT,DFF,TOL,ZLO,ZHI,F1,F2,PA,RMAT
      REAL(8) RREC,R1PLUS,U1,UM,UMAT(2),U1PLUS(2),MSYLD,TERM
      INTEGER J,IPNT, I,IERR,L
      LOGICAL GOOD
C
C     Set up tolerances
      DFF = 0.0001d0
      TOL = 0.0001d0

C     set RECF =fraction of unrecruited animals of age A which recruit
C     at age A+1, except RECF(0) = fraction recruited of age 0. Note
C     that selectvity equals maturity
      CALL SETO  (RECF,MSIG,MAT1,MAXAGE)

C     Set up maturity ogive FMATUR = proportion mature of age A
      CALL SETO (FMATUR,MSIG,MAT1,MAXAGE)
      DO 7 L = 0,MINMAT
       FMATUR(L) = 0.D0
       RECF(L) = 0.D0
7     CONTINUE

C     Set up natural mortality-at-age array SUR &  calculate the relative
C     population size starting with unity in the 1st age class (A=0)
C     Calculate the relative mature population RMAT (no mature of age 0)
      PA= 1.d0
      RMAT = 0.d0
      RREC = 0.d0
      R1PLUS = 0.d0
      DO 200 L = 0,MAXAGE-1
        SUR(L) = SURV(MORT,L)
        RMAT = RMAT + PA*FMATUR(L)
        RREC = RREC + PA*RECF(L)
        IF (L.GT.0) R1PLUS = R1PLUS + PA
c        PRINT *, 'MSYPAR PA:',PA
        PA = PA*SUR(L)
200   CONTINUE
      SUR(MAXAGE) = SURV(MORT,MAXAGE)
C
C     Adjust for last age class being pooled (and always fully mature)
      PA = PA/(1.d0 - SUR(MAXAGE))
      RMAT = RMAT + PA
      RREC = RREC + PA
      R1PLUS = R1PLUS + PA

C     Set the birth rate so as to give balance at equilibrium
      FEC = 1.d0/RMAT

C     Set up the recruitment ogive in transition form:
C     RECF(A+1) = fraction of unrecruited animals of age A which recruit
C     at age A+1, except RECF(0) = fraction recruited of age 0  {Eqn A4.2}
C     CALL TRFORM (RECF,MAXAGE)
      
      DO 100 J=1,NSTK

C       Estimate FEC at levels of survival which bound MSYR
C       Set two levels of survival after fishing to closely straddle MSY
        UF(1) = 1.d0 - MSYR(J) + DFF*.50d0
        UF(2) = UF(1) - DFF
        DO 30 I = 1,2

C         Set up unrecruited & recruited components relative to # of age 0.
C         Also sum number mature
          P(0) = 1.d0
          FMAT = 0.d0
          U1 = 0.d0
          IF (OPTF.EQ.0) THEN
C          Fishing pattern is 1+
           DO 101 L = 1,MAXAGE-1
            IF (L.NE.1) THEN
             P(L) = SUR(L-1)*UF(I)*P(L-1)
            ELSE
             P(L) = SUR(L-1)*P(L-1)
            ENDIF 
            FMAT = FMAT + P(L)*FMATUR(L)
            U1 = U1 + P(L)
  101      CONTINUE        
          ELSE
C          Fishing pattern is recruited
           DO 102 L = 1,MAXAGE-1
            P(L) = SUR(L-1)*(1-RECF(L-1)*(1-UF(I)))*P(L-1)
            FMAT = FMAT + P(L)*FMATUR(L)
            U1 = U1 + P(L)
  102      CONTINUE
          ENDIF

C         Pooled age class. NB last 2 age classes are all mature & all recruited
          P(MAXAGE) = SUR(MAXAGE-1) * P(MAXAGE-1)*UF(I) /
     +                                  (1.d0 - SUR(MAXAGE)*UF(I))
          FMAT = FMAT + P(MAXAGE)
          U1 = U1 + P(MAXAGE)

C         Save the birth rate which balances with this age-structure
          UF(I) = 1.d0/FMAT
          UMAT(I) = FMAT
          U1PLUS(I) = U1
   30   CONTINUE
        UM = (UMAT(1)+UMAT(2)) * 0.5d0
        U1 = (U1PLUS(1) + U1PLUS(2)) * 0.5

C       Setup MSYLD (JCRM 1:270 eqn 20)
        IF (OPTDD.EQ.OPMSYL) THEN
C        MSYL defined on density dependent component of population
         MSYLD = MSYL
        ELSE IF (OPMSYL.EQ.0 .AND. OPTDD.EQ.1) THEN
C        MSYL on 1+;  DD on mature
         MSYLD = MSYL * R1PLUS * UM * FEC / U1
        ELSE IF (OPMSYL.EQ.1 .AND. OPTDD.EQ.0) THEN
C        MSYL on mature; DD on 1+
         MSYLD = MSYL * U1 / (R1PLUS * UM * FEC)
        ENDIF

C       TERM (1st 2 terms in {} in eqn A7.3)
        IF (OPTDD.EQ.OPTF) THEN
C        MSYR defined on density dependent component of population
         TERM = 0.D0
        ELSE
         TERM = ((U1PLUS(2)-U1PLUS(1))/U1 - (UMAT(2)-UMAT(1))/UM) / DFF
         IF (OPTF.EQ.1 .AND. OPTDD.EQ.0) TERM = -TERM
        ENDIF

C       Place values to be passed into PARS.  PARS2 = Df(Fmsy)/DF / f(Fmsy)-f(0)
        PARS(1) = 1.D0 + TERM*MSYR(J)
        FECMSY = (UF(1) + UF(2)) * .5D0
        PARS(2) = MSYR(J) * (UF(2) - UF(1))/(DFF*(FECMSY - FEC))
        PARS(3) = MSYLD

C       Find density dependent exponent
        ZLO =  -5.d0
        ZHI =   5.d0
        CALL ZBRAC (ZCALC, PARS, 4, ZLO, ZHI, F1, F2, GOOD)
        IF (GOOD) THEN
          CALL ZBRENT (ZCALC, PARS, 4, ZLO, ZHI, F1, F2, TOL, IERR)
        ELSE
          IF (IPNT.GE.0) STOP '*** FAILURE TO SOLVE FOR Z'
        ENDIF

C       Comptute the density dependent exponent & the resilience parameter
        Z(J) = PARS(4)
        A(J) = (FECMSY/FEC - 1.d0)/(1.d0 - MSYL**Z(J))
  100 CONTINUE

      RETURN
      END

C     -----------------------------------------------------------------

      FUNCTION ZCALC (PARS,IERR)

C     The root of this function = density dependent exponent Z
C     PARS(1) = 1 + MSYR(J) (dP/dF/P - dPD/dF/PD)         
C     PARS(2) = MSYR(J) * df(Fmsy)/dF / (f(Fmsy)-f(0))
C     PARS(3) = MSYLD
C     PARS(4) = Z(J)

      DOUBLE PRECISION ZCALC,PARS(4)
      INTEGER IERR

      IERR = 0
      IF (PARS(4) .EQ. 0.D0) THEN
        ZCALC = PARS(1) - PARS(2)*LOG(PARS(3))
      ELSE
        ZCALC = PARS(1) - PARS(2)*(PARS(3)**(-PARS(4))-1.D0) / PARS(4)
      ENDIF
      END



