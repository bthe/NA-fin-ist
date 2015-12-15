C     Last change:  CA   25 Sep 2015    2:54 pm

C *********************************************************************
C     MANNAF.FOR
C
C     MANNAF.FOR applies the RMP to N.Atlantic fin whales
C
C  NAFpar.PAR gives the specifications of the management variant to apply
C  and includes Small Area definitions & whether to cap &/or cascade.
C  Details of the NAFin variants are given at the end of NAFpar.PAR
C
C  If OPTCLC=N the catch limit is set using the Norwegian CLC;
C  if OPTCLC=C the COOKE CLC is used
C
C *********************************************************************
C
      SUBROUTINE ICLC (OPTCLC)
C
C     Initialisation on start of new trial

      USE ARYPAR

      CHARACTER OPTCLC*1
      COMMON /STORE/ IASESS, NCAP, ATYPE, ITYPE
      INTEGER IASESS(I1SYR:IENDYR),NCAP(2),ATYPE(2,MXAR),ITYPE

      COMMON /MANPAR/ PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
      DOUBLE PRECISION PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL

C     Initialise Reporting statistics: NCAP stores no. of capped quotas
      NCAP = 0

C     Open the input file (only read once).

      IF (OPTCLC=='N') THEN
C       Read parameters for the Norwegian version of the CLC
        OPEN (UNIT=21,FILE='CLC-N.PAR')
        READ (21,'(T30,F10.0)') PPROB,PYMIN,PYMAX,DTMIN,DTMAX,
     +                          PBMIN,PBMAX,PSCALE,PHASET,
     +                          PHASEP,PCYCLE,PLEVEL,PSLOPE,ACCTOL
      ELSEIF (OPTCLC=='C') THEN
C       Read parameters for the Cooke version of the CLC
        OPEN (UNIT=21,FILE='CLC-C.PAR')
        READ (21,'(T30,F10.0)') PPROB,  PYMAX,  PNYSTP,PKSTEP,PDSTEP,
     +                          PBMIN,  PBMAX,  PNBSTP,PSCALE,PHASET,
     +                          PHASEP, PCYCLE, PLEVEL,PSLOPE
      ENDIF

      CLOSE (21)
      RETURN
      END
C
C *********************************************************************
C
      SUBROUTINE RCLC ()
C
C     Initialisation on start of each replicate
C
      RETURN
      END

C
C ************************************************************************
C
      SUBROUTINE RUNCLC(CATKA,SGHT,CVEST,IYR,INITYR,CATSM,
     +                  MAPC,OPTCLC,IPROC,NPROC)

C   This subroutine sets the catch quotas for all assessment areas.

C   Output:
C   CATSM(KA)   catch limits for each small area KA in year IYR.

C   Input parameters:
C   OPTCLC      The CLC version (Norwegian or Cooke) to set the catch limit
C   CATKA(IYR,KA) Reported Catch from area KA in IYR (males+females)
C   SGHT(IYR,KA)Absolute abundance estimate for area KA in year IYR.
C               = -1 if no survey performed ie no data available
C   CVEST       the CVs of the sightings estimates, organised in the
C               same way as SGHT.
C   KA          the assessment areas (possible Small, Medium or Combination areas)
C               Areas 1 to NSUBA are subareas; NSUBA+1 to MXAR area are groups of
C               subareas (as defined by MAP in MAIN).
C   IYR         the current year (for which will set catch quotas).
C   INITYR      the first year.
C   NSUBA       Number of sub-areas   (K =1,NSUBA)
C   MAPC(KA,KS) Map defining the Small Areas making up combination or medium areas
C               eg Assessment area 11 is a Combination Area made up of subareas 3 4 & 5
C                  so MAPC(11,KS)=0 0 1 1 1 0 0 0 0 0 0
C   IPROC       The number of the procedure to apply (=1,NPROC) where
C   NPROC       = No. of procedures = 1 in all initial trials.
C               = 2 to use one variant for 10yrs & then phase in to another variant
C   ITYPE       Defines the management variant.  TYPEDS is its description
C               -1: Set catch = 0 AND future bycatch=0
C               0 : Set Catch = 0 for all Small Areas
C               1-8 : Management variants tested in  NAFin trials
C   ATYPE(IPROC,KA) Defines area type & how the catches in Area KA are set
C       Note: AType is not used if ITYPE=-1 or 0
C       0 : Area not assessed. Catch quota = 0
C       1 : A Small Area. Set quota by Small Area - no cascading
C       2 : A Combination Area: set quota & allocate Small Areas by Catch-cascading
C       3 : A Medium Area. Set quota & use to cap Small Areas making up the area.
C       4 : Small area for which catch limit set by cascading
C       5 : Small area for which catch limit set by cascading, but catch not taken
C       9 : Area not assessed. Fixed catch set (used to test effect of given catch level)
C
C   Note: Small Areas are equal to subareas or groups of subareas.
C         MXAR different sub-area groupings are allowed.
C
C   Extract from the RMP specifications:
C   26B) The order in which catch limits are calculated is as follows:
C   (i) the CLA is applied to compute catch limits for Small Areas &/or Medium/
C       Large Areas and Combination Areas as required, with  the associated
C       abundance estimates utilised having the time stamps specified in annotation 20a;
C  (ii) when Catch-cascading is involved the associated catch limit for a
C       Combination Area is distributed amongst the constituent Small Areas;
C (iii) the Phaseout Rule (Section 3.4) is applied to catch limits for Small Areas;
C  (iv) the adjustment for recent sex ratios in the catch (see Section 3.5)
C       is applied to catch limits for Small Areas;
C   (v) the adjustment for other sources of humancaused mortality (Section 3.6)
C       is applied to the catch limits for each Management Area (Small,Medium, Large);
C  (vi) Catch-capping limitations relate to Small Area limits as evaluated at stage (v).
C
C **** GLOBAL PARAMETERS USED IN THE CATCH LIMIT CALCULATION (CLC):
C
      USE ARYPAR

      COMMON /MANPAR/ PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
      DOUBLE PRECISION PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
C
C     Global variables
      REAL(8) CATKA(I1YR:IENDYR,MXAR), CATSM(2,MXAR),
     +       SGHT(I1SYR:IENDYR,MXAR), CVEST(I1SYR:IENDYR,MXAR)
      INTEGER IYR,INITYR,MAPC(MXAR,MXAR),IPROC,NPROC,IL
      CHARACTER OPTCLC*1
C
C     Common data
      COMMON /STORE/ IASESS, NCAP, ATYPE, ITYPE
      INTEGER IASESS(I1SYR:IENDYR),NCAP(2),ATYPE(2,MXAR),ITYPE
C     IASESS(IYR) = 1 if an assessment of the region is done in year IYR
C     NCAP(1) = # of assessments which are capped out of NCAP(2) possibles

C     Local variables
      REAL(8) QUOTA(2,MXAR),S2(MXAR),S3(MXAR),CL,GG,CTOT
      INTEGER IY,KA,KS,ILAST(MXAR),IA
C     ILAST(KA)= Year of the most recent survey in area KA used in an assessment

      SAVE QUOTA
C
C     Reset all catch limits for IPROC to zero
      CATSM(IPROC,1:MXAR) = 0.d0

      IF (ITYPE<=0) RETURN

C     Catch limits for IPROC=1 are not used if IYR>15 & NPROC=2
      IF (NPROC==2 .AND. IPROC==1 .AND. IYR>15) RETURN

C     Set fixed future catch - test option or fixed aboriginal catch - not used here
C      DO KA=1,MXAR
C       IF (ATYPE(IPROC,KA)==9) CATSM(IPROC,KA) = 30.d0
C      ENDDO

C     Use the stored Catch Limit if an assessment is not run this year
      IF (IASESS(IYR)==0) THEN
       DO 5 KA=1,MXAR
        CATSM(IPROC,KA) = QUOTA(IPROC,KA)

C       Find year of last survey in area KA
        DO 4 IY = IYR-1, I1SYR, -1
          IF (SGHT(IY,KA)  /=  -1.0d0) THEN
            ILAST(KA) = IY
            GO TO 5
          ENDIF
    4   CONTINUE
    5  CONTINUE

C      Test for phaseout
       GO TO 600
      ENDIF


C     1st do an assessment for all Small, Medium and Combination Areas
C        except for Small Areas that are part of a Combination Area (ATYPE=4)

      DO 150 KA = 1,MXAR
        IF (ATYPE(IPROC,KA)==0 .OR. ATYPE(IPROC,KA)>=4) GOTO 150
C
C       Call SETCAT to set the catch limit CL and the last survey year IL
        IF (OPTCLC=='N') THEN
         CALL SETCN  (CL,CATKA(I1YR,KA),SGHT(I1SYR,KA),CVEST(I1SYR,KA),
     +               INITYR, IYR,IL)
        ELSEIF (OPTCLC=='C') THEN
         CALL SETCCK (CL,CATKA(I1YR,KA),SGHT(I1SYR,KA), CVEST(I1SYR,KA),
     +               INITYR, IYR,IL)
        ENDIF
        CATSM(IPROC,KA) = CL
        ILAST(KA) = IL
  150 CONTINUE


      DO 400 KA = 1,MXAR

        IF (ATYPE(IPROC,KA)==2 .AND. CATSM(IPROC,KA)>0.d0) THEN

C ******* Catch Cascading (ATYPE = 2) ******************************************
C         For each (non overlapping) Combination Area KA set the Catch Limits
C         for all the Small Areas within it.
C         MAPC(KA,KS) defines the Small Areas making up Combination Area KA
C         CATSM(IPROC,KA) = CL for Combination Area KA was set above (DO 150 loop)

C         Set the relative abundance of each Small Area within Combination Area KA.
          DO 240 KS = 1,MXAR
            S2(KS) = 0.d0
            S3(KS) = 0.d0
            DO 240 IY = I1SYR,IYR-1
              IF (MAPC(KA,KS)==1  .AND. SGHT(IY,KS) > 0.d0) THEN
                GG = 0.9d0**(IYR-IY)/CVEST(IY,KS)**2
                S2(KS) = S2(KS) + LOG(SGHT(IY,KS))*GG
                S3(KS) = S3(KS) + GG
              ENDIF
240       CONTINUE
C
C         Sum the relative abundance indices
          GG = 0.d0
          DO 250 KS = 1,MXAR
           IF (S3(KS)/=0.d0.AND.MAPC(KA,KS)==1) GG=GG+EXP(S2(KS)/S3(KS))
250       CONTINUE
C
C         Set the Small Area catch limits (i.e. do the catch cascading)
C         i.e. split CL between the Small areas making up Combination Area KA
C         (CATSM(IPROC,KS) was initialised to 0 on entry to this routine)
          DO 260 KS = 1,MXAR
            IF (MAPC(KA,KS)/=1 .OR. ATYPE(IPROC,KS)==5) GO TO 260
            IF (S3(KS) > 0.d0)
     +          CATSM(IPROC,KS) = CATSM(IPROC,KA)*EXP(S2(KS)/S3(KS))/GG

C           Set year of the last survey within combination area KA (used to set phaseout)
            ILAST(KS) = ILAST(KA)
260       CONTINUE
C

        ELSE IF (ATYPE(IPROC,KA)==3) THEN

C ******* Catch Capping (ATYPE=3) **********************************************
C         The quotas of the Small Areas within a Medium Area are 'Capped'
C         so the Medium Area limit is not exceeded.  CATSM(IPROC,KA) for
C         the Medium Area KA & its component Small Areas was set above (DO 150 loop)
C         Note: ATYPE=1 for these component Small Areas.

C         Sum the Small Area catch limits within the Medium Area KA
          CTOT = 0.d0
          DO 300 KS = 1,MXAR
            IF (MAPC(KA,KS)==1) CTOT = CTOT + CATSM(IPROC,KS)
  300     CONTINUE
C       
C         Apply capping if necessary (CL = Medium area limit)
          CL = CATSM(IPROC,KA)
          IF (CTOT > CL) THEN
            DO 380 KS = 1,MXAR
             IF (MAPC(KA,KS)==1) CATSM(IPROC,KS)=CATSM(IPROC,KS)*CL/CTOT
  380       CONTINUE
C           Increment the catch capping statistics
            NCAP(1) = NCAP(1) + 1
          ENDIF
          NCAP(2) = NCAP(2) + 1
C
        ENDIF

400   CONTINUE

C     Store all the raw quotas:
      QUOTA(IPROC,1:MXAR) = CATSM(IPROC,1:MXAR)


C *********************************************************************

600   CONTINUE

C *** PHASE OUT

C     Apply Phase out to Small Area catch limits if necessary: call PHOUT
C     ILAST = Year of most recent survey data (on IYR time scale)
      DO 650 KA=1,MXAR
       IA = ATYPE(IPROC,KA)
C  NO PHASEOUT IN THIS VERSION
C       IF ((IA==1 .OR. IA==4).AND. CATSM(IPROC,KA) > 0.d0) THEN
C         CALL PHOUT(CATSM(IPROC,KA),QUOTA(IPROC,KA),ILAST(KA),IYR)
C       ENDIF
       IF (IA>0 .AND. IA<9) WRITE (89,'(6I6,5F10.2)')IYR,IASESS(IYR),
     +      IPROC,IA,KA,ILAST(KA),CATSM(IPROC,KA),QUOTA(IPROC,KA),
     +      SGHT(ILAST(KA),KA)
 650  CONTINUE

      RETURN
      END

C ***********************************************************************
C ***********************************************************************
C
      SUBROUTINE PHOUT (CL,RAWCL,ILAST,IY)
C
      COMMON /MANPAR/ PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
      DOUBLE PRECISION PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
C     PHASET  Number of years without surveys before phaseout invoked 
C     PHASEP  Phaseout annual reduction proportion
C
      DOUBLE PRECISION CL,RAWCL
      INTEGER ILAST, IY
C
C     Phaseout: Reduce catch limit if there is no survey data in the 
C     last PHASET years
      IF (IY .GE. ILAST + PHASET) THEN
        CL = RAWCL * (1.D0 - PHASEP*(IY-ILAST-PHASET))
        IF (CL < 0.D0) CL = 0.D0
      ELSE
        CL = RAWCL
      ENDIF
C
      RETURN
      END
C
C ***********************************************************************
C ***********************************************************************
C
      SUBROUTINE SETCN(CL,CATKA,FSIGHT,CVEST,INITYR,IYR,IL)
C
C     Rearrange catch & sightings arrays into form required by the CLC
C     Call the CLC (catch limit calculation - Norwegian version).
C
      USE ARYPAR

      INTEGER INITYR,IYR,IL
      REAL(8) CL,NCLIMIT
      REAL(8) CATKA(I1YR:IENDYR),FSIGHT(I1SYR:IENDYR),
     +        CVEST(I1SYR:IENDYR)

      INTEGER MAXYR,MAXEST,MSIZE
      PARAMETER (MAXYR=1001, MAXEST=101)
      PARAMETER (MSIZE=(MAXEST*(MAXEST+1)/2))

      REAL(8) CATCH(1:MAXYR),SIGHT(1:MAXEST),FMATRX(1:MSIZE)
      INTEGER ISYR(1:MAXEST)
C     CATCH(I) = Catch in year I
C     SIGHT(N) = Nth Abundance estimate, N=1,NS where NS is the number of surveys
C     FMATRX Information matrix (H) of the log sightings estimates
C            (excluding zero estimates), stored as a lower triangle
C            FMATRX((i*(i+1))/2 + j = element H(i,j)
C     ISYR(N) Year of Nth abundance survey SIGHT(N). N=1,NS
C
C     MANDAT variables are defined on the CLC timescale.  In this case
C            year of first catch = 1
      COMMON   /MANDAT/  ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER            ISTART,IYEAR,NS,NZ
      REAL(8)             RKLO,RKHI
C     IYEAR  Current year (year for which to set catch limit)
C     NS     Number of non-zero estimates
C     NZ     Number of zero estimates
C     ISTART Year of first recorded catch (not used by the Norwegian CLC)
C     RKLO & RKHI Lower & Upper bounds used in integration over K (not used)
C
      COMMON /MANPAR/ PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
      DOUBLE PRECISION PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
C
C     Local variables
      INTEGER N,I, I1
      EXTERNAL NCLIMIT
C
C     Arrange catch and sightings arrays into form required by the CLC
C
C     Set I1 = 1st year to consider = 1st catch yr. Re-scale IYR on CLC timescale
      DO 30 I1 = INITYR,-1
        IF (CATKA(I1) > 0.0d0) GO TO 40
   30 CONTINUE
   40 IF (I1>I1SYR) STOP 'GENDAT ERROR: SIGHTINGS START BEFORE CATCHES'
      IYEAR = IYR - I1+1
C
C     Store catch history for area KA in CATCH.
      CATCH  = 0.d0
      DO 50 I = I1,IYR-1
        CATCH(I-I1+1) = CATKA(I)
   50 CONTINUE
C
C     Set up SIGHT & the information matrix FMATRX: diagonal with elements =1/CV**2
C     Set IL = Year of most recent survey data (on IYR time scale)
      FMATRX = 0.d0
      NS = 0
      NZ = 0
      IL = -1000
      DO 60 I = I1SYR, IYR-1
        IF (FSIGHT(I)  >  0.0d0) THEN
          NS = NS + 1
          SIGHT(NS) = FSIGHT(I)
          N = NS*(NS-1)/2 + NS
          FMATRX(N) = 1.d0/(CVEST(I)*CVEST(I))
          ISYR(NS) = I-I1+1
          IL = I
C
        ELSE IF (FSIGHT(I) .NE. -1.0d0) THEN
          NZ = NZ + 1
          IL = I
        ENDIF
   60 CONTINUE
C
C     Apply the CLC (Norwegian version)
      IF (NS==0) THEN
        CL=0.D0
      ELSE
        CL = NCLIMIT (CATCH,SIGHT,FMATRX,ISYR)
      ENDIF
      RETURN
      END

C ***********************************************************************
C ***********************************************************************
C
      SUBROUTINE SETCCK(CL,CATKA,FSIGHT,CVEST,INITYR,IYR,IL)
C
C     Rearrange catch & sightings arrays into form required by the CLC
C     Call the CLC (catch limit calculation - Cooke version).
C
      USE ARYPAR

      INTEGER INITYR,IYR,IL
      REAL(8) CL,CLIMIT
      REAL(8) CATKA(I1YR:IENDYR),FSIGHT(I1SYR:IENDYR),
     +        CVEST(I1SYR:IENDYR)

      INTEGER MAXYR,MAXEST,MSIZE
      PARAMETER (MAXYR=1001, MAXEST=101)
      PARAMETER (MSIZE=(MAXEST*(MAXEST+1)/2))

      REAL(8) CATCH(0:MAXYR),SIGHT(0:MAXEST),FMATRX(0:MSIZE),
     +       ZERO(0:MAXEST),POP(0:MAXYR), G(0:MAXYR)
      INTEGER ISYR(0:MAXEST), IZYR(0:MAXEST)
C     CATCH(I) = Catch in year I
C     SIGHT(N) = Nth Abundance estimate, N=0,NS-1 where NS = no. of surveys
C     FMATRX Information matrix (H) of the log sightings estimates
C            (excluding zero estimates), stored as a lower triangle
C            FMATRX((i*(i+1))/2 + j = element H(i,j)   NB i & j start at 0
C     ISYR(N) Year of Nth abundance estimate SIGHT(N). N=0,NS-1
C     ZERO(N) Poisson multiplier for Nth zero estimate
C     IZYR(N) Year of Nth zero estimate.  N=0,NZ-1
C     POP(I)  Modelled population size in year I (set in STKSIM)
C     G(I)    (set & used in DEVIAN)

C     MANDAT variables are defined on the CLC timescale.  In this case
C            year of first catch = 1
      COMMON   /MANDAT/  ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER            ISTART,IYEAR,NS,NZ
      REAL(8)             RKLO,RKHI
C     IYEAR  Current year (year for which to set catch limit)
C     NS     Number of non-zero estimates
C     NZ     Number of zero estimates
C     ISTART Year of first recorded catch
C     RKLO & RKHI Lower & Upper bounds used in integration over K
C
      COMMON /MANPAR/ PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
      DOUBLE PRECISION PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
C
C     Local variables
      INTEGER N,I1,ISCALE,I
cx      REAL(8) TOTCAT,PBAR,CV2,CVP

      EXTERNAL NCLIMIT

C     Arrange catch and sightings arrays into form required by the CLC

C     Set year of first catch (ISTART)
C     Set up the first year to consider. Re-scale IYR on CLC timescale (=MANDAT scale)
      ISCALE = 1-INITYR
      CATCH  = 0.d0
      DO 30 I1 = INITYR,-1
        IF (CATKA(I1)  >  0.0d0) GO TO 40
   30 CONTINUE
   40 ISTART = I1 + ISCALE
      IYEAR = IYR + ISCALE
C
C     Store catch history for area KA in CATCH & sum total catch to date
cx      CV2 = 0.d0
cx      PBAR = 0.d0
cx      TOTCAT = 0.D0
      DO 50 I = I1,IYR-1
        CATCH(I+ISCALE) = CATKA(I)
cx        TOTCAT = TOTCAT + CATCH(I+ISCALE)
   50 CONTINUE
C
C     Set up SIGHT & the information matrix FMATRX: diagonal with elements =1/CV**2
C     Set IL = Year of most recent survey data (on IYR time scale)
      FMATRX = 0.d0
      NS = 0
      NZ = 0
      IL = -1000
      DO 60 I = I1SYR, IYR-1
        IF (FSIGHT(I)  >  0.0d0) THEN
          SIGHT(NS) = FSIGHT(I)
          N = NS*(NS+1)/2 + NS
          FMATRX(N) = 1.d0/(CVEST(I)*CVEST(I))
cx          CV2 = CV2 + (CVEST(I)*CVEST(I))
cx          PBAR = PBAR + SIGHT(NS)
          ISYR(NS) = I + ISCALE
          NS = NS + 1
          IL = I
C
        ELSE IF (FSIGHT(I) .NE. -1.0d0) THEN
          STOP 'ERROR: ZERO ESTIMATES SHOULD NOT OCCUR'
          IZYR(NZ) = I + ISCALE
          ZERO(NZ) = CVEST(I)
          NZ = NZ + 1
          IL = I
        ENDIF
   60 CONTINUE

C     No zero estimates are passed to the CLC.
      NZ = 0
C
C     Bound the range for K (RKHI & RKLO)
      RKHI = 1.D7
      RKLO = 0.D0

Cx     Old method to bound the range for K. No longer used but code kept in case!
Cx     Normalise the mean abundance to PBMAX; Calculate the s.e. of the mean CV
cx      CVP = 0.d0
cx      IF (NS > 0)  THEN
cx        CVP  = PSCALE*SQRT(CV2)/DBLE(NS)
cx        PBAR = PBAR / (DBLE(NS)*PBMAX)
cx      ENDIF
cx      RKHI = PBAR*EXP(3.d0*CVP) + TOTCAT
cx      RKHI = MAX (RKHI, 200.d0*TOTCAT)
cx      RKLO = PBAR*EXP(-3.d0*CVP)
C
C     Apply the CLC (Cooke Version)
      CL = CLIMIT(CATCH,SIGHT,FMATRX,ISYR,IZYR,ZERO,POP,G)

      RETURN
      END                                                                                         v
