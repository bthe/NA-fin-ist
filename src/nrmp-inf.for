C     Last change:  CA   12 Jun 2001   10:17 am
C    INTERFACE

C ***********************************************************************
C
      SUBROUTINE RDPARS
C
C     Read the file of input parameters
C
      COMMON /MANPAR/ PPROB,PYMAX,PNYSTP,PKSTEP,PDSTEP,PNBSTP,PBMIN,
     1                PBMAX, PSCALE,PHASET, PHASEP, PCYCLE,PLEVEL,PSLOPE
      DOUBLE PRECISION PPROB,PYMAX,PNYSTP,PKSTEP,PDSTEP,PNBSTP,PBMIN,
     1       PBMAX,PSCALE,PHASET,PHASEP,PCYCLE,PLEVEL,PSLOPE
      COMMON /FILES/ copydat, randomf, clcfile
      character(len=30)  copydat, randomf, clcfile
      INTEGER IFILE
      DATA IFILE /21/
C
C     Open the input file (only read once)
      OPEN (UNIT=IFILE,FILE=clcfile)
C
      READ (IFILE,'(T30,F10.0)') PPROB,  PYMAX,
     +                           PBMIN,  PBMAX,  PSCALE, PHASET,
     +                           PHASEP, PCYCLE, PLEVEL
      CLOSE (IFILE)
      RETURN
      END
C
C
C ***********************************************************************
C ***********************************************************************
C
      SUBROUTINE PHOUT (CL,RAWCL,ILAST,IY)
C
      COMMON /MANPAR/ PPROB,PYMAX,PNYSTP,PKSTEP,PDSTEP,PBMIN,PBMAX,
     1             PNBSTP,PSCALE, PHASET, PHASEP, PCYCLE, PLEVEL, PSLOPE
      DOUBLE PRECISION PPROB,PYMAX,PNYSTP,PKSTEP,PDSTEP,PNBSTP,PBMIN,
     1       PBMAX,PSCALE,PHASET,PHASEP,PCYCLE,PLEVEL,PSLOPE
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
        IF (CL.LT.0.D0) CL = 0.D0
      ELSE
        CL = RAWCL
      ENDIF
C
      RETURN
      END

C *********************************************************************

      SUBROUTINE RSETSS

C     Subroutine initialises parameters at start of each trial

      COMMON /CLCPAR/ IASESS, LSURV, PASESS, POUT
      INTEGER IASESS, LSURV, PASESS, POUT

C     IASESS = Year of last assessment 
C     LSURV =  Year of the most recent survey used in an assessment
C              NB IASESS & LSURV are set to -1000 at the start of each trial 
C     PASESS = 1 is an option allowing assessments to be done if new 
C              survey data arrives.   Otherwise assessments are done 
C              at fixed intervals of PCYCLE years.
C     POUT    Parameter determining whether phase out may be applied if 
C             necessary.  Test for phase out if POUT=1.  (Phase out is
C             not applied to Medium or Large area nominal catch limits) 

      IASESS = -1000
      LSURV  = -1000
      POUT   = 1
      PASESS = 1
      END

C *********************************************************************

      SUBROUTINE MANSS (FCATCH, FSIGHT, CVX, IYR, INITYR, ITEST)

C ****
C **** GLOBAL PARAMETERS USED IN THE CATCH LIMIT CALCULATION (CLC):
C
      COMMON /MANPAR/ PPROB,PYMAX,PNYSTP,PKSTEP,PDSTEP,PNBSTP,PBMIN,
     1                PBMAX, PSCALE,PHASET,PHASEP,PCYCLE,PLEVEL,PSLOPE
      DOUBLE PRECISION PPROB,PYMAX,PNYSTP,PKSTEP,PDSTEP,PNBSTP,PBMIN,
     1       PBMAX,PSCALE,PHASET,PHASEP,PCYCLE,PLEVEL,PSLOPE
C     PPROB   Probability level
C     PYMAX   Maximum value of the productivity parameter (Y) for integration 
C     PBMIN   Minimum bias eg 0.5 means 50% downward bias
C     PBMAX   Maximum bias eg 1.5 means 50% upward bias    
C     PSCALE  Raw deviance scaling factor S = 1/PSCALE**2
C     PHASET  Number of years without surveys before phaseout invoked 
C     PHASEP  Phaseout annual reduction proportion
C     PCYCLE  Maximum number of years before a new CLC 
C             = number of years a CLC is valid 
C     PLEVEL  Internal protection level     
C
C
C     MANDAT: All years are scaled so 0 = a year prior to the 1st recorded catch
      COMMON   /MANDAT/  ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER            ISTART,IYEAR,NS,NZ
      DOUBLE PRECISION   RKLO,RKHI
C     ISTART Year of first recorded catch
C     IYEAR  Current year (year to set catch limit)
C     NS     Number of non-zero estimates
C     NZ     Number of zero estimates
C     RKLO   Lower bound used in integration over K
C     RKHI   Upper bound used in integration over K
C
C
      INTEGER MAXYR,MAXEST,MSIZE
      PARAMETER (MAXYR=1001, MAXEST=101)
      PARAMETER (MSIZE=(MAXEST*(MAXEST+1)/2))
      REAL*8 CATCH(1:MAXYR),SIGHT(1:MAXEST),FMATRX(1:MSIZE),
     +       ZERO(1:MAXEST)
      REAL*8 POP(0:MAXYR), DEVPOP(1:MAXEST),
     +       R8W(1:(MAXEST*(MAXEST+7)/2))
      INTEGER ISYR(1:MAXEST), IZYR(1:MAXEST)
C
C     CATCH  Catch array, indexed by year
C     SIGHT  Abundance estimate array, indexed by estimate no.
C     FMATRX Information matrix (H) of the log sightings estimates (excluding 
C            zero estimates, stored as a lower triangle
C            FMATRX((i*(i+1))/2 + j = element H(i,j)   NB i & j start at 0
C     ISYR(N) Year of Nth abundance estimate SIGHT(N). N=0,NS-1
C     ZERO(N) Poisson multiplier for Nth zero estimate
C     IZYR(N) Year of Nth zero estimate.  N=0,NZ-1
C     POP(I)  Modelled population size in year I (set in STKSIM)
C     G(I)    (set & used in DEVIAN)
C     RAWCL     Nominal catch limit i.e. output of the Catch Limit Calculation
C     CL      The catch limit for the area considered
C     POUT    Parameter determining whether phase out may be applied if 
C             necessary.  Test for phase out if POUT=1.  (Phase out is
C             not applied to Medium or Large area nominal catch limits) 
C
C ****
C **** Parameters passed to this subroutine:
      DOUBLE PRECISION FCATCH(-65:2000), FSIGHT(-1:2000), CVX(-1:2000)
      INTEGER IYR,INITYR,ITEST
C     FCATCH(IYR) Catch in IYR.
C     FSIGHT(IYR) Absolute abundance estimate in IYR1 eg sightings.
C                 Set to -1 if no data available.
C     CVX(IYR)    Estimated CV of the sightings estimate
C     IYR    Current year on scale with 0 = 1st year of management
C     INITYR Year in which first premanagement catch taken (-65 to 0)
C     ITEST  Set to 1 for test calls to MANSS: quota not stored
C
C
      COMMON /CLCPAR/ IASESS, LSURV, PASESS, POUT
      INTEGER IASESS, LSURV, PASESS, POUT
C
C ****
C **** Local variables  (defined later)
      REAL*8 IN_PPROB,IN_MU_MIN,IN_MU_MAX,IN_DT_MIN,IN_DT_MAX
      REAL*8 IN_B_MIN, IN_B_MAX, IN_PLEVEL, OUT_QUOTA
      REAL*8 accQuota, outDiff
      DOUBLE PRECISION RAWCL,STCL
      INTEGER ILAST,NUM,npRule,INIOUT,IFAIL,INFOLEVEL
      SAVE STCL
C
C     Re-scale IYR on MANDAT scale.  
      IYEAR = IYR + 1 - INITYR
C
C     Call GENDAT to rearrange variables into form required for the CLC
      CALL GENDAT (CATCH,SIGHT,FMATRX,ZERO,ISYR,IZYR,
     +             FCATCH,FSIGHT,CVX, INITYR, MAXYR, MSIZE)

C     Set ILAST = Year of most recent survey data
      IF (NS.GT.0) ILAST = ISYR(NS)+ISTART-1
      IF (NZ.GT.0) ILAST = MAX(ILAST, IZYR(NZ)+ISTART-1)

C     Do an assessment if one has not been done in the last PCYCLE-1 years 
C                   or if PASESS=1 & new survey data has arrived 
C     RAWCL    = Nominal Catch Limit for year IYEAR (before any phaseout reduction)
C     STCL   = Stored Nominal Catch Limit = RAWCL set on last assessment 
C     IASESS = Year of last assessment 
C     LSURV =  Year of the most recent survey used in an assessment
C              NB IASESS & LSURV are set to -100 at the start of each trial 
C     PASESS = 1 is an option allowing assessments to be done if new 
C              survey data arrives.   Otherwise assessments are done 
C              at fixed intervals of PCYCLE years.
      IF (IYEAR.GE.PCYCLE+IASESS .OR. ILAST.GE.LSURV+PCYCLE .OR. 
     +     (PASESS.EQ.1.0D0 .AND. ILAST.GT.LSURV) ) THEN
C     
         IF (NZ.GT.0) THEN
            RAWCL=0.0d0
            write(*,*) 'Zero-estimate, QUOTA=', RAWCL
         ELSE
            NUM = IYEAR-ISTART
            IN_PPROB = PPROB
            IN_MU_MIN = 0.0d0
            IN_MU_MAX = PYMAX
            IN_DT_MIN = 0.0d0
            IN_DT_MAX = 1.0d0
            IN_B_MIN = PBMIN
            IN_B_MAX = PBMAX
            IN_PLEVEL = PLEVEL
            accQuota=0.2d0
            INIOUT=50
            INFOLEVEL=0
            IFAIL=0
C-----------------------------------------------------------------------
            CALL CATCHLIMIT(NUM,NS,CATCH,SIGHT,FMATRX,R8W,ISYR,
     *           IN_PPROB,IN_MU_MIN,IN_MU_MAX,IN_DT_MIN,IN_DT_MAX,
     *           IN_B_MIN,IN_B_MAX,IN_PLEVEL,
     *           OUT_QUOTA,accQuota,outDiff,npRule,
     *           POP,DEVPOP,
     *           INFOLEVEL,INIOUT,IFAIL)
            RAWCL = OUT_QUOTA
c            write(*,*) 'No Zero-estimate, QUOTA=', RAWCL
         ENDIF
C     Store the new Nominal Catch Limit (RAWCL), except on test calls
         IF (ITEST.EQ.0) THEN
            STCL = RAWCL
            IASESS = IYEAR
            LSURV = ILAST
         ENDIF
      ELSE
C     No new assessment.  Reset ILAST to last survey used in an assessment
         RAWCL = STCL
         ILAST = LSURV
      ENDIF
      
C     If this catch limit may be subject to phaseout, call PHOUT
      IF (POUT.EQ.1) THEN
        CALL PHOUT (FCATCH(IYR),RAWCL,ILAST,IYEAR)
      ELSE
        FCATCH(IYR) = RAWCL
      ENDIF  

      RETURN
      END
C
C
C ***********************************************************************
C ***********************************************************************
C
      SUBROUTINE GENDAT (CATCH,SIGHT,FMATRX,ZERO,ISYR,IZYR,FCATCH,
     +                   FSIGHT,CVX, INITYR, MAXYR, MSIZE)
C
C     Arrange catch and sightings arrays into form required by the CLC,
C       set year of first catch (ISTART) & bound range for K (RKHI & RKLO)
C
      INTEGER IORGC, IORGS
C     Set IORGC & IORGS = 1st year of catch & of sightings arrays
      PARAMETER (IORGC=-65, IORGS=-1)

      REAL*8 CATCH(*),SIGHT(*),FMATRX(*),ZERO(*)
      DOUBLE PRECISION CV2,CVP,PBAR,FCATCH(IORGC:2000),
     + FSIGHT(IORGS:2000),CVX(IORGS:1000),TOTCAT
      INTEGER ISYR(1:*), IZYR(1:*), INITYR, MAXYR, MSIZE, ISCALE,I,N,IS
C
      COMMON /MANDAT/ ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER ISTART,IYEAR,NS,NZ
      DOUBLE PRECISION RKLO,RKHI
C     MANDAT variables are defined so ISTART = year of first catch
C 
      COMMON /MANPAR/ PPROB,PYMAX,PNYSTP,PKSTEP,PDSTEP,PNBSTP,PBMIN,
     1                PBMAX,PSCALE,PHASET,PHASEP,PCYCLE,PLEVEL,PSLOPE
      DOUBLE PRECISION PPROB,PYMAX,PNYSTP,PKSTEP,PDSTEP,PNBSTP,PBMIN,
     1       PBMAX,PSCALE,PHASET,PHASEP,PCYCLE,PLEVEL,PSLOPE
C
C     SIGHT(N) = Nth survey result, N=0,NS-1 where NS is the number of surveys
C     ISYR(N)  = Year of Nth survey  
C     CATCH(I) = Catch in year I  

      DO 10 N=1,MSIZE
         FMATRX(N) = 0.0D0
   10 CONTINUE
  
      DO 12 N=1,MAXYR
         CATCH(N) = 0.0D0
   12 CONTINUE
      NS = 0
      NZ = 0
C
      CV2 = 0.0D0
      PBAR = 0.0D0

      ISCALE = 1 - INITYR

C     Set up the first year to consider (year of 1st recorded catch)
      DO 30 IS = INITYR,-1
        IF (FCATCH(IS) .GT. 0.D0) GO TO 40
   30 CONTINUE
   40 ISTART = IS + ISCALE

C     Store the catches in CATCH & calculate the total catch to date
      TOTCAT = 0.D0
      DO 50 I = IS,IYEAR-ISCALE-1
        CATCH(I+ISCALE-ISTART+1) = FCATCH(I)
        TOTCAT = TOTCAT + FCATCH(I)
   50 CONTINUE

      DO 20 I = IORGS, IYEAR - (1 + ISCALE)
        IF (FSIGHT(I) .GT. 0.0D0) THEN
          NS = NS + 1
          SIGHT(NS) = FSIGHT(I)
C     Set up the information matrix FMATRX: diagonal with elements =1/CV**2
          N = NS*(NS-1)/2 + NS
          FMATRX(N) = 1./(CVX(I)*CVX(I))
          CV2 = CV2 + (CVX(I)*CVX(I))
          PBAR = PBAR + SIGHT(NS)
          ISYR(NS) = I+ISCALE-ISTART+1
C
        ELSE IF (FSIGHT(I) .NE. -1.0D0) THEN
          NZ = NZ + 1
          IZYR(NZ) = I + ISCALE-ISTART+1
          ZERO(NZ) = CVX(I)
        ENDIF
   20 CONTINUE

C     Normalise the sightings data to PBMAX.  
C     Calculate the s.e. of the mean C.V.
      CVP = 0.0D0
      IF (NS.GT.0)  THEN
        CVP  = PSCALE*SQRT(CV2)/DBLE(NS)
        PBAR = PBAR / (DBLE(NS)*PBMAX)
      ENDIF

C     Bound the range for K
      RKHI = PBAR*EXP(3.0D0*CVP) + TOTCAT
      RKHI = MAX (RKHI, 200.D0*TOTCAT)
      RKLO = PBAR*EXP(-3.0D0*CVP)

      RETURN
      END
