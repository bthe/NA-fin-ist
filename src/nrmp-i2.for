C     Last change:  CA   25 Sep 2015    2:54 pm
C    INTERFACE between trial control progam & NRMP (=CATCHLIMIT program)

C ***********************************************************************
C *********************************************************************

      DOUBLE PRECISION FUNCTION NCLIMIT (CATCH,SIGHT,FMATRX,ISYR)

C     NCLIMIT = the Nominal catch limit i.e. output of the Catch Limit Calculation

C ****
C **** GLOBAL PARAMETERS USED IN THE CATCH LIMIT CALCULATION (CLC):
C
      COMMON /MANPAR/ PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
      DOUBLE PRECISION PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
C     PPROB   Probability level
C     PYMAX   Maximum value of the productivity parameter (Y) for integration 
C     PBMIN   Minimum bias eg 0.5 means 50% downward bias
C     PBMAX   Maximum bias eg 1.5 means 50% upward bias    
C     PSCALE  Raw deviance scaling factor S = 1/PSCALE**2
C     PHASET  Number of years without surveys before phaseout invoked 
C     PHASEP  Phaseout annual reduction proportion
C     PCYCLE  Maximum number of years before a new CLC =no. of yrs a CLC is valid
C     PLEVEL  Internal protection level     
C
C
C     MANDAT variables are defined so 1 = year of first catch (=CLC timescale)
      COMMON   /MANDAT/  ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER            ISTART,IYEAR,NS,NZ
      REAL(8)             RKLO,RKHI
C     IYEAR  Current year (year for which to set catch limit)
C     NS     Number of non-zero estimates
C     NZ     Number of zero estimates
C     ISTART,RKLO & RKHI are unused - maintained for comparison with Cooke CLC
C
C
      INTEGER MAXYR,MAXEST,MSIZE
      PARAMETER (MAXYR=1001, MAXEST=101)
      PARAMETER (MSIZE=(MAXEST*(MAXEST+1)/2))
      REAL(8) CATCH(1:MAXYR),SIGHT(1:MAXEST),FMATRX(1:MSIZE)
      REAL(8) POP(0:MAXYR),DEVPOP(1:MAXEST),R8W(1:(MAXEST*(MAXEST+7)/2))
      INTEGER ISYR(1:MAXEST)
C
C     CATCH  Catch array, indexed by year
C     SIGHT  Abundance estimate array, indexed by estimate no.
C     FMATRX Information matrix (H) of the log sightings estimates (excluding 
C            zero estimates, stored as a lower triangle
C            FMATRX((i*(i+1))/2 + j = element H(i,j)   NB i & j start at 0
C     ISYR(N) Year of Nth abundance estimate SIGHT(N). N=0,NS-1
C
C ****
C **** Local variables  (defined later)
      REAL(8) IN_PPROB,IN_MU_MIN,IN_MU_MAX,IN_DT_MIN,IN_DT_MAX
      REAL(8) IN_B_MIN, IN_B_MAX, IN_PLEVEL, OUT_QUOTA, accQuota,outDiff
      INTEGER NUM,npRule,INIOUT,IFAIL,INFOLEVEL
C     If IFAIL = 0 on entry to CATCHLIMIT & the routine detects an error
C              it will output an error message and terminate.
C
      IF (NZ.GT.0) THEN
         NCLIMIT=0.0d0
c         write(*,*) 'Zero-estimate, QUOTA=', RAWCL
      ELSE
         NUM = IYEAR-1
         IN_PPROB = PPROB
         IN_MU_MIN = PYMIN
         IN_MU_MAX = PYMAX
         IN_DT_MIN = DTMIN
         IN_DT_MAX = DTMAX
         IN_B_MIN = PBMIN
         IN_B_MAX = PBMAX
         IN_PLEVEL = PLEVEL
         accQuota=ACCTOL
         INIOUT=50
         INFOLEVEL=0
         IFAIL=0
C--------------------------------------------------------------------
         CALL CATCHLIMIT(NUM,NS,CATCH,SIGHT,FMATRX,R8W,ISYR,
     *        IN_PPROB,IN_MU_MIN,IN_MU_MAX,IN_DT_MIN,IN_DT_MAX,
     *        IN_B_MIN,IN_B_MAX,IN_PLEVEL,
     *        OUT_QUOTA,accQuota,outDiff,npRule,
     *        POP,DEVPOP, INFOLEVEL,INIOUT,IFAIL)
         NCLIMIT = OUT_QUOTA
c         write(*,*) 'No Zero-estimate, QUOTA=', RAWCL
      ENDIF

      RETURN
      END

