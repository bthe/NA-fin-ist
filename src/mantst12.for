C     Last change:  CA    5 Apr 2011   12:57 pm
      PROGRAM MANTST                                                    

C   MANTST Version 12 (June 2005)  SINGLE STOCK CONTROL PROGRAM
C   Changes :

C   V12   Extended the final year to 2000
C         Variable depletion inputs
C         Added time-varying K and MSYL (age-structured model)
C   V11   Added standard Hitter Fitter population model calculation
C   V10   Tuning parameter calculated within program
C         See alternative line in DATGEN used to give Oman 98 results
C         with different compliers
C   V9:   V9 used for tuning results presented in Oman 98
C         Program changed to double precision
C   V8:   Calls CLC, the implementation program
C         In CLC bias mid points used - differs from COOKEV5
C         New way of specifying & generating CVX & SIGHT
C               (read in expected CVest + process error parameter)
C         Inc. new method of under-reporting historic catch (SC94)
C         Added options for varying K &/or MSYR  (K99,KYEAR etc)
C         Added option not to reseed random number generators
C         CVs set using initial KSURV (only affects varying K trials)
C            Ref: Punt/Butterworth Email 20/3/95
C   V7:   Calls COOKEV5 & sets CVX (replaces SGHT95)
C         Option allowing variable K & A added
C         Variable CVs allowed
C   V6:   Coastal trials removed
C         Unused recruitment & maturity parameters deleted

C     The program generates 'true' population trajectories using a
C     series of catches set by the management procedure.  'Measurements'
C     or estimates of the stock size are generated, incorporating random
C     noise, for input to the management procedure.
C     The management procedure (MANSS) is a self-contained module called
C     to set catch limits using the population estimates & past catches.

C     Note: Define IYEAR=0 as first year of management.  
C     NTRIAL simulations are run for each case, each lasting NYEAR years.


C PARAMETERS -----------------------------------------------------------

C  Population parameters, common to all population models:
C     PTRUE(IYR)  Mature population size of the true stock at start of year
C                 IYR (prior to CATCHM(IYR) being set & removed)
C     PSURV(IYR)  Size of surveyed population in IYR = 1+ population
C     PSURV1      Stores PSURV from year -1 = 1+ population in year -1
C     CATCH(IYR)  True Catch in year IYR 
C     CATCHM(IYR) Catch in IYR, as passed to MANSS. The premanagement 
C             catch is in error if OPTC>=1. Set by MANSS from year 0.
C     SIGHT(IYR)  Absolute abundance estimate, in IYR  eg sightings. 
C             Set to -1 if no data available that year.  Passed to MANSS.
C     CVX(IYR) Estimated CV of the sightings estimate
C            Set to -1 if no data available.  Passed to MANSS.
C     NTRIAL Number of trials. Read in. 1-400 if ITUNE=0 or 1-100000 if ITUNE=1
C     NYEAR  Number of years of management.  Read in. 1-100
C     NPCAT  Number of years of the (constant) premanagement catch.
C            Read in or if OPTRAN=1 a random value from U[15,45] used.
C     NPPROT Number of years of premanagement protection. Read in.
C     IPPROT First year of premanagement protection period = -NPPROT
C     INITYR Year in which first premanagement catch taken (-65 to 0)
C            A constant catch (=1 unit) is taken from INITYR to IPPROT-1
C            followed by 0 catch from years IPPROT to -1, if any. 
C     K1     Pristine mature population size. Set in SETK(A)
C            In cases when K varies during management, K1 = premanagement value
C     K1P    Pristine 1+ population size. Set in SETK(A)
C     KSURV  Pristine 1+ population in IYR. Set in SETK(A) 
C            Updated in DATGEN if K is time dependent ie if K99.ne.0
C     MSYL   MSY Level.  Read in.
C     MSYR1  Initial MSY rate used to set A. Read in  or  if OPTRAN=1 a
C            random value from U[.001,.05] is generated for each trial.
C            If OPTMOD=4 then the value read in = MSY/1+ population at 
C            MSYL, and the true MSYR is calculated in SMSYR.
C     Z      Density dependent exponent, calculated from MSYL.
C     A1     Initial resilience parameter, calculated from MSYR1.
C     DEPL   Depletion level in year 0, =PTRUE(0)/K. Used to set K.
C            Read in or if OPTRAN=1 a random value from U[.01,.99] used.
C     PROBE(IYR) Set if an epidemic occurs in year IYR.  Set in RESET.
C     ERATE  Rate at which epidemics occur.   AGE STRUCTURED MODEL ONLY
C            If epidemics occur, SUR is reduced after year 0 
C     CATERR Reported historic catch / true catch  = 1 in base case trials
C     ISEED1-5 Seeds for random number generators:
C            1: Set SIGHT    2: Set epidemic years if ERATE>0
C            3: Set CVX      4: Set random initial parameters if OPTRAN=1
C            5: Set survey interval if IFREQ=-1    

C P-T model parameters: (OPTMOD<2)
C     M      Mortality rate.  Read in.
C     S      Annual survival rate
C     TM     Age of maturity.  Read in
C     TLAG   Time delay in depletion density dependent response.  Usually = TM
C     K(IYR) Pristine mature population in IYR, set in SETK. = K1 if K99=0
C     K99    If K99>0 then there is a linear change in K
C                     from  K1 in year KYEAR  to  K1*K99 in year NYEAR-1  
C            If K99=-1: K cyclic, starting min;  K99=-2: K cyclic, starting max
C     KYEAR  If K99>0 then K varies linearly from years KYEAR to NYEAR-1  
C     A(IYR) Resilience parameter in IYR. = A1 unless MSYR99>0 or ISTEP>0
C     ISTEP  MSYR changes every ISTEP years, P-T model only. Usually=0.
C     MSYR99 If MSYR99 > 0 then there is a linear change in MSYR
C                      from  MSYR1 in year MSYRYR  to   MSYR1*MSYR99 in NYEAR-1 
C     MSYRYR If MSYR99>0 then MSYR varies linearly from year MSYRYR to NYEAR-1 
C     SREC(IYR) Stores number of recruits in year IYR+TM+1
C            = #of births in IYR surviving to age of recruitment/maturity 

C Survey estimate parameters:
C     BIAS   Bias in absolute abundance estimates 
C     BIAS0  Initial bias in absolute abundance estimates 
C     BINC   Increment added to BIAS every year if OPTB > 0
C     IFREQ  Frequency with which absolute abundance estimates made.
C            1st survey made in year -1. If IFREQ=5 the 2nd is in year 4
C            If IFREQ=-1, time to next estimate generated from U[1,9]
C            If OPTSUR=1 survey costs taken into account & IFREQ ignored
C     ISIGHT Set to year in which next sightings estimate is to be made
C     ENDSUR Year after which no surveys take place.  Read in.
C     CV1EST value of CV(est) when P=0.6K. CV(est) = CV of CVX.  Read in.
C     IYRCV  Year (if any) when CV(est) changes from CV1EST to CV2EST. Read in.
C     DOFMIN Minimum number of degrees of freedom in chi-square distribution
C     CV2EST 2nd value of CV(est): used if CV changes i.e. if IYRCV<NYEAR
C     ETA    Process error parameter.  If ETA = 0 there is no process error. 

C     IYR    Current year
C     N      Current trial number  (1,NTRIAL)
C     IOUT   Output file
C     IN     Input file
C     REF    Reference number of run

C OPTIONS --------------------------------------------------------------
C     OPTRAN Option determining use of random parameters and seeds:
C            0: No random parameters;  Use a new seed for each replicate
C            1: MSYR1, NPCAT (& INITYR) & DEPL generated for each trial.  
C            2: The random number generators are NOT reseeded (factorial expt.) 
C               IRAN defines the inital seeds.
C            3: New seeds are generated (tuning trial)
C     OPTB   Option controlling bias   0: Constant
C            1: Linear, doubles        2: Linear, decreases by 1/3 eg 1.5 -> 1.0
C     OPTC   Defines reported premanagement catch (i.e. as reported to 
C            the management procedure).
C            0: True catch (CATERR=1)     
C            1: Catch in error; true catch & depletion = base case values
C            2: Catch in error; reported catch & PTRUE(0) = base case values
C     OPTMOD Sets the population model 
C            0: Standard Pella Tomlinson
C            1: P-T with maximum recruitment limitation  
C            2: Tent Model  (constant MSYR only)         
C            3: Age structured, maturity = recruitment [base case]
C            4: Age structured, differing ages of maturity & recruitment
C               with MSY / 1+ population at MSYL = base case ratio
C            5: Age structured with MSYL, MSYR and density-dependent
C               components entered
C *** CARE: See different BIRTHS eqn used with different values of OPTMOD ********
C     OPTDPL Specifies whether depletion changes among simulations (not randomly)
C            0: Constant depletion
C            1: Depletion read-in
C     OPTDET Set on deterministic run    0: stochastic   1:deterministic
C     OPTSUR Survey control option       0: none
C            1: Survey costs taken into account
C     ITUNE  0: Standard run  1: Tuning run:

C DEFINITIONS ----------------------------------------------------------

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      COMMON /PTVRS/ K,A,SREC,S,TM,TLAG, K99,MSYR99,KYEAR,MSYRYR
      DOUBLE PRECISION K(-65:2000),A(-65:2000),SREC(-65:2000),S, K99,
     +  MSYR99
      INTEGER TM, TLAG, KYEAR, MSYRYR

      COMMON / TIMEVAR/ AT(-65:2000),KT(-65:2000),ZT(-65:2000),
     +                  K99AEP,KYEARAEP,MSYR99AEP,MSYRYRAEP,ISTEPAEP
      DOUBLE PRECISION AT,KT,ZT,K99AEP,MSYR99AEP
      INTEGER KYEARAEP,ISTEPAEP,MSYRYRAEP

      COMMON /AGEPAR/ MAT1,MSIG,REC1,RSIG,FEC,SUR(0:60),UNREC(0:60),
     +                REC(0:60),RECF(0:60),FMATUR(0:60),MAXAGE,MINMAT
      DOUBLE PRECISION MAT1,MSIG,REC1,RSIG,FEC,SUR,UNREC,REC,RECF,FMATUR
      INTEGER MAXAGE,MINMAT

      COMMON /DATVRS/ BIAS,CV1EST,CV2EST,DOFMIN,ETA,IFREQ,ISIGHT,
     +                ENDSUR,OPTDET,OPTSUR,IYRCV
      DOUBLE PRECISION BIAS,CV1EST,CV2EST,DOFMIN,ETA
      INTEGER IFREQ,ISIGHT,ENDSUR,OPTDET,OPTSUR,IYRCV

      COMMON /RANPAR/ ISEED1,ISEED2,ISEED3,ISEED4,ISEED5,SSEED1,SSEED3
      INTEGER ISEED1, ISEED2, ISEED3, ISEED4, ISEED5,SSEED1,SSEED3

      COMMON /HITFIT/ OPTDD,OPTF,OPTMSYL,FMSY,SMSYL,RFTOT0,RFEXP0,
     +       FMSYL,FECMSY,MSYLT1,MSYLT0,MSYLT2,AMSYR1,AMSYL1,
     +       AMSYR0,AMSYL0,AMSYR2,AMSYL2
      DOUBLE PRECISION FMSY,SMSYL,RFTOT0,RFEXP0,FMSYL,FECMSY
      DOUBLE PRECISION AMSYR0,AMSYL0,AMSYR1,AMSYL1,AMSYR2,AMSYL2
      DOUBLE PRECISION MSYLT0,MSYLT1,MSYLT2
      INTEGER OPTDD,OPTF,OPTMSYL

      DOUBLE PRECISION CATCHM(-65:2000),SIGHT(-1:2000),C,CVX(-1:2000),
     +    CATCH(-65:2000),M,DEPL,BINC,SETZ,MORT1,MORT2,BIAS0,PFIN,
     +    DEPLVALS(400),DEPS(6)
      INTEGER NTRIAL,INITYR,NYEAR,NPCAT,NPPROT,IPPROT,ISTEP,OPTRAN,OPTB,
     +     IYR,I,N,IOUT,IN,IN2,ISEEDS(400,5),IRAN,NZERO,ITUNE,
     +     MORTIP,OPTDEPL,DEPLCOL
      CHARACTER REF*10,DESC*50,PROG*8,DEPLFILE*12, DEPLREF*9
      DATA IOUT/8/, IN/7/, IN2/10/, NZERO/0/, IN3/11/
      
      COMMON /FILES/ copydat, randomf, clcfile, restest, manlog
      character(len=30)  copydat, randomf, clcfile, restest, manlog



C     read in file names from the command line
      integer :: num_args, ix
      character(len=30), dimension(:), allocatable :: args
      
      num_args = command_argument_count()
      allocate(args(num_args))  ! I've omitted checking the return status of the allocation 
      
c     Assign default values to 
      copydat = 'copy.dat'
      randomf = 'random.num'
      clcfile = 'CLC-N.PAR' 
      restest = 'restest'
      manlog = 'manlog'
      if (num_args>0) then
      ix = 1  
      do while(ix < num_args)
         call get_command_argument(ix,args(ix))
         ix = ix + 1
         if(ix == num_args) continue
         select case(adjustl(args(ix-1)))
         case("-main")
            call get_command_argument(ix,args(ix))
            copydat = args(ix)
         case("-rnd")
            call get_command_argument(ix,args(ix))
            randomf = args(ix)
         case("-clc")
            call get_command_argument(ix,args(ix))
            clcfile = args(ix)
         case("-res")
            call get_command_argument(ix,args(ix))
            restest = args(ix)
         case("-log")
            call get_command_argument(ix,args(ix))
            manlog = args(ix)
         end select
      end do
      end if

      OPEN (1,FILE=manlog)
      OPEN (IOUT,FILE=restest)
      OPEN (IN,FILE=copydat)
      OPEN (IN2,FILE=randomf)


C INITIALISATION -------------------------------------------------------
C     Read in the true population parameters and option settings 
C     Check parameters are within allowed range
      READ (IN,'(T37,A /A/)') REF,DESC      
      READ (IN,'(45X,I3,I4)') OPTRAN,IRAN
      READ (IN,'((45X,I3))') OPTB,OPTC,OPTMOD,OPTDET,OPTSUR
      READ (IN,'((43X,I5))') NTRIAL
      if (ntrial.eq.99000) ntrial = 600000
      READ (IN,'((44X,I4))') NYEAR
      READ (IN,'((45X,I3))') NPCAT, NPPROT
      IPPROT = -NPPROT
      IF (NYEAR.LT.0  .OR. NYEAR.GT.2000) STOP 'NYEAR out of range'
      IF (NPCAT.LT.2  .OR. NPCAT.GT.65)  STOP 'NPCAT out of range'
      IF (NPPROT.LT.0 .OR. NPPROT.GT.63) STOP 'NPPROT out of range'
      READ (IN,'((45X,F10.0))') MSYL,MSYR1,DEPL
      IF (MSYL .LT.0.D0 .OR. MSYL.GE.0.95D0) STOP 'MSYL out of range'
      IF (MSYR1.LT.0.D0 .OR. MSYR1.GE.1.D0) STOP 'MSYR1 out of range'
      READ (IN,'((45X,F8.0))') K99,MSYR99
      READ (IN,'((45X,I3))') ISTEP
      IF (K99.GT.0.D0) READ (IN,'((45X,I3))')  KYEAR
      IF (MSYR99.GT.0.D0) READ (IN,'((45X,I3))')  MSYRYR
      IF (KYEAR .LT. 0) STOP 'KYEAR out of range'
      IF (MSYRYR.LT. 0) STOP 'MSYRYR out of range'
      IF (OPTMOD.LT.3) THEN
        READ (IN,'((45X,I3))') TM,TLAG
        READ (IN,'((45X,F8.0))')  M
        ERATE = 0.D0
      ELSE
        READ (IN,'((45X,F8.0))') MAT1,MSIG,REC1,RSIG,MORT1,MORT2
        READ (IN,'(/(45X,I3))')  MORTIP,MAXAGE,MINMAT
        WRITE(*,*) MINMAT
        READ (IN,'((45X,F8.0))') ERATE
        IF (MAXAGE.GE.60) STOP 'MAXAGE out of range'
        READ (IN,*)
        READ (IN,'((45X,I3))') OPTF,OPTMSYL,OPTDD
      ENDIF
      K99AEP = K99
      MSYR99AEP = MSYR99
      KYEARAEP = KYEAR
      MSYRYRAEP = MSYRYR
      ISTEPAEP = ISTEP
      READ (IN,'((45X,I3))') IFREQ, ENDSUR,IYRCV
      READ (IN,'((45X,F8.2))') BIAS0, CV1EST, ETA, DOFMIN
      IF (IYRCV.LT.NYEAR) READ (IN,'((45X,F8.2))') CV2EST
      IF (OPTC.GE.1) READ (IN,'((45X,F8.2))') CATERR
      READ (IN,'(45X,I3)') OPTDEPL
      IF (OPTDEPL.EQ.1) THEN
       READ (IN,'(45x,A)') DEPLFILE
       OPEN (IN3,FILE=DEPLFILE)
       READ (IN,'(45X,I3)') DEPLCOL
       READ(IN3,'(A)') DEPLREF
       DO 109 N = 1,NTRIAL
        READ(IN3,*) (DEPS(I),I=1,6)
        DEPLVALS(N) = DEPS(DEPLCOL)
        write(*,*) N, DEPLCOL, DEPLVALS(N)
109    CONTINUE        
      ELSE
       DO 110 N = 1,NTRIAL
        DEPLVALS(N) = DEPL
110    CONTINUE        
      ENDIF 
      CLOSE (IN)

C     Read in random number seeds for all trials : store in ISEEDS
      READ (IN2,'((5I8))') ((ISEEDS(N,I),I=1,5),N=1,MIN(NTRIAL,400))
      CLOSE (IN2)

C     Call EXTRA to reset any variables if necessary, and open output files. 
C     In zero catch trials it sets NTRIAL=1 
      CALL EXTRA (NTRIAL,OPTRAN,ERATE,ENDSUR,IOUT,PROG,ITUNE,OPTDEPL)
      WRITE(1,'(/2(2X,A)//A)')REF,DESC
      IF (ITUNE.gt.0) OPEN  (31,FILE='TUNEFULL',status='unknown')

C     Initialise variables:

C     Don't do surveys in NYEAR-1 as they will never be used
      IF (ENDSUR.GT.NYEAR-2) ENDSUR = NYEAR-2

C     Set Z, the density dependent exponent using MSYL.
      IF (OPTMOD.LT.3) THEN
        S = EXP(-M)
        Z = SETZ(MSYL)
      ENDIF

C     Set BINC, the annual bias increment. (1): Doubled  (2): Reduced by 1/3 
      IF (OPTB.EQ.1) THEN
        BINC = BIAS0 / DBLE(NYEAR-1)        
      ELSE IF (OPTB.EQ.2) THEN
        BINC = -BIAS0 / DBLE(3*NYEAR-3)   
      ENDIF

C     Call RDPARS to initialise management routine
      CALL RDPARS


C TRIALS BEGIN ---------------------------------------------------------

      DO 200 N = 1,NTRIAL

C       Reseed the random number generators & print seeds, set BIAS & PROBE 
        CALL RESET (ISEEDS,N,NYEAR,BIAS0,BIAS,OPTRAN,IRAN)

C       Set up true population on 1st trial & on random parameter trials
        IF (OPTDEPL.EQ.1. OR. OPTRAN.EQ.1 .OR. N.EQ.1) THEN

C         SETPAR generates MSYR1,NPCAT & DEPL in random parameter trials
          IF (OPTRAN.EQ.1) CALL SETPAR (MSYR1,DEPLVALS(N),NPCAT)
          INITYR = IPPROT - NPCAT
          IF (INITYR.LT.-65 .OR. INITYR.GT.0) STOP 'INITYR out of range'

          IF (OPTMOD.LT.3) THEN
C           P-T or Tent model:
C           SETA sets the resilience array A using MSYR1, for all years
C           SETK finds the carrying capacity giving a depletion DEPL in 
C           year 0 & also sets PTRUE & CATCH prior to management.
            IF (OPTMOD.NE.2) CALL SETA (INITYR,NYEAR,ISTEP)
            CALL SETK (DEPLVALS(N),CATCH,IPPROT,INITYR,NYEAR)
            MSYLT1 = 0.6
            MSYLT0 = MSYL
            MSYLT1 = MSYL
            MSYLT2 = MSYL
            AMSYR0 = MSYR1*100
            AMSYR1 = MSYR1*100
            AMSYR2 = MSYR1*100
            AMSYL0 = MSYL
            AMSYL1 = MSYL
            AMSYL2 = MSYL
          ELSE
C           Age structured model.  
C           SETKA sets K,PTRUE & CATCH
            CALL SETKA (DEPLVALS(N),CATCH,IPPROT,INITYR,0,MORT1,MORT2,
     +       MORTIP,NYEAR)
          ENDIF

C         Modify array CATCHM if necessary for passing to MANSS
          CALL PMCAT (CATCHM,IPPROT,INITYR,CATCH,OPTC,CATERR)

C         Print parameter list & check values are within allowed range.
          IF (N.EQ.1) CALL PNTOUT(IOUT,M,ISTEP,NPCAT,NPPROT,DEPL,NTRIAL,
     +              NYEAR,REF,MORT1,MORT2,OPTRAN,OPTB,IRAN,DESC,PROG,
     +              OPTDEPL,INITYR,ITUNE,DEPLREF,DEPLCOL)

        ELSE IF (OPTMOD.GE.3) THEN
C         Reset the age structured population
          CALL SETKA (DEPLVALS(N),CATCH,IPPROT,INITYR,1,MORT1,MORT2,
     +       MORTIP,NYEAR)
        ENDIF

C       Print out new parameters
        IF (ITUNE.EQ.0) WRITE (IOUT,'(/A6,I4,A10,6I8)') 'Trial:',N,
     +            'Seeds:',ISEED1,ISEED2,ISEED3,ISEED4,ISEED5
        WRITE(IOUT,'(A15,20F8.5)') 'MSYL Outputs:',MSYLT1,MSYLT0,
     +        MSYLT2,AMSYR1,AMSYL1,AMSYR0,AMSYL0,AMSYR2,AMSYL2
        IF (OPTRAN.EQ.1.OR.OPTDEPL.EQ.1) THEN
          WRITE (IOUT,'(A20,I8,8F14.7)')  'New Parameters:',
     +         NPCAT,DEPLVALS(N),MSYR1,K1,A1,Z
        ELSE
          WRITE (IOUT,'()')
        ENDIF

C       Set the premanagement survey estimate, SIGHT(-1), & its CV
        PSURV(-1) = PSURV1
        IYR = -1
        ISIGHT = -1
        CALL DATGEN (SIGHT,CVX,IYR,CATCHM,INITYR,NYEAR,N,NZERO)

C       Reinitialise the management routine
        CALL RSETSS
C
        PRINT '('' '',A,2X,A,I8)', REF,DESC,N
        DO 100 IYR = 0,NYEAR-1

c        PRINT '(''+'',A,2X,A,2I5)', REF,DESC,N,IYR
C         Call MANSS (which contains the catch limit algorithm) to
C         set CATCHM(IYR).  SIGHT(IYR-1) was set in the last call to DATGEN
          CALL MANSS (CATCHM,SIGHT,CVX,IYR,INITYR,0)

C         Advance stock to next year, ie remove catch & set PTRUE(IYR+1)
C         If catch quota C > recruited population then STKUPD resets C
          IF (CATCHM(IYR).LT.0.D0) CATCHM(IYR) = 0.D0
          C = CATCHM(IYR)
          IF (OPTMOD.LT.3) THEN
            CALL STKUPD (IYR,INITYR,C)
          ELSE
            CALL STKUPA (IYR,C)
          ENDIF

C         Reset CATCHM to actual catch taken and also
C         store it in CATCH which is NOT available to MANAGE - for security!
          CATCHM(IYR) = C
          CATCH(IYR)  = C

C         Set estimates of abundance, SIGHT(IYR), & its CV + year of next survey
          CALL DATGEN (SIGHT,CVX,IYR,CATCHM,INITYR,NYEAR,N,NZERO)

C         Update BIAS if it is time dependent
          IF (OPTB.GT.0) BIAS = BIAS + BINC

  100   CONTINUE

C       Call REPORT to print out PTRUE & CATCH in year IYR
        IF (ITUNE.EQ.0) THEN
          CALL REPORT (IOUT,CATCH,IYR)
        ELSE
          PFIN = PTRUE(NYEAR)/K1
          OPEN  (31,FILE='TUNEFULL',ACCESS='APPEND')
          WRITE (31,*) N,PFIN
          CLOSE (31)
        END IF

C        IF (N.LE.2) WRITE(1,'(A3,I3,4(A8,F9.3),2(T80,''SGT & CV in yr'',
C     +      I3,2F9.3/))')'N:',N,'K:',K1,'P0:',PTRUE(0),'P99:',PTRUE(99),
C     +      'P0/K:',PTRUE(0)/K1,(I,SIGHT(I),CVX(I),I=-1,IFREQ-1,IFREQ)

  200 CONTINUE

C     Print out the number of zero survey estimates which have occurred
C     Zero estimates are rare & are treated as a special case by the CLC
      WRITE (IOUT,'(/ A,I5)') ' NZERO:',NZERO
      CLOSE (IOUT)
      CLOSE (1)
      STOP
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE DATGEN (SIGHT,CVX,IYR,CATCHM,INITYR,NYEAR,N,NZERO)

C     DATGEN sets absolute abundance estimate SIGHT(IYR) adding 
C     random noise to the surveyed population PSURV (see equations 
C     in Appendix 2 of IWC/45/4 Annex D).  Also set next survey year.

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      COMMON /DATVRS/ BIAS,CV1EST,CV2EST,DOFMIN,ETA,IFREQ,ISIGHT,
     +                ENDSUR,OPTDET,OPTSUR,IYRCV
      DOUBLE PRECISION BIAS,CV1EST,CV2EST,DOFMIN,ETA
      INTEGER IFREQ,ISIGHT,ENDSUR,OPTDET,OPTSUR,IYRCV

      COMMON /RANPAR/ ISEED1,ISEED2,ISEED3,ISEED4,ISEED5,SSEED1,SSEED3
      INTEGER ISEED1, ISEED2, ISEED3, ISEED4, ISEED5,SSEED1,SSEED3

      EXTERNAL RAN1,GAMMA,XNORM,POISSN
      DOUBLE PRECISION RAN1,GAMMA,XNORM,POISSN,SIGHT(-1:2000),RANNO,
     +   CVX(-1:2000),CATCHM(-65:2000),SIG2,EPS,CATCH1,CATCH2,W,CHISQ,
     +   CVESQ,BKCON,SIG,DOF,ACON,BCON,THTASQ,ALPHSQ,BETASQ,ASQ,BSQ,
     +   DCON,MA1(55),MA3(55),MA5(55),X
      INTEGER IYR,INITYR,NYEAR,NZERO,N,ISEED,
     +        INEXT1,INXTP1,INEXT3,INXTP3,INEXT5,INXTP5

      SAVE BKCON,SIG,DOF,ACON,BCON,
     +     MA1,MA3,MA5,INEXT1,INXTP1,INEXT3,INXTP3,INEXT5,INXTP5
      PARAMETER (ASQ=0.02D0,BSQ=0.012D0,DCON=2.9D0)

C     On 1st call to routine or if CV or K change, set the constants
      IF (IYR.LT.0 .OR. IYR.EQ.IYRCV) THEN   
        CVESQ = CV1EST*CV1EST
        IF (IYR.GE.IYRCV) CVESQ = CV2EST*CV2EST

C       CV(est)**2 = THTASQ (ASQ + BSQ K/P)    Use this eqn to give THTASQ
C       CV(true)**2= ALPHSQ + BETASQ K/P where ALPHSQ = THTASQ*ASQ + ETA*.1
C                                              BETASQ = THTASQ*BSQ + ETA*.013
        THTASQ = CVESQ / (ASQ + BSQ/0.6D0)
        ALPHSQ = THTASQ * ASQ + ETA*0.1D0
        BETASQ = THTASQ * BSQ + ETA*0.013D0
        BKCON = BETASQ * KSURV
        ACON = THTASQ * ASQ
        BCON = THTASQ * BSQ / BETASQ
C       SIG = sd of lognormal component: used to set EPS
        SIG = SQRT(LOG(1.D0 + ALPHSQ))

C       Set DOF = no. of degrees of freedom = 2.9 / CV(true)**2 - 1,
C                 with minimum value DOFMIN.  CV(true) is value at P=K 
        DOF = ANINT (DCON / (ALPHSQ+BETASQ)   - 1.D0)
        IF (DOF.LT.DOFMIN) DOF = DOFMIN

        IF (N.EQ.1 .AND. (IYR.EQ.-1.OR.IYR.EQ.IYRCV))  WRITE(1,1)
     +       SQRT(ALPHSQ+BETASQ/.6),SQRT(ACON+THTASQ*BSQ/.6),DOF,THTASQ
    1   FORMAT(' CV(tru) @P=.6K',F8.4,7X,'CV(est) @P=.6K',F8.4,2F10.3/)
      ENDIF


C     If OPTSUR=1 test to see if a survey will be performed this year.
      IF (OPTSUR.EQ.1 .AND. IYR.GE.0) THEN
C       Apply the management procedure without any new data
        SIGHT(IYR)  = -1.D0
        CVX(IYR) = -1.D0
        CALL MANSS (CATCHM,SIGHT,CVX,IYR+1,INITYR,1)
        CATCH1 = CATCHM(IYR+1)

C       Apply management procedure USING the true new data
        SIGHT(IYR) = BIAS * PSURV(IYR)
        CVX(IYR) = SQRT( LOG(1.D0 + ACON + BCON*BKCON/PSURV(IYR)))
        CALL MANSS (CATCHM,SIGHT,CVX,IYR+1,INITYR,1)
        CATCH2 = CATCHM(IYR+1)

C       If CATCH2 > CATCH1 then take the survey
        IF (CATCH2.GE.CATCH1) THEN
          ISIGHT = IYR
        ELSE
          ISIGHT = IYR+1
        ENDIF
      ENDIF


C     Set EPS, W & SEEDS to update random generators even if no survey
C     XNORM(SIG,MEAN,ISEEDi,...) produces a random number from a normal 
C     distribution N [MEAN,SIG**2] using the ith random generator.
C          (if MEAN = 0,  XNORM (SIG,0) = SIG * XNORM (1,0) )

      EPS = SIG * XNORM (1.D0,0.D0,ISEED1,MA1,INEXT1,INXTP1)
      RANNO = RAN1(ISEED1,MA1,INEXT1,INXTP1)
C     The .0000000001 ensures ISEED is the same for different compilers
C     for match with Oman 98 results
C      X  = -RAN1(ISEED3,MA3,INEXT3,INXTP3)*100000.D0 + 0.0000000001d0
      X  = -RAN1(ISEED3,MA3,INEXT3,INXTP3)*100000.D0
      ISEED  = INT(X)


C     Set abundance estimate SIGHT & CVX if a survey is done this year.
      IF (IYR.EQ.ISIGHT) THEN

        IF (OPTDET.EQ.0) THEN
C         Stochastic trial: Set SIGHT = Bias * K * BETASQ * EXP(EPS) * W
C             where W is the Poisson component.  
C             1st set W: W takes its expectation value if (P/K)/BETA**2 >70
          W = PSURV(IYR)/BKCON
          IF (W .LT. 70.D0) W = POISSN (W, RANNO)

          SIGHT(IYR) = BIAS * BKCON * EXP(EPS) * W

C         Set CVX, the estimate of the CV:
          IF (SIGHT(IYR).GT.0.D0) THEN
C           CHISQ= random no from chi square distribution (DOF=deg.of freedom)
C           CVX = THTASQ (ASQ + BSQ/(W*BETASQ))
            CHISQ = GAMMA (DOF, SQRT(DOF+DOF), ISEED)
            SIG2 = LOG(1.D0 + ACON + BCON/W)
            CVX(IYR) = SQRT(SIG2 * CHISQ/DOF)
          ELSE
C           Zero estimate. Store Z(i), the Poisson multiplier, in CVX
            CVX(IYR) = BIAS * BKCON * EXP(EPS)
          ENDIF

        ELSE
C         Deterministic trial
          SIGHT(IYR) = BIAS * PSURV(IYR)
          CVX(IYR) = SQRT( LOG(1.D0 + ACON + BCON*BKCON/PSURV(IYR)))

        ENDIF

C       Set year of next survey
        IF (IFREQ.GT.0) THEN
          ISIGHT = ISIGHT + IFREQ
          IF (ISIGHT.GT.ENDSUR) ISIGHT = NYEAR + 1
        ELSE
C         Time to next survey is taken randomly from U[1,9]
C         If ISIGHT>ENDSUR then no more surveys
          ISIGHT = ISIGHT + INT(RAN1(ISEED5,MA5,INEXT5,INXTP5)*9.D0) + 1
        ENDIF

C       Count number of zero estimates
        IF (SIGHT(IYR).LE.0.D0) NZERO = NZERO + 1

      ELSE
C       No survey this year.  Set SIGHT = -1
        SIGHT(IYR)  = -1.D0
        CVX(IYR)  = -1.D0
      ENDIF

      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE PMCAT (CATCHM,IPPROT,INITYR,CATCH,OPTC,CATERR)

C     PMCAT sets reported premanagement catch (CATCHM)

      DOUBLE PRECISION CATCHM(-65:2000),CATCH(-65:2000),CATERR
      INTEGER IPPROT,INITYR,OPTC,IYR

      DO 5 IYR = -65,-1
        CATCHM(IYR) = 0.D0
    5 CONTINUE

      IF (OPTC.EQ.0) THEN
C       Reported historic catch = true catch = 1 unless CATCH was reset by
C       STKUPD or STKUPA (this can happen if the required initial depletion
C       is small and there is a period of protection prior to management). 
        DO 10 IYR=INITYR,IPPROT-1
          CATCHM(IYR) = CATCH(IYR)
   10   CONTINUE

      ELSE
C       Management procedure is passed erroneous catch history
C       Reported catch = True catch * CATERR
        DO 60 IYR=INITYR,IPPROT-1
          CATCHM(IYR) = CATCH(IYR) * CATERR      
   60   CONTINUE
        WRITE (1,999) 'Reported',CATCHM(INITYR),CATCHM(IPPROT-1)
      ENDIF

      WRITE (1,999) 'True',CATCHM(INITYR),CATCHM(IPPROT-1)
  999 FORMAT (1X,A,' historic catch',T35,2F7.2)
      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE PNTOUT (IOUT,M,ISTEP,NPCAT,NPPROT,DEPL,NTRIAL,
     +              NYEAR,REF,MORT1,MORT2,OPTRAN,OPTB,IRAN,DESC,PROG,
     +              OPTDEPL,INITYR,ITUNE,DEPLREF,DEPLCOL)

C     Print parameter list & check values are within range.

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC,OPTDEPL

      COMMON /PTVRS/ K,A,SREC,S,TM,TLAG, K99,MSYR99,KYEAR,MSYRYR
      DOUBLE PRECISION K(-65:2000),A(-65:2000),SREC(-65:2000),S, K99,
     +  MSYR99
      INTEGER TM, TLAG, KYEAR, MSYRYR

      COMMON /AGEPAR/ MAT1,MSIG,REC1,RSIG,FEC,SUR(0:60),UNREC(0:60),
     +                REC(0:60),RECF(0:60),FMATUR(0:60),MAXAGE,MINMAT
      DOUBLE PRECISION MAT1,MSIG,REC1,RSIG,FEC,SUR,UNREC,REC,RECF,FMATUR
      INTEGER MAXAGE,MINMAT

      COMMON /DATVRS/ BIAS,CV1EST,CV2EST,DOFMIN,ETA,IFREQ,ISIGHT,
     +                ENDSUR,OPTDET,OPTSUR,IYRCV
      DOUBLE PRECISION BIAS,CV1EST,CV2EST,DOFMIN,ETA
      INTEGER IFREQ,ISIGHT,ENDSUR,OPTDET,OPTSUR,IYRCV

      COMMON /MANPAR/ PPROB, PYMAX, PNYSTP, PKSTEP, PDSTEP,PNBSTP,PBMIN,
     1                PBMAX, PSCALE,PHASET, PHASEP, PCYCLE,PLEVEL,PSLOPE
      DOUBLE PRECISION PPROB, PYMAX, PNYSTP, PKSTEP, PDSTEP,PNBSTP,
     +        PBMIN,PBMAX, PSCALE,PHASET, PHASEP, PCYCLE,PLEVEL,PSLOPE
      COMMON /CLCPAR/ IASESS, LSURV, PASESS, POUT
      INTEGER IASESS, LSURV, PASESS, POUT

      COMMON /HITFIT/ OPTDD,OPTF,OPTMSYL,FMSY,SMSYL,RFTOT0,RFEXP0,
     +       FMSYL,FECMSY,MSYLT1,MSYLT0,MSYLT2,AMSYR1,AMSYL1,
     +       AMSYR0,AMSYL0,AMSYR2,AMSYL2
      DOUBLE PRECISION FMSY,SMSYL,RFTOT0,RFEXP0,FMSYL,FECMSY
      DOUBLE PRECISION AMSYR0,AMSYL0,AMSYR1,AMSYL1,AMSYR2,AMSYL2
      DOUBLE PRECISION MSYLT0,MSYLT1,MSYLT2
      INTEGER OPTDD,OPTF,OPTMSYL
      
      DOUBLE PRECISION M, DEPL, MORT1,MORT2
      INTEGER IOUT,ISTEP,NPCAT,NPPROT,NTRIAL,NYEAR,OPTRAN,OPTB,I,IRAN
      INTEGER INITYR,DEPLCOL
      CHARACTER OPT(-1:5)*40, REF*10, DESC*50, PROG*8,DEPLREF*9
      LOGICAL ERR
      ERR = .FALSE.
      OPT(-1) = '**** ERROR: UNKNOWN VALUE'

C     Print out parameter list
      WRITE (IOUT,'(/2A,2A15/A10,2x,I3/)') 'CASE: ',REF,'RESULTS',PROG,
     +   DEPLREF,DEPLCOL
      OPT(1) = '1 true stock managed as 1'
      WRITE (IOUT,95) 'Option controlling numbers of stocks',0,OPT(1)

      OPT(0) = 'Constant values'
      OPT(1) = 'Random parameters'
      OPT(2) = 'Random number generators not reseeded'
      OPT(3) = 'Random number seeds are generated'
      I = OPTRAN
      IF (OPTRAN.LT.0.OR.OPTRAN.GT.4) THEN
        ERR = .TRUE.
        I = -1
      ENDIF
      WRITE (IOUT,95) 'Random parameters option',OPTRAN,OPT(I)

      OPT(0) = 'Single value'
      OPT(1) = 'Pre-specified set'
      I = OPTDEPL
      IF (OPTRAN.LT.0.OR.OPTRAN.GT.1) THEN
        ERR = .TRUE.
        I = -1
      ENDIF
      WRITE (IOUT,95) 'Depletion values',OPTDEPL,OPT(I)

      IF (K99.EQ.-1.) THEN
        WRITE (IOUT,97) 'K cyclic, starting minimum',K99
      ELSE IF (K99.EQ.-2.) THEN
        WRITE (IOUT,97) 'K cyclic, starting maximum',K99
      ELSE IF (K99.GT.0.) THEN 
        WRITE (IOUT,97) 'K varies linearly to K1*',K99,KYEAR
      ELSE
        WRITE (IOUT,97) 'K is constant',K99
      ENDIF

      OPT(0) = 'Constant bias'
      OPT(1) = 'Bias doubles linearly'
      OPT(2) = 'Bias decreases to 2/3 initial value'
      I = OPTB
      IF (OPTB.LT.0.OR.OPTB.GT.2) THEN
        ERR = .TRUE.
        I = -1
      ENDIF
      WRITE (IOUT,95) 'Variable bias option',OPTB,OPT(I)

      IF (MSYR99.GT.0.D0) THEN
        WRITE (IOUT,97) 'MSYR changes linearly to MSYR1*',MSYR99,MSYRYR
      ELSE
        WRITE (IOUT,97) 'MSYR99 = ',MSYR99
      ENDIF

      OPT(0) = 'True reported catch'
      OPT(1) = 'Error in reported catch; Depletion fixed'
      OPT(2) = 'Error in reported catch; PTRUE(0) fixed'
      I = OPTC
      IF (OPTC.LT.0.OR.OPTC.GT.2) THEN
        ERR = .TRUE.
        I = -1
      ENDIF
      WRITE (IOUT,95) 'Option controlling catch history',OPTC,OPT(I)

      OPT(0) = 'Standard P-T model'
      IF (OPTMOD.EQ.0) THEN
        IF (K99.NE.0. .OR. MSYR99.NE.0. .OR. ISTEP.NE.0)
     +    OPT(0) = 'P-T model with variable K &/or MSYR'
      END IF
      OPT(1) = 'P-T + Limit max recruitment'
      OPT(2) = 'Tent model'
      OPT(3) = 'Age-structured, Tm=Tr'
      OPT(4) = 'Age-structured, Tm.NE.Tr'
      OPT(5) = 'Age-structured, comps spec'
      I = OPTMOD
      IF (OPTMOD.LT.0.OR.OPTMOD.GT.5) THEN
        ERR = .TRUE.
        I = -1
      ENDIF
      WRITE (IOUT,95) 'Option defining population model',OPTMOD,OPT(I)

      OPT(0) = 'Stochastic run'
      OPT(1) = 'Deterministic run'
      I = OPTDET
      IF (OPTDET.LT.0.OR.OPTDET.GT.1) THEN
        ERR = .TRUE.
        I = -1
      ENDIF
      WRITE (IOUT,95) 'Option controlling stochasticity',OPTDET,OPT(I)

      OPT(0) = 'No'
      OPT(1) = 'Yes'
      I = OPTSUR
      IF (OPTSUR.LT.0.OR.OPTSUR.GT.1) THEN
        ERR = .TRUE.
        I = -1
      ENDIF
      WRITE (IOUT,95) 'Survey costs taken into account:',OPTSUR,OPT(I)
      WRITE (IOUT,'()')
      WRITE (IOUT,95) 'Number of trials',NTRIAL
      WRITE (IOUT,95) 'Number of years in simulation',NYEAR
      WRITE (IOUT,95) 'Number of years of premanagement catch',NPCAT
      WRITE (IOUT,95) 'Number of years of protection',NPPROT
      WRITE (IOUT,95) 'Number of true stocks',1
      WRITE (IOUT,95) 'Number of managed whaling grounds',1
      WRITE (IOUT,95) 'MSYR step length (in years)',ISTEP
      IF (OPTMOD.LT.3) THEN
        WRITE (IOUT,95) 'Age at recruitment (TM)',TM
        WRITE (IOUT,95) 'Time lag in recruitment',TLAG
        WRITE (IOUT,96) 'Mortality rate',M
      ELSE
        WRITE (IOUT,'(A,'' parameters'',T41,2(F9.2,9X),I3/
     +                A,'' parameters'',T41,2(F9.2,9X)/,
     +                A,'' parameters'',T41,F9.2,F6.2)') 
     +               ' Maturity',MAT1,MSIG,MINMAT,
     +               ' Recruitment',REC1,RSIG,
     +               ' Mortality',MORT1,MORT2
        WRITE (IOUT,95) 'Maximum age class',MAXAGE
C
C       Output Components
        OPT(0) = 'EXPLOITABLE'
        OPT(1) = 'TOTAL1+'
        OPT(2) = 'MATURE'
        I = OPTF
        IF (OPTF.LT.0.OR.OPTF.GT.2) THEN
          ERR = .TRUE.
          I = -1
        ENDIF
        WRITE (IOUT,95) 'MSYR component:',OPTF,OPT(I)
        I = OPTMSYL
        IF (OPTMSYL.LT.0.OR.OPTMSYL.GT.2) THEN
          ERR = .TRUE.
          I = -1
        ENDIF
        WRITE (IOUT,95) 'MSYL component:',OPTMSYL,OPT(I)
        I = OPTDD
        IF (OPTDD.LT.0.OR.OPTDD.GT.2) THEN
          ERR = .TRUE.
          I = -1
        ENDIF
        WRITE (IOUT,95) 'Density-dependent component:',OPTDD,OPT(I)
C      
      ENDIF
      WRITE (IOUT,'(A,(T16,2(A8,F10.7),3(A9,F6.1)))')' CLA parameters',
     1   'PPROB:',PPROB, 'PKSTEP:',PKSTEP,'PNYSTP:',PNYSTP,
     2   'PHASET:',PHASET,'PASESS:',REAL(PASESS), 'PLEVEL:',PLEVEL,
     3   'PDSTEP:',PDSTEP,'PNBSTP:',PNBSTP,'PCYCLE:',PCYCLE
      WRITE (1,'(A,(T20,4(A10,F10.7)))')
     + ' CLA parameters',
     1 'PPROB',PPROB, 'PNYSTP',PNYSTP, 'PKSTEP',PKSTEP, 'PDSTEP',PDSTEP, 
     2 'PNBSTP',PNBSTP,'PCYCLE',PCYCLE, 'PLEVEL',PLEVEL, 'PASESS',
     3  REAL(PASESS)
      PRINT '(//3(A,F10.7,3X),2(A,F7.1,3X)//)', ' PROB',PPROB,
     +      'KSTP',PKSTEP,'DSTP',PDSTEP,'YSTP',PNYSTP,'BSTP',PNBSTP
      IF (ITUNE.GT.0) THEN
       write (31,'(3(F10.7,A,3X,I6))') PPROB,' =PPROB;',IRAN
      ENDIF 
      WRITE (IOUT,96) 'MSYL', MSYL
      WRITE (IOUT,96) 'Initial MSY rate', MSYR1
      WRITE (IOUT,96) 'Density dependent exponent (Z)',Z
      WRITE (IOUT,96) 'Resilience parameter (A)',A1
      WRITE (IOUT,96) 'Carrying capacity (mature)', PTRUE(INITYR)
      WRITE (IOUT,96) 'Carrying capacity (1+)',K1P
      WRITE (IOUT,96) 'Initial depletion (P0/K)', DEPL
      IF (OPTC.GE.1) THEN
        WRITE (IOUT,96) 'Error in reported catch (CATERR)',CATERR 
      ELSE
        WRITE (IOUT,'()')
      ENDIF
      WRITE (IOUT,96) 'Bias in abs. abundance estimates (initial)',BIAS
      IF (IYRCV.GE.NYEAR) THEN
        WRITE (IOUT,98) 'CV1EST',CV1EST,'DOFMIN',DOFMIN,'ETA',ETA
      ELSE
        WRITE (IOUT,98) 'CV1EST',CV1EST,'DOFMIN',DOFMIN,'ETA',ETA,
     +                  'CV2EST',CV2EST,'IYRCV',IYRCV
      ENDIF
      WRITE (IOUT,96) 'Frequency of epidemics',ERATE
      WRITE (IOUT,95) 'Absolute abundance estimate frequency',IFREQ
      WRITE (IOUT,95) 'No surveys after year',ENDSUR

      IF (ERR) STOP
 
   95 FORMAT (' ',A,T41,I6,4X,A)
   96 FORMAT (' ',A,T41,F12.5,F14.5)
   97 FORMAT (' ',A,T41,F6.2, : 4X,'starting in year',I3) 
   98 FORMAT (' Survey parameters:',T20,4(A8,':',F6.3),A9,I5)

      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE REPORT (IOUT, CATCH, IYR)

C     REPORT prints out PTRUE & CATCH arrays from year 0

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      DOUBLE PRECISION CATCH(-65:2000)
      INTEGER IOUT,IYR,I

      WRITE (IOUT,'(A,A9,6A11)') 'Year','Stock 1','CM 1'
      DO 10 I = 0,IYR-1
        WRITE (IOUT,'(I3,6F11.5)') I,PTRUE(I), PSURV(I), CATCH(I)
   10 CONTINUE
      WRITE (IOUT,'(I3,6F11.5)') IYR, PTRUE(IYR),PSURV(IYR)
      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE RESET (ISEEDS,N,NYEAR,BIAS0,BIAS,OPTRAN,IRAN)

C     RESET reseeds the random number generators, resets BIAS to its 
C           initial value and sets PROBE

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      COMMON /RANPAR/ ISEED1,ISEED2,ISEED3,ISEED4,ISEED5,SSEED1,SSEED3
      INTEGER ISEED1, ISEED2,ISEED3,ISEED4,ISEED5,SSEED1,SSEED3

      EXTERNAL RAN1
      DOUBLE PRECISION BIAS0,BIAS,RAN1, MA1(55), MA2(55),MA3(55)
      INTEGER ISEEDS(400,5),N,NYEAR,OPTRAN,IRAN,IYR,ISEED
      INTEGER INEXT1,INXTP1,INEXT2,INXTP2,INEXT3,INXTP3
      SAVE MA1,MA2,MA3,INEXT1,INXTP1,INEXT2,INXTP2,INEXT3,INXTP3

C     Set new negative values for the random number generator seeds -
C     the generator is reset whenever a negative seed is used.
C     If OPTRAN=2 the generators are NOT reseeded and IRAN defines the 
C       inital seeds.  It will not be reproducible on different machines but
C       was added as a method of getting a different set of random numbers 
C       for each run in the factorial expt. - in the event it was not used. 
C     If OPTRAN=3 the generator is not reseeded. IRAN defines the inital seeds.
C       Used in tuning when 100000 replicates are needed
C     If OPTRAN=4 the a new seed is generated (1 & 3 only) instead of
C       using values which were read in. IRAN defines the inital seeds.
      ISEED = N
      IF (OPTRAN.LT.2 .OR. N.EQ.1) THEN
        IF (OPTRAN.GE.2) ISEED = IRAN
        ISEED1 = ISEEDS(ISEED,1)
        ISEED2 = ISEEDS(ISEED,2)
        ISEED3 = ISEEDS(ISEED,3)
        ISEED4 = ISEEDS(ISEED,4)
        ISEED5 = ISEEDS(ISEED,5)
      ELSE IF (OPTRAN.EQ.4) THEN
        ISEED1 = INT(-RAN1(SSEED1,MA1,INEXT1,INXTP1)*1000000.d0)
        ISEED3 = INT(-RAN1(SSEED3,MA3,INEXT3,INXTP3)*1000000.d0)
      ENDIF
C     Store the seeds used in the tuning trials
      SSEED1 = ISEED1
      SSEED3 = ISEED3

C     Reset the bias to its initial value
      BIAS = BIAS0

C     Set PROBE: Generate years in which epidemics occur. 1st initialise to 0.
      DO 6 IYR=0,NYEAR
        PROBE(IYR) = 0
    6 CONTINUE
      IF (ERATE.GT.0.D0) THEN
        DO 10 IYR=1,NYEAR
          IF (RAN1(ISEED2,MA2,INEXT2,INXTP2).LE.ERATE) PROBE(IYR) = 1
   10   CONTINUE
      ENDIF

      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE SETA (INITYR,NYEAR,ISTEP)

C     P-T model only.  SETA sets array A (the resilience) for each year 

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      COMMON /PTVRS/ K,A,SREC,S,TM,TLAG, K99,MSYR99,KYEAR,MSYRYR
      DOUBLE PRECISION K(-65:2000),A(-65:2000),SREC(-65:2000),S, K99,
     + MSYR99
      INTEGER TM, TLAG, KYEAR, MSYRYR

      DOUBLE PRECISION MSYR2,A2,MSYINC
      INTEGER INITYR,NYEAR,ISTEP,IYR

C     MSYR1 = initial MSY rate.  A1 = resilience 
      A1 = MSYR1 * S * (1.D0 + 1.D0/Z) / (1.D0-S)

      IF (MSYR99.GT.0.D0) THEN
C       MSYR varies linearly from MSYR1 to MSYR1*MSYR99 over the period 
C                                 MSYRYR to NYEAR-1
        MSYINC = MSYR1 * (MSYR99 - 1.D0) / DBLE(NYEAR-1-MSYRYR)
        DO 5 IYR=INITYR,MSYRYR
          A(IYR) = A1
    5   CONTINUE
        MSYR2 = MSYR1
        DO 6 IYR=MSYRYR+1,NYEAR
          MSYR2 = MSYR2 + MSYINC
          A(IYR) = MSYR2 * S * (1.D0 + 1.D0/Z) / (1.D0-S)
    6   CONTINUE
  
      ELSE IF (ISTEP.EQ.0) THEN
C       MSYR remains constant  
        DO 10 IYR=INITYR,NYEAR
          A(IYR) = A1
   10   CONTINUE

      ELSE
C       MSYR changes every ISTEP years from MSYR1 to .05-MSYR1 & back
        MSYR2 = .05 - MSYR1
        A2 = MSYR2 * S * (1.D0 + 1.D0/Z) / (1.D0-S)

        DO 20 IYR=INITYR,-1
          A(IYR) = A1
   20   CONTINUE
   25   DO 30 IYR = IYR, MIN (IYR+ISTEP-1, NYEAR)
          A(IYR) = A2
   30   CONTINUE
        DO 40 IYR = IYR, MIN (IYR+ISTEP-1, NYEAR)
          A(IYR) = A1
   40   CONTINUE
        IF (IYR.LE.NYEAR) GO TO 25
      ENDIF
      IF (MSYR99.GT.0.D0 .OR. ISTEP.NE.0) WRITE (1,'(2(A,6F8.2,9X))')
     +             ' A in yrs 0,1,2,25,50,99:',A(0),A(1),A(2),
     +             A(25),A(50),A(99),'MSYR = A / ', S*(1.+1./Z)/(1.-S)
      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE SETK (D,CATCH,IPPROT,INITYR,NYEAR)

C     P-T model only.
C     SETK sets K, CATCH(IYR) & PTRUE(IYR) for IYR=INITYR+1,0
C     Subroutine finds the carrying capacity K giving the required 
C     depletion D in year 0.  A catch of 1 unit is taken for NPCAT years
C     followed by a protection period of NPPROT years (if any).
C     Note: In cases with a protection period the catch may be reduced 
C           if PTRUE becomes too small to support it (STKUPD resets the 
C           catch in this case). CATCH stores the actual catch taken.
C     PTRUE is also set (by STKUPD) for years INITYR+1 to 0.
 
C     A catch of 1 unit is set so that management procedures under test
C     can not derive information they are not entitled to from the size of
C     the catch. (This was possible, in an early version of the program,
C     when K was fixed & the catch giving the required depletion found)

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      COMMON /PTVRS/ K,A,SREC,S,TM,TLAG, K99,MSYR99,KYEAR,MSYRYR
      DOUBLE PRECISION K(-65:2000),A(-65:2000),SREC(-65:2000),S, K99,
     +  MSYR99
      INTEGER TM, TLAG, KYEAR, MSYRYR

      DOUBLE PRECISION D,CATCH(-65:2000),KMAX,KMIN,KINC,C,CHIST,SS,
     +       PIFACT,P0
      INTEGER IPPROT,INITYR,NYEAR,ICOUNT,IYR

C     Find the mature carrying capacity K1 giving a depletion D after removing
C     a catch CHIST for NPCAT yrs. The historic catch CHIST=1 unless OPTC=1
      KMIN = 0.D0
      KMAX = 10000.D0
      ICOUNT = 0
      DO 25 IYR=IPPROT,-1
        CATCH(IYR) = 0.D0
   25 CONTINUE
      CHIST = 1.D0
      IF (OPTC.EQ.1) CHIST = CHIST / CATERR

C --- Start of loop.  Repeat until projected population is at required level
  100 K1 = (KMIN + KMAX) / 2.D0
      PTRUE(INITYR) = K1
      DO 105 IYR=INITYR,0
        K(IYR) = K1
  105 CONTINUE

C     Forward project the population as necessary
C     (CATCH(IYR) is reset each time as it may be reduced in STKUPD)
      DO 110 IYR = INITYR,IPPROT-1
        CATCH(IYR) = CHIST
        CALL STKUPD (IYR,INITYR,CATCH(IYR))
  110 CONTINUE
C     Protection period (CATCH(IYR) was set to 0 above)
      C = 0.D0
      DO 120 IYR = IPPROT,-1
        CALL STKUPD (IYR,INITYR,C)
  120 CONTINUE

C     Has convergence taken place?
      IF ((ABS(PTRUE(0)/K1 - D).GE.0.0001D0) .or. (ICOUNT.LT.500)) THEN
        IF ((PTRUE(0)/K1 - D).GT.0.D0) THEN
          KMAX = K1
        ELSE
          KMIN = K1
        ENDIF

C       Update loop count and check is not wasting time
        ICOUNT = ICOUNT + 1
C        IF (ICOUNT.GT.500) 
        
c STOP '**** ERROR: SETK HAS NOT CONVERGED'

      GO TO 100

C --- ELSE: Routine has converged so continue
      ENDIF
  

C     If OPTC=2 the reported catch, is erroneous.  The current (year 0) 
C       population size, PTRUE(0) is fixed at the same size as when 
C       the reported = true catch (=1). This was set above.
 101  IF (OPTC.EQ.2) THEN
C       Now find K (& hence also D) which gives a population PTRUE(0), 
C       in year 0, when the true catch = reported catch / CATERR.
        P0 = PTRUE(0)
        KMIN = P0
        KMAX = 10000.D0
        CHIST = 1.D0 / CATERR
        ICOUNT = 0

C --- Start of loop.  Repeat until projected population is at required level
  600   K1 = (KMIN + KMAX) / 2.D0
        PTRUE(INITYR) = K1
        DO 605 IYR=INITYR,0
          K(IYR) = K1
  605   CONTINUE

C       Forward project the population as necessary
        DO 610 IYR = INITYR,IPPROT-1
          CATCH(IYR) = CHIST
          CALL STKUPD (IYR,INITYR,CATCH(IYR))
  610   CONTINUE
        C = 0.D0
        DO 620 IYR = IPPROT,-1
          CALL STKUPD (IYR,INITYR,C)
  620   CONTINUE

C       Has convergence taken place?
        IF (ABS(PTRUE(0)/P0-1.D0) .GE.0.0001D0) THEN
          IF ((PTRUE(0) - P0) .GT.0.D0) THEN
            KMAX = K1
          ELSE
            KMIN = K1
          ENDIF
C         Update loop count and check is not wasting time
          ICOUNT = ICOUNT + 1
          IF (ICOUNT.GT.500) STOP '**** ERROR: SETK HAS NOT CONVERGED'

          GO TO 600

        ELSE
C ---     Routine has converged
          D = PTRUE(0)/K1
        ENDIF
      ENDIF

C     Set KSURV = pristine 1+ population
C               = mature population (K1) + whales of age 1 to TM
      KSURV = K1
      SS = 1
      DO 830 IYR = 1,TM
        SS = SS * S
        KSURV = KSURV + SREC(INITYR)/SS
  830 CONTINUE
C
C     AEP check this AEP!!!!
      K1P = KSURV

C     Store PSURV = survey population in year -1
      PSURV1 = PSURV(-1)


C     Set K array = Mature carrying capacity in IYR. K is constant if K99=0
      IF (K99 .EQ. 0.D0) THEN
        DO 850 IYR=1,NYEAR
          K(IYR) = K1
  850   CONTINUE
 
      ELSE IF (K99 .GT. 0.D0) THEN
C      K is constant up to year KYEAR, and then varies linearly 
C               from K1 in year KYEAR to K1*K99 in year NYEAR-1
        KINC = K1 * (K99 - 1.D0) / DBLE(NYEAR-1-KYEAR)
        DO 859 IYR=1,KYEAR
          K(IYR) = K1
  859   CONTINUE
        DO 860 IYR=KYEAR+1,NYEAR
          K(IYR) = K(IYR-1) + KINC
  860   CONTINUE

      ELSE
C       K cyclic.  K99=-1: minima in years 0 & 100; max = 3K0
C                  K99=-2: maxima in years 0 & 100; min = K0/3
        PIFACT = 2.D0 * 3.141593D0 / DBLE(NYEAR)
        IF (K99.EQ.-1.D0) THEN
          DO 870 IYR = 1,NYEAR
            K(IYR) = K1 * (2.D0 - COS(PIFACT * IYR)) 
  870     CONTINUE
        ELSE IF (K99.EQ.-2.D0) THEN
          DO 880 IYR = 1,NYEAR
            K(IYR) = K1/3.D0 * (2.D0 + COS(PIFACT * IYR))
  880     CONTINUE
        ENDIF      

      ENDIF

      WRITE (1,'(A,6F8.2)') ' K in yrs 0,1,2,25,50,99:',K(0),K(1),K(2),
     +                  K(25),K(50),K(99)
      
      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE SETPAR (MSYR1,DEPL,NPCAT)

C     Subroutine SETPAR sets MSYR1 to a random value from U[.001,.05]
C                            DEPL  to a random value from U[.01,.99]
C                            NPCAT to a random integer from U[15,40]

      COMMON /RANPAR/ ISEED1,ISEED2,ISEED3,ISEED4,ISEED5,SSEED1,SSEED3
      INTEGER ISEED1, ISEED2, ISEED3, ISEED4, ISEED5,SSEED1,SSEED3

      DOUBLE PRECISION MSYR1, DEPL, RAN1, MA4(55)
      INTEGER NPCAT, INEXT4,INXTP4
      EXTERNAL RAN1
      SAVE MA4,INEXT4,INXTP4

      MSYR1   = RAN1(ISEED4,MA4,INEXT4,INXTP4) * .049D0 + 0.001D0
      DEPL    = RAN1(ISEED4,MA4,INEXT4,INXTP4) * 0.98D0 + 0.01D0
      NPCAT = INT (RAN1(ISEED4,MA4,INEXT4,INXTP4) * 26.D0) + 15D0

      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      FUNCTION SETZ (MSYL)

C     Use bisection method to find Z.
C     NB maximum MSYL allowed = .95, corresponding to Z=87.4

      INTEGER ICOUNT
      DOUBLE PRECISION SETZ,MSYL,ZMIN,ZMAX,DERIV

      ICOUNT = 0
      ZMIN = 0.D0
      ZMAX = 90.D0

   10 SETZ = (ZMIN + ZMAX) / 2.D0
      DERIV = 1.D0 - (SETZ+1.D0)*(MSYL**SETZ)
      IF (ABS(DERIV).LT.0.00001D0) GO TO 12
      IF (DERIV.LT.0.D0) THEN
        ZMIN = SETZ
      ELSE
        ZMAX = SETZ
      ENDIF

      ICOUNT = ICOUNT + 1
      IF (ICOUNT.GT.500) STOP '**** ERROR: SETZ HAS NOT CONVERGED'
      GO TO 10

   12 CONTINUE

      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE STKUPD (IYR,INITYR,C)

C     P-T model only.
C     STKUPD updates the stock size at start of IYR+1 
C     ie sets PTRUE(IYR+1), removing C=catch in IYR.  PSURV is also set.

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      COMMON /PTVRS/ K,A,SREC,S,TM,TLAG, K99,MSYR99,KYEAR,MSYRYR
      DOUBLE PRECISION K(-65:2000),A(-65:2000),SREC(-65:2000),S, K99,
     + MSYR99
      INTEGER TM, TLAG, KYEAR, MSYRYR

      DOUBLE PRECISION C,REC,S1,CONST,SS
      INTEGER IYR,INITYR, ITM,ITLAG,I,I2

C     Set constants:
      S1 = 1.D0 - S

C     Store this year's births which survive to recruitment (in SREC)
C     They will enter the recruited population in IYR+TM. 
C     NB: for standard PT model TM=TLAG.  
      IF (OPTMOD.LE.1) THEN
        ITLAG = IYR-TLAG+TM
        IF (ITLAG.LT.INITYR) ITLAG = INITYR

        REC =S1*PTRUE(IYR)*(1.D0+A(IYR)*(1.-(PTRUE(ITLAG)/K(ITLAG))**Z))
        IF (OPTMOD.EQ.1) REC = MIN (REC, 
     +                     S1*PTRUE(IYR) * (1.D0+A(IYR)*(1.D0-MSYL**Z)))

      ELSE
C       Tent model:  (constant MSYR only)
        IF (PTRUE(IYR) .LE. K1*MSYL) THEN
          REC = PTRUE(IYR) * (S*MSYR1 + S1)
        ELSE
          CONST = S * MSYL * MSYR1 / (1.-MSYL)
          REC = PTRUE(IYR) * (-CONST + S1) + CONST*K1
        ENDIF
      ENDIF

      SREC(IYR) = MAX (REC,0.D0)

C     Add this years recruits (born in IYR-TM) to population & remove catches
      ITM = IYR-TM
      IF (ITM.LT.INITYR) ITM = INITYR

      IF (PTRUE(IYR).GE.C) THEN
        PTRUE(IYR+1) = S * (PTRUE(IYR) - C) + SREC(ITM)
      ELSE
        C = PTRUE(IYR)
        PTRUE(IYR+1) = SREC(ITM)
      ENDIF

      IF (IYR.GE.-1) THEN
C       Set PSURV = PTRUE + number unrecruited age 1 & over
C       Note: SREC(IYR) = # born in IYR which will survive to age TM+1
C       Another AEP change AEP !!!
        PSURV(IYR) = PTRUE(IYR)
        PSURV(IYR+1) = PTRUE(IYR+1)
        SS = 1
        DO 20 I2 = IYR-TM,IYR-1
          I = MAX(I2,INITYR)
          SS = SS * S
          PSURV(IYR) = PSURV(IYR) + SREC(I) / SS
          PSURV(IYR+1) = PSURV(IYR+1) + SREC(I+1) / SS
   20   CONTINUE

C       If an epidemic occurs the population is halved - NOT IN USE  C
C       BEFORE USE ADJUST THE MORTALITY RATE - FOR STANDARD P-T ONLY C
C        IF (PROBE(IYR+1).EQ.1) THEN                                 C
C         PTRUE(IYR+1) = 0.5*PTRUE(IYR+1)                            C
C         DO 15 I=IYR-TM+1,IYR                                       C
C           SREC(I) = SREC(I) * 0.5D0                                C
C  15     CONTINUE                                                   C
C       ENDIF                                                        C
      ENDIF

      RETURN
      END


C *********************************************************************
C *********************************************************************

      SUBROUTINE SETKA (D,CATCH,IPPROT,INITYR,NK,MORT1,MORT2,MORTIP,
     +    NYEAR)

C  Age structured model only.
C  Find the mature carrying capacity K giving the required depletion, D
C  and set up the premanagement population.  
C  K is already set if NK=1 (i.e. K was found on the 1st call to this routine)
C  NB the mature & exploitable populations are both set in this routine,
C   but are identical if their ogives are the same, as in the base case.

      IMPLICIT NONE

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      COMMON /AGEPAR/ MAT1,MSIG,REC1,RSIG,FEC,SUR(0:60),UNREC(0:60),
     +                REC(0:60),RECF(0:60),FMATUR(0:60),MAXAGE,MINMAT
      DOUBLE PRECISION MAT1,MSIG,REC1,RSIG,FEC,SUR,UNREC,REC,RECF,FMATUR
      INTEGER MAXAGE,MINMAT
      
      COMMON / TIMEVAR/ AT(-65:2000),KT(-65:2000),ZT(-65:2000),
     +                  K99AEP,KYEARAEP,MSYR99AEP,MSYRYRAEP,ISTEPAEP
      DOUBLE PRECISION AT,KT,ZT,K99AEP,MSYR99AEP
      INTEGER KYEARAEP,ISTEPAEP,MSYRYRAEP

      COMMON /HITFIT/ OPTDD,OPTF,OPTMSYL,FMSY,SMSYL,RFTOT0,RFEXP0,
     +       FMSYL,FECMSY,MSYLT1,MSYLT0,MSYLT2,AMSYR1,AMSYL1,
     +       AMSYR0,AMSYL0,AMSYR2,AMSYL2
      DOUBLE PRECISION FMSY,SMSYL,RFTOT0,RFEXP0,FMSYL,FECMSY
      DOUBLE PRECISION AMSYR0,AMSYL0,AMSYR1,AMSYL1,AMSYR2,AMSYL2
      DOUBLE PRECISION MSYLT0,MSYLT1,MSYLT2
      INTEGER OPTDD,OPTF,OPTMSYL
      
      DOUBLE PRECISION D, CATCH(-65:2000),KMAX,KMIN,C,CHIST,SURV,KREC,
     +       MORT1,MORT2,RINIT(0:60),UINIT(0:60),PA,RREC,RMAT,R1PLUS,P0,
     +       KINC,PIFACT,MSYR2,A1OLD,ZOLD,MSYINC,Kused,KMAT
      INTEGER INITYR,NK,IPPROT,L,ICOUNT,IYR,MORTIP,IOK,NYEAR
      EXTERNAL SURV
      SAVE P0
      
C     MSYPAR sets A1 & Z (initial MSYR for time-varying A
      CALL MSYPAR (IOK,MORT1,MORT2,MORTIP,MSYR1)
C
C     Set A and Z (defaults)
      DO 401 IYR = INITYR,NYEAR
       AT(IYR) = A1
       ZT(IYR) = Z
401   CONTINUE        
C
      IF (MSYR99AEP.EQ.1.0d0) THEN
C      
      ELSE IF  (MSYR99AEP.GT.0d0) THEN
C      MSYR varies linearly from MSYR1 to MSYR1*MSYR99AEP over the period
C                                       MSYRYR to NYEAR-1
       MSYINC = MSYR1 * (MSYR99AEP-1.0d0)/DBLE(NYEAR-1-MSYRYRAEP)
       DO 415 IYR = INITYR,MSYRYRAEP
        AT(IYR) = A1
        ZT(IYR) = Z
415    CONTINUE
       MSYR2 = MSYR1
       DO 416 IYR = MSYRYRAEP+1,NYEAR
        MSYR2 = MSYR2 + MSYINC
        CALL MSYPAR (IOK,MORT1,MORT2,MORTIP,MSYR2)
        AT(IYR) = A1
        ZT(IYR) = Z
416    CONTINUE        
C
      ELSE IF (ISTEPAEP.EQ.0) THEN
C      MSYR remains constant
       DO 425 IYR = INITYR,NYEAR
        AT(IYR) = A1
        ZT(IYR) = Z
425    CONTINUE        
C
      ELSE
C      MSYR changes every ISTEP years from MSYR1 to .05-MSYR1 & back
       MSYR2 = .05 - MSYR1
       A1OLD = A1
       ZOLD = Z
       CALL MSYPAR (IOK,MORT1,MORT2,MORTIP,MSYR2)

       DO 431 IYR = INITYR,-1
        AT(IYR) = A1OLD
        ZT(IYR) = ZOLD
431    CONTINUE
499    DO 432 IYR = IYR, MIN(IYR+ISTEPAEP-1,NYEAR)
        AT(IYR) = A1
        ZT(IYR) = Z
432    CONTINUE
       DO 433 IYR = IYR, MIN(IYR+ISTEPAEP-1,NYEAR)        
        AT(IYR) = A1OLD
        ZT(IYR) = ZOLD
433    CONTINUE
       IF (IYR.LE.NYEAR) GOTO 499
C               
      ENDIF
      CALL MSYPAR (IOK,MORT1,MORT2,MORTIP,MSYR1)
      
C     Set up the recruitment ogive
      CALL SETO (RECF,RSIG,REC1,MAXAGE)

C     Reset SUR if necessary i.e. in trials where epidemics occur after year 0
      IF (ERATE.GT.0.D0) THEN
        DO 4 L = 0,MAXAGE
          SUR(L) = SURV(MORT1,MORT2,L,MORTIP)
    4   CONTINUE
      ENDIF

C     Calculate the relative recruited population size starting with
C     unity in the zero age class (A=0)
      PA   = SUR(0)
      RMAT = FMATUR(0)
      RREC = RECF(0)
      DO 9 L = 1,MAXAGE-1
        RMAT = RMAT + PA*FMATUR(L)
        RREC = RREC + PA*RECF(L)
        PA   = PA*SUR(L)
    9 CONTINUE

C     Adjust for last age class being pooled (and always fully recruited)
      PA   = PA/(1.D0 - SUR(MAXAGE))
      RREC = RREC + PA
      RMAT = RMAT + PA

C     Scale the zero age class so the relative recruited population is unity. 
      PA  = 1.D0/RREC
C     Set up recruited & unrecruited age vectors (relative to recruited 
C     population) and RMAT = relative mature population
C     Also set relative 1+ population R1PLUS = 1 + # unrecruited -age 0 recruits
      UINIT(0) = PA*(1.D0-RECF(0))
      RINIT(0) = PA*RECF(0)
      R1PLUS = 1.D0 - RINIT(0)
      DO 20 L = 1,MAXAGE
        PA = PA*SUR(L-1)
        UINIT(L) = PA*(1.D0 - RECF(L))
        RINIT(L) = PA*RECF(L)
        R1PLUS = R1PLUS + UINIT(L)
   20 CONTINUE
      RINIT(MAXAGE) = RINIT(MAXAGE)/(1.D0 - SUR(MAXAGE))

C     Set the birth rate so as to give balance at equilibrium
      FEC = 1.D0/RMAT
      IF (OPTMOD.EQ.3) FEC = 1.0/RREC

C     Reset the recruitment ogive to transition form for use by STKUPA
      CALL TRFORM (RECF,MAXAGE)

      IF (NK.EQ.1 .AND. OPTC.EQ.2) GO TO 500
      IF (NK.EQ.0) THEN
        KMIN = 0.D0
        KMAX = 10000.D0
      ELSE 
C       Set up population using K found in last call to this subroutine
        KMIN = K1
        KMAX = K1
      ENDIF
      ICOUNT = 0
      DO 25 IYR=IPPROT,-1
        CATCH(IYR) = 0.D0
   25 CONTINUE
      CHIST = 1.D0
      IF (OPTC.EQ.1) CHIST = CHIST / CATERR

C --- Start of loop.  Repeat until projected population is at required level
  100 K1 = (KMIN + KMAX) / 2.D0
C
C     Set up the carrying capacity
      KREC = RREC * K1 / RMAT 
      KSURV  = KREC * R1PLUS
      KMAT = K1
      IF (OPTMOD.EQ.3) THEN
       Kused = KREC
       KMAT = KREC
      ELSEIF (OPTMOD.EQ.4) THEN
       Kused = K1
      ELSE
       IF (OPTDD.EQ.0) Kused = KREC
       IF (OPTDD.EQ.1) Kused = KSURV
       IF (OPTDD.EQ.2) Kused = K1
      ENDIF  
C
C     Set up the carrying capacity in each year (past)
      PTRUE(INITYR) = KMAT
      DO 305 IYR = INITYR,0
       KT(IYR) = Kused
305   CONTINUE
C
C     Set up the age-structure for this K1
      DO 105 L = 0,MAXAGE
        UNREC(L) = KREC*UINIT(L)
        REC(L)   = KREC*RINIT(L)
  105 CONTINUE

C     Forward project the population as necessary
C     (CATCH(IYR) is reset each time as it may be reduced in STKUPA)
      DO 110 IYR = INITYR,IPPROT-1
        CATCH(IYR) = CHIST
        CALL STKUPA (IYR,CATCH(IYR))
  110 CONTINUE
C     Protection period (CATCH(IYR) was set to 0 above)
      C = 0.D0
      DO 120 IYR = IPPROT,-1
        CALL STKUPA (IYR,C)
  120 CONTINUE

C     Has convergence taken place? ##### suggest changing this in future
c      IF (ABS(PTRUE(0)/K1 - D).GE.0.0000001D0) THEN        
      IF ((ABS(PTRUE(0)/K1 - D).GE.0.001D0) .or. (icount.lt.500)) THEN
        IF ((PTRUE(0)/K1 - D).GT.0.D0) THEN
          KMAX = K1
        ELSE
          KMIN = K1
        ENDIF
C       Update loop count and check is not wasting time
        ICOUNT = ICOUNT + 1
C        IF (ICOUNT.GT.500) STOP '**** ERROR: SETKA HAS NOT CONVERGED'

        GO TO 100

C --- ELSE: Routine has converged so continue
      ENDIF


C     If OPTC=2 the reported catch, is erroneous.  The current (year 0) 
C       population size, PTRUE(0) is fixed at the same size as when 
C       the reported = true catch (=1). This was set above.
  500 IF (OPTC.EQ.2) THEN
C       Now find K (& hence also D) which gives a population PTRUE(0), 
C       in year 0, when the true catch = reported catch / CATERR.
        P0 = PTRUE(0)
        IF (NK.EQ.0) THEN
          KMIN = P0
          KMAX = 10000.D0
        ELSE 
          KMIN = K1
          KMAX = K1
        ENDIF
        CHIST = 1.D0 / CATERR
        ICOUNT = 0
        DO 26 IYR=IPPROT,-1
         CATCH(IYR) = 0.D0
   26   CONTINUE

C --- Start of loop.  Repeat until projected population is at required level
  600   K1 = (KMIN + KMAX) / 2.D0
C
C       Set up the carrying capacity
        KREC = RREC * K1 / RMAT 
        KSURV  = KREC * R1PLUS
        KMAT = K1
        IF (OPTMOD.EQ.3) THEN
         Kused = KREC
         KMAT = KREC
        ELSEIF (OPTMOD.EQ.4) THEN
         Kused = K1
        ELSE
         IF (OPTDD.EQ.0) Kused = KREC
         IF (OPTDD.EQ.1) Kused = KSURV
         IF (OPTDD.EQ.2) Kused = K1
        ENDIF  
  
C       Set up the carrying capacity in each year (past)
        PTRUE(INITYR) = KMAT
        DO 306 IYR = INITYR,0
         KT(IYR) = Kused
306     CONTINUE

C       Set up the age-structure for this K1
        PTRUE(INITYR) = KMAT
        DO 605 L = 0,MAXAGE
          UNREC(L) = KREC*UINIT(L)
          REC(L)   = KREC*RINIT(L)
  605   CONTINUE

C       Forward project the population as necessary
        DO 610 IYR = INITYR,IPPROT-1
          CATCH(IYR) = CHIST
          CALL STKUPA (IYR,CATCH(IYR))
  610   CONTINUE
C       Protection period (CATCH(IYR) was set to 0 above)
        C = 0.D0
        DO 620 IYR = IPPROT,-1
          CALL STKUPA (IYR,C)
  620   CONTINUE

C       Has convergence taken place?
        IF ((ABS(PTRUE(0)/P0-1.D0).GE.0.0001D0).or.(icount.lt.500)) THEN
          IF ((PTRUE(0) - P0) .GT.0.D0) THEN
            KMAX = K1
          ELSE
            KMIN = K1
          ENDIF
C         Update loop count and check is not wasting time
          ICOUNT = ICOUNT + 1
c          IF (ICOUNT.GT.500) STOP '**** ERROR: SETKA HAS NOT CONVERGED'

          GO TO 600

        ELSE
C ---     Routine has converged
          D = PTRUE(0)/K1
        ENDIF
      ENDIF

C     Store the pristine 1+ population in KSURV, & year -1 size in PSURV1
      KSURV  = KREC * R1PLUS
      K1P = KSURV
      PSURV1 = PSURV(-1)
 
C     Reset SUR if necessary i.e. in trials where epidemics occur after year 0
C     The factor 0.017 is chosen such that if no catches are taken, the total
C     population will, on average, maintain the premanagement productivity.  
C     The value of this factor IS NOT GENERAL and is only correct for base 
C       case mortality rates and a 0.02 probability of epidemics.
      IF (ERATE.GT.0.D0) THEN
        DO 800 L = 0,MAXAGE
          SUR(L) = SURV(MORT1-0.017D0,MORT2-0.017D0,L,MORTIP)
  800   CONTINUE
      ENDIF
C
C     Set up the carrying capacity in each year (future)  
      IF (K99AEP .EQ. 0.d0) THEN
       DO 310 IYR = 1,NYEAR
        KT(IYR) = Kused
310    CONTINUE
C
      ELSE IF (K99AEP .GT. 0.d0) THEN
C      K changes linearly from Kused in year KYEAR to Kused*K99AEP in year NYEAR-1
       KINC = Kused * (K99AEP - 1.d0) / DBLE(NYEAR-1-KYEARAEP)
       DO 321 IYR = 1,KYEARAEP
        KT(IYR) = Kused
321    CONTINUE
       DO 322 IYR = KYEARAEP+1,NYEAR
        KT(IYR) = KT(IYR-1) + KINC
322    CONTINUE                 
       
      ELSE
C      K cyclic. K99AEP=-1; minima in years 0 & 100; max = 3K0
C                K99AEP=-2; maxima is years 0 & 100; min - K0/3
       PIFACT = 2.D0 * 3.141593D0 / DBLE(NYEAR)
       IF (K99AEP .EQ. -1.D0) THEN
        DO 331 IYR = 1,NYEAR
         KT(IYR) = Kused * (2.D0 - COS(PIFACT * IYR)) 
  331   CONTINUE
       ELSE IF (K99AEP .EQ. -2.D0) THEN
        DO 332 IYR = 1,NYEAR
         KT(IYR) = Kused/3.D0 * (2.D0 + COS(PIFACT * IYR))
  332   CONTINUE
       ENDIF      

      ENDIF
      
      RETURN
      END

C *********************************************************************
C *********************************************************************

      SUBROUTINE STKUPA (IYR,C)

C     Age structured model only.
C     STKUPA updates the stock size at the start of IYR+1, i.e. it
C     updates the REC & UNREC arrays and sets PTRUE(IYR+1) & PSURV
C     where PSURV is the 1+ population in year IYR (& not IYR+1).

      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC

      COMMON /AGEPAR/ MAT1,MSIG,REC1,RSIG,FEC,SUR(0:60),UNREC(0:60),
     +                REC(0:60),RECF(0:60),FMATUR(0:60),MAXAGE,MINMAT
      DOUBLE PRECISION MAT1,MSIG,REC1,RSIG,FEC,SUR,UNREC,REC,RECF,FMATUR
      INTEGER MAXAGE,MINMAT
      
      COMMON / TIMEVAR/ AT(-65:2000),KT(-65:2000),ZT(-65:2000),
     +                  K99AEP,KYEARAEP,MSYR99AEP,MSYRYRAEP,ISTEPAEP
      DOUBLE PRECISION AT,KT,ZT,K99AEP,MSYR99AEP
      INTEGER KYEARAEP,ISTEPAEP,MSYRYRAEP

      COMMON /HITFIT/ OPTDD,OPTF,OPTMSYL,FMSY,SMSYL,RFTOT0,RFEXP0,
     +       FMSYL,FECMSY,MSYLT1,MSYLT0,MSYLT2,AMSYR1,AMSYL1,
     +       AMSYR0,AMSYL0,AMSYR2,AMSYL2
      DOUBLE PRECISION FMSY,SMSYL,RFTOT0,RFEXP0,FMSYL,FECMSY
      DOUBLE PRECISION AMSYR0,AMSYL0,AMSYR1,AMSYL1,AMSYR2,AMSYL2
      DOUBLE PRECISION MSYLT0,MSYLT1,MSYLT2
      INTEGER OPTDD,OPTF,OPTMSYL
      
      DOUBLE PRECISION C,BIRTHS,PMAT,PEXP,URM,DENST,P1P,PREC
      INTEGER IYR,L

C     Sum the numbers of exploitable and 1+ whales
      PEXP = REC(0)
      PSURV(IYR) = 0.D0
      DO 10 L = 1,MAXAGE
        PEXP = PEXP + REC(L)
        PSURV(IYR) = PSURV(IYR) + REC(L) + UNREC(L)
   10 CONTINUE

C     Remove the catches for this year by calculating average 
C     survivorship after catching and multiplying over each age. 
      IF (C .LE. PEXP) THEN
        URM = 1.D0 - C/PEXP
      ELSE
        URM = 0.D0
        C = PEXP
      ENDIF

C     First advance the pooled age class (gets different treatment)
      REC(MAXAGE) = URM * (REC(MAXAGE-1)*SUR(MAXAGE-1) +
     1                     REC(MAXAGE)*SUR(MAXAGE) )

C     Loop advances the remaining age classes, working from oldest down.
C     Accumulate PTRUE(IYR+1)= No. mature (defined as no. recruited if OPTMOD=3) 
      PMAT  =  REC(MAXAGE)
      PREC  =  REC(MAXAGE)
      P1P   =  REC(MAXAGE) 
      IF (OPTMOD.EQ.3) THEN
        DO 40 L = MAXAGE-1,1,-1
          UNREC(L) = UNREC(L-1)*SUR(L-1)*(1.D0 - RECF(L))
          REC(L)   = SUR(L-1)*(REC(L-1)*URM+UNREC(L-1)*RECF(L))
          PMAT = PMAT + REC(L)
          P1P = P1P + UNREC(L) + REC(L)
          PREC = PREC + REC(L)
   40   CONTINUE
      ELSE
        DO 42 L = MAXAGE-1,1,-1
          UNREC(L) = UNREC(L-1)*SUR(L-1)*(1.D0 - RECF(L))
          REC(L)   = SUR(L-1)*(REC(L-1)*URM+UNREC(L-1)*RECF(L))
          PMAT = PMAT + (REC(L)+UNREC(L))*FMATUR(L)
          P1P = P1P + UNREC(L) + REC(L)
          PREC = PREC + REC(L)
   42   CONTINUE
      ENDIF
      PTRUE(IYR+1) = PMAT

C     Add new births to zero age-class (all defined to be unrecruited)
C     A different equation for BIRTHS is used unless OPTMOD=5 for
C     consistency with previous runs of the program (tuning in Oman 1998,
C     tuning in 2001 and the factorial expt).
C     (see correspondence for this decision)

      IF (OPTMOD.NE.5) THEN
        DENST = AT(IYR)*(1.D0 - (PTRUE(IYR)/KT(IYR))**ZT(IYR))
        BIRTHS = PTRUE(IYR)*FEC*(1.D0 + DENST)
      ELSE
        IF (OPTDD.EQ.0) THEN
         DENST = AT(IYR+1)*(1.D0 - (PREC/KT(IYR+1))**ZT(IYR+1))
        ELSEIF (OPTDD.EQ.1) THEN 
         DENST = AT(IYR+1)*(1.D0 - (P1P/KT(IYR+1))**ZT(IYR+1))
        ELSEIF (OPTDD.EQ.2) THEN
         DENST = AT(IYR+1)*(1.D0 - (PMAT/KT(IYR+1))**ZT(IYR+1))
        ELSE
         STOP
        ENDIF 
        BIRTHS = PTRUE(IYR+1)*FEC*(1.D0 + DENST)
      END IF
C
C ***
      IF (BIRTHS.LT.0.D0) BIRTHS = 0.D0
      UNREC(0) = BIRTHS*(1.D0-RECF(0))
      REC(0) = BIRTHS * RECF(0)
C
C     Update 1+ numbers
      PSURV(IYR+1) = 0.D0
      DO 53 L = 1,MAXAGE
        PSURV(IYR+1) = PSURV(IYR+1) + REC(L) + UNREC(L)
   53 CONTINUE

C     Halve the population if an epidemic occurs - only from IYR=0
      IF (IYR.GE.0 .AND. PROBE(IYR+1).EQ.1) THEN
        PTRUE(IYR+1) = 0.5D0*PTRUE(IYR+1)
        DO 50 L=0,MAXAGE
          UNREC(L) = UNREC(L) * 0.5D0
          REC(L) = REC(L) * 0.5D0
   50   CONTINUE
      ENDIF

      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE SETO (V,SIG,MEAN,MAXAGE)

C     SETO sets the proportion of each age class mature or recruited
C          into vector V

      DOUBLE PRECISION V(0:60),SIG,MEAN
      INTEGER MAXAGE,L

C     Loop over all ages
      DO 10 L = 1, MAXAGE-2
        IF (SIG.GT.0) THEN
         IF ((DBLE(L)-MEAN)/SIG.GT.10.D0) THEN
          V(L) = 1.D0
         ELSE IF ((DBLE(L)-MEAN)/SIG.LT.-10.D0) THEN
          V(L) = 0.D0
         ELSE
          V(L) = 1.D0/(1.D0 + EXP(-(DBLE(L)-MEAN)/SIG))
         ENDIF
        ELSE
         IF (DBLE(L).LT.MEAN) THEN
          V(L) = 0.D0
         ELSE
          V(L) = 1.D0
         ENDIF  
        ENDIF 
   10 CONTINUE
      V(0) = 0.D0
      V(MAXAGE-1) = 1.D0
      V(MAXAGE) = 1.D0

      RETURN
      END

C     ------------------------------------------------------------------

      FUNCTION SURV (MEAN0,MEAN1,A,MORTIP)

C     This function computes survival as function of age

      DOUBLE PRECISION SURV,MEAN0,MEAN1,BETA,ALPHA
      INTEGER A,MORTIP

      IF (MORTIP.LT.0) THEN
       BETA = (MEAN1-MEAN0)/16.D0
       ALPHA = MEAN0 - 4.D0*BETA
       IF (A.LE.4) THEN
          SURV = EXP(-MEAN0)
       ELSE
          SURV = EXP(-(ALPHA+BETA*A))
       ENDIF
      ELSE
       IF (A.LE.MORTIP) THEN
        SURV = MEAN0
       ELSE
        SURV = MEAN1
       ENDIF  
      ENDIF
      
      RETURN
      END

C     ------------------------------------------------------------------

      SUBROUTINE TRFORM (V,MAXAGE)

C     Adjust an ogive to transition form, that is so that V(L) =
C     the proportion of animals in a given class at age A-1 which make 
C     the transition to a different class age A

      DOUBLE PRECISION V(0:60),RM,D
      INTEGER MAXAGE,L

      RM   = V(0)
      DO 90 L = 1,MAXAGE
        IF (RM .LT. 1.D0) THEN
          D = RM
          RM = V(L)
          V(L) = (RM - D)/(1.D0 - D)
        ELSE
          RM = V(L)
          V(L) = 1.D0
        ENDIF
   90 CONTINUE

      RETURN
      END

C --------------------------------------------------------------------------
C --------------------------------------------------------------------------
C
      SUBROUTINE MSYPAR (OK,MORT1,MORT2,MORTIP,MSYR2)
C
C    Subroutine finds the resilience A1 & density dependent exponent Z
C    which give the required MSY at the specified MSYL.  It equates the
C    slope in the balancing per capita birthrates at two levels of
C    fishing mortality just above and below MSY, with the derivative of
C    per capita birth rate at MSYL. A root finding routine using Brent's
C    method is used to solve for the density dependent exponent.

      IMPLICIT NONE
      COMMON /HITFIT/ OPTDD,OPTF,OPTMSYL,FMSY,SMSYL,RFTOT0,RFEXP0,
     +       FMSYL,FECMSY,MSYLT1,MSYLT0,MSYLT2,AMSYR1,AMSYL1,
     +       AMSYR0,AMSYL0,AMSYR2,AMSYL2
      DOUBLE PRECISION FMSY,SMSYL,RFTOT0,RFEXP0,FMSYL,FECMSY
      DOUBLE PRECISION AMSYR0,AMSYL0,AMSYR1,AMSYL1,AMSYR2,AMSYL2
      DOUBLE PRECISION MSYLT0,MSYLT1,MSYLT2
      INTEGER OPTDD,OPTF,OPTMSYL
      
      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC
      
      COMMON /AGEPAR/ MAT1,MSIG,REC1,RSIG,FEC,SUR(0:60),UNREC(0:60),
     +                REC(0:60),RECF(0:60),FMATUR(0:60),MAXAGE,MINMAT
      DOUBLE PRECISION MAT1,MSIG,REC1,RSIG,FEC,SUR,UNREC,REC,RECF,FMATUR
      INTEGER MAXAGE,MINMAT
      
      EXTERNAL ZCALC,SURV,SMSYR
      DOUBLE PRECISION UF(2), UT(2), UE(2), UNRC(0:2000), RC(0:2000), 
     + PARS(4),DF,TOL,FM,URMSY,DRDF,UDMSY,DDDF,ZLO,ZHI,F1,F2,
     + MORT1,MORT2,R1PLUS,RMAT,RREC,PA,SURV,SMSYR,MSYR2
      INTEGER I,J,IERR,OK,L,MORTIP
      LOGICAL GOOD
C
C     Set tolerances
      DF = 0.0001D0
      TOL = 0.0001D0
C
C     Small (but nasty) conversion from Operating Model variables to Hitter-
C     Fitter variables
      FMSY = MSYR2
      SMSYL = MSYL
C
C     Note SUR array set in RESET or in conditioning (in GENPAR & HITPEG)
C     Set up maturity ogive: FMATUR =proportion of age A which are mature {A5}
C     Knife edge at age MAT50 in this version
      CALL SETO (FMATUR,MSIG,MAT1,MAXAGE)
      DO 7 L = 1,MINMAT
       FMATUR(L) = 0  
7     CONTINUE
C       
C     Set up recruitment ogive RECF: knife edge at age 1  {Eqn A4.1}
      CALL SETO (RECF,RSIG,REC1,MAXAGE)
C
C     Set up the survival vector
      DO 8 L = 0,MAXAGE
        SUR(L) = SURV(MORT1,MORT2,L,MORTIP)
    8 CONTINUE
C    
C *** Add the relative mature & recruited pristine population sizes starting
C     with unity in the zero age class (L=0)  (PA=relative no.in Lth age class)
      RMAT = FMATUR(0)
      RREC = RECF(0)
      R1PLUS = 0.D0
      PA   = SUR(0)
      RFTOT0 = 0.D0
      RFEXP0 = RECF(0)
      DO 9 L = 1,MAXAGE-1
        RMAT = RMAT + PA*FMATUR(L)
        RREC = RREC + PA*RECF(L)
        R1PLUS = R1PLUS + PA
        RFTOT0 = RFTOT0 + PA
        RFEXP0 = RFEXP0 + PA*RECF(L)
        PA   = PA*SUR(L)
    9 CONTINUE

C     Adjust for last age class being pooled (and fully recruited / mature)
      PA   = PA/(1.D0 - SUR(MAXAGE))
      RMAT = RMAT + PA
      RREC = RREC + PA
      R1PLUS = R1PLUS + PA
      RFTOT0 = (RFTOT0+PA)
      RFEXP0 = (RFEXP0+PA)

C     Set the birth rate FEC so as to give balance at equilibrium
C     FEC = # of age 0 / # mature in pristine population
      FEC = 1.d0 / RMAT
C
C     Set up the recruitment ogive in transition form:
C     set RECF =fraction of unrecruited animals of age A which recruit
C     at age A+1, except RECF(0) = fraction recruited of age 0
      CALL SETO  (RECF,RSIG,REC1,MAXAGE)
      CALL TRFORM(RECF,MAXAGE)
C      
C     If OPTMOD=4 MSYR2 contains MSY/1+population at MSYL. 
C     SMSYR sets the true MSYR (= MSY/Recruited population at MSYL)
      IF (OPTMOD.EQ.4) FMSY = SMSYR(FMSY,RECF,SUR,MAXAGE)
C
C  Set two levels of survival after fishing to closely straddle MSY
      UF(1) = 1.D0 - FMSY + DF*.5D0
      UF(2) = UF(1) - DF
C
C  Set up equilibrium population age structure under each F
      UNRC(0) = 1.D0 - RECF(0)
      RC(0)   = RECF(0)
C
      DO 50 I = 1,2
C
        IF (OPTF.EQ.1) THEN
          UNRC(0) = 1
          RC(0) = 0
          UNRC(1) = 0
          RC(1)   = SUR(0)*UNRC(0)
          DO 10 J = 2,MAXAGE
            UNRC(J) = 0
            RC(J)   = SUR(J-1)*UF(I)*RC(J-1)
 10       CONTINUE
        ELSE IF (OPTF.EQ.0) THEN
          DO 20 J = 1,MAXAGE
            UNRC(J) = SUR(J-1)*UNRC(J-1)*(1.D0 - RECF(J))
            RC(J)   = SUR(J-1)*(RC(J-1)*UF(I) + UNRC(J-1)*RECF(J))
 20       CONTINUE
        ELSE IF (OPTF.EQ.2) THEN
          DO 30 J = 1,MAXAGE
            FM = 1.D0 - (1.D0 - UF(I))*FMATUR(J-1)
            UNRC(J) = SUR(J-1)*FM*UNRC(J-1)*(1.D0 - RECF(J))
            RC(J)   = SUR(J-1)*FM*(RC(J-1) + UNRC(J-1)*RECF(J))
 30       CONTINUE
        ENDIF
        RC(MAXAGE) = RC(MAXAGE)/(1. - SUR(MAXAGE)*UF(I))
C
C       Calculate mature, total1+ and recruited totals
        UF(I) = FMATUR(0)
        UT(I) = 0.D0
        UE(I) = RC(0)
        DO 40 J = 1,MAXAGE
          UT(I) = UT(I) + RC(J) + UNRC(J)
          UE(I) = UE(I) + RC(J)
          UF(I) = UF(I) + (RC(J) + UNRC(J))*FMATUR(J)
 40     CONTINUE
        IF (OPTMOD.EQ.3) UF(I) = UE(I)
C
 50   CONTINUE
C
C     Calculate the 'standard recruited' relative population and derivative
      IF (OPTF.EQ.1) THEN
        URMSY = .5D0*(UT(1) + UT(2))
        DRDF  = (UT(2) - UT(1))/DF
      ELSE IF (OPTF.EQ.0) THEN
        URMSY = .5D0*(UE(1) + UE(2))
        DRDF  = (UE(2) - UE(1))/DF
      ELSE IF (OPTF.EQ.2) THEN
        URMSY = .5D0*(UF(1) + UF(2))
        DRDF  = (UF(2) - UF(1))/DF
      ENDIF
C
C     Calculate the density dependent relative population and derivative
      IF (OPTDD.EQ.1) THEN
        UDMSY = .5D0*(UT(1) + UT(2))
        DDDF  = (UT(2) - UT(1))/DF
      ELSE IF (OPTDD.EQ.0) THEN
        UDMSY = .5D0*(UE(1) + UE(2))
        DDDF  = (UE(2) - UE(1))/DF
      ELSE IF (OPTDD.EQ.2) THEN
        UDMSY = .5D0*(UF(1) + UF(2))
        DDDF  = (UF(2) - UF(1))/DF
      ENDIF
C
C  Calculate the equilibrium MSYL for the density dependent component 
C    which corresponds to the MSYL input for the specified component
C    (The values of UT, UE and UF from the two levels of fishing
C     mortality are implicitly averaged in the following expressions,
C     to give the required values at MSYL.  This saves doing all this
C     again for the actual FMSY)
      IF (OPTMSYL.EQ.1) THEN
        IF (OPTDD.EQ.1) THEN
          FMSYL = SMSYL
        ELSE IF (OPTDD.EQ.0) THEN
          FMSYL = SMSYL*(UE(1)+UE(2))*RFTOT0/(RFEXP0*(UT(1)+UT(2)))
        ELSE IF (OPTDD.EQ.2) THEN
          FMSYL = SMSYL*(UF(1)+UF(2))*RFTOT0*FEC/(UT(1)+UT(2))
        ENDIF   
        MSYLT1 = SMSYL
        MSYLT0 = SMSYL*(UE(1)+UE(2))*RFTOT0/(RFEXP0*(UT(1)+UT(2)))
        MSYLT2 = SMSYL*(UF(1)+UF(2))*RFTOT0*FEC/(UT(1)+UT(2))
      ELSE IF (OPTMSYL.EQ.0) THEN
        IF (OPTDD.EQ.1) THEN
          FMSYL = SMSYL*(UT(1)+UT(2))*RFEXP0/(RFTOT0*(UE(1)+UE(2)))
        ELSE IF (OPTDD.EQ.0) THEN
          FMSYL = SMSYL
        ELSE IF (OPTDD.EQ.2) THEN
          FMSYL = SMSYL*(UF(1)+UF(2))*RFEXP0*FEC/(UE(1)+UE(2))
        ENDIF
        MSYLT1 = SMSYL*(UT(1)+UT(2))*RFEXP0/(RFTOT0*(UE(1)+UE(2)))
        MSYLT0 = SMSYL
        MSYLT2 = SMSYL*(UF(1)+UF(2))*RFEXP0*FEC/(UE(1)+UE(2))
      ELSE IF (OPTMSYL.EQ.2) THEN
        IF (OPTDD.EQ.1) THEN
          FMSYL = SMSYL*(UT(1)+UT(2))/(FEC*RFTOT0*(UF(1)+UF(2)))
        ELSE IF (OPTDD.EQ.0) THEN
          FMSYL = SMSYL*(UE(1)+UE(2))/(FEC*RFEXP0*(UF(1)+UF(2)))
        ELSE IF (OPTDD.EQ.2) THEN
          FMSYL = SMSYL
        ENDIF   
        MSYLT1 = SMSYL*(UT(1)+UT(2))/(FEC*RFTOT0*(UF(1)+UF(2)))
        MSYLT0 = SMSYL*(UE(1)+UE(2))/(FEC*RFEXP0*(UF(1)+UF(2)))
        MSYLT2 = SMSYL
      ENDIF
C
      PARS(1) = FMSY
      UF(1) = 1.D0/UF(1)
      UF(2) = 1.D0/UF(2)
      FECMSY = (UF(1) + UF(2))*.5D0
      PARS(2) = (UF(2) - UF(1))/(DF*(FECMSY - FEC))
      PARS(3) = DRDF/URMSY - DDDF/UDMSY
C
C  Find exponent
      ZLO =  -5.D0
      ZHI =   5.D0
      CALL ZBRAC (ZCALC, PARS, 4, ZLO, ZHI, F1, F2, GOOD)
      IF (GOOD) THEN
        CALL ZBRENT (ZCALC, PARS, 4, ZLO, ZHI, F1, F2, TOL, IERR)
      ENDIF
C
C     Density dependent exponent and Resilience
      Z = PARS(4)
      A1 = (FECMSY/FEC - 1.D0)/(1.D0 - FMSYL**Z)
C
C     Convert among the various MSYL, MSYRs 
      CALL ALTMSY (1, IERR, AMSYR1, AMSYL1)
      CALL ALTMSY (2, IERR, AMSYR0, AMSYL0)
      CALL ALTMSY (3, IERR, AMSYR2, AMSYL2)
C      
      RETURN
      END
C
C  ********************************************************************
C
      FUNCTION ZCALC(PARS,IERR)
C
C  The root of this function is the density dependent exponent
C
      IMPLICIT NONE
      
      COMMON /HITFIT/ OPTDD,OPTF,OPTMSYL,FMSY,SMSYL,RFTOT0,RFEXP0,
     +       FMSYL,FECMSY,MSYLT1,MSYLT0,MSYLT2,AMSYR1,AMSYL1,
     +       AMSYR0,AMSYL0,AMSYR2,AMSYL2
      DOUBLE PRECISION FMSY,SMSYL,RFTOT0,RFEXP0,FMSYL,FECMSY
      DOUBLE PRECISION AMSYR0,AMSYL0,AMSYR1,AMSYL1,AMSYR2,AMSYL2
      DOUBLE PRECISION MSYLT0,MSYLT1,MSYLT2
      INTEGER OPTDD,OPTF,OPTMSYL
      
      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC
      
      COMMON /AGEPAR/ MAT1,MSIG,REC1,RSIG,FEC,SUR(0:60),UNREC(0:60),
     +                REC(0:60),RECF(0:60),FMATUR(0:60),MAXAGE,MINMAT
      DOUBLE PRECISION MAT1,MSIG,REC1,RSIG,FEC,SUR,UNREC,REC,RECF,FMATUR
      INTEGER MAXAGE,MINMAT
      
      DOUBLE PRECISION PARS(4),ZCALC
      INTEGER IERR

      IERR = 0
      IF (PARS(4) .EQ. 0.D0) THEN
        ZCALC = 1.D0 + PARS(1)*(PARS(3) - PARS(2)*DLOG(FMSYL))
      ELSE
        ZCALC = 1.D0 + PARS(1)*(PARS(3) - PARS(2)*
     1                           (FMSYL**(-PARS(4)) - 1.D0)/PARS(4))
      ENDIF
      RETURN
      END
C
C **********************************************************************
C
      SUBROUTINE ALTMSY (ICOMP, IERR, AMSYR, AMSYL)
C
C  Subroutine finds the MSYR and MSYL for a population component
C    for specified values of the resilience and exponent
C
      IMPLICIT NONE

      COMMON /HITFIT/ OPTDD,OPTF,OPTMSYL,FMSY,SMSYL,RFTOT0,RFEXP0,
     +       FMSYL,FECMSY,MSYLT1,MSYLT0,MSYLT2,AMSYR1,AMSYL1,
     +       AMSYR0,AMSYL0,AMSYR2,AMSYL2
      DOUBLE PRECISION FMSY,SMSYL,RFTOT0,RFEXP0,FMSYL,FECMSY
      DOUBLE PRECISION AMSYR0,AMSYL0,AMSYR1,AMSYL1,AMSYR2,AMSYL2
      DOUBLE PRECISION MSYLT0,MSYLT1,MSYLT2
      INTEGER OPTDD,OPTF,OPTMSYL
      
      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC
      
      COMMON /AGEPAR/ MAT1,MSIG,REC1,RSIG,FEC,SUR(0:60),UNREC(0:60),
     +                REC(0:60),RECF(0:60),FMATUR(0:60),MAXAGE,MINMAT
      DOUBLE PRECISION MAT1,MSIG,REC1,RSIG,FEC,SUR,UNREC,REC,RECF,FMATUR
      INTEGER MAXAGE,MINMAT
      
      EXTERNAL OMSY
      DOUBLE PRECISION DF,TOL,PARS(4),FLO,FHI,AMSYR,AMSYL,F1,F2
      INTEGER ICOMP,IERR
      LOGICAL GOOD
C
      DF = 0.1D-3
      TOL = 0.1D-6
C
C  Find exponent
      IF (FMSY .GE. 1.D-8) THEN
        PARS(2) = A1
        PARS(3) = ICOMP
        FLO =  FMSY*.5D0
        FHI =  FMSY*1.5D0
        CALL ZBRAC (OMSY, PARS, 1, FLO, FHI, F1, F2, GOOD)
        IF (GOOD) THEN
          CALL ZBRENT (OMSY, PARS, 1, FLO, FHI, F1, F2, TOL, IERR)
        ENDIF
      ELSE
        PARS(1) = 0.D0
        PARS(4) = SMSYL
      ENDIF
      AMSYR = PARS(1)*100.0
      AMSYL = PARS(4)
C
      RETURN
      END
C
C  ********************************************************************
C
      FUNCTION OMSY (PARS, IERR)
C
C  Function calculates used in calculating MSY and MSYL for a designated 
C    component of the population, given values for resilience A and 
C    exponent beta.
C
      IMPLICIT NONE

      COMMON /HITFIT/ OPTDD,OPTF,OPTMSYL,FMSY,SMSYL,RFTOT0,RFEXP0,
     +       FMSYL,FECMSY,MSYLT1,MSYLT0,MSYLT2,AMSYR1,AMSYL1,
     +       AMSYR0,AMSYL0,AMSYR2,AMSYL2
      DOUBLE PRECISION FMSY,SMSYL,RFTOT0,RFEXP0,FMSYL,FECMSY
      DOUBLE PRECISION AMSYR0,AMSYL0,AMSYR1,AMSYL1,AMSYR2,AMSYL2
      DOUBLE PRECISION MSYLT0,MSYLT1,MSYLT2
      INTEGER OPTDD,OPTF,OPTMSYL
      
      COMMON /STKVRS/ PTRUE(-65:2000),MSYL,MSYR1,Z,KSURV,
     +         PSURV(-65:2000),
     +         PSURV1,ERATE,CATERR,A1,K1,K1P,PROBE(0:2000),OPTMOD,OPTC
      DOUBLE PRECISION PTRUE,MSYL,MSYR1,Z,KSURV,PSURV,PSURV1,ERATE,
     +        CATERR,A1,K1,K1P
      INTEGER PROBE,OPTMOD,OPTC
      
      COMMON /AGEPAR/ MAT1,MSIG,REC1,RSIG,FEC,SUR(0:60),UNREC(0:60),
     +                REC(0:60),RECF(0:60),FMATUR(0:60),MAXAGE,MINMAT
      DOUBLE PRECISION MAT1,MSIG,REC1,RSIG,FEC,SUR,UNREC,REC,RECF,FMATUR
      INTEGER MAXAGE,MINMAT
      
      DOUBLE PRECISION UF(2), UT(2), UE(2), UNRC(0:2000), RC(0:2000), 
     +    PARS(4),C1,URMSY,DRDF,UDMSY,DDDF,C2,DF,TOL,F,A,OMSY,OMSYL,FM
      INTEGER IERR,TYPE,I,J
C
C     Tolerances 
      DF = 0.1D-5
      TOL = 0.1D-5
C
      IERR = 0
      F = PARS(1)
      IF (F.LE.0) F = 0
      A = PARS(2)
      TYPE = PARS(3)
C
C  PARS(1) = F, is the trial value of MSY
C  PARS(2) = A, is the given value of resilience
C  PARS(3) = COMP, determines which MSY is calculated 
C            (1. = Total, 2. = recruited, 3. = mature)
C  PARS(4) = MSYL for exploitation of the designated component
C
C  Set two levels of survival after fishing to closely straddle MSY
      UF(1) = 1.D0 - F + DF*.5D0
      UF(2) = UF(1) - DF
C
C  Set up equilibrium population age structure under each F
      UNRC(0) = 1.D0 - RECF(0)
      RC(0)   = RECF(0)
C
      DO 50 I = 1,2
C
        IF (TYPE .EQ. 1.D0) THEN
          UNRC(0) = 1
          RC(0) = 0
          UNRC(1) = 0
          RC(1)   = SUR(0)*UNRC(0)
          DO 10 J = 2,MAXAGE
            UNRC(J) = 0
            RC(J)   = SUR(J-1)*UF(I)*RC(J-1)
 10       CONTINUE
        ELSE IF (TYPE .EQ. 2.D0) THEN
          DO 20 J = 1,MAXAGE
            UNRC(J) = SUR(J-1)*UNRC(J-1)*(1.D0 - RECF(J))
            RC(J)   = SUR(J-1)*(RC(J-1)*UF(I) + UNRC(J-1)*RECF(J))
 20       CONTINUE
        ELSE IF (TYPE .EQ. 3.D0) THEN
          DO 30 J = 1,MAXAGE
            FM = 1.D0 - (1.D0 - UF(I))*FMATUR(J-1)
            UNRC(J) = SUR(J-1)*FM*UNRC(J-1)*(1.D0 - RECF(J))
            RC(J)   = SUR(J-1)*FM*(RC(J-1) + UNRC(J-1)*RECF(J))
 30       CONTINUE
        ENDIF
C
C  Adjust the oldest class for effect of pooling
        RC(MAXAGE) = RC(MAXAGE)/(1. - SUR(MAXAGE)*UF(I))
C
C  Calculate mature, total1+ and recruited totals
        UF(I) = 0.D0
        UT(I) = 0.D0
        UE(I) = RC(0)
        DO 40 J = 1,MAXAGE
          UT(I) = UT(I) + RC(J) + UNRC(J)
          UE(I) = UE(I) + RC(J)
          UF(I) = UF(I) + (RC(J) + UNRC(J))*FMATUR(J)
 40     CONTINUE
C
 50   CONTINUE
C
C     Calculate the 'standard recruited' relative population and derivative
      IF (TYPE .EQ. 1.D0) THEN
        URMSY = .5D0*(UT(1) + UT(2))
        DRDF  = (UT(2) - UT(1))/DF
      ELSE IF (TYPE .EQ. 2.D0) THEN
        URMSY = .5D0*(UE(1) + UE(2))
        DRDF  = (UE(2) - UE(1))/DF
      ELSE IF (TYPE .EQ. 3.D0) THEN
        URMSY = .5D0*(UF(1) + UF(2))
        DRDF  = (UF(2) - UF(1))/DF
      ENDIF
C
C     Calculate the density dependent relative population and derivative
      IF (OPTDD.EQ.1) THEN
        UDMSY = .5D0*(UT(1) + UT(2))
        DDDF  = (UT(2) - UT(1))/DF
      ELSE IF (OPTDD.EQ.0) THEN
        UDMSY = .5D0*(UE(1) + UE(2))
        DDDF  = (UE(2) - UE(1))/DF
      ELSE IF (OPTDD.EQ.2) THEN
        UDMSY = .5D0*(UF(1) + UF(2))
        DDDF  = (UF(2) - UF(1))/DF
      ENDIF
C
      UF(1) = 1.D0/UF(1)
      UF(2) = 1.D0/UF(2)
      FECMSY = (UF(1) + UF(2))*.5D0
      FMSYL = (1.D0 - (FECMSY - FEC)/(A*FEC))
      IF (FMSYL .GT. 0.D0) THEN
        FMSYL = FMSYL**(1.D0/Z)
      ELSE
        FMSYL = 1.D-2
      ENDIF
      C1 = DRDF/(URMSY+0.00001) - DDDF/(UDMSY+0.00001)
      C2 = (UF(2) - UF(1))/(DF*(FECMSY - FEC))*(FMSYL**(-Z) - 1.D0)/
     1                                 Z
      OMSY = 1.D0 + F*(C1 - C2)
C
      UF(1) = 1.D0/UF(1)
      UF(2) = 1.D0/UF(2)
C
      IF (TYPE .EQ. 1.D0) THEN
        IF (OPTDD.EQ.1) THEN
          OMSYL = FMSYL
        ELSE IF (OPTDD.EQ.0) THEN
          OMSYL = FMSYL/((UE(1)+UE(2))*RFTOT0/(RFEXP0*(UT(1)+UT(2))))
        ELSE IF (OPTDD.EQ.2) THEN
          OMSYL = FMSYL/((UF(1)+UF(2))*RFTOT0*FEC/(UT(1)+UT(2)))
        ENDIF
      ELSE IF (TYPE .EQ. 2.D0) THEN
        IF (OPTDD.EQ.1) THEN
          OMSYL = FMSYL/((UT(1)+UT(2))*RFEXP0/(RFTOT0*(UE(1)+UE(2))))
        ELSE IF (OPTDD.EQ.0) THEN
          OMSYL = FMSYL
        ELSE IF (OPTDD.EQ.2) THEN
          OMSYL = FMSYL/((UF(1)+UF(2))*RFEXP0*FEC/(UE(1)+UE(2)))
        ENDIF   
      ELSE IF (TYPE .EQ. 3.D0) THEN
        IF (OPTDD.EQ.1) THEN
          OMSYL = FMSYL/((UT(1)+UT(2))/(FEC*RFTOT0*(UF(1)+UF(2))))
        ELSE IF (OPTDD.EQ.0) THEN
          OMSYL = FMSYL/((UE(1)+UE(2))/(FEC*RFEXP0*(UF(1)+UF(2))))
        ELSE IF (OPTDD.EQ.2) THEN
          OMSYL = FMSYL
        ENDIF   
      ENDIF
C
      PARS(4) = OMSYL
C
      RETURN
      END

C *********************************************************************

      FUNCTION SMSYR (MSYR1P,RECF,SUR,MAXAGE)

C     Calculate MSYR given MSYR1P = MSY/1+population at MSYL

      DOUBLE PRECISION SMSYR,MSYR1P,RECF(0:60),SUR(0:60),UF,DIF,UFMIN,
     +       UFMAX,R1PLUS,R(0:60),UNR(0:60),U1PLUS
      INTEGER MAXAGE,ICOUNT,L

C     Set up unrecruited component relative to # of age 0 & store sum
      UNR(0) = 1.D0 - RECF(0)
      U1PLUS = 0.D0
      DO 5 L = 1,MAXAGE
        UNR(L) = UNR(L-1)*SUR(L-1)*(1.D0 - RECF(L))
        U1PLUS = U1PLUS + UNR(L)
   5  CONTINUE
      R(0)   = RECF(0)

      ICOUNT = 0
      UFMAX = 1.D0 - MSYR1P
      UFMIN = 0.8D0

   10 UF = (UFMIN + UFMAX) / 2.D0
      R1PLUS = 0.D0
      DO 8 L = 1,MAXAGE-1
        R(L)   = SUR(L-1)*(R(L-1)*UF+UNR(L-1)*RECF(L))
        R1PLUS = R1PLUS + R(L)
   8  CONTINUE
      R(MAXAGE) = SUR(MAXAGE-1) * (R(MAXAGE-1)*UF + UNR(MAXAGE-1))
      R(MAXAGE) = R(MAXAGE)/(1. - SUR(MAXAGE)*UF)
      R1PLUS = R1PLUS + R(MAXAGE)
      DIF = (1.D0-UF)*(R1PLUS+R(0))/(R1PLUS+U1PLUS) - MSYR1P
      IF (ABS(DIF).LT.0.00001D0) GO TO 100
      IF (DIF.GT.0.D0) THEN
        UFMIN = UF
      ELSE
        UFMAX = UF
      ENDIF

      ICOUNT = ICOUNT + 1
      IF (ICOUNT.GT.500) STOP ' **** ERROR: MSYR NOT FOUND'
      GO TO 10

  100 CONTINUE
C     MSYR found (where MSYR = MSY / Recruited population size at MSYL)
      SMSYR = 1.D0 - UF
      END
