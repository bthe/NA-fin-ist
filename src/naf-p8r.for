C     Last change:  CA    4 Dec 2015   10:26 am
C
C   NAF-TST   P8                              Dec 2014
C
C     To do: Generalise the CPUE code
C
C     Changes
C       CA  Density dependent Dispersal options removed 10/4/15
C       CA  Allow Subareas used in Management (CATKA & SIGHT arrays)
C           to differ from those used in the population model
C       CA  Put minimisation stuff into RUNMIN
C       CA  Changed management definitions & DIVCL to match NPMi (28/1/2015)
C       AEP Moved declining selex code to HITFUN (1/20/2014)
C       AEP Moved MSYPAR to HITFUN because selectivity is estimated (1/20/2014)
C
      MODULE ARYPAR
      IMPLICIT NONE

      CHARACTER(3) :: PROG="P8 "

      INTEGER,PARAMETER ::MXSTK=6,  MXAGE2=100,MXSUBA=7,   MXAR=12,
     +                    I1YR=-151,IENDYR=100,I1SYR =-30, MXMIX=2,
     +                    MXAB=20,  MXOP=40,   MXGAM =10,  MXFIT=18000

C     REF      Reference number of run: = Trial, variant & program version
C     MXSTK    Maximum number of stocks/substocks      = 6   J
C     MXAGE2   Maximum allowed pooled age-class        =100  L
C     MAXAGE   Pooled age-class, read in.
C     MXSUBA   Maximum number of subareas              = 7   K
C     MXAR     Maximum number of group areas           = 11 KA
C              i.e.assessment areas (= subareas or unions thereof)
C     I1YR     first year of CATCHM & CATCHF arrays    =-151 (1864)
C     IENDYR   last year of catch arrays               = 100
C     I1SYR    1st year of survey arrays (SGHT & CV)   =-30  (1985)
C     MXMIX    Maximum number of mixing matrix options =2

C *** CONDITIONING arrays sizes:
C     MXAB     Maximum no. of HISTORICAL abundance estimates =20
C     MXOP     Maximum no. of parameters = MXGAM+MXSTK+5   =29
C     MXGAM    Maximum no. of gamma parameters             =10
C     MXFIT    Maximum no. of fit components :MXABS+others =180

C     SUB-AREA NUMBERS:
C     1   2   3   4   5   6   7
C     EC, WG, EG, WI, EI, No, Sp  Hyp 1-5, read in data & CLA arrays
C     EC, WG, G+I,EI, No, Sp      Hyp 6-7
      CHARACTER*5 ANAM(0:MXSUBA),ANAM2(0:MXSUBA),STKNAM(MXSTK)*2,REF*15
      DATA ANAM /' ','EC','WG','EG', 'WI','EI','No','Sp'/
      DATA ANAM2/' ','EC','WG','EG+WI','EI','No','Sp','*'/

      END MODULE ARYPAR

C --------------------------------------------------------------------

      MODULE DECLARECOM

      USE ARYPAR

      REAL(8) PMATF(I1YR:IENDYR,MXSTK), PP(I1YR:IENDYR,MXSTK,MXSUBA),
     +     RM(0:MXAGE2,MXSTK),RF(0:MXAGE2,MXSTK), KMAT(MXSTK),FEC,G0,
     +     SUR(0:MXAGE2),
     +     RECF(0:MXAGE2),FMATUR(0:MXAGE2), A(MXSTK),Z(MXSTK),
     +     JVE(MXSUBA,MXSTK,MXMIX), JV(MXSUBA,MXSTK,MXMIX),
     +     V (MXSUBA,MXSTK,MXAGE2), VE(MXSUBA,MXSTK),
     +     SELN(I1YR:IENDYR,MXSUBA),D(MXSTK,MXSTK),AVDISP(MXSTK),
     +     INITDEP(MXSTK),CVADD(MXSUBA),SSAMPL,CPENAL,SELDEC
      INTEGER NSTK,NSUBA,INITYR,MAXAGE,OPTMIX,
     +     NMIX(MXSTK),MIXNO(I1YR:IENDYR+1,MXSTK),
     +     JVEGAM(MXSUBA,MXSTK,MXMIX),OPTDSP(MXSTK),IMIX1,IMIX2,IMIX3,
     +     FIRSTYR,OPTEARLY,NODISP
C
      COMMON /STKVRS/ PMATF,A,Z,KMAT,SUR,FEC,FMATUR,PP,
     +             RM,RF,RECF,V,VE,JV,JVE,G0,CVADD,
     +             SSAMPL,CPENAL,D,AVDISP,SELN,INITDEP,SELDEC,
     +             MIXNO,NMIX,INITYR,NSTK,NSUBA,OPTDSP,NODISP,MAXAGE,
     +             OPTMIX,IMIX1,IMIX2,IMIX3,FIRSTYR,OPTEARLY,JVEGAM
C
C     PMATF(IYR,J) Mature female population size of stock J at start
C                  of year IYR (prior to catch being set & removed).
C                  Pristine values are read from NAFCON.DAT or are
C                  calculated when applying the conditioning algorithm.
C     RM(L,J) / RF No. of males/females of age L, stock J in current year.
C                  Note: this version uses a selectivity vector SELM/SELF
C                        instead of keeping recruited & unrecruited arrays
C     RECF(L)      Fraction of animals of age L-1 which recruit at age L
C     A(J)         Resilience parameter for stock J. Set in MSYPAR.
C     Z(J)         Density-dependence exponent for stock J. Set in MSYPAR.
C     KMAT(J)      Pristine no. of mature females in stock J =PMATF(INITYR,J)
C                  Set in SETKA
C     SUR(L)       Annual survival rate at age L. Set in MSYPAR.
C                  Sex/stock independent.
C     FEC          "Fecundity" when the resource was pristine.
C     FMATUR(L)    Fraction of age L females that are mature. Set in MSYPAR
C                  MINMAT = minimum age of maturity.
C                  Maturity is knife edge @ MAT1 (age of 1st parturition) if MSIG=0
C     PP(IYR,J,K)  No. of 1+ of stock J in IYR, subarea K. Set in STKUPA
C                  (Subarea is defined using TRUE stock boundaries)
C     SELN(IYR,K)  Selected population in subarea K in year IYR (relevant to CPUE)
C     G0           G(0) Correction factor. Read in.
C     JVE(K,J,N)   The mixing matrices excluding the gamma factors.
C                  Read in, K=1,NSUBA, J=1,NSTK, N=1,NMIX(J)
C     JVEGAM(K,J,N)The gamma factors for the mixing matrices. Read in.
C                  If JVEGAM(K,J,N)=i then mixing matrix element
C                  JV(K,J,N) = JVE(K,J,N)*gamma(i).
C     JV(K,J,N)    The mixing matrices set using JVE & JVEGAM once the
C                  gamma values are set (in conditioning).
C     VE(K,J)      The expected mixing matrix = mean of JV (over N), unless
C                  OPTMIX>1 when VE = JV(K,J,N=1) = initial value. Set in HITFUN
C     V(K,J,L)     Mixing matrix in current year = fraction of stock J
C                  of age L in each subarea  (so sum V(K,J,L) over K = 1)
C     MIXNO(IYR,J) Mixing Random number for stock J in year IYR.
C                  Set in SETMIX; used in SETV.
C     CVADD(K)     Additional process error.  Read in (=0 for most trials)
C     SSAMPL       Increase in survey sampling intensity. =1 in basecase.
C                  In NF-P intensity is tripled (SSAMPL=3) by changing n, alpha**2
C                  & beta**2 to 3n, alpha**2/3 & beta**2/3. (CVadd unchanged).
C                  This is implemented by setting TAU2 = TAU2/3 & DF to 3DF
C                  RUNSUR is reset so EG and EI/F are not surveyed in future
C                  Not used if EG & WI subareas are combined (Hyp. 7 & 8)
C     OPTMIX       =0 for 1 mixing matrix, =1 for average of 2,
C                  =2 for linear change from matrix M1 in yr IMIX1 to M2 in IMIX2
C                   then constant (NF-G)
C                  =3 as for OPTMIX=2 until IMIX2; changes back to M1 in IMIX3 (NF-F)
C     SELDEC       trial NF-U: selectivity decreases by SELDEC /yr after age 8
C     NODISP       1: dispersal; =0 no dispersal (hyp. 4 & 8). Used to set up gammas
C     OPTDSP(J)    Estimate Average dispersal between stocks J & J+1. [J/=NSTK]
C     AVDISP(J)    Mean dispersal rate, set in HITFUN if OPTDSP(J)=1
C                  If OPTDSP=0 & AVDISP>0 (read in), dispersal is fixed - not in NAF
C                  Dispersal ONLY occurs between adjacent central substocks
C     D(J1,J2)     Dispersal rate from stock J1 to J2. Set in SETKA
C     INITDEP(J)   Initial depletion
C     CPENAL       historical catch that is not taken. Set during conditioning (Devin
C     NSTK         No. of stocks   (Read in) J=1,NSTK (substocks are treated as stocks)
C     NSUBA=6 or 7 No. of subareas (Read in) K=1,NSUBA
C     NMIX(J)      No. of mixing matrix options (=1 or 2)
C     INITYR       1st year of model when stocks are pristine
C     FIRSTYR      Year when model begins. Read in.
C     OPTEARLY     0:start model in 1864; 1:start in later year (1915 or 50) (NF-I & -J)

C Parameters used by MSYPAR (all read in)

      REAL(8) MSYR(MXSTK),MSYL,MORT(3),MAT1,MSIG
      INTEGER OPTF,OPMSYL,OPTDD,MINMAT
      COMMON /MSYVAR/ MSYR,MSYL,MORT,MAT1,MSIG,MINMAT,OPTF,OPMSYL,OPTDD

C     MSYR(J)  MSYR of stock J. Read in. J=1,NSTK
C     MSYL     MSYL (of all stocks).  Read in.
C     MAT1,MSIG,MINMAT Used to set maturity ogive. Read in.
C     MORT(II) Used in setting up survivorship vectors.  Read in.
C     OPTDD    Density dependence on: 0: 1+ population (not used) 1: mature
C     OPMSYL   MSYL is defined on:    0: 1+ population (not used) 1: mature
C     OPTF     Fishing pattern defining MSYR selectivity on: 0:1+ 1: mat.pop (default)

      INTEGER        PRDIAG,DOCON,OPTDET,OPTCPE,SELYR,NSEL,OPTSEL
      COMMON /OPTNS/ PRDIAG,DOCON,OPTDET,OPTCPE,SELYR,NSEL,OPTSEL
C
C     PRDIAG   Print diagnostics:  0=No,  1=Yes
C     OPTDET   Option controlling stochasticity.  (Read from MANAGE.DAT)
C              0 : Stochastic   1 : Deterministic
C     DOCON    Option controlling conditioning i.e. whether to estimate or
C              read in K, GAMMA & BKCON parameters  (Read from MANAGE.DAT)
C              -1: Check conditioning: Read data from NAFCON.DAT
C              0 : Read data from NAFCON.DAT
C              1 : Conditioning: start with parameters in COPYNP.DAT
C              2 : Conditioning: start deterministic run with values in NAFCON.DAT
C              3 : Conditioning: start each replicate with values in NAFCON.DAT
C                  beginning at trial NSCON
C              4 : Read from NAFCON.DAT. Print quota & tuning level, in real applcn
C     OPTCPE   =1 to use CPUE data (NF-C); =0 otherwise (Read in)
C     SELYR    If >0 selectivity changes in year SYEAR (trial NF-S)
C     NSEL     Defines which selectivity values are used; =1 unless
C              SELYR>0 in which case NSEL=2 from SELYR on (NF-S trial).
C     OPTSEL   =1 to estimate selectivity (baseline); =0 otherwise (Read in)
C
C
      REAL(8) CATT(I1YR:IENDYR,MXSTK),CATKA(I1YR:IENDYR,MXAR),
     +       SELM(0:MXAGE2,2),SELF(0:MXAGE2,2),TOTCAT,SIGSQ
      REAL(8), DIMENSION(I1YR:IENDYR,MXSUBA):: CATCHM,CATCHF,OCAM,OCAF
      REAL(8), DIMENSION(I1YR:IENDYR,1:MXAGE2,MXSUBA):: OCBYAM,OCBYAF,
     +                                                  PCBYAM,PCBYAF,
C     True observed values
     +                                                  TOCBYAM, TOCBYAF

      INTEGER ICBYA1,ICBYA2
      COMMON /CTCG/ CATCHM,CATCHF,CATT,CATKA,SELM,SELF,TOTCAT,
     +     OCBYAM,OCBYAF,PCBYAM,PCBYAF,OCAM,OCAF,ICBYA1,ICBYA2
C
C     CATCHM(IYR,K) True catch of males in IYR in Subarea K. Read in.
C                   CATCHF=females
C     CATT(IYR,J)   Catch from stock J in IYR. Set in STKUPA
C     CATKA(IYR,KA) Catch from area KA in IYR, 'best' series for use by CLC
C                   If NSUBA/=MXSUBA (Hyp. 7 & 8) the area definitions
C                   differ from those used in the population model.
C     TOTCAT        Total true historical catch.  Set in INIT.
C     SELM(L,NSEL)  Selectivity of males age L. Estimated (unless OPTSEL=0).
C                   NSEL=1 except if SELYR>IYR then NSEL=2 (NF-S trial)
C     OCAM(IYR,K)   Observed catch of males in subarea K in IYR known by age
C     OCBYAM(IYR,L,K) Observed proportion of catch of age L males in subarea K in IYR. Read in.
C     PCBYAM(IYR,L,K) Modelled catch of age L males in subarea K in IYR.
C     SELF,OCAF,OCBYAF,PCBYAF = female versions of above arrays
C     ICBYA1 & 2    1st & last years to fit catch at age data
C
      REAL(8) CONEST(MXAB),TAREST(MXAB),TARVAR(MXAB),CONVAR(MXAB),
     +        TAGM(40,0:MXAGE2,MXSTK),TAGF(40,0:MXAGE2,MXSTK)
      REAL(8) RECTAR(40,I1YR:0,MXSUBA),TAGU(40,I1YR:0,MXSUBA),NREL(40),
     +        NREL2(40),PSI(MXSUBA),TLOSS1,TLOSS2,TAGWT,CPEWT
      REAL(8) CPUE(6,30),CPUEV(4,4),CPUEsig(5:6),qCPUE(6)
      INTEGER TARYR(MXAB),TARA(MXAB),RECOBS(40,I1YR:0,MXSUBA)
      INTEGER IGenTag(40,90),KGenTag(40,90),TAGOP(I1YR:0,0:MXSUBA),
     +        TAGA(40),TAGIYR(40),TAGREV(40),
     +        TAGSS(40),TTYPE,NOP,NAB,NROP,ITAG1,OPTPSI,
     +        NGAMMA,NDELTA,NCPUE,Np(6),CPUEYR(6,30),CPUEK(6)
      COMMON /CONDIT/ CONEST,TAREST,TARVAR,CONVAR,TAGM,TAGF,RECTAR,TAGU,
     +     CPUE,CPUEV,CPUEsig,qCPUE,NREL2,TARYR,TARA,RECOBS,TAGSS,
     +     TAGOP,NREL,PSI,TLOSS1,TLOSS2,TAGWT,TAGA,TAGIYR,TTYPE,NOP,NAB,
     +     IGenTag,KGenTag,TAGREV,
     +     NROP,ITAG1,OPTPSI,NCPUE,Np,CPUEYR,CPUEK,NGAMMA,NDELTA,TAGEXC1
C
C *** CONDITIONING variables: Abundances: (all read in EXCEPT TAREST)
C     NOP        Number of conditioning parameters (max=MXOP):
C                =NGAMMA + NDELTA + 3 (tagging) + NSTK*2 + 4 or 8 (selectivity)
C     NGAMMA     No. of Gamma parameters (for mix matrix), normally 6 (Read in)
C     NDELTA     No. of Delta (dispersal) parameters (Read in)
C     CONEST(NA) historical abundance about which TAREST(NA) is
C                lognormally distributed.  NA=1,NAB
C     CONVAR(NA) Variance of CONEST(NA)
C     TAREST(NA) Target abundance estimate to fit to when conditioning.{Ref Ga}
C                In Hyp. 7&8 data for EG & WI is combined
C     TARVAR(NA) Target variance = CV**2.  CV is read in.
C     TARYR(NA)  Target Year corresponding to TAREST(NA)
C     TARA(NA)   Target subarea   "       "    "      - use model subarea
C                Is reset to -1 to skip the estimate in likelihood calculation
C                (because estimate is combined with WI value, Hyp 7&8).
C     TTYPE      Target type 1=standard  In hyp. 7 & 8, if TTYPE==-3 the
C                estimate is combined with previous one (which has TTYPE=-2)
C *** Tagging data: A release operation is a year/area when marking was done
C                eg WI/1965, EC/1966, EC/1967 etc
C                The likelihood includes a term for every release operation
C     ITAG1      1st year whales were tagged (Read in)
C     NROP       Number of release operations (Read In) NR=1,NROP
C     NREL(NR)   No. of marks released in release operation NR (Read in)
C     TAGA(NR)   Area of release NR (Read in) - used for printing
C     TAGIYR(NR) Year of release NR (Read in).
C       E.g. Release NR=10 is 5 tags in WG in 1967 so
C                NREL(10)=5, TAGA(10)=2, TAGIYR(10)=1967+ISCALE=IT & TAGOP(IT,2)=10
C *** The following hold tag data by MODEL area
C                ie with subareas 3 (EG) & 4 (WI) combined if Hyp = 7 or 8
C     NR2        = release operation no. after combining areas if necessary
C     TAGREV(NR) = revised release operation # for original operation NR
C                TAGREV(NR) = NR except for combined releases in Hyp 7 & 8
C     NREL2(NR2) = No. of marks released in revised release operation NR2
C                Set to 0 if operation NR2 is included with another one.
C     TAGSS(NR2) = Number of same season recoveries from release NR2
C     TAGOP(IT,K)= Release operation no. (NR2) for releases in area K year IT
C                = 0 if there were no releases in year IT & area K
C                TAGUPD uses NREL2 & TAGOP to set up the marked popln TAGM & TAGF
C     TAGM(NR2,L,J)/TAGF Current no. of stock J, age L m/f tagged in release NR2
C *** Tag recoveries
C     RECOBS(NR,IYR,K) No.of observed recoveries in area K in IYR for release NR
C                They were marked in area TAGA(NR) & year TAGIYR(NR). Read in
C     RECTAR(NR2,IYR,K) Target no. of recoveries in model area K in IYR, release NR2
C                In the deterministic run RECTAR = RECOBS unless Hyp=7 or 8
C     TAGU(NR2,IYR,K) Modelled no. of recoveries in model area K in IYR, release NR2
C                for comparison with target. (No same season recoveries).
C     IGenTag & KGenTag(NR,N) = Recovery Year & Area (0 if not recovered) of mark N
C                of release NR for use in generating tag recovery targets.
C     NOTREC	 Number of tags not recovered
C     NREC       Number of tags recovered
C     PSI(K)     Tag reporting rate parameter for subarea K
C     OPTPSI     =0 if PSI fixed; =1 if PSI-C estimated; =2 if both PSI estimated
C                Baseline: OPTPSI=1; old trials NF27 used OPTPSI=0
C     TLOSS1 & 2 =tag loss rate in year 1 & from year 2 on
C     TAGWT      Weighting of tag data in likelihood. =1 (basecase) or 10 (NF-W)
C *** CPUE data
C     NCPUE      Number of CPUE series
C     Np(NCPUE)  Number of data points in each CPUE series
C     CPUEyr(NCPUE,30) CPUE years for each series
C     CPUE(NCPUE,30)   CPUE data for each series
C     CPUEK(NCPUE)     which area for CPUE series
C     CPUEV(4,4)       variance-covariance matrix for CPUE series 1 to 4
C     CPUEsig(5:6)     variance for CPUE series 5 and 6
C     TAGEXC1    1=exlude all tags after the first year, otherwise include all tags except same season recoveries
C     CPEWT      Weighting of cpue data in likelihood. =1e-3 
      REAL(8) BIAS,KSGT(MXSUBA),CV2TAU(I1SYR:-1,MXSUBA),TAU2(MXSUBA),
     +        SRATIO(MXSUBA),VRATIO(MXSUBA),DF,EGWIP,EGWIV
      INTEGER RUNSUR(I1SYR:IENDYR,MXAR),MAP(MXAR,MXSUBA),
     +        MAPC(MXAR,MXAR),I2SYR
      COMMON /DATVRS/BIAS,KSGT,CV2TAU,TAU2,SRATIO,VRATIO,DF,EGWIP,
     +               EGWIV,RUNSUR,MAP,MAPC,I2SYR
C     BIAS          The bias of the abundance estimates.  (Read in)
C     KSGT(IYR,KM)  Pristine TRUE mean 1+ abundance in Management subarea KM.
C                   Set in RESET. Used in SURVEY to generate SIGHT
C     RUNSUR(IYR,KA) Set to conduct a survey in area KA in IYR. Read in.
C                   For KA>MXSUBA the surveys are combinations of subarea surveys
C                   Reset if SSAMPL>1 so EG and EI/F are not surveyed in future
C                   KA index is by management subarea
C     TAU(KM)       Used to generate abundance estimates KM=1,MXSUBA. Set in RESET
C     DF            No.of degrees of freedom in CHISQ: used to set CVs. Read in
C     CV2TAU(IYR,K) historical CV**2 used in calculating TAU2 (in RESET)
C                   Read from SURVEYN file. Very large CVs should be EXCLUDED:
C                   in NP minke CVs > 0.9 were replaced by 0.9 to set TAU.
C     SRATIO(K)     Used to set proportion of Future survey estimates that are N 60N
C                   in trial NF-Q & Variant V4. Basecase value = 1. Read in
C     VRATIO(K)     SE of future survey estimate in trial NF-Q & Variant V4. Read in
C     EGWIP & EGWIV Historical mean EG proportion in combined EG-WI area + variance.
C                   Used in Hyp 7&8 to split EG-WI area in future. Read in
C     I2SYR         Year of first generated survey
C     MAP(KA,K)     Definition of subareas making up Small, Medium & Combination
C                   areas (which are also survey areas). Read from NAFpar.PAR.
C                   Catches are set for subarea groupings defined by MAP. KA=1,MXAR
C     CSUB(IP,KA,K) specifies subareas from which to take the catch for Small Area KA
C     MAPC(KA,KS)   Map of Small Areas making up combination & medium areas
C     Note: Small Areas are equal to subareas or groups of subareas.
C           MXAR different subarea groupings are allowed.
C
      INTEGER         ISEED1,ISEED2,ISEED3,ISEED4,ISEED5,ISEED6
      COMMON /RANPAR/ ISEED1,ISEED2,ISEED3,ISEED4,ISEED5,ISEED6

C     Uses of the random numbers
C      1 - Abundance estimates (lognormal, Poisson, Gamma)
C      2 - Used in SURVEY in beta Generator for N60N estimates & EG-WI split
C      3 - CV of abundance estimates
C      4 - Conditioning: set abundance & tagging targets
C      5 - Mixing matrices
C      6 - Conditioning: set catch at age targets

      REAL(8) FITDAT(MXFIT),FITRSD(MXFIT),FITRSS(MXFIT),GBEST(MXOP),
     +        BESTLL,INITV(MXOP),WtAge
      REAL(8) LikeAb,LikeCat,LikeTag,LikeAge
      INTEGER PHASE,ISCALE,NFIT,COND1,ISEST(MXOP)
      COMMON /FITTING/ FITDAT,FITRSD,FITRSS,GBEST,BESTLL,INITV,
     +              WtAge,LikeAb,LikeCat,LikeTag,LikeAge,
     +              PHASE,ISCALE,NFIT,COND1,ISEST

C     FITDAT   Values achieved in conditioning. Init to 0 in MAIN
C     FITRSD   Residuals = contributions to likelihood
C     FITRSS   Residuals
C     GBEST    Stores the parameters from the best fit so far
C     BESTLL   Best fit identified so far (to keep things down)
C     PHASE    Estimation phase (fit is adjusted if PHASE<=2)
C     ISCALE   Scale factor = -(IPEND+1) used to scale years
C     NFIT     Total number of conditioning targets (=NAB+NPS+extras)
C     COND1    =1 to initialise conditioning (w. mean mix matrices); Otherwise=0
C     ISEST    =1 is the parameter is estimated
C     INITV    Parameter values at start of the current call
C     WTAGE    Weight for catch at age data in likelihood
C     LikeAb etc Likelihood components, for printing

C     Filenames
      CHARACTER(LEN=30)  COPYNA, RANDOMF, MANAGE, NAFCON, NAFPAR,SURVEYn
      END MODULE DECLARECOM
C
C **********************************************************************
C **********************************************************************
C **********************************************************************

      PROGRAM NAFTST

C
C   Program implementing western North Atlantic Fin whale trials.
C   The equation numbers and references are taken from the trial
C   specifications (Jan 2013 Workshop report, Annex C)
C
C   The program generates 'true' population trajectories using a
C   series of catches set by the management procedure.  'Measurements'
C   of the stock sizes are generated, incorporating random noise, for
C   input to the management procedure (or CLC).
C
C   The management procedure is a self-contained module called by this
C   program to set catch limits using the population 'measurements'.
C
C   Note: Define IYEAR=0 as first year of management.  NTRIAL (=100)
C   replicates are run for each case, each lasting NYEAR years.
C
C PARAMETERS -----------------------------------------------------------
C   CATCH(K)     Catch of males+females in current year from subarea K
C                set in DIVCL (using values set by RUNCLC) from IYR=0.
C   SIGHT(IYR,KA)Absolute abundance estimate for area KA in year IYR.
C                Set to -1 if no survey performed ie no data available
C   CV(IYR,KA)   The estimate of the CV of SIGHT(IYR,KA)
C
C   CATSM(IP,KA) Catch quota for Small Area(KA) & procedure IP as set by RUNCLC
C   NSCAT(IP,KA) No. of subareas over which quota for Small Area KA is split
C   NTRIAL       Number of trials.  Read in.  1-100
C   NYEAR        Number of years of management.  Read in. 1-100
C   F,STEP,V0    Vectors used when applying the conditioning algorithm.
C                F: Test values for gammas + initial stock sizes.
C   F1ST         Initial values (guesses) of F: NGAMMA gamma values +
C                NDELTA dispersals, initial stock sizes.   Read in.
C   ISTART,IPEND First & final years of historical catches.
C
C   IYR          Current year
C   N            Current trial number  (1,NTRIAL)
C   IPNT         Output file
C   IN           Parameter input file
C   IN2 - 4      Input files: 2: Random numbers  3: historical catch data
C                4: Future survey details
C
C
C DEFINITIONS ----------------------------------------------------------
C
      USE DECLARECOM
C
      COMMON /STORE/ IASESS, NCAP, ATYPE, ITYPE
      INTEGER IASESS(I1SYR:IENDYR),NCAP(2),ATYPE(2,MXAR),ITYPE
C
      REAL(8) CATCH(MXSUBA),CATSM(2,MXAR),
     +     RM50,RMSIG,RF50,RFSIG,R,CVIN,FIT,OLDFIT,XNORM,
     +     MA4(55),MA6(55),MA(55),CCM(MXSUBA),CCF(MXSUBA),CCK(MXSUBA),
     +     CV(I1SYR:IENDYR,MXAR),SIGHT(I1SYR:IENDYR,MXAR),STEP(MXOP),
     +     F(MXOP),F1ST(MXOP),FBEST(MXOP),V0(MXOP),RAN1,NOSEX,
     +     SE, CREG(20),STEP2(MXOP),S,CVCLA,ABPR,SCV,SS,WTAGEI
      INTEGER NTRIAL,NYEAR,IYR,ISTART,IPEND,N,J,K,I,I2,KA,KS,
     +     L,NA,IPNT,IN,IN2,IN4,ISEEDS(400,6),CSUB(2,MXAR,MXSUBA),
     +     INEXT4,INXTP4,INEXT6,INXTP6,INEXT,INXTP,ISEED,NTOT,
     +     NS,KT,N2,IP,NPROC,NSCON,CHK(MXSUBA),KA2,KM,TSS,NSCAT(2,MXAR),
     +     NCOL, NREG(20),NR,IT,N3,NREC,NOTREC,R2,II,
     +     IAGE,CMULT,OPTSGT,NR2
      REAL(8) TMPF(1:MXAGE2),TMPM(1:MXAGE2),CMRATIO,CFRATIO,MCHAMMER,
     +     SIF,SIM, SIGSQEST
      LOGICAL USEGAM(MXGAM)
      CHARACTER DESC*50,TYPEDS*45,CTYPE*1,VARNAM*3,SCLA*6,OPTCLC*1
      
      EXTERNAL HITFUN,XNORM,RAN1
C
      DATA IPNT/8/,IN/5/,IN2/10/,IN4/12/

C     read in file names from the command line
      INTEGER :: NUM_ARGS, IX
      CHARACTER(LEN=30), DIMENSION(:), ALLOCATABLE :: ARGS
      
      NUM_ARGS = COMMAND_ARGUMENT_COUNT()
      ALLOCATE(ARGS(NUM_ARGS))  ! I'VE OMITTED CHECKING THE RETURN STATUS OF THE ALLOCATION 
      
C     ASSIGN DEFAULT VALUES TO 
      COPYNA = 'copyna.dat'
      RANDOMF = 'random.num'
      SURVEYn = 'surveyn.dat'
      MANAGE = 'manage.dat'
C      CATBYSEX = 'CATBYSEX.DAT'
C      CATCHBYAGE = 'CATCHBYAGE.DAT'
      NAFCON = 'nafcon.dat'
      NAFPAR = 'nafpar.par'

      IF (NUM_ARGS>0) THEN
      IX = 1  
      DO WHILE(IX < NUM_ARGS)
         CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
         IX = IX + 1
         IF(IX == NUM_ARGS) CONTINUE
         SELECT CASE(ADJUSTL(ARGS(IX-1)))
         CASE("-main")
            CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
            COPYNA = ARGS(IX)
         CASE("-rnd")
            CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
            RANDOMF = ARGS(IX)
         CASE("-man")
            CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
            MANAGE = ARGS(IX)
         CASE("-con")
            CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
            NAFCON = ARGS(IX)
         CASE("-par")
            CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
            NAFPAR = ARGS(IX)
         CASE("-survey")
            CALL GET_COMMAND_ARGUMENT(IX,ARGS(IX))
            surveyn = ARGS(IX)
         END SELECT
      END DO
      END IF

C
C     OPEN INPUT FILES
      OPEN (IN, FILE=COPYNA)
      OPEN (IN2,FILE=RANDOMF)
      OPEN (IN4,FILE=surveyn)                                           # Survey data & timetable
      OPEN (16, FILE=MANAGE,STATUS='OLD')                               # Conditioning opt.
C xxx UNIT  17, FILE=tag data     (Name set below)                      # Tags
C xxx UNIT  18, FILE=catch data   (Name set below)                      # Catch data
      OPEN (19, FILE='catbysex.dat',STATUS='OLD')                       # CATCHES KNOWN BY SEX
C xxx OPEN (20, FILE='CPUE.DAT')  (Opened below if OPTCPE>0)            # CPUE DATA
C xxx OPEN (21, FILE=CLC parameters (Opened in INITNA)                  # CLC parameters
      OPEN (22, FILE=NAFPAR)                                            # Management option
C xxx UNIT  23, FILE='NAFCON.DAT'= Conditioning parameters opened below
      OPEN (24, FILE='catchbyage.dat')
      
C     Read name of catch file
C     read in the first line from copyna.dat to name all the "naf" files
      READ (IN,'(A9,1X,A /)') REF,DESC
c restest files opened when varnam is read in
C      OPEN(IPNT,FILE=ref(1:7)//'.restest')
      OPEN (95,FILE=ref(1:7)//'.all')                                           # Conditioning parameters output
C     NAF.ALL is written when conditioning; it is read in as NAFCON.DAT for projections
      OPEN (94,FILE=ref(1:7)//'.age')                                           # age data
C ### OPEN (96,FILE='ref(1:7)//.CPE')                                           # CPUE (if OPTCPE>0)
      OPEN (97,FILE=ref(1:7)//'.mix')                                           # mix matrix
      OPEN (99,FILE=ref(1:7)//'.pop')                                           # matF by stock, 1+ by subarea
      OPEN (98,FILE=ref(1:7)//'.cat')                                           # catch by stck, 1+ by subarea
      OPEN (34,FILE=ref(1:7)//'.his')                                           # Values achieved in conditioning
      OPEN (35,FILE=ref(1:7)//'.rsd')                                           # Residuals =likelihood contributions
      OPEN (36,FILE=ref(1:7)//'.rss')                                           # Residuals
      OPEN (37,FILE=ref(1:7)//'.tar')                                           # Target values
      OPEN (38,FILE=ref(1:7)//'.sum')                                           # deterministic summary
      OPEN (39,FILE=ref(1:7)//'.tag')                                           # deterministic tag summary
      OPEN (41,FILE=ref(1:7)//'.sur')                                           # SURVEY VALUES
      OPEN (42,FILE=ref(1:7)//'.chk')                                           # Minimisation Check file.
      OPEN (89,FILE=ref(1:7)//'.cla')                                           # QUOTA OUTPUT FOR CHECKING
C ### OPEN (40,FILE='NAF.APP') Used in real application (SAG 2010) if DOCON=4

C     Set PRDIAG =1 here to print input catch & abundance data summaries.
      PRDIAG = 0

C **********************************************************************
C
C First read in the various files                                          # BASIC INPUT
C
C *********************************************************************
C
C     Read in true population parameters and option settings

      READ (IN,'((44X,I4))') NSTK,NTRIAL,NYEAR,NSUBA,ISTART,IPEND,
     +         FIRSTYR,OPTMIX,OPTCPE,OPTF,OPMSYL,OPTDD,MINMAT,MAXAGE

      IF (MAXAGE>MXAGE2) STOP 'ERROR: MXAGE2 needs increasing'
C
C     Premanagement catches were taken between actual years ISTART to IPEND.
C     Define IPEND+1 as year 0 = 1st year of management. Set scaling factor
      OPTEARLY = 0
      ISCALE = - (IPEND + 1)
      INITYR = ISTART + ISCALE
      FIRSTYR = FIRSTYR+ISCALE
      IF (INITYR.NE.FIRSTYR) OPTEARLY = 1
      IF (INITYR.LT.I1YR .OR. INITYR.GE.0) STOP 'ERROR IN INITYR'
C
C     Read in biological parameters
      READ (IN,'(45X,F8.2)') MSYL
      READ (IN,'((45X,F8.2))') (MSYR(J),J=1,NSTK)
      READ (IN,'((45X,F8.2))') (MORT(I),I=1,3),MAT1,MSIG

C     Read inital recruitment parameters (M & Fem)
      READ (IN,'((45X,F8.2))') RM50,RMSIG,RF50,RFSIG,SELDEC
C     Read SELYR. If >0 the selectivity changes in year SELYR (trial NF-S)
      READ (IN,'((44X,I4))') SELYR,OPTSEL
      IF (SELYR>0.AND.SELDEC>0.0D0) STOP'ERROR: SELYR & SELDEC MISMATCH'
      IF (SELYR>0) SELYR = SELYR + ISCALE
C
C     Read in abundance-estimation-related parameters
      READ (IN,'((45X,F6.3))') BIAS,DF,G0, SSAMPL
      READ (IN,'((45X,10F6.3))') (CVADD(L),L=1,MXSUBA)
C     If SSAMPL > 1 change DF to DF*SSAMPL
      DF = DF * SSAMPL

C     Read in tag loss rates and likelihood weights
      READ (IN,'((45X,F6.3))') TLOSS1,TLOSS2,TAGWT,WTAGEI
      TLOSS1 = EXP(-TLOSS1)
      TLOSS2 = EXP(-TLOSS2)
C     Read years defining linear change between 2 mixing matrices (used if OPTMIX>=2)
      READ (IN,'((44X,I4))') IMIX1,IMIX2,IMIX3
      IMIX1 = IMIX1+ISCALE
      IMIX2 = IMIX2+ISCALE
      IMIX3 = IMIX3+ISCALE
C
C     Read catch data type L (low), B (best) or H (high) series
      READ (IN,'(/A / I2)') CTYPE , NCOL
C     Read the subareas to which the catches (by column) are allocated
      READ (IN,*) (NREG(N),N=1,NCOL)
C     Read added pre 1915 lost rate & which data in SURVEY.DAT to use
      READ (IN,'(//(I2))') CMULT,OPTSGT

C     Hyp. 7 & 8: Read parameters used to split future EG & WI abundances
C                 & check values
      READ (IN,*) EGWIP,EGWIV
      IF (NSUBA/=MXSUBA) THEN
        IF (NSTK/=5) STOP 'ERROR Hyp 7/8 but not 5 stks'
        IF (EGWIP<0.d0 .OR. EGWIP>=1.d0 .OR. EGWIV<0.d0)
     +     STOP 'ERROR in EGWIP or EGWIV'
 
      ELSE
C       Set ANAM2 = names of the MODELLED areas. (ANAM=Management areas.)
        ANAM2 = ANAM
      END IF

C     Read 1st & last years of catch at age data to use in conditioning
      READ (IN,*) ICBYA1,ICBYA2
      ICBYA1 = ICBYA1+ISCALE
      ICBYA2 = ICBYA2+ISCALE

C     Should the second year (+) be excluded
      READ(IN,*) TAGEXC1 

C     Skip 1 spare parameters

      READ (IN,'(//)')

C     Read SRATIO(K): used to set N60N proportion of future abundance
C     =1 except in NF-Q (& Variant V4, for which this value is over-written below)
      READ(IN,*) (SRATIO(K),K=1,MXSUBA)
      READ(IN,*) (VRATIO(K),K=1,MXSUBA)
C
C     Read the numbers of conditioning parameters & their initial values
C     Parameters: NGAMMA gammas (mixing matrix parms),
C             Ndelta dispersal parms,
C             3 tagging parms: lambda, psiCan & psiWI+
C             NSTK*2 pristine stock sizes & initial depletions
C             4 selectivity parameters: RMSIG, RM50, RFSIG,RF50)
C           + 4 extra selectivity parameters if SELYR (as input)>0 (i.e. SELYR + ISCALE<0)
C     Set NOP= no. of parameters used in conditioning
      F1ST=0.1d0
      READ (IN,*) OPTPSI
      READ (IN,*) NGAMMA, F1ST(1:NGAMMA)
      READ (IN,*) NDELTA, F1ST(NGAMMA+1:NGAMMA+NDELTA)
C     Read initial tagging parameters: LAMBDA, psiC, psiWI
      N2 = NGAMMA+NDELTA+3
      READ (IN,*) F1ST(N2-2:N2)
C
      NOP = NGAMMA + NDELTA + 3 + 2*NSTK + 4
      IF (SELYR<0) NOP = NOP + 4
      IF (NOP > MXOP) STOP ' *** ERROR: MXOP IS TOO SMALL'
C
C     Read Carrying capacities
      READ (IN,*) (F1ST(N), N=N2+1,N2+NSTK)
C     Read Initial depletions (not used if OPTEARLY = 0)
      READ (IN,*) (F1ST(N), N=N2+NSTK+1,N2+NSTK*2)
      F1ST(NOP-3) = RMSIG
      F1ST(NOP-2) = RM50
      F1ST(NOP-1) = RFSIG
      F1ST(NOP)   = RF50
      IF (SELYR<0) THEN
        F1ST(NOP-7) = RMSIG
        F1ST(NOP-6) = RM50
        F1ST(NOP-5) = RFSIG
        F1ST(NOP-4) = RF50
      ENDIF
C
C     Read average dispersal rates from stock J to J+1
C          (applies to C substocks only so J=1,NSTK is sufficient)
      READ (IN,*) (OPTDSP(J), J=1,NSTK)
      READ (IN,*) (AVDISP(J), J=1,NSTK)
      NODISP = 0
      IF (NDELTA>0 .OR. SUM(AVDISP(1:NSTK))>0.d0) NODISP = 1

C     Read in mixing matrices. Set USEGAM(I)=T if parameter gamma(i) is used.
C     If JVEGAM(K,J,LP,NS,N)=i then matrix element JVE(K,J,LP,NS,N) is * gamma(i)
      READ (IN,*) NMIX(1:NSTK)
      USEGAM = .FALSE.
      READ (IN,'(///)')
      DO 30 J = 1,NSTK
       IF (OPTMIX>0 .AND. NMIX(J)<1) STOP 'ERROR IN OPTMIX AND NMIX'
       DO 30 N = 1,NMIX(J)
        READ(IN,*) I,I2,(JVE(K,J,N),K=1,NSUBA),(JVEGAM(K,J,N),K=1,NSUBA)
        DO 30 K = 1,NSUBA
          I = JVEGAM(K,J,N)
          IF (I>NGAMMA) STOP 'ERROR: JVEGAM out of range'
          IF (I>0.AND.JVE(K,J,N)==0) STOP'ERROR: GAMMA MULTPLR=0'
          IF (I>0) USEGAM(I) = .TRUE.
   30 CONTINUE
      PRINT *, 'USEGAM:',USEGAM(1:NGAMMA)
C
      CLOSE(IN)
C      
C     Read MANAGE.DAT. If DOCON>0 this is a conditioning run              # INPUT: RUN TYPE
      READ (16,*) DOCON
      IF (DOCON==3) READ (16,*) NSCON
      READ (16,*) OPTDET
C     Read the CLC version (Norwegian or Cooke) to set the catch limit    # INPUT: CLC
      READ (16,*) OPTCLC
      CLOSE(16)
      IF (OPTDET==1) NTRIAL = 0
      IF (DOCON/=0)  NYEAR = 0
C     Open file of conditioning parameter estimates
      IF (DOCON/=1) THEN
        OPEN (23,FILE=nafcon)
        READ (23,'(19X,I5/)') N
        IF (NSTK/=N) STOP 'ERROR IN NAFCON.DAT (No.of stocks)'
      ELSEIF (DOCON == 4)  THEN
C       Real application: was used in 2010 for SAG report
        NTRIAL = 2
        OPEN (40,FILE=ref(1:7)//'.app')
      ENDIF

C     Read the random number generator seeds                              # INPUT: RANDOM NO.
      READ (IN2,'(6I8)') ((ISEEDS(N,I), I=1,6), N=1,400)
      CLOSE(IN2)
C
C     Read tagging data from TAGS.DAT                                     # INPUT: TAGGING
C     Store release operation no.(=NR), no. of releases, area & year
C     RECOBS(NR,IYR,K) No.of observed recoveries in IYR, area K from release NR
C           They were marked in area TAGA(NR), year TAGIYR(NR).
      TAGOP = 0
      RECOBS = 0
      RECTAR = 0.D0
      IGenTag = 0
      KGenTag = 0
      TAGSS = 0
      OPEN (17,FILE="tags.dat")
      READ (17,*) ITAG1
      READ (17,*) NROP
      READ (17,*)
      IF (NAB+NROP>MXFIT) STOP '*** ERROR: NAB+NROP out of range'           ####### change to inc. age data
      READ (17,'(/)')
      ITAG1 = ITAG1 + ISCALE

C     Read: Op,Area,Yr,No Released,Not Recoverd,Same Season Rec,Recoverd
C     NROP= Number of release operations;
C     TSS = # of same season recoveries in release NR (recovered in release area)
      NREL2 = 0.d0
      DO 40 NR=1,NROP
        READ (17,*) N,KT,IT,NREL(NR),NOTREC,TSS,NREC
        IF (N/=NR) STOP 'ERROR 1 READING TAGGING DATA'
        IF (NREL(NR)/=NOTREC+TSS+NREC) STOP 'ERROR IN TAG DATA'
        IT = IT + ISCALE
        TAGIYR(NR) = IT
        K = KT
        IF (NSUBA/=MXSUBA .AND. KT>3) K=KT-1

        IF (TAGOP(IT,K)==0) THEN
C         No other releases so far in year IT, model subarea K
          TAGOP(IT,K)= NR
          TAGA(NR)  = K
          NREL2(NR) = NREL(NR)
          TAGSS(NR)=  TSS
          NR2 = NR
        ELSE
C         This is a combined operation (shown by NREL2(NR) = 0):
C              add these releases to operation NR2
          IF(K==KT)STOP 'TAG ERROR: 1 release only /yr&K except Hyp 7&8'
          NR2 = TAGOP(IT,K)
          IF(NR2==0 .OR. NREL(NR2)<=0) STOP 'TAG DATA ERRROR'
          NREL2(NR2) = NREL(NR2)  + NREL(NR)
          TAGSS(NR2) = TAGSS(NR2) + TSS
        END IF
        TAGREV(NR) =  NR2

C       Read Tag Recovery Data
C       IGenTag & KGenTag(NR,N) = Recovery Year & Area (0 if no recov) for mark N
C       of release NR for use in generating psuedo tag data, so use original subarea
C       1st add same season recoveries (in release area) & then the later recoveries
        DO N2=NOTREC+1,NOTREC+TSS
          KGenTag(NR,N2) = KT
          IGenTag(NR,N2) = IT
        ENDDO
        NTOT = 0
        N2 = NOTREC+TSS
        DO 39 N=1,NREC
          READ (17,*) R2,KM,I
          I=I + ISCALE
          IF (I<=IT) STOP 'ERROR 2 READING TAGGING DATA'
          DO N3 = N2+NTOT+1, N2+NTOT+R2
            KGenTag(NR,N3) = KM
            IGenTag(NR,N3) = I
          ENDDO
C         Add recoveries (other than same season ones) to RECOBS & RECTAR
C         RECOBS is used to generate targets so refers to original subarea KM & release NR
C         RECTAR is compared with model data so uses model subarea K & new release # NR2
          K = KM
          IF (NSUBA/=MXSUBA .AND. KM>3) K=KM-1
          RECOBS(NR,I,KM) = R2
          RECTAR(NR2,I,K) = RECTAR(NR2,I,K) + R2
          NTOT = NTOT+R2
          IF (NTOT==NREC) GO TO 40
   39   CONTINUE
   40 CONTINUE
      CLOSE (17)

      IF (SUM(RECOBS)/=SUM(RECTAR) .OR. SUM(NREL)/=SUM(NREL2))
     +                   STOP 'ERROR 3 IN TAG DATA'

      IF (PRDIAG==1) THEN
        PRINT '(/A/A,F10.1/9G6.0)',' INPUT TAG DATA',
     +     ' TOTAL DIFF SEASON RECOVERIES',SUM(RECTAR),
     +     'NR','Yr','K','K','Rel','SS','TAGREV','TagOp','Rel2'
        DO 44 NR=1,NROP
         I = TAGIYR(NR)
         K = TAGA(NR)
         PRINT '(9G6.0)',NR,I-ISCALE,K,ANAM2(K),NINT(NREL(NR)),
     +          TAGSS(NR),TAGREV(NR),TAGOP(I,K),NINT(NREL2(NR))
         DO 44 N=1,90
          K = KGenTag(NR,N)
          IF (K>0) PRINT '(4G6.0)', 'Rec',IGenTag(NR,N)-ISCALE,K,ANAM(K)
   44   CONTINUE
        PRINT *,' TARGET RECOVERIES Trial 0'
        PRINT '(12X,9A6)','Rel',ANAM2(1:NSUBA)
        DO 45 NR2=1,NROP
         DO 45 I=ITAG1,-1
          IF (SUM(RECTAR(NR2,I,1:NSUBA))>0)PRINT'(12G6.1)',NR2,I-ISCALE,
     +             ANAM2(TAGA(NR2)),NINT(RECTAR(NR2,I,1:NSUBA))
   45   CONTINUE
      END IF
C
C     Read the Variant & the area definitions for use in the RMP      # INPUT: RMP VARIANT
      READ (22,'(I1/1X,I2,T1,2A)') NPROC,ITYPE,VARNAM,TYPEDS
      OPEN(IPNT,FILE=ref(1:7)//'-'//VARNAM(1:2)//'.restest')
      IF (NPROC<1 .OR. NPROC >2) STOP 'ERROR: NPROC not 1 or 2'
C     SCLA determines if the CLA is given different survey data (Variant 4)
      READ (22,*) SCLA
      READ (22,'(///)')

C     Read MAP = area definitions: it is the same for all Variants / procedures
      DO KA=1,MXAR
        READ(22,*) MAP(KA,1:MXSUBA)
      END DO
C     Read the management variant details: Small Area definitions, whether
C     to cap &/or cascade & CSUB = the subareas from which to take the catch
      CSUB = 0
      DO 68 IP=1,NPROC

        CHK  = 0
        READ (22,'(//)')
        DO 67 KA = 1,MXAR
          READ(22,*) KA2,I,NS,CSUB(IP,KA,1:NS)
          IF (KA/=KA2) STOP 'ERROR: READING ATYPE & NSCAT'
          ATYPE(IP,KA) = I
          NSCAT(IP,KA) = NS
          IF ((I==1.OR.I==4.OR.I==9).AND.NS==0) STOP 'ERROR WHEN NS=0'
          IF (I/=1.AND.I/=4.AND.I/=9.AND.NS>0)  STOP 'ERROR WHEN NS>0'
C         Check that a subarea catch is only set once
          DO N=1,NS
            KM=CSUB(IP,KA,N)
            IF (KM==0) STOP 'ERROR IN CSUB'
            IF (CHK(KM)==1) STOP 'ERROR: KM CATCH SET TWICE'
            CHK(KM) = 1
          END DO
   67   END DO
   68 CONTINUE

C     Read MAPC (the same for all variants). Only used when cascading
C          (when ATYPE(IP,KA)==2) or capping (not used in NAF trials)
      READ (22,'(//)')
      DO 70 KA = 1,MXAR
        READ(22,*) (MAPC(KA,KS),KS=1,MXAR)

        IF (ATYPE(IP,KA)==2) THEN
         DO 69 IP = 1,NPROC
          DO 69 K=1,MXAR
           IF (MAPC(KA,K)>0 .AND. ATYPE(IP,K)/=4.AND.ATYPE(IP,K)/=5)
     +                     STOP'ERROR:ATYPE=2 but Small Areas/=4 or 5'
   69    CONTINUE
        END IF
   70 CONTINUE

C     REF holds the Trial no., variant, program version and CLC version
C     Add variant & program version to REF
      REF = REF(1:8) // VARNAM//' '// PROG // OPTCLC
      PRINT '(/1X,4A)', REF, TYPEDS, DESC

C     Initialise & read the catch data                                    # INPUT: CATCHES

C     Read catch by age data, accumulating those in MAXAGE class
C     For hyp. 7 & 8 (NSUBA/=MXSUBA) reorder the catches & sum EG & WI areas
      OCBYAM = 0.d0
      OCBYAF = 0.d0
      READ (24,*)
      DO WHILE (.TRUE.)
       READ (24,*,END=669) I2,IAGE,(CCM(KM),CCF(KM),KM=1,MXSUBA)
       I2 = I2 + ISCALE
       DO 667 KM=1,MXSUBA
C        Set subarea for the current hypothesis
         K = KM
         IF (NSUBA/=MXSUBA .AND. KM>3) K=KM-1
         IF (IAGE > MAXAGE) IAGE = MAXAGE
         OCBYAM(I2,IAGE,K) = CCM(KM) + OCBYAM(I2,IAGE,K)
         OCBYAF(I2,IAGE,K) = CCF(KM) + OCBYAF(I2,IAGE,K)
 667   CONTINUE
      END DO
 669  CLOSE(24)


      IF (PRDIAG==1) THEN
        PRINT '(/A)', ' CATCH AT AGE'
        DO 71 I=ICBYA1,ICBYA2
          DO 71 L=1,MAXAGE
            IF (SUM(OCBYAM(I,L,1:NSUBA))+SUM(OCBYAF(I,L,1:NSUBA))>0)
     +        PRINT'(20I5)',I,L,NINT(OCBYAM(I,L,1:NSUBA)),
     +                          NINT(OCBYAF(I,L,1:NSUBA))
   71   CONTINUE
      END IF

C     Convert to proportions by age and store for conditioning
      DO 72 I=ICBYA1,ICBYA2
        DO 72 K=1,NSUBA
C         Sum over all ages
          OCAM(I,K) = SUM(OCBYAM(I,1:MAXAGE,K))
          OCAF(I,K) = SUM(OCBYAF(I,1:MAXAGE,K))
          IF (OCAM(I,K)>0.d0) THEN 
             OCBYAM(I,1:MAXAGE,K) = OCBYAM(I,1:MAXAGE,K)/OCAM(I,K)
             TOCBYAM(I,1:MAXAGE,K) = OCBYAM(I,1:MAXAGE,K)
          END IF
          IF (OCAF(I,K)>0.d0) THEN
             OCBYAF(I,1:MAXAGE,K) = OCBYAF(I,1:MAXAGE,K)/OCAF(I,K)
             TOCBYAF(I,1:MAXAGE,K) = OCBYAF(I,1:MAXAGE,K)
          END IF 
   72 CONTINUE
C
      CATCH  = 0.d0
      CATCHM = 0.d0
      CATCHF = 0.d0
C     Note: Catches are read by Management area.  They are combined
      OPEN (18,FILE="cat-"//CTYPE//".dat")
      READ (18,'(/)')
      READ (19,*)
      TOTCAT = 0.d0
      DO 53 IYR=INITYR,-1
C
C       1st read Catch by region/operation & add to area defined by NAREA
        READ (18,*) I, (CREG(N),N=1,NCOL)
        IF (IYR/=I+ISCALE) STOP 'ERROR READING CATCH FILE1'
C       If Region NREG(N) = 0 ignore these catches in this run
C       TOTCAT = total historical catch (used to check all taken).
        CCK = 0.d0
        DO 50 N=1,NCOL
          KM = NREG(N)
          IF (KM<0 .OR. KM>MXSUBA) STOP 'ERROR IN NREG'
          IF (KM>0) CCK(KM) = CCK(KM) + CREG(N)
   50   CONTINUE

C       Now read catches known by sex.  Split the rest (=NOSEX) 50:50 males:fem
        READ (19,*) I,(CCM(KM),CCF(KM),KM=1,MXSUBA)
        IF (IYR/=I+ISCALE) STOP 'ERROR READING CATCH FILE2'
C       CMULT = option (used Feb 2015) to increase proportion of lost up to 1915
        IF (I>1915) CMULT = 1.D0

C       KM = management area (7 in all for NAF). K = population model area (6 or 7)
C       For hyp. 7 & 8 (NSUBA/=MXSUBA) reorder the catches & sum EG & WI areas
        DO 51 KM=1,MXSUBA
          K = KM
          IF (NSUBA/=MXSUBA .AND. KM>3) K=KM-1
          IF (CCM(KM)+CCF(KM) > CCK(KM)) STOP 'ERROR IN CATCHES BY SEX'
          NOSEX = CCK(KM)-CCM(KM)-CCF(KM)
          CATCHM(IYR,K) = CATCHM(IYR,K) + CMULT*(CCM(KM) + NOSEX*0.5d0)
          CATCHF(IYR,K) = CATCHF(IYR,K) + CMULT*(CCF(KM) + NOSEX*0.5d0)
          IF (IYR>=FIRSTYR) TOTCAT = TOTCAT + CMULT*CCK(KM)
   51   CONTINUE
   53 CONTINUE
      CLOSE(18)
      CLOSE(19)

      IF (PRDIAG==1) THEN
        PRINT '(/A)', ' CATCHES BY SEX'
        DO IYR=INITYR,-1
          PRINT'(20I5)',IYR-ISCALE,NINT(CATCHM(IYR,1:NSUBA)),
     +                             NINT(CATCHF(IYR,1:NSUBA))
        ENDDO
      END IF
      PRINT *,"cat-"//CTYPE//".dat",'  Total catch:',TOTCAT
C     
C     Read Best catch series for use by CLC (may be the same as read above)
      CATKA  = 0.d0
      OPEN (18,FILE="cat-B.dat")
      READ (18,'(/)')
      DO 57 IYR=INITYR,-1
         READ (18,*) I, (CREG(N),N=1,NCOL)
         IF (IYR/=I+ISCALE) STOP 'ERROR READING CATCH FILE3'
         DO 54 N=1,NCOL
            KM = NREG(N)
            IF (KM>0) CATKA(IYR,KM) = CATKA(IYR,KM)+CREG(N)
 54      CONTINUE
C        Sum catches in assessment areas composed of groups of subareas
         DO 56 KM = 1,MXSUBA
          DO 56 KA = MXSUBA+1,MXAR
            IF (MAP(KA,KM)==1) CATKA(IYR,KA)=CATKA(IYR,KA)+CATKA(IYR,KM)
 56       CONTINUE
C       PRINT '(20I5)',IYR,NINT(CATKA(IYR,1:MXAR))
 57   CONTINUE
      CLOSE(18)
C
C     Read historical abundances for conditioning & CLA (unless variant=V4)   # INPUT: ABUNDANCES
C     CONEST & VAR=observed values used later to reset TAREST for each replicate
C     If OPTSGT=1 use the pro-rated abundance estimates (ABPR) as targets (NF-A3)
      SIGHT  = -1.d0
      CV     = -1.d0
      CV2TAU = -1.d0
      READ(IN4,*) NAB
      READ(IN4,*) I2SYR
      I2SYR = I2SYR+ISCALE
      READ(IN4,*)
      IF (NAB.GT.MXAB) STOP '*** ERROR: NAB out of range'
      IF (PRDIAG==1) PRINT '(/A)', ' ABUNDANCE INPUT: N,K,YR,EST,VAR'
      DO 20 NA = 1,NAB
        READ(IN4,*) IYR,KM,CONEST(NA),CVIN,S,CVCLA,ABPR,TTYPE
        IYR = IYR + ISCALE
        IF (IYR.GE.I2SYR) STOP '*** ERROR IN SURVEY YEAR'
        TARYR(NA) = IYR
        CONVAR(NA) = CVIN * CVIN
        IF (OPTSGT==1) CONEST(NA) = ABPR

C       Add survey data for CLA to SIGHT (is replaced below for Variant 4)
        SIGHT(IYR,KM) = S
        CV(IYR,KM)   = CVCLA
C       Set TAREST = CONEST for the initial (deterministic) run
        TAREST(NA) = CONEST(NA)
        TARVAR(NA) = CONVAR(NA)
C       Store CVs for use in setting TAU
        CV2TAU(IYR,KM)= CVCLA*CVCLA

C       If NSUBA/=MXSUBA (hyp. 7 & 8) reset the model area numbers
        K = KM
        IF (NSUBA/=MXSUBA .AND. KM>3) K=KM-1
        TARA(NA) = K

C       If NSUBA/=MXSUBA (hyp. 7 & 8) combine estimates for EG & WI (KM=3 & 4)
        IF (NSUBA/=MXSUBA .AND. TTYPE==-3) THEN
C         TTYPE=-3: combine this estimate with the previous one (TTYPE=-2)
C         & omit TAREST(NA-1) from the likelihood calculation (shown by TARA=-1)
          IF (KM/=4 .AND. TARA(NA-1)/=3) STOP 'ERROR 1 IN SURVEY DATA'
          TAREST(NA) = CONEST(NA-1)+CONEST(NA)
          TARVAR(NA) = ((CONEST(NA-1)**2)*CONVAR(NA-1)+
     +                  (CONEST(NA)**2)  *CONVAR(NA) ) / TAREST(NA)**2
          TARA(NA-1) = -1
        ENDIF
        IF (PRDIAG==1) PRINT'(3I4,2(F8.0,F8.4),F8.4)',NA,TARA(NA),
     +       TARYR(NA),CONEST(NA),
     +       SQRT(CONVAR(NA)),TAREST(NA),SQRT(TARVAR(NA)),CV2TAU(IYR,K)
 20   CONTINUE

C     Read survey info: RUNSUR = list of areas to be surveyed each year
      READ(IN4,'(////)')
      READ(IN4,*)(IYR,(RUNSUR(I,KA),KA=1,MXAR),IASESS(I),I=I1SYR,IENDYR)
      IF (IYR+ISCALE/=IENDYR) STOP '*** ERROR READING SURVEY YR'
      IF (IASESS(0)/=1) STOP ' ERROR: NO ASSESSMENT IN YEAR 0'
      CLOSE(IN4)
C     Reset RUNSUR if SSAMPL>1 so EG & EI/F are not surveyed in future
      IF (SSAMPL>1.D0) THEN
        RUNSUR(0:NYEAR,3)=0
        RUNSUR(0:NYEAR,5)=0
      END IF

C     If Variant=V4 read replacement historical survey data for use by CLA
C     & replacement values for SRATIO & VRATIO.
      IF (SCLA/='Trial ') THEN
        OPEN (IN4,FILE=SCLA//".DAT")
        READ(IN4,'(/)')
        DO 22 NA = 1,NAB
          READ(IN4,*) IYR,K,S,CVIN
          IYR = IYR + ISCALE
          IF (IYR.GE.I2SYR) STOP '*** ERROR IN SURVEY YEAR'
          SIGHT(IYR,K) = S
          CV(IYR,K) = CVIN
   22   CONTINUE
        READ(IN4,'(/)')
        READ(IN4,*) SRATIO(1:MXSUBA)
        READ(IN4,*) VRATIO(1:MXSUBA)
      ENDIF
      DO K=1,MXSUBA
       IF (SRATIO(K)<=0.0001d0.OR.SRATIO(K)>1.d0) STOP 'ERROR IN SRATIO'
      END DO

C
C     Recheck CPUE.DAT file ##############################                    INPUT: CPUE
C     Read data from CPUE.DAT:
      CPEWT = 0.001
      IF (OPTCPE>0) THEN
        OPEN (20,  FILE='CPUE.dat')
        READ (20,*) NCPUE
        DO 58 I=1,NCPUE
          READ (20,*) I2, Np(I),CPUEK(I)
          IF (I2.NE.I) STOP '*** ERROR READING CPUE SERIES'
          DO 58 II=1,Np(I)
            READ (20,*) CPUEYR(I,II),CPUE(I,II)
            CPUEYR(I,II) = CPUEYR(I,II) + ISCALE
   58   CONTINUE
        READ (20,'(/)')
        DO 59 I=1,4
          READ (20,*) (CPUEV(I,I2),I2=1,4)
   59   CONTINUE
        READ (20,'(/)')
        READ (20,*) (CPUEsig(I),I=5,6)
        CLOSE(20)
        OPEN (96, FILE=ref(1:7)//'.cpe')
      ENDIF

C
c *********************************************************************
C
C     Initialisation & DETERMINISTIC CONDITIONING / RUN                  # SET UP
C
c *********************************************************************
C
C     Call ICLC to initialise management routine
      CALL ICLC(OPTCLC)
C
C     Set up various arrays
      CALL INIT ()
C     Previously called MSYPAR here to set A & Z - now called by HITFUN
C     Set MIXNO which defines the mixing matrices used
      COND1=1
      CALL SETMIX (NYEAR)
C
C     Print list of parameters used, and check values are within allowed range.
      CALL PNTOUT (NTRIAL,NYEAR,IPNT,DESC,TYPEDS)
      PRINT *, 'Initial deterministic setup'

C     Set the initial minimisation step sizes (0 if parameter not used)
      STEP = 0.d0
      F1ST = LOG(F1ST)
C     Set STEP for Initial stock sizes & depletions
      I2 = NGAMMA+NDELTA+3
      STEP(I2+1:I2+NSTK*2) = -0.1d0*F1ST(I2+1:I2+NSTK*2)
      IF (OPTEARLY==0) STEP(I2+NSTK+1:I2+NSTK*2) = 0.d0
C
C     Set STEP2 for next phase: gammas,selectivity dispersion & tag reporting
      STEP2 = STEP
      DO I = 1,NGAMMA
        IF (USEGAM(I)) STEP2(I) = -0.4d0*F1ST(I)
      ENDDO
      STEP2(NGAMMA+1:NGAMMA+NDELTA) =-0.4d0*F1ST(NGAMMA+1:NGAMMA+NDELTA)
      IF (OPTSEL==1) STEP2(I2+NSTK*2+1:NOP)=-0.1d0*F1ST(I2+NSTK*2+1:NOP)

C     Tag reporting rate parameters PSI (NGAMMA+NDELTA+2 & +3 = I2-1 & I2)
C     STEP2(NGAMMA+NDELTA+1) = 0.D0  Set above                           # LAMBDA FIXED
      IF (OPTPSI> 0) STEP2(I2-1) = -0.4d0*F1ST(I2-1)                     # PSI-C  est
      IF (OPTPSI==2) STEP2(I2)   = -0.4d0*F1ST(I2)                       # PSI-WI est

C     Cases when gammas are fixed (i.e. are functions of the other gammas)
C     In these trials only Gamma(1) and Gamma(3) are estimated
      STEP2(2)   = 0.d0
      STEP2(4:6) = 0.d0

      WRITE (42,'(9X,3A3,A13,60(A8,I1.0))') 'N','NS','Ph','Best',
     +     ('G',N,N=1,NGAMMA),('D',N,N=1,2),
     +     'lambda',0,'psiCan',0,'psiWI+',0, ('K',N,N=1,NSTK),
     +     ('Dpln',N,N=1,NSTK),'RMSIG',0,'RM50',0,'RFSIG',0,'RF50'
      WRITE (42,6671) ' ','STEP   ', 0,0,0,0.0,STEP(1:NOP)
      WRITE (42,6671) ' ','STEP2  ', 0,0,0,0.0,STEP2(1:NOP)
      WRITE (42,6671) ' ','F1ST   ', 0,0,0,0.0,EXP(F1ST(1:NOP))
C
      IF (DOCON.EQ.1) THEN
C       Deterministic Conditioning: use results to initialise each
C       stochastic trial. Abund. & propn. targets = best values (set above)
C       NOP fit Parameters = (NGAMMA) gammas (mixing matrix parms),
C             delta (dispersal parms),
C             lambda, psiCan & psiWI+ (tagging parms),
C             (NSTK) pristine stock sizes, (NSTK) initial depletions
C             & 4 (or 8) selectivity parameters
C       F1ST = log of guesses for parameters (read in /set above)

        PRDIAG = 0
        BESTLL = 1.d20
        PHASE = 0
C     
C       Call HITFUN to calculate fit using the initial parameters
        CPENAL = 0.d0
        N=1
   80   CALL HITFUN (F1ST,FIT)
C       The minimum allowed stock depletion is .05
        DO 85 J = 1,NSTK
          PRINT '(A,I5,2F12.4)',' Call 1',J,PMATF(0,J)/KMAT(J),CPENAL
          IF (PMATF(0,J)/KMAT(J) .LT.0.05d0 .OR. CPENAL.GT.1) THEN
            I = NGAMMA+NDELTA+3+J
            F1ST(I) = F1ST(I) + 1.4d0
            STEP(I) = -0.1d0*F1ST(I)
          ENDIF
          N=N+1
  85    CONTINUE
        IF (CPENAL>1 .AND. n<10) GOTO 80
        IF (N>=10) STOP 'ERROR: INIITIAL VALUES NEED ADJUSTING'
C
C       Run non-linear minimisation, phase 1: population sizes & selectivity
C       WtAge downweights the age data to start with
        PHASE = 1
        WtAge = 0.0001d0
        CALL MINIM (F1ST,STEP,NOP,FIT,2000,-1,0.1d0,5,0,1.d0,V0,
     +              HITFUN,6,I2)
        PRINT    6671,REF,'Min01: ',0,0,PHASE,BESTLL
        WRITE(42,6671)' ','Min01: ',0,0,PHASE,BESTLL,EXP(GBEST(1:NOP))
6671    FORMAT(1X,A,A7,3I3,F13.4,60F9.3)

        F1ST = GBEST
        WtAge = WTAGEI
        CALL HITFUN (F1ST,FIT)
        WRITE(*,*) "Best with reweight",FIT
C       
C       Rerun the non-linear minimisation, phase 1: population sizes
        PHASE = 2
        BESTLL = 1.d20
        CALL MINIM (F1ST,STEP,NOP,FIT,2000,-1,0.1d0,5,0,1.d0,V0,
     +              HITFUN,6,I2)

        PRINT    6671,Ref,'Min02: ',0,0,PHASE,BESTLL
        WRITE(42,6671)' ','Min02: ',0,0,PHASE,BESTLL,EXP(GBEST(1:NOP))

C       RUNMIN runs the minimisations. F=input parameters; F1ST=output
        F = GBEST

c        CALL RUNMIN(F,F1ST,STEP2,1,3,0.06d0,0.01d0,0)
        CALL RUNMIN(F,F1ST,STEP2,2,18,0.06d0,0.01d0,0)
        WRITE(42,6671)' ','Min:   ',0,3,PHASE,BESTLL,EXP(GBEST(1:NOP))

      ELSE
C       Read the parameters from NAFCON.DAT
        F = 0.1D0
        READ (23,*) N,OLDFIT,(F(I),I=1,NOP)
        F1ST = LOG(F)
        WtAge = WTAGEI
        IF (DOCON/=2 .AND. DOCON/=3) GO TO 98

C       Deterministic conditioning using initial values from NAFCON.DAT
        PRDIAG= 0
        PHASE = 3
        BESTLL = 1.d20
        CALL HITFUN (F1ST,FIT)
        PRINT '(/A,F18.4)','******  Init fit ' ,FIT
C       Do initial fit with GAMMAS fixed (use STEP)
        CALL MINIM(F1ST,STEP,NOP,FIT,2000,-1,0.05d0,5,0,1.d0,V0,
     +            HITFUN,6,I2)
        PRINT *,'minim',PHASE, BESTLL

C       RUNMIN runs the minimisations. F=input, F1ST=output parameters
        F(1:NOP) = GBEST(1:NOP)
        CALL RUNMIN(F,F1ST,STEP2,2,14,0.04d0,0.001D0,0)
      ENDIF

C     Recall HITFUN to print diagnostics
   98 PHASE = 25
      PRDIAG = 1
      CALL HITFUN (F1ST,FIT)

C     Call PCOND & REPORT to print Deterministic output (inc. best parameters),
C                         PMATF & Catch
      IF (DOCON/=0) CALL PCOND(0,FIT)
      CALL REPORT (IPNT,-1)
      PRINT *, 'End of deterministic setup:',FIT
      IF (DOCON<1) THEN
        PRINT'(/2(A,I3),9F12.4)',' Trial:',0,'/',NTRIAL,OLDFIT,FIT
cc     +                 ,LikeAb,LikeCat,LikeTag,LikeAge
        IF (ABS(OLDFIT-FIT)>0.01d0) STOP ' ERROR IN FIT'
      ENDIF

C     Combine historical surveys as necessary
      PRINT '(/A/A4,12A8)',' CLA HISTORICAL SURVEY DATA','Yr',ANAM
      DO 100 I=I1SYR,-2
        CALL CONSUR (SIGHT,CV,I)
        IF (SUM(SIGHT(I,1:MXSUBA))>0.)
     +    PRINT '(I4,12F8.1/4X,20F8.4)',I,SIGHT(I,1:MXAR),CV(I,1:MXAR)
  100 CONTINUE

      COND1=0

C Remember the estimated sigma parameter form the C@age likelihood
      SIGSQEST = SIGSQ
C
C TRIALS BEGIN ---------------------------------------------------------
C
      DO 500 N = 1,NTRIAL
C **********************************************************************
C
C       START OF TRIAL SET UP - CONDITIONING
C
C **********************************************************************
C
C       If not conditioning, read in the parameters for replicate N
        NS = N
        IF (DOCON<=0 .or. DOCON==3) READ(23,*) NS,OLDFIT,(F(I),I=1,NOP)

C       Reseed random # generators (use NS as some condit.fits may be discarded)
        CALL RESEED (ISEEDS,NS)
        WRITE (IPNT,'(/ A6,I4,A10,6I8/)') 'Trial:',N,
     +         'Seeds:',ISEED1,ISEED2,ISEED3,ISEED4,ISEED5,ISEED6
C
C       Set up mixing (MIXNO array) so rates are fixed for this N
C           (& don't change during conditioning). Random seed 5
c        CALL SETMIX(NYEAR)
C
C       Generate data to use in conditioning or to confirm correct file was read.
C       Both abundance & tag data are generated for all 7 sub-areas and then
C       combined for Hyp 7 & 8, to ensure the targets are independent of hypothesis.

C       Generate abundance targets for management subareas using ISEED4.  {Ref G(a), eqn G.1}
C       Always generate a new (non zero) seed - ensures the same sequence in all trials
C       (-1 used when setting ISEED4 to ensure it is non-zero)
        PRINT '(/A,I3)',' ABUNDANCE TARGETS TRIAL',N
        DO 106 NA = 1,NAB
          ISEED4 = -NINT(RAN1(ISEED4,MA4,INEXT4,INXTP4)*100000.d0)-1
          ISEED  = ISEED4
          SE = SQRT (CONVAR(NA))
  104     R = XNORM(SE,0.d0,ISEED,MA,INEXT,INXTP) - 0.5d0*CONVAR(NA)
          TAREST(NA) = CONEST(NA)*EXP(R)
          IF (TAREST(NA).LT.50.d0) GOTO 104

C         In Hyp. 7 & 8 combine estimates for EG & WI: TARA=-1 for the EG estimate
C                     so it is not used in the likelihood calculation
          IF (TARA(NA)==3) THEN
            IF (TARA(NA-1)==-1) TAREST(NA) = TAREST(NA-1)+TAREST(NA)
          ENDIF
          PRINT '(3I4,2F10.2,F9.5)',NA,TARA(NA),TARYR(NA),CONEST(NA),
     +                                         TAREST(NA),TARVAR(NA)
  106   CONTINUE
C
C       Generate psuedo tag release data: for each release operation (NR)
C       of NREL(NR) releases select the fate of each (with replacement).
C       TAGREV defines the operation to which the recoveries belong.
        RECTAR = 0.d0
        TAGSS  = 0
        DO 107 NR=1,NROP
         NR2 = TAGREV(NR)
         DO 107 N2=1,NINT(NREL(NR))
          N3 = INT(1.0d0+RAN1(ISEED4,MA4,INEXT4,INXTP4) * NREL(NR))
          KM= KGenTag(NR,N3)
          I = IGenTag(NR,N3)
          IF (I.EQ.TAGIYR(NR)) THEN
C           Same season recovery selected (not inc. in likelihood)
            TAGSS(NR2) = TAGSS(NR2) + 1
          ELSE IF (KM>0) THEN
C           Recovery in original area KM year I selected; add to RECTAR
            K = KM
            IF (NSUBA/=MXSUBA .AND. KM>3) K=KM-1
            RECTAR(NR2,I,K) = RECTAR(NR2,I,K) + 1.D0
          ENDIF
C         Otherwise mark is unrecovered, so not inc. in target
  107   CONTINUE
        PRINT '(A,I4)',' TARGET RECOVERIES Trial',N
        PRINT '(12X,12A6)','Rel',ANAM2(1:NSUBA)
        DO 108 NR2=1,NROP
         DO 108 I=ITAG1,-1
          IF(SUM(RECTAR(NR2,I,1:NSUBA))>0) PRINT'(12G6.1)',NR2,I-ISCALE,
     +           ANAM2(TAGA(NR2)),NINT(RECTAR(NR2,I,1:NSUBA))
 108   CONTINUE

C       To come: ##############################################################
C          Generate Catch at age data targets using ISEED6
       ISEED6 = -NINT(RAN1(ISEED6,MA6,INEXT6,INXTP6)*100000.d0)-1
       ISEED  = ISEED6
       DO 109 IYR=ICBYA1,ICBYA2
          DO 109 KM=1,NSUBA
             MCHAMMER = SUM(TOCBYAM(IYR,1:MAXAGE,KM)) + 
     +            SUM(TOCBYAF(IYR,1:MAXAGE,KM))
C     OCAF = SUM(OCBYAF(IYR,1:MAXAGE,K))
             TMPF = 0.D0
             TMPM = 0.D0
             IF(MCHAMMER>0.D0) THEN
                DO 1091 L=1,MAXAGE
                   IF(TOCBYAF(IYR,L,KM) > 0.d0) THEN 
                   SIF = sqrt(SIGSQEST/TOCBYAF(IYR,L,KM))
                   TMPF(L) = EXP(LOG(TOCBYAF(IYR,L,KM)) +  
     +                  XNORM(SIF,0.d0,ISEED,MA4,INEXT4,INXTP))
                   END IF
                   
                   IF(TOCBYAM(IYR,L,KM) > 0.d0) THEN 
                   SIM = sqrt(SIGSQEST/TOCBYAM(IYR,L,KM))
                   TMPM(L) = EXP(LOG(TOCBYAM(IYR,L,KM)) +  
     +                  XNORM(SIM,0.d0,ISEED,MA4,INEXT4,INXTP))
                   END IF
 1091           CONTINUE
                CMRATIO = 1.D0/SUM(TMPM(1:MAXAGE))
                CFRATIO = 1.D0/SUM(TMPF(1:MAXAGE))
                DO 1092 L=1,MAXAGE
                   OCBYAM(IYR,L,KM) = TMPM(L)*CMRATIO
                   OCBYAF(IYR,L,KM) = TMPF(L)*CFRATIO
 1092           CONTINUE
                ENDIF
 109   CONTINUE
C
C       Is this a condition test?
        IF (DOCON>0) THEN
          PRINT *, 'CONDITIONING, TRIAL ',N

C         RUNMIN runs the minimisations. F1ST=input parameters; FBEST=output
          PRDIAG = 0
          IF (DOCON==3) THEN
C           For DOCON==3 skip trials until NS=NSCON
            IF (NS<NSCON) GO TO 500
            F1ST = LOG(F)
          ENDIF
          WRITE(42,6671)' ','F1ST:  ',N,0,0,0.0,EXP(F1ST(1:NOP))

c          CALL RUNMIN(F1ST,FBEST,STEP2,1,3,0.06d0,0.01d0,N)
          CALL RUNMIN(F1ST,FBEST,STEP2,2,20,0.06d0,0.01d0,N)

C         Recall HITFUN to write parameters (so the program can be re-run
C         without re-conditioning) to NAF.ALL  & the diagnostics
  142     PRDIAG = 1
          BESTLL = 1.d20
          PHASE = 25
          CALL HITFUN (FBEST,FIT)
          PRINT *, 'Trial',N,fit
          CALL PCOND (N,FIT)
          CALL REPORT (IPNT,-1)
          PRDIAG = 0
C
C         Conditioning only.  Go on to next replicate
          GO TO 500
C
        ELSE
C
C         Parameters were read above. Check the expected fit is achieved  # CHECK EXPECTED FIT
C         Call HITFUN to run the model forward from INITYR to start of
C           management (year 0), using parameters from conditioning run.
C         (Note: SETMIX was called above to set MIXNO)

          F1ST = LOG(F)
          PRDIAG = 0
          PHASE = 25
C
C         Call HITFUN to setup population at year 0 & calculate fit to data
          CALL HITFUN (F1ST,FIT)
          PRINT'(2(1X,A),I3,A,9F12.4)',REF,'Trial:',N,' FIT:',OLDFIT,FIT
cc     +                 ,LikeAb,LikeCat,LikeTag,LikeAge
          IF (ABS(OLDFIT-FIT) > 0.05d0) STOP ' ERROR IN FIT'
          IF (DOCON==-1) THEN
            CALL PCOND(N,FIT)
            CALL REPORT (IPNT,-1)
            GO TO 500
          ENDIF
        ENDIF

        PRDIAG = 0
C
C
C **********************************************************************
C
C       END OF TRIAL SET UP
C
C **********************************************************************
C
C       RESET resets SIGHT=0 for IYR>0.  Also sets TAU2.
        CALL RESET (SIGHT,CV)

C       Reinitialise management routine
        CALL RCLC ()
C
C       Apply the procedure for NYEARs.
        DO 300 IYR = 0,NYEAR-1
C
C         Call the CLA to set catch limit CATSM for each Small Area & procedure

          DO 200 IP = 1,NPROC
            CALL RUNCLC (CATKA,SIGHT,CV,IYR,INITYR,CATSM,MAPC,
     +                   OPTCLC,IP,NPROC)

C           Set fixed aboriginal catch in WG (subarea 2)
            CATSM(IP,2) = 19.0d0

 200      CONTINUE
          

C         DIVCL allocates the Small Area catch CATSM to subareas (CATCH)
C               If NPROC=2 the value from each procedure is combined
          CALL DIVCL(CATSM,CATCH,NSCAT,CSUB,SIGHT,CV,IYR,NPROC)

C     
C         Allocate catch to sex (50:50) and combine EG+WI in Hyp 7&8
          DO 250 KM=1,MXSUBA
            K = KM
            IF (NSUBA/=MXSUBA .AND. KM>3) K=KM-1
            CATCHM(IYR,K) = CATCHM(IYR,K) + CATCH(K)*0.5D0
            CATCHF(IYR,K) = CATCHF(IYR,K) + CATCH(K)*0.5D0
  250     CONTINUE
          

C         Call SURVEY to set survey estimates SIGHT(IYR-1,K)
C         This allows a 2-year delay between a survey being performed
C         and the results being used to set a catch limit
          CALL SURVEY (SIGHT,CV,IYR-1)
C         Construct estimates for areas larger than subareas
          CALL CONSUR (SIGHT,CV,IYR-1)
C     
C         Project population forward to IYR+1. Set new V matrix using SETV
          CALL STKUPA (IYR)

  300   CONTINUE
C
C       Call REPORT to print PMATF & Catch
        CALL REPORT (IPNT,NYEAR)
C
  500 CONTINUE
C
C     Use PSTATS to print out final statistics
      IF (DOCON.NE.1) CALL PSTATS (IPNT,SIGHT)
      CLOSE (38)
C
      STOP
      END
C
C -----------------------------------------------------------------------
C     CATCH ALLOCATION SUBROUTINE
C       DIVCL  - Allocates the catch limits to subarea. Called by MAIN
C -----------------------------------------------------------------------
C
      SUBROUTINE DIVCL(CATSM,CATCH,NSCAT,CSUB,SIGHT,CV,IYR,NPROC)
C
C     DIVCL takes the catch limits by Small Area set by the RMP & maps
C     them on to the subareas from which the catch will be taken.
C
      USE DECLARECOM
C
      INTEGER NSCAT(2,MXAR),CSUB(2,MXAR,NSUBA),  IYR,NPROC
      REAL(8) CATCH(MXSUBA),CATSM(2,MXAR),CPROC(2,MXSUBA),
     +       SIGHT(I1SYR:IENDYR,MXAR),CV(I1SYR:IENDYR,MXAR)
C
C     Local variables
      INTEGER KA,K,IK,IY,IP
      REAL(8) S2(3),S3(3),GG,Tot,WEIGHT
C
C     CATSM(IP,KA)=Quota for Assessment Area(KA), procedure IP
C     Assessment areas are Small, Combination or Medium areas

C     1st set CPROC(IP,K) = catch in subarea K under procedure IP
C     For Procedure IP
C     If NSCAT(IP,KA)=1 catch is taken from subarea for which CSUB(IP,KA,K)=1
C     If NSCAT>1, the catch from Small Area KA is divided
C                        between the subareas listed in CSUB

C     Clear Catch array
      CPROC = 0.d0
C
      DO 100 IP=1,NPROC
        DO 100 KA = 1,MXAR

         IF (NSCAT(IP,KA)==1) THEN
           K = CSUB(IP,KA,1)
           CPROC(IP,K) = CATSM(IP,KA)

C        Else divide catch using cascading method [not used in NAFin]
         ELSE IF (NSCAT(IP,KA)>1) THEN
          TOT = 0.d0
          S2  = 0.d0
          S3  = 0.d0
          DO 20 IK=1,NSCAT(IP,KA)
           K = CSUB(IP,KA,IK)
           DO 10 IY = I1SYR,IYR-1
             IF (SIGHT(IY,K)>0.0d0) THEN
              GG = 0.9d0**(IYR-IY)/CV(IY,K)**2
              S2(IK) = S2(IK) + LOG(SIGHT(IY,K))*GG
              S3(IK) = S3(IK) + GG
             END IF
   10      CONTINUE
           IF (S2(IK)/=0.0d0) TOT = TOT + EXP(S2(IK)/S3(IK))
   20     CONTINUE
          DO 30 IK=1,NSCAT(IP,KA)
            K = CSUB(IP,KA,IK)
            GG = EXP(S2(IK)/S3(IK))/TOT
            IF(S2(IK)/=0.d0) CPROC(IP,K) = CATSM(IP,KA)*GG
 30    CONTINUE
         ENDIF
100    CONTINUE

C     Set CATCH(K) = catch in subarea
C     If NPROC=2 the Catch is set using variant 1 in years 1-10
C                & then phased in to variant 2 by year 15.
         IF (NPROC==1 .OR. IYR<=9) THEN
            DO 120 K = 1,MXSUBA
               CATCH(K) = CPROC(1,K)
 120        CONTINUE

C     ELSE create weighted sum in yrs 10-14. From yr 15 on use CPROC(2,K)
      ELSE
         WEIGHT = 0.D0
         IF (IYR<15) WEIGHT = 1.D0-(IYR-9.D0)/6.D0
         DO 130 K = 1,MXSUBA
            CATCH(K) = WEIGHT*CPROC(1,K) + (1.D0-WEIGHT)*CPROC(2,K)
 130     CONTINUE
      ENDIF
C
C     CATKA stores catch by assessment area for use by the CLA.
C     1st store catches for the NSUBA areas which are single subareas
      DO 200 K=1,MXSUBA
        CATKA(IYR,K) = CATCH(K)

C       Sum catches in assessment areas composed of groups of subareas
        DO 200 KA = MXSUBA+1,MXAR
          IF (MAP(KA,K).EQ.1) CATKA(IYR,KA) = CATKA(IYR,KA) + CATCH(K)
 200   CONTINUE

      RETURN
      END
C
C     ------------------------------------------------------------------
C     INITIALISATION SUBROUTINES
C       INIT    - Initialise arrays.  Called by MAIN.
C     ------------------------------------------------------------------
C
      SUBROUTINE INIT ()
C
C     Set up tasks
C
      USE DECLARECOM
C
      REAL(8) F
C
C     Call PCOND to write titles
      CALL PCOND(-1,F)
C
      RETURN
      END
C
C     -------------------------------------------------------------------
C     MIXING-RELATED SUBROUTINES
C       SETV    - Sets up the mixing matrix.
C                 Called by MAIN, SETKA
C       SETMIX  - Sets up the mixing options by year
C                 Called by MAIN at the start of each trial
C
C     -------------------------------------------------------------------
C
      SUBROUTINE SETV (IYR)
C
C     SETV sets the V matrix for year IYR
C     This version assumes mixing matrix is independent of age and sex
C
      USE DECLARECOM
C
      INTEGER IYR,J,K,L,N
      REAL(8) TOT,VVAL(MXSUBA)
C
C     At start of each replicate & for deterministic runs (unless OPTMIX>=2)
C       set the mixing matrix V to its expected (normalised) value
      IF (IYR.EQ.INITYR .OR.
     +   (OPTMIX<2 .AND. (OPTDET==1 .OR. COND1.EQ.1)) ) THEN
        DO 10 J = 1,NSTK
        DO 10 K = 1,NSUBA
        DO 10 L  = 1,MAXAGE
          V(K,J,L) = VE(K,J)
  10    CONTINUE

      ELSE
        DO 50 J = 1,NSTK
          TOT = 0.d0

          IF (OPTMIX<2 .OR. NMIX(J).EQ.1) THEN
C           Mixing matrix depends on values in MIXNO (set in SETMIX)
            N = MIXNO(IYR,J)
            DO 21 K = 1,NSUBA
C             Select JV(K,J,N) where N = 1 or 2
              VVAL(K)  = JV(K,J,N)
              TOT = TOT + VVAL(K)
  21        CONTINUE

          ELSE
C         Else OPTMIX=2 or 3 & NMIX=2. Value depends on year. NF-F & NF-G trials
            DO 22 K = 1,NSUBA

              VVAL(K) = JV(K,J,1)
              IF (IYR>IMIX1 .AND. IYR<IMIX2) THEN
C               Linear change from JV(K,J,1) in IMIX1 to JV(K,J,2) in IMIX2
                VVAL(K)= JV(K,J,1) + (JV(K,J,2)-JV(K,J,1))*(IYR-IMIX1)
     +                                                 / (IMIX2-IMIX1)
              ELSE IF (IYR>=IMIX2 .AND. OPTMIX.EQ.2) THEN
C               OPTMIX=2 so matrix remains constant from year IMIX2
                VVAL(K) = JV(K,J,2)
              ELSE IF (IYR>=IMIX2 .AND. IYR<IMIX3) THEN
C               OPTMIX=3. Change from JV(K,J,2) in IMIX2 to JV(K,J,1) in IMIX3
                VVAL(K)= JV(K,J,2) + (JV(K,J,1)-JV(K,J,2))*(IYR-IMIX2)
     +                                                 / (IMIX3-IMIX2)
              ENDIF
              TOT = TOT + VVAL(K)
  22        CONTINUE

          ENDIF

C         Ensure elements of V matrix sum to 1 (=fraction of J stock in area K)
          DO 50 K = 1,NSUBA
            IF (TOT>0.d0) VVAL(K) = VVAL(K) / TOT
            DO 50 L = 1,MAXAGE
              V(K,J,L) = VVAL(K)
  50    CONTINUE
      ENDIF
C
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE SETMIX (NYEAR)
C
C     Called at start of each replicate to set the catch mixing matrices
C     Uses random seed 5. ISEED5 is always -ve on entry
C
      USE DECLARECOM
C
      INTEGER I,NYEAR,INEXT5,INXTP5,J
      REAL(8) RANDOM, RAN1, MA5(55), PMIX(MXSTK)
      EXTERNAL RAN1
C
C     If NMIX(J)=1 there is only 1 mixing matrix for stock J
C     If NMIX(J)=2 for more than 1 stock, they both move in phase (use same random no.)
C     Do for every possible year & stock, so MIXNO is independent of 1st catch yr
C     PMIX = proportion of times matrix 1 is used. Not used here - read in if used ####
      PMIX = 0.5d0
      MIXNO = 1
      IF (OPTDET==1 .OR. COND1>0) RETURN

      DO 50 I = I1YR,NYEAR
       RANDOM = RAN1(ISEED5,MA5,INEXT5,INXTP5)

       DO 50 J = 1,NSTK
         IF (NMIX(J)>1 .AND. RANDOM>PMIX(J)) MIXNO(I,J) = 2
  50  CONTINUE
C
      RETURN
      END
C
C     -------------------------------------------------------------------
C     -------------------------------------------------------------------
C     STOCK-DYNAMICS SUBROUTINES
C       SETKA  - Set up the initial age-structure & project the model
C                forward to the start of year 0.
C                Called by HITFUN
C       INITP  - Called from SETKA to set up the initial age-distribution
C       STKUPA - Update the population vector, allocate catches to stock
C                & calculate the number of births.
C                Called by MAIN, SETKA
C       TAGUPD - Sets TAGM and TAGU, number tagged and recovered
C				 Called by SETKA
C       TAGUPD2- Updates TAGM & TAGU to subsequent year & adjusts for movement
C				 Called by SETKA
C
C     -------------------------------------------------------------------
C
      SUBROUTINE SETKA ()
C
C     Subroutine sets up the pristine population age vectors.
C     Note: survival & birth rates are assumed to be independent of sex.
C
      USE DECLARECOM
C
      REAL(8) P,KTOT(MXSTK),PPRJ
      INTEGER J,L,IYR,K,J2
C
C     Initialise annual population vectors:
      PMATF = 0.D0
      PP    = 0.D0
      KTOT  = 0.D0
C     Initialize the catch at age matrix
      PCBYAM = 0.d0
      PCBYAF = 0.d0
C
C     Call SETV to initialise mixing matrix
      CALL SETV(INITYR)

      DO 501 J = 1,NSTK
C
C      Set mature female carrying capacity
       PMATF(INITYR,J) = KMAT(J)
C
C       Scale the 1st age class so KMAT(J) = # of mature females
C       (=# mature males) in pristine stock J
C       FEC = # female births / # mature females = #births/#mature
       RM(0,J) = KMAT(J) * FEC
       RF(0,J) = RM(0,J)
C
C      Set up pristine nos. by age (male & female)      {Eqn B.1a with C=0}
       DO 20 L = 1,MAXAGE-1           
        RM(L,J) = SUR(L-1) * RM(L-1,J)
        RF(L,J) = RM(L,J)
 20    CONTINUE
       RM(MAXAGE,J) = SUR(MAXAGE-1)*RM(MAXAGE-1,J)/(1.d0-SUR(MAXAGE))
       RF(MAXAGE,J) = RM(MAXAGE,J)
C
C      Set PP(IYR,J,K) = 1+ size of stock J in IYR, subarea K
C      using mean (expectation) mixing matrix.
       DO 50 K = 1,NSUBA
        P = 0.d0
        DO 45 L = 1,MAXAGE
         P = P + RM(L,J)*V(K,J,L)
   45   CONTINUE
        PP(INITYR,J,K) = P+P
        KTOT(J) = KTOT(J) + PP(INITYR,J,K)
   50  CONTINUE
C
C      If Population is not pristine when projection begins
C      set up the initial age-distribution for 1+ population=PPRJ in IYRPRJ
       IF (OPTEARLY.EQ.1) THEN
        NSEL = 1
        PPRJ = INITDEP(J)*KTOT(J)
        CALL INITP(PPRJ,J)
       ENDIF 
C   
  501 CONTINUE   
C
C     Set up the dispersal rates to ensure equilibrium when pristine
      D = 0.d0
      DO 60 J = 1,NSTK
        IF (AVDISP(J)>0.d0) THEN
         J2 = J+1
         D(J,J2) =MIN(0.99D0,AVDISP(J)*(KTOT(J)+KTOT(J2))*0.5d0/KTOT(J))
         D(J2,J) =MIN(0.99D0,D(J,J2) * KTOT(J) /KTOT(J2))
        END IF
   60 CONTINUE
C
      IF (PRDIAG.EQ.1) THEN
        PRINT '(A,10F10.1)',' ktot ',KTOT
        PRINT '(A,10F10.4/ (6X,10F10.4))',' SELM ',SELM(1:MAXAGE,1)
        PRINT '(A,10F10.4/ (6X,10F10.4))',' SELF ',SELF(1:MAXAGE,1)
        PRINT *,'dispersal'
        PRINT '(10F10.7)',(D(1,J2),J2=1,NSTK)
        PRINT '(10F10.7)',(D(2,J2),J2=1,NSTK)
        PRINT '(10F10.7)',(D(3,J2),J2=1,NSTK)
        PRINT '(10F10.7)',(D(4,J2),J2=1,NSTK)
        PRINT '(10F10.7)',(D(5,J2),J2=1,NSTK)
      ENDIF
C
      TAGM = 0.D0
      TAGF = 0.D0
      TAGU = 0.D0
C     Call STKUPA to advance population vectors to start of year 0
C          (and set up the new V matrix, by calling SETV)
C     TAGUPD updates the number of tags during the tagging data period
C     Set NSEL. Reset NSEL if IYR>SELYR to change the selectivity (NF-S trial)
      NSEL = 1 
      DO 200 IYR=FIRSTYR,-1
        IF (IYR.GE.SELYR) NSEL=2
        IF (IYR.GE.ITAG1) CALL TAGUPD(IYR)
        CALL STKUPA (IYR)
        IF (IYR.GE.ITAG1) CALL TAGUPD2()
  200 CONTINUE
C
      RETURN
      END

C     -------------------------------------------------------------------
C
      SUBROUTINE INITP(PPRJ,J)
C
C    Called from SETKA 
C    INITP sets up the initial age-distribution in IYRPRJ by
C    solving the population equations to "hit" PPRJ = relative population size
C
      USE DECLARECOM
C
C     Global variables
      INTEGER J
      REAL(8) PPRJ
C
      REAL(8) P1PL,ROI,FMIN,FMAX,B
      INTEGER II,L
C
      ROI = 0.d0
      CALL FUN2(ROI,P1PL,B,J)
      IF (ABS(P1PL-PPRJ).LT.0.0000001d0*PPRJ) GOTO 1001
C
C     Use a Bisection approach to find the equilibrium age-structure
      FMIN =  0.d0
      FMAX =  1.d0
      DO 1000 II = 1,50
C
C      Bisect F    
       ROI = (FMIN + FMAX) * 0.5d0
       Call FUN2(ROI,P1PL,B,J)
C
C      Check for convergence
       IF (ABS(P1PL-PPRJ).LT.0.0000001d0*PPRJ) GOTO 1001
       IF (P1PL.GT.PPRJ) THEN
         FMIN = ROI
       ELSE
         FMAX = ROI
       ENDIF
C
1000  CONTINUE
1001  CONTINUE
C
C     Rescale the numbers at age to actual values:
C     (FUN2 set RF & UNRF to sizes scaled by B)
      DO 2000 L = 0,MAXAGE
        RM(L,J) =   RM(L,J) * B
        RF(L,J) =   RF(L,J) * B
2000  CONTINUE      
C
      RETURN
      END
C
C     -------------------------------------------------------------------
C
      SUBROUTINE FUN2(ROI,P1PL,B,J)
C
C     FUN2 sets the 1+ population for the given value of ROI (rate of inc.)
C     Called by INITP (only called if IYRPRJ.NE.INITYR)
C
      USE DECLARECOM
C
      INTEGER J
      REAL(8) ROI,P1PL
C
C     Local variables
      REAL(8) PMAT,B,B2
      INTEGER L
      
C
C     Calculate the numbers at age relative to the number of 0 yr olds {C.4}
C     NB Population is in equilibrium so males=females as SUR is sex independent
      RF(0,J) = 1.0d0
      RM(0,J) = 1.0d0
      DO 70 L = 0,MAXAGE-1
       RF(L+1,J) = SUR(L)*RF(L,J)*(1.d0-ROI*SELF(L,NSEL))
       RM(L+1,J) = SUR(L)*RM(L,J)*(1.d0-ROI*SELM(L,NSEL))
 70   CONTINUE
C   
C     Adjust for last age class being pooled (always fully recruited so SELF(MAXAGE)=1 )
      RF(MAXAGE,J) = RF(MAXAGE,J) / (1.d0-SUR(MAXAGE)*(1.d0-ROI))
      RM(MAXAGE,J) = RM(MAXAGE,J) / (1.d0-SUR(MAXAGE)*(1.d0-ROI))

C     Add 1+ & mature female totals
      P1PL = 0.d0
      PMAT = 0.d0
      DO 90 L = 1,MAXAGE
       P1PL = P1PL + (RM(L,J) + RF(L,J))
       PMAT = PMAT + RF(L,J)*FMATUR(L)
 90   CONTINUE
C
C     Find no. of births B {C.5} and hence 1+ population (P1PL)
      B = 0.d0
      B2 = 1.d0 - (1.d0/(FEC*PMAT)-1.d0)/A(J)
      IF (B2.GT. 0.d0) THEN                                      
       B2 = KMAT(J) * B2**(1.d0/Z(J))
       B = B2 / PMAT
      ENDIF
      P1PL = P1PL*B
C
      RETURN
      END
C
C     -------------------------------------------------------------------
C
      SUBROUTINE STKUPA (IYR)
C
C     STKUPA updates the stock sizes at the start of IYR+1, i.e. it
C     updates the RM & F arrays and sets PMATF & F(IYR+1,J).
C
      USE DECLARECOM
C
      REAL(8) FM(MXSUBA),FF(MXSUBA),SUMVRM,SUMVRF,UCRM,UCRF,
     +       FMAT, CCM,CCF,SM(MXSTK),SF(MXSTK),
     +       RMT(0:MXAGE2,MXSTK),RFT(0:MXAGE2,MXSTK),VSFM,VSFF
      INTEGER IYR,L,J,K,I1,J2
C
C     STKUPA sets population at start of year I1
      I1 = IYR + 1
C
C     Zero catch array (PMATF & PP were set to 0 in SETKA)
      CATT(IYR,1:NSTK) = 0.d0
C
C     First calculate fishing mortality by area
      DO 50 K = 1,NSUBA
         
       FM(K) = 0.d0
       FF(K) = 0.d0
       CCM = CATCHM(IYR,K)
       CCF = CATCHF(IYR,K)
       IF (CCM+CCF==0.D0) GO TO 50

C      Set SUMVR = no. of whales available to catching in subarea K  {Eqn D.1}
C      Mixing matrix V(K,J,L) = fraction of age L in stock J that are in K
       SUMVRM = 0.d0
       SUMVRF = 0.d0
       DO 45 J=1,NSTK
        SM(J) = 0.d0
        SF(J) = 0.d0
        DO 40 L=1,MAXAGE
         SM(J) = SM(J) + V(K,J,L)*RM(L,J)*SELM(L,NSEL)
         SF(J) = SF(J) + V(K,J,L)*RF(L,J)*SELF(L,NSEL)
 40     CONTINUE
        SUMVRM = SUMVRM + SM(J)
        SUMVRF = SUMVRF + SF(J)
 45    CONTINUE
C
C      Set the fishing mortality in subarea K         {Eqn D.2}
       IF (SUMVRM.GT.0.d0) FM(K) = MIN(1.d0, CCM / SUMVRM)
       IF (SUMVRF.GT.0.d0) FF(K) = MIN(1.d0, CCF / SUMVRF)
C     
C      Save selected numbers (relevant to CPUE)       {Eqn G.10}
       SELN(IYR,K) = SUMVRM+SUMVRF
C     
C      Set the total catch from stock (assuming stock is large enough
C      to support catch)
       DO 46 J=1,NSTK
        CATT(IYR,J) = CATT(IYR,J) + SM(J)*FM(K) + SF(J)*FF(K)
 46    CONTINUE
C
 50   CONTINUE
C
C     Now update dynamics (computes the post-survival numbers and hence total abundance)
      DO 100 J = 1,NSTK
       DO 52 L = 1,MAXAGE
C       Set UCRM = mean survivorship of stock J of age L after catching {EqnD.1}
C       & accumulate the catch @ age matrix PCBYAM (sum over stocks)
        UCRM = 1.d0
        UCRF = 1.d0
        DO 51 K = 1,NSUBA
          VSFM = V(K,J,L) * SELM(L,NSEL) * FM(K)
          VSFF = V(K,J,L) * SELF(L,NSEL) * FF(K)
          UCRM = UCRM - VSFM
          UCRF = UCRF - VSFF
          PCBYAM(IYR,L,K) = PCBYAM(IYR,L,K) + RM(L,J)*VSFM
          PCBYAF(IYR,L,K) = PCBYAF(IYR,L,K) + RF(L,J)*VSFF
 51     CONTINUE
C     
C       Remove catches & natural deaths from all age classes      {Eqn B.1}
        RM(L,J) = MAX(0.d0, RM(L,J) * UCRM * SUR(L))
        RF(L,J) = MAX(0.d0, RF(L,J) * UCRF * SUR(L))
 52    CONTINUE
       RM(0,J) = RM(0,J) * SUR(0)
       RF(0,J) = RF(0,J) * SUR(0)
         
C      Pool last 2 age classes (prior to advancing to next year)
       RM(MAXAGE-1,J) = RM(MAXAGE,J) + RM(MAXAGE-1,J)
       RF(MAXAGE-1,J) = RF(MAXAGE,J) + RF(MAXAGE-1,J)
         
C      Update dynamics: Advance age-classes to next year
       DO 100 L=MAXAGE,1,-1
        RM(L,J) = RM(L-1,J)
        RF(L,J) = RF(L-1,J)
        RMT(L,J)= RM(L,J)
        RFT(L,J)= RF(L,J)
  100 CONTINUE
C
C     Update dynamics: Add in dispersal (Density independent dispersal between stocks)
      IF (NODISP>0) THEN
        DO 1022 J=1,NSTK
         DO 1022 J2=1,NSTK
          IF (J2/=J .AND. (D(J,J2)+D(J2,J))>0.D0) THEN
           DO L=1,MAXAGE
             RM(L,J) = RM(L,J) - RMT(L,J)*D(J,J2) + RMT(L,J2)*D(J2,J)
             RF(L,J) = RF(L,J) - RFT(L,J)*D(J,J2) + RFT(L,J2)*D(J2,J)
           ENDDO
          ENDIF
1022    CONTINUE
      ENDIF

      DO 220 J = 1,NSTK

C       Accumulate number mature                         {Eqn C.2}
        DO 210 L=1,MAXAGE
          PMATF(I1,J) = PMATF(I1,J) + RF(L,J)*FMATUR(L)
  210   CONTINUE

C       Put new births into age-class 0                  {Eqn C.1}
C       -ve Births occur in conditioning if dispersed stocks have v.different Ks
        FMAT = PMATF(I1,J)
        RM(0,J) = FEC*FMAT*(1.d0+A(J)*(1.d0-(FMAT/KMAT(J))**Z(J)))
        RM(0,J) = MAX(0.d0, RM(0,J))
        RF(0,J) = RM(0,J)
  220 CONTINUE

C     Set up the mixing matrix for year I1.
      CALL SETV (I1)

C     Accumulate 1+ of stock J in subarea K = PP(I1,J,K)
      DO 250 J = 1,NSTK
       DO 250 L=1,MAXAGE
        DO 250 K = 1,NSUBA
          PP(I1,J,K) = PP(I1,J,K) + (RM(L,J)+RF(L,J))*V(K,J,L)
C          if (rm(L,J)<-.01 .or. rf(L,J)<-.01) STOP 'STKUPA error: -ve P'
  250 CONTINUE
C
      RETURN
      END
C
C     -------------------------------------------------------------------
C
      SUBROUTINE TAGUPD (IYR)
C
C     Update TAGM(NR,L,J)/TAGF =no. of stock J, age L m/f tagged during release NR
C     Set TAGU(NR,IYR,K) = no. recovered in area K in IYR for release NR
C
      USE DECLARECOM
C
      REAL(8) FM(MXSUBA),FF(MXSUBA),TM(MXSUBA),TF(MXSUBA),SUMVRM,SUMVRF
      REAL(8) CCM,CCF,UCRM,UCRF,SUMM,SUMF,SL,TNEW
      INTEGER IYR,L,J,K,NR2
C
C     TAGUPD is called before STKUPA.
C     1st set catches in the current tagged population (if any) & calculate recoveries
C       (before adding New tags to the tagged population)
      DO 50 K = 1,NSUBA
        CCM = CATCHM(IYR,K)
        CCF = CATCHF(IYR,K)

C       SUMVR = whales available to catching in subarea K  {bottom of Eqn D.2}
C       SUM   = whales available for tagging in subarea K  {end of Eqn G.5}
C       Mixing matrix V(K,J,L) = fraction of age L in stock J that are in K
        SUMVRM = 0.d0
        SUMVRF = 0.d0
        SUMM = 0.d0
        SUMF = 0.d0
        DO 45 J=1,NSTK
          DO 45 L=1,MAXAGE
            SUMVRM = SUMVRM + V(K,J,L)*RM(L,J)*SELM(L,NSEL)
            SUMVRF = SUMVRF + V(K,J,L)*RF(L,J)*SELF(L,NSEL)
            SUMM   = SUMM   + V(K,J,L)*RM(L,J)
            SUMF   = SUMF   + V(K,J,L)*RF(L,J)
   45   CONTINUE
C
C       Set the fishing mortality in subarea K         {Eqn D.2}
        FM(K) = 0.d0
        FF(K) = 0.d0
        IF (SUMVRM.GT.0.d0) FM(K) = MIN(1.d0, CCM / SUMVRM)
        IF (SUMVRF.GT.0.d0) FF(K) = MIN(1.d0, CCF / SUMVRF)
C
C       Set TAGU(NR,IYR,K) = whales reported recovered this year in area K
C       that were marked in release NR2 (ignore same season recoveries)  {Eqn G.6}
C       TAGU was initalised to 0 in SETKA
        IF (CCM+CCF.GT.0.001D0) THEN
          DO 48 J=1,NSTK
           DO 48 L=1,MAXAGE
            DO 48 NR2 = 1,NROP
             TAGU(NR2,IYR,K) = TAGU(NR2,IYR,K) +
     +           PSI(K)*V(K,J,L)*(TAGM(NR2,L,J)*FM(K)*SELM(L,NSEL) +
     +                            TAGF(NR2,L,J)*FF(K)*SELF(L,NSEL) )
   48     CONTINUE
C       ELSE No catch in K this year
        ENDIF

C       Set TM=proportion of males tagged this year in K  {part of Eqn G.5}
c       If no catch taken in subarea then assume sex ratio 50:50 instead
C       Release operation = NR2 as given by TAGOP
        NR2 = TAGOP(IYR,K)
        TM(K) = 0.0d0
        TF(K) = 0.0d0
        IF (NR2>0 .AND. SUMM+SUMF> 0.D0)THEN
          TNEW = NREL2(NR2) - DBLE(TAGSS(NR2))/PSI(K)
          TNEW = MAX(0.D0,TNEW)
          IF (CCM<1.D0 .AND. CCF<1.D0) THEN
           TM(K) = TNEW * 0.5D0 / SUMM
           TF(K) = TNEW * 0.5D0 / SUMF
          ELSE
           TM(K) = TNEW * CCM / ((CCM+CCF) * SUMM)
           TF(K) = TNEW * CCF / ((CCM+CCF) * SUMF)
          ENDIF
        ENDIF

   50 CONTINUE
C
C
C     Update TAGM & TAGF(L,J,K): remove deaths & add new tags. K=release area
C     TAGM(NR,L,J) =no. of stock J males, age L tagged in release NR

      DO 100 J = 1,NSTK

        DO 56 L = 1,MAXAGE
C         Set UCRM = mean survivorship of stock J of age L after catching {EqnG.5a}
          UCRM = 1.d0
          UCRF = 1.d0
          DO 51 K = 1,NSUBA
            UCRM = UCRM - V(K,J,L) * SELM(L,NSEL) * FM(K)
            UCRF = UCRF - V(K,J,L) * SELF(L,NSEL) * FF(K)
  51      CONTINUE

C         Update TAGM & TAGF: remove catches & add new tags  {Eqn G.5a}
          DO 52 NR2 = 1,NROP
C           First remove catches and natural deaths from old tagged whales
            TAGM(NR2,L,J) = MAX(0.D0, TAGM(NR2,L,J) *UCRM*SUR(L)*TLOSS2)
            TAGF(NR2,L,J) = MAX(0.D0, TAGF(NR2,L,J) *UCRF*SUR(L)*TLOSS2)
  52      CONTINUE

C         Add new tags this released this year   {Eqn G.5a}
          SL =  SQRT(SUR(L)*TLOSS1)
          DO 54 K=1,NSUBA
            NR2 = TAGOP(IYR,K)
C           If NR2=0 there was no tagging in IYR, area K; skip to next area
            IF (NR2.EQ.0) GO TO 54
            TAGM(NR2,L,J) = TAGM(NR2,L,J) + TM(K)*V(K,J,L)*RM(L,J)*SL
            TAGF(NR2,L,J) = TAGF(NR2,L,J) + TF(K)*V(K,J,L)*RF(L,J)*SL
  54      CONTINUE
  56    CONTINUE

C       Pool last 2 age classes (prior to advancing to next year)
        DO 60 NR2 = 1,NROP
          TAGM(NR2,MAXAGE-1,J) = TAGM(NR2,MAXAGE,J)+TAGM(NR2,MAXAGE-1,J)
          TAGF(NR2,MAXAGE-1,J) = TAGF(NR2,MAXAGE,J)+TAGF(NR2,MAXAGE-1,J)
  60    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C
C     -------------------------------------------------------------------
C
      SUBROUTINE TAGUPD2 ()
C
C     Only for movement of tags (has to come after STUKPA {Eqn G.5b)
C     i.e. apply dispersal to the tagged population

      USE DECLARECOM
C
C     Local variables
      INTEGER L,J,NR,J2      
      REAL(8) RMT(0:MXAGE2,MXSTK),RFT(0:MXAGE2,MXSTK)
C
C     Update dynamics: Advance age-classes to next year
C     (TAGM always =0 for L=0 age class (set in SETKA & never changed)
      DO 220 NR = 1,NROP

       DO 100 J = 1,NSTK
        DO 100 L=MAXAGE,1,-1
         TAGM(NR,L,J) = TAGM(NR,L-1,J)
         TAGF(NR,L,J) = TAGF(NR,L-1,J)
         RMT(L,J)  = TAGM(NR,L,J)
         RFT(L,J)  = TAGF(NR,L,J)
  100  CONTINUE
C      
C      Update dynamics: Add in dispersal (not stock J to stock J)
       IF (NODISP>0) THEN
        DO 211 J=1,NSTK
         DO 211 J2=1,NSTK
          IF (J2/=J .AND. (D(J,J2)+D(J2,J))>0.D0) THEN
           DO L=1,MAXAGE
            TAGM(NR,L,J)=TAGM(NR,L,J)-RMT(L,J)*D(J,J2)+RMT(L,J2)*D(J2,J)
            TAGF(NR,L,J)=TAGF(NR,L,J)-RFT(L,J)*D(J,J2)+RFT(L,J2)*D(J2,J)
           ENDDO
          ENDIF
  211   CONTINUE
       ENDIF
  220 CONTINUE
C
      RETURN
      END
C
C     -------------------------------------------------------------------
C     CONDITIONING SUBROUTINES
C      RUNMIN - controls the calling of the minimisation routines
C      HITFUN - Calculates the deviance as a function of the model parameters
C             - & performs most of the housekeeping.  Called by MAIN
C      DEVIN  - Calculates the deviance               Called by HITFUN
C      HITFUN2- called by the HOOKE minimiser to set FF & call HITFUN
C      FUNK   - called by the FIT22 minimiser to set FF & call HITFUN
C
C ----------------------------------------------------------------------
C
      SUBROUTINE RUNMIN(FINP,FOUT,STEP1,MAXS,MAXP,STOP1,STOP2,NT)
C
C     RUNMIN controls the calling of the minimisation routines
C
      USE DECLARECOM
C
      REAL(8) FINP(MXOP),FOUT(MXOP), STEP1(MXOP),STEP(MXOP),STOP1,STOPC,
     +       STOP2,BLL,XVALS(100),SSBEST,F(MXOP),V0(MXOP),FIT,F2(100)
      REAL(8) XNORM,RAN1,HITFUN2,FUNK
      INTEGER MAXS,MAXP,NS,NT,IPNT,I,IFAIL,HOOKE,NN

      EXTERNAL HITFUN,XNORM,RAN1,HITFUN2,HOOKE,FUNK

    1 FORMAT (1x,2A,3I3,F13.4,60F9.3)
      BESTLL = 1.d20
      GBEST = FINP
      STEP = STEP1

C     Now estimate ALL parameters
      DO 99 NS=MAXS,1,-1

C       Call the non-linear minimisation routine
        STOPC = STOP1
        DO 95 PHASE=3,MAXP
          SSBEST = BESTLL
          BESTLL = 1.d20
C
C         Hooke minimisation, all parameters
          INITV = GBEST
          IPNT = 0
          ISEST = 0
          DO 86 I = 1,NOP
           IF (STEP(I)/=0.d0) THEN
             ISEST(I) = 1
             IPNT = IPNT + 1
             XVALS(IPNT) = GBEST(I)
            ENDIF
 86       CONTINUE
          NN = hooke (IPNT,XVALS,F2,0.5d0,STOPC,10000,HITFUN2)
          PRINT 1,Ref,'Hooke: ',NT,NS,PHASE,BESTLL

C         MINIM Minimisation, all parameters
          F = GBEST
          CALL MINIM (F,STEP,NOP,FIT,2000,-1,STOPC,5,0,1.d0,V0,
     +               HITFUN,6,IFAIL)
          PRINT 1,Ref,'MINIM: ',NT,NS,PHASE,BESTLL

C         BFGS minimisation, all parameters
          INITV = GBEST
          IPNT = 0
          ISEST = 0
          DO 88 I = 1,NOP
           IF (STEP(I)/=0.d0) THEN
             ISEST(I) = 1
             IPNT = IPNT + 1
             XVALS(IPNT) = GBEST(I)
            ENDIF
 88       CONTINUE
          CALL BFGS(IPNT,XVALS)
          PRINT 1,Ref,'BFGS:  ',NT,NS,PHASE,BESTLL

C         FIT22 Minimisation, Gammas & dispersion parameters only
          INITV = GBEST
          IPNT = 0
          ISEST = 0
          DO 90 I = 1,NGAMMA+NDELTA
           IF (STEP(I)/=0.d0) THEN
             ISEST(I) = 1
             IPNT = IPNT + 1
             XVALS(IPNT) = GBEST(I)
            ENDIF
 90       CONTINUE
          CALL FIT22(XVALS,BLL,IPNT)
          PRINT 1,Ref,'FIT22a:',NT,NS,PHASE,BESTLL

C         MINIM Minimisation, all parameters
          F = GBEST
          CALL MINIM (F,STEP,NOP,FIT,2000,-1,STOPC,5,0,1.d0,V0,
     +               HITFUN,6,IFAIL)
          PRINT 1,Ref,'MINIM: ',NT,NS,PHASE,BESTLL,
     +                 LikeAb,LikeCat,LikeTag,LikeAge
          WRITE (42,1)' ','RUNMIN:',NT,NS,PHASE,BESTLL,
     +          EXP(GBEST(1:NOP)),LikeAb,LikeCat,LikeTag,LikeAge
C
          IF (PHASE<6) STOPC = STOPC * .05d0

          IF (ABS(BESTLL-SSBEST)<STOP2 .AND.PHASE>8) GOTO  98
   95   CONTINUE

C       Reset the minimisation step sizes
   98   STEP = STEP1 * DBLE(NS)/10.d0
   99 CONTINUE
C
      FOUT = GBEST
      RETURN
      END

C -----------------------------------------------------------------------
C
      REAL(8) FUNCTION HITFUN2(F,N)
C
C     HITFUN2 is called by the HOOKE minimiser to set up FF & call HITFUN
C
      USE DECLARECOM
C
      REAL(8) F(*),FOUT, FF(MXOP)
      INTEGER II,IPNT2,N
C
C     Move through all parameters
      IPNT2 = 0
      DO 1000 II = 1,MXOP
       IF (ISEST(II).EQ.1) THEN
        IPNT2 = IPNT2 + 1
        FF(II) = F(IPNT2)
       ELSE
        FF(II) = INITV(II)
       ENDIF
1000  CONTINUE
C      
      CALL HITFUN(FF,FOUT)
      HITFUN2 = FOUT
      
      RETURN
      END

      REAL(8) FUNCTION FUNK(F)
C
      USE DECLARECOM
C
      REAL(8) F(*),FOUT,FF(MXOP)
      INTEGER II,IPNT2
C
C     Move through all parameters
      IPNT2 = 0
      DO 1000 II = 1,MXOP
       IF (ISEST(II).EQ.1) THEN
        IPNT2 = IPNT2 + 1
        FF(II) = F(IPNT2)
       ELSE
        FF(II) = INITV(II)
       ENDIF
1000  CONTINUE
C
      CALL HITFUN(FF,FOUT)
      FUNK = FOUT

      RETURN
      END

C ----------------------------------------------------------------------

      SUBROUTINE HITFUN (F,FUNCT)
C
C     Function calculates the likelihood to be minimised.
C
      USE DECLARECOM
C
      INTEGER I,J,K,N,L
      REAL(8) F(*),FUNCT, DEVIN,TOT,GG(0:MXOP),XAGE,DPLMEAN,RM50,RF50
      EXTERNAL devin
C
C     Reset FUNCT
      FUNCT = 0.d0
C
C     F array contains log values for the gamma parameters + pristine stock sizes
C     On first call in conditioning runs F = log of guessed values (read in)
C     First anti-log the parameters. 
      DO 10 I = 1,NOP
       IF (F(I).GT.100) THEN
         FUNCT = FUNCT + 100000.d0 + F(I)
       ELSE
         GG(I) = EXP(F(I))
         IF (GG(I)>200000.d0.OR.GG(I)<0.999d-6) FUNCT=FUNCT+100000.d0
       ENDIF
  10  CONTINUE

C     Gammas 1 & 3 are proportions so constrain to be < 1
      IF (GG(1) > 1.D0) FUNCT = FUNCT + 100000.d0 + GG(1)
      IF (GG(3) > 1.D0) FUNCT = FUNCT + 100000.d0 + GG(3)
C      if(ngamma .gt. 6) then
C         IF (GG(8) > 1.D0) FUNCT = FUNCT + 100000.d0 + GG(8)
C      endif
C     Cases when gammas are fixed (i.e. are functions of the other gammas)
      GG(4) = 1.D0 - GG(1)
      IF (NODISP==0) THEN
C       Hypotheses 4 & 8: NODISP = 0 (no dispersal). Ensure GG(5) is +ve
        GG(2) = GG(1) * GG(3)
        IF (GG(3) > 0.4999D0) FUNCT = FUNCT + 100000.d0 + GG(3)
        GG(5) = 1.D0 - 2*GG(3)
        GG(6) = GG(4) * GG(3)
      ELSE 
CIF (NSTK==5) THEN 
C       Hypotheses 6 & 7
        GG(2) = 1.D0 - GG(3)
        GG(5) = GG(1) * GG(2)
        GG(6) = GG(4) * GG(2)
C        GG(7) = 1.D0 - GG(8)
      ENDIF
      IF (FUNCT>=100000.D0 ) RETURN
C
C     Set up the JV matrix = JVE matrix * GAMMA parameters
C     JVEGAM gives index of gamma parameter, or else is 0
C     ie if JVEGAM(K,J,N)=i then JV(k,j,n) = JVE(K,J,N) * gamma(i).
C     Set GG(0)=1 so JV values are fixed @ JVE(K,J,N) if JVEGAM(K,J,N)=0
      GG(0) = 1.d0
      DO 200 J = 1,NSTK
      DO 200 N = 1,NMIX(J)
        TOT = 0.d0
        DO 150 K = 1,NSUBA
          JV(K,J,N)  = JVE(K,J,N) * GG(JVEGAM(K,J,N))
          TOT = TOT + JV(K,J,N)
 150    CONTINUE
C
C       Ensure elements of JV matrix sum to 1 (=fraction of J stock in area K)
        DO 180 K = 1,NSUBA
          IF (TOT.NE.0.d0) JV(K,J,N) = JV(K,J,N) / TOT
 180    CONTINUE
 200  CONTINUE
C
C     Set up the expected mixing matrix (VE) =mean of JV matrices
      DO 500 J = 1,NSTK
        DO 250 K = 1,NSUBA
C         1 mixing matrix (=norm for NAFin trials) or OPTMIX>=2 (when VE=initial matrix)
          VE(K,J) = JV(K,J,1)
C         If NMIX=2 take the average of 2 mixing matrices ###### not used here
C          IF (NMIX(J).EQ.2.AND.OPTMIX<2) VE(K,J) =(VE(K,J)+JV(K,J,2))*.5d0
 250    CONTINUE
C
C       Set pristine mature stock sizes (stored in GG after gamma & other parameters)
        I = NGAMMA+NDELTA+3+J
        KMAT(J)   = GG(I)
        INITDEP(J)= GG(I+NSTK)
        IF (KMAT(J)<300.d0) FUNCT = FUNCT+100000.d0+(300-KMAT(J))*100.d0
        IF (INITDEP(J)>1.d0)FUNCT = FUNCT+100000.d0+ INITDEP(J)*100.d0
  500 CONTINUE
C     Penalise significant deviations in initial depletion
C     Note: if this is used in future it must be inc. in specifications #########
C     NOTE: is adding to likelihood in all cases - removed as not using INITDEP now
c      DPLMEAN = SUM(INITDEP(1:NSTK)) / DBLE(NSTK)
c      DO 502 J = 1,NSTK
c        FUNCT = FUNCT + 100*(INITDEP(J) - DPLMEAN)**2
c 502  CONTINUE
C
C     Set the average dispersals using GG(NGAMMA+1 & +2)
      I = NGAMMA
      DO 600 J = 1,NSTK
       IF (OPTDSP(J)==1) THEN
         I = I + 1
         AVDISP(J) = GG(I)
         IF (AVDISP(J)>0.4d0) FUNCT = FUNCT + 100000.0d0
       ENDIF
 600  CONTINUE
C
C     Set reporting rate parameters PSI (MAX VALUE=1).
C     PSI(K) Tag reporting rate parameter for subarea K
C     OPTPSI =0 if all PSI fixed (Old trial NF27);
C            =1 to estimate PSI in EC (& WG??) (the Baseline)
C            =2 to estimate PSI in EC (& WG??) & PSI for other areas (not used now)
C  we are using the same PSI for EC & WG but specs just say C - which should be used
C  I think this arose because in the  1st set of trials EC/WG were combined
C  ####### check - specifications seem to say PSI(WG) = PSI(other)
      PSI(1:2)    = GG(NGAMMA+NDELTA+2)
      PSI(3:NSUBA)= GG(NGAMMA+NDELTA+3)
      IF (PSI(1)>1.D0) FUNCT = FUNCT + 100000.d0 + PSI(1)*10.d0
      IF (PSI(3)>1.D0) FUNCT = FUNCT + 100000.d0 + PSI(3)*10.d0
C
C     Set Selectivity SETO (SELM,RMSIG,RM50,MAXAGE).
C     Constrain RM50=GG(NOP-2) to be between 1 & MAXAGE
      XAGE = DBLE(MAXAGE)
      RM50 = GG(NOP-2)
      IF (RM50>XAGE) FUNCT = FUNCT + 100000.d0 * (1.D0+RM50-XAGE)
      IF (RM50<1.d0) FUNCT = FUNCT + 100000.d0 * (2.d0-RM50)
      RF50 = GG(NOP)
      IF (RF50>XAGE) FUNCT = FUNCT + 100000.d0 * (1.d0+RF50-XAGE)
      IF (RF50<1.d0) FUNCT = FUNCT + 100000.d0 * (2.d0-RF50)

      IF (FUNCT>=100000.D0 ) RETURN

      CALL SETO (SELM(:,1),GG(NOP-3),RM50,MAXAGE)
      CALL SETO (SELF(:,1),GG(NOP-1),RF50,MAXAGE)
C
C     NF-S trials: set the Selectivity for use from SELYR on (NSEL=2)
      IF (SELYR<0) THEN
        RM50 = GG(NOP-6)
        IF (RM50>XAGE) FUNCT = FUNCT + 100000.d0 * (1.D0+RM50-XAGE)
        IF (RM50<1.d0) FUNCT = FUNCT + 100000.d0 * (2.d0-RM50)
        RF50 = GG(NOP-4)
        IF (RF50>XAGE) FUNCT = FUNCT + 100000.d0 * (1.d0+RF50-XAGE)
        IF (RF50<1.d0) FUNCT = FUNCT + 100000.d0 * (2.d0-RF50)
        IF (FUNCT>=100000.D0 ) RETURN
        CALL SETO (SELM(:,2),GG(NOP-7),RM50,MAXAGE)
        CALL SETO (SELF(:,2),GG(NOP-5),RF50,MAXAGE)
C
C     TRIAL NF-U: selectivity decreases by SELDEC /yr after age 8
      ELSE IF (SELDEC>0.0D0) THEN
       DO 15 L=0,8
        SELM(L,1) = SELM(L,1)/SELM(8,1)
        SELF(L,1) = SELF(L,1)/SELF(8,1)
   15  CONTINUE
       DO 16 L=9,MAXAGE
        SELM(L,1) = SELM(L-1,1)*(1.D0-SELDEC)
        SELF(L,1) = SELF(L-1,1)*(1.D0-SELDEC)
   16  CONTINUE
      ENDIF    
C
C     Call MSYPAR to set A & Z 
      CALL MSYPAR (8)
C
C     SETKA sets up the pristine population vectors & advances them to year 0.
      CALL SETKA ()
C
      DO J = 1,NSTK
       IF (PMATF(0,J)/KMAT(J)<0.05d0)
     +    FUNCT = FUNCT + 100000.d0*(1.05d0-PMATF(0,J)/KMAT(J))
      ENDDO
      IF (FUNCT>=100000.D0 ) RETURN
C
C     Calculate Deviance using DEVIN(LAMBDA) so LAMBDA=GG(NGAMMA+NDELTA+1)
      FUNCT = FUNCT + DEVIN(GG(NGAMMA+NDELTA+1))

      IF (FUNCT.LT.BESTLL) THEN
        IF (PRDIAG>1) WRITE(*,801) FUNCT,GG(1:NOP)
        BESTLL = FUNCT
        DO I = 1,NOP
         GBEST(I) = F(I)
        ENDDO
C        PRINT '(A,I2,30F10.4)',' BEST:',PHASE,BESTLL,GG(1:NOP)
      ENDIF
      RETURN
 801  FORMAT("FUNCTION VALUE: ",F20.7 / (100F13.5))
      END

C ----------------------------------------------------------------------

      REAL(8) FUNCTION DEVIN (LAMBDA)
C
C     Calculate twice negative of the log-likelihood. This subroutine is used
C     when fitting the operating model to the abundance & other data
C
      USE DECLARECOM
C
      INTEGER IYR,J,K,N,I,NR2,Ns(6),II,I2,II2,NF1,L,NF,LASTTIYR
      REAL(8) X,Y,S,TOT1,ALPHA,TERM1,TERM2,TERM3,GAMMLN,LAMBDA,
     +       Utarg,Upred
      REAL(8) V3(3,3),VI3(3,3),V4(4,4),VI4(4,4),V2(2,2),
     +     VI2(2,2), ETA(6,-47:-22),E3(3),E4(4),E2(2),C,
     +     PRDM,PRDF,PRDRM,PRDRF,OC(MXFIT)
C     NF1 used to indicate what model fits are written to NAF.HIS
C
C     Zero fitting arrays
      FITDAT = 0.d0
      FITRSD = 0.d0
      FITRSS = 0.d0
      
      LikeAb = 0.d0
      LikeTag = 0.d0
      LikeAge = 0.d0
C      
C *** ABUNDANCES *****************
C
      DO 50 N = 1,NAB
C
C       PP = 1+ abundance of stock J in (true) Subarea K
C       Set Y = 1+ abundance in K (and apply correction)
C       If TARA=-1 skip this estimate: in Hyp. 7 & 8 the estimate has been
C                  combined with that for WI.
        IF (TARA(N)==-1) GO TO 50
        K   = TARA(N)
        IYR = TARYR(N)
        Y = 0.d0
        DO 10 J=1,NSTK
          Y = Y + PP(IYR,J,K)
   10   CONTINUE
        Y = Y*G0

        X = TAREST(N)
C       Compute log-likelihood contribution.     {Ref G(a)}
        IF (Y>1.d0) THEN
          TOT1 = LOG(Y/X)
        ELSE
          TOT1 = 1000.d0*(Y - 1.d0 - LOG(X))
        ENDIF
        S = TOT1*TOT1/TARVAR(N)
C
        FITDAT(N) = Y
        FITRSD(N) = S
        IF (Y<0.D0) STOP 'DEVIN ERROR'
        FITRSS(N) = log(X/(Y+1.d0))/SQRT(TARVAR(N))
        LikeAb = LikeAb + FITRSD(N)
        if (prdiag==1) PRINT '(a,3i4,3f12.3,f20.1)',
     +                       ' devin',n,k,iyr,y,x,tarvar(N),s
  50  CONTINUE
      NF = NAB
C
C *** HISTORICAL CATCH *****************
C     Check all historical catch is taken (CATT=actual catch modelled to be taken)
      TOT1 = 0.d0
      DO 251 J=1,NSTK
       DO 251 IYR = INITYR,-1
        TOT1 = TOT1 + CATT(IYR,J)
 251  CONTINUE
      CPENAL = ABS(TOTCAT-TOT1)
      NF=NF+1
      FITDAT(NF) = CPENAL
      FITRSD(NF) = CPENAL*10.d0
      FITRSS(NF) = TOTCAT-TOT1
      LikeCat  = FITRSD(NF)
C
C
C *** Tagging likelihood *****************
C
C     RECTAR(NR2,IYR,K)= Target no. recovered in area K in IYR for release NR
C     TAGU(NR2,IYR,K)  = Model  no. recovered
C     In deterministic case RECTAR=RECOBS(NR2,IYR,K)
C     Exc. same season recoveries; incorporates PSI=recovery reporting rate parameter

      DO 350 NR2=1,NROP
       NF=NF + 1
C      Skip out of loop if 0 releases in this op (ie in Hyp 7&8 when EG&WI are combined)
       IF (NREL2(NR2)==0) GO TO 350

C      Start Year AFTER tags were released (as NO same season recov in likelihood)
       IF(TAGEXC1 .EQ. 1) THEN
          LASTTIYR = TAGIYR(NR2)+2
       ELSE
          LASTTIYR = -1
       ENDIF
       
       DO 340 IYR = TAGIYR(NR2)+1,LASTTIYR
          DO 340 K = 1,NSUBA

         Utarg = RECTAR(NR2,IYR,K)
         Upred = TAGU(NR2,IYR,K)
C        Skip out of loop if no predicted or observed tags
         IF (Utarg.EQ.0.D0 .AND. Upred.EQ.0.D0) GOTO 340

         FITDAT(NF) = FITDAT(NF) + Upred
         FITRSS(NF) = FITRSS(NF) + ABS(Upred - Utarg)
         ALPHA = Upred * LAMBDA
         TERM1 = GAMMLN(Utarg+ALPHA) - GAMMLN(Utarg+1.D0) -GAMMLN(ALPHA)
         TERM2 = ALPHA * LOG(LAMBDA/(LAMBDA+1.D0))
         TERM3 = Utarg * LOG(1.D0/(LAMBDA+1.D0))
         FITRSD(NF) = FITRSD(NF) + TERM1 + TERM2 + TERM3
         
 340   CONTINUE
       FITRSD(NF) = -2.D0 * FITRSD(NF) * TAGWT
       LikeTag = LikeTag + FITRSD(NF)

 350  CONTINUE

C *** Catch at age likelihood ************
C
      NF1 = NF
      SIGSQ = 0.d0

      DO 660 K=1,NSUBA 
       DO 660 IYR=ICBYA1,ICBYA2
        IF (OCAM(IYR,K)+OCAF(IYR,K) > 0.d0) THEN
          PRDM = SUM(PCBYAM(IYR,1:MAXAGE,K))
          PRDF = SUM(PCBYAF(IYR,1:MAXAGE,K))
          DO 680 L=1,MAXAGE
            NF = NF + 1
            PRDRM = MAX(PCBYAM(IYR,L,K)/MAX(PRDM,1d0),1d-10)
            FITDAT(NF) = PRDRM
            OC(NF) = OCBYAM(IYR,L,K)
            IF (OC(NF) > 0.d0) THEN
              IF(PRDM < 3d-10) THEN
C               Observed catch at age but no model whales in K in IYR!
                FITRSS(NF) = 1d10
              ELSE
                FITRSS(NF) = (LOG(OC(NF))-LOG(PRDRM))**2
              END IF
C           ELSE FITRSS(NF) = 0  (initial value)
            END IF

            NF = NF + 1
            PRDRF = MAX(PCBYAF(IYR,L,K)/MAX(PRDF,1d0),1d-10)
            FITDAT(NF) = PRDRF
            OC(NF) = OCBYAF(IYR,L,K)
            IF(OC(NF) >  0.d0) THEN
              IF(PRDF < 3d-10) THEN
                FITRSS(NF) = 1d10
              ELSE
                FITRSS(NF) = (LOG(OC(NF))-LOG(PRDRF))**2
              END IF
            END IF

            SIGSQ = SIGSQ + OC(NF-1)*FITRSS(NF-1) + OC(NF)*FITRSS(NF)
 680      CONTINUE
         ENDIF
 660  CONTINUE
      SIGSQ = SIGSQ/(NF-NF1)

      DO N = NF1+1,NF
       IF (OC(N)>0.D0) THEN
         TERM1 = 0.5*LOG(SIGSQ/OC(N))
         TERM2 = OC(N)*FITRSS(N)/(2.d0*SIGSQ)
         FITRSD(N) = WtAge*(TERM1 + TERM2)
         LikeAge = LikeAge + FITRSD(N)
C         IF (PRDIAG==1) PRINT '(I5,9F10.4)', N,WtAge,OC(N),
C     +      FITDAT(N),FITRSS(N), FITRSD(N),TERM1,TERM2,LIKEAGE
C      ELSE FITRSD(N) = 0 (initial value)
       ENDIF
      ENDDO

C *** End C@age ***

C *** CPUE likelihood *****************
      IF (OPTCPE.GT.0) THEN
        Ns(1)=NF
        DO I=2,6
          Ns(I)=Ns(I-1)+Np(I-1)
        ENDDO
C       Compute the q's {Eqn G.11}
        DO 400 I=1,NCPUE
          qCPUE(I) = 0.d0
          DO 410 II=1,Np(I)
            X = CPUE(I,II)
            IYR = CPUEYR(I,II)
            K = CPUEK(I)
            Y = SELN(IYR,K)
            qCPUE(I) = qCPUE(I) + LOG(X/Y)
  410     CONTINUE
          qCPUE(I) = EXP(qCPUE(I)/Np(I))
  400  CONTINUE

C      Likelihood for the later series (1-4)
C      ****Note: hardwired years throughout this CPUE section**** ######
C      ****      Including the array size (years) of ETA
C
       DO 420 IYR=1962+ISCALE,1987+ISCALE
        DO 420 I = 1,4
         DO 420 II = 1,Np(I)
           IF (CPUEYR(I,II).EQ.IYR) THEN
             X = CPUE(I,II)
             K = CPUEK(I)
             Y = SELN(IYR,K)
             IF (X>0.D0) ETA(I,IYR) = LOG(X/(qCPUE(I)*Y))
           ENDIF
  420  CONTINUE
C      Split it into the 4 different periods
C      Period 1: 1962 to 1965: only series 1 to 3
       DO 425 I=1,3
        DO 425 II=1,3
          V3(I,II) = CPUEV(I,II)
  425  CONTINUE
       CALL INVM(V3,VI3,3,3)
       DO 430 IYR=1962+ISCALE,1966+ISCALE
         DO I=1,3
           K = CPUEK(I)
           Ns(I) =Ns(I)+1
           NF = Ns(I)
           E3(I) = ETA(I,IYR)
           FITDAT(NF) = qCPUE(I)*SELN(IYR,K)
           FITRSS(NF) = ETA(I,IYR)
         ENDDO
         CALL MultEVE(E3,VI3,C,3)
         FITRSD(NF) = CPEWT*C
  430  CONTINUE
C      Period 2: 1966 to 1982: series 1 to 4
       DO 525 I=1,4
        DO 525 II=1,4
          V4(I,II) = CPUEV(I,II)
  525  CONTINUE
       CALL INVM(V4,VI4,4,4)
       DO 530 IYR=1966+ISCALE,1982+ISCALE
         DO I=1,4
           K = CPUEK(I)
           Ns(I) =Ns(I)+1
           NF = Ns(I)
           E4(I) = ETA(I,IYR)
           FITDAT(NF) = qCPUE(I)*SELN(IYR,K)
           FITRSS(NF) = ETA(I,IYR)
         ENDDO
         CALL MultEVE(E4,VI4,C,4)
         FITRSD(NF) = CPEWT*C
  530  CONTINUE
C      Period 3: 1983 to 1985: series 1, 3 and 4
       I = 0
       DO 620 I2=1,4
         II = 0
         IF (I2/=2) THEN
           I = I+1
           DO II2=1,4
             IF (II2/=2) THEN
               II = II+1
               V3(I,II) = CPUEV(I2,II2)
             ENDIF
           ENDDO
         ENDIF
  620  CONTINUE
       CALL INVM(V3,VI3,3,3)
       DO 630 IYR=1983+ISCALE,1985+ISCALE
         I = 0
         DO I2=1,4
           IF (I2/=2) THEN
             K = CPUEK(I2)
             I = I+1
             Ns(I2) = Ns(I2)+1
             NF = Ns(I2)
             E3(I) = ETA(I2,IYR)
             FITDAT(NF) = qCPUE(I2)*SELN(IYR,K)
             FITRSS(NF) = ETA(I2,IYR)
           ENDIF
          ENDDO
          CALL MultEVE(E3,VI3,C,3)
          FITRSD(NF) = CPEWT*C
  630  CONTINUE
C      Period 4: 1986 to 1987: series 3 to 4
       I = 0
       DO 725 I2=3,4
         I=I+1
         II = 0
         DO 725 II2=3,4
           II = II+1
           V2(I,II) = CPUEV(I2,II2)
  725    CONTINUE
         CALL INVM(V2,VI2,2,2)
         DO 730 IYR=1986+ISCALE,1987+ISCALE
         I=0
         DO I2=3,4
           K = CPUEK(I2)
           I=I+1
           Ns(I2)= Ns(I2)+1
           NF    = Ns(I2)
           E2(I) = ETA(I2,IYR)
           FITDAT(NF) = qCPUE(I2)*SELN(IYR,K)
           FITRSS(NF) = ETA(I2,IYR)
         ENDDO
         CALL MultEVE(E2,VI2,C,2)
         FITRSD(NF) = CPEWT*C
  730  CONTINUE
C
C      Now the older CPUE series {Eqn G.14}
       DO 800 I=5,6
        DO 810 II=1,Np(I)
         NF = NF+1
         IYR=CPUEYR(I,II)
         K=CPUEK(I)
         Y = qCPUE(I)*SELN(IYR,K)
         TOT1 = LOG(CPUE(I,II)/Y)
         S=TOT1*TOT1/(2.0d0*CPUEsig(I)*CPUEsig(I))
         FITDAT(NF) = Y
         FITRSS(NF) = TOT1
         FITRSD(NF) = CPEWT*2.D0 * S
  810   CONTINUE
  800  CONTINUE
      ENDIF
C *** End CPUE likelihood *****************


C     Sum deviance for all conditioning parameters
      DEVIN = SUM(FITRSD(1:NF))
C     Store NF (used in printing conditioning output)

C     Temp: need to inc. catch @ age & cpue in NFIT ##############################
      NFIT = NF1
C
      RETURN
      END
C
C
C     -------------------------------------------------------------------
C     SURVEY-RELATED SUBROUTINES
C       RESET   - Resets the survey arrays.  Called by MAIN
C                 Also resets some other bits
C       SURVEY  - Performs a survey in subareas defined by RUNSUR.
C                 Called by MAIN
C       CONSUR  - Sets sightings estimates for groups of subareas
C                 Called by MAIN
C
C ------------------------------------------------------------------------------
C
      SUBROUTINE RESET (SIGHT,CV)
C
C     Reset is called in year 0 (after the premanagement population is setup)
C     Set TAU2 & KSGT arrays. Reset SIGHT(IYR,K), CV & catch arrays
C
      USE DECLARECOM
C
      REAL(8) SIGHT(I1SYR:IENDYR,MXAR),CV(I1SYR:IENDYR,MXAR),
     +        SUM2,P,PSGT(I1SYR:0),PROP
      INTEGER IYR,K,J,NCV,KM
C
C     Initialise FUTURE SIGHT & CATCH arrays (don't reset historical values)
      DO 5 IYR = I2SYR,IENDYR
       DO 5 K = 1,MXAR
        SIGHT(IYR,K) = -1.d0
        CV(IYR,K) = -1.d0
    5 CONTINUE
      DO 9 IYR = 0,IENDYR
       DO 8 K = 1,MXSUBA
        CATCHM(IYR,K) = 0.d0
        CATCHF(IYR,K) = 0.d0
    8  CONTINUE
       CATKA(IYR,1:MXAR) = 0.D0
    9 CONTINUE

C     Set KSGT & PSGT = 1+ population in MANAGEMENT subareas KM           {Eqn F.2}
C     PP= 1+ popln in MODEL area K
C     Hyp 7&8 (NSUBA/=MXSUBA): split the combined EG+WI area (K=3) into
C             EG (KM=3) & WI (KM=4) using EGWIP, the mean historical EG proportion
      KSGT = 0.d0
      TAU2 = 0.d0
      DO 25 KM=1,MXSUBA
        K = KM
        PROP = 1.D0

        IF (NSUBA/=MXSUBA) THEN
          IF (KM>3) K=KM-1
          IF (KM==3) THEN
            PROP = EGWIP
          ELSE IF (KM==4) THEN
            PROP = 1.D0 - EGWIP
          ENDIF
        ENDIF

        KSGT(KM)= SUM(PP(INITYR,1:NSTK,K)) * PROP
        DO IYR = I1SYR,0
         PSGT(IYR) = SUM(PP(IYR,1:NSTK,K)) * PROP
        END DO

C       Set TAU2=TAU**2 using mean CV**2 of the conditioning abundance estimates
C       & the mean depletion in the corresponding years. {Eqn F.4}
        SUM2 = 0.d0
        NCV = 0
        P = 0.d0
        DO 20 IYR = I1SYR,I2SYR-1
          IF (CV2TAU(IYR,KM) /= -1.d0) THEN
            SUM2 = SUM2 + CV2TAU(IYR,KM)
            NCV = NCV+1
            P = P + PSGT(IYR)
          ENDIF
   20   CONTINUE

C       Set Tau**2                                                 {Eqn F.3}
        IF (NCV>0) TAU2(KM) = SUM2/(NCV*(0.12d0+0.025d0*KSGT(KM)*NCV/P))
C        PRINT'(A,2I4,6F10.3)',' TAU2:',KM,NCV,TAU2(KM),P,SUM2,KSGT(KM),SSAMPL
C
   25 CONTINUE

C     SSAMPL = Increase in survey sampling intensity (=1 in basecase).
C     In trial NF-P survey intensity is tripled in WI by setting TAU2=TAU2/3 & DF=3DF
      TAU2(4) = TAU2(4) / SSAMPL
C
      RETURN
      END
C
C ------------------------------------------------------------------------
C
      SUBROUTINE SURVEY (SIGHT,CV,IY)
C
C     SURVEY generates survey data : sets abundance estimate (SIGHT) &
C            CV in each subarea for year IY (=IYR-1).

      USE DECLARECOM
C
      REAL(8) SIGHT(I1SYR:IENDYR,MXAR),CV(I1SYR:IENDYR,MXAR)
      REAL(8) ABUND,BETA2K,SIGMA,Y,W,SIG2,CHISQ,RANNO,
     +     GAMMA,POISSN,RAN1,XNORM,MA1(55),MA2(55),MA3(55),SGT,SCV,PSGT,
     +     ALPHA2,ALPHABETA,BETABETA,MU,VAR,PROP,GENBET,P,PROP2
      INTEGER IY,K,INEXT1,INXTP1,INEXT2,INXTP2,INEXT3,INXTP3,J,
     +     I2SEED,I3SEED,IPSEED,KM
      EXTERNAL GAMMA,POISSN,RAN1,XNORM,GENBET
      SAVE INEXT1,INXTP1,INEXT2,INXTP2,INEXT3,INXTP3,MA1,MA2,MA3
C
C     Hyp 7&8 (NSUBA/=MXSUBA): set propn of EG+WI whales that are in EG using
C         a BETA generator. EGWIP = mean historical EG propn with variance EGWIV
    
      IPSEED = INT(-RAN1(ISEED2,MA2,INEXT2,INXTP2)*100000.d0)
      IF (NSUBA/=MXSUBA) THEN
        ALPHABETA = EGWIP*EGWIP*(1.d0-EGWIP)/EGWIV - EGWIP
        BETABETA = ALPHABETA*(1.D0-EGWIP)/EGWIP
        PROP2 = GENBET(ALPHABETA,BETABETA,IPSEED)
      ENDIF

      

C     Generate the survey results
      PROP2 =0.5D0
      DO 100 KM=1,MXSUBA
C       First increment all the random number generators for every area
        Y = XNORM (1.d0,0.d0,ISEED1,MA1,INEXT1,INXTP1)
        RANNO = RAN1(ISEED1,MA1,INEXT1,INXTP1)
        I2SEED = INT(-RAN1(ISEED2,MA2,INEXT2,INXTP2)*100000.d0)
        I3SEED = INT(-RAN1(ISEED3,MA3,INEXT3,INXTP3)*100000.d0)
C       Leave loop if no survey this year/area (SIGHT & CV initialised to -1 in RESET)
        IF (RUNSUR(IY,KM)==0) GOTO 100

C       Set PSGT = Sum of 1+ population in MANAGEMENT subarea KM         {Eqn F.2}
C       PP= 1+ popln in MODEL area K
        K = KM
        PROP = 1.D0

        IF (NSUBA/=MXSUBA) THEN
C         Hyp 7&8: split combined EG+WI area (K=3) into EG & WI (KM=3 & 4)
          IF (KM>3) K=KM-1
          IF (KM==3) PROP = PROP2
          IF (KM==4) PROP = 1.D0 - PROP2
        ENDIF
        
       
        
        PSGT = SUM(PP(IY,1:NSTK,K)) * PROP
C       Adjust abundance for bias (whales missing from the survey area)
        ABUND = PSGT*G0*BIAS
C     
C       Specify parameters used to generate estimates. TAU2 = TAU**2  {Eqns F.5 & F.6}
        ALPHA2 = 0.120d0*TAU2(KM)
        SIGMA  = SQRT(LOG(ALPHA2 + CVADD(KM)**2 +1.d0))
        BETA2K = 0.025d0*TAU2(KM)*KSGT(KM)*G0*BIAS
C
C       Generate Poisson component
        W = ABUND/BETA2K
        IF (W < 70.d0) W = POISSN (W, RANNO)

C
C       Construct estimate {Eqn F.1} (unless OPTDET = 1)
        SGT = BETA2K * EXP(Y*SIGMA) * W
        IF (OPTDET==1) THEN
          SGT = ABUND
          SCV = SQRT (ALPHA2 + BETA2K/ABUND)
C
C       CV Generation {Eqn F.5} (unless estimate is non-zero)
        ELSE IF (SGT.GT.0.d0) THEN
          SIG2  = LOG(1.d0 + ALPHA2 + BETA2K/SGT)
          CHISQ = GAMMA (DF,SQRT(DF+DF),I3SEED)
          SCV = SQRT(SIG2 * CHISQ/DF)
        ELSE
C         Zero estimate. Store Z(i), the Poisson multiplier, in SCV
          SCV =  BETA2K * EXP(Y)
        ENDIF
C
C       Now generate proportion N of 60N for NF-Q trials & Variant 4 (in EI/F & WI)
C           from a beta distribution (with Mean Mu & Var based on obs.data).
C       For all other trials/areas/variants  Prop = Mu = SRATIO = 1.
        MU = SRATIO(KM)
        VAR = VRATIO(KM)**2
        IF (VAR<=0.D0 .OR. MU==1.D0) THEN
         PROP = MU
        ELSE
         ALPHABETA = MU*MU*(1.D0-MU)/VAR-MU
         BETABETA = ALPHABETA*(1.D0-MU)/MU
         PROP = GENBET(ALPHABETA,BETABETA,I2SEED)
        ENDIF
C
C       Store results.
        SIGHT(IY,KM) = SGT*PROP
        CV(IY,KM) = SCV
C       PRINT '(A,3I5,7F12.4)',' SURVEY',IY,KM,K,PSGT,SGT,SGT*PROP,SCV,
C     +       prop,ALPHABETA,BETABETA,prop2
C
100   CONTINUE

      RETURN
      END
C
C     -------------------------------------------------------------------
C
      SUBROUTINE CONSUR (SIGHT,CV,IY)
C
C     This subroutine constructs survey results for areas larger than subareas
C     for use by CLA (so uses Management areas through out).
C     SIGHT & CV arrays were initialised to -1 in RESET
C
      USE DECLARECOM
C
      COMMON /STORE/ IASESS, NCAP, ATYPE, ITYPE
      INTEGER IASESS(I1SYR:IENDYR),NCAP(2),ATYPE(2,MXAR),ITYPE
C
      REAL(8) SIGHT(I1SYR:IENDYR,MXAR),CV(I1SYR:IENDYR,MXAR),
     +       SGT,SCV,SUM1,SUM2,SK
      INTEGER IY,K,KA,I,IY2,N
C
C     Construct surveys for groups of subareas as defined by MAP
      DO 100 KA = MXSUBA+1,MXAR
C
C       Skip loop if not constructing a survey this year (RUNSUR=0)
        IF (RUNSUR(IY,KA).EQ.0) GO TO 100

C       Calculate the estimate (SGT) and the CV (SCV)
        SGT = 0.d0
        SCV = 0.d0
        N = 0

C       If RUNSUR>1 then combine surveys over previous RUNSUR years
        IY2 = IY-RUNSUR(IY,KA)+1

        DO 80 K = 1,MXSUBA
C         Skip loop if subarea K is not part of assessment area KA (MAP=0)
          IF (MAP(KA,K).EQ.0) GO TO 80

C           There may be more than one historical survey/subarea in assessment area
C           1st find inverse variance weighted average of surveys in same subarea
            SUM1 = 0.d0
            SUM2 = 0.d0
            DO 70 I=IY2,IY
              IF (SIGHT(I,K).GE.0.000000001d0) THEN
                N = N+1
                SUM1 = SUM1 + LOG(SIGHT(I,K))/CV(I,K)**2
                SUM2 = SUM2 + 1.d0 /CV(I,K)**2
              ENDIF
  70        CONTINUE

C           Add result to assessment area totals (1/CVK**2 = SUM2)
            IF (SUM1.GT.0.d0) THEN
              SK = EXP (SUM1/SUM2)
              SGT = SGT + SK
              SCV = SCV + (SK**2)/SUM2
            ENDIF

  80    CONTINUE
C
C       Store result in mean year (rounded up). Assume no 0s in combined area
C          (if 0 combined estimates do occur, Poisson multipler must be added here) #####
        I = (IY2+IY)/2
        SIGHT(I,KA) = SGT
        IF (N==0 .or. SGT==0.d0) then
          WRITE(41,'(A,5I5,2F12.4)')' CONSUR ERROR',I,KA,N,IY2,IY,SGT
          STOP ' CONSUR ERROR: no surveys in area KA'
        ENDIF
        CV(I,KA) = SQRT(SCV)/SGT
        WRITE(41,'(A,5I5,2F12.4)')' CONSUR',I,KA,N,IY2,IY,SGT,CV(I,KA)
100   CONTINUE
C
      RETURN
      END
C
C     ------------------------------------------------------------------
C     OUTPUT-RELATED SUBROUTINES
C       PNTOUT  - Called by MAIN
C       REPORT  - CAlled by MAIN
C       PCOND   - Called by MAIN
C       PSTATS  - Called by MAIN
C ----------------------------------------------------------------------


      SUBROUTINE PNTOUT (NTRIAL,NYEAR,IPNT, DESC,TYPEDS)
C     Print list of parameters & options. Check values are within range

      USE DECLARECOM

      COMMON /MANPAR/ PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL
      REAL(8) PPROB,PYMIN,PYMAX,PNYSTP,PKSTEP,PDSTEP,
     1     DTMIN,DTMAX,PBMIN,PBMAX,PNBSTP,PSCALE,PHASET,PHASEP,
     1     PCYCLE,PLEVEL,PSLOPE,ACCTOL

      INTEGER NTRIAL,NYEAR,IPNT,J,I
      CHARACTER DESC*50, TYPEDS*45, OPT(-1:6)*50
      LOGICAL ERR

C     Check parameters are within allowed ranges. If not print error message
      IF (NYEAR<0 .OR. NYEAR>IENDYR)STOP' *** ERROR IN SIMULATION YEARS'
      IF (DOCON<-1.OR. DOCON> 3)    STOP' *** ERROR IN DOCON'
      IF (MSYL <0.d0.OR.MSYL>=1.d0) STOP' *** ERROR IN MSYL'
      DO 10 J=1,NSTK
        IF (MSYR(J)<0.d0 .OR.MSYR(J)>=1.d0) STOP' *** ERROR IN MSY rate'
   10 CONTINUE

      ERR = .FALSE.
      OPT(-1) = '**** ERROR: UNKNOWN VALUE'
C     Print parameter list
      WRITE (IPNT,'(1X,A,3X,A/19X,A/)') REF,DESC,TYPEDS
      WRITE (IPNT,'(A,T41,8F7.3)') ' CLA parameters',
     1       PPROB,PNYSTP,PKSTEP,PDSTEP,PNBSTP,PCYCLE,PLEVEL

      OPT(0) = 'Stochastic run'
      OPT(1) = 'Deterministic run'
      I = OPTDET
      IF (I.LT.0 .OR. I.GT.1) CALL SETERR (ERR,I)
      WRITE (IPNT,100) 'Option controlling stochasticity',I,OPT(I)

      WRITE (IPNT,120)
      WRITE (IPNT,120) 'Constant assessment interval, in years',PCYCLE
      WRITE (IPNT,100) 'Number of trials',NTRIAL
      WRITE (IPNT,100) 'Number of years in simulation',NYEAR
      WRITE (IPNT,100) 'Number of years of premanagement catch',-INITYR
      WRITE (IPNT,100) '1st year of premanagement catch', INITYR-ISCALE
      WRITE (IPNT,100) 'Number of stocks',NSTK
      WRITE (IPNT,100) 'Number of subareas',NSUBA
      WRITE (IPNT,'(A,T41,F9.3, F7.3)') ' Maturity parameters',MAT1,MSIG
      WRITE (IPNT,'(A,T41,F9.3,2F7.3)') ' Mortality parameters',MORT
      WRITE (IPNT,100) 'Maximum age class',MAXAGE
      WRITE (IPNT,'(/ T41,10(A9,I3,:2X))') ('Stock',J,J=1,NSTK)
      WRITE (IPNT,120) 'MSYL', MSYL
      WRITE (IPNT,110) 'MSY rate', MSYR(1:NSTK)
      WRITE (IPNT,110) 'Density dependent exponent (Z)',Z(1:NSTK)
      WRITE (IPNT,110) 'Resilience parameter (A)', A(1:NSTK)
      WRITE (IPNT,100)
      WRITE (IPNT,110) 'Bias in absolute abundance estimates',BIAS
      WRITE (IPNT,120) 'Degrees of Freedom',DF
      WRITE (IPNT,120) 'g(0) factor',G0
      IF (ERR) STOP

  100 FORMAT (1X,A,T41,I6,4X,A)
  110 FORMAT (1X,A,T41,F12.5,10F14.5)
  120 FORMAT (1X,A,T41,F9.2,14F8.0)
  200 FORMAT (/' ***** ERROR: ',A,' OUT OF RANGE')

      RETURN
      END

      SUBROUTINE SETERR (ERR,I)
      INTEGER I
      LOGICAL ERR
      ERR = .TRUE.
      I = -1
      RETURN
      END

C ------------------------------------------------------------------------

      SUBROUTINE REPORT (IP,NYEAR)

C     Subroutine REPORT prints PMAT & Catches by year

      USE DECLARECOM

      INTEGER PK(MXSUBA),K1(MXSTK),P1(MXSTK),CK(MXSUBA),IP,NYEAR,K,I,J,L
     +        ,CPK(MXSUBA)

    1 FORMAT (I4,100I7)
    2 FORMAT (A4,100(A5,I2))


C     Set pristine & yr 0 1+ population size by stock
      DO 5 J=1,NSTK
        K1(J)= NINT(SUM(PP(INITYR,J,1:NSUBA)))
        P1(J)= NINT(SUM(PP(0,J,1:NSUBA)))
   5  CONTINUE
      WRITE (IP,1)INITYR,K1(1:NSTK)
      WRITE (IP,1)  0,   P1(1:NSTK)

      WRITE (IP,2) 'Yr',('Fem:',J,J=1,NSTK),('PK:',K,K=1,NSUBA),
     +                  ('CJ:', J,J=1,NSTK),('CK:',K,K=1,NSUBA)
      OPEN (99, FILE=ref(1:7)//'.pop',ACCESS='APPEND')                          # matF by stock, 1+ by subarea
      OPEN (98, FILE=ref(1:7)//'.cat',ACCESS='APPEND')                          # catch by stck, 1+ by subarea

      OPEN (94, FILE=ref(1:7)//'.age',ACCESS='APPEND')
      WRITE(94,*) 'YEAR ','AREA ','AGE ','OBS-FEM ','PRD-FEM ',
     +     'OBS-M  ','PRD-M '
      DO 8 K=1,NSUBA
        IF (SUM(OCAM(ICBYA1:ICBYA2,K))>0) THEN
          DO 7 I = ICBYA1,ICBYA2
           DO 7 L=1,MAXAGE
            WRITE(94,*) I-ISCALE,K,L, OCBYAM(I,L,K),
     +               PCBYAM(I,L,K), OCBYAF(I,L,K),PCBYAF(I,L,K)
    7     CONTINUE
        END IF
    8 CONTINUE

      DO 20 I = INITYR,NYEAR
        DO 10 K=1,NSUBA
C         SUM 1+ no. & catch in subarea
          CK(K) = NINT(CATCHF(I,K)+CATCHM(I,K))
          PK(K) = NINT(SUM(PP(I,1:NSTK,K)))
          CPK(K)= NINT(SUM(PCBYAM(I,1:MAXAGE,K))
     +                +SUM(PCBYAF(I,1:MAXAGE,K)))
   10   CONTINUE

        WRITE (IP,1) I-ISCALE,NINT(PMATF(I,1:NSTK)),PK(1:NSUBA),
     +                        NINT(CATT(I,1:NSTK)), CK(1:NSUBA)
C       NAF.STS(99): write Mat females & 1+ pop in subarea K
        WRITE (99,1) I-ISCALE,NINT(PMATF(I,1:NSTK)), PK(1:NSUBA)
        IF (I<100) WRITE(98,1) I-ISCALE,NINT(CATT(I,1:NSTK)),CK(1:NSUBA)
     +              ,CPK(1:NSUBA)

C       OUTPUT 1st year catches only for real application test
        IF (DOCON==3.AND.I==0) WRITE (40,1) I-ISCALE,CK(1:NSUBA)
   20 CONTINUE
      CLOSE (98)
      CLOSE (99)
      CLOSE (94)

      RETURN
      END
C
C ------------------------------------------------------------------------
C
      SUBROUTINE PCOND (NT,FIT)

C     PCOND prints the conditioning diagnostics

      USE DECLARECOM
      REAL(8) FIT,FA(MXSUBA),FP(MXSUBA),SA(MXSUBA,MXSUBA),
     +       SP(MXSUBA,MXSUBA),SREL(MXSUBA),
     +       YA(I1YR:0,MXSUBA),YP(I1YR:0,MXSUBA),TARG(40)
      INTEGER NT,N,J,I,K,KT,TARYR2(MXAB),NA,NR,IYR,I1,I2
      CHARACTER C2(34:37)*15,FMT2*22
      CHARACTER(2) :: COLS(18:20) = (/'18','19','20'/)
      DATA C2 /'VALUES ACHIEVED','RESIDUALS','RESIDUALS','TARGETS'/

    1 FORMAT (A,I5,4X,A)
    2 FORMAT (2I5,1x,99F7.2)
    4 FORMAT (A,99G6.0)
    5 FORMAT (A,99F8.3)
    6 FORMAT (A,99F6.2)
    7 FORMAT (A3,A14,99(A4,I2.2))
    8 FORMAT (I3,F14.4,99F6.2)
    9 FORMAT (2X,A,99(1X,A,G2.0))

C     Initialisation NT = -1: PRINT OUTPUT HEADINGS
      IF (NT==-1) THEN

        WRITE(95,1) REF,NSTK,'Trial Constants'
        WRITE(95,'(A3,A11,99(A12,I2.0))')'N','Fit',
     +      ('Gamma',I,I=1,NGAMMA),('Disp',I,I=1,2),'Lambda',0,'PsiC',0,
     +       'PsiWI+',0,('Kmat',J,J=1,NSTK),('Depn',J,J=1,NSTK)
        WRITE (97,1) Ref,NSTK,'Mixing Matrices'
        WRITE (99,1) Ref,NSTK,'Mature F by stock, mean 1+ by subarea'
        WRITE (99,9) 'Yr',(' Pf:',J,J=1,NSTK),(' P1:',K,K=1,NSUBA)
        WRITE (98,1) Ref,NSTK,'Catches by stock and by subarea'
        WRITE (98,9) 'Yr',(' Cf:',J,J=1,NSTK),('Cat:',K,K=1,NSUBA)
        WRITE (94,1) REF,NSTK,'Catch at age'
        IF (OPTCPE>0) WRITE (96,'(2A/I2,A/6I3,A)') Ref,
     +     ' Observed & model-predicted CPUE',NCPUE,'#No CPUE series',
     +      Np(1:NCPUE),'#No data point per series'
        DO I=94,99
          CLOSE(I)
        END DO

        DO I=1,NAB
          TARYR2(I) = TARYR(I)-ISCALE-1900
          IF (TARYR2(I)>=100) TARYR2(I) = TARYR2(I)-100
        END DO
        DO 10 I=34,37
         WRITE (I,1) REF,NSTK,C2(I)
         WRITE (I,7) 'N','Fit',(ANAM2(TARA(NA)),TARYR2(NA),NA=1,NAB),
     +          'NoC',0,('Tg',N,N=1,NROP)
         CLOSE(I)
   10   CONTINUE
        RETURN
      END IF

      OPEN (95, FILE=ref(1:7)//'.all',ACCESS='APPEND')
      WRITE(95,'(I3,F13.6,1x,99G14.8e1)') NT,FIT,EXP(GBEST(1:NOP))
      CLOSE(95)

C     Print results of conditioning trial NT (inc. NT=0=deterministic run)
      OPEN (34,FILE=ref(1:7)//'.his',ACCESS='APPEND')
      OPEN (35,FILE=ref(1:7)//'.rsd',ACCESS='APPEND')
      OPEN (36,FILE=ref(1:7)//'.rss',ACCESS='APPEND')
      OPEN (37,FILE=ref(1:7)//'.tar',ACCESS='APPEND')
      FMT2 = '(I3,F14.4,' // COLS(NAB) // 'I6,99F6.3)'
C     Needs to be changed to print more of age fit & cpue
      WRITE(34,*) NT,FIT,NINT(FITDAT(1:NAB)),FITDAT(NAB+1:NFIT),
     +               LikeAge
      WRITE(35,8) NT,FIT,(FITRSD(I),I=1,NFIT)
      WRITE(36,8) NT,FIT,(FITRSS(I),I=1,NFIT)

      TARG = 0.d0
      DO 20 NR=1,NROP
       DO 20 IYR = ITAG1,-1
        DO 20 K = 1,NSUBA
         TARG(NR) = TARG(NR) + RECTAR(NR,IYR,K)
  20   CONTINUE
      WRITE(37,*) NT,FIT,NINT(TAREST(1:NAB)),0.,TARG(1:NROP)
      DO I=34,37
        CLOSE(I)
      END DO

C     Write 97: NAF.MIX = Mixing matrices (N=1 for all these trials)
      OPEN (97, FILE=ref(1:7)//'.mix',ACCESS='APPEND')
      WRITE (97,'(I3,'' Fit:'',F12.5,2F8.4)')NT,FIT
      DO J = 1,NSTK
        WRITE(97,'(2I5,10F7.4)')J,INITYR-ISCALE,JV(1:NSUBA,J,1)
      END DO
      IF (OPTMIX>1) THEN
        WRITE (97,*)
        DO 52 IYR=IMIX1-1,0
          CALL SETV (IYR)
          DO 52 J = 1,NSTK
           WRITE(97,'(2I5,1x,10F7.4)') J,IYR-ISCALE,V(1:NSUBA,J,1)
   52   CONTINUE
      ENDIF
      CLOSE (97)

      IF (NT.GE.1) RETURN

C     Print deterministic fit diagnostics (NT=0) to NAF.SUM (unit 38)   # deterministic summary
      WRITE(38,1) REF,NSTK,'  Summary'
      WRITE(38,'(A3,A13,50(A11,I2.0))')'N','Fit',('Gamma',I,I=1,NGAMMA),
     + ('Disp',I,I=1,2),'Lambda',0,'PsiC',0,'PsiWI+',0,
     + ('Kmat',J,J=1,NSTK)
      WRITE(38,'(I3,1x,99G13.6e1)') 0,FIT,EXP(GBEST(1:NOP))
      WRITE(38,'(/7X,6A8)')  'C/WG','EG','WI','EI/F','N','S'
      WRITE(38,5) 'DEPL 0',(PMATF(0,J)/KMAT(J),J=1,NSTK)
      WRITE(38,5) '  1915',(PMATF(1915+ISCALE,J)/KMAT(J),J=1,NSTK)
      WRITE(38,'(A,6I8)') ' KMAT  ', NINT(KMAT(1:NSTK))
      WRITE(38,'(A,6I8)') ' P 0   ', NINT(PMATF(0,1:NSTK))
      WRITE(38,'(/3(A,3F6.3)/)')' DISPL 1,2:',AVDISP(1),D(1,2),D(2,1),
     +                      '     DISPL 2,3:',AVDISP(2),D(2,3),D(3,2),
     +                      '     DISPL 3,4:',AVDISP(3),D(3,4),D(4,3)
C     2. Abundances
      WRITE(38,9) '  Ab:',(ANAM2(TARA(NA)),TARYR2(NA),NA=1,NAB),'NoC',0
      WRITE(38,4) ' Target', NINT(TAREST(1:NAB)), 0
      WRITE(38,4) ' Fitdat', NINT(FITDAT(1:NAB+1))
      WRITE(38,6) ' FitRSS',(FITRSS(I),I=1,NAB+1)
      WRITE(38,6) ' Resid', (FITRSD(I),I=1,NAB+1)
      WRITE(38,*)
C     3. Tags
      WRITE(38,9) ' Tags',('Tg',NR,NR=1,NROP)
      I1 = ISCALE-1900
      WRITE(38,9) '     ',(ANAM2(TAGA(NR)),TAGIYR(NR)-I1,NR=1,NROP)
      WRITE(38,4) ' Target', NINT(TARG(1:NROP))
      WRITE(38,6) ' Fitdat', FITDAT(NAB+2:NFIT)
      WRITE(38,6) ' FitRSS', FITRSS(NAB+2:NFIT)
      WRITE(38,6) ' Resid',  FITRSD(NAB+2:NFIT)

C     Check output & print tag data
      WRITE(38,'(/A)') ' Observed + Predicted Tag Recoveries'
      WRITE(38,'(3A5," No.  Recov:",A,9A8)') 'Rel:','Area','Yr',ANAM2
      SA=0.D0
      SP=0.D0
      SREL = 0.D0
      DO 350 NR =1,NROP
        FA = 0.D0
        FP = 0.D0
        KT = TAGA(NR)
        DO 320 K = 1,NSUBA
          DO 310 I = ITAG1,IYR
            FA(K) = FA(K) + RECTAR(NR,I,K)
            FP(K) = FP(K) + TAGU(NR,I,K)
  310     CONTINUE
          SP(KT,K) = SP(KT,K) + FP(K)
          SA(KT,K) = SA(KT,K) + FA(K)
  320   CONTINUE
        WRITE (38,'(3I5,9F8.2)') NR,KT,TAGIYR(NR)-ISCALE,NREL(NR),FA
        WRITE (38,'(23X,9F8.2)') FP
        SREL(KT) = SREL(KT) + NREL(NR)
  350 CONTINUE
      WRITE(38,'(/A)')' Summary: Obs + Predict Recovery by release area'
      DO 360 KT = 1,NSUBA
        WRITE (38,'(5X,I5,5X,9F8.2)') KT,SREL(KT),SA(KT,1:NSUBA)
        WRITE (38,'(23X,9F8.2)')      SP(KT,1:NSUBA)
  360 CONTINUE

      DO 365 I=38,39
       WRITE(I,1)Ref,NSTK,'Tag Recoveries by release area, Obs+Predict'
       WRITE(I,'(2A5," Obs:",A,A4,5A5,3x,"Pred:",A,6A7)') 'RelA','Yr',
     +           ANAM2,ANAM2
  365 CONTINUE
      DO 390 KT = 1,NSUBA
        YA = 0.D0
        YP = 0.D0
        I1 = -3
        DO 380 NR =1,NROP
C         Find whales released in KT; I1 = year of 1st recoveries in area
C         (ignore same season recoveries so I1= year of 1st release + 1)
          IF (KT.EQ.TAGA(NR)) THEN
            I1 = MIN(I1,TAGIYR(NR)+1)
            DO 370 K = 1,NSUBA
             DO 370 I = TAGIYR(NR)+1,-3
              YA(I,K) = YA(I,K) + RECTAR(NR,I,K)
              YP(I,K) = YP(I,K) + TAGU(NR,I,K)
  370       CONTINUE
          ENDIF  
  380   CONTINUE
        DO 390 I = I1,-3
         I2 = I-ISCALE
         WRITE(38,2) KT,I2,YA(I,1:NSUBA),YP(I,1:NSUBA)
         WRITE(39,2) KT,I2,YA(I,1:NSUBA),YP(I,1:NSUBA)
  390 CONTINUE
      CLOSE (38)
      CLOSE (39)
      OPEN (38,FILE=ref(1:7)//'.sum',ACCESS='APPEND')

      IF (OPTCPE>0) THEN
        OPEN (96,FILE=ref(1:7)//'.cpe',ACCESS='APPEND')
        DO 391 I=1,NCPUE
         K=CPUEK(I)
         WRITE(96,'(/A4,(A7,I1),4X,A4)') 'Year', ' CPUE:',I,ANAM2(K)
         DO 391 N=1,Np(I)
          I2 = CPUEYR(I,N)
          WRITE(96,'(I4,20F8.3)')I2-ISCALE,CPUE(I,N),qCPUE(I)*SELN(I2,K)
  391   CONTINUE
        CLOSE (96)
      ENDIF

      RETURN
      END
C
C ------------------------------------------------------------------------
C
      SUBROUTINE PSTATS (IPNT,SIGHT)

C     PSTATS prints the final statistics

      USE DECLARECOM

      COMMON /STORE/ IASESS, NCAP, ATYPE, ITYPE
      INTEGER IASESS(I1SYR:IENDYR),NCAP(2),ATYPE(2,MXAR),ITYPE
      INTEGER IPNT,NZ(MXAR),NS(MXAR),NZTOT,NSTOT,I,KA
      REAL(8) SIGHT(I1SYR:IENDYR,MXAR)
C
C     Write no. of zero & nonzero estimates & times capping was applied (if any)
      NZTOT = 0
      NSTOT = 0
      NZ = 0
      NS = 0
      DO 10 KA=1,MXAR
        DO 9 I=0,IENDYR-1
          IF (SIGHT(I,KA) .NE. -1.d0) THEN
            NS(KA) = NS(KA) + 1
            IF (SIGHT(I,KA)<=0.d0) NZ(KA) = NZ(KA) + 1
          ENDIF
    9   CONTINUE
        NSTOT = NSTOT + NS(KA)
        NZTOT = NZTOT + NZ(KA)
   10 CONTINUE
      WRITE (IPNT,98) ' Zero estimates:',NZ
      WRITE (IPNT,98) ' Surveys:',       NS
      WRITE (IPNT,99) ' Total zero estimates used',NZTOT,NSTOT
      IF (NCAP(2)>0) WRITE (IPNT,99)' Capping was invoked',NCAP
   98 FORMAT (A,T17,18I5)
   99 FORMAT (A,I6,'  out of ',I6 :' times')

      RETURN
      END
C
C     ------------------------------------------------------------------
C     FUNCTIONS
C       GAMMLN  - Used by DEVIN to calculate log of Gamma function
C       MultEVE - Used by DEVIN to multiply matrices
C
C ----------------------------------------------------------------------
C
      REAL(8) FUNCTION GAMMLN(XX)
      REAL(8) COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX,XY
      INTEGER J
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      XY = XX
      IF (XY.LT.1.0D-10) XY = 1.0D-10
      X=XY-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C
C     -------------------------------------------------------------------

      SUBROUTINE MultEVE(A,B,C,NN)

C     MultEVE multiplies matrices: A(1,NN)*B(NN,NN)*AT(NN,1) to give C
 
      INTEGER NN,TT,kk
      REAL(8) A(NN),B(NN,NN),Ctmp(NN),C

C     Mutliply A*B
      DO 10 TT = 1,NN
        Ctmp(TT) = 0.d0
        DO 10 kk = 1,NN
          Ctmp(TT) = Ctmp(TT) + A(KK)*B(KK,TT)
   10 CONTINUE

C     Multiply B*AT
      C=0.d0
      DO kk=1,NN
        C = C + Ctmp(KK)*A(KK)
      ENDDO

      RETURN
      END
C
C     -------------------------------------------------------------------
C
      INCLUDE 'mannaf4.for'
      INCLUDE 'nrmp-i2.for'
      INCLUDE 'nrmp.for'
      INCLUDE 'clc-d2.for'
      INCLUDE 'pelagic1.for'
      INCLUDE 'pelagic2.for'
      INCLUDE 'toms178.for'
      INCLUDE 'bfgs.for'
