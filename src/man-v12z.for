C     Last change:  CA   29 Apr 2006    2:51 am
C     Program to run MANTST12 with ZERO catches
      INCLUDE 'mantst12.for'
      INCLUDE 'manextrd.for'

C *********************************************************************

      SUBROUTINE EXTRA (NTRIAL, OPTRAN, ERATE, ENDSUR, IOUT, PROG,ITUNE,
     +                  OPTDEPL) 

C     To avoid wasting time EXTRA sets the number of trials to 1 (except in 
C     random parameter or epidemic trials) and the year of last survey to 0 
C     Routine also used to open output files.

      INTEGER NTRIAL, OPTRAN,ENDSUR,IOUT,ITUNE,OPTDEPL
      DOUBLE PRECISION ERATE
      CHARACTER PROG*8

      
     
      PROG = 'MAN-V12Z'
      ITUNE = 0

      IF (OPTDEPL.EQ.0.AND.OPTRAN.NE.1 .AND. ERATE.EQ.0.0) NTRIAL = 1
      ENDSUR = 0

      RETURN
      END

C *********************************************************************

      SUBROUTINE RSETSS ()

C     Subroutine initialises parameters at start of each trial

      COMMON /CLCPAR/ IASESS, LSURV, PASESS, POUT
      INTEGER IASESS, LSURV, PASESS, POUT
      PASESS = 0
      END

C *********************************************************************

      SUBROUTINE MANSS (FCATCH, FSIGHT, CVEST, IYR, INITYR, ITEST)

      DOUBLE PRECISION FCATCH(-65:2000), FSIGHT(-1:2000), CVEST(-1:2000)
      INTEGER IYR,INITYR,ITEST

C     FCATCH(IYR) = Catch in IYR.
      FCATCH(IYR) = 0.D0

      RETURN
      END

C *********************************************************************

      SUBROUTINE RDPARS
C
      RETURN
      END
