C     Last change:  CA   29 Apr 2006    1:56 am
C     MAN-V12 : Single Stock Control program
C     Program links MANTST11 modules with catches set by the CLC

      INCLUDE 'mantst12.for'
      INCLUDE 'manextrd.for'
      INCLUDE 'nrmp-inf.for'
      INCLUDE 'nrmp.for'


      SUBROUTINE EXTRA (NTRIAL, OPTRAN, ERATE, ENDSUR, IOUT,PROG,ITUNE,
     +                  OPTDEPL) 

C     Subroutine EXTRA is used to open output files.
C     and set global parameters.
C
      INTEGER NTRIAL, OPTRAN, ENDSUR,IOUT,ITUNE,OPTDEPL
      DOUBLE PRECISION ERATE
      CHARACTER PROG*8

C      OPEN (IOUT,FILE='RESTEST')
C      OPEN (1,FILE='MANLOG')
      PROG = 'MANV12'
      ITUNE = 0

C     EXTRA may also be used to avoid wasting time by setting the NTRIAL=1 
C     for trials when the catch quota is set to zero
C     IF (OPTDEPL.EQ.0.AND.OPTRAN.NE.1 .AND. ERATE.EQ.0.0) NTRIAL = 1
C     ENDSUR = 0

      RETURN
      END
