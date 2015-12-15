C     Last change:  CA   25 Sep 2015    2:56 pm
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C     xrmpSub_inc.f
C
C     Include file for xrmpSub.f
C     Norwegian Computing Center, december 1992, Jon Helgeland
C     Modified, April 1999, Ragnar Bang Huseby
C
C
C     Error message constants
C     -----------------------
C     On exit of the subroutine CATCHLIMIT, IFAIL is equal to one of
C     these constants.
C
      integer OK,NUMERR,ABDIMERR,CATCHERR,AERR,VERR,
     *     ABYEARERR,NMAXERR,
     *     PPROBERR,MUERR,DTERR,BERR,PLEVELERR,ACCQUOTAERR,PTERR1,
     *     PTERR2,JACOBIERR,SETSPLITERR,HALFGAUSSERR,PRODGAUSSERR,LTERR,
     *     RIDGEERR,QUOTAERR,ACCERR,IFAILERR
      parameter(OK=0,NUMERR=2,ABDIMERR=3,CATCHERR=4,AERR=5,VERR=6,
     *     ABYEARERR=7,PPROBERR=8,
     *     MUERR=9,DTERR=10,BERR=11,PLEVELERR=12,ACCQUOTAERR=13,
     *     NMAXERR=14,PTERR1=15,PTERR2=16,JACOBIERR=17,SETSPLITERR=18,
     *     HALFGAUSSERR=19,PRODGAUSSERR=20,
     *     LTERR=21,RIDGEERR=22,QUOTAERR=23,ACCERR=24,IFAILERR=-2)


C     Variables used by Gauss-Legendre integration rules
C     --------------------------------------------------
C
C     nmax specifies the size of the arrays xg, wg, xg1, xg2, wg1, wg12,
C     and intSplit. nmax is a constant. nmax should be at least as large 
C     as the maximum number of abscissas used in the integration rules.
C
C     xg and wg contain the abscissas and the weights, respectively, 
C     in the approximation of the b-integral.
C
C     xg1 and wg1 contain the abscissas and the weights, respectively, 
C     in the approximation of the mu-integral.
C
C     xg2 contains the abscissas in the approximation of the p0-integral.
C
C     wg12 contains the product of the weights used in the approximation of
C     the mu-integral and p0-integrals.
C
C     intSplit contains p0,split(mu) for the abscissas used in the
C     approximation of the mu-integral.
C
C     ng1 is the number of points used in the approximation of
C     the mu-integral.
C     
C     ng2 is the number of points used in the approximation of
C     the p0-integral.
C     
C     ngauss is the number of points used in the approximation of
C     the b-integral.
C
      integer nmax
      parameter(nmax=1200)
      REAL(8) xg(nmax),wg(nmax)
      common /gauss1/ xg,wg
      integer ngauss
      common /gauss2/ ngauss

      REAL(8) xg1(nmax),xg2(nmax,2*nmax),wg1(nmax),
     *                  wg12(nmax,2*nmax),intSplit(nmax)
      common /pgauss1/ xg1,xg2,wg1,wg12,intSplit
      integer ng1,ng2
      common /pgauss2/ ng1,ng2





