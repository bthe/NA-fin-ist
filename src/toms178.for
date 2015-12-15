C     Last change:  CA   25 Sep 2015    2:57 pm
      function best_nearby ( delta, point, prevbest, nvars, f, 
     &  funevals )

c*********************************************************************72
c
cc BEST_NEARBY looks for a better nearby point, one coordinate at a time.
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    The ALGOL original is by Arthur Kaupe.
c    C version by Mark Johnson
c    FORTRAN77 version by John Burkardt
c
c  Reference:
c
c    M Bell, Malcolm Pike,
c    Remark on Algorithm 178: Direct Search,
c    Communications of the ACM,
c    Volume 9, Number 9, September 1966, page 684.
c
c    Robert Hooke, Terry Jeeves,
c    Direct Search Solution of Numerical and Statistical Problems,
c    Journal of the ACM,
c    Volume 8, Number 2, April 1961, pages 212-229.
c
c    Arthur Kaupe,
c    Algorithm 178:
c    Direct Search,
c    Communications of the ACM,
c    Volume 6, Number 6, June 1963, page 313.
c
c    FK Tomlin, LB Smith,
c    Remark on Algorithm 178: Direct Search,
c    Communications of the ACM,
c    Volume 12, Number 11, November 1969, page 637-638.
c
c  Parameters:
c
c    Input, double precision DELTA(NVARS), the size of a step in each direction.
c
c    Input/output, double precision POINT(NVARS); on input, the current candidate.
c    On output, the value of POINT may have been updated.
c
c    Input, double precision PREVBEST, the minimum value of the function seen
c    so far.
c
c    Input, integer NVARS, the number of variables.
c    NVARS should not be more than 20, or whatever the internal
c    value of NVARS_MAX has been set to in HOOKE and BEST_NEARBY.
c
c    Input, external double precision F, the name of the function routine,
c    which should have the form:
c      function f ( x, n )
c      integer n
c      double precision f
c      double precision x(n)
c
c    Input/output, integer FUNEVALS, the number of function evaluations.
c
c    Output, double precision BEST_NEARBY, the minimum value of the function seen
c    after checking the nearby neighbors.
c
      implicit none

      integer nvars
      integer nvars_max
      parameter ( nvars_max = 50 )

      double precision best_nearby
      double precision delta(nvars)
      double precision f
      external f
      double precision ftmp
      integer funevals
      integer i
      double precision minf
      double precision point(nvars)
      double precision prevbest
      double precision z(nvars_max)

      minf = prevbest

      do i = 1, nvars
       z(i) = point(i)
      end do

      do i = 1, nvars

        z(i) = point(i) + delta(i)

        ftmp = f ( z, nvars )
        funevals = funevals + 1

        if ( ftmp .lt. minf ) then

          minf = ftmp

        else

          delta(i) = - delta(i)
          z(i) = point(i) + delta(i)
          ftmp = f ( z, nvars )
          funevals = funevals + 1

          if ( ftmp .lt. minf ) then
            minf = ftmp
          else
            z(i) = point(i)
          end if

        end if

      end do

      do i = 1, nvars
        point(i) = z(i)
      end do

      best_nearby = minf

      return
      end
      function hooke ( nvars, startpt, endpt, rho, eps, 
     &  itermax, f )   

c*********************************************************************72
c
c  Discussion:
c
c    This routine find a point X where the nonlinear objective function 
c    F(X) has a local minimum.  X is an N-vector and F(X) is a scalar.  
c    The objective function F(X) is not required to be differentiable
c    or even continuous.  The program does not use or require derivatives 
c    of the objective function. 
c
c    The user supplies three things: 
c    1) a subroutine that computes F(X), 
c    2) an initial "starting guess" of the minimum point X, 
c    3) values for the algorithm convergence parameters.  
c
c    The program searches for a local minimum, beginning from the    
c    starting guess, using the Direct Search algorithm of Hooke and  
c    Jeeves.
c
c    This program is adapted from the Algol pseudocode found in the
c    paper by Kaupe, and includes improvements suggested by Bell and Pike,
c    and by Tomlin and Smith.
c
c    The algorithm works by taking "steps" from one estimate of
c    a minimum, to another (hopefully better) estimate.  Taking 
c    big steps gets to the minimum more quickly, at the risk of 
c    "stepping right over" an excellent point.  The stepsize is 
c    controlled by a user supplied parameter called RHO.  At each 
c    iteration, the stepsize is multiplied by RHO  (0 < RHO < 1), 
c    so the stepsize is successively reduced. 
c
c    Small values of rho correspond to big stepsize changes, 
c    which make the algorithm run more quickly.  However, there 
c    is a chance (especially with highly nonlinear functions) 
c    that these big changes will accidentally overlook a 
c    promising search vector, leading to nonconvergence. 
c
c    Large values of RHO correspond to small stepsize changes, 
c    which force the algorithm to carefully examine nearby points 
c    instead of optimistically forging ahead.  This improves the 
c    probability of convergence. 
c
c    The stepsize is reduced until it is equal to (or smaller 
c    than) EPS.  So the number of iterations performed by 
c    Hooke-Jeeves is determined by RHO and EPS:
c 
c      RHO^(number_of_iterations) = EPS
c
c    In general it is a good idea to set RHO to an aggressively 
c    small value like 0.5 (hoping for fast convergence).  Then, 
c    if the user suspects that the reported minimum is incorrect 
c    (or perhaps not accurate enough), the program can be run 
c    again with a larger value of RHO such as 0.85, using the 
c    result of the first minimization as the starting guess to 
c    begin the second minimization.
c
c    Normal use: 
c    (1) Code your function F() in the C language;
c    (2) Install your starting guess;
c    (3) Run the program.
c
c    If there are doubts about the result, the computed minimizer 
c    can be used as the starting point for a second minimization attempt.
c
c    To apply this method to data fitting, code your function F() to be 
c    the sum of the squares of the errors (differences) between the 
c    computed values and the measured values.  Then minimize F() 
c    using Hooke-Jeeves. 
c
c    For example, you have 20 datapoints (T(i), Y(i)) and you want to
c    find A, B and C so that:
c
c      A*t*t + B*exp(t) + C*tan(t)
c
c    fits the data as closely as possible.  Then the objective function
c    F() to be minimized is just
c
c      F(A,B,C) = sum ( 1 <= i <= 20 )
c        ( y(i) - A*t(i)*t(i) - B*exp(t(i)) - C*tan(t(i)) )^2.
c
c    Several internal arrays, both in HOOKE and in BEST_NEARBY, are
c    dimensioned by the parameter NVARS_MAX, which is set to 20.
c    The user cannot run a problem for which NVARS is greater than
c    NVARS_MAX.  If NVARS is greater than the current value of
c    NVARS_MAX, the user must alter the program source code so that
c    NVARS_MAX is large enough, and recompile.
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    The ALGOL original is by Arthur Kaupe.
c    The C version was by Mark Johnson.
c    The FORTRAN77 version is by John Burkardt.
c
c  Reference:
c
c    M Bell, Malcolm Pike,
c    Remark on Algorithm 178: Direct Search,
c    Communications of the ACM,
c    Volume 9, Number 9, September 1966, page 684.
c
c    Robert Hooke, Terry Jeeves,
c    Direct Search Solution of Numerical and Statistical Problems,
c    Journal of the ACM,
c    Volume 8, Number 2, April 1961, pages 212-229.
c
c    Arthur Kaupe,
c    Algorithm 178:
c    Direct Search,
c    Communications of the ACM,
c    Volume 6, Number 6, June 1963, page 313.
c
c    FK Tomlin, LB Smith,
c    Remark on Algorithm 178: Direct Search,
c    Communications of the ACM,
c    Volume 12, Number 11, November 1969, page 637-638.
c
c  Parameters:
c
c    Input, integer NVARS, the number of spatial dimensions.
c    NVARS should not be more than 20, or whatever the internal
c    value of NVARS_MAX has been set to in HOOKE and BEST_NEARBY.
c
c    Input, double precision STARTPT(NVARS), the user-supplied
c    initial estimate for the minimizer. 
c
c    Output, double precision ENDPT(NVARS), the estimate for the
c    minimizer, as calculated by the program.
c
c    Input, double precision RHO, a user-supplied convergence parameter
c    which should be set to a value between 0.0 and 1.0.  Larger values 
c    of RHO give greater probability of convergence on highly nonlinear 
c    functions, at a cost of more function evaluations.  Smaller     
c    values of RHO reduce the number of evaluations and the program 
c    running time, but increases the risk of nonconvergence. 
c
c    Input, double precision EPS, the criterion for halting   
c    the search for a minimum.  When the algorithm   
c    begins to make less and less progress on each   
c    iteration, it checks the halting criterion: if  
c    the stepsize is below EPS, terminate the    
c    iteration and return the current best estimate  
c    of the minimum.  Larger values of EPS (such 
c    as 1.0e-4) give quicker running time, but a     
c    less accurate estimate of the minimum.  Smaller 
c    values of EPS (such as 1.0e-7) give longer  
c    running time, but a more accurate estimate of   
c    the minimum.            
c
c    Input, integer ITERMAX, a limit on the number of iterations.
c
c    Input, external double precision F, the name of the function routine,
c    which should have the form:
c      function f ( x, n )
c      integer n
c      double precision f
c      double precision x(n)
c
c    Output, integer HOOKE, the number of iterations taken.
c
      implicit none

      integer nvars

      integer nvars_max
      parameter ( nvars_max = 50 )

      double precision best_nearby
      double precision delta(nvars_max)
      double precision endpt(nvars)
      double precision eps
      double precision f
      external f
      double precision fbefore
      integer funevals
      integer hooke
      integer i
      integer itermax
      integer iters
      integer j
      integer keep
      double precision newf
      double precision newx(nvars_max)
      double precision rho
      double precision startpt(nvars)
      double precision steplength
      double precision tmp
      logical verbose
      parameter ( verbose = .false. )
      double precision xbefore(nvars_max)

      if ( nvars_max .lt. nvars ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HOOKE - Fatal error!'
        write ( *, '(i8,a,i8)' ) 
     &  nvars_max, ' =  NVARS_MAX < NVARS = ', nvars
        stop
      end if

      do i = 1, nvars
        newx(i) = startpt(i)
      end do

      do i = 1, nvars
        xbefore(i) = startpt(i)
      end do

      do i = 1, nvars
        if ( startpt(i) .eq. 0.0D+00 ) then
          delta(i) = rho
        else
          delta(i) = rho * abs ( startpt(i) )
        end if
      end do

      funevals = 0
      steplength = rho
      iters = 0
      fbefore = f ( newx, nvars )
      funevals = funevals + 1
      newf = fbefore

10    continue

      if ( iters .lt. itermax .and. eps .lt. steplength ) then

        iters = iters + 1

        if ( verbose ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8,a,g14.6)' ) 
     &    '  FUNEVALS, = ', funevals, '  F(X) = ', fbefore

          do j = 1, nvars
            write ( *, '(2x,i8,2x,g14.6)' ) j, xbefore(j)
          end do

        end if
c
c  Find best new point, one coordinate at a time.
c
        do i = 1, nvars
          newx(i) = xbefore(i)
        end do

        newf = best_nearby ( delta, newx, fbefore, nvars, f, 
     &    funevals )
        if (funevals .GT. itermax) go to 99
c
c  If we made some improvements, pursue that direction.
c
        keep = 1

20      continue

        if ( newf .lt. fbefore .and. keep .eq. 1 ) then

          do i = 1, nvars
c
c  Arrange the sign of DELTA.
c
            if ( newx(i) .le. xbefore(i) ) then
              delta(i) = - abs ( delta(i) )
            else
              delta(i) = abs ( delta(i) )
            end if
c
c  Now, move further in this direction.
c
            tmp = xbefore(i)
            xbefore(i) = newx(i)
            newx(i) = newx(i) + newx(i) - tmp
          end do

          fbefore = newf
          newf = best_nearby ( delta, newx, fbefore, nvars, f, 
     &      funevals )
          if (funevals .GT. itermax) go to 99
c
c  If the further (optimistic) move was bad...
c
          if ( fbefore .le. newf ) then
            go to 30
          end if
c
c  Make sure that the differences between the new and the old points 
c  are due to actual displacements; beware of roundoff errors that 
c  might cause NEWF < FBEFORE.
c
          keep = 0

          do i = 1, nvars
            if ( 0.5D+00 * abs ( delta(i) ) .lt.
     &        abs ( newx(i) - xbefore(i) ) ) then
              keep = 1
              go to 20
            end if
          end do

          go to 20

        end if

30      continue

        if ( eps .le. steplength .and. fbefore .le. newf ) then
          steplength = steplength * rho
          do i = 1, nvars
            delta(i) = delta(i) * rho
          end do
        end if

        go to 10

      end if
99    continue
      do i = 1, nvars
        endpt(i) = xbefore(i)
      end do

      hooke = iters

      return
      end
C
C ===========================================================================
C
      SUBROUTINE FIT22(X,SS,NDIM)

C     Set up the parameters for a fit

C     IMPLICIT NONE
      REAL(8) P(101,100),Y(101),X(100),SS,FUNK,TOL,GRD
      INTEGER(4) NDIM,I,J,ITER
      EXTERNAL FUNK

C     SET UP TOLERANCES AND GRIDDING
      TOL = 0.0001
      GRD = 1.2

      DO 1 I=1,NDIM+1
         DO 2 J=1,NDIM
            P(I,J)=X(J)
            IF((I-1).EQ.J) P(I,J) = GRD*P(I,J)
2        CONTINUE
1     CONTINUE

      DO 3 I=1,NDIM+1
         DO 4 J=1,NDIM
            X(J)=P(I,J)
4        CONTINUE
         Y(I)=FUNK(X)
3     CONTINUE

      CALL AMOEBA(P,Y,101,100,NDIM,TOL,ITER)

      DO 5 J = 1,NDIM
         X(J) = P(1,J)
5     CONTINUE
      SS = Y(1)

      RETURN
      END

      SUBROUTINE AMOEBA(P,Y,MP,NP,NDIM,FTOL,ITER)
C
C     MULTIDIMENSIONAL MINIMISATION OF THE FUNCTION FUNK(X) WHERE X IS
C     AN NDIM-DIMENSIONAL VECTOR, BY THE DOWNHILL SIMPLEX METHOD OF
C     NELDER AND MEAD. INPUT IS A MATRIX P WHOSE NDIM+1 ROWS ARE THE
C     NDIM-DIMENSIONAL VECTORS WHICH ARE THE VERTICES OF THE STARTING
C     SIMPLEX. [LOGICAL DIMENSIONS OF P ARE P(NDIM+1,NDIM); PHYSICAL
C     DIMENSIONS ARE INPUT AS P(MP,NP).] ALSO INPUT IS THE VECTOR Y
C     OF LENGTH NDIM+1, WHOSE COMPONENTS MUST BE PRE-INITIALISED TO
C     THE VALUES OF FUNK EVALUATED AT THE NDIM+1 VERTICES (ROWS) OF P;
C     AND FTOL IS THE FRACTIONAL CONVERGENCE TOLERANCE TO BE ACHIEVED
C     IN THE FUNCTION VALUE (N.B.!). ON OUTPUT, P AND Y WILL HAVE BEEN
C     RESET TO NDIM+1 NEW POINTS ALL WITHIN FTOL OF A MINIMUM FUNCTION
C     VALUE, AND ITER GIVES THE NUMBER OF ITERATIONS TAKEN.
C
C     FROM: NUMERICAL RECIPES - THE ART OF SCIENTIFIC COMPUTING
C           BY W. H. PRESS ET AL, CAMBRIDGE UNIVERSITY PRESS
C           ISBN 0-251-30811-9
C
C     ********************************************************************
C

C     SPECIFY THE MAXIMUM NUMBER OF DIMENSIONS, THREE PARAMETERS WHICH
C     DEFINE THE EXPANSIONS AND CONTRACTIONS, AND THE MAXIMUM NUMBER OF
C     ITERATIONS ALLOWED

C     IMPLICIT NONE
      REAL(8) ALPHA,BETA,GAMMA,EPS
      INTEGER(4) ITMAX
      PARAMETER (ALPHA=1.0,BETA=0.5,GAMMA=2.0,ITMAX=500,EPS=1.0E-10)

C     Global Data
      INTEGER(4) MP,NP,NDIM,ITER
      REAL(8) P(MP,NP),Y(MP),PR(100),PRR(100),PBAR(100),FTOL,FUNK
      EXTERNAL FUNK

C     Local Data
      REAL(8) YPR,YPRR
      INTEGER(4) I,ILO,IHI,INHI,J,MPTS

C     NOTE THAT MP IS THE PHYSICAL DIMENSION CORRESPONDING TO THE LOGICAL
C     DIMENSION MPTS, NP TO NDIM.

      MPTS = NDIM+1
      ITER=0

C     FIRST WE MUST DETERMINE WHICH POINT IS THE HIGHEST (WORST), NEXT
C     HIGHEST, AND LOWEST (BEST).

1     ILO=1
      IF(Y(1).GT.Y(2)) THEN
         IHI=1
         INHI=2
      ELSE
         IHI=2
         INHI=1
      ENDIF
      DO 11 I=1,MPTS
         IF(Y(I).LT.Y(ILO)) ILO=I
         IF(Y(I).GT.Y(IHI)) THEN
            INHI=IHI
            IHI=I
         ELSE IF(Y(I).GT.Y(INHI)) THEN
            IF(I.NE.IHI) INHI=I
         ENDIF
11    CONTINUE

C     COMPUTE THE FRACTIONAL RANGE FROM THE HIGHEST TO THE LOWEST AND
C     RETURN IF SATISFACTORY

      IF(2.0*ABS(Y(IHI)-Y(ILO)).LT.FTOL*(ABS(Y(IHI))+ABS(Y(ILO))+EPS)) 
     +   RETURN
      IF(ITER.EQ.ITMAX) THEN
c        WRITE(6,200)
200      FORMAT(1H ,'AMOEBA EXCEEDING MAXIMUM ITERATIONS')
         RETURN
      ENDIF
      ITER=ITER+1

      DO 12 J=1,NDIM
         PBAR(J)=0.
12    CONTINUE

C     BEGIN A NEW ITERATION. COMPUTE THE VECTOR AVERAGE OF ALL POINTS
C     EXCEPT THE HIGHEST, I.E. THE CENTRE OF THE "FACE" OF THE SIMPLEX
C     ACROSS FROM THE HIGH POINT. WE WILL SUBSEQUENTLY EXPLORE ALONG
C     THE RAY FROM THE HIGH POINT THROUGH THE CENTRE.

      DO 14 I=1,MPTS
         IF(I.NE.IHI) THEN
            DO 13 J=1,NDIM
               PBAR(J)=PBAR(J)+P(I,J)
13          CONTINUE
         ENDIF
14    CONTINUE

C     EXTRAPOLATE BY A FACTOR ALPHA THROUGH THE FACE, I.E. REFLECT THE
C     SIMPLEX FROM THE HIGH POINT

      DO 15 J=1,NDIM
         PBAR(J)=PBAR(J)/NDIM
         PR(J)=(1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15    CONTINUE

C     EVALUATE THE FUNCTION AT THE REFLECTED POINT

      YPR=FUNK(PR)

C     GIVES A BETTER RESULT THAN THE BEST POINT, SO TRY AN ADDITIONAL
C     EXTRAPOLATION BY A FACTOR GAMMA

      IF(YPR.LE.Y(ILO)) THEN
         DO 16 J=1,NDIM
            PRR(J)=GAMMA*PR(J)+(1.-GAMMA)*PBAR(J)
16       CONTINUE

C        CHECK THE FUNCTION THERE

         YPRR=FUNK(PRR)

C        THE ADDITIONAL EXTRAPOLATION SUCCEEDED, AND REPLACES THE
C        HIGHEST POINT

         IF(YPRR.LT.Y(ILO)) THEN
            DO 17 J=1,NDIM
               P(IHI,J)=PRR(J)
17          CONTINUE
            Y(IHI)=YPRR
         ELSE

C        THE ADDITIONAL EXTRAPOLATION FAILED, BUT WE CAN STILL USE THE
C        REFLECTED POINT

         DO 18 J=1,NDIM
            P(IHI,J)=PR(J)
18       CONTINUE
         Y(IHI)=YPR
      ENDIF

C     THE REFLECTED POINT IS WORSE THAN THE SECOND HIGHEST

      ELSE IF(YPR.GE.Y(INHI)) THEN

C        IF IT'S BETTER THAN THE HIGHEST, THEN REPLACE THE HIGHEST

         IF(YPR.LT.Y(IHI)) THEN
            DO 19 J=1,NDIM
               P(IHI,J)=PR(J)
19          CONTINUE
            Y(IHI)=YPR
         ENDIF

C        BUT LOOK FOR AN INTERMEDIATE LOWER POINT; IN OTHER WORDS
C        PERFORM A CONTRACTION OF THE SIMPLEX ALONG ONE DIMENSION
C        AND THEN EVALUATE THE FUNCTION

         DO 21 J=1,NDIM
            PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
21       CONTINUE
         YPRR=FUNK(PRR)

C        CONTRACTION GIVES AN IMPROVEMENT, SO ACCEPT IT

         IF(YPRR.LT.Y(IHI)) THEN
            DO 22 J=1,NDIM
               P(IHI,J)=PRR(J)
22          CONTINUE
            Y(IHI)=YPRR
         ELSE

C           CAN'T SEEM TO GET RID OF THAT HIGH POINT. BETTER CONTRACT
C           AROUND THE LOWEST (BEST) POINT

            DO 24 I=1,MPTS
               IF(I.NE.ILO) THEN
                  DO 23 J=1,NDIM
                     PR(J)=0.5*(P(I,J)+P(ILO,J))
                     P(I,J)=PR(J)
23                CONTINUE
                  Y(I)=FUNK(PR)
               ENDIF
24          CONTINUE
         ENDIF

      ELSE

C        WE ARRIVE HERE IF THE ORIGINAL REFLECTION GIVES A MIDDLING
C        POINT. REPLACE THE OLD HIGH POINT AND CONTINUE

         DO 25 J=1,NDIM
            P(IHI,J)=PR(J)
25       CONTINUE
         Y(IHI)=YPR

      ENDIF

      GO TO 1

      END

C -----------------------------------------------------------------------------


