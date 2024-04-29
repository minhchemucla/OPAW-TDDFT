      subroutine cg_driver
!      program cg_driver
      implicit none
      external fcn_rosen
      integer i;      integer, parameter :: n= 2
      real*8 :: x(n), toll
      toll = 1.d-8
      do i=1,n
         x(i) = -2.d+00
      end do
      call cg_top_sub(n, x, toll, fcn_rosen)
      end

      module cg_module
      implicit none
      save
      integer mp,lp
      end
c
c       --------------------------------------------------------------------
c       Main program for running the conjugate gradient methods described in 
c       the paper:
c
c       Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence Properties 
c       of Conjugate Gradient Methods", SIAM Journal on Optimization, Vol. 2,
c       pp. 21-42.
c
c       A web-based Server which solves unconstrained nonlinear optimization
c       problems using this Conjugate Gradient code can be found at:
c
c       http://www-neos.mcs.anl.gov/neos/solvers/UCO:CGPLUS/
c
c       Written by G. Liu, J. Nocedal and R. Waltz
c       October 1998
c       --------------------------------------------------------------------
c 
      subroutine CG_top_sub(n, x, eps, fcn_sub)
      use cg_module
      !
      ! n: dimension; x(n): vector;eps: tolerance
      ! fcn_sub: a user supplied function declared 
      ! external in calling program, of the form:
      !    call fcn_sub(n,x,f,g)
      ! where f is the function to be minimized, g is the gradient.    
      !
c
c     Change the maximum size of the problem dimension here
c     
      implicit none
      external fcn_sub
      integer n
      double precision x(n),g(n),d(n),gold(n),w(n)
      double precision f,eps,tlev
      double precision time1,time2,tottime
      
      logical          finish
      integer          iprint(2),iflag,icall,method,i
      integer          iter,nfun
      common /runinf/  iter,nfun
      real*8 one
      data one/1.0D+0/
      integer irest

      FINISH= .FALSE.
c     
c     Read problem input information
c     
      method =    3
      irest =     1   
      iprint(1) = 1 
      iprint(2) = 0  
c     
c     Check for correct dimension value n
c     
      if (n .lt. 0) then
         iflag = -3
         write(*,850)
         go to 50
      end if
      if (n .gt. n) then
         iflag = -3
         write(*,860)
         go to 50
      end if
c     
c     Print parameters
c     
      if (iprint(1) .ge. 0) then
         write (*,820)
         write (*,840) n, method, irest
      end if
      
      ICALL=0
c     
c     This is the convergence constant 
c     
      
c     IFLAG=0 indicates an initial entry to program
      
      IFLAG=0
c     
c     Begin counting CPU time. 
c     (Note: This function may not work on all operating systems.)
c     
      call timer_dn(time1)
      
 20   CONTINUE
c     
c     Calculate the function and gradient values here
c     
      call fcn_sub(n,x,f,g)
      
 30   CONTINUE
c     
c     Call the main optimization code
c     
      CALL CGFAM(N,X,F,G,D,GOLD,IPRINT,EPS,W,
     *     IFLAG,IREST,METHOD,FINISH )
c     
c     IFLAG=
c     0 : successful termination
c     1 : return to evaluate F and G
c     2 : return with a new iterate, try termination test
c     -i : error
c     
      IF(IFLAG.LE.0.OR.ICALL.GT.10000) GO TO 50
      IF(IFLAG.EQ.1) THEN
         ICALL=ICALL + 1   
         write(6,*)' function (in cgfam) ',f
         GO TO 20
      ENDIF 
      IF(IFLAG.EQ.2) THEN
c     
c     Termination Test.  The user may replace it by some other test. However, 
c     the parameter 'FINISH' must be set to 'TRUE' when the test is satisfied.
c     
         TLEV= EPS*(ONE + DABS(F))
         I=0
 40      I=I+1
         IF(I.GT.N) THEN
            FINISH = .TRUE.
            GO TO 30
         ENDIF
         IF(DABS(G(I)).GT.TLEV) THEN
            GO TO 30
         ELSE
            GO TO 40
         ENDIF
         
      ENDIF
      
 50   continue
c     
c     End CPU counting
c     
      call timer_dn(time2)
c     
c     Calculate the elapsed time
c     
      tottime = time2-time1
c     
c     Code has terminated; print final results
c     
      if (iprint(1).ge.0.and.iflag.ge.0) then
         write (*,890) f
         write (*,900) tottime
      end if
      write(6,*)' returning '
c     
c     Formatting
c     
 800  format (12x, i3)
 820  format (//, ' Conjugate Gradient Minimization Routine', /)
 840  format (/, ' n  =', i6, /,
     *     ' method   =', i6,/,
     *     ' irest    =', i6,/)
 850  format (/,'  Error: negative N value'/)
 860  format (/,'  Error: N too large, incrs param ndim'/)
 890  format (/,' f(x*) =', 1pd16.8)
 900  format (' It took ',1pd13.6,' CPU seconds'/)
      
      end
      























