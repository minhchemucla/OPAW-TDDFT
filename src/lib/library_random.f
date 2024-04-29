!      call driver_ran_normal
!      end
!
!     <x^2 ex(-x^2/2) dx = d/dbeta log(exp(-x^2 beta)) (beta=0.5)
!     =d/dbeta log(sqrt(1/beta))=1/2beta= 1
!
      subroutine driver_ran_normal
      implicit none
      real*8 sum_r, sum_r2, ran_normal, x
      integer i
     
      do i=1, 1000
         x =  ran_normal()
         sum_r  = sum_r + x
         sum_r2 = sum_r2 + x**2
      enddo
      
      write(6,*)' x, sum_r ',x, sum_r, sum_r2/1000.d0
      end

      function ran_normal_explct_seed(iseed)
      implicit none
      integer iseed
      real*8 box_muller, gasdev_explct, term, ran_normal_old
      real*8 ran_normal_explct_seed
      integer, parameter :: method_normal = 2
      select case (method_normal)
      case(2)
         term = gasdev_explct(iseed)
      case default
         write(6,*)' method_normal ',method_normal
         stop
      end select
      ran_normal_explct_seed = term
      end function

      function ran_normal()
      implicit none
      real*8 box_muller, gasdev, term, ran_normal_old, ran_normal
      integer, parameter :: method_normal = 2
      select case (method_normal)
      !case(1)
      !   term = box_muller()
      case(2)  ! only case adapated for mpi
         term = gasdev()
      !case(3)
      !   term = ran_normal_old()
      case default
         write(6,*)' method_normal ',method_normal
         stop
      end select
      ran_normal = term
      end function

      function box_muller()
      ! Box Muller, a-la
      ! http://www.taygeta.com/random/gaussian.html

      implicit none
      real*8 box_muller, ran_ps_neg
      real*8 x1, x2, w, y1, y2;
 
      w = 2.d0
      do while(w.ge.1.d0)
         x1 = ran_ps_neg()
         x2 = ran_ps_neg()
         w = x1 * x1 + x2 * x2;
      enddo

      w = sqrt( (-2.d00 * log( w ) ) / w );
      y1 = x1 * w;
      y2 = x2 * w;
      box_muller = y1
      end function


      function ran_normal_old()   ! returns normal a-la exp(-x**2/2.d0)
      implicit none
      real*8 ran_normal_old, x, ran_ps_neg
      integer irep, nrep

      nrep=100
      ran_normal_old = 0.d0

      do irep= 1, nrep
         ran_normal_old = ran_normal_old + ran_ps_neg()
      enddo

      ran_normal_old = ran_normal_old / dsqrt( nrep / 3.d0 )

      end

      function ran_ps_neg()   ! dist. equall between -1, 1
      implicit none
      real*8 ran_ps_neg, ran_ps
      ran_ps_neg = ran_ps()*2d0 - 1d0
      end

!+++++++!(99.10.5) instead of ran_normal, get from recipe and modified by dxn
      FUNCTION gasdev()
      implicit none
      INTEGER idum

      integer ifirst
      save idum, ifirst
       data ifirst / 1 /


      REAL*8 gasdev !ran2
      real*8 ran_ps; external ran_ps
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/

      if(ifirst==1) then
	 ifirst = -1
	 idum = -99
      endif

      if (iset.eq.0) then
! 1      v1=2.*ran2(idum)-1.d0
!        v2=2.*ran2(idum)-1.d0
 1       v1 = 2d0*ran_ps()-1d0
        v2 = 2d0*ran_ps()-1d0
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

!+++++++!(99.10.5) instead of ran_normal, get from recipe and modified by dxn
      FUNCTION gasdev_explct(iseed)
      implicit none
      INTEGER iseed

      REAL*8 gasdev_explct, ran2
!     USES ran2 here; use ran1 if not working.
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/

      if (iset.eq.0) then
1       v1=2.*ran2(iseed)-1.d0
        v2=2.*ran2(iseed)-1.d0
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev_explct=v2*fac
        iset=1
      else
        gasdev_explct=gset
        iset=0
      endif
      return
      END


      FUNCTION RAN2_notreally(IDUM)
      implicit none
      real*8   ran2_notreally
      integer, PARAMETER:: M=714025,IA=1366,IC=150889
      real*8, parameter :: RM=1.4005112d-6
      integer IR(97), iff, idum, iy, j
      save IR, iY
      save iff
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1)then; write(6,*)' j_r',j;stop;endif
      IY=IR(J)
      RAN2_notreally=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END 

      FUNCTION ran2(idum)
!      use mpi_lib_ours, only : rank
      implicit none;       save
      integer, parameter :: rank=0  !erase
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      real*8 ran2
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      write(6,*)' seed and rank in ran2 begining ',idum,rank
      stop ' should not use ran2 anymore, use ran_ps() '
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
 11    continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      write(6,*)' seed,rank ran2 at the end ',idum,rank,ran2
      return
      END

      subroutine ran2_sub(idum, ran_p)
!      use mpi_lib_ours, only : rank
      implicit none
      integer, parameter :: rank=0  !erase
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      real*8 ran_p
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      write(6,*)' seed and rank in ran2 begining ',idum,rank
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
 11    continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran_p=min(AM*iy,RNMX)
      write(6,*)' seed,rank ran_p at the end ',idum,rank,ran_p
      return
      END












