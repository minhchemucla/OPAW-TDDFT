!      call drvr_dummy_t_dvr1
!      end

      subroutine drvr_dummy_t_dvr1
      implicit none
      integer i, j
      integer, parameter :: n=200
      real*8 dx, x_grd(n), wf(n), wf1(n), wf1_an(n), t_dvr_1drv
      real*8, parameter :: xmin=-3, xmax=5, xmid = (xmin+xmax)/2 + 0.1

      dx = (xmax-xmin)/(n+1)

      do i=1, n
         x_grd(i) = xmin+ i*dx
      enddo
      
      wf     = exp(-(x_grd-xmid)**2/0.3**2)
      wf1_an =   -2.d0/0.3**2 *(x_grd - xmid) * wf

      wf1 = 0.d0
      do i=1, n
         do j=1, n
            wf1(j) = wf1(j) + t_dvr_1drv(j-i)/dx*wf(i)
         enddo
      enddo

      do i=1, n
         write(6,88) i, x_grd(i), wf(i), wf1(i), wf1_an(i)
      enddo
 88   format(' ',i5,4f16.8)

      call r_1drv_dvr(wf,wf1,n,dx,1) 
      write(6,*)' new function '

      do i=1, n
         write(6,88) i, x_grd(i), wf(i), wf1(i), wf1_an(i)
      enddo

      end

      function t_dvr_2drv(j)
      implicit none
      integer j
      real*8 t_dvr_2drv
      real*8 pi

      if(j.eq.0) then
         pi = dacos(-1.d0)
         t_dvr_2drv = pi**2/3.d0
      else
         t_dvr_2drv = 2* (-1)**j /dble(abs(j))**2.d0
      endif

      end

c-----
c call this as follows:
c
c      wf1 = 0.d0
c      do i=1, n
c         do j=1, n
c            wf1(j) = wf1(j) + t_dvr_1drv(j-i)/dx*wf(i)
c         enddo
c      enddo
c
c-----
      function t_dvr_1drv(j)
      implicit none
      integer j
      real*8 t_dvr_1drv
      real*8 pi

      if(j==0) then 
         t_dvr_1drv = 0.d0
      else
         t_dvr_1drv = (-1)**j/dble(j)
      endif

      end

      subroutine r_1drv_dvr(phi,phi_1,ng,dx,lot) 
      implicit none
      integer ng, lot, j,k 

      real*8 phi(  ng, lot)
      real*8 phi_1(ng, lot)

      real*8 t_dvr_1drv, dx

      integer ifirst
      save ifirst
      data ifirst / 1 /

      integer,parameter :: ntop=2100
      real*8 tt(-ntop:ntop); save tt

      if(ng>ntop) then
         write(6,*)' increase ntop in l*t_dvr_f; ntop,ng ',ntop,ng
      endif

      if(ifirst==1) then
         ifirst = -1
         do j=-ntop, ntop
            tt(j) = t_dvr_1drv(j)
         enddo
      endif

      phi_1 = 0.d0
      do j=1, ng
         do k=1, ng
            phi_1(j,:) = phi_1(j,:) + tt(j-k)/dx*phi(k,:)
         enddo
      enddo
      end

      subroutine r_2drv_dvr(phi,phi_1,ng,dx,lot) 
      implicit none
      integer ng, lot, j,k 

      real*8 phi(  ng, lot)
      real*8 phi_1(ng, lot)

      real*8 t_dvr_2drv, dx

      integer ifirst
      save ifirst
      data ifirst / 1 /

      integer,parameter :: ntop=2100
      real*8 tt(-ntop:ntop); save tt

      if(ng>ntop) then
         write(6,*)' increase ntop in l*t_dvr_f; ntop,ng ',ntop,ng
      endif

      if(ifirst==1) then
         ifirst = -1
         do j=-ntop, ntop
            tt(j) = t_dvr_2drv(j)
         enddo
      endif

      phi_1 = 0.d0
      do j=1, ng
         do k=1, ng
            phi_1(j,:) = phi_1(j,:) + tt(j-k)/dx**2*phi(k,:)
         enddo
      enddo
      end


      subroutine c_1drv_dvr(phi,phi_1,ng,dx,lot) 
      implicit none
      integer ng, lot, j,k 

      complex*16 phi(  ng, lot)
      complex*16 phi_1(ng, lot)

      real*8 t_dvr_1drv, dx

      integer ifirst
      save ifirst
      data ifirst / 1 /

      integer,parameter :: ntop=2100
      real*8 tt(-ntop:ntop); save tt

      if(ng>ntop) then
         write(6,*)' increase ntop in l*t_dvr_f; ntop,ng ',ntop,ng
      endif

      if(ifirst==1) then
         ifirst = -1
         do j=-ntop, ntop
            tt(j) = t_dvr_1drv(j)
         enddo
      endif

      phi_1 = 0.d0
      do j=1, ng
         do k=1, ng
            phi_1(j,:) = phi_1(j,:) + tt(j-k)/dx*phi(k,:)
         enddo
      enddo
      end

      subroutine c_2drv_dvr(phi,phi_1,ng,dx,lot) 
      implicit none
      integer ng, lot, j,k 

      complex*16 phi(  ng, lot)
      complex*16 phi_1(ng, lot)

      real*8 t_dvr_2drv, dx

      integer ifirst
      save ifirst
      data ifirst / 1 /

      integer,parameter :: ntop=2100
      real*8 tt(-ntop:ntop); save tt

      if(ng>ntop) then
         write(6,*)' increase ntop in l*t_dvr_f; ntop,ng ',ntop,ng
      endif

      if(ifirst==1) then
         ifirst = -1
         do j=-ntop, ntop
            tt(j) = t_dvr_2drv(j)
         enddo
      endif

      phi_1 = 0.d0
      do j=1, ng
         do k=1, ng
            phi_1(j,:) = phi_1(j,:) + tt(j-k)/dx**2*phi(k,:)
         enddo
      enddo
      end


