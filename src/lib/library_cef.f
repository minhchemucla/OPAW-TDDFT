      subroutine  c_get_cef(cef, cevl, nx, nch, xm_x, cvd, dx)
      implicit none
      integer nx, nch, i, j, l, ix, nbas
      complex*16 cef(nx, nch), cevl(nch), cvd(nx)
      real*8 xm_x,  dx,pi
      real*8 Lx_sin,  ovrlp, delta_k; external delta_k
      
      complex*16, allocatable, dimension(:, :) :: cbas  ! dont use real
      complex*16, allocatable, dimension(:, :) :: h_cevc, h_cbas
      complex*16, dimension(:), allocatable    :: h_cevl

      write(6,*)' dx,xm_x ',dx,xm_x
      write(6,*)' nch, nx ',nch,nx
      write(6,*)' cvd ',cvd
      
      if(nch>nx.or.nch<1.or.nx<1) then
         write(6,*)' nch, nx ',nch, nx
         stop
      endif

      nbas = nx
      
      allocate(cbas(nx, nbas), h_cbas(nbas, nbas), 
     x        h_cevc(nbas, nbas), h_cevl(nbas), stat=j);if(j/=0)stop
      
      if(Nx == 1) then
         cef(1,1) = 1; cevl(1)=cvd(1)
         return
      endif

c     mnake sin cbasis
      pi = dacos(-1.d0)
      Lx_sin = (nx+1)*dx
      
      do ix=1, nx
         do j=1, nbas
            cbas(ix, j) = dsqrt(2.d0/(nx+1)/dx)*
     x           sin(ix*j*pi/(nx+1.d0))
         enddo
      enddo

      do j=1, nbas
         do l=1, nbas
            ovrlp = sum(cbas(:,j)*cbas(:,l)*dx)
            if(abs(ovrlp-delta_k(j,l))>1.d-8) then
               write(6,*)' ovrlp ',j,l,ovrlp
               stop
            endif
         enddo
      enddo

c     make hamiltonian eigenvalues
      
      do j=1, nbas
         do l=1, nbas
            h_cbas(j, l) = 0.d0
            if(j==l) h_cbas(j, l) = j**2*pi**2/2/xm_x/Lx_sin**2
            
            h_cbas(j, l) = h_cbas(j, l) + 
     x           sum(cvd(:)*cbas(:,j)*cbas(:, l))*dx
         enddo
      enddo

      call cs_diag_normlz(h_cbas, h_cevc, h_cevl, nbas, 1, 1)  ! flags on evc,
                ! and method


      cbas = matmul(cbas, h_cevc) ! convert to grd
      
      do i=1, nch
         cevl(i) = h_cevl(i)
         cef(:, i) = cbas(:, i)  
      enddo
      
c     print, check
      
      do i =1, nch
         write(6,*)' i ',i , ' cevl(c_get) ', cevl(i)
      enddo

      do j=1, nch
         do l=1, nch
            ovrlp = sum(cef(:,j)*cef(:,l)*dx)
            if(abs(ovrlp-delta_k(j,l))>1.d-5) then
               write(6,*)' cef, ovrlp ',j,l,ovrlp
               stop
            endif
            write(10,*)' cef, ovrlp',ovrlp, j, l
         enddo
      enddo


      deallocate(cbas, h_cbas , h_cevl, h_cevc)
      end

      function delta_k(i, j)
      implicit none
      real*8 delta_k
      integer i, j
      
      if( i == j) then
         delta_k = 1.d0
      else
         delta_k = 0.d0
      endif
      
      end
