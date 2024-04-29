      subroutine  get_ef(ef, evl, nx, nch, xm_x, vd, dx)
      implicit none
      integer nx, nch, nbas
      real*8 ef(nx, nch), evl(nch), xm_x, vd(nx), dx

      nbas = min(3*nch+8,nx)
      call get_ef_nbasin(ef, evl, nx, nch, nbas, xm_x, vd, dx)

      end


      subroutine  get_ef_nbasin(ef, evl, nx, nch, nbas, xm_x, vd, dx)
      implicit none
      integer nx, nch, i, j, l, ix, nbas

      
      real*8 ef(nx, nch), evl(nch), xm_x, vd(nx), dx,pi
      real*8 Lx_sin,  ovrlp, delta_l; external delta_l
      
      real*8, allocatable, dimension(:, :) :: bas, h_bas, h_evc
      real*8, dimension(:), allocatable :: h_evl

      write(6,*)'vd(1),vd(nx), nx ',vd(1),vd(nx),nx
      write(6,*)' dx,xm_x ',dx,xm_x
      write(6,*)' nch, nx ',nch,nx
      write(6,*)' min_vd, max_vd, min_abs_vd ',
     x  	minval(vd),maxval(vd),minval(abs(vd))
      

      if(nch>nx.or.nch<1.or.nx<1.or.nx<nbas.or.nbas<nch) then
         write(6,*)' nch, nx,nbas ',nch, nx,nbas
         stop
      endif
      
      allocate(bas(nx, nbas), h_bas(nbas, nbas), 
     x        h_evc(nbas, nbas), h_evl(nbas))
      
      if(Nx == 1) then
         ef(1,1) = 1; evl(1)=vd(1)
         return
      endif

c     mnake sin basis
      pi = dacos(-1.d0)
      Lx_sin = (nx+1)*dx
      
      do ix=1, nx
         do j=1, nbas
            bas(ix, j) = dsqrt(2.d0/(nx+1)/dx)*
     x           sin(ix*j*pi/(nx+1.d0))
         enddo
      enddo


      do j=1, nbas
         do l=1, nbas
            
            ovrlp = sum(bas(:,j)*bas(:,l)*dx)
            if(abs(ovrlp-delta_l(j,l))>1.d-8) then
               write(6,*)' ovrlp ',j,l,ovrlp
               stop
            endif
!!            write(10,*)' ovrlp ',ovrlp, j, l
         enddo
      enddo

c     make hamiltonian eigenvalues
      
      write(6,*)' step 152 '
      do j=1, nbas
         do l=1, nbas
            h_bas(j, l) = 0.d0
            if(j==l) h_bas(j, l) = j**2*pi**2/2/xm_x/Lx_sin**2
            
            h_bas(j, l) = h_bas(j, l) + 
     x           sum(vd(:)*bas(:,j)*bas(:, l))*dx
         enddo
      enddo
      
      call r_diag_normlz(h_bas, h_evc, h_evl, nbas)

      bas = matmul(bas, h_evc) ! convert to grd

      do i=1, nch
         evl(i) = h_evl(i)
         ef(:, i) = bas(:, i)  
      enddo

c     print, check
      
      do i =1, nch
         write(6,*)' i ',i , ' evl ', evl(i)
      enddo

      do j=1, nch
         do l=1, nch
            ovrlp = sum(ef(:,j)*ef(:,l)*dx)
            if(abs(ovrlp-delta_l(j,l))>1.d-8) then
               write(6,*)' ef, ovrlp ',j,l,ovrlp
               stop
            endif
!!            write(10,*)' ef, ovrlp',ovrlp, j, l
         enddo
      enddo
      deallocate(bas, h_bas, h_evc, h_evl)

      end

      function delta_l(i, j)
      implicit none
      real*8 delta_l
      integer i, j
      
      if( i == j) then
         delta_l = 1.d0
      else
         delta_l = 0.d0
      endif
      
      end
