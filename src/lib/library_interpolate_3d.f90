!
! This routine does a 3d interpolation of a function,
! The grid used is with xmin,nx,dx,etc., and the output grid has
! rmin, dr, nr.
!
! The interpolation uses sinc.

subroutine interpolate_plot_3d(  &
     rho, &
     xmin, dx, nx,   &
     ymin, dy, ny,   &
     zmin, dz, nz,    &
     rmin, dr, nr,   &
     ar) 

  implicit none

  complex*16, parameter :: ci=(0.d0,1.d0)
  integer nx, ny, nz
  integer ix, iy, iz
  integer ikx, iky, ikz
  real*8 rho(nx, ny, nz)
  real*8     dx, dy, dz
  real*8     dkx, dky, dkz
  real*8      xx, yy, zz
  real*8      kx, ky, kz
  real*8      xmin, ymin, zmin

  integer ir, nr
  real*8  ar(nr), dr, rmin, rmax, kk, rr, tt
  complex*16 ck(nx, ny, nz)
  
  real*8 pi

  pi = dacos(-1.d0)

  if(nx<2.or.ny<2.or.nz<2) then
     write(6,*)' nx,ny,nz problem  ',nx,ny,nz 
     stop
  endif
  
  !
  ! ck(ikx,iky,ikz) =
  ! sum_ix_iy_iz &
  !  (exp(-ci*sum(kvec(:,ikx,iky,ikz)*rvec(:,ix,iy,iz) ))*  &
  !                                                      rho(ix,iy,iz)
  ! = 
  !  sum_ix_iy_iz
  !  (exp(-ci*sum(kvec(:,ikx,iky,ikz)*rvec(:,ix,iy,iz)-rvec(:,1,1,1)))*  &
  !                                                      rho(ix,iy,iz)
  !   * exp(-ci*sum(kvec(:,ikx,iky,ikz)*rvec(:,1,1,1)))

  ar = 0.d0
  ck = rho
  call fft_regular_3d(ck,Nx,Ny,Nz,-1)
  ck = ck / (nx*ny*nz)
  
  dkx = 2*pi/nx/dx
  dky = 2*pi/ny/dy
  dkz = 2*pi/nz/dz

  do ikx=1, nx
     do iky=1, ny
        do ikz=1, nz

           kx = (ikx-1)*dkx; if(ikx>nx/2+1) kx = kx - 2*pi/dx
           ky = (iky-1)*dky; if(iky>ny/2+1) ky = ky - 2*pi/dy
           kz = (ikz-1)*dkz; if(ikz>nz/2+1) kz = kz - 2*pi/dz
         
           ck(ikx, iky, ikz) = &
           ck(ikx, iky, ikz) * &
           exp(-ci*(kx*xmin+ky*ymin+kz*zmin))

           kk = sqrt( kx**2 + ky**2 +kz**2)

           do ir=1, nr
              rr = rmin+dr*(ir-1)              
              tt = kk*rr
              if(tt>1.e-8) then
                 ar(ir) = &
                 ar(ir) + sin(tt)/tt * ck(ikx, iky, ikz)
              else
                 ar(ir) = &
                 ar(ir) + ck(ikx, iky, ikz)
              end if
           end do
        end do
     end do
  end do

  write(6,*) ' sum(ck) ',sum(ck)
end subroutine interpolate_plot_3d

!
! comment out
!

!call interpolate_plot_3d_driver
! end

subroutine interpolate_plot_3d_driver

  implicit none
  integer, parameter ::  nx=8
  integer, parameter ::  ny=nx
  integer, parameter ::  nz=nx
  real*8 rho(nx, ny, nz)
  real*8, parameter :: dx = 2d0
  real*8, parameter :: dy = dx
  real*8, parameter :: dz = dx
  integer   ix, iy, iz, ir
  real*8      xx, yy, zz
  real*8      xmin, ymin, zmin
  
  integer, parameter ::  nr = 100

  real*8  ar(nr), rr, dr, rmin, rmax
  
  xmin = -(nx-1)*dx/2.d0
  ymin = -(ny-1)*dy/2.d0
  zmin = -(nz-1)*dz/2.d0

  rmin = 0.d0
  rmax = xmin+(nx-1)*dx
  dr = rmax/(nr-1)
  do ir=1, nr
     rr = (ir-1)*dr
  end do

  do ix=1, nx
     do iy=1, ny
        do iz=1, nz
           xx = xmin + (ix-1)*dx
           yy = ymin + (iy-1)*dy
           zz = zmin + (iz-1)*dz
           rr = sqrt(xx**2 + yy**2 + zz**2)
           rho(ix, iy, iz) = ff(rr)
        enddo
     enddo
  enddo
  
  call interpolate_plot_3d(  &
     rho, &
     xmin, dx, nx,   &
     ymin, dy, ny,   &
     zmin, dz, nz,   &
     rmin, dr, nr,   &
     ar) 

  do ir=1, nr
     rr = (ir-1)*dr
     write(6,*)rr,ar(ir),ff(rr)
  enddo

contains
  function ff(qq)
    implicit none
    real*8 ff, qq
    ff = 0.5d0*exp(-qq**2/2.d0/3.d0**2)
  end function ff
  
end subroutine interpolate_plot_3d_driver
