module electrostatic_param   ! erase it and repalce by your favorite module that holds nx,ny,nz
  implicit none              ! and dx,dy,dz
  save
  integer, parameter :: nx=32,ny=nx,nz=nx
  real*8,  parameter :: dx=2d0, dy=dx,dz=dx
end module electrostatic_param

subroutine rho_to_pot(rho, pot)
  use electrostatic_param
  implicit none

  real*8 pi
  real*8        rho(nx,ny,nz)
  real*8        pot(nx,ny,nz)
  complex*16  ctemp(nx,ny,nz)
  real*8  :: k2_inv(nx,ny,nz)
  save k2_inv
  integer :: ifirst=1 
  save    :: ifirst

  pi = dacos(-1.d0)

  if(ifirst==1) then
     call set_k2_inv(nx,ny,nz,dx,dy,dz,k2_inv)
     ifirst=-1
  end if

  ctemp = rho
  call fftk_3d(ctemp,1, Nx,Ny,Nz,1)
  
  ctemp = ctemp/(nx*ny*nz)

  ctemp = 4*pi*k2_inv * ctemp
 
  call fftk_3d(ctemp, 1, Nx,Ny,Nz,-1)
  call check_real(ctemp, size(ctemp))
  
  pot = ctemp

end subroutine rho_to_pot

subroutine set_k2_inv(nx,ny,nz,dx,dy,dz,k2_inv)
  implicit none
  integer ::  nx,ny,nz
  real*8  :: dx,dy,dz,k2_inv(nx,ny,nz)

  integer              :: ix, iy, iz
  real*8               :: kxmin, kymin, kzmin
  real*8               :: kx, ky, kz
  real*8               :: dkx, dky, dkz
  real*8               :: k2
  real*8               :: pi
  real*8, parameter    :: zerothresh = 1d-8

  pi = dacos(-1.d0)
  
  kxmin = 0d0
  kymin = 0d0
  kzmin = 0d0
  
  dkx = 2d0*pi/(nx*dx)
  dky = 2d0*pi/(ny*dy)
  dkz = 2d0*pi/(nz*dz)
  
  do iz = 1, nz
     do iy = 1, ny
        do ix = 1, nx
           
           ! reciprocal space grid:
           kx = kxmin + (ix-1)*dkx
           ky = kymin + (iy-1)*dky
           kz = kzmin + (iz-1)*dkz
           
           if(kx>pi/dx)     kx = kx - 2d0*pi/dx
           if(ky>pi/dy)     ky = ky - 2d0*pi/dy
           if(kz>pi/dz)     kz = kz - 2d0*pi/dz
           k2 = kx**2 + ky**2 + kz**2
           
           ! we dont want to divide by 0
           if (k2 > zerothresh) then
              k2_inv(ix,iy,iz) = 1d0 / k2
           else
              k2_inv(ix,iy,iz) = 0d0
           end if
           
        end do
     end do
  end do
  
end subroutine set_k2_inv
  
