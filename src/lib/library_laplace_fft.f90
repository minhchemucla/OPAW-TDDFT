subroutine laplace_fft(nx, ny, nz, dx, dy, dz, dens, lap)
  !
  ! inefficient now, later improve
  !
  implicit none
  integer       st
  integer       ix, iy, iz
  integer       nx, ny, nz
  real*8        dx, dy, dz
  real*8   dens(nx, ny, nz)
  real*8    lap(nx, ny, nz )
  
  real*8,      allocatable :: k2(:, :, :)
  complex*16,  allocatable :: ck(:, :, :)

  allocate( k2(nx, ny, nz), stat=st); call check0(st,' k2  ')
  allocate( ck(nx, ny, nz), stat=st); call check0(st,' ck  ')

  call make_k2

  ck = dens
  call fft3d_general(ck, nx, ny, nz, 1)

  do ix= 1, nx
     ck(ix, iy, iz) = k2(ix, iy, iz) * ck(ix, iy, iz)
  enddo

  call fft3d_general(ck, nx, ny, nz, -1)

  call check_real(ck, nx*ny*nz)
  lap = ck

  deallocate(ck, k2)
contains
  subroutine make_k2
    implicit none
    real*8 dkx, dky, dkz, pi
    real*8 kx, ky, kz

    pi = dacos(-1d0)

    dkx = 2d0*pi/dx
    dky = 2d0*pi/dy
    dkz = 2d0*pi/dz

    do iz=1, nz
       do iy=1, ny
          do ix=1, nx

             kx = (ix-1)*dkx
             if(ix.gt.nx/2+1) kx = kx - dkx* nx 
             
             ky = (iy-1)*dky
             if(iy.gt.ny/2+1) ky = ky - dky* ny 
             
             kz = (iz-1)*dkz
             if(iz.gt.nz/2+1) kz = kx - dkz* nz 
             
             k2( ix, iy, iz) = kx**2 + ky**2 + kz**2 
          enddo
       enddo
    enddo
  end subroutine make_k2

end subroutine laplace_fft
