subroutine grad_fft(nx, ny, nz, dx, dy, dz, dens, grad)
  !
  ! inefficient now, later improve
  !
  implicit none
  integer       st
  integer       ix, iy, iz
  integer       nx, ny, nz
  real*8        dx, dy, dz
  real*8   dens(nx, ny, nz)
  real*8   grad(nx, ny, nz, 3)
  
  real*8,      allocatable :: kx(:), ky(:), kz(:)
  complex*16,  allocatable :: bk(:, :, :)
  complex*16,  allocatable :: ck(:, :, :)
  complex*16,  allocatable :: cr(:, :, :)
  complex*16,  parameter   :: ci=(0d0,1d0)

  allocate( kx(nx),         stat=st); call check0(st,' kx  ')
  allocate( ky(ny),         stat=st); call check0(st,' ky  ')
  allocate( kz(nz),         stat=st); call check0(st,' kz  ')
  allocate( bk(nx, ny, nz), stat=st); call check0(st,' bk  ')
  allocate( ck(nx, ny, nz), stat=st); call check0(st,' ck  ')
  allocate( cr(nx, ny, nz), stat=st); call check0(st,' cr  ')

  call make_k
  bk = dens
  call fft3d_general(bk, nx, ny, nz, 1)
  !
  ! grad_x
  !
  ck = bk
  do ix= 1, nx
     ck(ix, iy, iz) = ci* kx(ix) * ck(ix, iy, iz)
  enddo
  call fft3d_general(ck, nx, ny, nz, -1)
  call check_real(ck, nx*ny*nz)
  grad( :, :, :, 1) = ck

  !
  ! grad_y
  !
  ck = bk
  do iy= 1, ny
     ck(ix, iy, iz) = ci* ky(iy) * ck(ix, iy, iz)
  enddo
  call fft3d_general(ck, nx, ny, nz, -1)
  call check_real(ck, nx*ny*nz)
  grad( :, :, :, 2) = ck

  !
  ! grad_z
  !
  ck = bk
  do iz= 1, nz
     ck(ix, iy, iz) = ci* kz(iz) * ck(ix, iy, iz)
  enddo
  call fft3d_general(ck, nx, ny, nz, -1)
  call check_real(ck, nx*ny*nz)
  grad( :, :, :, 3) = ck
  !
  ! 
  !
  deallocate(bk, ck, cr)
  deallocate(kx, ky, kz)

contains
  subroutine make_k
    implicit none
    call make_kx
    call make_ky
    call make_kz
  end subroutine make_k

  subroutine make_kx
    implicit none
    real*8 dkx, pi
    pi = dacos(-1d0)
    dkx = 2d0*pi/dx

    do ix=1,nx
       if(ix.le.nx/2) then
          kx(ix) = (ix-1)*dkx
       elseif (ix==nx/2+1) then
          kx(ix) = 0d0
       else
          kx(ix) = (ix-1-nx)*dkx
       end if
    end do
  end subroutine make_kx

  subroutine make_ky
    implicit none
    real*8 dky, pi
    pi = dacos(-1d0)
    dky = 2d0*pi/dy

    do iy=1,ny
       if(iy.le.ny/2) then
          ky(iy) = (iy-1)*dky
       elseif (iy==ny/2+1) then
          ky(iy) = 0d0
       else
          ky(iy) = (iy-1-ny)*dky
       end if
    end do
  end subroutine make_ky

  subroutine make_kz
    implicit none
    real*8 dkz, pi
    pi = dacos(-1d0)
    dkz = 2d0*pi/dz

    do iz=1,nz
       if(iz.le.nz/2) then
          kz(iz) = (iz-1)*dkz
       elseif (iz==nz/2+1) then
          kz(iz) = 0d0
       else
          kz(iz) = (iz-1-nz)*dkz
       end if
    end do
  end subroutine make_kz
end subroutine grad_fft
