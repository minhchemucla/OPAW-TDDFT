subroutine ek_prep
  use main_mod
  use tddft_mod, only : dt, tddft_flag
  use paw_mod
  use mpi_lib_ours
  implicit none
  integer i
  integer ix, iy, iz, ik, jk
  real*8  kx, ky, kz, kx1,ky1,kz1
  real*8  dkx, dky, dkz
  real*8  pi
  complex*16, parameter :: ci = (0d0,1d0)

  allocate(ek(nx,ny,nz,nk_loc))
  if(tddft_flag > -1) allocate(expk(nx,ny,nz,nk_loc))
  pi = dacos(-1d0)
  
  dkx = 2d0*pi/(nx*dx)
  dky = 2d0*pi/(ny*dy)
  dkz = 2d0*pi/(nz*dz)

  do iz=1,nz
     do iy=1,ny
        do ix=1,nx

           kx = (ix-1)*dkx
           ky = (iy-1)*dky
           kz = (iz-1)*dkz
           
           if(kx>pi/dx) kx = kx - 2d0*pi/dx
           if(ky>pi/dy) ky = ky - 2d0*pi/dy
           if(kz>pi/dz) kz = kz - 2d0*pi/dz

           do ik=1,nk_loc
                jk=(ik-1)*nodes+rank+1
                kx1 = kx+kpt(1,jk)
                ky1 = ky+kpt(2,jk)
                kz1 = kz+kpt(3,jk)

                ek(ix,iy,iz,ik) = (kx1**2+ky1**2+kz1**2)/2d0

                if(ek(ix,iy,iz,ik)>ekcut) then
                    ek(ix,iy,iz,ik)=ekcut !0d0
                endif
           enddo
        enddo
     enddo
  enddo
  if(tddft_flag > -1) expk=exp(-ci*dt*ek)/dble(nn)
  ek=ek/(dble(nx)*dble(ny)*dble(nz))
end subroutine ek_prep
