subroutine excite_psi_bar_pert
  use main_mod,  only : nx, ny, nz, dx, dy, dz
  use tddft_mod
  implicit none
  integer ix,iy,iz
  real*8  rr(3)
  complex*16, parameter :: ci=(0d0,1d0)
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           rr = (/ (ix-1-nx/2)*dx, (iy-1-ny/2)*dy, (iz-1-nz/2)*dz /)
           phi_bar_tot_pert(ix,iy,iz,:,:) = phi_bar_tot(ix,iy,iz,:,:)* exp(-ci*sm*rr(ipol))
        enddo
     enddo
  enddo
end subroutine excite_psi_bar_pert


subroutine excite_psi_bar_pert_homo
  use main_mod,  only : nx, ny, nz, dx, dy, dz, nocc
  use tddft_mod
  implicit none
  integer ix,iy,iz
  real*8  rr(3)
  complex*16, parameter :: ci=(0d0,1d0)

  phi_bar_tot_pert = phi_bar_tot
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           rr = (/ (ix-1-nx/2)*dx, (iy-1-ny/2)*dy, (iz-1-nz/2)*dz /)
           phi_bar_tot_pert(ix,iy,iz,nocc,:) = phi_bar_tot(ix,iy,iz,nocc,:)* exp(-ci*sm*rr(ipol))
        enddo
     enddo
  enddo
end subroutine excite_psi_bar_pert_homo
