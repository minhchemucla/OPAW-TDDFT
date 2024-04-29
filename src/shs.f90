!apply S^-1/2 H S^-1/2 to a wavefunction
!S^-1/2 (K+Vloc) S^-1/2 + S^-1/2 D S^-1/2
subroutine shs(pin,pout,ik)
    use mpi_lib_ours
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    implicit none

    complex*16,intent(in) :: pin(nx,ny,nz)
    complex*16 :: pout(nx,ny,nz)
    complex*16 :: tmp(nx,ny,nz)

    integer :: ib,jb,ia,is,js,it,ik
    integer :: ix,iy,iz,ix1,iy1,iz1,jx,jy,jz
    integer :: ixs,iys,izs,ixe,iye,ize

    call sn_phi(pin,tmp,ik,-0.5d0)
    call h_phi(tmp,pout,ik)
    tmp=pout
    call sn_phi(tmp,pout,ik,-0.5d0)
end subroutine shs      

subroutine shs_debug(pin,pout,ik)
    use mpi_lib_ours
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    implicit none

    complex*16,intent(in) :: pin(nx,ny,nz)
    complex*16 :: pout(nx,ny,nz)
    complex*16 :: tmp(nx,ny,nz)

    integer :: ib,jb,ia,is,js,it,ik
    integer :: ix,iy,iz,ix1,iy1,iz1,jx,jy,jz
    integer :: ixs,iys,izs,ixe,iye,ize

    call sn_phi(pin,tmp,ik,-0.5d0)
    write(6,*) 's', sum(abs(tmp)**3d0)
    call h_phi(tmp,pout,ik)
    write(6,*) 'sh', sum(abs(pout)**3d0)
    tmp=pout
    call sn_phi(tmp,pout,ik,-0.5d0)
    write(6,*) 'shs', sum(abs(tmp)**3d0)
end subroutine shs_debug

subroutine shs_pert(pin,pout,ik)
    use mpi_lib_ours
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    implicit none

    complex*16,intent(in) :: pin(nx,ny,nz)
    complex*16 :: pout(nx,ny,nz)
    complex*16 :: tmp(nx,ny,nz)

    integer :: ib,jb,ia,is,js,it,ik
    integer :: ix,iy,iz,ix1,iy1,iz1,jx,jy,jz
    integer :: ixs,iys,izs,ixe,iye,ize

    call sn_phi(pin,tmp,ik,-0.5d0)
    call h_phi_pert(tmp,pout,ik)
    tmp=pout
    call sn_phi(tmp,pout,ik,-0.5d0)
end subroutine shs_pert
