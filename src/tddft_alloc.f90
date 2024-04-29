subroutine alloc_tddft
  use param_mod, only : k_flg, iscf_hminmax, flg_bin
  use libpaw_mod
  use main_mod
  use mpi_lib_ours
  use paw_mod
  use tddft_mod 
  use atom_mod
  implicit none
  integer :: st
  allocate(nhat_pert(nx,ny,nz), dens_pert(nn), dens0(nn))
  allocate(vks_pert(nx,ny,nz),vxc_pert(nn,1),vh_pert(nx,ny,nz))
  allocate(vks_pert_old(nx,ny,nz))

  allocate(phi_bar(nx,ny,nz,nb,nk_loc))
  allocate(phi_bar_pert(nx,ny,nz,nb,nk_loc))
  if(rank==0) then
    !allocate(phi_bar_tot(nx,ny,nz,nocc,nk_loc))
    !allocate(phi_bar_tot_pert(nx,ny,nz,nocc,nk_loc))
    allocate(phi_bar_tot(nx,ny,nz,nb*nodes,nk_loc))
    allocate(phi_bar_tot_pert(nx,ny,nz,nb*nodes,nk_loc))
    !for parallelization
  else
     !phit_tot, phi_bar_tot, and phi_bar_tot_pert are stored only on rank=0
    allocate(phi_bar_tot(0,0,0,0,0))
    allocate(phi_bar_tot_pert(0,0,0,0,0))
  endif
end subroutine

