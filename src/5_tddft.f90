module tddft_mod
  implicit none
  save
 
  integer :: tddft_flag ! determines if tddft is being used or not
  integer :: nt !number of time steps
  integer :: it !iterator
  integer :: ipol !direction of perturbation
  integer :: nstates !total number of states in wf.txt
  integer :: dyn ! 1=static hamiltonian, 2=RPA/TDHF, 3=TDDFT
  integer :: prop_type=2! 1 = debug, 2 = rk4, 3 = split_op !LEGACY
  integer :: debug_prop_flag!LEGACY
  integer :: tbar_taylor_power = 50!LEGACY
  integer :: vbar_taylor_power = 50!LEGACY
  integer :: na_node ! divide number of atoms amongst nodes
  !integer :: nfrozen !number of frozen states
  !integer :: nprop !number of propagated states
  integer :: nvirtual !number of unoccupied states to include for time propagation 
  integer :: nvirtual_read
  integer :: n_restart = -1 !every nt step write out wfs so can restart programp

  integer :: nc_init
  integer :: nc_vbar
  integer :: nb_init
  integer :: it_restart
  integer ,allocatable :: state_map(:) !assigns which state to each rank

  real*8  :: dt
  real*8  :: sm
  !real*8  :: pfrozen !percent frozen
  !dens0 i.e no perturbation is defined in main_mod as dens
  !nhat0 i.e no perturbation is defined in main_mod as nhat
  real*8, allocatable :: dens0(:)
  real*8, allocatable :: dens_pert(:) !pert = pertubed
  real*8, allocatable :: nhat_pert(:,:,:) 
  real*8, allocatable :: vxc_pert(:,:),vks_pert(:,:,:),vh_pert(:,:,:)
  real*8, allocatable :: vks_pert_old(:,:,:)
  complex*16, allocatable :: phi_bar(:,:,:,:,:) !phibar=s^1/2 phit
  complex*16, allocatable :: phi_bar_tot(:,:,:,:,:) !phibar=s^1/2 phit
  complex*16, allocatable :: phi_bar_pert(:,:,:,:,:) !perturbed states
  complex*16, allocatable :: phi_bar_tot_pert(:,:,:,:,:) !perturbed states
  logical :: prop_bar!LEGACY
  logical :: tddft_stop_print = .false.!LEGACY
end module tddft_mod
