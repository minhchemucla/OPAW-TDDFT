subroutine propdt
  use main_mod
  use tddft_mod
  use mpi_lib_ours
  implicit none
  integer :: flg
  integer :: is,jb
  !real*8  :: wtime

  !call scatter_c16(phi_bar_tot,phi_bar,size(phi_bar),0)
  !call scatter_c16(phi_bar_tot_pert,phi_bar_pert,size(phi_bar_pert),0)
  !do is=1,nb
    !jb=state_map(is)
    !if(jb<=nocc) then

  if (h_type .eq. 1) then
    if (prop_type .eq. 1) then
      stop 'prop_type 1 debugging broke now'
      !call prop_bar_dt_debug_1(phi_bar(:,:,:,is,1))
      !call prop_bar_dt_debug_1_pert(phi_bar_pert(:,:,:,is,1))
    else if (prop_type .eq. 2) then
      !if(is==1) then
      !  call rk4_prop_shs_debug(phi_bar(:,:,:,is,1), phi_bar_pert(:,:,:,is,1))
      !else
        call rk4_prop_shs(phi_bar(:,:,:,:,1), phi_bar_pert(:,:,:,:,1), nn,nb)
      !endif
    else
      write(6,*) 'only implemented prop_type=1(debug) and 2(RK4)'
      stop
    endif
  else
    if (prop_type .eq. 1) then
      stop 'prop_type 1 debugging broke now'
      !call prop_dt_sh_debug(debug_prop_flag, phi_bar(:,:,:,is,1), phi_bar_pert(:,:,:,is,1))
    else if (prop_type .eq. 2) then
      stop 'need to fix prop_sh '
      !wtime = mpi_wtime()
      call rk4_prop_sh
      !wtime = mpi_wtime() - wtime
      !write(6,*) 'rank, wall time rk4', rank, wtime
    else 
      write(6,*) 'only implemented prop_type=1(debug) and 2(RK4)'
      stop
    endif
  endif
    !endif
  !enddo
  !call gather_c16(phi_bar,phi_bar_tot,size(phi_bar),0)
  !call gather_c16(phi_bar_pert,phi_bar_tot_pert,size(phi_bar_pert),0)
end subroutine
