subroutine ncpaw_make_hamiltonian
    use main_mod, only : dens,iscf, nx, ny, nz, nn
    use main_mod, only : nk_loc, h_type, nb, phit_tot, phit
    use mpi_lib_ours
    implicit none
    integer :: st, i
    complex*16, allocatable :: tmp(:,:,:,:,:), sp(:,:,:,:)

    if(h_type .eq. 1) then
      !for hbar=shs, need to form phi_tilde=S^-1/2 phi_bar
      allocate(tmp(nx,ny,nz,nb*nodes, nk_loc), stat=st)
      if(st/=0) stop 'allocate phi_tilde'
      allocate(sp(nx,ny,nz,nb), stat=st)
      if(st/=0) stop 'allocate phi_tilde'

      tmp = phit_tot

      call scatter_c16(phit_tot,phit,size(phit),0)
      do i=1,nb
         call sn_phi(phit(:,:,:,i,1),sp(:,:,:,i),1,-0.5d0)
      enddo
      call gather_c16(sp,phit_tot,size(phit),0)
    endif


    if(iscf>1) call update_dens
    call get_pot(dens)
    call get_dij

    if(h_type .eq. 1) then
      phit_tot = tmp
      deallocate(tmp, sp)
    endif
end subroutine

subroutine ncpaw_make_hamiltonian_tddft_init
    use main_mod, only : dens,iscf, nx, ny, nz, nn
    use main_mod, only : nk_loc, h_type, nb, phit_tot, phit
    use tddft_mod, only  : nstates, nc_init, nb_init
    use mpi_lib_ours
    implicit none
    integer :: st, i, j
    complex*16, allocatable :: tmp(:,:,:,:,:)

    if(rank < nc_init .and. h_type .eq. 1) then
      !for hbar=shs, need to form phi_tilde=S^-1/2 phi_bar
      !allocate(tmp(nx,ny,nz,nstates, nk_loc), stat=st)
      
      !if(allocated(phit))deallocate(phit)
      !allocate(phit(nx,ny,nz,nb_init,nk_loc))
      do i=1,nc_init-1
        do j=1,nb_init
          if(rank==0)then
            phit(:,:,:,j,1) = phit_tot(:,:,:,i*nb_init+j,1)
          endif
          call send_receive_c16_array_in_mpi(phit(:,:,:,j,1), size(phit(:,:,:,j,1)), 0, i, i)
        enddo
      enddo
      if(rank==0) then
        do i=1,nb_init
          phit(:,:,:,i,1) = phit_tot(:,:,:,i,1)
        enddo
      endif

      allocate(tmp(nx,ny,nz,nb_init, nk_loc), stat=st)
      if(st/=0) stop 'allocate tmp '
      tmp = phit

      do i=1,nb_init
         !call sn_phi(tmp(:,:,:,i,1),phit_tot(:,:,:,i,1),1,-0.5d0)
         call sn_phi(tmp(:,:,:,i,1),phit(:,:,:,i,1),1,-0.5d0)
         !write(6,*) 'rank, i, snphit', rank, i, sum(abs(tmp(:,:,:,i,1)-phit(:,:,:,i,1)))
      enddo

      if(rank==0) then
        do i=1,nb_init
          phit_tot(:,:,:,i,1) = phit(:,:,:,i,1) 
        enddo
      endif
      do i=1,nc_init-1
        do j=1,nb_init
          call send_receive_c16_array_in_mpi(phit(:,:,:,j,1), size(phit(:,:,:,j,1)), i, 0, i)
          if(rank==0)then
             phit_tot(:,:,:,i*nb_init+j,1) = phit(:,:,:,j,1)
          endif
        enddo
      enddo
    endif

    call sync_mpi


    call update_dens
    call get_pot(dens)
    !vks will be different from full scf because of no DIIS
    call get_dij


    if(rank<nc_init .and. h_type .eq. 1) then
      if(rank==0) then
        do i=1,nb_init
          phit_tot(:,:,:,i,1) = tmp(:,:,:,i,1) 
        enddo
      endif
      do i=1,nc_init-1
        do j=1,nb_init
          call send_receive_c16_array_in_mpi(tmp(:,:,:,j,1), size(tmp(:,:,:,j,1)), i, 0, i)
          if(rank==0)then
             phit_tot(:,:,:,i*nb_init+j,1) = tmp(:,:,:,j,1)
          endif
        enddo
      enddo
      deallocate(tmp)
    endif
    call sync_mpi
end subroutine

subroutine ncpaw_make_hamiltonian_tddft
    use main_mod
    use tddft_mod
    use mpi_lib_ours
    use atom_mod, only : at => atominfo
    implicit none
    integer :: st, i
    complex*16, allocatable :: tmp(:,:,:,:,:), sp(:,:,:,:)
    complex*16, allocatable :: tmp_pert(:,:,:,:,:)

    if (h_type .eq. 1) then
      allocate(tmp(nx,ny,nz,nb, nk_loc), stat=st)
      if(st/=0) stop 'allocate tmp'
      allocate(tmp_pert(nx,ny,nz,nb, nk_loc), stat=st)
      if(st/=0) stop 'allocate tmp_pert'
      allocate(sp(nx,ny,nz,nb), stat=st)
      if(st/=0) stop 'allocate sp'

      tmp = phi_bar
      tmp_pert = phi_bar_pert
      do i=1,nb
         call sn_phi(tmp(:,:,:,i,1),phi_bar(:,:,:,i,1),1,-0.5d0)
         call sn_phi(tmp_pert(:,:,:,i,1),phi_bar_pert(:,:,:,i,1),1,-0.5d0)
      enddo
    endif

    !call tddft_rhoij
    call tddft_rhoij_2(.false.) !unpert
    call tddft_rhoij_2(.true.)  !pert
    !call tddft_dens 
    call tddft_dens_2(.false.)
    call tddft_dens_2(.true.)
    call get_nhat
    call get_nhat_pert
    !call get_pot_tddft(dens, dens_pert)
    call get_pot_tddft_2(.false.)
    call get_pot_tddft_2(.true.)
    call get_dij
    call get_dij_pert
    !call interpolate_vks_dij

    if (h_type .eq. 1) then
      phi_bar= tmp
      phi_bar_pert = tmp_pert
      deallocate(tmp, tmp_pert, sp)
    endif

    call tddft_dens_2(.false.)
    call tddft_dens_2(.true.)
end subroutine

subroutine ncpaw_make_hamiltonian_rpa
    use main_mod
    use tddft_mod
    use mpi_lib_ours
    use atom_mod, only : at=>atominfo
    implicit none
    integer :: st, i
    complex*16, allocatable :: tmp(:,:,:,:,:), sp(:,:,:,:)
    complex*16, allocatable :: tmp_pert(:,:,:,:,:)

    if (h_type .eq. 1) then
      allocate(tmp(nx,ny,nz,nb, nk_loc), stat=st)
      if(st/=0) stop 'allocate tmp'
      allocate(tmp_pert(nx,ny,nz,nb, nk_loc), stat=st)
      if(st/=0) stop 'allocate tmp_pert'
      allocate(sp(nx,ny,nz,nb), stat=st)
      if(st/=0) stop 'allocate sp'

      tmp = phi_bar
      tmp_pert = phi_bar_pert
      do i=1,nb
         call sn_phi(tmp(:,:,:,i,1),phi_bar(:,:,:,i,1),1,-0.5d0)
         call sn_phi(tmp_pert(:,:,:,i,1),phi_bar_pert(:,:,:,i,1),1,-0.5d0)
      enddo
    endif

    call tddft_rhoij
    call tddft_dens 
    call get_nhat
    call get_nhat_pert
    !call get_pot_tddft(dens, dens_pert)
    call get_vh_tddft(dens,dens_pert)
    call get_dij
    call get_dij_pert

    if (h_type .eq. 1) then
      phi_bar = tmp
      phi_bar_pert = tmp_pert
      deallocate(tmp, tmp_pert, sp)
    endif

    call tddft_dens 
end subroutine

subroutine ncpaw_make_hamiltonian_tddft_unpert
    use main_mod
    use tddft_mod
    use mpi_lib_ours
    use atom_mod, only : at => atominfo
    implicit none
    integer :: st, i
    complex*16, allocatable :: tmp(:,:,:,:,:), sp(:,:,:,:)

    if (h_type .eq. 1) then
      allocate(tmp(nx,ny,nz,nb, nk_loc), stat=st)
      if(st/=0) stop 'allocate tmp'
      allocate(sp(nx,ny,nz,nb), stat=st)
      if(st/=0) stop 'allocate sp'

      tmp = phi_bar
      do i=1,nb
         call sn_phi(tmp(:,:,:,i,1),phi_bar(:,:,:,i,1),1,-0.5d0)
      enddo
    endif

    call tddft_rhoij_2(.false.) !unpert
    call tddft_dens_2(.false.)
    call get_nhat
    call get_pot_tddft_2(.false.)
    call get_dij

    if (h_type .eq. 1) then
      phi_bar= tmp
      deallocate(tmp, sp)
    endif

    call tddft_dens_2(.false.)
end subroutine

subroutine ncpaw_make_hamiltonian_tddft_pert
    use main_mod
    use tddft_mod
    use mpi_lib_ours
    use atom_mod, only : at => atominfo
    implicit none
    integer :: st, i
    complex*16, allocatable :: sp(:,:,:,:)
    complex*16, allocatable :: tmp_pert(:,:,:,:,:)

    if (h_type .eq. 1) then
      allocate(tmp_pert(nx,ny,nz,nb, nk_loc), stat=st)
      if(st/=0) stop 'allocate tmp_pert'
      allocate(sp(nx,ny,nz,nb), stat=st)
      if(st/=0) stop 'allocate sp'

      tmp_pert = phi_bar_pert
      do i=1,nb
         call sn_phi(tmp_pert(:,:,:,i,1),phi_bar_pert(:,:,:,i,1),1,-0.5d0)
      enddo
    endif

    call tddft_rhoij_2(.true.)  !pert
    call tddft_dens_2(.true.)
    call get_nhat_pert
    call get_pot_tddft_2(.true.)
    call get_dij_pert

    if (h_type .eq. 1) then
      phi_bar_pert = tmp_pert
      deallocate(tmp_pert, sp)
    endif
    call tddft_dens_2(.true.)
end subroutine
