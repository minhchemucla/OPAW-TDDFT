subroutine opaw_tddft
  use param_mod, only : k_flg, iscf_hminmax, flg_bin
  use libpaw_mod
  use main_mod
  use mpi_lib_ours
  use paw_mod
  use tddft_mod 
  use atom_mod
  implicit none
  integer :: is, st, i,nds, istart
  !real*8 :: eig(nb*nodes)
  real*8, allocatable :: eig(:)
  complex*16, allocatable :: tmp(:,:,:,:,:), sp(:,:,:,:)
  complex*16, allocatable :: tmp_pert(:,:,:,:,:)
  real*8 :: wtime, wtime2, wtime3
  logical fex

  !call test_T
  !call test_V
  !call test_debug_sgw_prop

  nds = max(nodes,0)

  fex=.false.
  if(n_restart>0) inquire(file='./tddft_save',exist=fex)
  if(fex) then
    call read_phi_bar
    istart = it_restart + 1
  else
    if(tddft_flag .eq. 1) then
      call read_wf_tddft
    else
      nstates = nb*nds
    endif
    !call calc_nvirtual !unecessary 
    call tddft_get_nb
    call alloc_tddft
    if(tddft_flag .eq. 1) call tddft_single_scf
    !call write_wf_for_sgw_debug(eig,nstates)
    if(rank==0) phi_bar_tot = phit_tot(:,:,:,1:nocc+nvirtual,:)
    !if(rank==0) phi_bar_tot = phit_tot(:,:,:,nfrozen+1:nocc,:)
    if(rank==0) call excite_psi_bar_pert
    !if(rank==0) call excite_psi_bar_pert_homo 
    !call pert_transform
    !call prep_frozen
    call sync_mpi

    !wtime = mpi_wtime()
    !call prep_tbar_sh
    !call prep_vbar_sh
    !wtime = mpi_wtime() - wtime
    !if(rank==0) write(6,*) 'prep_vbar_time', wtime
    !call finalize_mpi_lib
    !stop

    if(dyn==1)call prep_static_H !static

    !call scatter_c16(phi_bar_tot,phi_bar,size(phi_bar),0)
    !call scatter_c16(phi_bar_tot_pert,phi_bar_pert,size(phi_bar_pert),0)
    !tddft_stop_print= .true.
    call scatter_phi_bar
    if(dyn.eq.2) call tddft_rpa_xc_prep
    istart = 1
    it_restart = 0
  endif
  !if(rank==0) write(6,*) 'sum(phi_bar_tot) pre prop', sum(abs(phi_bar_tot(:,:,:,1:2,1)))
  do it=istart,nt
    wtime = mpi_wtime()
    select case (dyn)
      case (1) !STATIC
        call tddft_dens
      case (2) !RPA
        call ncpaw_make_hamiltonian_rpa
      case (3) !TDDFT
        call ncpaw_make_hamiltonian_tddft
    end select
    wtime2 = mpi_wtime() - wtime
    if(rank==0) write(800,*) 'it, ham_time', it, wtime2
    wtime2 = mpi_wtime() 
    call propdt
    wtime2 = mpi_wtime() - wtime2
    if(rank==0) write(800,*) 'it, prop time', it, wtime2
    wtime2 = mpi_wtime() 
    call plot_dip
    wtime2 = mpi_wtime() - wtime2
    if(rank==0) write(800,*) 'it, plot time', it, wtime2
    if(n_restart>0 .and. mod(it,n_restart) .eq. 0) call write_to_restart
    wtime = mpi_wtime() - wtime
    if(rank==0) write(800,*) 'it, tot time', it, wtime
    !call gather_phi_bar
  enddo
contains
  subroutine tddft_single_scf
    implicit none
    integer :: i, j
    integer :: a, b, c
    complex*16 :: nrm
    !real*8,allocatable  :: eigs(:)
    !perform a single SCF cycle with as many cores as possible
    !because there is no DIIS the results will differ slightly from continuing an SCF
    iscf_hminmax = -1
    iscf = 1
    do nc_init=nodes,1,-1
      if(mod(nstates, nc_init)==0) exit
    enddo
    nb_init = nstates/nc_init
    if(rank==0)write(6,*) 'using ', nc_init,' cores for initial SCF with ', nb_init, ' states per core'

    if(allocated(phit))deallocate(phit)
    allocate(phit(nx,ny,nz,nb_init,nk_loc))
    call ncpaw_make_hamiltonian_tddft_init
    call set_hminmax

    close(441)
    if(h_type .eq. 0) then
      if(flg_bin)then; open (unit=441,file='wf.bin',status='old',form='unformatted');
      else;           open (unit=441,file='wf.txt',status='old');  endif
    else if (h_type .eq. 1) then
      if(flg_bin)then; open (unit=441,file='wf_bar.bin',status='old',form='unformatted');
      else;           open (unit=441,file='wf_bar.txt',status='old'); endif
    endif

    rewind(441)

    do i=1,9
     if(flg_bin)then; read(441)
     else;            read(441,*); endif
    enddo

    allocate(eig(nstates))
    if(flg_bin) then; read(441) eig(:)
    else; read(441,*) eig(:); endif
    mu=eig(nocc+1)+dmua
    !if(rank==0) write(6,*) 'eigs ', eig
    !if(rank==0) write(6,*) 'mu   ', mu
    !if(rank==0) write(6,*) 'dmua',  dmua
    !if(rank==0) write(6,*) 'nocc ', nocc
    !if(rank==0) write(6,*) 'eig(nocc+1) ', eig(nocc+1)

    !write(*,*) 'rank,hmax,mu',rank,hmax,mu
    a = nn
    b = nstates
    c = b
    !Euler's algorithm for GCF
    do while (mod(a,b) .ne. 0) 
      c = b
      b = mod(a,b)
      a = c
    enddo
    if(rank==0) write(6,*) 'Greatest common factor of nstates,nn', c
    if (c>nodes) then
      do nc_init=nodes,1,-1
        if(mod(b,nc_init).eq.0) exit
      enddo
    else 
      nc_init = c
    endif
    if(mod(nstates,nc_init).ne.0) stop 'there is a bug in calculating nc_init 1'
    if(mod(nn,nc_init).ne.0) stop 'there is a bug in calculating nc_init 2'
    if(rank==0) write(6,*) 'Greatest common factor of nstates,nn <= nodes', nc_init
    if(rank==0) write(6,*) 'nb_init', nb_init
    nb_init = nstates/nc_init
    
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
    call sync_mpi
    call chebyr_m_tddft_init(hmax,mu)
    call orthog_phi_tddft_init
    call diagh_tddft_init
    !vks will be different from full scf because of no DIIS thus eigenvalues will not be exactly the same
    call sync_mpi
  end subroutine

  subroutine test_completeness
    use atom_mod, only : p=> pawinfo, at=> atominfo
    implicit none
    complex*16, allocatable :: test_wf(:), final_wf(:)
    real*8, allocatable :: print_wf(:), p3d(:,:,:), p1d(:)
    complex*16, allocatable :: test_wf3d(:,:,:), final_wf3d(:,:,:)
    integer :: i, ii, ia, t
    integer :: ixs,iys,izs,ixe,iye,ize,jx,jy,jz
    integer :: ix,iy,iz

    close(139)
    open(unit=139, file='complete_diff')
    allocate(test_wf(nn), final_wf(nn), print_wf(nn))
    allocate(test_wf3d(nx,ny,nz), final_wf3d(nx,ny,nz))
    allocate(p3d(nx,ny,nz), p1d(nn))
    

    p3d = 0d0
    do ia=1,natom
      it=atom_map(ia)
      ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
      izs=at(ia)%ir_start(3)
      ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
      ize=izs+p(it)%nrough(3)-1
      if(stat/=0) stop 'ca alloc problem in sphi'
        if(.not.at(ia)%edge) then
          do is=1,p(it)%mstates
            p3d(ixs:ixe,iys:iye,izs:ize)=&
             p3d(ixs:ixe,iys:iye,izs:ize)+&
             at(ia)%local_p3d1_c(:,:,:,is,1)
          enddo
        else
          do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
            jx=mod(ix+nx-1,nx)+1
            jy=mod(iy+ny-1,ny)+1
            jz=mod(iz+nz-1,nz)+1
            p3d(jx,jy,jz)=p3d(jx,jy,jz)+&
              sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,1))
          enddo;enddo;enddo
        endif
    enddo

    test_wf3d = phit_tot(:,:,:,nb,1)

    call proj_phi(test_wf3d, final_wf3d, 1)

    test_wf = reshape(test_wf3d, (/nn/))
    final_wf = reshape(final_wf3d, (/nn/))
    print_wf = real(final_wf-test_wf)
    p1d = reshape(p3d, (/nn/))
    rewind(139)
    do i=1,nn
      write(139,*) print_wf(i), real(final_wf(i)), real(test_wf(i)), p1d(i)
    enddo
    call flush(139)
    !write(139,*) maxval(real(final_wf-test_wf))
    close(139)
    deallocate(test_wf, final_wf)

  end subroutine

  subroutine scatter_phi_bar
    implicit none
    integer :: i,j

    do i=1,nodes-1
      do j=1,nb
        if (rank==0) then
          phi_bar(:,:,:,j,:) = phi_bar_tot(:,:,:,i*nb+j,:)
          phi_bar_pert(:,:,:,j,:) = phi_bar_tot_pert(:,:,:,i*nb+j,:)
        endif
      enddo
      call send_receive_c16_array_in_mpi(phi_bar, size(phi_bar), 0, i, i)
      call send_receive_c16_array_in_mpi(phi_bar_pert, size(phi_bar_pert), 0, i, i)
    enddo
    if (rank==0) then
      do i=1,nb
        phi_bar(:,:,:,i,:) = phi_bar_tot(:,:,:,i,:)
        phi_bar_pert(:,:,:,i,:) = phi_bar_tot_pert(:,:,:,i,:)
      enddo
    endif
    call sync_mpi
  end subroutine

  subroutine gather_phi_bar
    implicit none
    integer :: i,j

    if(rank==0) then
      do i=1,nb
        phi_bar_tot(:,:,:,i,1)=phi_bar(:,:,:,i,1)
        phi_bar_tot_pert(:,:,:,i,1)=phi_bar_pert(:,:,:,i,1)
      enddo
    endif
    do i=1,nc_init-1
      do j=1,nb
        call send_receive_c16_array_in_mpi(phi_bar(:,:,:,j,1), size(phi_bar(:,:,:,j,1)), i, 0, i)
        call send_receive_c16_array_in_mpi(phi_bar_pert(:,:,:,j,1), size(phi_bar_pert(:,:,:,j,1)), i, 0, i)
        if(rank==0)then
          phi_bar_tot(:,:,:,i*nb+j,1)=phi_bar(:,:,:,j,1)
          phi_bar_tot_pert(:,:,:,i*nb+j,1)=phi_bar_pert(:,:,:,j,1)
        endif
      enddo
    enddo
    call sync_mpi
  end subroutine

  subroutine prep_static_H
    use atom_mod, only : at => atominfo
    implicit none
    integer :: ia

    allocate(tmp(nx,ny,nz,nocc, nk_loc), stat=st)
    if(st/=0) stop 'allocate phi_tilde'
    allocate(tmp_pert(nx,ny,nz,nocc, nk_loc), stat=st)
    if(st/=0) stop 'allocate phi_tilde'
    allocate(sp(nx,ny,nz,nb), stat=st)
    if(st/=0) stop 'allocate phi_tilde'

    tmp = phi_bar_tot

    if (h_type .eq. 1) call make_psi_tilde

    call tddft_rhoij
    do ia=1,natom
      at(ia)%rhoij_pert =  at(ia)%rhoij
    enddo
    call tddft_dens 
    call get_pot_tddft(dens, dens_pert)
    call get_nhat
    vh_pert=vh
    vxc_pert=vxc
    vks_pert=vks
    nhat_pert=nhat
    !  call get_nhat_pert

    call get_dij
    !call get_dij_pert
    do ia=1,natom
      at(ia)%dij_pert =  at(ia)%dij
    enddo
    
    if (h_type .eq. 1) call close_psi_tilde
  end subroutine

  subroutine make_psi_tilde
    implicit none
    integer :: i


    call scatter_c16(phi_bar_tot,phi_bar,size(phi_bar),0)
    do i=1,nb
       call sn_phi(phi_bar(:,:,:,i,1),sp(:,:,:,i),1,-0.5d0)
    enddo
    call gather_c16(sp,phi_bar_tot,size(phi_bar),0)

    tmp_pert = phi_bar_tot_pert

    call scatter_c16(phi_bar_tot_pert,phi_bar_pert,size(phi_bar_pert),0)
    do i=1,nb
       call sn_phi(phi_bar_pert(:,:,:,i,1),sp(:,:,:,i),1,-0.5d0)
    enddo
    call gather_c16(sp,phi_bar_tot_pert,size(phi_bar_pert),0)

  end subroutine

  subroutine close_psi_tilde
    implicit none
    phi_bar_tot = tmp
    phi_bar_tot_pert = tmp_pert

  end subroutine
end subroutine opaw_tddft
