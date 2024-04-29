subroutine tddft_rpa_xc_prep
  use main_mod
  use tddft_mod
  implicit none
  integer :: st,i
  real*8, allocatable :: ncoret1d(:)
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

  if(funct==0) then
      funct_x = 1
      funct_c = 12
  else if(funct==1) then
      funct_x = 101
      funct_c = 130
  else
      write(*,*) 'funct should be 0 or 1'
      stop
  endif

  allocate(ncoret1d(nn),stat=st); if(st/=0) stop 'ncoret1d rpa_xc_prep'

  ncoret1d = reshape(ncoret,(/nn/))

  write(6,*) 'pre vxc_libxc dens, ncoret1d', maxval(dens), maxval(ncoret1d)
  call vxc_libxc(dens+ncoret1d,vxc(:,1),nn,1)
  call vxc_libxc(dens_pert+ncoret1d,vxc_pert(:,1),nn,1)
  deallocate(ncoret1d)
  write(6,*) 'sum(vxc) rpa_xc_prep', sum(vxc**3d0)

  if (h_type .eq. 1) then
    phi_bar= tmp
    phi_bar_pert = tmp_pert
    deallocate(tmp, tmp_pert, sp)
  endif
end subroutine
