subroutine calc_exx_homo
  use main_mod
  use mpi_lib_ours
  implicit none
  integer :: i,st
  real*8  :: exx,exx1, exx2, tmp2(nn), pot(nn), pot3d(nx,ny,nz)
  complex*16  :: psi_i(nx,ny,nz), psi_j(nx,ny,nz), nhatij(nx,ny,nz)
  complex*16, allocatable :: tmp(:,:,:,:,:), sp(:,:,:,:)
  

  call read_wfs

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

  if(rank==0) psi_i = phit_tot(:,:,:,nocc,1)
  if(rank==0) call exx_expect(psi_i,psi_i,phit_tot(:,:,:,1:nocc,1),nx,ny,nz,nocc,exx)

  if(h_type .eq. 1) then
    phit_tot = tmp
    deallocate(tmp)
  endif

  call sync_mpi
  stop 'exx debug'
  contains 
    subroutine read_wfs
      implicit none

      if(rank==0) then
        close(413)
        if(h_type.eq.0) then
          open(unit=413,file='wf.txt',status='old')
        else if(h_type.eq.1) then
          open(unit=413,file='wf_bar.txt',status='old')
        endif
        do i=1,11;read(413,*);enddo
        do i=1,nocc
          read(413,*)
          read(413,*) tmp2
          phit_tot(:,:,:,i,1) = reshape(tmp2,(/nx,ny,nz/))
        enddo
      endif

      call update_dens
    end subroutine read_wfs

end subroutine
