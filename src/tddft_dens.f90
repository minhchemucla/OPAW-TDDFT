subroutine tddft_dens
  use main_mod
  use tddft_mod 
  use mpi_lib_ours
  use tddft_mod, only : tddft_stop_print
  implicit none
  integer :: ik, jk, i, jb
  real*8  :: rhok(nx,ny,nz),rhok_pert(nx,ny,nz)

  ik = 1
  dens = 0d0
  dens_pert = 0d0
  rhok= 0d0
  rhok_pert = 0d0
  jk=(ik-1)*nodes+rank+1
  !do i=1,nocc
  !if (rank .le. nodes-2) then
  do i=1,nb
    jb = state_map(i)
    if(jb<=nocc) then
      rhok=rhok+abs(phi_bar(:,:,:,i,ik))**2
      rhok_pert=rhok_pert+abs(phi_bar_pert(:,:,:,i,ik))**2
    endif
  enddo
  !else
    !
  !  do i=1,nb-nvirtual
  !      rhok=rhok+abs(phi_bar(:,:,:,i,ik))**2
  !      rhok_pert=rhok_pert+abs(phi_bar_pert(:,:,:,i,ik))**2
  !  enddo
  !endif
  dens=dens+2d0*reshape(rhok*wk(jk), (/nn/))
  dens_pert=dens_pert+2d0*reshape(rhok_pert*wk(jk), (/nn/))
  call allsum_r8(dens,size(dens))
  call allsum_r8(dens_pert,size(dens_pert))
  !dens=dens+dens0
  !dens_pert=dens_pert+dens0
    
  !if(rank==0 .and. .not. tddft_stop_print) &
  !  write(6,*) 'sum(dens), sum(dens_pert)', sum(dens)*dv, sum(dens_pert)*dv

end subroutine

subroutine tddft_dens_2(pert)
  use main_mod
  use tddft_mod 
  use mpi_lib_ours
  use tddft_mod, only : tddft_stop_print
  implicit none
  integer :: ik, jk, i, jb
  real*8  :: rhok(nx,ny,nz),rhok_pert(nx,ny,nz)
  logical :: pert

  ik = 1
  if(.not. pert) then
    dens = 0d0
    rhok= 0d0
  else
    dens_pert = 0d0
    rhok_pert = 0d0
  endif
  jk=(ik-1)*nodes+rank+1
  !do i=1,nocc
  !if (rank .le. nodes-2) then
  do i=1,nb
    jb = state_map(i)
    if(jb<=nocc) then
      if(.not. pert) rhok=rhok+abs(phi_bar(:,:,:,i,ik))**2
      if(pert)       rhok_pert=rhok_pert+abs(phi_bar_pert(:,:,:,i,ik))**2
    endif
  enddo
  if(.not. pert) then
    dens=dens+2d0*reshape(rhok*wk(jk), (/nn/))
    call allsum_r8(dens,size(dens))
  else
    dens_pert=dens_pert+2d0*reshape(rhok_pert*wk(jk), (/nn/))
    call allsum_r8(dens_pert,size(dens_pert))
  endif
  !dens=dens+dens0
  !dens_pert=dens_pert+dens0
    
  !if(rank==0 .and. .not. tddft_stop_print) &
  !  write(6,*) 'sum(dens), sum(dens_pert)', sum(dens)*dv, sum(dens_pert)*dv

end subroutine
