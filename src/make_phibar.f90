subroutine make_phi_bar
  use main_mod
  use tddft_mod
  use mpi_lib_ours
  implicit none
  complex*16 :: sp(nx,ny,nz,nb)
  real*8, allocatable :: overlap(:,:)
  integer :: i, j, nstates_r
  logical :: file_exists

!Block 1
  call scatter_c16(phit_tot,phit,size(phit),0)
  do i=1,nb
     call sn_phi(phit(:,:,:,i,1),sp(:,:,:,i),1,0.5d0)
  enddo
  call gather_c16(sp,phi_bar_tot,size(phit),0)

  if (rank==0) then
    inquire(file='overlap_psibar', exist=file_exists)
    if (file_exists) then
      open(unit=123, file='overlap_psibar', status='replace')
    else
      open(unit=123, file='overlap_psibar', status='new')
    endif
  !End Block 1

  !Block 2
    !phi_bar_tot = phit_tot
    !inquire(file='overlap_psi', exist=file_exists)
    !if (file_exists) then
    !  open(unit=123, file='overlap_psi', status='replace')
    !else
    !  open(unit=123, file='overlap_psi', status='new')
    !endif
  !End Block 2

    nstates_r = nb*nodes
    allocate(overlap(nstates_r, nstates_r))

    rewind(123)
    if (rank==0) then
    write(6,*) 'checking orthogonality'
    write(6,*) 'nstates', nstates_r
    do i=1,nstates_r
      do j=1,nstates_r
         overlap(i,j) = sum(conjg(phi_bar_tot(:,:,:,i,1))*phi_bar_tot(:,:,:,j,1))*dv
         if (i .eq. j) then
           if (abs(overlap(i,i)-1d0) > 1e-8) then
             write(6,*) 'warning: abs(<psi_bar_i|psi_bar_j>-1d0) larger than 1e-8'
             write(6,*) 'i, overlap(i,i)', i, overlap(i,i)
             !stop 
           endif
         else
           if (abs(overlap(i,j)) > 1e-8) then
             write(6,*) 'warning: <psi_bar_i|psi_bar_j> larger than 1e-8'
             write(6,*) 'i, j, overlap(i,j)', i, j, overlap(i,j)
             !stop
           endif
         endif
      enddo
      write(123,*) overlap(i,:)
    enddo
    endif

    deallocate(overlap)
  endif

end subroutine

