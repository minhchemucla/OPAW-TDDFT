subroutine read_wf_debug
  use main_mod
  use tddft_mod
  use param_mod
  use mpi_lib_ours

  implicit none
  
  character*9 ch
  integer :: i,j, ix,iy,iz
  integer :: nx_r, ny_r, nz_r
  integer :: nstates_r
  integer :: i_orb, nsp_r
  real*8  :: dx_r, dy_r, dz_r, norm
  real*8, allocatable :: overlap(:,:)
  real*8 :: tmp(nx*ny*nz)
   complex*16 :: sp(nx,ny,nz,nb)

  if(k_flg) stop 'not printing wf for kpoint'


  close(441)
  !if(flg_bin)then; open (441,file='wfbar.bin',status='old',form='unformatted');
  !else;           open (441,file='wfbar.txt',status='old');  
  if(flg_bin)then; open (441,file='wf_debug.bin',status='old',form='unformatted');
  else;           open (441,file='wf_debug.txt',status='old');  
  end if
  rewind(441)
  
  if(flg_bin) then
     read(441) ch, nx_r
     read(441) ch, ny_r
     read(441) ch, nz_r
     read(441) ch, dx_r
     read(441) ch, dy_r
     read(441) ch, dz_r
     read(441) ch, nsp_r
     read(441) ch, nstates_r
  else
     read(441,*) ch, nx_r
     read(441,*) ch, ny_r
     read(441,*) ch, nz_r
     read(441,*) ch, dx_r
     read(441,*) ch, dy_r
     read(441,*) ch, dz_r
     read(441,*) ch, nsp_r
     read(441,*) ch, nstates_r
  endif

  call check(nx,nx_r,' nx, nx_r')
  call check(ny,ny_r,' ny, ny_r')
  call check(nz,nz_r,' nz, nz_r')
  call check_r(dx,dx_r,' dx, dx_r')
  call check_r(dy,dy_r,' dy, dy_r')
  call check_r(dz,dz_r,' dz, dz_r')
  call check(1,nsp_r,' nsp, nsp_r')
  call check(nb*nodes,nstates_r,' nb*nodes,nstates_r')

  read(441,*) ch !skip eigenvalues
  read(441,*) !skip eigenvalues
  read(441,*) ch !skip orbitals

  
  !debug
  !close(331)
  !open(unit=331, file='wf_debug.dat')
  
  if (rank==0) then
  do i=1,nstates_r
    read(441,*) i_orb, nsp_r
    read(441,*) tmp
    !norm = sqrt(sum(tmp*tmp)*dv)
    !write(6,*) norm
    phit_tot(:,:,:,i,1) = reshape(tmp,(/nx,ny,nz/))!/norm
    !phi_bar_tot(:,:,:,i,1) = reshape(tmp,(/nx,ny,nz/))
      !debug
      !write(331,*) i_orb, nsp_r
      !write(331, *) real(phit_tot(:,:,:,i,1))
  enddo
  endif

  !call scatter_c16(phit_tot,phit,size(phit),0)
  !do i=1,nb
  !   call sn_phi(phit(:,:,:,i,1),sp(:,:,:,i),1,0.5d0)
  !enddo
  !call gather_c16(sp,phi_bar_tot,size(phit),0)

  !!call scatter_c16(phi_bar_tot,phi_bar,size(phi_bar),0)
  !!do i=1,nb
  !!   call sn_phi(phi_bar(:,:,:,i,1),sp(:,:,:,i),1,0.5d0)
  !!enddo
  !!call gather_c16(sp,phi_bar_tot,size(phi_bar),0)
  !
  !if (rank==0) then
  !write(6,*) 'checking orthogonality'
  !allocate(overlap(nstates_r, nstates_r))
  !do i=nstates_r-5,nstates_r
  !  do j=nstates_r-5,nstates_r
  !     overlap(i,j) = sum(conjg(phi_bar_tot(:,:,:,i,1))*phi_bar_tot(:,:,:,j,1))*dv
  !     if (i .eq. j) then
  !       if (abs(overlap(i,i)-1d0) > 1e-8) then
  !         write(6,*) 'abs(norm-1d0) larger than 1e-8'
  !         stop 
  !       endif
  !     else
  !       if (abs(overlap(i,j)) > 1e-8) then
  !         write(6,*) '<psi_bar_i|psi_bar_j> larger than 1e-8'
  !         stop
  !       endif
  !     endif
  !  enddo
  !  !write(*,*) overlap(i,nstates_r-5:nstates_r)
  !enddo
  !endif

  !deallocate(overlap)

    !do i=1,nstates_r
    !   write(6,*) i, sum(conjg(phit_tot(:,:,:,i,1))*phit_tot(:,:,:,i,1))*dv
    !enddo

end subroutine read_wf_debug
