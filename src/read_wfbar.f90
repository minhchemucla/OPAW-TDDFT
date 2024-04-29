subroutine read_wfbar
   use main_mod
   use param_mod
   use mpi_lib_ours

   implicit none
  
   character*9 ch
   integer :: i,j, ix,iy,iz
   integer :: nx_r, ny_r, nz_r
   integer :: nstates_r
   integer :: i_orb, nsp_r
   real*8  :: dx_r, dy_r, dz_r
   real*8, allocatable :: eig(:)
   real*8 :: tmp(nx*ny*nz)

   if(k_flg) stop 'not printing wf for kpoint'

   close(441)
   if(flg_bin)then; open (441,file='wf.bin',status='old',form='unformatted');
   else;           open (441,file='wf.txt',status='old');  
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
   call check(nb*nodes,nstates_r,' nb,nstates_r')

   read(441,*) ch !skip eigenvalues
   read(441,*) !skip eigenvalues
   read(441,*) ch !skip orbitals

   
   !debug
   !close(331)
   !open(unit=331, file='wf_debug.dat')
   
   do i=1,nb
     read(441,*) i_orb, nsp_r
     read(441,*) tmp
     phit_tot(:,:,:,i,1) = reshape(tmp,(/nx,ny,nz/))
       !debug
       !write(331,*) i_orb, nsp_r
       !write(331, *) real(phit_tot(:,:,:,i,1))
   enddo
   
end subroutine read_wfbar

