subroutine prop_tbar_split_sh(pin,pout, ik)
  use main_mod
  use tddft_mod, only : tbar_taylor_power, dt
  use mpi_lib_ours
  use atom_mod, only : p=> pawinfo, at=> atominfo, natom, atom_map
  implicit none
  integer :: ik
  complex*16 :: pin(nx,ny,nz), pout(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
  complex*16 :: tmp(nx,ny,nz)
  complex*16 :: ci = (0d0,1d0)
  real*8  :: nrm0, nrm1
  integer :: is, js, ms, i 
  real*8  :: wtime

  !nrm0 = sqrt(dv*sum(abs(pin)**2))
  nrm0 = 1d0
  call add_ctaylor_phi(pin,cin,ik)
  wtime = mpi_wtime()
  call fft3d_forward(nx,ny,nz,cin,cout)
  cout=cout*expk(:,:,:,ik)
  call fft3d_backward(nx,ny,nz,cout,tmp)
  wtime = mpi_wtime() - wtime
  write(6,*) 'rank, exp(T) wtime', rank, wtime
  !call add_ctaylor_phi1(tmp,pout,ik)
  call add_ctaylor_phi(tmp,pout,ik)
  !nrm1 = sqrt(dv*sum(abs(tmp)**2))
  nrm1 = 1d0
  pout = pout*nrm0/nrm1
end subroutine
