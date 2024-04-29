subroutine prop_vbar_split_sh(pin,pout, ik, dt)
  use main_mod
  use tddft_mod, only : tbar_taylor_power
  use mpi_lib_ours
  use atom_mod, only : p=> pawinfo, at=> atominfo, natom, atom_map
  implicit none
  integer :: ik
  real*8  :: dt
  complex*16 :: pin(nx,ny,nz), pout(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
  complex*16 :: tmp(nx,ny,nz)
  complex*16 :: ci = (0d0,1d0)
  real*8  :: nrm0, nrm1
  integer :: is, js, ms, i 

  !nrm0 = sqrt(dv*sum(abs(pin)**2))
  nrm0 = 1d0
  call add_cvtaylor_phi(pin,cin,ik)
  tmp = exp(-ci*dt*vks)*cin
  call add_cvtaylor_phi(tmp,pout,ik)
  !call add_cvtaylor_phi1(tmp,pout,ik)
  !nrm1 = sqrt(dv*sum(abs(tmp)**2))
  nrm1 = 1d0
  pout = pout*nrm0/nrm1
end subroutine
