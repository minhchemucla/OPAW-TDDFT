subroutine h_phi_debug_2(pin,hp)
  use main_mod
  use atom_mod
  use atom_mod, only : p=> pawinfo, at=> atominfo
  use tddft_mod 
  implicit none

  
  complex*16 :: pin(nx,ny,nz)
  complex*16 :: hp(nx,ny,nz),tmp(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)


end subroutine
