subroutine st_phi(pin,hp,ik)
  use main_mod
  use atom_mod
  use atom_mod, only : p=> pawinfo, at=> atominfo
  implicit none

  integer :: ik
  complex*16 :: pin(nx,ny,nz)
  complex*16 :: hp(nx,ny,nz),tmp(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
  
  call fft3d_forward(nx,ny,nz,pin,cout)
  cout=cout*ek(:,:,:,1)
  call fft3d_backward(nx,ny,nz,cout,tmp)
  call sn_phi(tmp,hp,1,-1d0)
  !call cmat_phi(pin,cout,ik)
  !call c_phi(pin,cout,ik)
  !hp = tmp + cout
end subroutine

subroutine t_phi(pin,tp,ik)
  use main_mod
  use atom_mod
  use atom_mod, only : p=> pawinfo, at=> atominfo
  implicit none

  integer :: ik
  complex*16 :: pin(nx,ny,nz)
  complex*16 :: tp(nx,ny,nz),tmp(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
  
  call fft3d_forward(nx,ny,nz,pin,cout)
  cout=cout*ek(:,:,:,ik)
  call fft3d_backward(nx,ny,nz,cout,tp)
end subroutine

subroutine sum_inv_T_phi(pin,stp,ik)
  use main_mod
  use atom_mod
  use atom_mod, only : p=> pawinfo, at=> atominfo
  implicit none
  integer :: ixs,iys,izs,ixe,iye,ize,jx,jy,jz
  integer :: ib,jb,ia,ix,iy,iz,is,js,it,ix1,iy1,iz1,ig,ik
  complex*16 :: pin(nx,ny,nz)
  complex*16 :: stp(nx,ny,nz),tmp(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
  complex*16, allocatable :: ca(:),sn(:)
  
  call fft3d_forward(nx,ny,nz,pin,cout)
  cout=cout*ek(:,:,:,1)
  call fft3d_backward(nx,ny,nz,cout,tmp)
  call sn_phi(tmp,stp,1,-1d0)
  stp = stp - tmp
end subroutine
