subroutine fft_fff1(nx, ny, nz, p, pout,rank)  
  implicit none
  integer nx,ny,nz,rank
  complex*16 p(nx,ny,nz),pout(nx,ny,nz)
  call fft3d_forward(nx, ny, nz, p, pout)
end subroutine fft_fff1

subroutine get_vcon( dens, nx, ny, nz, vkn,vh, rank)
  implicit none
  integer nx, ny, nz,rank,i
  logical, parameter :: use_r2c=.true.
  real*8  dens(nx,ny,nz)
  real*8   vkn(nx,ny,nz)
  real*8    vh(nx,ny,nz)
  
  complex*16, allocatable :: ca(:,:,:)

  if(use_r2c) then
     call rconv(dens,vkn,vh,nx,ny,nz)
     vh = vh*(dble(nx)*dble(ny)*dble(nz))
  else
     allocate(ca(nx,ny,nz),stat=i); if(i/=0) stop ' ca '
     ca = dens
     call fft3d_forward(nx, ny, nz, ca, ca)
     ca = ca * vkn
     call fft3d_backward(nx, ny, nz, ca, ca)
     vh = ca
     deallocate(ca)
  end if
     
end subroutine get_vcon

subroutine get_cvco( cdens, nx, ny, nz, vkn,cvh, rank)
  implicit none
  integer nx, ny, nz,rank,i
  complex*16 cdens(nx,ny,nz)
  real*8   vkn(nx,ny,nz)
  complex*16       cvh(nx,ny,nz)
  
  complex*16, allocatable :: ca(:,:,:)
  allocate(ca(nx,ny,nz),stat=i); if(i/=0) stop ' ca '
  ca = cdens
  call fft3d_forward(nx, ny, nz, ca, ca)
  ca = ca * vkn
  call fft3d_backward(nx, ny, nz, ca, ca)
  cvh = ca

  deallocate(ca)
end subroutine get_cvco

subroutine fft3d_forward_many(nx, ny, nz, ns, pt)
  implicit none
  integer nx,ny,nz,ns,is
  complex*16 pt(nx,ny,nz,ns)
  do is=1,ns
     call fft3d_forward(nx,ny,nz,pt(1,1,1,is),pt(1,1,1,is))
  enddo
end subroutine fft3d_forward_many

subroutine fft3d_backward_many(nx, ny, nz, ns, pt)
  implicit none
  integer nx,ny,nz,ns,is
  complex*16 pt(nx,ny,nz,ns)
  do is=1,ns
     call fft3d_backward(nx,ny,nz,pt(1,1,1,is),pt(1,1,1,is))
  enddo
end subroutine fft3d_backward_many
