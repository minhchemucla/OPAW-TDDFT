subroutine test_T
  use main_mod
  use mpi_lib_ours
  implicit none
  complex*16 :: tmp(nx,ny,nz)
  complex*16 :: p1(nx,ny,nz), p2(nx,ny,nz)
  complex*16 :: tp1(nx,ny,nz), tp2(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz), cout(nx,ny,nz)
  complex*16 :: a, b

  call rand_c(p1,nn,dv)
  cin=p1
  call fft3d_forward(nx,ny,nz,cin,cout)
  cout=cout*ek(:,:,:,1)
  call fft3d_backward(nx,ny,nz,cout,tp1)

  call rand_c(p2,nn,dv)
  cin=p2
  call fft3d_forward(nx,ny,nz,cin,cout)
  cout=cout*ek(:,:,:,1)
  call fft3d_backward(nx,ny,nz,cout,tp2)

  a = sum(conjg(p1)*tp2)*dv
  b = sum(conjg(tp1)*p2)*dv
  if(rank==0)write(6,*) '<1|T2>, <1T|2>, <1|T2>-<1T|2>', a, b, a-b

  call sync_mpi

  stop
  
end subroutine 
