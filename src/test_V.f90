subroutine test_V
  use main_mod
  use mpi_lib_ours
  implicit none
  complex*16 :: tmp(nx,ny,nz)
  complex*16 :: p1(nx,ny,nz), p2(nx,ny,nz)
  real*8 :: vtest(nx,ny,nz)
  complex*16 :: vp1(nx,ny,nz), vp2(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz), cout(nx,ny,nz)
  complex*16 :: a, b

  call rand_c(p1,nn,dv)
  call rand_c(p2,nn,dv)
  call rand_c(vtest,nn,dv)

  
  vp1 = vtest*p1
  vp2 = vtest*p2

  a = sum(conjg(p1)*vp2)*dv
  b = sum(conjg(vp1)*p2)*dv
  if(rank==0)write(6,*) '<1|V2>, <1V|2>, <1|V2>-<1V|2>', a, b, a-b

  call sync_mpi

  stop
  
end subroutine 
