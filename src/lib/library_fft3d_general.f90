subroutine fft3d_general(p, nx, ny, nz,sign)
  use mpi_lib_ours, only : rank
  implicit none
  integer nx, ny, nz,sign
  complex*16 p(nx,ny,nz)

  if(abs(sign).ne.1) stop ' sign in fft3d_general '

  if(sign==1) p = conjg(p)
  call fft_fff(nx, ny, nz, p, p,rank)  ! for-fortran-forward
  if(sign==1) p = conjg(p)

end subroutine fft3d_general

