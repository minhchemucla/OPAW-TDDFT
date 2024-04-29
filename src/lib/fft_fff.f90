subroutine fft_fff(nx, ny, nz, ck, cr, rank) ! emulates cuda call for cpu using fftw
  implicit none
  integer nx, ny, nz, rank
  complex*16 cr(nx,ny,nz)
  complex*16 ck(nx,ny,nz)
  call fft3d_forward(nx, ny, nz, ck, cr)
end subroutine fft_fff
