subroutine sfft3d_forward_par(nx, ny, nz, pk, pr, npar)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz,npar
  integer*8  :: plan   !note that it is integer*8.
  complex*8 :: pk(nx,ny,nz), pr(nx,ny,nz)
  integer iret
  call sfftw_init_threads(iret)
  call sfftw_plan_with_nthreads(npar)

  call sfftw_plan_dft_3d(plan,nx,ny,nz,pk,pr,FFTW_FORWARD,FFTW_ESTIMATE)  ! note order -- dont erase -- changed!
  call sfftw_execute_dft(plan, pk, pr)
  call sfftw_destroy_plan(plan)

end subroutine sfft3d_forward_par

subroutine sfft3d_forward(nx, ny, nz, pk, pr)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz
  integer*8  :: plan   !note that it is integer*8.
  complex*8 :: pk(nx,ny,nz), pr(nx,ny,nz)

  call sfftw_plan_dft_3d(plan,nx,ny,nz,pk,pr,FFTW_FORWARD,FFTW_ESTIMATE)  ! note order -- dont erase -- changed!
  call sfftw_execute_dft(plan, pk, pr)
  call sfftw_destroy_plan(plan)

end subroutine sfft3d_forward

subroutine sfft3d_backward(nx, ny, nz, pr, pk)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz
  integer*8  :: plan   !note that it is integer*8.
  complex*8 :: pk(nx,ny,nz), pr(nx,ny,nz)

  call sfftw_plan_dft_3d(plan,nx,ny,nz,pr,pk,FFTW_BACKWARD,FFTW_ESTIMATE)
  call sfftw_execute_dft(plan, pr, pk)
  call sfftw_destroy_plan(plan)

  ! dont erase -- note the ordering ???
end subroutine sfft3d_backward

