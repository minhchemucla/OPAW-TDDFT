subroutine fft2d_forward(nx, ny, pk, pr)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny
  integer*8  :: plan   !note that it is integer*8.
  complex*16 :: pk(nx,ny), pr(nx,ny)

  call dfftw_plan_dft_2d(plan,nx,ny,pk,pr,FFTW_FORWARD,FFTW_ESTIMATE)  ! note order -- dont erase -- changed!
  call dfftw_execute_dft(plan, pk, pr)
  call dfftw_destroy_plan(plan)

end subroutine fft2d_forward

subroutine fft2d_backward(nx, ny, pr, pk)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny
  integer*8  :: plan   !note that it is integer*8.
  complex*16 :: pk(nx,ny), pr(nx,ny)

  call dfftw_plan_dft_2d(plan,nx,ny,pr,pk,FFTW_BACKWARD,FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, pr, pk)
  call dfftw_destroy_plan(plan)

! dont erase -- note the ordering ???
end subroutine fft2d_backward



subroutine fft3d_backward_realcast(nx, ny, nz, pk, pr)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz
  integer*8  :: plan   !note that it is integer*8.
  real*8 :: pk(2,nx,ny,nz), pr(2,nx,ny,nz)
  call dfftw_plan_dft_3d(plan,nx,ny,nz,pk,pr,FFTW_BACKWARD,FFTW_ESTIMATE)  ! note order -- dont erase -- changed!
  call dfftw_execute_dft(plan, pk, pr)
  call dfftw_destroy_plan(plan)
end subroutine fft3d_backward_realcast

subroutine fft3d_forward_realcast(nx, ny, nz, pk, pr)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz
  integer*8  :: plan   !note that it is integer*8.
  real*8 :: pk(2,nx,ny,nz), pr(2,nx,ny,nz)
  call dfftw_plan_dft_3d(plan,nx,ny,nz,pk,pr,FFTW_FORWARD,FFTW_ESTIMATE)  ! note order -- dont erase -- changed!
  call dfftw_execute_dft(plan, pk, pr)
  call dfftw_destroy_plan(plan)
end subroutine fft3d_forward_realcast

subroutine fft3d_forward(nx, ny, nz, pk, pr)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz
  integer*8  :: plan   !note that it is integer*8.
  complex*16 :: pk(nx,ny,nz), pr(nx,ny,nz)

  call dfftw_plan_dft_3d(plan,nx,ny,nz,pk,pr,FFTW_FORWARD,FFTW_ESTIMATE)  ! note order -- dont erase -- changed!
  call dfftw_execute_dft(plan, pk, pr)
  call dfftw_destroy_plan(plan)

end subroutine fft3d_forward

subroutine fft3d_backward(nx, ny, nz, pr, pk)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz
  integer*8  :: plan   !note that it is integer*8.
  complex*16 :: pk(nx,ny,nz), pr(nx,ny,nz)

  call dfftw_plan_dft_3d(plan,nx,ny,nz,pr,pk,FFTW_BACKWARD,FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, pr, pk)
  call dfftw_destroy_plan(plan)

! dont erase -- note the ordering ???
end subroutine fft3d_backward

subroutine fft3d_r2cinplace(nx, ny, nz, ar)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz
  integer*8  :: plan   !note that it is integer*8.
  real*8     :: ar(2*(nx/2+1),ny,nz)
  
  call dfftw_plan_dft_r2c_3d(plan,nx,ny,nz,ar,ar,FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan, ar, ar)
  call dfftw_destroy_plan(plan)
end subroutine fft3d_r2cinplace

subroutine fft3d_c2rinplace(nx, ny, nz, ar)
  use mpi_lib_ours, only : rank
  implicit none
  include "fftw3.f_include"
  integer nx,ny,nz
  integer*8  :: plan   !note that it is integer*8.
  real*8     :: ar(2*(nx/2+1),ny,nz)

  call dfftw_plan_dft_c2r_3d(plan,nx,ny,nz,ar,ar,FFTW_ESTIMATE)
  call dfftw_execute_dft_c2r(plan, ar, ar)
  call dfftw_destroy_plan(plan)
end subroutine fft3d_c2rinplace

