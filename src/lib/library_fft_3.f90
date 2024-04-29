

! fft easy: use as call fft_easy(zt, cw, N, tmin, dt)
!
! input: zt, N, tmin, dt.  
!    N should be 2**power(i.e., 2,4,8, etc.)
!    zt(N) can be real or complex
!
! Output: cw (complex).  
!
! this routine returns:
!  cw(j) = 1/N * sum_k exp(-ci*wj*tk) zt(tk)
!  where tk = tmin+(k-1)*dt
!  and   wj = wmin+(j-1)*dw 
!
! It calls fft_multiple, which does the same job, except that it (fft_multiple)
! assumes that tmin and wmin are zero.

module lib_interface
  interface fft_easy
     module procedure  fft_easy_cinput
     module procedure  fft_easy_rinput 
  end interface

  interface fft_inv_easy
     module procedure  fft_inv_easy_cinput
     module procedure  fft_inv_easy_rinput
  end interface

contains

subroutine fft_easy_cinput(ct, cw, N, tmin, dt)
  implicit none

  integer  :: N
  real*8              :: tmin, dt
  complex*16          :: cw(N), ct(N)
  integer             :: is_it_power2, st
  external             is_it_power2  ! This routine checks whether N=2**power

  integer                          :: j
  complex*16 , parameter           :: ci=(0.d0,1.d0) 
  real*8                           :: wmin, dw, pi
  
  !
  ! first we need to verify that nout is a product of 2.
  !  
  
  if(is_it_power2(N) /=1) then
     write(6,*)' problem with N = ',N,' not power of 2 '
     stop
  end if
  
  !
  ! Now define some variables.
  !
  pi   = dacos(-1.d0)
  dw   = 2*pi/dble(N)/dt
  wmin = -dw*N/2
  
  cw = ct
  !
  ! cw(j) = 1/N * sum_k exp(-ci*wj*tk) zt(tk)
  !       = 1/N * sum_k exp(-ci*(wmin+(j-1)*dw)*(tmin+(k-1)*dt))
  !       = 1/N * sum_k exp(-ci*wmin*tmin) * exp(-ci*wmin*dt*(k-1)) * &
  !                     exp(-ci*tmin*dw*(j-1))* exp(-ci*(k-1)*(j-1)*dt*dw)
  
  wmin = -(N/2)*dw
  
  do j=1, N
     cw(j) = cw(j)* exp(-ci*wmin*dt*(j-1))
  enddo
  
  call fft_multiple(cw, N, 1) ! fft calculation; 
                              ! cw(j)=sum_k cw_old(k) exp(-i(k-1)(j-1)2pi/N)
  do j=1,N
     cw(j) = cw(j) * exp(-ci*tmin*dw*(j-1))
  enddo
  
  cw = cw* exp(-ci*tmin*wmin)/N
end subroutine fft_easy_cinput

subroutine fft_easy_rinput(rt, cw, N, tmin, dt)
  implicit none
  
  integer    N
  real*8     rt(N)
  complex*16 cw(N)
  real*8     tmin,  dt
  complex*16, allocatable :: ct(:)

  integer st
  
  allocate(ct(N), stat=st)
  if(st/=0) then
     write(6,*)' problem in ct_fft ',st
     stop
  endif
  
  ct = rt
  call fft_easy_cinput(ct, cw, N, tmin, dt)
  deallocate(ct)
end subroutine fft_easy_rinput

subroutine fft_inv_easy_rinput(rt, cw, N, wmin, dw, tmin, dt)
  implicit none
  
  real*8 dt, tmin
  integer    N
  real*8     rt(N)
  complex*16 cw(N)
  real*8     wmin,  dw
  complex*16, allocatable :: ct(:)
  integer st
  
  allocate(ct(N), stat=st)
  if(st/=0) then
     write(6,*)' problem in ct_fft ',st
     stop
  endif
  
  ct = rt
  call fft_easy(ct, cw, N, wmin, dw)
  cw = conjg(cw)*N 
  deallocate(ct)
end subroutine fft_inv_easy_rinput

subroutine fft_inv_easy_cinput(ct, cw, N, wmin, dw)
  implicit none
  
  integer    N
  complex*16 ct(N)
  complex*16 cw(N)
  real*8     wmin,  dw
  
  ct = conjg(ct)
  call fft_easy(ct, cw, N, wmin, dw)
  cw = conjg(cw)*N 
  ct = conjg(ct)

end subroutine fft_inv_easy_cinput
end module lib_interface








