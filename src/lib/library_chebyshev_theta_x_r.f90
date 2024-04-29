subroutine cheb_coeff_theta_power_r(nchb, Temperature, mu, havg, DL,toll, co, nchbtop, power_r)
  ! chebyshev coeff. for representating f(x)= theta_smooth(x) = erfc(-x/Temperature)/2.d0 as
  ! f(x)=sum_k=0,1,2...nchb  co(k) T_n( x/DL), for -DL<x<DL.
  ! nchbtop is  the maxumum index subject to a tolerance on the values of toll.
  use mpi_lib_ours, only : rank
  implicit none

  integer nchb, ip, nchbtop ! nchbtop: output
  real*8  temperature, mu, havg, dl, toll, co(0:nchb) ! co is output

  integer i,m,N
  real*8  pi,w,D,power_r
  complex*16, allocatable :: f(:)

  pi = dacos(-1.d0)

! N  is about 4 or 8 times biger than nchb, and is 2**m

  m = 3 + nint( log(dble(nchb))/log(2.d0))
  N = 2**m  
  
  allocate(f(0:N-1), stat=i); call check(i,0,' falloc  ')

  do i=0,N-1
     w = 2d0*pi/dble(N) * i
     D = DL * cos(w)+havg
     f(i) = erfc((D-mu)/Temperature)/2.d0
  enddo

  f = f**power_r

  call fftsa(N,f,m)
  f = 2.d0*f/dble(N)
  f(0) = f(0)/2.d0

  co   = f(0:nchb)

  maxneed : do i=nchb,10,-1
     if(abs(co(i))>toll) exit maxneed 
  enddo maxneed

  nchbtop =i
  if(nchbtop>nchb-5.or.nchbtop<2) then
     write(6,*)' problem in nchbtop,nchb ',nchbtop,nchb; stop
  end if
  
  deallocate(f,stat=i); call check(i,0,' fdealloc  ')
end subroutine cheb_coeff_theta_power_r


