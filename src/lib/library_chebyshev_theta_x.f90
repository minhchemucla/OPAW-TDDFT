subroutine cheb_coeff_theta_power(nchb, Temperature, mu, havg, DL,toll, co, nchbtop, power)
  ! chebyshev coeff. for representating f(x)= theta_smooth(x) = erfc(-x/Temperature)/2.d0 as
  ! f(x)=sum_k=0,1,2...nchb  co(k) T_n( x/DL), for -DL<x<DL.
  ! nchbtop is  the maxumum index subject to a tolerance on the values of toll.
  use mpi_lib_ours, only : rank
  implicit none

  integer nchb, power,ip, nchbtop ! nchbtop: output
  real*8  temperature, mu, havg, dl, toll, co(0:nchb) ! co is output

  integer i,m,N
  real*8  pi,w,D
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

  select case (power)
  case(1)
     f = f
  case(2)
     f = f*f
  case default
     write(6,*)' power should be 1 or 2, not ', power
     stop
  end select

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
end subroutine cheb_coeff_theta_power


subroutine cheb_coeff_theta2_x(nchb, Temperature, mu, havg, DL,toll, co, nchbtop)
  use mpi_lib_ours, only : rank
  implicit none

  integer nchb, nchbtop ! nchbtop: output
  real*8  temperature, mu, havg, dl, toll, co(0:nchb) ! co is output

  integer i,m,N
  real*8  pi,w,D
  complex*16, allocatable :: f(:)

  pi = dacos(-1.d0)

  m = 3 + nint( log(dble(nchb))/log(2.d0))
  N = 2**m  
  
  allocate(f(0:N-1), stat=i); call check(i,0,' falloc  ')
  do i=0,N-1
     w = 2d0*pi/dble(N) * i
     D = DL * cos(w)+havg
     f(i) = erfc((D-mu)/Temperature)/2.d0
     f(i) = f(i)*f(i)*D
  enddo

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
end subroutine cheb_coeff_theta2_x



subroutine cheb_coeff_theta_driver
  implicit none
  integer, parameter  :: nchb = 5000
  integer                nchbtop
  real*8 , parameter  :: Temperature =  0.02, DL = 2d0, mu=0.05, havg=0.2d0
  real*8              :: co(0:nchb)

  integer             :: k,i
  real*8              :: D, tk,tk1,tk2,S,h,ddl
  call cheb_coeff_theta_power(nchb, Temperature, mu, havg, DL,1.d-8, co, nchbtop,1)
  ddl = temperature/8d0
  do i=nint(-DL/ddl)+1, nint(DL/ddl)-1
     D = i*ddl 
     S =co(0)
     S =S+ co(1)*D/DL
     tk2= 1
     tk1= D/DL
     do k=2,nchbtop
        tk = 2*D/DL*tk1 - tk2
        S = S+ co(k)*tk
        tk2 = tk1
        tk1 = tk
     enddo
     h = D+havg
     write(778,*)h,D,S, erfc(-(mu-h)/Temperature)/2.d0,' h,d,s,1/2* erfc(-(mu-h)/T)  '
  enddo
  call flush(778)

end subroutine cheb_coeff_theta_driver

!program dummy
! call cheb_coeff_theta_driver
!end program dummy
