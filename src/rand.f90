subroutine rand_r(g,n,dv)
  implicit none
  integer               :: n
  integer               :: i,j
  real*8                :: g(n)
  real*8                :: dv,r, b
  real*8, external      :: ran_ps

  b = dsqrt(3d0/dv)
  do i=1,n;  
     g(i) = b* ( 2d0*ran_ps() - 1d0 ) 
  enddo
end subroutine rand_r

subroutine rand_c(g,n,dv)
  implicit none
  integer                  :: n
  integer                  :: i
  complex*16               :: g(n)
  complex*16, parameter    :: ci=(0d0,1d0)
  real*8                   :: dv, b, pi
  real*8, external         :: ran_ps

  pi = dacos(-1d0)
  b = 1/dsqrt(dv)

  do i=1,n;  
     g(i) = b * exp(2d0*pi*ci*ran_ps())
  enddo
end subroutine rand_c

complex*16 function crand()
  implicit none
  real*8, external         :: ran_ps
  complex*16, parameter    :: ci=(0d0,1d0)
  !  real*8                   :: pi
  ! pi = dacos(-1d0)
  !crand = exp(2d0*pi*ci*ran_ps())
  ! modified  -- note

  if(ran_ps()<0.5d0) then
     crand = 1d0
  else
     crand = -1d0
  endif
end function crand

complex*16 function crand_unitcirc()
  implicit none
  real*8, external         :: ran_ps
  complex*16, parameter    :: ci=(0d0,1d0)
  real*8                   :: pi
  pi = dacos(-1d0)
  crand_unitcirc = exp(2d0*pi*ci*ran_ps())
end function crand_unitcirc

real*8 function rrand()
  implicit none
  real*8, external         :: ran_ps

  if(ran_ps()<0.5d0) then
     rrand = 1d0
  else
     rrand = -1d0
  endif
end function rrand
