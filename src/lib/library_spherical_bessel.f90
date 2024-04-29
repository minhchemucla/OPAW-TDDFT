! spehrical bessel function program
! calculates jn(z) = sqrt(pi/2/z) * J_(n+0.5)(z) for n=0,ntop, for x,ntop<1000000
! Using rercusion formula (abramowitz and stegun, 10.1.10)    f_(n-1)(z)+f_(n+1)(z) = (2n+1) f_n(z)/z
subroutine spherical_bessel( f, ntop, x)
  implicit none
  
  integer ntop
  integer i, nhigh
  real*8  f(0:ntop), x
  real*8  bm, b, bp
  real*8  s0, s1
  real*8, parameter :: xm = 1d6
  
  if(ntop>xm.or.abs(x)>xm.or.ntop<0) then
     write(6,*)' stopping, problem in spherical bessel:x,ntop,xm ',x,ntop,xm
     stop
  end if
  
  if(ntop==0) then
     f(0) = sb0(x)
     return
  endif
   
  if(abs(x)<1d-8) then
     f = 0d0
     f(0) = sb0(x)
     if(ntop.ge.1) f(1) = sb1(x)
     if(ntop.ge.2) f(2) = x**2d0/15d0-x**4d0/210d0
     if(ntop.ge.3) f(3) = x**3d0/105d0-x**5d0/1890d0
     return
  end if
  
  nhigh = 1.1* max(abs(x),dble(ntop))+500

  bp = 0d0
  b  = 1d0
  do i=nhigh,ntop+2,-1
     bm = dble(2*i+1)/x*b - bp
     bp = b
     b  = bm
     if(abs(b)>1d10) then
        bp = bp/b
        b  = 1d0
     end if
  end do

  ! now b= (scaled)  func(ntop+1);    bp= (scaled) func(ntop+2)
  
  do i=ntop+1,1,-1
     bm = dble(2*i+1)/x*b - bp
     f(i-1) = bm

     bp = b
     b  = bm

     if(abs(b)>1d10) then
        f(i-1:ntop) = f(i-1:ntop)/b
        bp = bp/b
        b  = 1d0
     end if
  enddo

  s0 = sb0(x)
  s1 = sb1(x)

  if(abs(f(0))>abs(f(1))) then
     f(0:ntop) = f(0:ntop)* s0/f(0)
  else
     f(0:ntop) = f(0:ntop)* s1/f(1)
  end if
  
contains
  real*8 function sb0(x) ! = sin(x)/x, abramowtiz-stegun 10.1.11
    implicit none
    real*8 x
    if(abs(x)>1d-8) then
       sb0 = sin(x)/x
    else
       sb0 = 1d0 - x**2/6d0 + x**4d0/120d0
    end if
  end function sb0

  real*8 function sb1(x) ! = sin(x)/x**2-cos(x)/x, abramowtiz-stegun 10.1.11
    implicit none
    real*8 x
    if(abs(x)>1d-4) then
       sb1 = (sin(x)-x*cos(x))/x**2
    else ! from wolphram alpha...
       sb1 = x/3d0- x**3d0/30d0 + x**5d0/840d0
    end if
  end function sb1
end subroutine spherical_bessel
  
