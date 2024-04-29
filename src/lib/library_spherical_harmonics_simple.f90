subroutine spherical_harmonics_simple(theta,l,m,yh)
  implicit none
  real*8 theta, yh, ct, st, pi
  integer l,m

  if(abs(m)>l) then
     write(6,*)' abs(m) bigger than l ',m,l
     stop
  endif

  pi = acos(-1.d0)
  
  ct = cos(theta)
  st = sin(theta)

  select case(l)
  case(0)
     yh = 1.d0/(2.d0*sqrt(pi))
  case(1)
     select case(m)
     case(-1)
        yh =  1.d0/2.d0 * sqrt(3.d0/2.d0/pi) * st
     case(0)
        yh =  1.d0*sqrt(3.d0/pi)* ct
     case(1)
        yh = -1.d0/2.d0 * sqrt(3.d0/2.d0/pi) * st
     end select
  case(2)
     select case(m)
     case(-2)
        yh =   1.d0/4.d0 * sqrt(15.d0/2.d0/pi) * st**2
     case(-1)
        yh =   1.d0/2.d0 * sqrt(15.d0/2.d0/pi) * ct* st
     case(0)
        yh =   1.d0/4.d0 * sqrt( 5.d0     /pi) * (3*ct**2 - 1) 
     case(1)
        yh = - 1.d0/2.d0 * sqrt(15.d0/2.d0/pi) * ct* st
     case(2)
        yh =   1.d0/4.d0 * sqrt(15.d0/2.d0/pi) * st**2
     end select
  case default
     write(6,*)' problem, l= ',l
     stop
  end select
end subroutine spherical_harmonics_simple
