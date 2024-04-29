!-------------------------------------
!
! runge kutta. fixed step size.
!
!  input: n, x (x is a vector: x(n))  --- x and xout are complex here
!         t
!         dt  (small)
!         cderivs : the name of a subroutine (decalred external in main)
!               call cderivs(t, x, dxdt, n)
!
!
!
!   output: xout          
!-------------------------------------

subroutine real_rk4(x,n,t,dt,xout,cderivs)
  implicit none

  integer n, ierr
  real*8 dt, t 
  real*8 x(n), xout(n), dxdt(n),dxm(n),dxt(n),xt(n)

  external cderivs

  integer i

  real*8  dt6, dth, th, tf 

  if(n>10000101.or.n<1) then
     write(6,*)' in rk4 problem with n; recent value : ' ,n
     stop
  endif

  call cderivs(t, x, dxdt, n)

  dth = dt*0.5d0
  dt6  = dt/6.d0
  th = t + dth

  xt = x + dth*dxdt

  call cderivs(th, xt, dxt, n)

  xt = x + dth*dxt

  call cderivs(th, xt, dxm, n)

  xt  = x + dt*dxm
  dxm = dxt + dxm

  tf = t+dt
  call cderivs(tf, xt, dxt, n)
  
  xout = x + dt6*(dxdt + dxt + 2.d0* dxm)

end subroutine

!--------------
!  the following is a simple program to check.
!--------------
module simple_r_r_data
  implicit none
  save
  integer, parameter :: n=2
  real*8 aa(n)
  data aa / 0.5d0, 0.3d0 /
end module simple_r_r_data

!program drv_driver_real_rk4
!  implicit none
!  call driver_real_rk4
!end program drv_driver_real_rk4

subroutine driver_real_rk4
  use simple_r_r_data
  implicit none
  integer i,it
  real*8 :: dt, t, tt
  real*8 x(n),x0(n),xout(n)
  external r_derivs_simple

  x0 =5

  x = x0

  dt = 0.02d0
  do it=0,int(0.2/dt)
     t = it*dt
     call real_rk4(x,n,t,dt,xout,r_derivs_simple)
     x = xout

     tt = t+dt
     do i=1, n
        write(6,*) ' x',x(i),x0(i)*exp(-aa(i)*tt)
     enddo

  end do
end subroutine driver_real_rk4
 
subroutine r_derivs_simple(tt, x, dx, nin)
  use simple_r_r_data
  implicit none
  real*8 tt
  integer nin

  real*8 x(n),dx(n)

  call check(nin,n,' nin n ')

  dx = -aa*x
end subroutine r_derivs_simple
