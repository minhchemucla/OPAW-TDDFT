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

subroutine complex_rk4(x,n,t,dt,xout,cderivs)
  implicit none

  integer n, ierr
  real*8 dt, t 
  complex*16 x(n), xout(n), dxdt(n),dxm(n),dxt(n),xt(n)

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





!-------------------------------------
!
! runge kutta. fixed step size, mixed complex and real.
!
!  input: n, x (x is a vector: x(n))  --- x and xout are complex here
!         m, y (y is real here)
!         t
!         dt  (small)
!         cr_derivs : the name of a subroutine (decalred external in main)
!               call cr_derivs(t, x, y, dxdt, dydt, n, m)
!
!
!
!   output: x,y
!-------------------------------------

subroutine complex_real_rk4(x,y,n,m,t,dt,cr_derivs)
  implicit none

  integer n, m,ierr
  real*8 dt, t 
  complex*16 x(n), dxdt(n),dxm(n),dxt(n),xt(n)
  real*8     y(m), dydt(m),dym(m),dyt(m),yt(m)

  external cr_derivs

  integer i

  real*8  dt6, dth, th, tf 

  if(n>10000101.or.n<1) then
     write(6,*)' in rk4 problem with n; recent value : ' ,n
     stop
  endif

  call cr_derivs(t, x, dxdt, y, dydt, n, m)


  dth = dt*0.5d0
  dt6  = dt/6.d0
  th = t + dth

  xt = x + dth*dxdt
  yt = y + dth*dydt

  call cr_derivs(th, xt, dxt, yt, dyt, n,m)

  xt = x + dth*dxt
  yt = y + dth*dyt

  call cr_derivs(th, xt, dxm, yt, dym, n)

  xt  = x + dt*dxm
  dxm = dxt + dxm
  yt  = y + dt*dym
  dym = dyt + dym

  tf = t+dt
  call cr_derivs(tf, xt, dxt, yt, dyt, n,m)
  
  x = x + dt6*(dxdt + dxt + 2.d0* dxm)
  y= y + dt6*(dydt + dyt + 2.d0* dym)

end subroutine

!--------------
!  the following is a simple program to check.
!--------------
module simple_c_r_data
  implicit none
  save
  integer, parameter :: n=2,m=3
  complex*16, parameter :: ci = (0.d0,1.d0)
  real*8 aa(n), bb(m)
  data aa / 0.5d0, 0.3d0 /
  data bb /0.1d0,0.2d0,0.25d0/
end module simple_c_r_data

subroutine driver_complex_real_rk4
  use simple_c_r_data
  implicit none
  integer i,it
  real*8 :: dt, t, y(m),y0(m), tt
  complex*16 x(n),x0(n)
  external cr_derivs_simple

  x0 =1+ci*2
  y0 = 5

  x = x0
  y = y0 

  dt = 0.02d0
  write(6,*)' y = ',y
  do it=0,int(0.2/dt)
     t = it*dt

     call complex_real_rk4(x,y,n,m,t,dt,cr_derivs_simple)
     
     tt = t+dt
     do i=1, n
        write(6,*) ' x',x(i),x0(i)*exp(ci*aa(i)*tt)
     enddo

     do i=1, m
        write(6,*) ' y',y(i),y0(i)*exp(-bb(i)*tt)
     enddo
  enddo
end subroutine driver_complex_real_rk4
 
subroutine cr_derivs_simple(tt, x, dx, y, dy, nin,min)
  use simple_c_r_data
  implicit none
  real*8 tt
  integer nin, min

  complex*16 x(n),dx(n)
  real*8     y(m), dy(m)

  call check(nin,n,' nin n ')
  call check(min,m,' min m ')

  dx = ci*aa*x
  dy = -bb*y
end subroutine cr_derivs_simple
