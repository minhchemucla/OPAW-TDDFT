!-------------------------------------
!
! runge kutta. fixed step size.
!
!  input: n, x (x is a vector: x(n))
!         t
!         dt  (small)
!         derivs : the name of a subroutine (decalred external in main)
!               call derivs(t, x, dxdt, n)
!
!
!
!   output: xout          
!-------------------------------------

subroutine real_rk4_new(x,n,t,dt,xout,derivs)
  implicit none

  integer n, ierr

  real*8 dt, t, x(n), xout(n)   ! could be the same -- check

  external derivs

  integer i

  real*8  dt6, dth, th 

  real*8, allocatable, dimension(:) :: dxdt, dxm,dxt, xt

  if(n>10000101.or.n<1) then
     write(6,*)' in rk4 problem with n; recent value : ' ,n
     stop
  endif

  allocate( dxdt(n), stat=ierr); call check(ierr,0,' dxdt ')
  allocate( dxm(n),  stat=ierr); call check(ierr,0,' dxm ')
  allocate( dxt(n),  stat=ierr); call check(ierr,0,' dxt ')
  allocate( xt(n),   stat=ierr); call check(ierr,0,' xt ')

  call derivs(t, x, dxdt, n)

  dth = dt*0.5d0
  dt6  = dt/6.d0
  th = t + dth

  xt = x + dth*dxdt
  call derivs(th, xt, dxt, n)

  xt = x + dth*dxt

  call derivs(th, xt, dxm, n)

  xt  = x + dt*dxm
  dxm = dxt + dxm

  call derivs(t+dt, xt, dxt, n)
  
  xout = x + dt6*(dxdt + dxt + 2.d0* dxm)

  deallocate(dxdt, dxm, dxt, xt)

end subroutine real_rk4_new


!--------------
!  the following is a simple program to check.
!--------------


module rk_trial_simple
  implicit none
  save
  real*8,  parameter :: omega2 = 1.0
end module rk_trial_simple

subroutine derivs_rk_simple(t, x, dxdt, n)
  use rk_trial_simple
  implicit none

  integer n
  real*8 t, x(n), dxdt(n)
  
  dxdt(1) = x(2)
  dxdt(2) = -omega2 * x(1)

end subroutine derivs_rk_simple

!program rk_trial_simple_program
subroutine rk_trial_simple_program
  use rk_trial_simple
  implicit none
  external derivs_rk_simple
  integer it
  integer, parameter :: n=2
  real*8 t, x(n), dt, tmax

  tmax = 15
  t = 0
  dt = 0.001
  x(1) = 0.5
  x(2) = 0


  do it=1,nint(tmax/dt)+1
     t = (it-1)*dt
     write(10,*)'t,x,v ',t,x
     write(11,*) t, x(1)
     write(12,*) t, x(2)
     call  real_rk4_new(x,n,t,dt,x,derivs_rk_simple)     
  enddo
     
end








