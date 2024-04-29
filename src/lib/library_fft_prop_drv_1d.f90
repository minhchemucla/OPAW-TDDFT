
subroutine fft_1drv_multiple(cain, ca, n, dx, lot)  ! ca can be same as cain
  implicit none

  real*8                  dx
  integer                 lot
  integer                n, j
  complex*16   ::  cain(n, lot), ca(n, lot)

  real*8                 pi
  real*8                 dp
  real*8                 pmom

  complex*16,parameter ::           eye = (0.d0,1.d0)

  pi = dacos(-1.d0)
  dp = 2.d0* pi / dx/n

  ca = cain
  call fft_multiple(ca, n, lot)

  do j=1, n
     if(j <= n/2) then
        pmom = dp*(j-1)
     else if(j == n/2+1 ) then
	pmom = 0
     else
        pmom = dp*(j-1-n)
     endif
     ca(j,:) = eye*pmom * ca(j,:)
  end do

  ca = conjg(ca)
  call fft_multiple(ca, n, lot)
  ca = conjg(ca)
  
  ca = ca / n
end subroutine fft_1drv_multiple


subroutine fft_2drv_multiple(cain, ca, n, dx, lot)  ! ca can be same as cain
  implicit none

  real*8                  dx
  integer                 lot
  integer                n, j
  complex*16   ::  cain(n, lot), ca(n, lot)

  real*8                 pi
  real*8                 dp
  real*8                 pmom

  pi = dacos(-1.d0)
  dp = 2.d0* pi / dx/n

  ca = cain
  call fft_multiple(ca, n, lot)

  do j=1, n
     if(j <= n/2) then
        pmom = dp*(j-1)
     else
        pmom = dp*(j-1-n)
     endif
     ca(j,:) = -pmom**2.d0 * ca(j,:)
  end do

  ca = conjg(ca)
  call fft_multiple(ca, n, lot)
  ca = conjg(ca)
  
  ca = ca / n
end subroutine fft_2drv_multiple


subroutine fft_prop_1d_good(cwf, n, mass, spacing, dt) ! out could=in.
  implicit none
  
  real*8     mass, dt, spacing, dp, pi, pmom
  integer    n, j
  complex*16 cwf(n)
  complex*16, parameter :: ci=(0.d0, 1.d0)


  call fft_good(cwf, n)

  pi = dacos(-1.d0)
  dp = 2.d0* pi / spacing/n

  do j=1, n
     if(j <= n/2) then
        pmom = dp*(j-1)
     else
        pmom = dp*(j-1-n)
     endif
     cwf(j) = exp(-ci*pmom**2.d0/2.d0/mass *dt)* cwf(j)
  end do

  cwf = conjg(cwf)
  call fft_good(cwf, n)
  cwf = conjg(cwf)
  
  cwf = cwf / n
end subroutine fft_prop_1d_good

subroutine fft_prop_imag_1d_good(cwf, n, mass, spacing, dt) ! out could=in.
  implicit none
  
  real*8     mass, dt, spacing, dp, pi, pmom
  integer    n, j
  complex*16 cwf(n)
  complex*16, parameter :: ci=(0.d0, 1.d0)


  call fft_good(cwf, n)

  pi = dacos(-1.d0)
  dp = 2.d0* pi / spacing/n

  do j=1, n
     if(j <= n/2) then
        pmom = dp*(j-1)
     else
        pmom = dp*(j-1-n)
     endif
     cwf(j) = exp(-pmom**2.d0/2.d0/mass *dt)* cwf(j)
  end do

  cwf = conjg(cwf)
  call fft_good(cwf, n)
  cwf = conjg(cwf)
  
  cwf = cwf / n
end subroutine 


subroutine dummy_driver_fft_2drv
  implicit none
  
  integer, parameter :: n = 128
  complex*16         :: ci = (0.d0,1.d0)

  complex*16         :: cwf(n,3),cwf2drv(n,3),cwf1drv(n,3), cwf2drv_anly(n,3)

  real*8             :: grid(n), xmin, xmax, dx, p0(3), x0, width

  integer            :: ix, j

  xmin = -4
  xmax = 6
  p0(1)   = 0.5; p0(2) = 0.53; p0(3) = -0.3
  x0   = 1.3
  width = 0.6

  dx = (xmax-xmin)/(n-1)
  do ix=1, n
     grid(ix) = xmin+(ix-1)*dx
  enddo
  
  do j=1,3
  cwf(:,j)     = exp(-(grid-x0)**2/2.d0/width**2)* 0.3* exp(ci*p0(j)*grid)
  cwf1drv(:,j) = cwf(:,j)*(-(grid-x0)/width**2 + ci*p0(j))
  cwf2drv(:,j) = cwf1drv(:,j)*(-(grid-x0)/width**2 + ci*p0(j)) + &
                     cwf(:,j)*(-1.d0/width**2)
  enddo

  call fft_2drv_multiple(cwf, cwf2drv_anly, n, dx, 3) !ca can be same as cain
  do j=1,3
     do ix=1, n
        write(6,888)ix,cwf2drv(ix,j),cwf2drv_anly(ix,j)
888     format(' ',i8,2f14.6,3x,2f14.6)
     enddo
  enddo

  write(6,*)' abs_diff ',sum(abs(cwf2drv-cwf2drv_anly))

end




