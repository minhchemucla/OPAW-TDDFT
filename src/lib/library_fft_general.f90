!
! fft_general_multiple is actually less general but does same as 
! fft_easy_multiple
!
subroutine fft_general_multiple(ca, caout, nx, nk, lot, dx, xmin,kmin, dk)
  implicit none

  complex*16, parameter :: ci = (0.d0, 1.d0)

  integer nx, nk, lot, ilot, st, is_it_power2, ix, ik
  real*8  dx, xmin, kmin, dk, x, kk
  complex*16 ca(nx, lot), caout(nk, lot)
  complex*16, allocatable :: cx(:), ck(:)
  
  if((nx/=nk).or.is_it_power2(nx)/=1) then
     write(6,*)' nx, nk, not ready for diff. ones or non 2** power', &
          nx, nk
     stop
  endif

  allocate(cx(nx), ck(nk), stat=st); call check(st,0,' stkx0   ')

  do ilot=1, lot
     cx = ca(:, ilot)

     !
     ! what we need is 
     !      ck(ik) = sum_x exp(-i*k(ik)*x(ix)) cx(ix)
     !
     !  what we have is
     !      ck_have(ik) = sum_x exp(-i*(k(ik)-kmin)*x(ix)-xmin)cx(ix)

     !
     !  So first define
     !      cx_tilde(ix) = exp(-i*(kmin)*(x(ix)-xmin))*cx(ix)
     !
     !
     !  Next fft ck_tilde, to get
     !      ck_tilde = (post fft) 
     !               sum_x exp(-i*(k(ik))*(x(ix)-xmin)*cx(ix)
     !    and finally
     !      ck(ik) = ck_tilde * exp(-i*k(ik)*xmin)*ck_tilde
     
     do ix=1, nx
        x = xmin+(ix-1)*dx
        cx(ix) = exp(-ci*kmin*(x-xmin))*cx(ix)
     enddo

     call fft_good(cx, nx)
     ck = cx

     do ik=1,nk
        kk = kmin + (ik-1)*dk
        ck(ik) = exp(-ci*kk*xmin)* ck(ik)
     enddo

     caout(:, ilot) = ck
  enddo
  deallocate(cx, ck)
end subroutine fft_general_multiple
  
subroutine drive_fft_general
  implicit none

  integer, parameter :: nx=8, nk=nx
  integer ix, ik
  
  complex*16 cx(nx), ck(nk)
  real*8     dx, dk, pi, xmin, kmin
  
  dx = 0.2d0
  pi = dacos(-1.d0)
  xmin = -dx*(nx/2.d0-0.5d0)
  dk = 2*pi/(nx*dx)
  kmin = -dk*(nk/2.d0-0.5)

  do ix=1,nx
     cx(ix) = cos(0.5578*ix)+cos(2.598*ix**2)
  enddo

  call fft_general_multiple(cx, ck, nx, nk, 1, dx, xmin,kmin, dk)

  do ik=1,nk
     write(6,*)ik, ck(ik)
  enddo
end subroutine drive_fft_general

!program drive_fft_general_prog
!call drive_fft_general
!end
