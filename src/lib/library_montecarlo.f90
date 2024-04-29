subroutine get_new_pos_gen_pot(xv, xa, wa, LL)
  
  ! this subroutine chooses a new x such that w is distributed randomly between
  ! 0 and 1.  Inputs are x, xa and wa. It is assumed that wa(1)=0 & wa(LL)=1.d0
  ! (These features are checked).  It is also assumed that wa is monotonic --
  ! a check is optionally done in the called routine.

  implicit none
  integer  LL
  real*8   wa(LL)
  real*8   xa(LL)
  real*8   ran_ps
  external ran_ps
  real*8   wv, xv

! careful on the interpolation -- wa is the "x","xa" is the"dependent variable"
 
  if((abs(wa(1))>1.d-8).or.abs(wa(LL)-1)>1.d-8) then
     write(6,*)' wa(1), wa(LL) '
     write(6,*)  wa(1)
     write(6,*)  wa(LL)
     stop
  endif

  wv = ran_ps()
  call interpol8_gen(wa, xa, wv, xv, LL)   ! Note order.

end subroutine get_new_pos_gen_pot


!
! this subroutine integrates, using trapezoidal rule, ya. so that
! w(x) = integral y(x') dx' from xmin.
! LLater make a more sophisticated function.
! For now assume that ya is vanishing at the bottom end, and x
! is equispaced.

subroutine wintg_trpz(dx, ya, wa, LL)
  implicit none
  
  integer LL, i
  real*8  ya(LL), wa(LL), dx
  
  wa(1) = 0.d0
  wa(2) = (ya(1)+ya(2))*dx/2.d0
  
  if (wa(2)>1.d-5) then
     write(6,*)' wa(2) is quite large.  Are you sure theres no mistake ? '
     stop
  endif
  
  do i=1+1,LL-1
     wa(i+1) = wa(i-1) + 2*dx*ya(i)
  enddo
  
end subroutine wintg_trpz

subroutine interpol8_gen(ra, sa, r, s, LL)   ! outputs s gien r and sa,ra)
  implicit none
  integer LL, jlo, jhi, j
  real*8  ra(LL)
  real*8  sa(LL)
  real*8  r
  real*8  s
  real*8  h, a, b
  integer  :: icheck, i

  icheck= 0
!
! this subroutine assumes that the order is monotonic increasing in ra;
! you can check it explicitly

  if(icheck==1) then
     do i=1, LL-1
        call check_real_le(ra(i), ra(i+1), 'rai      ' )
     enddo
  end if
  
  if(r>ra(LL).or.r<ra(1)) then
     write(6,*)' problem ; ra(1),r,ra(LL) '
     write(6,*)ra(1)
     write(6,*)r
     write(6,*)ra(LL)
     stop
  end if

  jlo=1
  jhi=LL

  do while (jhi-jlo>1) 
     j=(jhi+jlo)/2
     if(ra(j).gt.r)then
        jhi=j
     else
        jlo=j
     endif
  end do

  h=ra(jhi)-ra(jlo)
  
  if (h.eq.0.) then
     write(6,*) 'bad ra input.'
     stop
  endif
  a=(ra(jhi)-r)/h
  b=(r-ra(jlo))/h
  s=a*sa(jlo)+b*sa(jhi)

end subroutine interpol8_gen
  

subroutine interpol8_smpl_equi(xmin, dx, x, w, wa, LL)
  implicit none
  integer LL
  real*8 x, xj,xjh
  real*8 xmin
  real*8 dx
  real*8 w       ! output
  real*8 wa(LL)
  
  integer j, jh

  call check_real_le(0.d0,dx,' pl:dx ')

  j=int((x-xmin)/dx)+1
 
  if(j<1) then
     w=wa(1)
     return
  endif

  if(j>LL) then 
     w=wa(LL)
     return
  endif

  jh = min(j+1,LL)
  
  xj = xmin+(j-1)*dx
  xjh = xj+dx

  if(x>xjh.or.x<xj) then
     write(6,*)' problem '
     write(6,*)' xj, x, xjh '
     write(6,*) xj
     write(6,*) x
     write(6,*) xjh
     write(6,*)' j ',j
     stop
  endif
  w = (x-xj)/dx*wa(j) + (xjh-x)/dx*wa(jh)
end subroutine interpol8_smpl_equi

