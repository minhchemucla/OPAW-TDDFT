!
!  dummy routine to check crank nicolson.
!
subroutine drvr_crank()
  implicit none
  integer, parameter :: n=20, lot=5

  integer l, j
  complex*16 ket(n,lot),ket_new(n,lot), diag(n), off_diag(n)

  real*8 dt, ran_ps_neg; external ran_ps_neg
  dt = 0.03

  do l=1, lot
     do j=1, n
        ket(j, l) = ran_ps_neg()
        diag(j) =   ran_ps_neg()
        off_diag(j) = ran_ps_neg()
     enddo
  enddo
  off_diag(:) = off_diag(1)
  diag(:)     = diag(1)

  call c_cranck_nicolson(ket, diag, off_diag, dt, n, lot, ket_new)
  write(6,*)sum(abs(ket_new-ket))/sum(abs(ket)), ' diff  -shouldnt be zero'

  !  write(6,*)' dt       ',dt
  !  write(6,*)' diag     ',diag
  !  write(6,*)' off_diag ',off_diag

  !  write(6,*)' final '
  !  write(6,*)' ket ',     ket
  !  write(6,*)' ket_new  ',ket_new

end subroutine drvr_crank

!
!  calculates (1-H*dt/2)/(1+H*dt/2)*ket --> ket_new.
!  here all is complex, except dt. H is tridiagonal (the off-diag
!  start from 1).  But note that you will need slight modifications
!  to deal with dt --> i*dt. (Just multiply H by i*dt rahter than i)
!
module cranck_nicolson_work
  implicit none
  save
  complex*16, allocatable, dimension(:)   :: diag_temp,off_diag_temp
  complex*16, allocatable, dimension(:,:) :: ket_mid
end module cranck_nicolson_work

subroutine c_cranck_nicolson(ket, diag, off_diag, dt, n, lot, ket_new)
  use cranck_nicolson_work
  implicit none
  integer n, lot, j
  complex*16 off_diag( n-1), diag(n)
  complex*16  ket(n, lot), ket_new(n, lot)
  real*8  dt

  
  if(allocated(diag_temp)) then
     if(size(diag_temp).ne.n)then
        deallocate(diag_temp, off_diag_temp, stat=j); call check(j,0,'  deall diag_temp ')
     endif
  endif


  if(allocated(ket_mid)) then
     if(size(ket_mid,1).ne.n.or.size(ket_mid,2).ne.lot)then
        deallocate(ket_mid,stat=j); call check(j,0,' ket_mid, dealloc ')
     end if
  endif


  if(.not.allocated(diag_temp)) then
     allocate( diag_temp(n), off_diag_temp(n-1),stat=j);if(j/=0)stop
  end if

  if(.not.allocated(ket_mid)) then
     allocate( ket_mid(n,lot),stat=j);if(j/=0)stop
  endif

  !
  ! this subroutine calculates ket_new = 1/(1+h*dt/2)  * (1-h*dt/2)*ket
  !

  diag_temp =      1 - diag* dt/2.d0
  off_diag_temp = -off_diag* dt/2.d0
  
  ket_mid(1,:) = diag_temp(1)*ket(1,:) + off_diag_temp(1)*ket(2,:)
  do j=1+1,n-1
     ket_mid(j,:) =     diag_temp(j)  *    ket(j,:) + &
                    off_diag_temp(j-1)*    ket(j-1,:) + &
                    off_diag_temp(j)  *    ket(j+1,:) 

  enddo
  ket_mid(n,:) = diag_temp(n)*ket(n,:) + off_diag_temp(n-1)*ket(n-1,:)

  !  now multip. by 1/(1+H*dt/2)

  diag_temp =     1 + diag* dt/2.d0
  off_diag_temp = off_diag* dt/2.d0  ! change the minuses

  call c_symm_tridiag_lot(off_diag_temp, diag_temp, ket_mid, ket_new, n, lot)
  
end subroutine c_cranck_nicolson


subroutine c_cranck_nicolson_1(ket, diag, off_diag, dt, n, lot)  ! (1-ci*dtH/2)/(..)
  use cranck_nicolson_work
  implicit none
  integer n, lot, j
  complex*16 off_diag( n-1), diag(n), ket(n, lot)
  complex*16, parameter :: ci=(0.d0,1.d0)
  real*8  dt

  if(allocated(diag_temp)) then
     if(size(diag_temp).ne.n)then
        deallocate(diag_temp, off_diag_temp, stat=j); call check(j,0,'  deall diag_temp ')
     endif
  endif

  if(allocated(ket_mid)) then
     if(size(ket_mid,1).ne.n.or.size(ket_mid,2).ne.lot)then
        deallocate(ket_mid,stat=j); call check(j,0,' ket_mid, dealloc ')
     end if
  endif

  if(.not.allocated(diag_temp)) then
     allocate( diag_temp(n), off_diag_temp(n-1),stat=j);call check(j,0,' cnico b ')
  end if

  if(.not.allocated(ket_mid)) then
     allocate( ket_mid(n,lot),stat=j);call check(j,0,' cnico c ')
  endif

  !
  ! this subroutine calculates ket = 1/(1+h*dt/2)  * (1-h*dt/2)*ket
  !

  diag_temp =      1 - ci*diag* dt/2.d0
  off_diag_temp = -ci*off_diag* dt/2.d0
  
  ket_mid(1,:) = diag_temp(1)*ket(1,:) + off_diag_temp(1)*ket(2,:)
  do j=1+1,n-1
     ket_mid(j,:) =     diag_temp(j)  *    ket(j,:) + &
                    off_diag_temp(j-1)*    ket(j-1,:) + &
                    off_diag_temp(j)  *    ket(j+1,:) 

  enddo
  ket_mid(n,:) = diag_temp(n)*ket(n,:) + off_diag_temp(n-1)*ket(n-1,:)

  !  now multip. by 1/(1+H*dt/2)

  diag_temp =     1 + ci*diag* dt/2.d0
  off_diag_temp = ci*off_diag* dt/2.d0  ! change the minuses

  call c_symm_tridiag_lot(off_diag_temp, diag_temp, ket_mid, ket, n, lot)
  
end subroutine c_cranck_nicolson_1

