subroutine invcb(rhoarr,rspts,npts)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invcb'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npts
!arrays
 real*8,intent(in) :: rhoarr(npts)
 real*8,intent(out) :: rspts(npts)

!Local variables-------------------------------
!scalars
 integer :: ii,ipts
 real*8,parameter :: c2_27=2.0d0/27.0d0,c5_9=5.0d0/9.0d0
 real*8,parameter :: c8_9=8.0d0/9.0d0,m1thrd=-1.0/3.0
 real*8 :: del,prod,rho,rhom1,rhomtrd
 logical :: test
!character(len=500) :: message

! *************************************************************************

!Loop over points : here, brute force algorithm
!do ipts=1,npts
!rspts(ipts)=sign( (abs(rhoarr(ipts)))**m1thrd,rhoarr(ipts))
!end do
!

 rhomtrd=sign( (abs(rhoarr(1)))**m1thrd, rhoarr(1) )
 rhom1=1d0/rhoarr(1)
 rspts(1)=rhomtrd
 do ipts=2,npts
   rho=rhoarr(ipts)
   prod=rho*rhom1
!  If the previous point is too far ...
   if(prod < 0.01d0 .or. prod > 10.0d0 )then
     rhomtrd=sign( (abs(rho))**m1thrd , rho )
     rhom1=1d0/rho
   else
     del=prod-1d0
     do ii=1,5
!      Choose one of the two next lines, the last one is more accurate
!      rhomtrd=((one+third*del)/(one+two_thirds*del))*rhomtrd
       rhomtrd=((1d0+c5_9*del)/(1d0+del*(c8_9+c2_27*del)))*rhomtrd
       rhom1=rhomtrd*rhomtrd*rhomtrd
       del=rho*rhom1-1d0
!      write(std_out,*)rhomtrd,del
       test = del*del < 1.0d-24
       if(test) exit
     end do
     if( .not. test) then
       rhomtrd=sign( (abs(rho))**m1thrd , rho )
     end if
   end if
   rspts(ipts)=rhomtrd
 end do

end subroutine invcb
!!***
