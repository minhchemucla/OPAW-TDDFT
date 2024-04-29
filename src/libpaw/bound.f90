subroutine bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section
 use m_libpaw_defs

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: plane
 real(dp),intent(out) :: dsqmax,dsqmin
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(out) :: gbound(3)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer :: i1,i1min,i2,i2min,i3,i3min
 real(dp) :: dsm,dsp
 character(len=500) :: message

! *************************************************************************

! dsq(i1,i2,i3)=gmet(1,1)*(kpt(1)+dble(i1))**2&
!& +gmet(2,2)*(kpt(2)+dble(i2))**2&
!& +gmet(3,3)*(kpt(3)+dble(i3))**2&
!& +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
!& +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
!& +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))

!Set plane to impossible value
 plane=0

!look at +/- g1 planes:
 dsqmax=zero
 dsqmin=dsq(ngfft(1)/2,-ngfft(2)/2,-ngfft(3)/2,gmet,kpt)+0.01_dp
 do i2=-ngfft(2)/2,ngfft(2)/2
   do i3=-ngfft(3)/2,ngfft(3)/2
     dsp = dsq(ngfft(1)/2, i2, i3,gmet,kpt)
     dsm = dsq( - ngfft(1)/2, i2, i3,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = ngfft(1)/2
       i2min = i2
       i3min = i3
       plane=1
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min =  - ngfft(1)/2
       i2min = i2
       i3min = i3
       plane=1
     end if
   end do
 end do
!
!+/- g2 planes:
 do i1=-ngfft(1)/2,ngfft(1)/2
   do i3=-ngfft(3)/2,ngfft(3)/2
     dsp = dsq(i1,ngfft(2)/2,i3,gmet,kpt)
     dsm = dsq(i1,-ngfft(2)/2,i3,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = i1
       i2min = ngfft(2)/2
       i3min = i3
       plane=2
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min = i1
       i2min =  - ngfft(2)/2
       i3min = i3
       plane=2
     end if
   end do
 end do
!
!+/- g3 planes:
 do i1=-ngfft(1)/2,ngfft(1)/2
   do i2=-ngfft(2)/2,ngfft(2)/2
     dsp = dsq(i1,i2,ngfft(3)/2,gmet,kpt)
     dsm = dsq(i1,i2,-ngfft(3)/2,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = i1
       i2min = i2
       i3min = ngfft(3)/2
       plane=3
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min = i1
       i2min = i2
       i3min =  - ngfft(3)/2
       plane=3
     end if
   end do
 end do

 if (plane==0) then
!  Trouble: missed boundary somehow
   write(message, '(a,a,a,3f9.4,a,3i5,a,a,a,a,a)' )&
&   'Trouble finding boundary of G sphere for',ch10,&
&   'kpt=',kpt(:),' and ng=',ngfft(1:3),ch10,&
&   'Action : check that kpt lies',&
&   'reasonably within first Brillouin zone; ',ch10,&
&   'else code bug, contact ABINIT group.'
   write(*,*) message
   stop
 end if

 gbound(1)=i1min
 gbound(2)=i2min
 gbound(3)=i3min

 contains

   function dsq(i1,i2,i3,gmet,kpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

     integer :: i1,i2,i3
     real(dp) :: dsq
     real(dp) :: kpt(3),gmet(3,3)

     dsq=gmet(1,1)*(kpt(1)+dble(i1))**2&
&      +gmet(2,2)*(kpt(2)+dble(i2))**2&
&      +gmet(3,3)*(kpt(3)+dble(i3))**2&
&      +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
&      +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
&      +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))
   end function dsq

end subroutine bound

