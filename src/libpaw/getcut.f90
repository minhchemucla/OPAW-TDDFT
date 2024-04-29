subroutine getcut(boxcut,ecut,gmet,gsqcut,iboxcut,iout,kpt,ngfft,ek_factor)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section
 use m_libpaw_defs
 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iboxcut,iout
 real(dp),intent(inout) :: ecut
 real(dp),intent(out) :: boxcut,gsqcut
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer :: plane
 real(dp) :: boxsq,cutrad,ecut_pw,effcut,largesq,sphsq,ek_factor
 character(len=500) :: message
!arrays
 integer :: gbound(3)
 real*8  :: ten=10d0
! *************************************************************************

!This is to treat the case where ecut has not been initialized,
!for wavelet computations. The default for ecut is -1.0 , allowed
!only for wavelets calculations
 ecut_pw=ecut
 if(ecut<-tol8)ecut_pw=ten

!gcut(box)**2=boxsq; gcut(sphere)**2=sphsq
!get min. d**2 to boundary of fft box:
!(gmet sets dimensions: bohr**-2)
!ecut(sphere)=0.5*(2 pi)**2 * sphsq:
 call bound(largesq,boxsq,gbound,gmet,kpt,ngfft,plane)
 effcut=0.5_dp * (two_pi)**2 * boxsq
 ecut=effcut*(ek_factor-1d-10)
 ecut_pw=ecut
 write(*,*) 'ecut adjusted to:',ecut
 sphsq=2._dp*ecut_pw/two_pi**2

 if (iboxcut/=0) then
   boxcut=10._dp
   gsqcut=(largesq/sphsq)*(2.0_dp*ecut)/two_pi**2

   write(message, '(a,a,3f8.4,a,3i4,a,a,f11.3,a,a)' ) ch10,&
&   ' getcut: wavevector=',kpt,'  ngfft=',ngfft(1:3),ch10,&
&   '         ecut(hartree)=',ecut_pw+tol8,ch10,'=> whole FFT box selected'
   if(iout/=std_out) then
     write(iout,*)message
   end if
   write(*,*) message

 else

!  Get G^2 cutoff for sphere of double radius of basis sphere
!  for selecting G s for rho(G), V_Hartree(G), and V_psp(G)--
!  cut off at fft box boundary or double basis sphere radius, whichever
!  is smaller.  If boxcut were 2, then relation would be
!$ecut_eff = (1/2) * (2 Pi Gsmall)^2 and gsqcut=4*Gsmall^2$.
   boxcut = sqrt(boxsq/sphsq)
   cutrad = min(2.0_dp,boxcut)
   gsqcut = (cutrad**2)*(2.0_dp*ecut_pw)/two_pi**2

   if(ecut>-tol8)then

     write(message, '(a,a,3f8.4,a,3i4,a,a,f11.3,3x,a,f10.5)' ) ch10,&
&     ' getcut: wavevector=',kpt,'  ngfft=',ngfft(1:3),ch10,&
&     '         ecut(hartree)=',ecut+tol8,'=> boxcut(ratio)=',boxcut+tol8
     if(iout/=std_out) then
       write(iout,*)message
     end if
     write(*,*) message

     if (boxcut<1.0_dp) then
       write(message, '(a,a,a,a,a,a,a,a,a,f12.6)' )&
&       '  Choice of acell, ngfft, and ecut',ch10,&
&       '  ===> basis sphere extends BEYOND fft box !',ch10,&
&       '  Recall that boxcut=Gcut(box)/Gcut(sphere)  must be > 1.',ch10,&
&       '  Actio: try larger ngfft or smaller ecut.',ch10,&
&       '  Note that ecut=effcut/boxcut**2 and effcut=',effcut+tol8
       if(iout/=std_out) then
           write(iout,*) message
       end if
       write(*,*) message
     end if

     if (boxcut>2.2_dp) then
       write(message, '(a,a,a,a,a,a,a,a,a,a,a,f12.6,a,a)' ) ch10,&
&       ' getcut : COMMENT -',ch10,&
&       '  Note that boxcut > 2.2 ; recall that',' boxcut=Gcut(box)/Gcut(sphere) = 2',ch10,&
&       '  is sufficient for exact treatment of convolution.',ch10,&
&       '  Such a large boxcut is a waste : you could raise ecut',ch10,&
&       '  e.g. ecut=',effcut*0.25_dp+tol8,' Hartrees makes boxcut=2',ch10
       if(iout/=std_out) then
           write(iout,*) message
       end if
       write(*,*) message
     end if

     if (boxcut<1.5_dp) then
       write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&       ' getcut : WARNING -',ch10,&
&       '  Note that boxcut < 1.5; this usually means',ch10,&
&       '  that the forces are being fairly strongly affected by','  the smallness of the fft box.',ch10,&
&       '  Be sure to test with larger ngfft(1:3) values.',ch10
       if(iout/=std_out) then
           write(iout,*) message
       end if
    write(*,*) message  
   end if

   end if

 end if  ! iboxcut

 write(74,*) 'gsqcut',gsqcut
end subroutine getcut

