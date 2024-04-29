!!****f* m_kg/getph
!!
!! NAME
!! getph
!!
!! FUNCTION
!! Compute three factors of one-dimensional structure factor phase
!! for input atomic coordinates, for all planewaves which fit in fft box.
!! The storage of these atomic factors is made according to the
!! values provided by the index table atindx. This will save time in nonlop.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  natom=number of atoms in cell.
!!  n1,n2,n3=dimensions of fft box (ngfft(3)).
!!  xred(3,natom)=reduced atomic coordinates.
!!
!! OUTPUT
!!  ph1d(2,(2*n1+1)*natom+(2*n2+1)*natom+(2*n3+1)*natom)=exp(2Pi i G.xred) for
!!   integer vector G with components ranging from -nj <= G <= nj.
!!   Real and imag given in usual Fortran convention.
!!
!! PARENTS
!!      m_afterscfloop,m_berryphase_new,m_bethe_salpeter,m_cgtk,m_cut3d
!!      m_dfpt_looppert,m_epjdos,m_extraprho,m_fock,m_gkk,m_gstate
!!      m_hamiltonian,m_inwffil,m_nonlinear,m_orbmag,m_pead_nl_loop,m_phgamma
!!      m_phpi,m_prcref,m_respfn_driver,m_scfcv_core,m_screening_driver
!!      m_sigma_driver,m_sigmaph,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine getph(atindx,natom,n1,n2,n3,ph1d,xred,dim1,dim2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,natom
!arrays
 integer,intent(in) :: atindx(natom)
 real*8,intent(in) :: xred(3,natom)
 integer :: dim1, dim2
 real*8 :: ph1d(dim1,dim2)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ii,ph1d_size1,ph1d_size2,ph1d_sizemin
 !character(len=500) :: msg
 real*8 :: arg

! *************************************************************************
! write(979,*) atindx,natom,n1,n2,n3,xred

 ph1d_size1=size(ph1d,1);ph1d_size2=size(ph1d,2)
! write(*,*) ph1d_size2
 ph1d_sizemin=(2*n1+1+2*n2+1+2*n3+1)*natom
 if (ph1d_size1/=2.or.ph1d_size2<ph1d_sizemin) then
   stop 'Wrong ph1d sizes!'
 end if

 do ia=1,natom

!  Store the phase factor of atom number ia in place atindx(ia)
   i1=(atindx(ia)-1)*(2*n1+1)
   i2=(atindx(ia)-1)*(2*n2+1)+natom*(2*n1+1)
   i3=(atindx(ia)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)

   do ii=1,2*n1+1
     arg=6.28318530717959d0*dble(ii-1-n1)*xred(1,ia)
     ph1d(re,ii+i1)=dcos(arg)
     ph1d(im,ii+i1)=dsin(arg)
   end do

   do ii=1,2*n2+1
     arg=6.28318530717959d0*dble(ii-1-n2)*xred(2,ia)
     ph1d(re,ii+i2)=dcos(arg)
     ph1d(im,ii+i2)=dsin(arg)
   end do

   do ii=1,2*n3+1
     arg=6.28318530717959d0*dble(ii-1-n3)*xred(3,ia)
     ph1d(re,ii+i3)=dcos(arg)
     ph1d(im,ii+i3)=dsin(arg)
   end do

 end do

!This is to avoid uninitialized ph1d values
 if (ph1d_sizemin<ph1d_size2) then
   ph1d(:,ph1d_sizemin+1:ph1d_size2)=0d0
 end if

end subroutine getph

