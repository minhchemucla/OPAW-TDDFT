!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_atm2fft
!! NAME
!!  m_atm2fft
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2018 ABINIT group (FJ, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#include "libpaw.h"

module m_atm2fft

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_pawtab,      only : pawtab_type
! use m_fft,         only : zerosym, fourdp

 implicit none

 private
!!***

 public :: atm2fft
!!***

contains
!!***

!!****f* ABINIT/atm2fft
!! NAME
!! atm2fft
!!
!! FUNCTION
!! This routine sums atomic functions (density or potential) defined
!! (in rec. space) on a radial grid to get global quantities on the
!! fine FFT grid. It can also compute contribution to energy derivatives
!! of these atomic functions.
!!
!! Possible options:
!!   optn=1: compute a sum of local atomic densities or contrib. to energy derivatives
!!   optv=1: compute a sum of local atomic potentials or contrib. to energy derivatives
!!
!!   optatm =1: computes sum of atomic potentials/densities
!!   optgr  =1: computes contribution of atomic pot./dens. to forces
!!   optstr =1: computes contribution of atomic pot./dens. to stress tensor
!!   optdyfr=1: computes contribution of atomic pot./dens. to frozen part of dyn. matrix
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  distribfft<type(distribfft_type)>=--optional-- contains infos related to FFT parallelism
!!  gauss(2,ntypat)= params for gaussian atm density (optn2=3) for each atom type
!!  gmet(3,3)=reciprocal space metric
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on |G|^2: see setup1 for definition (doubled sphere)
!!  kxc(2,nfft)=exchange and correlation kernel
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when comm_fft is present)
!!  mgfft=maximum size of 1D FFTs
!!  comm_fft=--optional-- MPI communicator over FFT components
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  ntypat=number of types of atoms.
!!  optatm,optdyfr,optgr,optn,optn2,optstr,optv= (see NOTES below)
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qprtrb(3)= integer wavevector of possible perturbing potential
!!             in basis of reciprocal lattice translations
!!  rhog(2,nfft)=electron density rho(G) in reciprocal space
!!               (used only if optv=1 and (optgr=1 or optstr=1 or optdyfr=1))
!!  ucvol=unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vspl(mqgrid,2,ntypat)=q^2 v(q) spline of an atomic potential
!!                        (used only if optv=1)
!!  vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!            perturbing potential is added of the form V(G)=(vprtrb(1)+I*vprtrb(2))/2
!!            at the values G=qprtrb and (vprtrb(1)-I*vprtrb(2))/2 at G=-qprtrb
!!  vg(2,nfft)= potential V(G) in reciprocal space
!!              (used only if optn=1 and (optgr=1 or optstr=1 or optdyfr=1))
!!  vg1(2,nfft)= 1st-order potential V(G) in reciprocal space
!!               (used only if opteltfr==1)
!!  vg1_core(2,nfft)= 1-order potential V(G) in reciprocal space with only core contribution
!!                    (used only if opteltfr==1)
!! OUTPUT
!!  ======= if optv==1 =======
!!  ============================
!!   --- if optatm==1
!!    atmvloc(nfft)=sum of local atomic potentials in real space
!!   --- if optgr==1
!!    grv(3,natom)=contribution of atomic potentials to forces
!!   --- if optstr==1
!!    strv(6)=contribution of atomic potentials to stress tensor
!!            cart. coordinates, symmetric tensor, 6 comp. in order 11,22,33,32,31,21
!!   --- if optdyfr==1
!!    dyfrv(3,3,natom)=contribution of atomic potentials to frozen part of dyn. matrix
!!
!!  ======= if optn==1 =======
!!  ============================
!!   --- if optatm==1
!!    atmrho(nfft)=sum of atomic densities in real space
!!   --- if optgr==1
!!    grn(3,natom)=contribution of atomic densities to forces
!!   --- if optstr==1
!!    strn(6)=contribution of atomic densities to stress tensor
!!            cart. coordinates, symmetric tensor, 6 comp. in order 11,22,33,32,31,21
!!   --- if optdyfr==1
!!    dyfrn(3,3,natom)=contribution of atomic densities to frozen part of dyn. matrix
!!   --- if opteltfr==1
!!    eltfrn(6+3*natom,6)=contribution of atomic density to frozen part of stress tensor
!!
!! NOTES
!! Details on possible options:
!! ============================
!! optv: controls the computation of a local potential as sum of atomic potentials
!!          Vloc(r)=Sum_R[V^AT(r-R)]
!!       or its contributions to energy derivatives, i.e. derivatives of Int[Vloc(r).rho(r).dr]
!!          rho(r) is stored in reciprocal space in array rhog()
!!          V^AT is stored in reciprocal space in array vspl (in practice vspl(q)=q^2.V^AT(q))
!!
!! optn: controls the computation of a density as sum of atomic densities
!!          n(r)=Sum_R[n^AT(r-R)]
!!       or its contributions to energy derivatives, i.e. derivatives of Int[n(r).V(r).dr]
!!          V(r) is stored in reciprocal space in array vg()
!!          n^AT is stored in reciprocal space:
!!          if optn2=1: n^AT is the atomic PAW PS core density stored in array pawtab%tcorespl()
!!                   2: n^AT is the atomic PAW PS valence density stored in array pawtab%tvalespl()
!!                   3: n^AT is a gaussian density: n(g)=gauss(1,ityp)*exp[-(gauss(2,ityp)*G)^2]
!! Note: optv and optn can be activated together
!!
!! Options controlling which contrib. to Etot derivatives are computed:
!!   optatm  =1: computes Vloc(r) or n(r) as sum of atomic potentials/densities
!!   optgr   =1: computes contribution of atomic Vloc(r) or n(r) to forces
!!   optstr  =1: computes contribution of atomic Vloc(r) or n(r) to stress tensor
!!   optdyfr =1: computes contribution of atomic Vloc(r) or n(r) to fr part of dyn. matrix
!!   opteltfr=1: computes contribution of atomic Vloc(r) or n(r) to elastic tensor
!! Note: optatm, optgr, optstr, optelfr and optdyfr can be activated together
!!
!! Typical uses:
!! =============
!! Computation of:
!!  - local potential: optv=1, optatm=1
!!  - contrib. of local potential to Etot derivatives: optv=1, rhog=total valence density
!!                                                     optgr=1 or optstr=1 or optdyfr=1
!!  - PS core density: optn=1, optn2=1, optatm=1
!!  - contrib. of NLCC to Etot derivatives: optn=1, optn2=1, vg=XC potential
!!                                          optgr=1 or optstr=1 or optdyfr=1
!!  - sum of atomic valence densities: optn=1, optn2=2 or 3, optatm=1
!!  - correction of forces due to potential residual: optn=1, optn2=2 or 3, optgr=1
!!                                                    vg=potential residual
!!    etc...
!!
!! PARENTS
!!      m_dfpt_elt,m_extraprho,m_forces,m_nonlinear,m_prcref,m_respfn_driver
!!      m_setvtr,m_stress
!!
!! CHILDREN
!!      destroy_distribfft,fourdp,init_distribfft_seq,initmpi_seq
!!      set_mpi_enreg_fft,unset_mpi_enreg_fft,zerosym
!!
!! SOURCE

subroutine atm2fft(atindx1,atmrho,atmvloc,gmet,gprimd,&
&                  gsqcut,mgfft,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&
&                  pawtab,ph1d,qgrid,ucvol,vspl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,natom,nfft,ntypat
 real(dp),intent(in) :: gsqcut,ucvol
! type(distribfft_type),optional,intent(in),target :: distribfft
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(3)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid)
 real(dp),intent(in) :: vspl(mqgrid,2,ntypat)
 real(dp),intent(out) :: atmrho(nfft)
 real(dp),intent(inout) :: atmvloc(nfft)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

 complex*16 :: cin(nfft),cout(nfft)
!Local variables ------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ierr,ig1,ig1_,ig2,ig2_,ig3,ig3_,ii,is1,is2
 integer :: itypat,jj,js,ka,kb,kd,kg,ndir,n1,n2,n3,nproc_fft,paral_kgb_fft
 integer :: shift1,shift2,shift3
 logical :: have_g0
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,alf2pi2,bb,cc,cutoff,dbl_ig1,dbl_ig2,dbl_ig3,dd,dg1,dg2,d2g,diff
 real(dp) :: dn_at,d2n_at,d2n_at2,dq,dq2div6,dqdiv6,dqm1,dv_at,ee,ff,gauss1,gauss2,gg,gmag,gsquar,n_at
 real(dp) :: ph12i,ph12r,ph1i,ph1r,ph2i,ph2r,ph3i,ph3r,sfi,sfr,term,term1,term2,tmpni,tmpnr
 real(dp) :: tmpvi,tmpvr,v_at,xnorm
 character(len=500) :: message
! type(distribfft_type),pointer :: my_distribfft
! type(mpi_type) :: mpi_enreg_fft
!arrays
! integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 real(dp), CONTIGUOUS, pointer :: tvalespl(:,:),tcorespl(:,:)
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer  :: delta(6)=(/1,1,1,0,0,0/)
 real(dp) :: dgm(3,3,6),d2gm(3,3,6,6),gcart(3),tsec(2)
 real(dp),allocatable :: phim_igia(:),phre_igia(:),workn(:,:)
 real(dp),allocatable :: workv(:,:)

! *************************************************************************
! write(982,*) vspl
! write(984,*) qgrid
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Zero out arrays to permit accumulation over atom types
 ALLOCATE(workv(2,nfft))
 workv(:,:)=zero
 ALLOCATE(workn(2,nfft))
 workn(:,:)=zero

 dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
 dqm1=1.0_dp/dq
 dqdiv6=dq/6.0_dp
 dq2div6=dq**2/6.0_dp
 cutoff=gsqcut*tolfix
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 ALLOCATE(phre_igia(natom))
 ALLOCATE(phim_igia(natom))

 ia1=1
 do itypat=1,ntypat
!  ia1,ia2 sets range of loop over atoms:
   ia2=ia1+nattyp(itypat)-1
   ii=0

   tcorespl => pawtab(itypat)%tcorespl
   tvalespl => pawtab(itypat)%tvalespl

!   write(10030+itypat,*) tcorespl
!   write(10050+itypat,*) tvalespl

   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     ig3_=ig3;if (ig3_==(n3/2+1)) ig3_=0
     do i2=1,n2
       ig2=i2-(i2/id2)*n2-1
       ig2_=ig2;if (ig2_==(n2/2+1)) ig2_=0
         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1
           ig1_=ig1;if (ig1_==(n1/2+1)) ig1_=0
           ii=ii+1
           gsquar=gsq_atm(ig1,ig2,ig3)

!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then

             gmag=sqrt(gsquar)
             have_g0=(ig1==0.and.ig2==0.and.ig3==0)

             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Compute structure factor for all atoms of given type:
             do ia=ia1,ia2
               shift1=1+n1+(ia-1)*(2*n1+1)
               shift2=1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1)
               shift3=1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
               ph1r=ph1d(1,ig1+shift1);ph1i=ph1d(2,ig1+shift1)
               ph2r=ph1d(1,ig2+shift2);ph2i=ph1d(2,ig2+shift2)
               ph3r=ph1d(1,ig3+shift3);ph3i=ph1d(2,ig3+shift3)
               ph12r=ph1r*ph2r-ph1i*ph2i
               ph12i=ph1r*ph2i+ph1i*ph2r
               phre_igia(ia)=ph12r*ph3r-ph12i*ph3i
               phim_igia(ia)=ph12r*ph3i+ph12i*ph3r
             end do

!            Assemble structure factors for this type of atom= sum[exp(-i.piG.R)]
             sfr=zero;sfi=zero
             do ia=ia1,ia2
               sfr=sfr+phre_igia(ia)
               sfi=sfi-phim_igia(ia)
             end do

!            Compute V^AT(G) and/or n^AT(G) for given type of atom
!            Evaluate spline fit: p. 86 Numerical Recipes, Press et al;
!            NOTE: error in book for sign of "aa" term in derivative;
!            !           also see splfit routine.
             bb = diff*dqm1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*dq2div6
             dd = bb*(bb**2-1.0_dp)*dq2div6
             if (have_g0) then
               v_at=zero
             else
               v_at=(aa*vspl(jj,1,itypat)+bb*vspl(jj+1,1,itypat)+&
&               cc*vspl(jj,2,itypat)+dd*vspl(jj+1,2,itypat)) &
&               /gsquar
             end if
             n_at=(aa*tcorespl(jj,1)+bb*tcorespl(jj+1,1)+cc*tcorespl(jj,2)+dd*tcorespl(jj+1,2))

!            Compute sum of local atomic potentials or densities
!            ---------------------------------------------------
!              Accumulate V^AT(G)*SF(G) or n^AT(G)*SF(G)
             workv(re,ii)=workv(re,ii)+sfr*v_at
             workv(im,ii)=workv(im,ii)+sfi*v_at
             workn(re,ii)=workn(re,ii)+sfr*n_at
             workn(im,ii)=workn(im,ii)+sfi*n_at


           end if
!          End loop on n1, n2, n3
         end do
     end do
   end do

   ia1=ia2+1

!  End loop on type of atoms
 end do

 DEALLOCATE(phre_igia)
 DEALLOCATE(phim_igia)

 xnorm=one/ucvol

! write(990,*) workv(re,:)
! write(993,*) workv(im,:)
 cin=cmplx(workv(re,:),workv(im,:))
 call fft3d_backward(n1,n2,n3,cin,cout)
 atmvloc=dble(cout)
! write(992,*) atmvloc
 atmvloc=atmvloc*xnorm 

 open(file='myncoret_k',unit=994)
 write(994,*) workn(re,:) 
 close(994)
 
 cin=cmplx(workn(re,:),workn(im,:))
 call fft3d_backward(n1,n2,n3,cin,cout)
 atmrho=dble(cout)
 atmrho=atmrho*xnorm
! write(995,*) atmrho
 deallocate(workn,workv)
 contains

   function gsq_atm(i1,i2,i3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsq_atm'
!End of the abilint section

   real(dp) :: gsq_atm
   integer,intent(in) :: i1,i2,i3
   gsq_atm=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+dble(i3*i3)*gmet(3,3) &
&   +two*(dble(i1*i2)*gmet(1,2)+dble(i2*i3)*gmet(2,3)+dble(i3*i1)*gmet(3,1))
 end function gsq_atm

end subroutine atm2fft
!!***
end module m_atm2fft
!!***
