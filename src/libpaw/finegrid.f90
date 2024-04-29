module m_paw_nhat
    implicit none

    contains
!!****f* m_paw_nhat/nhatgrid
!! NAME
!! nhatgrid
!!
!! FUNCTION
!! Determine parts of the rectangular (fine) grid that are contained
!! inside spheres around atoms (used to compute n_hat density).
!! If corresponding option is selected, compute also g_l(r)*Y_lm(r)
!! (and derivatives) on this grid (g_l=radial shape function).
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  distribfft<type(distribfft_type)>=--optional-- contains all the informations related
!!                                    to the FFT parallelism and plane sharing
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  comm_fft=--optional-- MPI communicator over FFT components
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nattyp(ntypat)= # atoms of each type.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in unit cell
!!  optcut= option for the cut-off radius of spheres:
!!          if optcut=0, cut-off radius=pawtab%rshp=cut-off radius of compensation charge
!!          if optcut=1, cut-off radius=pawtab%rpaw=radius of PAW augmentation regions
!!  optgr0= 1 if g_l(r)*Y_lm(r) are computed
!!  optgr1= 1 if first derivatives of g_l(r)*Y_lm(r) are computed
!!  optgr2= 1 if second derivatives of g_l(r)*Y_lm(r) are computed
!!  optrad= 1 if vectors (r-r_atom) on the fine grid around atoms have to be stored
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type (integer) for each atom
!!  typord=1 if the output is ordered by type of atoms, 0 otherwise
!!  ucvol=unit cell volume in bohr**3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  pawfgrtab(natom)%ifftsph(nfgd)=FFT index (fine grid) of a points in paw spheres around each atom
!!  pawfgrtab(natom)%nfgd= number of (fine grid) FFT points in paw spheres around atoms
!!  if (optgr0==1)
!!    pawfgrtab(natom)%gylm(nfgd,l_size**2)= g_l(r)*Y_lm(r) around each atom
!!  if (optgr1==1)
!!    pawfgrtab(natom)%gylmgr(3,nfgd,l_size**2)= derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optgr2==1)
!!    pawfgrtab(natom)%gylmgr2(6,nfgd,l_size**2)= second derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optrad==1)
!!    pawfgrtab(natom)%rfgd(3,nfgd)= coordinates of r-r_atom around each atom
!!
!! PARENTS
!!      m_afterscfloop,m_bethe_salpeter,m_classify_bands,m_exc_analyze
!!      m_nonlinear,m_paw_mkaewf,m_paw_mkrho,m_respfn_driver,m_scfcv_core
!!      m_screening_driver,m_sigma_driver,m_wfd,m_wfk_analyze
!!
!! CHILDREN
!!      pawgylm,pawrfgd_wvl,timab,xred2xcart
!!
!! SOURCE

subroutine nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfft,ntypat,optcut,&
        & optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,typat,ucvol,xred,&
        & mpi_atmtab,comm_atom,comm_fft,typord) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'nhatgrid'
!End of the abilint section
 use m_pawtab
 use m_pawfgrtab
 use m_libpaw_defs
 use m_libpaw_mpi
 use m_paral_atom
 use m_paw_finegrid
 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat,optcut,optgr0,optgr1,optgr2,optrad
 integer,optional,intent(in) :: comm_atom,comm_fft,typord
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 integer,intent(in),target :: atindx1(natom),nattyp(ntypat)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gmet(3,3),rprimd(3,3),xred(3,natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ------------------------------
!scalars
 integer :: i3,iat,iatm,iatom,iatom_,iatom_tot,itypat,lm_size,me_fft,my_comm_atom,n1,n2,n3,nfgd
 logical :: grid_found,my_atmtab_allocated,paral_atom
 real(dp) :: rcut
 character(len=500) :: msg
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 integer,pointer :: my_atindx1(:),my_atmtab(:),my_nattyp(:)
! integer, pointer :: fftn3_distrib(:),ffti3_local(:)
! integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: rfgd_tmp(:,:)

! write(*,*) 'flag 0'
! write(180,*) 'atindx1',atindx1
! write(180,*) 'gmet',gmet
! write(180,*) 'my_natom',my_natom
! write(180,*) 'natom',natom
! write(180,*) 'nattyp',nattyp
! write(180,*) 'ngfft',ngfft
! write(180,*) 'ntypat',ntypat
! write(180,*) 'optcut',optcut
! write(180,*) 'optgr0',optgr0
! write(180,*) 'optgr1',optgr1
! write(180,*) 'optgr2',optgr2 
! write(180,*) 'optrad',optrad
! write(180,*) 'rprimd',rprimd
! write(180,*) 'typat',typat
! write(180,*) 'ucvol',ucvol
! write(180,*) 'xred',xred

! *************************************************************************

! DBG_ENTER("COLL")

!write(*,*) 'flag 1'
! call timab(559,1,tsec)
 if (my_natom==0) return

!Set up parallelism over FFT
 me_fft=0
 if (present(comm_fft)) then
   me_fft=xpaw_mpi_comm_rank(comm_fft)
 end if
!write(*,*) 'flag 2'

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
!write(*,*) 'flag 2_1'
 nullify(my_atmtab);
!write(*,*) 'flag 2_2'
 if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
!write(*,*) 'flag 2_3'
 my_comm_atom=xpaw_mpi_comm_self
!write(*,*) 'flag 2_4'
! if (present(comm_atom)) my_comm_atom=comm_atom
!write(*,*) 'flag 2_5'

 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

!write(*,*) 'flag 3'
 if (paral_atom) then
   allocate(my_atindx1(natom))
   allocate(my_nattyp(ntypat))
   my_atindx1(:)=0;my_nattyp(:)=0
   iat=1
   do itypat=1,ntypat
     if (my_natom>0) then
       do iatom=1,my_natom
         if(typat(my_atmtab(iatom))==itypat)then
           my_nattyp(itypat)=my_nattyp(itypat)+1
           my_atindx1(iat)=iatom
           iat=iat+1
         end if
       end do
     end if
   end do
 else
   my_atindx1 => atindx1
   my_nattyp => nattyp
 end if
!write(*,*) 'flag 4'

!Get the distrib associated with this fft_grid
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
! if (present(distribfft)) then
!   grid_found=.false.
!   if (n2 == distribfft%n2_coarse) then
!     if (n3== size(distribfft%tab_fftdp3_distrib)) then
!       fftn3_distrib => distribfft%tab_fftdp3_distrib
!       ffti3_local => distribfft%tab_fftdp3_local
!       grid_found=.true.
!     end if
!   end if
!   if (n2 == distribfft%n2_fine) then
!     if (n3 == size(distribfft%tab_fftdp3dg_distrib)) then
!       fftn3_distrib => distribfft%tab_fftdp3dg_distrib
!       ffti3_local => distribfft%tab_fftdp3dg_local
!       grid_found = .true.
!     end if
!   end if
!   if (.not.(grid_found)) then
!     msg='Unable to find an allocated distrib for this fft grid!'
!     MSG_BUG(msg)
!   end if
! else
!   LIBPAW_ALLOCATE(fftn3_distrib,(n3))
!   LIBPAW_ALLOCATE(ffti3_local,(n3))
!   fftn3_distrib=0;ffti3_local=(/(i3,i3=1,n3)/)
! end if

!Loop over types of atom
!-------------------------------------------
 iatm=0
 do itypat=1,ntypat

   if (optcut==1) then
     rcut=pawtab(itypat)%rpaw
   else
     rcut=pawtab(itypat)%rshp
   end if
!write(*,*) 'flag 5'

!  Loop over atoms
!  -------------------------------------------
   do iat=1,my_nattyp(itypat)
     iatm=iatm+1;iatom=my_atindx1(iatm)
     iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
     iatom_=iatom!if(present(typord)) iatom_=merge(iatm,iatom,typord==1)
     lm_size=pawfgrtab(iatom_)%l_size**2

!    ------------------------------------------------------------------
!    A-Determine FFT points and r-R vectors around the atom
!    ------------------------------------------------------------------
!write(*,*) 'flag 6'

     call pawrfgd_fft(ifftsph_tmp,gmet,n1,n2,n3,nfgd,rcut,rfgd_tmp,rprimd,ucvol,&
&     xred(:,iatom_tot),me_fft=me_fft)
!     call pawrfgd_fft(ifftsph_tmp,gmet,n1,n2,n3,nfgd,rcut,rfgd_tmp,rprimd,ucvol,&
!&     xred(:,iatom_tot),fft_distrib=fftn3_distrib,fft_index=ffti3_local,me_fft=me_fft)

!    Allocate arrays defining sphere (and related data) around current atom
     if (allocated(pawfgrtab(iatom_)%ifftsph)) then
       deallocate(pawfgrtab(iatom_)%ifftsph)
     end if
     allocate(pawfgrtab(iatom_)%ifftsph(nfgd))
     pawfgrtab(iatom_)%nfgd=nfgd
!     write(180,*) 'iat,nfgd',nfgd
     pawfgrtab(iatom_)%ifftsph(1:nfgd)=ifftsph_tmp(1:nfgd)
!write(*,*) 'flag 7'

     if (optrad==1) then
       if (allocated(pawfgrtab(iatom_)%rfgd))  then
         deallocate(pawfgrtab(iatom_)%rfgd)
       end if
       allocate(pawfgrtab(iatom_)%rfgd(3,nfgd))
       pawfgrtab(iatom_)%rfgd_allocated=1
       pawfgrtab(iatom_)%rfgd(1:3,1:nfgd)=rfgd_tmp(1:3,1:nfgd)
     end if
!write(*,*) 'flag 8'

     if (optgr0==1) then
       if (allocated(pawfgrtab(iatom_)%gylm))  then
         deallocate(pawfgrtab(iatom_)%gylm)
       end if
       ALLOCATE(pawfgrtab(iatom_)%gylm(nfgd,lm_size))
       pawfgrtab(iatom_)%gylm_allocated=1
     end if
!write(*,*) 'flag 9'

     if (optgr1==1) then
       if (allocated(pawfgrtab(iatom_)%gylmgr))  then
         DEALLOCATE(pawfgrtab(iatom_)%gylmgr)
       end if
       ALLOCATE(pawfgrtab(iatom_)%gylmgr(3,nfgd,lm_size))
       pawfgrtab(iatom_)%gylmgr_allocated=1
     end if
!write(*,*) 'flag 10'

     if (optgr2==1) then
       if (allocated(pawfgrtab(iatom_)%gylmgr2))  then
         DEALLOCATE(pawfgrtab(iatom_)%gylmgr2)
       end if
       ALLOCATE(pawfgrtab(iatom_)%gylmgr2(6,nfgd,lm_size))
       pawfgrtab(iatom_)%gylmgr2_allocated=1
     end if
!write(*,*) 'flag 11'

!    ------------------------------------------------------------------
!    B-Calculate g_l(r-R)*Y_lm(r-R) for each r around the atom R
!    ------------------------------------------------------------------
!     if (optgr0+optgr1+optgr2>0) then
!       call pawgylm(pawfgrtab(iatom_)%gylm,pawfgrtab(iatom_)%gylmgr,pawfgrtab(iatom_)%gylmgr2,&
!&       lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),rfgd_tmp(:,1:nfgd))
!     end if

!    End loops over types/atoms
!    -------------------------------------------
     DEALLOCATE(ifftsph_tmp)
     DEALLOCATE(rfgd_tmp)
   end do
 end do
!write(*,*) 'flag 12'

!Destroy atom tables used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 if (paral_atom) then
   DEALLOCATE(my_atindx1)
   DEALLOCATE(my_nattyp)
 end if
!write(*,*) 'flag 13'

! if (.not.present(distribfft)) then
!   LIBPAW_DEALLOCATE(fftn3_distrib)
!   LIBPAW_DEALLOCATE(ffti3_local)
! end if

! call timab(559,2,tsec)

! DBG_EXIT("COLL")

end subroutine nhatgrid
end module m_paw_nhat
