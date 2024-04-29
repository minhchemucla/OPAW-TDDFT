!!****f* m_paw_occupancies/initrhoij
!! NAME
!! initrhoij
!!
!! FUNCTION
!! Initialize PAW rhoij occupancies (in packed storage)
!! from atomic ones
!!
!! INPUTS
!!  cplex=1 if rhoij are REAL, 2 if they are complex
!!  lexexch(ntypat)=l on which local exact-exchange is applied for a given type of atom
!!  lpawu(ntypat)=l on which U is applied for a given type of atom (PAW+U)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of atom types
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated in PAW augmentation regions
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!                                     (containing initial rhoij)
!!  spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!!  typat(natom)=type of each atom
!!  === Optional arguments
!!    ngrhoij=number of gradients to be allocated (OPTIONAL, default=0)
!!    nlmnmix=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
!!    use_rhoij_=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
!!    use_rhoijres=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)

!!
!! OUTPUT
!!  pawrhoij(natom) <type(pawrhoij_type)>=rhoij quantities for each atom
!!                                        in packed storage
!!
!! PARENTS
!!      m_gstate,m_nonlinear,m_positron,m_respfn_driver
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawrhoij_alloc
!!
!! SOURCE




subroutine initrhoij(cplex,lexexch,lpawu,my_natom,natom,&
&                   nspden,nspinor,nsppol,ntypat,pawrhoij,pawspnorb,pawtab,spinat,typat,&
                    ngrhoij,nlmnmix,use_rhoij_,use_rhoijres,& ! optional arguments
                    comm_atom) ! optional arguments (parallelism)
!                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!#undef ABI_FUNC
!#define ABI_FUNC 'initrhoij'
!End of the abilint section
 use m_pawtab
 use m_libpaw_defs
 use m_libpaw_mpi
 use m_paral_atom
 use m_pawrhoij
 use m_libpaw_tools, only : libpaw_msg_hndl

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,my_natom,natom,nspden,nspinor,nsppol,ntypat,pawspnorb
 integer,intent(in),optional :: comm_atom,ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
 character(len=500) :: message
!arrays
 integer,intent(in) :: lexexch(ntypat),lpawu(ntypat)
 integer,intent(in) :: typat(natom)
! integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: spinat(3,natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!Arrays
!scalars
 integer :: iatom,iatom_rhoij,ilmn,ispden,itypat,j0lmn,jl,jlmn,jspden,klmn,klmn1,ln,lnspinat0,my_comm_atom
 integer :: ngrhoij0,nlmnmix0,nselect,nselect1,nspden_rhoij,use_rhoij_0,use_rhoijres0
 real(dp) :: ratio,ro,roshift,zratio,zz
 logical :: my_atmtab_allocated,paral_atom,spinat_zero,test_exexch,test_pawu,test_lnspinat
!arrays
 integer,pointer :: my_atmtab(:),lnspinat(:)
 real(dp),allocatable :: occ(:)
!************************************************************************
 
!write(*,*) 'flag 1'
! write(174,*) 'cplex',cplex
! write(174,*) 'natom',natom
! write(174,*) 'my_natom',my_natom
! write(174,*) 'nspden',nspden
! write(174,*) 'nspinor',nspinor
! write(174,*) 'nsppol',nsppol
! write(174,*) 'ntypat',ntypat
! write(174,*) 'pawspnorb',pawspnorb
! write(174,*) 'lexexch',lexexch
! write(174,*) 'lpawu',lpawu
! write(174,*) 'typat',typat
 write(174,*) 'spinat',spinat
! DBG_ENTER("COLL")
!write(*,*) 'flag 2'

!PAW+U and local exact-exchange restriction
 do itypat=1,ntypat
   if (lpawu(itypat)/=lexexch(itypat).and. lpawu(itypat)/=-1.and.lexexch(itypat)/=-1) then
     message = ' lpawu must be equal to lexexch !'
     call libpaw_msg_hndl(message,'ERROR','PERS') 
   end if
 end do
!write(*,*) 'flag 3'

!Set up parallelism over atoms
! paral_atom=(present(comm_atom).and.(my_natom/=natom))
! nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
! my_comm_atom=xpaw_mpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
! call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)
!write(*,*) 'flag 4'

 nspden_rhoij=pawrhoij_get_nspden(nspden,nspinor,pawspnorb)
 ratio=one;if (nspden_rhoij==2) ratio=half
 spinat_zero=all(abs(spinat(:,:))<tol10)
 write(*,*) 'spinatzero',spinat_zero
!write(*,*) 'flag 5'

 if (my_natom>0) then
!   ngrhoij0=0;if (present(ngrhoij)) ngrhoij0=ngrhoij
!   nlmnmix0=0;if (present(nlmnmix)) nlmnmix0=nlmnmix
!   use_rhoij_0=0;if (present(use_rhoij_)) use_rhoij_0=use_rhoij_
!   use_rhoijres0=0;if (present(use_rhoijres)) use_rhoijres0=use_rhoijres
!   if (paral_atom) then
!     call pawrhoij_alloc(pawrhoij,cplex,nspden_rhoij,nspinor,nsppol,typat,&
!&     ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0,&
!&     pawtab=pawtab,comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
!   else
     call pawrhoij_alloc(pawrhoij,cplex,nspden_rhoij,nspinor,nsppol,typat,pawtab=pawtab,&
&     ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0)
!   end if
 end if
!write(*,*) 'flag 6'

 do iatom_rhoij=1,my_natom
   iatom=iatom_rhoij;!if (paral_atom) iatom=my_atmtab(iatom_rhoij)
   itypat=typat(iatom)
   nselect=0
   allocate(lnspinat(pawtab(itypat)%basis_size))
   lnspinat=-1
!write(*,*) 'flag 7'
! Determine occupancies of each orbital
   if (nspden_rhoij==2) then
     ALLOCATE(occ(pawtab(itypat)%basis_size))
     occ=zero
     do jlmn=1,pawtab(itypat)%lmn_size
       ln=pawtab(itypat)%indlmn(5,jlmn)
       klmn=jlmn*(jlmn+1)/2
       occ(ln)=occ(ln)+pawtab(itypat)%rhoij0(klmn)
     end do
     do ln=1,pawtab(itypat)%basis_size
       if(pawtab(itypat)%orbitals(ln)==0.and.occ(ln)==1) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==1.and.(occ(ln)>=1.and.occ(ln)<=5)) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==2.and.(occ(ln)>=1.and.occ(ln)<=9)) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==3.and.(occ(ln)>=1.and.occ(ln)<=13)) lnspinat(ln)=ln
     end do
     DEALLOCATE(occ)
   end if
   lnspinat0=maxval(lnspinat)
   lnspinat0=-1
!write(*,*) 'flag 8'

!  Determine Z (trace of rhoij0 or part of it)
   zz=zero
   do jlmn=1,pawtab(itypat)%lmn_size
     jl=pawtab(itypat)%indlmn(1,jlmn)
     ln=pawtab(itypat)%indlmn(5,jlmn)
     j0lmn=jlmn*(jlmn-1)/2
     test_lnspinat=(lnspinat0==-1.or.lnspinat(ln)==ln)
     test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
     test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       if ((ilmn==jlmn).and.test_pawu.and.test_exexch.and.test_lnspinat) &
&       zz=zz+pawtab(itypat)%rhoij0(klmn)
     end do
   end do
!write(*,*) 'flag 9'

!  Compute rhoij from tabulated value and magnetization
   do ispden=1,nspden_rhoij

     zratio=zero
     roshift=one
     ratio=one
     if (nspden_rhoij==2) then
       ratio=half
       if ((spinat(3,iatom)>zero.and.ispden==1).or.&
&       (spinat(3,iatom)<zero.and.ispden==2)) then
         if(abs(zz)>tol12)then
           zratio=two*abs(spinat(3,iatom))/zz
         else
           zratio=zero
         end if
       end if
     else if (nspden_rhoij==4.and.ispden>=2) then
       roshift=zero
       if(abs(zz)>tol12)then
         zratio=spinat(ispden-1,iatom)/zz
       else
         zratio=zero
       end if
     end if
!write(*,*) 'flag 10'

     nselect=0;nselect1=1-cplex
     do jlmn=1,pawtab(itypat)%lmn_size
       jl=pawtab(itypat)%indlmn(1,jlmn)
       ln=pawtab(itypat)%indlmn(5,jlmn)
       j0lmn=jlmn*(jlmn-1)/2
       test_lnspinat=(lnspinat0==-1.or.lnspinat(ln)==ln)
       test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
       test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         ro=pawtab(itypat)%rhoij0(klmn)
         if ((ilmn==jlmn).and.test_pawu.and.test_exexch.and.test_lnspinat) then
           ro=ro*ratio*(roshift+zratio)
         else
           ro=ro*ratio*roshift
         end if
!write(*,*) 'flag 11'

         klmn1=cplex*(klmn-1)+1
         if (abs(ro)>tol10) then
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=ro
         else
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=zero
         end if

         if (ispden==nspden_rhoij) then
           if (any(abs(pawrhoij(iatom_rhoij)%rhoijp(klmn1,:))>tol10)) then
             nselect=nselect+1;nselect1=nselect1+cplex
             pawrhoij(iatom_rhoij)%rhoijselect(nselect)=klmn
             do jspden=1,nspden_rhoij
               pawrhoij(iatom_rhoij)%rhoijp(nselect1,jspden)=pawrhoij(iatom_rhoij)%rhoijp(klmn1,jspden)
             end do
           end if
         end if
!write(*,*) 'flag 12'

       end do
     end do

   end do
   pawrhoij(iatom_rhoij)%nrhoijsel=nselect
!write(*,*) 'flag 13'

!  Non-collinear magnetism: avoid zero magnetization, because it produces numerical instabilities
!    Add a small real to the magnetization ; not yet activated => must be tested.
!   if (pawrhoij(iatom_rhoij)%nspden==4.and.spinat_zero) then
!     pawrhoij(iatom_rhoij)%rhoijp(:,4)=pawrhoij(iatom_rhoij)%rhoijp(:,4)+tol10
!   end if
   DEALLOCATE(lnspinat)
 end do ! iatom_rhoij

!Destroy atom table used for parallelism
! call free_my_atmtab(my_atmtab,my_atmtab_allocated)
!write(*,*) 'flag 14'

! DBG_EXIT("COLL")

end subroutine initrhoij
