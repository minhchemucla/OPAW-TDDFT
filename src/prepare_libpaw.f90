subroutine prep_libpaw
    use param_mod, only : ek_factor
    use main_mod
    use paw_mod
    use libpaw_mod
    use atom_mod
    use mpi_lib_ours, only : rank, bcast_scalar_r8
    
    implicit none

    integer :: ia,it

!    allocate(vxc_r(nn,1),vtrial_r(nn))
    call set_alloc_libpaw
    call set_param_alloc

    if (rank==0) write(*,*) 'setting up libpaw'
    
    do it=1,ntype
!        if(rank==0) then
            filename=pawinfo(it)%filename
        
            if(rank==0) write(*,*) 'it,filename',it,filename

        !read information from xml file and calculate some stuff
            call rdpawpsxml(filename,pawsetup)
            call rdpawpsxml(filename,paw_setuploc)
            call pawpsp_read_header_xml(lloc,lmax,pspcod,pspxc,&
                & pawsetup,r2well,zion,znucl(it))
            call pawpsp_read_pawheader(pawpsp_header%basis_size,&
                &   lmax,pawpsp_header%lmn_size,&
                &   pawpsp_header%l_size,pawpsp_header%mesh_size,&
                &   pawpsp_header%pawver,pawsetup,&
                &   pawpsp_header%rpaw,pawpsp_header%rshp,pawpsp_header%shape_type)
                

            pawtab(it)%has_vhnzc=1
            pawtab(it)%has_vhtnzc=1
!    call libxc_functionals_init(ixc,nspden)
!    call print_input

            call pawpsp_17in(epsatm,ffspl,icoulomb,ipsp,ixc,lmax,&
                &     lnmax,pawpsp_header%mesh_size,mqgrid_ff,mqgrid_vl,pawpsp_header,pawrad(it),pawtab(it),&
                &     xcdev,qgrid_ff,qgrid_vl,usewvl,usexcnhat,vlspl(:,:,it),xcccrc,&
                &     xclevel,denpos,zion,znucl(it))
        
            call paw_setup_free(pawsetup)
            call paw_setup_free(paw_setuploc)    
!        endif
!        call pawrad_bcast(pawrad(it),0)
!        call pawtab_bcast(pawtab(it),0,only_from_file=.true.)
   enddo        
   

   do ia=1,natom
        it=atom_map(ia)
        l_size_atm(ia)=pawtab(it)%l_size
   enddo

    mpsang=1+lmax

    !initialize data structures and calculate some stuff
    call pawinit(gnt_option,gsqcut,hyb_range_fock,lcutdens,lmix,mpsang,nphi,nsym,ntheta,&
        &     pawang,pawrad,pawspnorb,pawtab,xcdev,xclevel,usepotzero)

    call pawfgrtab_init(pawfgrtab,cplex,l_size_atm,nspden,typat)

    call paw_ij_init(paw_ij,cplex,nspinor,nsppol,nspden,pawspnorb,natom,ntypat,typat,pawtab,&
        has_dij=1,has_dijhartree=1,has_dijso=1,has_pawu_occ=1,has_exexch_pot=1)
    call paw_an_init(paw_an,natom,ntypat,0,0,nspden,cplex,xcdev,typat,pawang,pawtab,&
        has_vxc=1,has_vxc_ex=1)
    call initrhoij(cplex,lexexch,lpawu,natom,natom,nspden,nspinor,&
        nsppol,ntypat,pawrhoij,pawspnorb,pawtab,spinat,typat)

    call nhatgrid(atindx1,gmet,natom,natom,nattyp,ngfft,ntypat,&
        &optcut=0,optgr0=0,optgr1=0,optgr2=0,optrad=1,&
        pawfgrtab=pawfgrtab,pawtab=pawtab,rprimd=rprimd,typat=typat,&
        ucvol=ucvol,xred=xred)


    if (rank==0) then
      write(*,*) 'finish preparing libpaw'
      write(*,*) '========================'
      write(*,*)
     endif
!    write(170,*) 'pawrhoij%rhoijp',pawrhoij(1)%rhoijp
!    call prep_rhoij(natom,pawrhoij,pawtab)

!    write(2201,*) 'pawrhoij%nrhoijsel',pawrhoij(1)%nrhoijsel
!    write(2201,*) 'pawrhoij%rhoijselect',pawrhoij(1)%rhoijselect

contains

    subroutine set_alloc_libpaw
        implicit none

        integer :: st

        allocate(pawrad(ntype),pawtab(ntype),pawrhoij(natom),paw_ij(natom),paw_an(natom),pawfgrtab(natom),typat(natom),stat=st)
        if(st/=0) stop 'libpaw alloc problem 1'
        allocate(znucl(ntype),nattyp(ntype),lexexch(ntype),lpawu(ntype),l_size_atm(natom),stat=st)
        if(st/=0) stop 'libpaw alloc problem 2'
        allocate(xred(3,natom),stat=st)
        if(st/=0) stop 'libpaw alloc problem 3'

        typat=atom_map
        nattyp=natom_type
        znucl=diff_z
        lexexch=-1
        lpawu=-1

    end subroutine set_alloc_libpaw

    subroutine set_param_alloc
        implicit none

        integer :: ia
        real*8  :: boxcut,effcut
        real*8  :: kpt(3)=0d0
        real*8  :: ek_target, ek_high, ek_low

        if(funct==0) then
            ixc=7
            xclevel=1
        else if(funct==1) then
            ixc=11
            xclevel=2
        else
            write(*,*) 'funct should be 0 or 1'
            stop
        endif
        rprimd=0d0
        rprimd(1,1)=xmax*2d0
        rprimd(2,2)=ymax*2d0
        rprimd(3,3)=zmax*2d0

        gmet=0d0
        gmet(1,1)=1d0/(xmax*2d0)**2
        gmet(2,2)=1d0/(ymax*2d0)**2
        gmet(3,3)=1d0/(zmax*2d0)**2

        ngfft=0
        ngfft(1)=nx
        ngfft(2)=ny
        ngfft(3)=nz
        ngfft(4)=2*(ngfft(1)/2)+1
        ngfft(5)=2*(ngfft(2)/2)+1
        ngfft(6)=ngfft(3)
        ngfft(7)=112
        ngfft(8)=16
        ngfft(9)=1
        ngfft(10)=1
        ngfft(12)=ngfft(2)
        ngfft(13)=ngfft(3)

        nfft=nn
        nfftot=nn

        ucvol=xmax*ymax*zmax*8d0
        
        spinat=0d0
        nucdipmom=0d0

        ntypat=ntype

        allocate(ffspl(mqgrid_ff,2,lnmax),vlspl(mqgrid_vl,2,ntype))
        allocate(qgrid_ff(mqgrid_ff),qgrid_vl(mqgrid_vl))

    !ekcut = 100d0
    if(rank==0) then
        if(ekread) then
          ek_target= ekcut
          ek_low = 0d0
          ek_high = 2d0*ek_factor
          do while (ek_high .gt. ek_low) 
            call getcut(boxcut,ecutpaw,gmet,gsqcut,0,74,kpt,ngfft,ek_factor)
            gsqcutdg=gsqcut
            ekcut=ecutpaw
            if ((ekcut-ek_target)> 0.5e-3) then
              !ek_low = (ek_high+ek_low)/2d0
              ek_high = ek_factor
            else if ((ekcut-ek_target) < 0.5e-3) then
              !ek_high = (ek_high+ek_low)/2d0
              ek_low = ek_factor
            else 
              exit 
            endif
            ek_factor=(ek_high+ek_low)/2d0
          enddo
        else 
            call getcut(boxcut,ecutpaw,gmet,gsqcut,0,74,kpt,ngfft,ek_factor)
            gsqcutdg=gsqcut
            ekcut=ecutpaw
        endif
        !if(.not.ekread) ekcut=ecutpaw
        write(*,*) 'ekcut',ekcut
        write(*,*) 'ek_factor',ek_factor
    endif

        call bcast_scalar_r8(ecutpaw)
        call bcast_scalar_r8(ekcut)
        call bcast_scalar_r8(gsqcutdg)
        call bcast_scalar_r8(gsqcut)

        qmax    =1.2d0*sqrt(gsqcut)
        dq=qmax/(1.0*(mqgrid_ff-1))
        do ii=1,mqgrid_ff
            qgrid_ff(ii)=(ii-1)*dq
        end do
    
        qmax    =1.2d0*sqrt(gsqcutdg)
        dq=qmax/(1.0*(mqgrid_vl-1))
        do ii=1,mqgrid_vl
            qgrid_vl(ii)=(ii-1)*dq
        end do

        do ia=1,natom
            xred(1,ia)=atominfo(ia)%coord(1)/xmax/2d0+0.5d0
            xred(2,ia)=atominfo(ia)%coord(2)/ymax/2d0+0.5d0
            xred(3,ia)=atominfo(ia)%coord(3)/zmax/2d0+0.5d0
        enddo
    end subroutine set_param_alloc
    
end subroutine

