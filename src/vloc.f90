subroutine get_vloc_ncoret
    use main_mod
    use libpaw_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use m_atm2fft
    use mpi_lib_ours
    implicit none

    real*8, allocatable :: ph1d(:,:)
    integer :: mgrid

    allocate(vloc_tot(nx,ny,nz),ncoret(nx,ny,nz),stat=stat)
    if (stat/=0) stop 'vloc_tot alloc problem'

    if(rank==0) then
        call get_mgrid
        allocate(ph1d(2,3*(2*mgrid+1)*natom),stat=stat)
        if(stat/=0) stop 'ph1d alloc problem'

        ph1d=0d0
        call getph(atindx,natom,nx,ny,nz,ph1d,xred,size(ph1d,1),size(ph1d,2))

        !Get the PS valence density (optn2=2)
        call atm2fft(atindx1,dens,vloc_tot,gmet,gprimd,gsqcut,mgrid,mqgrid_vl,&
            natom,nattyp,nn,ngfft,ntypat,pawtab,ph1d,&
            qgrid_vl,ucvol,vlspl,2)
        !Get the PS core density (optn2=1)
        call atm2fft(atindx1,ncoret,vloc_tot,gmet,gprimd,gsqcut,mgrid,mqgrid_vl,&
            natom,nattyp,nn,ngfft,ntypat,pawtab,ph1d,&
            qgrid_vl,ucvol,vlspl,1)

        !Minh
        !if(.not. periodic) then
        !    call vloc_tuma_prep
        !endif

        deallocate(ph1d)

    endif

    !Minh
    if(.not. periodic) then
        call vloc_tuma_prep
    endif
    call bcast_r8(ncoret,size(ncoret),0)
    call bcast_r8(vloc_tot,size(vloc_tot),0)
    call bcast_r8(dens,size(ncoret),0)
    
    if(rank==0) then
        write(*,*) 'sum dens0+nhat0',sum(dens)*dv
        write(*,*) 'finish preparing vloc,ncoret'
        write(*,*) '========================'
        write(*,*)
    endif

!write(977,*) ph1d

!    open(unit=1,file='vloc_tot')
!    write(1,*) vloc_tot
!    close(1)
!    open(unit=1,file='myncoret')
!    write(1,*) ncoret
!    close(1)
contains
    subroutine get_mgrid
        implicit none

        mgrid=nx
        if (ny>mgrid) then
            mgrid=ny
        endif
        if (nz>mgrid) then
            mgrid=nz
        endif

    end subroutine
end subroutine      
