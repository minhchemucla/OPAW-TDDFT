subroutine ncpaw_libpaw_prepare
    use main_mod, only : nx,ny,nz,nn,nhat,dens,vks,vxc,vh,ek
    implicit none

    call prep_kpt
    call prep_paw
    call prep_libpaw
    allocate(dens(nn))
    allocate(nhat(nx,ny,nz),vks(nx,ny,nz),vxc(nn,1),vh(nx,ny,nz))
    !nhat = compensation charge
    call prep_vks
    call get_vloc_ncoret
contains
    subroutine prep_vks
        implicit none

        call vk_prep
        call ek_prep
    end subroutine prep_vks    
end subroutine ncpaw_libpaw_prepare
