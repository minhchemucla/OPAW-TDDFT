subroutine get_pot1(dens)
    use main_mod, only :  nx,ny,nz,vh,vks,vloc_tot,vxc,nn,&
        nhat,ncoret,scale_vh,nn,funct,periodic,funct_c,funct_x
    use param_mod, only : k_flg
    use mpi_lib_ours
    implicit none

    real*8  :: dens(nx,ny,nz)
    real*8  :: vxc3d(nx,ny,nz)

    real*8, allocatable :: vks_tot(:)


!    call debug_rho
    if(k_flg) then
      if(rank==0) then
        call get_vxc(vxc3d)
        call vh_sub(dens+nhat,vh,scale_vh)
        if(periodic) then
            vh=vh-sum(vh)/size(vh)
        endif
        vks=vh+vxc3d+vloc_tot
      endif
      call bcast_r8(vks,size(vks),0)
      call bcast_r8(vxc,size(vxc),0)
      call bcast_r8(vh,size(vh),0)
    else
      if(rank==0) then
        call get_vxc(vxc3d)    
        call vh_sub(dens+nhat,vh,scale_vh)
        if(periodic) then
          vh=vh-sum(vh)/size(vh)
        endif
        write(*,*) 'vh',sum(vh**3)
        write(*,*) 'vloc_tot',sum(vloc_tot**3)
        vks=vh+vxc3d+vloc_tot
      endif
      call bcast_r8(vks,size(vks),0)
    endif
!    open(unit=1,file='vtot')
!    write(1,*) vks
!    close(1)
!    open(unit=1,file='vxc')
!    write(1,*) vxc
!    close(1)
!    open(unit=1,file='vh')
!    write(1,*) vh
!    close(1)
contains
    subroutine debug_rho
        implicit none

        open(unit=1,file='nhat')
        read(1,*) nhat
        close(1)
        open(unit=1,file='ncore')
        read(1,*) ncoret
        close(1)
    end subroutine

    subroutine get_vxc(v)
        implicit none

        real*8 :: vx(nn)
        real*8 :: vc(nn)
        real*8 :: v(nn)

        if(funct==0) then
            !call get_vx_lda(dens+ncoret,nn,vx)
            !call get_vcn_spn_lda(dens+ncoret,nn,1,vc)
            !vxc(:,1)=vx+vc
            !v=vxc(:,1)

            !call vxc_lda_libxc(dens+ncoret,v,nn,1)
            !vxc(:,1)=v

            funct_x = 1
            funct_c = 12
        else if(funct==1) then
            !call vxc_pbe_new(dens+ncoret,v,nn,1)
            !vxc(:,1)=v

            funct_x = 101
            funct_c = 130
        else
            write(*,*) 'funct should be 0 or 1'
            stop
        endif
        call vxc_libxc(dens+ncoret,v,nn,1)
        vxc(:,1)=v

!        vxc=vx
    end subroutine get_vxc
end subroutine get_pot1
