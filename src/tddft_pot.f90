subroutine get_pot_tddft(dens,dens_pert)
    use main_mod, only :  nx,ny,nz,vh,vks,vloc_tot,vxc,nn,&
        nhat,ncoret,scale_vh,nn,funct,periodic, h_type, &
        funct_x, funct_c
    use param_mod, only : k_flg
    use mpi_lib_ours
    use tddft_mod, only : vh_pert, vxc_pert, vks_pert, nhat_pert, prop_bar
    implicit none

    real*8  :: dens(nx,ny,nz)
    real*8  :: dens_pert(nx,ny,nz)
    real*8  :: vxc3d(nx,ny,nz)
    real*8  :: vxc3d_pert(nx,ny,nz)

    real*8, allocatable :: vks_tot(:)
    real*8, allocatable :: vks_tot_pert(:)


!    call debug_rho
    if(k_flg) then
      !k_flg is always false for now
      !if(rank==0) then
      !  call get_vxc(vxc3d)
      !  call vh_sub(dens+nhat,vh,scale_vh)
      !  if(periodic) then
      !      vh=vh-sum(vh)/size(vh)
      !  endif
      !  vks=vh+vxc3d+vloc_tot
      !  call diis_general(vks,nn)
      !endif
      !call bcast_r8(vks,size(vks),0)
      !call bcast_r8(vxc,size(vxc),0)
      !call bcast_r8(vh,size(vh),0)
    else
      if(rank==0) then
        call get_vxc_tddft(vxc3d, vxc3d_pert)    
        !if (prop_bar) then
        !if (h_type .eq. 1) then
        !  call vh_sub(dens,vh,scale_vh)
        !  call vh_sub(dens_pert,vh_pert,scale_vh)
        !else
          call vh_sub(dens+nhat,vh,scale_vh)
          call vh_sub(dens_pert+nhat_pert,vh_pert,scale_vh)
        !endif
        if(periodic) then
          vh=vh-sum(vh)/size(vh)
          vh_pert=vh_pert-sum(vh_pert)/size(vh_pert)
        endif
        vks=vh+vxc3d+vloc_tot
        vks_pert=vh_pert+vxc3d_pert+vloc_tot
        !write(6,*)  'get_pot vxc', sum(vxc**3d0), sum(vxc)
        !write(6,*)  'get_pot vh', sum(vh**3d0), sum(vh)
        !write(6,*)  'get_pot vloc_tot', sum(vloc_tot**3d0), sum(vloc_tot)
        !write(6,*)  'get_pot vks', sum(vks**3d0), sum(vks)
      endif
      call bcast_r8(vks,size(vks),0)
      call bcast_r8(vks_pert,size(vks),0)
    endif
contains

    subroutine get_vxc_tddft(v,v_p)
        implicit none

        real*8 :: vx(nn)
        real*8 :: vc(nn)
        real*8 :: v(nn)
        real*8 :: vx_p(nn)
        real*8 :: vc_p(nn)
        real*8 :: v_p(nn)

        if(funct==0) then
            !call get_vx_lda(dens+ncoret,nn,vx)
            !call get_vx_lda(dens_pert+ncoret,nn,vx_p)
            !call get_vcn_spn_lda(dens+ncoret,nn,1,vc)
            !call get_vcn_spn_lda(dens_pert+ncoret,nn,1,vc_p)
            !vxc(:,1)=vx+vc
            !vxc_pert(:,1)=vx_p+vc_p
            !v=vxc(:,1)
            !v_p=vxc_pert(:,1)

            !call vxc_lda_libxc(dens+ncoret,v,nn,1)
            !call vxc_lda_libxc(dens_pert+ncoret,v_p,nn,1)
            !vxc(:,1)=v
            !vxc_pert(:,1)=v_p

            funct_x = 1
            funct_c = 12
        else if(funct==1) then
            !call vxc_pbe_new(dens+ncoret,v,nn,1)
            !call vxc_pbe_new(dens_pert+ncoret,v_p,nn,1)
            !vxc(:,1)=v
            !vxc_pert(:,1)=v_p

            funct_x = 101
            funct_c = 130
        else
            write(*,*) 'funct should be 0 or 1'
            stop
        endif

        call vxc_libxc(dens+ncoret,v,nn,1)
        call vxc_libxc(dens_pert+ncoret,v_p,nn,1)
        vxc(:,1)=v
        vxc_pert(:,1)=v_p

!        vxc=vx
    end subroutine get_vxc_tddft
end subroutine get_pot_tddft

subroutine get_pot_tddft_2(pert)
    use main_mod, only :  nx,ny,nz,vh,vks,vloc_tot,vxc,nn,&
        nhat,ncoret,scale_vh,nn,funct,periodic, h_type, &
        funct_x, funct_c, dens
    use param_mod, only : k_flg
    use mpi_lib_ours
    use tddft_mod, only : vh_pert, vxc_pert, vks_pert, nhat_pert, prop_bar
    use tddft_mod, only : dens_pert
    implicit none

    real*8  :: dens3d(nx,ny,nz)
    real*8  :: dens3d_pert(nx,ny,nz)
    real*8  :: vxc3d(nx,ny,nz)
    real*8  :: vxc3d_pert(nx,ny,nz)
    logical :: pert

    real*8, allocatable :: vks_tot(:)
    real*8, allocatable :: vks_tot_pert(:)


    if(.not.pert) dens3d = reshape(dens,(/nx,ny,nz/))
    if(     pert) dens3d_pert = reshape(dens_pert,(/nx,ny,nz/))
!    call debug_rho
    if(k_flg) then
      !k_flg is always false for now
    else
      if(rank==0) then
        call get_vxc_tddft(vxc3d, vxc3d_pert)    
        if(.not. pert) call vh_sub(dens3d+nhat,vh,scale_vh)
        if(pert)       call vh_sub(dens3d_pert+nhat_pert,vh_pert,scale_vh)
        !endif
        if(periodic) then
          if(.not. pert) vh=vh-sum(vh)/size(vh)
          if(pert)       vh_pert=vh_pert-sum(vh_pert)/size(vh_pert)
        endif
        if(.not. pert) vks=vh+vxc3d+vloc_tot
        if(pert)       vks_pert=vh_pert+vxc3d_pert+vloc_tot
      endif
      if(.not. pert) call bcast_r8(vks,size(vks),0)
      if(pert)       call bcast_r8(vks_pert,size(vks),0)
    endif
contains

    subroutine get_vxc_tddft(v,v_p)
        implicit none

        real*8 :: vx(nn)
        real*8 :: vc(nn)
        real*8 :: v(nn)
        real*8 :: vx_p(nn)
        real*8 :: vc_p(nn)
        real*8 :: v_p(nn)

        if(funct==0) then
            funct_x = 1
            funct_c = 12
        else if(funct==1) then
            funct_x = 101
            funct_c = 130
        else
            write(*,*) 'funct should be 0 or 1'
            stop
        endif

        if(.not. pert) call vxc_libxc(dens3d+ncoret,v,nn,1)
        if(.not. pert) vxc(:,1)=v
        if(pert)       call vxc_libxc(dens3d_pert+ncoret,v_p,nn,1)
        if(pert)       vxc_pert(:,1)=v_p

!        vxc=vx
    end subroutine get_vxc_tddft
end subroutine get_pot_tddft_2
