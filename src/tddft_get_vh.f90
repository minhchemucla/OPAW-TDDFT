subroutine get_vh_tddft(dens,dens_pert)
    use main_mod, only :  nx,ny,nz,vh,vks,vloc_tot,vxc,nn,&
        nhat,ncoret,scale_vh,nn,funct,periodic
    use param_mod, only : k_flg
    use mpi_lib_ours
    use tddft_mod, only : vh_pert, vxc_pert, vks_pert, nhat_pert, prop_bar
    implicit none

    real*8  :: dens(nx,ny,nz)
    real*8  :: dens_pert(nx,ny,nz)
    real*8  :: vxc3d(nx,ny,nz)
    real*8  :: vxc3d_pert(nx,ny,nz)


    if(rank==0) then
      call vh_sub(dens+nhat,vh,scale_vh)
      call vh_sub(dens_pert+nhat_pert,vh_pert,scale_vh)
      if(periodic) then
        vh=vh-sum(vh)/size(vh)
        vh_pert=vh_pert-sum(vh_pert)/size(vh_pert)
      endif
      vxc3d = reshape(vxc, (/nx,ny,nz/))
      vxc3d_pert = reshape(vxc_pert, (/nx,ny,nz/))
      vks=vh+vxc3d+vloc_tot
      vks_pert=vh_pert+vxc3d_pert+vloc_tot
!        write(6,*)  'get_vh vxc', sum(vxc**3d0), sum(vxc)
!        write(6,*)  'get_vh vh', sum(vh**3d0), sum(vh)
!        write(6,*)  'get_vh vloc_tot', sum(vloc_tot**3d0), sum(vloc_tot)
!        write(6,*)  'get_vh vks', sum(vks**3d0), sum(vks)
    endif
    call bcast_r8(vks,size(vks),0)
    call bcast_r8(vks_pert,size(vks),0)
end subroutine get_vh_tddft
