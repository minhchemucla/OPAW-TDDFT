subroutine prep_phi0
      use param_mod, only : k_flg
!      use mat_module, only : mat_diag
      use main_mod
      use mpi_lib_ours
      implicit none

      integer :: ib,i,j
      real*8  :: tmp(nx,ny,nz)

      if(k_flg) then
        do ib=1,nb
          do i=1,nk_loc
            call rand_r(tmp,nn,dv)     
            phit(:,:,:,ib,i)=tmp
          enddo
        enddo
      else
        do ib=1,nb*nodes
          do i=1,nk_loc
            call rand_r(tmp,nn,dv)
            if(rank==0) then
              phit_tot(:,:,:,ib,i)=tmp
            endif
          enddo
        enddo
        call scatter_c16(phit_tot,phit,size(phit),0)
      endif
!      if(.not.k_flg) call bcast_c16(phit,size(phit),0)
end subroutine prep_phi0      
