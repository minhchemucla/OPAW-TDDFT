subroutine get_sij
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use mpi_lib_ours, only : rank
    implicit none

    integer :: ib,jb,ig,ia,ix,iy,iz,is,js,it,ik
    complex*16  :: tmp(nx,ny,nz)
    complex*16  :: sij(nb,nb)
!    real*8  :: fac

!    fac=sqrt(dx*dy*dz)*sqrt(dble(nx*ny*nz))

    write(15000+rank+1,*) 'sij'
    do ik=1,nk_loc
        write(15000+rank+1,*) 'ik',ik
        do ib=1,nb
            call sn_phi(phit(:,:,:,ib,ik),tmp,ik,1d0)
            do jb=1,nb
                sij(ib,jb)=sum(conjg(phit(:,:,:,jb,ik))&
                    *tmp)*dv
            enddo
            write(15000+rank+1,'(100(1x,f22.15))') dble(sij(ib,:))
        enddo
        write(15000+rank+1,*) 'imag',maxval(dimag(sij))
    enddo
end subroutine      
