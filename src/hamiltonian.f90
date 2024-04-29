subroutine get_h
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    implicit none
    integer :: ib,jb,ik
    real*8  :: hphi(nx,ny,nz)
    real*8  :: hm(nb,nb)

    write(*,*) 'hamiltonian'
    do ik=1,nk_loc
        do ib=1,nb
                call h_phi(phit(:,:,:,ib,ik),hphi,ik)
                do jb=1,nb
                    hm(ib,jb)=sum(hphi*phit(:,:,:,jb,ik))*dv
                enddo
            write(*,'(100(1x,f9.5))') hm(ib,:)
        enddo        
    enddo

!    call get_h_debug
end subroutine get_h      
