subroutine acc_coeff1
    use main_mod, only : nocc, phit_tot,wk, nk_loc, nb, iscf
    use paw_mod
    use atom_mod, only : at => atominfo, p => pawinfo, atom_map, ngrid, natom
    use mpi_lib_ours
    implicit none

    integer :: ia,ig,ix,iy,iz,is,ms,it,js,ib,ik,jk,iib
    complex*16, allocatable :: ca(:),tmp(:,:)

    do ia=1,natom
        at(ia)%rhoij=0d0
    enddo

    do ik=1,nk_loc
        jk=(ik-1)*nodes+rank+1
        do ib=1,nocc
!          iib=rank*nb+ib
!          write(9000+rank,*) 'ib,rank,iib',ib,rank,iib
!          if(iib.le.nocc) then
            do ia=1,natom
!                write(116,*) 'ca for atom:',ia,'band',ib
                it=atom_map(ia)
                ms=p(it)%mstates
                allocate(ca(ms),tmp(ms,ms),stat=stat)
                tmp=0d0
                if(stat/=0) stop 'ca alloc problem in coeff'

                call proj(ia,phit_tot(:,:,:,ib,ik),ca,ms,ik) 
                do is=1,ms
                    do js=1,ms
                        tmp(is,js)=tmp(is,js)+conjg(ca(is))*ca(js)*2d0 !Note
                    enddo
                enddo
                at(ia)%rhoij=at(ia)%rhoij+dble(tmp)*wk(jk)
                deallocate(ca,tmp)
            enddo
!          endif
        enddo
    enddo

    !write(117,*) phit_tot(1:3,1:3,1:3,1,1)
    !do ia=1,natom
    !    write(117,*) at(ia)%rhoij
    !enddo
    !call debug_rhoij
contains
    subroutine debug_rhoij
        implicit none
        open(unit=1,file='rhoij')
        do ia=1,natom
            read(1,*) at(ia)%rhoij 
        enddo
        close(1)
    end subroutine debug_rhoij
end subroutine acc_coeff1

