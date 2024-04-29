subroutine get_nhat_pert
    use main_mod
    use tddft_mod, only : nhat_pert, dens_pert
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use mpi_lib_ours, only : rank
    implicit none

    real*8 :: a,b
    integer :: ia,ig,it,is,js,igl,ix,iy,iz

    nhat_pert=0d0

    do ia=1,natom
        it=atom_map(ia)
        do ig=1,ngrid(ia)
            ix=at(ia)%local_grid(ig,1)
            iy=at(ia)%local_grid(ig,2)
            iz=at(ia)%local_grid(ig,3)
            do is=1,p(it)%mstates
                do js=1,p(it)%mstates
                    do igl=1,(2*p(it)%nl-1)**2
                        nhat_pert(ix,iy,iz)=nhat_pert(ix,iy,iz)+at(ia)%rhoij_pert(is,js) &
                            *p(it)%qijlm(igl,is,js)*at(ia)%local_g3d(ig,igl)
                    enddo
                enddo
            enddo
        enddo
    enddo

    !if(rank==0)write(*,*) 'nhat_pert',sum(nhat_pert**3)
    a=sum(nhat_pert)*dv
    b=sum(dens_pert)*dv
    !if(rank==0)write(*,*) 'sum dens_pert',b
    !if(rank==0)write(*,*) 'sum dens_pert+nhat_pert',a+b
    nhat_pert=nhat_pert*(rnel)/(a+b)
    dens_pert=dens_pert*(rnel)/(a+b)
    a=sum(nhat_pert)*dv
    b=sum(dens_pert)*dv
    !if(rank==0)write(*,*) 'adjusted sum dens_pert+nhat_pert',a+b

end subroutine get_nhat_pert
