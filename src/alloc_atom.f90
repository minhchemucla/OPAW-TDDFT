subroutine alloc_atompaw
    use main_mod, only : nk_loc
    use paw_mod
    use atom_mod, only : at => atominfo, p => pawinfo,&
            natom, atom_map!, at_old => atominfo_old
    implicit none

    integer :: ms,nl,ia,it,nmp,is,ik

    do ia=1,natom
        it=atom_map(ia)
        ms=p(it)%mstates
        nl=p(it)%nl

        allocate(at(ia)%rhoij(ms,ms))
        allocate(at(ia)%rhoij_pert(ms,ms))
        allocate(at(ia)%rho0ij(ms,ms))
        if(stat/=0) stop 'problem allocated rhoij' 
!        allocate(at_old(ia)%rhoij(ms,ms))
!        if(stat/=0) stop 'problem allocated rhoij old' 
        allocate(at(ia)%dij(ms,ms),at(ia)%dij_old(ms,ms))
        allocate(at(ia)%dij_pert(ms,ms),at(ia)%dij_pert_old(ms,ms))
        if(stat/=0) stop 'problem allocated dij' 
        at(ia)%rhoij=0d0
        at(ia)%rhoij_pert=0d0
        at(ia)%rho0ij=0d0
        at(ia)%dij=0d0
        at(ia)%dij_pert=0d0
!        at_old(ia)%rhoij=0d0
!        do is=1,ms
!            at(ia)%rhoij(is,is)=1d0
!            at_old(ia)%rhoij(is,is)=1d0
!        enddo
        allocate(at(ia)%s(ms,nk_loc),at(ia)%sinv(ms,nk_loc)&
            ,at(ia)%ssqinv(ms,nk_loc),stat=stat)
        if(stat/=0) stop 'alloc s,sinv,ssqinv'
    enddo 

    write(*,*) 'finished allocating atompaw for rank0'
end subroutine alloc_atompaw

