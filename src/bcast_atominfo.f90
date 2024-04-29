subroutine bcast_atominfo
    use main_mod, only : nk_loc
    use mpi_lib_ours
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo

    implicit none

    integer :: nrough(3),ia,it,ms,is
    integer :: qmmax
    call bcast_scalar_i(natom)
    call bcast_scalar_i(ntype)
    if(rank/=0) then
        allocate(atom_map(natom),at(natom),p(ntype),atom_z(natom),stat=stat)
        if(stat/=0) stop 'atom alloc nonzero rank'
        allocate(ngrid(natom),atindx1(natom),atindx(natom),diff_z(ntype),stat=stat)
        if(stat/=0) stop 'atom alloc nonzero rank 1'
        allocate(natom_type(ntype),stat=stat)
        if(stat/=0) stop 'atom alloc nonzero rank 2'
    endif

    call bcast_i (atom_map,natom,0)
    call bcast_i (atom_z,natom,0)
    call bcast_i (atindx1,natom,0)
    call bcast_i (atindx,natom,0)
    call bcast_i (diff_z,ntype,0)
    call bcast_i (natom_type,ntype,0)
    call bcast_i (ngrid,natom,0)

    do it=1,ntype
        call bcast_scalar_i(p(it)%nl)
        call bcast_scalar_i(p(it)%mstates)
        call bcast_i(p(it)%nrough,3,0)
        call bcast_i(p(it)%nfine,3,0)
        call bcast_r8(p(it)%dxf,3,0)
        call bcast_char264(p(it)%filename,0)

        if(rank==0) qmmax=size(p(it)%qijlm,1)
        call bcast_scalar_i(qmmax)
        ms=p(it)%mstates
        if(rank/=0) allocate(p(it)%qijlm(qmmax,ms,ms),stat=stat)
        if(stat/=0) stop 'qmmax alloc problem rank/=0'
        call bcast_r8(p(it)%qijlm,qmmax*ms*ms,0)
    enddo

    do ia=1,natom
        call bcast_r8(at(ia)%coord,   3,0)
        call bcast_i (at(ia)%ir_start,3,0)
        call bcast_scalar_l(at(ia)%edge)
        it=atom_map(ia)
        ms=p(it)%mstates
        if(rank/=0) then
            allocate(at(ia)%rhoij(ms,ms),stat=stat)
            if(stat/=0) stop 'problem allocated rhoij, rank/=0' 
            allocate(at(ia)%rhoij_pert(ms,ms),stat=stat)
            if(stat/=0) stop 'problem allocated rhoij_pert, rank/=0' 
            allocate(at(ia)%rho0ij(ms,ms),stat=stat)
            if(stat/=0) stop 'problem allocated rho0ij, rank/=0' 
            allocate(at(ia)%dij(ms,ms),stat=stat)
            if(stat/=0) stop 'problem allocated dij, rank/=0' 
            allocate(at(ia)%dij_pert(ms,ms),stat=stat)
            if(stat/=0) stop 'problem allocated dij_pert, rank/=0' 
            allocate(at(ia)%dij_old(ms,ms),stat=stat)
            if(stat/=0) stop 'problem allocated dij_old, rank/=0' 
            allocate(at(ia)%dij_pert_old(ms,ms),stat=stat)
            if(stat/=0) stop 'problem allocated dij_pert_old, rank/=0' 
            at(ia)%rhoij=0d0
            at(ia)%dij=0d0
            at(ia)%dij_old=0d0
            at(ia)%rhoij_pert=0d0
            at(ia)%dij_pert=0d0
            at(ia)%dij_pert_old=0d0
            allocate(at(ia)%s(ms,nk_loc),at(ia)%sinv(ms,nk_loc)&
                ,at(ia)%ssqinv(ms,nk_loc),stat=stat)
            if(stat/=0) stop 'alloc s,sinv,ssqinv'
        endif

        if(rank==0) then
            do is=1,ms
                at(ia)%rhoij(is,is)=p(it)%occ_l(p(it)%ms_ls(is))
            enddo
            write(*,*) 'ia,rhoij0',ia,(at(ia)%rhoij(is,is), is=1,ms)
        endif

        call bcast_r8(at(ia)%rhoij,ms*ms,0)
        call bcast_r8(at(ia)%rhoij_pert,ms*ms,0)
        call bcast_r8(at(ia)%rho0ij,ms*ms,0)
        
        call bcast_r8(at(ia)%s,size(at(ia)%s),0)
        call bcast_r8(at(ia)%sinv,size(at(ia)%sinv),0)
        call bcast_r8(at(ia)%ssqinv,size(at(ia)%ssqinv),0)

    enddo
    if(rank==0) then
        write(*,*) 'finished bcast atompaw'
    endif
end subroutine      
