subroutine set_rhoij
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use libpaw_mod
    implicit none

    integer :: ia
    integer :: ndim,ndim1,i,j,ind,nselect

    do ia=1,natom
        !pawrhoij only stores the bottom lower triangle of rhoij
        pawrhoij(ia)%rhoijp=0d0
        pawrhoij(ia)%rhoijselect=0d0
        pawrhoij(ia)%nrhoijsel=0
        ndim=size(at(ia)%rhoij,1)
        ndim1=size(pawrhoij(ia)%rhoijp,1)
        !sum of arithmetric series from 1,...,ndim
        ! = n/2(a_1+a_n) = ndim*(1+ndim)/2
        if (ndim*(ndim+1)/2/=ndim1) then
            write(*,*) 'dim of rhoij not match,ia,ndim,ndim1',ia,ndim,ndim1
            stop
        endif
        
        ind=0

        do j=1,ndim
            do i=1,j
                ind=ind+1
                pawrhoij(ia)%rhoijp(ind,1)=at(ia)%rhoij(i,j)
            enddo
        enddo

        nselect=0
        do i=1,ndim1
            if(abs(pawrhoij(ia)%rhoijp(i,1))>1d-10) then
                nselect=nselect+1
                pawrhoij(ia)%rhoijselect(nselect)=i
                pawrhoij(ia)%rhoijp(nselect,1)=pawrhoij(ia)%rhoijp(i,1)
            endif
        enddo

        pawrhoij(ia)%nrhoijsel=nselect
!        write(600,*) pawrhoij(ia)%rhoijp
    enddo

end subroutine

