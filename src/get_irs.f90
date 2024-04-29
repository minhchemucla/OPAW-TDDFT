! Note : see
! https://www.abinit.org/sites/default/files/infos/7.8/tutorial/lesson_paw1.html
! for overlap of paw spheres
! getting starting indices and implicitly the ending indices for the rough grid
subroutine get_irs 
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use mpi_lib_ours 
    
    implicit none

    integer :: ia,it

    integer :: irs(3),ire(3)

    write(*,*) 'rough grid, start and end'
    do ia=1,natom
        write(*,*) 'atom number:',ia
        it=atom_map(ia)
        at(ia)%ir_start(1)=ceiling((at(ia)%coord(1)+xmax)/dx)-p(it)%nrough(1)/2+1

        at(ia)%ir_start(2)=ceiling((at(ia)%coord(2)+ymax)/dy)-p(it)%nrough(2)/2+1

        at(ia)%ir_start(3)=ceiling((at(ia)%coord(3)+zmax)/dz)-p(it)%nrough(3)/2+1
        irs=at(ia)%ir_start
        ire=irs+p(it)%nrough-1
        write(*,*) irs,ire
        if(minval(irs)<1 .or. ire(1)>nx .or. &
            ire(2)>ny .or. ire(3)>nz) then
            if(.not.periodic) then
                write(*,*) 'non periodic, local grid cannot got outside box'
                write(*,*) 'ia,irs,ire',ia,irs,ire
                stop
            else
                at(ia)%edge=.true.
            endif
        endif
    enddo
end subroutine get_irs 
