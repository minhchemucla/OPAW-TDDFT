subroutine map_atominfo
    use atom_mod
    use paw_mod
    use mpi_lib_ours
    implicit none

    integer :: iatom
    
    call count_types
    call   map_types
    call  read_filenames

    if(rank==0) write(*,*) 'finish reading atom info'
    if(rank==0) write(*,*) '========================'
    if(rank==0) write(*,*)
contains
    subroutine read_filenames
        implicit none

        integer :: i,z

        logical :: found(ntype)

        character*2  :: ch
        character*20 :: fn !filename

        open(unit=17,file='pawfiles',action='read')

        found=.false.

        do while (.true.)
        read(17,*,end=11) ch, fn
            call findtype(z,ch)
            do i=1,ntype
                if(z==diff_z(i)) then
                    if(.not.found(i)) then
                        pawinfo(i)%filename=trim(fn)
                        found(i)=.true.
                    else
                        write(*,*) 'error: paw filename for atom of type: ',z, ' already read'
                        stop
                    endif
                endif
            enddo
        enddo

11      do i=1,ntype
            if(.not.found(i)) then
                write(*,*) 'error: paw filename not given for atom of type: ',diff_z(i)
                stop
            endif
            if(rank==0) write(*,*) 'paw filename for atom of type: ',diff_z(i), ' is:'
            if(rank==0) write(*,*) pawinfo(i)%filename
        enddo
            
        close(17)

    end subroutine read_filenames

    subroutine map_types
        implicit none

        integer :: i,j,indx
        logical :: found

        natom_type=0

        do i=1,natom
            found=.false.
            do j=1,ntype
                if(atom_z(i)==diff_z(j)) then
                    atom_map(i)=j
                    found=.true.
                    natom_type(j)=natom_type(j)+1
                endif
            enddo
            if(.not.found) then
                write(*,*) 'error: type of atom #.', i, ' is not mapped'
                write(*,*) 'its atomic number is: ', atom_z(i)
                stop
            endif
        enddo
       
        indx=0
        do i=1,ntype
            do j=1,natom
                if(atom_z(j)==diff_z(i))then
                    indx=indx+1
                    atindx(j)=indx
                    atindx1(indx)=j
                endif
            enddo
        enddo

        if(rank==0) write(*,*) 'atoms types mapped: ', atom_map
        if(rank==0) write(*,*) 'number of atoms for each type: ',natom_type
        if(rank==0) write(*,*) 'atindx1',atindx1
    end subroutine map_types

    subroutine count_types
        implicit none

        integer :: i=0, min_val, max_val
        integer, allocatable :: unique(:)

        allocate(unique(natom),stat=stat)
        if(stat/=0) stop 'unique for counting types of atoms'

        min_val = minval(atom_z)-1
        max_val = maxval(atom_z)

        do while (min_val<max_val)
            i = i+1
            min_val = minval(atom_z, mask=atom_z>min_val)
            unique(i) = min_val
        enddo

        ntype=i

        allocate(diff_z(ntype),stat=stat,source=unique(1:ntype))
        if(stat/=0) stop 'diff_z alloc problem'
        allocate(pawinfo(ntype),stat=stat)
        if(stat/=0) stop 'pawinfo alloc problem'
        allocate(natom_type(ntype),stat=stat)
        if(stat/=0) stop 'natom_type alloc prblem'
        if(rank==0) write(*,*) 'number of different elements: ', ntype
        if(rank==0) write(*,*) 'z of the different elements: ', diff_z

        deallocate(unique)

    end subroutine count_types
end subroutine map_atominfo
