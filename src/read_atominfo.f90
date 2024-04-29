subroutine read_atominfo
    use atom_mod
    use main_mod, only : xmax,ymax,zmax
    implicit none

    integer :: i

    open(unit=19,file='cnt.ini',status='old',action='read')

    call count_atoms
    call alloc_atoms
    call  read_atoms

    close(19)
contains
    subroutine read_atoms
        implicit none

        character(len=2)  :: chz
        integer           :: z
        real*8            :: xx,yy,zz
       
        do i=1,natom
            read(19,*) chz,xx,yy,zz
            if(abs(xx)>xmax .or. abs(yy)>ymax .or. abs(zz)>zmax) then
                write(*,*) 'for now put atoms between -L/2 to L/2'
                write(*,*) 'ia,xx,yy,zz',i,xx,yy,zz
                stop
            endif
            call findtype(z,chz)
            atom_z(i)=z
            atominfo(i)%coord=(/xx,yy,zz/)
        enddo

        write(*,*) 'Number of atoms: ',natom
        write(*,*) 'And their atomic numbers: ',atom_z
    end subroutine read_atoms

    subroutine alloc_atoms
        implicit none

        allocate(atom_z(natom),atom_map(natom),atominfo(natom),&
            ngrid(natom),atindx1(natom),atindx(natom),&
            stat=stat)
!            atominfo_old(natom),stat=stat)
        if(stat/=0) stop 'problem allocate atoms'
    end subroutine alloc_atoms

    subroutine count_atoms
        implicit none
        !count number of atoms

        character(len=2)  :: chz
        i=0
        do while (.true.)
            read(19,*,end=11) chz
            i=i+1
        enddo
11      rewind(19)
            
        natom=i
    end subroutine count_atoms
end subroutine read_atominfo  
