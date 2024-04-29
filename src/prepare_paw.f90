subroutine prep_paw
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use mpi_lib_ours, only :rank, sync_mpi, bcast_r8

    integer :: i
    logical :: phase_rough = .true.

    !read from cnt.ini and store information
    if (rank==0) then
        call read_atominfo
        call  map_atominfo
    
        do i=1,ntype
            !read paw files
            call readpaw(i)
!            call check_mat(i)
        enddo

    !find local grid around each atom
        call get_irs
        call map_localgrid
        call alloc_atompaw 
    endif

    call bcast_atominfo

!    call writef(' post atominfo ')

    if (phase_rough) then !phase_rough is always true because using rough grid
!        call writef(' in phaserough ')
        call get_local_p3d !projectors
        call get_localg3d !shape functions for compensation charge
        !call get_local_pp3d !pseudo partial waves
        !call get_local_aep3d !ae partial waves
    else
         stop 'use phase_rough'
    endif

    call bcast_s
    if(rank==0) then
        write(6,*) 'finished preparing paw'
        write(6,*) '======================'
        write(6,*)
        call flush(6)
    endif
    call sync_mpi
end subroutine prep_paw      
