subroutine bcast_s
    use mpi_lib_ours,only : bcast_r8
    use atom_mod,only : at=>atominfo,natom
    implicit none

    integer :: ia

    do ia=1,natom

        call bcast_r8(at(ia)%s,size(at(ia)%s),0)
        call bcast_r8(at(ia)%sinv,size(at(ia)%sinv),0)
        call bcast_r8(at(ia)%ssqinv,size(at(ia)%ssqinv),0)

    enddo
end subroutine
