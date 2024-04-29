!subroutine transform_c(ia)
!    ! orthogonalize projectors on rough grid
!    use main_mod
!    use atom_mod
!    use atom_mod, only : p => pawinfo, at => atominfo
!    use mpi_lib_ours
!    implicit none
!
!    integer :: ia,it,ms
!    real*8, allocatable :: c3d(:,:,:,:)
!    integer :: nrough(3)
!
!    it=atom_map(ia)
!    nrough=p(it)%nrough
!
!    if(rank==0) then
!        ms=p(it)%mstates
!        allocate(c3d(nrough(1),nrough(2),nrough(3),ms),stat=stat)
!        if(stat/=0) stop 'alloc c3d for transform'
!    
!        call transform1 !equation A2 and A3 and equation after A3
!        call transform2  !equation after A4
!    endif
!contains
!    subroutine transform2
!        implicit none
!        integer :: is,js
!        c3d=at(ia)%local_c
!        at(ia)%local_c=0d0
!        do is=1,ms
!            do js=1,ms
!                at(ia)%local_c(:,:,:,is)=at(ia)%local_c(:,:,:,is)+&
!                    c3d(:,:,:,js)*at(ia)%vec_transform_c(js,is)
!            enddo
!        enddo
!
!    end subroutine transform2
!
!    subroutine transform1
!        implicit none
!        integer :: i,j,k
!        c3d=at(ia)%local_c
!        at(ia)%local_c=0d0
!
!        do i=1,ms
!            do j=1,ms
!                at(ia)%local_c(:,:,:,i)=at(ia)%local_c(:,:,:,i)+&
!                    c3d(:,:,:,j)*at(ia)%mat_transform_c(i,j)
!            enddo
!        enddo
!    end subroutine
!end subroutine      
