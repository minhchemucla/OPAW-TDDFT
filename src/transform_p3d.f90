subroutine transform_p3d(ia)
    ! orthogonalize projectors on rough grid
    use mat_module, only : mat_diag
    use main_mod
    use atom_mod
    use tddft_mod, only : tddft_flag
    use atom_mod, only : p => pawinfo, at => atominfo
    use mpi_lib_ours
    implicit none

    integer :: ia,it,ms
    real*8, allocatable :: mat(:,:),vec(:,:),eig(:)
    real*8, allocatable :: mat1(:,:), mat2(:,:), mat3(:,:)
    real*8, allocatable :: p3d(:,:,:,:)
    integer :: nrough(3)

    it=atom_map(ia)
    nrough=p(it)%nrough

    if(rank==0) then
        ms=p(it)%mstates
        allocate(mat(ms,ms),vec(ms,ms),eig(ms),&
            mat1(ms,ms),mat2(ms,ms),mat3(ms,ms),stat=stat)
        if(stat/=0) stop 'alloc mat,vec,eig'
        allocate(p3d(nrough(1),nrough(2),nrough(3),ms),stat=stat)
        if(stat/=0) stop 'alloc p3d for transform'

        !if(tddft_flag>-1) then
        !  allocate(at(ia)%mat_transform_c(ms,ms), at(ia)%vec_transform_c(ms,ms))
        !  if(stat/=0) stop 'alloc tddft mat,vec for transform'
        !endif
    
        call check_mat1 !equation A1. L matrix
        call mat_diag(mat,vec,eig) 
        call transform1 !equation A2 and A3 and equation after A3
        call mat_diag(mat,vec,eig)
        call transform2  !equation after A4
        
        deallocate(mat,vec,eig,p3d,mat1,mat2,mat3)
    endif
contains
    subroutine transform2
        implicit none
        integer :: is,js
        !if(tddft_flag>-1) at(ia)%vec_transform_c = vec
        p3d=at(ia)%local_p3d1
        at(ia)%local_p3d1=0d0
        do is=1,ms
            do js=1,ms
                at(ia)%local_p3d1(:,:,:,is)=at(ia)%local_p3d1(:,:,:,is)+&
                    p3d(:,:,:,js)*vec(js,is)
            enddo
            at(ia)%s(is,:)=eig(is)
            at(ia)%sinv(is,:)=1d0/(1d0+eig(is))-1d0
            at(ia)%ssqinv(is,:)=1d0/sqrt(1d0+eig(is))-1d0
        enddo


    end subroutine transform2

    subroutine transform1
        implicit none
        integer :: i,j,k
        p3d=at(ia)%local_p3d1
        at(ia)%local_p3d1=0d0
        mat3=0d0

        do j=1,ms
            mat3=p(it)%sij
            do k=1,ms
                mat1(j,k)=sum(vec(j,:)*vec(k,:)/sign(sqrt(abs(eig)),eig))
                mat2(j,k)=sum(vec(j,:)*vec(k,:)*sign(sqrt(abs(eig)),eig))
            enddo
        enddo

        do i=1,ms
            do j=1,ms
                at(ia)%local_p3d1(:,:,:,i)=at(ia)%local_p3d1(:,:,:,i)+&
                    p3d(:,:,:,j)*mat1(i,j)
            enddo
        enddo

        !if(tddft_flag>-1) at(ia)%mat_transform_c = mat1

        mat=matmul(matmul(transpose(mat2),mat3),mat2)
    end subroutine

    subroutine check_mat1
        implicit none
        integer :: is,js
        do is=1,ms
            do js=1,ms
                mat(is,js)=sum(at(ia)%local_p3d1(:,:,:,is)*&
                    at(ia)%local_p3d1(:,:,:,js))*dv
            enddo
!            write(*,'(100f9.5)') mat(is,:)
        enddo
    end subroutine    
end subroutine      
