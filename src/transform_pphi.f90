subroutine transform_pphi(it)
! This subroutine does the orthogonalization of the projectors and etc..
! that have the same l values on the radial grid .
    use mat_module, only : mat_diag
    use main_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use mpi_lib_ours, only : rank
    implicit none

    integer :: it,ms
    integer :: i,j,counter,nr
    integer, allocatable :: index_state(:)
    real*8,  allocatable :: mat(:,:),vec(:,:),eig(:)
    real*8,  allocatable :: mat1(:,:),mat2(:,:),mat3(:,:)

    nr=p(it)%nr_int

    if(rank==0) then
        write(*,*) 'transform_p'
        write(*,*) 's before:'
        do i=1,p(it)%nstates
            write(*,*) p(it)%mat_t(:,i)
        enddo
        write(*,*) 'overlap before:'
        do i=1,p(it)%nstates
            write(*,*) p(it)%mat_sp(:,i)
        enddo
        
        do i=1,nr
            write(106,*) p(it)%rr(i),p(it)%ptilde(:,i)
            write(107,*) p(it)%rr(i),p(it)%phitilde(:,i)
            write(108,*) p(it)%rr(i),p(it)%phi(:,i)
        enddo
        do i=1,p(it)%nl
            call get_index
            call alloc_submats
    
            call get_submat
            call mat_diag(mat,vec,eig)
!            call checkdiasym  
            call transform_pphi1
            call mat_diag(mat,vec,eig)
            call transform_pphi2
    
!            if(flg==1) call checkss12
            call dealloc_submats
        enddo

        do i=1,nr
            write(103,*) p(it)%rr(i),p(it)%ptilde1(:,i)
            write(104,*) p(it)%rr(i),p(it)%phitilde1(:,i)
            write(105,*) p(it)%rr(i),p(it)%phi1(:,i)
        enddo

        write(*,*) 'p(it)%s after:'
        do i=1,p(it)%nstates
            do j=1,p(it)%nstates
               p(it)%mat_t (i,j)     = &
                     sum(p(it)%rr(1:nr)*p(it)%rr(1:nr)*&
                     p(it)%dr(1:nr)*p(it)%phi1(i,1:nr)&
                     *p(it)%phi1(j,1:nr)) -&
                     sum(p(it)%rr(1:nr)*p(it)%rr(1:nr)*&
                     p(it)%dr(1:nr)*p(it)%phitilde1(i,1:nr)&
                     *p(it)%phitilde1(j,1:nr))
            enddo
            write(*,*) p(it)%mat_t(i,:)
        enddo

        write(*,*) 'overlap of projectors'
        do i=1,p(it)%nstates
            do j=1,p(it)%nstates
                p(it)%mat_sp(i,j)     =sum(p(it)%rr(1:nr)*&
                    p(it)%rr(1:nr)*p(it)%dr(1:nr)&  
                    *p(it)%ptilde1(i,1:nr)*p(it)%ptilde1(j,1:nr))
            enddo
            write(*,*) p(it)%mat_sp(i,:)
        enddo
    endif

contains

    subroutine transform_pphi2
        implicit none
        
        p(it)%ptilde1  (index_state,:)=matmul(transpose(vec), p(it)%ptilde1  (index_state,:)) !Eq after Eq. (A.4)
        p(it)%phi1     (index_state,:)=matmul(transpose(vec), p(it)%phi1     (index_state,:))
        p(it)%phitilde1(index_state,:)=matmul(transpose(vec), p(it)%phitilde1(index_state,:))
        p(it)%mat_t(index_state,index_state)=matmul(matmul(transpose(vec),&  !o_i^(a) in Eq. (A.4)
            mat),vec)
    end subroutine transform_pphi2

    subroutine transform_pphi1
        implicit none
        
        integer :: j,k

        mat3=p(it)%mat_t(index_state,index_state)

        do j=1,counter
            do k=1,counter
                mat1(j,k)=sum(vec(j,:)*vec(k,:)/sign(sqrt(abs(eig)),eig))
                mat2(j,k)=sum(vec(j,:)*vec(k,:)*sign(sqrt(abs(eig)),eig))
            enddo
        enddo
        
        p(it)%ptilde1  (index_state,:)=matmul(mat1,           p(it)%ptilde  (index_state,:)) !Eq. (A.2)
        p(it)%phi1     (index_state,:)=matmul(transpose(mat2),p(it)%phi     (index_state,:))
        p(it)%phitilde1(index_state,:)=matmul(transpose(mat2),p(it)%phitilde(index_state,:))

        !write(*,*) 'check pphi1'
        !do j=1,counter
        !    do k=1,counter
        !        mat(j,k)     =sum(p(it)%rr(1:nr)*&
        !            p(it)%rr(1:nr)*p(it)%dr(1:nr)&  
        !            *p(it)%ptilde1(index_state(j),1:nr)*p(it)%ptilde1(index_state(k),1:nr))
        !    enddo
        !    write(*,*) mat(j,:)
        !enddo


        mat=matmul(matmul(transpose(mat2),mat3),mat2) !O^(a) under Eq. (A.3)
        write(*,*) 'the second sub matrix for l= ',p(it)%lstate_diff(i),' is:'
        do j=1,counter
            write(*,*) mat(:,j)
        enddo
    end subroutine transform_pphi1

    subroutine get_index
        implicit none
        integer :: j,k

        counter=count(p(it)%lstate.eq.p(it)%lstate_diff(i))
        write(*,*) 'number of states for l=',p(it)%lstate_diff(i),' is: ',counter
        if(allocated(index_state)) deallocate(index_state)
        allocate(index_state(counter),stat=stat)
        if(stat/=0) stop 'index_state alloc problem'

        k=0
        do j=1,p(it)%nstates
            if(p(it)%lstate(j).eq.p(it)%lstate_diff(i)) then
                k=k+1
                index_state(k)=j
            endif
        enddo
        write(*,*) 'the states with l=',p(it)%lstate_diff(i),' are: ',index_state

    end subroutine get_index

    subroutine dealloc_submats
        implicit none

        deallocate(mat,mat1,mat2,mat3,vec,eig,stat=stat)
        if(stat/=0) stop 'submats dealloc problem'
    end subroutine dealloc_submats

    subroutine alloc_submats
        implicit none

        allocate(mat(counter,counter),vec(counter,counter),&
            mat1(counter,counter),mat2(counter,counter),&
            mat3(counter,counter),eig(counter),stat=stat)
        if(stat/=0) stop 'submats alloc problem'
    end subroutine alloc_submats

    subroutine get_submat
        implicit none
        integer :: j,k

        write(*,*) 'the sub matrix for l= ',p(it)%lstate_diff(i),' is:'

        do k=1,counter
            do j=1,counter
                mat(k,j)=p(it)%mat_sp(index_state(k),index_state(j))
            enddo
            write(*,*) mat(k,:)
        enddo
        
    end subroutine get_submat
    
!    subroutine checkdiasym
!        implicit none
!        integer :: j
!
!           10 format(I3,'   ',f14.8)
!           20 format(i3,'   ',10f14.8)
!           30 format(10f14.8)
!
!        open(10,file='dia.dat',status='unknown',access='append')
!        
!        write(10,*)'Original matrix:'
!        do j=1,counter
!           write(10,20)j,ss0(:,j)
!        enddo
!
!        write(10,*)
!
!        write(10,*)'Eigenvalues:'
!        do j=1,counter
!           write(10,10)j,eig(j)
!        enddo
!
!        write(10,*)
!
!        write(10,*)'Eigenvectors:'
!        do j=1,counter
!           write(10,20)j,ss(:,j)
!        enddo
!        write(10,*)
!        
!        ss1=matmul(transpose(ss),ss0)
!        ss2=matmul(ss1,ss)
!
!        write(10,*)'Transformed matrix (check):'
!        do j=1,counter
!           write(10,30)ss2(:,j)
!        enddo
!        write(10,*)
!
!        ss1=matmul(transpose(ss),ss)
!
!        write(10,*)'Check if eigenvector matrix is unitary:'
!        do j=1,counter
!           write(10,30)ss1(:,j)
!        enddo
!        write(10,*)
!
!        close(10)
!    end subroutine checkdiasym
!
!    subroutine checkss12
!        implicit none
!        integer :: j
!        real*8  :: tmp(counter,counter)
!
!           30 format(10f14.8)
!
!        open(10,file='dia.dat',status='unknown',access='append')
!        
!        write(10,*)'Original matrix:'
!        do j=1,counter
!           write(10,30)ss0(:,j)
!        enddo
!
!        write(10,*)
!
!        write(10,*)'S^(-1/2):'
!        do j=1,counter
!           write(10,30)ss1(:,j)
!        enddo
!
!        write(10,*)
!
!        write(10,*)'S^(1/2):'
!        do j=1,counter
!           write(10,30)ss2(:,j)
!        enddo
!        write(10,*)
!        
!        tmp=matmul(ss1,ss2)
!
!        write(10,*)'S^(-1/2)*S^(1/2)'
!        do j=1,counter
!           write(10,30)tmp(:,j)
!        enddo
!        write(10,*)
!
!        tmp=matmul(ss1,matmul(ss0,ss1))
!        write(10,*)'S^(-1/2)*S*S^(-1/2)'
!        do j=1,counter
!           write(10,30)tmp(:,j)
!        enddo
!        write(10,*)
!        
!        close(10)
!    end subroutine checkss12

end subroutine transform_pphi

