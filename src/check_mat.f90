subroutine check_mat(q)
    use paw_mod
    use atom_mod, only : p => pawinfo
    implicit none

    integer :: i,j,q

    do i=1,p(q)%nstates
        do j=1,p(q)%nstates
            if(p(q)%lstate(i)==p(q)%lstate(j)) then
!               p(q)%mat_s(i,j)      =sum(p(q)%rr*p(q)%rr*p(q)%dr*&
!                   p(q)%phi(i,:)*p(q)%phi(j,:))
!               p(q)%mat_stilde(i,j) =sum(p(q)%rr*p(q)%rr*p(q)%dr*&
!                   p(q)%phitilde(i,:)*p(q)%phitilde(j,:))
               p(q)%mat_sp(i,j)     =sum(p(q)%rr*p(q)%rr*p(q)%dr*&
                   p(q)%ptilde(i,:)*p(q)%ptilde(i,:))
               p(q)%mat_t(i,j)      =sum(p(q)%rr*p(q)%rr*p(q)%dr*&
                   p(q)%phitilde(i,:)*p(q)%ptilde(:,j))
            else
!                p(q)%mat_s(i,j)      =0d0
!                p(q)%mat_stilde(i,j) =0d0
                p(q)%mat_sp(i,j)     =0d0
                p(q)%mat_t(i,j)      =0d0
            endif
        enddo
    enddo

!    call printmat1(p(q)%mat_delta, 'Delta_ij, old')
    
!    p(q)%mat_delta=p(q)%mat_s-p(q)%mat_stilde

!    call printmat1(p(q)%mat_s,     '<phi|phi>')
!    call printmat1(p(q)%mat_stilde,'<phit|phit>')
    call printmat1(p(q)%mat_sp,    '<pt|pt>')
    call printmat1(p(q)%mat_t,     '<phit|pt>')
!    call printmat1(p(q)%mat_delta, 'Delta_ij')
    
contains
    subroutine printmat1(mat,matname)
        implicit none
        
        real*8 :: mat(p(q)%nstates,p(q)%nstates)
        character(len=*) :: matname
        
        write(*,*) 'print matrix: ',matname
        do i=1,p(q)%nstates
            write(*,*) mat(i,:)
        enddo
    end subroutine printmat1
end subroutine check_mat      
