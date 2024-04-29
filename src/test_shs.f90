subroutine test_shs
    use main_mod
    use mat_module
    use rand_seed_mod
    use paw_mod
    implicit none

    real*8 :: eta(nn),xi(nn)
    real*8 :: tmp(nn),tmp1(nx,ny,nz)
    real*8 :: hm(nb,nb)
    real*8 :: pp1(nn),pp2(nn),hpp1(nn),ham(nn,nn),vec(nn,nn)
    real*8 :: eig(nn)

    integer :: ib,jb,ik
    
    write(*,*) 'test_shs'
    write(*,*) '1. the matrix'
    do ik=1,nk_loc
    do ib=1,nb
        call shs(phit(:,:,:,ib,ik),tmp1,ik)
        do jb=1,nb
            hm(ib,jb)=sum(tmp1*phit(:,:,:,jb,ik))*dv
        enddo
        write(*,'(100(1x,f9.5))') hm(ib,:)
    enddo        
    enddo

!    write(*,*) '2. check if shs is Hermitian'
!    call read_initiate_seeds
!    call rand_r(eta,nn,dv)
!    call rand_r(xi,nn,dv)
!
!    call shs(eta,tmp)
!    write(*,*) '<xi|shs|eta>'
!    write(*,*) sum(xi*tmp)*dv
!    call shs(xi,tmp)
!    write(*,*) '<eta|shs|xi>'
!    write(*,*) sum(eta*tmp)*dv

!    write(*,*) '3. check eigenvalues'
!    do ib=1,nn
!        pp1=0d0
!        pp1(ib)=1d0/sqrt(dv)
!        call shs(pp1,hpp1)
!        do jb=1,nn
!            pp2=0d0
!            pp2(jb)=1d0/sqrt(dv)
!            ham(ib,jb)=sum(pp2*hpp1)*dv
!        enddo
!    enddo
!
!    call mat_diag(ham,vec,eig)
!    write(15,*) eig
end subroutine test_shs
