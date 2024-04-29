subroutine cg_coeff(l,m,l1,m1,l2,m2,cg)
!https://en.wikipedia.org/wiki/Table_of_Clebsch%E2%80%93Gordan_coefficients#_j1_=_2,_j2_=_1
    implicit none

    integer :: l,m,l1,m1,l2,m2,lmax,lmin,i
    logical :: zero,legit,flg_fact
    real*8  :: fact,srNum,srDen,cg

    zero=.false.
    legit=.true.
    
    call check_legit

    if(legit) call check_zero

    if(legit .and. .not.zero) then

        srNum = dble(2*l+1)*fact(l+l1-l2)*fact(l-l1+l2)*fact(l1+l2-l)*fact(l+m)*fact(l-m)
        srDen = fact(l1+l2+l+1)*fact(l1-m1)*fact(l1+m1)*fact(l2-m2)*fact(l2+m2)
        cg=0d0
        do i=0,l+m
            flg_fact=.true.
            call check_fact
            if(flg_fact) cg=cg+((((-1d0)**(i+l2+m2))*fact(l2+l+m1-i)*fact(l1-m1+i))/&
                (fact(i)*fact(l-l1+l2-i)*fact(l+m-i)*fact(i+l1-l2-m)))
!            write(*,*) i,flg_fact
        enddo
        cg=sqrt(srNum/srDen)*cg
    endif

    if(legit) then
    endif
contains   
    subroutine check_fact
        implicit none

        if(l-l1+l2-i<0) flg_fact=.false.
        if(l+m-i    <0) flg_fact=.false.
        if(i+l1-l2-m<0) flg_fact=.false.
    end subroutine check_fact

    subroutine check_zero
        implicit none

        if(abs(l1-l2).gt.l .or. l1+l2.lt.l) then
            cg=0d0
            zero=.true.
        endif

        if(m/=m1+m2) then
            cg=0d0
            zero=.true.
        endif

    end subroutine check_zero

    subroutine check_legit
        implicit none

        if(m<-l .or. m>l .or. m1<-l1 .or. m1>l1 .or. m2<-l2 .or. m2>l2) then
            write(*,*) 'error in cg: m,l,m1,l1,m2,l2',m,l,m1,l1,m2,l2
            stop
        endif
!        if(l1>3 .and. l2>3) then
!            write(*,*) 'error in cg: l1,l2 too large:',l1,l2
!            stop
!        endif
    end subroutine check_legit
end subroutine cg_coeff

subroutine cg3(l1,m1,l2,m2,l3,m3,cg)
!integral of product of three real spherical harmonics    
    implicit none

    integer :: l1,m1,l2,m2,l3,m3,k
    real*8  :: cg,pf,cg0,cgg(2,2,2)
    complex*16 :: tmp,s(2,2,2)

    k=0
    if(m1<0) k=k+1
    if(m2<0) k=k+1
    if(m3<0) k=k+1

    if(k==1 .or. k==3) then
        cg=0d0
    else

        pf=sqrt((2d0*l1+1d0)*(2d0*l2+1d0)/(2d0*l3+1d0)/12.5663706144)
        call get_cg
        call get_s
        tmp=sum(cgg*s)*cg0
        if(dimag(tmp)>1d-8) then
            write(*,*) 'tmp not real',tmp
        endif
        cg=pf*dble(tmp)

    endif
contains    
    subroutine get_s
        implicit none

        integer :: i,j,k
        complex*16 :: a(2),b(2),c(2)
        call get_aa(m1,a)
        call get_aa(m2,b)
        call get_aa(m3,c)
        a=conjg(a);b=conjg(b)
!        write(*,*) a
!        write(*,*) b
!        write(*,*) c

        do i=1,2;do j=1,2;do k=1,2
            s(i,j,k)=a(i)*b(j)*c(k)
        enddo;enddo;enddo
!        write(*,*) s(1,1,1),s(2,1,1),s(1,2,1),s(2,2,1),&
!            s(1,1,2),s(2,1,2),s(1,2,2),s(2,2,2)
    end subroutine get_s
    
    subroutine get_cg
        implicit none
        call cg_coeff(l3,0,l1,0,l2,0,cg0)
        call cg_coeff(l3, m3,l1, m1,l2, m2,cgg(1,1,1))
        call cg_coeff(l3, m3,l1, m1,l2,-m2,cgg(1,2,1))
        call cg_coeff(l3, m3,l1,-m1,l2, m2,cgg(2,1,1))
        call cg_coeff(l3, m3,l1,-m1,l2,-m2,cgg(2,2,1))
        call cg_coeff(l3,-m3,l1, m1,l2, m2,cgg(1,1,2))
        call cg_coeff(l3,-m3,l1, m1,l2,-m2,cgg(1,2,2))
        call cg_coeff(l3,-m3,l1,-m1,l2, m2,cgg(2,1,2))
        call cg_coeff(l3,-m3,l1,-m1,l2,-m2,cgg(2,2,2))
!        write(*,*) l3,m3,l1,m1,l2,m2
!        write(*,*) cgg(1,1,1),cgg(2,1,1),cgg(1,2,1),cgg(2,2,1),&
!            cgg(1,1,2),cgg(2,1,2),cgg(1,2,2),cgg(2,2,2)
    end subroutine get_cg
end subroutine cg3

subroutine get_aa(m,a)
    implicit none

    integer :: m
    complex*16 :: a(2),ci
    ci=(0d0,1d0)

    if(m<0) then
        a(1)=ci/sqrt(2d0)
        a(2)=ci/sqrt(2d0)*(-1)**(m+1)
    else if(m==0) then
        a(1)=1d0
        a(2)=0d0
    else if(m>0) then
        a(1)=1/sqrt(2d0)*(-1)**m
        a(2)=1/sqrt(2d0)
    endif
end subroutine get_aa

subroutine w3j_coeff(l1,m1,l2,m2,l3,m3,cg)
!Wigner 3j symbol
    implicit none

    integer :: l1,m1,l2,m2,l3,m3
    real*8  :: cg,fact

    call cg_coeff(l3,-m3,l1,m1,l2,m2,cg)
    cg=cg/sqrt(2d0*l3+1d0)
    if(abs(mod(l1-l2-m3,2))==1) then
        cg=-cg
    endif
end subroutine w3j_coeff    

function fact(n)
    integer :: n
    real*8  :: fact

    fact = 1
    do i = 1, n
        fact = fact *dble(i)
    end do
end function fact

!program test_cg
!    implicit none
!
!    integer :: l,m,l1,l2,m1,m2
!    real*8  :: cg, cg_read, error
!    integer :: i(6)
!
!    real*8  :: y1, y2, y3
!    real*8  :: dtheta, dphi, theta, phi
!    real*8  :: pi=3.14159265359
!    real*8  :: cross
!    integer :: itheta, iphi
!    integer :: ntheta=1000, nphi=1000
!
!!    character(len=5) :: inp
!!    call getarg(1,inp)
!!    read(inp,*) l
!!    call getarg(2,inp)
!!    read(inp,*) m
!!    call getarg(3,inp)
!!    read(inp,*) l1
!!    call getarg(4,inp)
!!    read(inp,*) m1
!!    call getarg(5,inp)
!!    read(inp,*) l2
!!    call getarg(6,inp)
!!    read(inp,*) m2
!open(unit=1,file='out')
!error=0d0
!    do l=0,4;do m=-l,l;do l1=0,4;do m1=-l1,l1;do l2=0,4; do m2=-l2,l2
!    call cg3(l1,m1,l2,m2,l,m,cg)
!
!    dtheta=    pi/ntheta
!    dphi  =2d0*pi/nphi
!
!    cross=0d0
!    do itheta=1,ntheta
!        do iphi=1,nphi
!            theta=itheta*dtheta
!            phi  =iphi  *dphi
!            call sphericalharmonics(l1,m1,theta,phi,y1)
!            call sphericalharmonics(l2,m2,theta,phi,y2)
!            call sphericalharmonics(l,m,theta,phi,y3)
!            cross=cross+y1*y2*y3*sin(theta)
!        enddo
!    enddo
!
!    cross=cross*dtheta*dphi
!
!    write(*,*) l,m,l1,m1,l2,m2,cg,cross
!    error=max(error,abs(cg-cross))
!    enddo;enddo;enddo;enddo;enddo;enddo
!write(*,*) error
!close(1)
!!    open(unit=1,file='output2')
!!
!!    error=0d0
!!
!!    do l=0,4
!!    do m=-5,5
!!    do l1=0,4
!!    do m1=-5,5
!!    do l2=0,4
!!    do m2=-5,5
!!    call cg_coeff(l,m,l1,m1,l2,m2,cg)
!!    read(1,*) i,cg_read
!!    if(abs(cg-cg_read)>1d-8) write(*,*) l,m,l1,m1,l2,m2,abs(cg-cg_read),cg,cg_read,srnum,srden
!!    error=max(error,abs(cg-cg_read))
!!    enddo
!!    enddo
!!    enddo
!!    enddo
!!    enddo
!!    enddo
!!
!!    close(1)
!!    write(*,*) error
!end program test_cg    
