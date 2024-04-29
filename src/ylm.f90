!return in y Ylm(theta,phi)
!Y normalized such dint YlmYlm d(costheta)dphi =1
!subroutine sphericalharmonics(l,m,theta,phi,y) !for testing 
subroutine sphericalharmonics(l,m,xx,yy,zz,y)
    implicit none

    real*8  :: y !output

    real*8  :: theta,phi
    real*8  :: xx,yy,zz !x/r,y/r,z/r
    integer :: l,m

    real*8  :: pi=3.14159265359
    real*8  :: pref(3)

    pref(1)=1d0/(2d0*sqrt(pi))
    pref(2)=pref(1)*sqrt(3d0)
    pref(3)=sqrt(15d0)/(2d0*sqrt(pi))

!    xx=sin(theta)*cos(phi)
!    yy=sin(theta)*sin(phi)
!    zz=cos(theta)

    sl : select case (l) 
    case (0)
        if(m/=0) then
            write(*,*) 'm should be 0 for l=0, instead m= ',m
            stop
        endif
        y=pref(1)
    case (1)
        sm1 : select case (m)
        case(-1)
            y=pref(2)*yy
        case(0)
            y=pref(2)*zz
        case(1)
            y=pref(2)*xx
        case default
            write(*,*) 'm should be -1 to 1 for l=1, instead m= ',m
            stop
        end select sm1
    case(2)
        sm2 : select case (m)
        case(-2)
            y=pref(3)*xx*yy
        case(-1)
            y=pref(3)*yy*zz
        case(0)
            y=pref(3)*(2d0*zz**2-xx**2-yy**2)/(2d0*sqrt(3d0))
        case(1)
            y=pref(3)*xx*zz
        case(2)
            y=pref(3)*(xx**2-yy**2)/2d0
        case default
            write(*,*) 'm should be -2 to 2 for l=2, instead m= ',m
            stop
        end select sm2
    case(3)
        sm3 : select case(m)
        case(-3)
            y=1d0/4d0*sqrt(35d0/2d0/pi)*yy*(3d0*xx*xx-yy*yy)
        case(-2)
            y=1d0/2d0*sqrt(105d0/pi)*xx*yy*zz
        case(-1)
            y=1d0/4d0*sqrt(21d0/2d0/pi)*yy*(4d0*zz*zz-xx*xx-yy*yy)
        case(0)
            y=1d0/4d0*sqrt(7d0/pi)*zz*(2d0*zz*zz-3d0*xx*xx-3d0*yy*yy)
        case(1)
            y=1d0/4d0*sqrt(21d0/2d0/pi)*xx*(4d0*zz*zz-xx*xx-yy*yy)
        case(2)
            y=1d0/4d0*sqrt(105d0/pi)*(xx*xx-yy*yy)*zz
        case(3)
            y=1d0/4d0*sqrt(35d0/2d0/pi)*xx*(xx*xx-3d0*yy*yy)
        case default
            write(*,*) 'm should be -3 to 3 for l=3, instead m= ',m
            stop
        end select sm3
    case(4)
        sm4 : select case(m)
        case(-4)
            y=3d0/4d0*sqrt(35d0/pi)*xx*yy*(xx*xx-yy*yy)
        case(-3)
            y=3d0/4d0*sqrt(35d0/2d0/pi)*(3d0*xx*xx-yy*yy)*yy*zz
        case(-2)
            y=3d0/4d0*sqrt(5d0/pi)*xx*yy*(7d0*zz*zz-1d0)
        case(-1)
            y=3d0/4d0*sqrt(5d0/2d0/pi)*yy*zz*(7d0*zz*zz-3d0)
        case(0)
            y=3d0/16d0/sqrt(pi)*(35d0*zz**4-30d0*zz*zz+3d0)
        case(1)
            y=3d0/4d0*sqrt(5d0/2d0/pi)*xx*zz*(7d0*zz*zz-3d0)
        case(2)
            y=3d0/8d0*sqrt(5d0/pi)*(xx*xx-yy*yy)*(7d0*zz*zz-1d0)
        case(3)
            y=3d0/4d0*sqrt(35d0/2d0/pi)*(xx*xx-3d0*yy*yy)*xx*zz
        case(4)
            y=3d0/16d0*sqrt(35d0/pi)*(xx*xx*(xx*xx-3d0*yy*yy)-yy*yy*(3*xx*xx-yy*yy))
        case default
            write(*,*) 'm should be -4 to 4 for l=4, instead m= ',m
            stop
        end select sm4
    case default
        stop 'l>4 not implemented yet'
    end select sl        

end subroutine sphericalharmonics

!program test_ylm
!    implicit none
!
!    real*8  :: y1, y2
!    real*8  :: dtheta, dphi, theta, phi
!    real*8  :: pi=3.14159265359
!    real*8  :: cross,error
!    
!    integer :: itheta, iphi, l1, m1, l2, m2
!    integer :: ntheta=1000, nphi=1000
!
!    dtheta=    pi/ntheta
!    dphi  =2d0*pi/nphi
!
!    do l1=0,4
!    do m1=-l1,l1
!    do l2=0,4
!    do m2=-l2,l2
!
!        cross=0d0
!        do itheta=1,ntheta
!            do iphi=1,nphi
!                theta=itheta*dtheta
!                phi  =iphi  *dphi
!                call sphericalharmonics(l1,m1,theta,phi,y1)
!                call sphericalharmonics(l2,m2,theta,phi,y2)
!                cross=cross+y1*y2*sin(theta)
!            enddo
!        enddo
!
!        cross=cross*dtheta*dphi
!        if(l1==l2.and.m1==m2) then
!            error=abs(cross-1d0)
!        else
!            error=abs(cross)
!        endif
!
!        if(error.gt.1d-6) then
!            write(*,*) 'cross terms with l1,m1,l2,m2: ',l1,m1,l2,m2,cross
!        endif
!    enddo
!    enddo
!    enddo
!    enddo
!            
!end program test_ylm
!
