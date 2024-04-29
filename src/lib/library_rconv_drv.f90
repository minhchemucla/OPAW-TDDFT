!program check_rconv
!  call rconv_drv
!end program check_rconv

subroutine rconv_drv
  implicit none
  integer, parameter :: nx=8,ny=6,nz=3
  integer            :: ix, iy, iz
   real*8             :: gin( nx,ny,nz)
  real*8             :: bk(  nx,ny,nz)
  real*8             :: gout(nx,ny,nz)
  real*8             :: gana(nx,ny,nz)
  
  call gin_bk_set
  call rconv_ana
  call rconv(gin,bk,gout,nx,ny,nz)
   
  write(6,*)' gout, gana, dif ',sum(abs(gout)),sum(abs(gana)),sum(abs(gout-gana))

contains

  subroutine gin_bk_set
    implicit none
    real*8 x, y, z
    do ix=1,nx
       do iy=1,ny
          do iz=1,nz
             x = ix-dble(nx)/2
             y = iy-dble(ny)/2
             z = iz-dble(nz)/2
             gin(ix,iy,iz)=exp(-(x+0.4)**2*0.5-(x-1.9)*(y-0.7)*0.4+(x-y)*(x-0.5)*0.2-(z+y)**2*0.4)
          enddo
       enddo
    enddo
    
    do ix=1,nx
       do iy=1,ny
          do iz=1,nz
             x = ix-1d0
             y = iy-1d0
             z = iz-1d0
             if(ix>nx/2+1)x=x-dble(nx)
             if(iy>ny/2+1)y=y-dble(ny)
             if(iz>nz/2+1)z=z-dble(nz)
             bk(ix,iy,iz)= 1/max(x**2+0.8*y**2+0.5*z**2,0.1)
          enddo
       enddo
    enddo

  end subroutine gin_bk_set

  subroutine rconv_ana
    implicit none
    complex*16 gk(nx,ny,nz)
    gk = gin
    call fft_ana(gk)
    gk = gk*bk 
    gk = conjg(gk)
    call fft_ana(gk)
    gk = conjg(gk)
    gana = gk/dble(nx)/dble(ny)/dble(nz)
  end subroutine rconv_ana

  subroutine fft_ana(gk)
    implicit none
    integer jx, jy, jz
    real*8 pi
    complex*16 gk(nx,ny,nz)
    complex*16 hk(nx,ny,nz)
    complex*16, parameter :: ci=(0d0,1d0)

    pi = dacos(-1d0)

    hk = 0d0
    do ix=1,nx
       do iy=1,ny
          do iz=1,nz
             do jx=1,nx
                do jy=1,ny
                   do jz=1,nz
                      hk(     ix,iy,iz) = &
                           hk(ix,iy,iz)+&
                           gk(jx,jy,jz)*&
                          exp(-2d0*pi*ci/dble(nx)*(ix-1)*(jx-1))*&
                          exp(-2d0*pi*ci/dble(ny)*(iy-1)*(jy-1))*&
                          exp(-2d0*pi*ci/dble(nz)*(iz-1)*(jz-1))
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    gk = hk
  end subroutine fft_ana
end subroutine rconv_drv
