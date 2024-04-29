! maps grd to grda; m's must be multiples of n's.
! grd assumed (not checked!!!!) equispaced.

subroutine increase_3dgrid(gr,nd,nx,ny,nz,gra,mx,my,mz)
  implicit none
  integer, intent(in) :: nd,nx,ny,nz,mx,my,mz
  real*8,  intent(in) :: gr( nx,ny,nz, nd)
  real*8,  intent(out):: gra(mx,my,mz, nd)

  real*8              :: gb( mx,ny,nz,nd)
  real*8              :: gc( mx,my,nz,nd)
  
  real*8              :: dd
  integer             :: pp, p,i


  if(mod(mx,nx)/=0.or.mod(my,ny)/=0.or.mod(mz,nz)/=0.or.min(nx,ny,nz).le.1) &
       then
     write(6,*)' sizes mismatches in cgrdirefine; smaller grid ',nx,ny,nz
     write(6,*)' is not a divisor of the larger grid ',mx,my,mz
     write(6,*)' or grid sizes have 1s, which is a pain '
     stop
  endif

  dd=gr(2,1,1,1)-gr(1,1,1,1)
  pp=mx/nx
  do i=1,nx
     do p=1,pp
        gb((i-1)*pp+p,:,:,:) = gr(i,:,:,:)
        gb((i-1)*pp+p,:,:,1) = gr(i,:,:,1)+dd/pp * (p-1)
     enddo
  enddo

  dd = gr(1,2,1,2)-gr(1,1,1,2)
  pp = my/ny
  do i=1,ny
     do p=1,pp
        gc(:,(i-1)*pp+p,:,:) = gb(:,i,:,:)
        gc(:,(i-1)*pp+p,:,2) = gb(:,i,:,2)+dd/pp*(p-1)
     enddo
  enddo

  dd = gr(1,1,2,3)-gr(1,1,1,3)
  pp = mz/nz
  do i=1,nz
     do p=1,pp
        gra(:,:,(i-1)*pp+p,:) = gc(:,:,i,:)
        gra(:,:,(i-1)*pp+p,3) = gc(:,:,i,3)+dd/pp*(p-1)
     enddo
  enddo

end subroutine increase_3dgrid
     


subroutine increase_3dwf(cf,nx,ny,nz,cz,mx,my,mz,lot)
  !
  ! increase cf(nx,ny,nz) to cz(nx,ny,nz); the m's must be integ. multip. of n
  !
  implicit none
  integer nx,ny,nz
  integer mx,my,mz
  integer lot
  complex*16 cf(nx,ny,nz,lot)
  complex*16 cz(mx,my,mz,lot)

  complex*16 ca(mx,ny,nz,lot)
  complex*16 cb(mx,my,nz,lot)

  if(mod(mx,nx)/=0.or.mod(my,ny)/=0.or.mod(mz,nz)/=0.or.min(nx,ny,nz)<1) then
     write(6,*)' sizes mismatches in cgrdirefine; smaller grid ',nx,ny,nz
     write(6,*)' is not a divisor of the larger grid ',mx,my,mz
     stop
  endif

  call increase_singledirecwf(cf, ca, 1,     nx, mx/nx, ny*nz*lot) 
  call increase_singledirecwf(ca, cb, mx,    ny, my/ny,    nz*lot)
  call increase_singledirecwf(cb, cz, mx*my, nz, mz/nz,       lot)
end subroutine increase_3dwf

subroutine increase_singledirecwf(cs, cr, nl, n, k, nr)
  implicit none
  integer nl, n, k, nr,i,j
  complex*16 cs(nl, n, nr)
  complex*16 cn(nl, n, nr)
  complex*16 ct(nl, n, nr)
  complex*16 cr(nl, k, n, nr)
  complex*16, parameter :: ci = (0.d0,1.d0)
  real*8 pi
  complex*16 cterm


  cn=cs
  pi = dacos(-1.d0)

  write(6,*)' nl,n,k,nr ',nl,n,k,nr

  if(k==1) then
     cr(:,1,:,:)=cs
     return
  end if

  
  write(6,*)' step 2 '
  write(6,*)' sum(cn**2) ',sum(abs(cn**2))

  call fftka( nl,n,  nr, cn)
  
  write(6,*)' post sum(cn**2) ',sum(abs(cn**2))

  do j=1,k
     ct = cn
     do i=1,n
        if(i.le.n/2) then
           cterm = exp(ci*2*pi* (j-1)/dble(n*k)* (i-1))
        else if(i>n/2+1) then
           cterm = exp(ci*2*pi* (j-1)/dble(n*k)* (i-1-n)) 
        else
           cterm = 0d0
        end if
        ct(:,i,:)= cterm* ct(:,i,:)
     enddo
     ct= conjg(ct)
     call fftka(nl,n,nr, ct)
     ct = conjg(ct)
     cr(:,j,:,:) = ct/n
  enddo

  write(6,*)' post sum(cr**2)/k ',sum(abs(cr**2))/k

end subroutine increase_singledirecwf

!program increase_grid_drv_prg
!  call increase_grid_drv
!end program increase_grid_drv_prg

subroutine increase_grid_drv
  implicit none
  integer, parameter :: nx=8,ny=16,nz=32,mx=nx,my=2*ny,mz=2*nz
  real*8 y,dy,ymin,y0,suma,sumd,sumg,x,dx,x0,xmin,z,dz,z0,zmin
  integer iy,ix,iz
  complex*16 cin(nx,ny,nz),cout(mx,my,mz)
  real*8 gr(nx,ny,nz,3), grout(mx,my,mz,3)
  complex*16, parameter :: ci=(0.d0,1.d0)

  xmin=-2.3
  ymin=0.7
  zmin=0.5
  dx=0.35
  dy=0.4
  dz=0.43
  x0 = xmin+(nx/3)*dx
  y0 = ymin+(ny/4+2)*dy
  z0 = zmin+((nz*5)/8+1)*dz

  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           x=xmin+(ix-1)*dx
           y=ymin+(iy-1)*dy
           z=zmin+(iz-1)*dz

           gr(ix,iy,iz,:)=(/x,y,z/)
           
           cin(ix,iy,iz)= fc(x,y,z)
           write(600+(iz-1)*mz/nz,*)x,y,dble( cin(ix,iy,iz))
           write(700+(iz-1)*mz/nz,*)x,y,aimag(cin(ix,iy,iz))
        enddo
     enddo
  enddo
  call increase_3dwf(cin,nx,ny,nz,   cout, mx,my,mz,1)
  call increase_3dgrid(gr,3,nx,ny,nz,grout, mx,my,mz)
  suma = 0d0
  sumd = 0d0
  sumg = 0d0
  do iz=1,mz
     do iy=1,my
        do ix=1,mx
           x=xmin+(ix-1)*dx*(nx/dble(mx))
           y=ymin+(iy-1)*dy*(ny/dble(my))
           z=zmin+(iz-1)*dz*(nz/dble(mz))
           
           write(800+iz-1,*)x,y,dble( cout(ix,iy,iz))
           write(900+iz-1,*)x,y,aimag(cout(ix,iy,iz))
           
           suma = suma + abs(fc(x,y,z))**2
           sumd = sumd + abs(fc(x,y,z)-cout(ix,iy,iz))**2
           sumg = sumg + sum((grout(ix,iy,iz,:)-(/x,y,z/))**2)
        enddo
     enddo
  enddo
  write(6,*)' suma, sumd , sumg',suma,sumd,sumg
  
contains
  complex*16 function fc(x,y,z)
    implicit none
    real*8 x,y,z
    fc = exp(-0.5*(x-x0)**2-0.6*(y-y0)**2-0.52*(z-z0)**2)*exp(ci*(0.4*x+0.8*y+0.3*z))
  end function fc

end subroutine increase_grid_drv
  
     
  
