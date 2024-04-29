! v(r)= sum_j exp(-(r-r_j)**2 b(r_j) = sum_j g(j) b(j)
! va(rk) = sum_j s_kj b(r_j)
! b = s^-1 va
! 


module gauss_fit_modu
  implicit none
  save
  integer, parameter   :: md = 3
  real*8,  parameter   :: sigma_i = 0.4d0, toll=1d-6
  real*8,  allocatable :: S_inv(:,:), g(:,:,:), b(:), F(:,:)
end module gauss_fit_modu

subroutine gauss_fit( xmin, dx, nx, &
                      ymin, dy, ny, &
                      zmin, dz, nz, &
                      va, &
                      x, y, z, &
                      v)
  use gauss_fit_modu
  implicit none
  integer nx, ny, nz, M,st
  real*8  xmin, ymin, zmin
  real*8  dx,   dy,   dz
  real*8  va(nx, ny, nz)
  real*8  x, y, z
  real*8  v ! output

  integer ix, iy, iz, ibz, iby, ibx,i
  integer, save :: i1=1
  real*8 :: rd(3), rd_mx, rd_mn
  real*8, allocatable :: w(:,:,:)

  M= 8*md*md*md
  if(i1==1) then
     allocate(S_inv(M,M), F(M,M), g(2*md,2*md,2*md),b(M), stat=st ); 
     if(st/=0) stop ' S_inv '
  endif
  allocate(w(2*md,2*md,2*md), stat=st );  if(st/=0) stop ' w(gauss) '
  !
  !
  !
  
  ibz = int((z-zmin)/dz)+1; call check_lele(md,ibz,nz-md,' md,ibz,nd-md ')
  iby = int((y-ymin)/dy)+1; call check_lele(md,iby,ny-md,' md,iby,nd-md ')
  ibx = int((x-xmin)/dx)+1; call check_lele(md,ibx,nx-md,' md,ibx,nd-md ')

  rd_mx = 0d0
  rd_mn = 0d0

  !write(6,*)' ibx, iby, ibz ',ibx,iby, ibz

  do iz= ibz-md+1, ibz+md
     do iy= iby-md+1, iby+md
        do ix= ibx-md+1, ibx+md

           rd = (/ x-(ix-1)*dx-xmin , y-(iy-1)*dy-ymin, z-(iz-1)*dz-zmin /)
           g(ix-ibx+md,iy-iby+md,iz-ibz+md)= ffit(rd)
           w(ix-ibx+md,iy-iby+md,iz-ibz+md)= wovrlp(rd)

           rd_mx = max(rd_mx, maxval(rd))
           rd_mn = min(rd_mn, minval(rd))
        enddo
     enddo
  enddo
  !write(6,*)' rd_mx, rd_mn ',rd_mx, rd_mn

  call make_gauss_inv_fit(w)
  if(i1==1) then
     i1=-1
  end if
  
  b = matmul(S_inv, reshape( va(ibx-md+1:ibx+md, iby-md+1:iby+md, ibz-md+1:ibz+md ), (/ M /) ))   
  v = sum(b* reshape(g, (/M/) )) 

  deallocate(w)
contains
  subroutine  make_gauss_inv_fit(w)
    use mat_module
    implicit none

    integer ii, jj, jx, jy, jz,st
    real*8,allocatable ::  S(:,:), s_evc(:,:),s_evl(:)
    real*8 w(M)
    
    allocate(S(M,M), S_evc(M,M), s_evl(M), stat=st ); if(st/=0) stop ' S '

    if(i1==1) then
       ii= 0
       do iz=1,2*md
          do iy=1,2*md
             do ix=1,2*md
                ii=ii+1
                
                jj = 0
                do jz=1,2*md
                   do jy=1,2*md
                      do jx=1,2*md
                         jj = jj+1
                         
                         rd = (/ (ix-jx)*dx , (iy-jy)*dy, (iz-jz)*dz /)
                         F(ii, jj) = ffit(rd)
                         
                      enddo
                   enddo
                enddo
                
             enddo
          enddo
       enddo
    end if

    S= F
    do i=1,M
       S(:,i) = S(:,i)* w(i)
    enddo
    S=matmul(S,F)
    do i=1,M
       S(i,i)=S(i,i)+toll
    enddo
    
    !call mat_diag(S, S_evc, S_evl)
    !write(6,*)' s_evl ',s_evl 

    !S_inv = S ! erase
    call mat_inv(S,S_inv)
    S_inv = matmul(S_inv, F)
    do i=1,M
       S_inv(:,i) = S_inv(:,i)*w(i)
    enddo

    deallocate(s, s_evc, s_evl)
  end subroutine make_gauss_inv_fit

  real*8 function ffit(rd)
    implicit none
    real*8 rd(3)
    ffit = exp(-sum(rd*rd)/2d0/sigma_i**2)
  end function ffit

  real*8 function wovrlp(rd)
    implicit none
    real*8 rd(3)
    wovrlp = exp(-1d0*sum(rd*rd)/2d0/sigma_i**2)
  end function wovrlp
end subroutine gauss_fit
  


subroutine gauss_fit_old( xmin, dx, nx, &
                      ymin, dy, ny, &
                      zmin, dz, nz, &
                      va, &
                      x, y, z, &
                      v)
  use gauss_fit_modu
  implicit none
  integer nx, ny, nz
  real*8  xmin, ymin, zmin
  real*8  dx,   dy,   dz
  real*8  va(nx, ny, nz)
  real*8  x, y, z
  real*8  v ! output

  integer ix, iy, iz, ibz, iby, ibx
  integer, save :: i1=1
  real*8 :: rd(3), rd_mx, rd_mn
  !
  !
  !
  if(i1==1) then
     i1=-1
     call make_gauss_inv_fit
  end if
  
  ibz = int((z-zmin)/dz)+1; call check_lele(md,ibz,nz-md,' md,ibz,nd-md ')
  iby = int((y-ymin)/dy)+1; call check_lele(md,iby,ny-md,' md,iby,nd-md ')
  ibx = int((x-xmin)/dx)+1; call check_lele(md,ibx,nx-md,' md,ibx,nd-md ')

  rd_mx = 0d0
  rd_mn = 0d0

  !write(6,*)' ibx, iby, ibz ',ibx,iby, ibz

  do iz= ibz-md+1, ibz+md
     do iy= iby-md+1, iby+md
        do ix= ibx-md+1, ibx+md

           rd = (/ x-(ix-1)*dx-xmin , y-(iy-1)*dy-ymin, z-(iz-1)*dz-zmin /)
           g(ix-ibx+md,iy-iby+md,iz-ibz+md)= exp(-sum(rd*rd)/2d0/sigma_i**2)

           rd_mx = max(rd_mx, maxval(rd))
           rd_mn = min(rd_mn, minval(rd))
        enddo
     enddo
  enddo

  !write(6,*)' rd_mx, rd_mn ',rd_mx, rd_mn
  
  b = matmul(S_inv, reshape( va(ibx-md+1:ibx+md, iby-md+1:iby+md, ibz-md+1:ibz+md ), (/ 8*md*md*md /) ))   
  v = sum(b* reshape(g, (/8*md*md*md/) )) 

contains
  subroutine  make_gauss_inv_fit
    use mat_module
    implicit none

    integer ii, jj, jx, jy, jz,st
    real*8,allocatable :: S(:,:), S_evc(:,:),S_evl(:)

    allocate(S_inv(8*md*md*md,8*md*md*md), g(2*md,2*md,2*md), b(8*md*md*md), stat=st ); 
    if(st/=0) stop ' S_inv '

    allocate(S(8*md*md*md,8*md*md*md), stat=st ); if(st/=0) stop ' S '
    allocate(S_evc(8*md*md*md,8*md*md*md), stat=st ); if(st/=0) stop ' S '
    allocate(S_evl(8*md*md*md), stat=st ); if(st/=0) stop ' S '

    ii= 0
    do iz=1,2*md
       do iy=1,2*md
          do ix=1,2*md
             ii=ii+1
             
             jj = 0
             do jz=1,2*md
                do jy=1,2*md
                   do jx=1,2*md
                      jj = jj+1
                      
                      rd = (/ (ix-jx)*dx , (iy-jy)*dy, (iz-jz)*dz /)
                      S(ii, jj) = exp(-sum(rd*rd)/2d0/sigma_i**2)
                      
                   enddo
                enddo
             enddo
             
          enddo
       enddo
    enddo
    
    call mat_diag(S, S_evc, S_evl)
    write(6,*)' S_evl_min_max ',real(minval(S_evl)),real(maxval(S_evl))
    
    S_inv = S_evc
    do ii=1,size(S_inv,2)
       if(S_evl(ii)>toll) then
          S_inv(:,ii) = S_inv(:,ii)/S_evl(ii)
       else
          S_inv(:,ii) = 0d0
       end if
    enddo
    S_inv = matmul(S_inv, transpose(S_evc))

    deallocate(s, s_evc, s_evl)
  end subroutine make_gauss_inv_fit
end subroutine gauss_fit_old
  
