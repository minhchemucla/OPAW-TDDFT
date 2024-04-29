! v(r)= sum_j exp(-(r-r_j)**2 b(r_j) = sum_j g(j) b(j)
! va(rk) = sum_j s_kj b(r_j)
! b = s^-1 va
! 


module lagrange_fit_modu
  implicit none
  save
  integer, parameter   :: md = 6
  real*8               :: fctrl(0:2*md)
end module lagrange_fit_modu

subroutine lagrange_fit( xmin, dx, nx, &
                         ymin, dy, ny, &
                         zmin, dz, nz, &
                         va, &
                         x, y, z, &
                         v)
  use lagrange_fit_modu
  implicit none
  integer nx, ny, nz
  real*8  xmin, ymin, zmin
  real*8  dx,   dy,   dz
  real*8  va(nx, ny, nz)
  real*8  x, y, z
  real*8  v ! output

  integer ix, iy, iz, ibz, iby, ibx, id, i, k
  integer, save :: i1=1
  real*8 :: A(-md+1:md,3),pv(3),s(3),p,t,tp,g(-md+1:md,-md+1:md,-md+1:md)
  !
  !
  !
  if(i1==1) then
     i1=-1
     fctrl(0)=1d0
     do i=1,2*md
        fctrl(i) = fctrl(i-1)*dble(i)
     enddo
  end if
  
  ibx = int((x-xmin)/dx)+1; call check_lele(md,ibx,nx-md,' md,ibx,nd-md ')
  iby = int((y-ymin)/dy)+1; call check_lele(md,iby,ny-md,' md,iby,nd-md ')
  ibz = int((z-zmin)/dz)+1; call check_lele(md,ibz,nz-md,' md,ibz,nd-md ')

  s(1) = x- (ibx-1)*dx-xmin
  s(2) = y- (iby-1)*dy-ymin
  s(3) = z- (ibz-1)*dz-zmin

  pv(1) = s(1)/dx
  pv(2) = s(2)/dy
  pv(3) = s(3)/dz

  idloop: do id=1,3
     p = pv(id)
     if(p<0d0.or.p>1d0) stop ' p wrong '
     if(abs(p)<1d-10) then
        A(:,id)=0d0
        A(0,id)=1d0
     elseif(abs(1-p)<1d-10) then
        A(:,id)=0d0
        A(1,id)=1d0
     else
        tp = 1d0
        do k=1,2*md
           tp = tp*(p+dble(md-k))
        enddo
        do k= -md+1,md
           A(k,id)= (-1)**(mod(md+k,2))* tp /(fctrl(md-1+k)*fctrl(md-k)*(p-k) )
        enddo
     endif
  end do idloop


  do iz= -md+1, md
     do iy= -md+1, md
        do ix= -md+1, md
           g(ix,iy,iz)= A(ix,1)*A(iy,2)*A(iz,3)
        enddo
     enddo
  enddo

  v = sum( va(ibx-md+1:ibx+md, iby-md+1:iby+md, ibz-md+1:ibz+md )* g)
end subroutine lagrange_fit
  
