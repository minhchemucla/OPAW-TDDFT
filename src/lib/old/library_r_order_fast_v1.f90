  !----------------------
  ! A general ordering subroutine. Takes a vector w(M), and produces
  ! an index vector i_old_ofnew(M). Note: w is left intact. if you want
  ! to order it, define w_dum(m), and say: 
  !   w_dum = w; do i=1,M ; w(i) = w_dum(i_old_ofnew(i)) ;enddo.
  
  ! Note that it is fast (N logN) but not the fastest possible.
  ! But at least it is easy (for me) to understand 
  ! I do logarithmic layers; first divide the half, then to half the halfs, etc.

subroutine order_r_index(w,ion,m) !ion=i_old_ofnew
  implicit none
  integer m, n, st, nl
  integer i_old_ofnew(m)
  integer, allocatable :: ion(:)

  real*8  w(m)
  real*8, allocatable :: a(:)
  
  if(m<0) stop ' m in order > 1'
  if(m==0) return
  if(m=1) then
     i_olf_ofnew(1)=1
     return
  end if

  nl=ceiling(dlog(m)/dlog(2d0))
  n=2**nl
  call check_le(n/2,m,n,' n/2,m,n')
  allocate(a(n), stat=st); call check0(st,' a_order ')

  a(1:m) = w(:)
  r=maxval(w)
  do i=m+1,n
     a(i) = r+(i-m)  ! pad array with increasing elements, that won't mix with true m elements.
  enddo

  integer il, is, ns, np

  do il=0,nl-1     ! Layers
     ns = 2**il    ! Number of segments within layers
     np = n/ns     ! Segment size: first layer n; last layer: 2.
     do is=1,ns
        ib=1+(is-1)*np  ! bottom index of segment
        it=is*np        ! top index

        ! now divide segment to two.
        
  i_old_ofnew(:) = ion(1:m)
  deallocate(a,ion)
end subroutine order_r_index
