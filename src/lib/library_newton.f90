!* prepare f(0), ... f(n)  (array of length n+1)
!* g(0)=(f(1)-f(0))/(x1-x0), g(1)=(f(2)-f(1))/(x1-x0).
!* f=g 
!* next stage 
!* g(0)=(f(1)-f(0))/(x(1+1)-x(0))

module newton
  interface newton_a
     module procedure newton_a_r
     module procedure newton_a_cr
     module procedure newton_a_cc
  end interface
contains

  subroutine newton_a_r(n, r1, r, scale, rf)
    implicit none
    integer k,l,n,st
    
    real*8 r1(0:n-1), r(0:n-1)
    real*8, external :: rf
    real*8 scale
    real*8, allocatable:: ra(:), rb(:)
    
    allocate(ra(0:n-1), rb(0:n-1), stat=st); if(st/=0)stop ' newton r '
    
    do l=0,n-1
       do k=0,n-1-l
          if(l==0) then
              rb(k) = rf(r(k)*scale)
          else
             rb(k) = (ra(k+1)-ra(k))/(r(k+l)-r(k))
          end if
       enddo
       ra(0:n-1-l) = rb(0:n-1-l)
       r1(l)     = ra(0)
    end do
    deallocate(ra,rb)
  end subroutine newton_a_r

  subroutine newton_a_cr(n, c1, x, scale, rf)
    implicit none
    integer k,l,n,st
    
    complex*16 c1(0:n-1)
    real*8 x(0:n-1)
    complex*16, external :: rf
    real*8 scale
    complex*16, allocatable:: ca(:), cb(:)
    
    allocate(ca(0:n-1), cb(0:n-1), stat=st); if(st/=0)stop ' newton r '
    
    do l=0,n-1
       do k=0,n-l
          if(l==0) then
             cb(k) = rf(x(k)*scale)
          else
             cb(k) = (ca(k+1)-ca(k))/(x(k+l)-x(k))
          end if
       enddo
       ca(0:n-1-l) = cb(0:n-1-l)
       c1(l)     = ca(0)
    end do
    deallocate(ca,cb)
  end subroutine newton_a_cr

  subroutine newton_a_cc(n, c1, z, scale, rf)
    implicit none
    integer k,l,n,st
    
    complex*16 c1(0:n-1)
    complex*16 z(0:n-1)
    complex*16, external :: rf
    real*8 scale
    complex*16, allocatable:: ca(:), cb(:)
    
    allocate(ca(0:n-1), cb(0:n-1), stat=st); if(st/=0)stop ' newton r '
    
    do l=0,n-1
       do k=0,n-l
          if(l==0) then
             cb(k) = rf(z(k)*scale)
          else
             cb(k) = (ca(k+1)-ca(k))/(z(k+l)-z(k))
          end if
       enddo
       ca(0:n-1-l) = cb(0:n-1-l)
       c1(l)     = ca(0)
    end do
    deallocate(ca,cb)
  end subroutine newton_a_cc
  
  subroutine newton_points_r(n, r)  
    implicit none
    
    integer,parameter :: mm =500   ! ntot = mm* npoints
    real*8, parameter :: dist2min = 1.d-30
    
    integer              :: k,n,nb,m,j,jp,st
    integer, external    :: is_it_power2
    integer, allocatable :: avail(:)
    real*8               :: r(1:n)
    real*8               :: amax,dist2,rho,pi
    real*8,  allocatable :: rp(:),a(:)
    real*8,     external :: ran_ps

    nb = mm*n
    allocate(rp(1:nb),avail(1:nb),a(1:nb),stat=st); call check(st,0,' rp  ')
    
    pi = dacos(-1d0)
    do j=1,nb
       rp(j) = cos(dble(2*j-1)/dble(2*nb+1)* pi)
    enddo
    
    avail= 1
    r(1) = 1d0
    avail(1) = 0
    A  = 0d0
    do k=1,n-1
       if(is_it_power2(k)==1)write(6,*)' k, r(k) ',k,r(k)
       do j=1,nb
          dist2 = (rp(j)-r(k))**2
          if(avail(j)==1.and.dist2>dist2min) then
             A(j) = A(j) + log(dist2)
          else
             avail(j)=0
          end if
       enddo
       jp=1; Amax = -1e30; 
       do j=1,nb
          if(avail(j)==1.and.A(j)>Amax) then
             Amax = A(j)
             jp = j
          endif
       end do
       r(k+1) = rp(jp)
       avail(jp) = 0
    enddo
    rho = exp(sum(log(abs(r)))/dble(n+1))
    write(6,*)' points scaled (divided) by rho ', rho; call flush(6)
    r = r/rho
    
    open(550,file='pointnewton.txt',status='replace')
    do j=1,n
       write(550,*)j,r(j), ' j, newton_point_j   '
    enddo
    close(550)  
    deallocate(rp,avail,a)
  end subroutine newton_points_r
  
  subroutine newton_read_points_r(n, r)  
    implicit none
    
    integer :: n,j,i
    real*8  :: r(n)
    
    open(550,file='pointnewton.txt',status='old')
    do j=1,n
       read(550,*,end=99)i,r(j)
       if(i/=j) stop ' ij pointnewton '
    enddo
    close(550)  
    return

99  continue
    write(6,*)' stopping; could not reach line ',j,' in pointnewton.txt '
    stop
    
  end subroutine newton_read_points_r
  


  subroutine test_newton
    implicit none
    integer,  parameter :: n=1000
    integer   st
    real*8              :: scale
    real*8, allocatable :: r(:),a(:)
    real*8, external    :: theta_simple
    allocate(r(n),a(n),stat=st)
    call check(st,0,'st newton ')
    scale = 1d0
    call newton_points_r(n,r)
    call newton_a_r(n, a, r, scale, theta_simple)
    call print_newton_coef(a, r, n)
  end subroutine test_newton
  
  subroutine print_newton_coef(a, r, n)
    implicit none
    integer i, n
    real*8 a(n), r(n)
    open(553,file='newton_coef_points.txt')
    do i=1,n
       write(553,*)i,a(i),r(i),' i, a(i), r(i) '
    enddo
    call flush(553)
    close(553)
  end subroutine print_newton_coef
end module newton

function theta_simple(x)
  implicit none
  real*8 theta_simple, x
  theta_simple = 0.5d0 * erfc((x-0.7)/0.01d0)
end function theta_simple
  
!program newtonprog
!  use newton
!  call test_newton
!end program newtonprog
