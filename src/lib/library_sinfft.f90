subroutine sinfft_1d(g, n) 
  ! g(r)-->sum(sin(kr) g(r)), k=0,dk,..,(nf-1)*dk, and dk= pi/(nf*dr). Not 2pi!
  ! r and k start from 0
  implicit none
  integer n,st,i,m
  real*8 g(0:n-1), d
  complex*16, allocatable :: ca(:)
  
  allocate(ca(0:n*2-1),stat=st); if(st/=0) stop ' ca in sinfft '
  ca(0:n-1)=g
  ca(n)=0d0
  do i=1,n-1   
     ca(2*n-i) = -ca(i)   !  0 1 2 3 4->zero 5=-3 6=-2 7=-1 
  enddo

  m=nint(dlog(dble(2*n))/dlog(2d0))
  call check(2**m,2*n,' 2^m,2*n ')

  call fftsa(2*n,ca,m)
  d = maxval(abs(dble(ca)))
  if(d>1d-9) stop ' problem in ca sinfft '

  g = -aimag(ca(0:n-1))/2d0
  deallocate(ca)
end subroutine sinfft_1d

