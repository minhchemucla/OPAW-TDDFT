!
! u(1:n,1:lot) = inv(matrix) r(1:n,1:lot), 
!  where matrixtridiag, with off_diag and diag, an symmetric.  note that off_diag has dim. n-1.
!
      subroutine c_symm_tridiag_lot(off_diag, diag, r, u, n, lot)
      implicit none

      integer n,j,lot
      
      complex*16 bet
      complex*16 off_diag(n-1),diag(n),r(n,lot),u(n,lot)!could be off_diag(n)
      complex*16 ,allocatable :: gam(:)

      allocate(gam(n),stat=j); if(j/=0) stop ' tridiag problems '
      if(abs(diag(1))==0.d0)then
         write(6,*)' problem, diag(1) ',diag(1)
         stop
      endif

      bet=diag(1)
      u(1,:)=r(1,:)/bet
      do 11 j=2,n
        gam(j)=off_diag(j-1)/bet
        bet=diag(j)-off_diag(j-1)*gam(j)
        if(bet.eq.0.)then
           write(6,*)' bet ',bet
           stop
        endif
           
        u(j,:)=(r(j,:)-off_diag(j-1)*u(j-1,:))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j,:)=u(j,:)-gam(j+1)*u(j+1,:)
12    continue

      deallocate(gam)
      end

      subroutine c_tridag(a,b,c,r,u,n)
      implicit none

      integer n,j
      integer, parameter :: nmax=100
      
      complex*16 bet
      complex*16 gam(nmax),a(n),b(n),c(n),r(n),u(n)

      if(b(1)==0.d0.or.n>nmax)then
         write(6,*)' problem, b(1),n,nmax ',b(1),n,nmax
         stop
      endif

      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.0) then
           write(6,*)' bet  ',bet
           stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      end

c     seems that a,b,c are the tridiag coeff.; 
c     a below, b diag, c above
c     Note that a(1) and c(n) are not used; i.e., matrix is a(j) = M(j,j-1),
c     b(j) = M(j,j), c(j) = M(j,j+1)

      subroutine r_tridag(a,b,c,r,u,n)
      implicit none

      integer n,j
      
      real*8 bet
      real*8 gam(n),a(n),b(n),c(n),r(n),u(n)

      if(b(1)==0.d0)then
         write(6,*)' problem, b(1) ',b(1)
         stop
      endif

      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)then
           write(6,*)' bet '
           stop
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      end

c     Correced for a, c being of order n-1!!!
c     a,b,c are the tridiag coeff.; 
c     a below, b diag, c above
c     Note that a(1) and c(n) are not used; i.e., matrix is a(j) = M(j,j-1),
c     b(j) = M(j,j), c(j) = M(j,j+1)
c
c     Calculates u = M_inv r.

      subroutine r_tridag_correct(a,b,c,r,u,n)
      implicit none
      integer n,j
      real*8 bet,a(2:n),b(n),c(1:n-1),r(n),u(n)
      real*8, allocatable :: gam(:)
      allocate(gam(2:n),stat=j); if(j/=0)stop ' r_tridiag '
      if(b(1)==0.d0)then
         write(6,*)' problem, b(1) ',b(1)
         stop
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)

        if(abs(bet)<1e-12) then
           write(6,*)' bet in tridiag is too small ',bet
           stop
        endif

        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      deallocate(gam)
      return
      end

c     Correced for a, c being of order n-1!!!
c     a,b,c are the tridiag coeff.; 
c     a below, b diag, c above
c     Note that a(1) and c(n) are not used; i.e., matrix is a(j) = M(j,j-1),
c     b(j) = M(j,j), c(j) = M(j,j+1)
c
c     Calculates u = M_inv r.

      subroutine c_tridag_correct(a,b,c,r,u,n)
      implicit none
      integer n,j
      complex*16 bet, a(2:n),b(n),c(1:n-1),r(n),u(n)
      complex*16, allocatable :: gam(:)
      allocate(gam(2:n),stat=j); if(j/=0) stop ' tridiag correct '

      if(b(1)==0.d0)then
         write(6,*)' problem, b(1) ',b(1)
         stop
      endif

      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)

        if(abs(bet)<1e-12) then
           write(6,*)' bet in tridiag is too small ',bet
           stop
        endif

        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      deallocate(gam)
      return
      end











