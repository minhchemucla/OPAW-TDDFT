subroutine rg_simple(aa, n, cz, cw)
  implicit none
  integer, parameter :: flag_rg_method =1
  integer nm, n, i, ierr, matz
  real*8 AA(n, n),norm2
  real*8, allocatable :: ain(:,:), wr(:), wi(:), z(:,:), fv1(:)
  integer, allocatable :: iv1(:)
  complex*16 cZ(n, n), cw(n)
  complex*16, parameter :: ci=(0.d0,1.d0)
  
!  if(flag_rg_method==1) then
!     call dgeev_simple(aa, n, cz, cw)
!     return
!  end if

  write(6,*)' n= ',n
  allocate(ain(n, n), wr(n), wi(n), z(n, n), fv1(n), iv1(n),stat=i); call check(i,0,' rgalloc  ')

  ain = AA

  nm= n
  matz = 1
  
  wr=0;wi=0
  write(6,*) ' pre rg,n= ',n; call flush(6)
  call rg(nm,n,ain,wr,wi,matz,z,iv1,fv1,ierr)
  write(6,*)' post rg  -- now checks: min,max wr ',minval(wr),maxval(wr),' minmax wi ',minval(wi),maxval(wi)
  

  call check(ierr,0,' ierr   ')
  
  do i=1, n
     if(wi(i)<0.d0) cycle
     if(wi(i)==0.d0) then
        cZ(:, i) = z(:, i)
     else
        call check_le(i, n-1, ' inm1    ')
        cZ(:, i)     = z(:, i) + ci* z(:, i+1)
        cZ(:, i+1)   = z(:, i) - ci* z(:, i+1)
     endif
  end do

  write(6,*)' formulation ok '

  do i=1, n
     norm2 = sum(abs(Cz(:,i))**2)
     if(norm2 < 1.d-14) then
        write(6,*)' cZ:, norm2 ',i,norm2
        stop
     endif
     cZ(:, i)= cZ(:, i)/sqrt(norm2)
  enddo

  cW = wr + ci*wi
  
!  call rg_check  ! can uncheck

  deallocate( ain, wr, wi, z, fv1, iv1)
contains
  subroutine rg_check
    implicit none

    complex*16 cdiff(N, N), cW_mat(N, N)
    real*8     unity(N, N)

    unity  = 0.d0
    cW_mat = 0.d0
    do i=1,N
       unity( i, i) = 1.d0
       cW_mat(i, i) = cw(i)
    enddo

    cdiff = matmul(AA, Cz)- matmul(Cz, Cw_mat)

    if(sum(abs(cdiff))>1.d-8*sum(abs(AA)))then
       write(6,*)' problem in cdiff '
       write(6,*)' sum(abs(cdiff))  ', sum(abs(cdiff))  
       write(6,*)' sum(abs(AA))      ', sum(abs(AA))  
       stop
    endif
  end subroutine rg_check

end subroutine rg_simple

!call rg_simple_driver
!end
subroutine rg_simple_driver
  implicit none
  integer, parameter :: N=3
  real*8 A(N, N)
  complex*16 Cw(N), Cz(n, n)
  integer i

  A = 0.d0
  A(1,1) = 0.3
  A(1,2) = 1.9
  A(2,1) = -1.78
  A(3,3) = 2.50

  call rg_simple(A, n, cZ, cW)

  do i=1,N
     write(6,*)' i, cW ',i,cw(i)
     write(6,*)' i, cZ ',i,cZ(:, i)
  enddo
end subroutine

