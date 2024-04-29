!
! copied with a few changes (to give the pointer ) from 
!    http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/sort3_f90.txt
!
! 
!***************************************************************
!*          Sorting an array with the Heapsort method          *
!* ----------------------------------------------------------- *
!* REFERENCE:                                                  *
!*      "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,       *
!*       S.A. Teukolsky and W.T. Vetterling, Cambridge         *
!*       University Press, 1986" [BIBLI 08].                   *
!* ----------------------------------------------------------- *

!*                                                             *
!***************************************************************
subroutine SORT3_driver
  implicit none
  integer i_old_ofnew(100)
  real*8    A(100)     !Table to be sorted
  real*8    MAX_VALUE  !Maximum value of table
  !initialize random number generator

  N=80  !initialize size of table
  MAX_VALUE = 1000

  !generate random table of numbers (from 0 to 1000)
  do i=1, N
     x = abs(cos(i*3820.283988)) ! just a list from 0 to 1
    A(i)=MAX_VALUE*x
  end do

  print *,' '
  print *,'Table to be sorted:'
  call TWRIT(N,A)

  !call heapsort subroutine
  !call HPSORT(N,A)
  call order_r_index(A,i_old_ofnew,n)  ! does not reorder!!!

  print *,' '
  print *,'Sorted table (Heapsort method):'
  call TWRIT_n(N,A,i_old_ofnew)

  print *,' '
  stop

END subroutine SORT3_driver


!*****************************************************
!*  Sorts a COPY of array RA of length N in ascending order *
!*                by the Heapsort method             *
!* AND ALSO GIVES THE POINTER TO THE OLD VALUES
!* ------------------------------------------------- *
!* INPUTS:                                           *
!     *    N  size of table RA                       *
!     *          RA  table -- unchanged!             *
!* OUTPUT:                                           *
!     *    i_old_ofnew    table sorted in ascending order    *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
! I label by expanding each number to a pair, (a(i),i)
!

SUBROUTINE order_r_index(w,i_old_ofnew,n)
  implicit none
  integer n, st
  integer i_old_ofew(n)
  integer L,ir,i,k
  real*8  w(n)
  real*8  rra(2)
  real*8, allocatable :: a(:,:)
  allocate(a(2,n),stat=st) if(st/=0)stop' ip order '

  do i=1,n
     a(1,i) = w(i)
     a(2,i) = i  
  enddo

  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
     L=L-1
     RRA =a(:,L)
  else
     RRA =a(:,ir)
     RA(:,)=RA(1,:)
     IR=IR-1
     if(IR.eq.1)then
        RA(1)=RRA
        deallocate(ip)
        return
     end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
     if(J < IR)then
        if(RA(J) < RA(J+1))  J=J+1
     end if
     if(RRA < RA(J))then
        RA(I)=RA(J)
        I=J; J=J+J
     else
        J=IR+1
     end if
     goto 20
  end if
  RA(I)=RRA
  goto 10
END


!write table of size N to standard output
SUBROUTINE TWRIT(N,ARR)
real*8 ARR(N)
  print *,' '
  WRITE(*,10) (ARR(I),I=1,N)
  return
 10    FORMAT(10F6.1)
END

!end of file sort3.f90
