!call check_mat_inv_diag
!end

subroutine check_mat_inv_diag
  use mat_module
  implicit none
  integer, parameter :: N=6
  real*8 A(N, N)
  integer i,j, iplot
  real*8 ran_ps_neg   ! a function that makes random numbers between -1 and 1.
  external ran_ps_neg
  
  do i=1, N
     do j=1, N
        A(i, j) = ran_ps_neg()
     enddo
  enddo
  
  iplot =6  ! plot to terminal-- see plot_mat later

  call check_mat_inv

  A = (A+transpose(A))/2.d0    ! symmetrizes A.  Our diag. routines only
                              ! work for symmetrical A.
  write(6,*)' Symm(A) '
  call plot_mat(A, iplot)

  call check_mat_diag

contains
  subroutine check_mat_inv
    implicit none
    real*8 A_inv(N, N)
    real*8 B(N, N)  ! work array
    call mat_inv(A,A_inv)
    write(6,*)' A '              ; call plot_mat(A,iplot)
    write(6,*)' A_inv '          ; call plot_mat(A_inv, iplot)
    write(6,*)' matmul(A,A_inv) '; B=matmul(A, A_inv); call plot_mat(B, iplot)
    !
    ! Note that I dont call plot_mat(matmul(A,A_inv)); it is safer to 
    ! call subroutines with actual matrices, not with results of operations.
    !
    !
  end subroutine check_mat_inv

  subroutine check_mat_diag
    implicit none
    real*8 A_vc(N, N), A_evl(N)
    real*8 B(N, N)  ! work array
    call mat_diag(A, A_vc, A_evl)

    B= matmul(A, A_vc)
    do i=1, N
       write(6,*)' i ',i,' matmul(A, A_vc(:,i) ) , A_vc(:,i)*A_evl(i) '
       do j=1,N
          write(6,*)B(j,i),A_vc(j,i)*A_evl(i)
       enddo
       write(6,*)
    enddo
  end subroutine check_mat_diag
end subroutine check_mat_inv_diag
  







