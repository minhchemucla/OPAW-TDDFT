!Calculates the first derivative of unequally spaced data using Lagrange 
! interpolating polynomials.
! input
! ----------------------------------------------
!   f    : a 1d array of points
!   x    : a 1d array of grid points
!   n    : the number of grid points
! output
! ----------------------------------------------
!   dfdx : the first derivative of f
subroutine first_derivative_lagrange(f, x, n, dfdx)
  implicit none
  integer :: n
  integer :: i
  real*8  :: f(n), x(n), dfdx(n)

  !x=x_1=x_i
  dfdx(1) = (2d0*x(1)-x(2)-x(3))/((x(1)-x(2))*(x(1)-x(3)))*f(1)+&
            (x(1)-x(3))/((x(2)-x(1))*(x(2)-x(3)))*f(2)+&
            (x(1)-x(2))/((x(3)-x(1))*(x(3)-x(2)))*f(3)
  !x=x_n=x_{i+2}
  dfdx(n) = (x(n)-x(n-1))/((x(n-2)-x(n-1))*(x(n-2)-x(n)))*f(n-2)+&
            (x(n)-x(n-2))/((x(n-1)-x(n-2))*(x(n-1)-x(n)))*f(n-1)+&
            (2d0*x(n)-x(n-2)-x(n-1))/((x(n)-x(n-2))*(x(n)-x(n-1)))*f(n)

  do i=2,n-1
    !x=x_i
    dfdx(i) = (x(i)-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)))*f(i-1)+&
              (2d0*x(i)-x(i-1)-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)))*f(i)+&
              (x(i)-x(i-1))/((x(i+1)-x(i-1))*(x(i+1)-x(i)))*f(i+1)
  enddo
end subroutine

