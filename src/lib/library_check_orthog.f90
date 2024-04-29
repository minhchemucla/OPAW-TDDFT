subroutine check_c_orthog(cvec, N)
  implicit none
  integer N, j, i
  complex*16 cvec(N, N)
  real*8 delta; external delta

  complex*16, allocatable, dimension(:, :) :: ctemp
  allocate(ctemp(N, N), stat=j); call check(j,0,' check_c_orth ')

  ctemp = matmul(cvec, transpose(cvec))
  do j=1, N
     do i=1, N
        if(abs(ctemp(i, j)-delta(i, j)) > 1.d-6) then
           write(6,*)' i, j, ctemp ',i, j, ctemp(i, j)
           stop
        endif
     enddo
  enddo

  deallocate(ctemp)

end subroutine check_c_orthog

