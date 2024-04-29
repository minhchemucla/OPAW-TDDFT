subroutine vxc_lda_libxc(dens, vxc, n, nspin) 
  implicit none
  integer n, nspin, st
  real*8  dens(n), vxc(n)
  real*8, allocatable :: exc_lda(:)
  allocate(exc_lda(n), stat=st); call check0(st,' exc_lda ')
  call lda_libxc(dens, n, nspin, vxc, exc_lda)
  deallocate(exc_lda)
end subroutine vxc_lda_libxc
