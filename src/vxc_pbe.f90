subroutine vxc_pbe_new(dens, vxc, n, nspin) ! remove
  implicit none
  integer n, nspin, st
  real*8  dens(n), vxc(n)
  real*8, allocatable :: exc_pbe(:)
  allocate(exc_pbe(n), stat=st); call check0(st,' exc_pbe ')
  call pbe_libxc(dens, n, nspin, vxc, exc_pbe)
  deallocate(exc_pbe)
end subroutine vxc_pbe_new
