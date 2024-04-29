subroutine calc_nvirtual
  use main_mod
  use tddft_mod
  use mpi_lib_ours
  implicit none
  integer :: i, j


  if (nvirtual_read < 0) then
    if (rank==0) then
      write(6,*) 'calculating possible choices for nodes'
      write(6,*) 'nvirtual, nodes, states/nodes'
      do nvirtual=0,nstates-nocc
        do j=2,nocc
          if (mod(nocc+nvirtual,j) .eq. 0 .and. (nvirtual < (nocc+nvirtual)/j)) then
            write(6,*) nvirtual, j, (nocc+nvirtual)/j
          endif
        enddo
      enddo
    endif
    call sync_mpi
    call finalize_mpi_lib
    stop
  else
    if (rank==0) then
      if(mod(nocc+nvirtual_read,nodes) .ne. 0) stop  'mod(nocc+nvirtual_read,nodes) /= 0'
    endif
    nvirtual = nvirtual_read
  endif
end subroutine
