subroutine writeflush(j)
  use mpi_lib_ours, only : rank
  implicit none
  integer j
  if(rank==0) then
     write(6,*)' stage: ',j
     call flush(6)
  end if
end subroutine writeflush
