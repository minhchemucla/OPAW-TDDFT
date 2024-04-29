subroutine writef(a)
  use mpi_lib_ours, only : rank
  implicit none
  integer i
  character(*) a
  if(rank==0) then 
     i=6
  else
     i=86000+rank
  endif
  write(i,*)a
  call flush(i)
end subroutine writef
