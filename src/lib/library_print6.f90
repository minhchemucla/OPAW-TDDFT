subroutine print6(a)
  use mpi_lib_ours, only : rank
  implicit none
  character(*) a

  integer ip
  real*8, save :: tstart, t0, t1
  logical, save :: first=.true.

  if(rank/=0) ip = 20000+rank
  if(rank==0) ip = 6

  if(first) then;
     first = .false.
     call cpu_time(tstart)
     t0 = tstart
  endif

  call cpu_time(t1)
  write(ip,*)a,real(t1-t0),' acc. ',real(t1-tstart); call flush(ip)
  t0=t1

end subroutine print6


subroutine prnt(j,a)
  use mpi_lib_ours, only : rank,nodes
  implicit none
  character(*) a
  integer j,ip
  real*8,  save :: tstart, t0, t1
  logical, save :: first=.true.
  integer, save :: cnt(20)

  if(rank/=0.and.rank.ne.nodes/4.and.rank.ne.nodes-1) return
  
  if(rank.ge.10000)stop ' rank bigger than 10^4 '
  if(j>20.or.j<1) stop ' prnt: 1st arguments should be 1-20 '
  ip = 20000+rank
  if(j==1.and.rank==0)ip=6

  if(first) then;
     cnt(1:20) = 0
     first = .false.
     call cpu_time(tstart)
     t0 = tstart
  endif

  if(cnt(j)>100.and.j/=1) return

  cnt(j) = cnt(j)+1
  call cpu_time(t1)
  write(ip,*)'stage: ',j,' ',a,real(t1-t0),' acc. ',real(t1-tstart); 
  call flush(ip)
  t0=t1

end subroutine prnt
