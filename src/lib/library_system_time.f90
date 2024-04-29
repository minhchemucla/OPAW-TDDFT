subroutine system_time(ts)
  implicit none
  real*8 ts
  integer(kind=8) count, count_rate
  call system_clock(count, count_rate)
  ts= dble(count)/dble(count_rate)
end subroutine system_time

subroutine system_time_stamp(a)
  implicit none
  character(*) a
  integer, save :: i1=1
  real*8,  save :: t0, t1

  if(i1==1) then
     call system_time(t0)
     i1=-1
  end if

  call system_time(t1)
  write(6,*)' system time for ',a,' was ',t1-t0
  t0=t1
end subroutine system_time_stamp
