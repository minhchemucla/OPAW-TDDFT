
subroutine hsub_cheby_dummy(ps,hps,nt)
  implicit none
  integer nt
  complex*16 ps(nt), hps(nt)
  write(6,*)' dummy -- should not get here '
  stop
end
