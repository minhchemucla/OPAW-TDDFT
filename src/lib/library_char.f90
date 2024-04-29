subroutine char3(i,ch)
  implicit none
  integer i
  character(len=1)c1
  character(len=2)c2
  character(len=3)ch
  select case(i)
  case(0:9)
     write(c1,"(I1)")i
     ch='00'//c1
  case(10:99)
     write(c2,"(I2)")i
     ch='0'//c2
  case(100:999)
     write(ch,"(I3)")i
  case default
     write(6,*)' problem, i in char3 =',i,' stopping '
     stop
  end select
end subroutine char3

subroutine char4(i,ch)
  implicit none
  integer i
  character(len=3)c3
  character(len=4)ch
  select case(i)
  case(0:999)
     call char3(i,c3)
     ch='0'//c3
  case(1000:9999)
     write(ch,"(I4)")i
  case default
     write(6,*)' problem, i in char4 =',i,' stopping '
     stop
  end select
end subroutine char4
