
subroutine r_first_drv_5_pnt_lot( ain, a1,  N, lot, dy)
  implicit none

  integer N, lot
  real*8  dy
  real*8 ain(N, lot), a1(N, lot)
  integer j

  real*8 aaa
  real*8 bbb


  aaa =  2.d0/3.d0
  bbb = -1.d0/12.d0

  if(N<4) then
     write(6,*)' problem, N too small for 1drv 5 points ',N
     stop
  endif

  a1 = 0.d0
  do j=1, N-1
     a1(j,:) =  a1(j,:) + aaa*ain(j+1, :)
  enddo

  do j=1+1, N
     a1(j,:) =  a1(j,:) - aaa*ain(j-1, :)
  enddo

  do j=1, N-2
     a1(j,:) =  a1(j,:) + bbb* ain(j+2, :) 
  enddo

  do j=1+2, N
     a1(j,:) =  a1(j,:) - bbb* ain(j-2, :) 
  enddo

  a1 = a1/dy
end

subroutine c_first_drv_5_pnt_lot( ain, a1,  N, lot, dy)
  implicit none

  integer N, lot
  real*8  dy
  complex*16 ain(N, lot), a1(N, lot)
  integer j

  real*8 aaa
  real*8 bbb

  aaa =  2.d0/3.d0
  bbb = -1.d0/12.d0

  if(N<4) then
     write(6,*)' problem, N too small for 1drv 5 points ',N
     stop
  endif

  a1 = 0.d0
  do j=1, N-1
     a1(j,:) =  a1(j,:) + aaa*ain(j+1, :)
  enddo

  do j=1+1, N
     a1(j,:) =  a1(j,:) - aaa*ain(j-1, :)
  enddo

  do j=1, N-2
     a1(j,:) =  a1(j,:) + bbb* ain(j+2, :) 
  enddo

  do j=1+2, N
     a1(j,:) =  a1(j,:) - bbb* ain(j-2, :) 
  enddo

  a1 = a1/dy

end


subroutine r_second_drv_5_pnt_lot( ain, a1,  N, lot, dy)
  implicit none

  integer N, lot
  real*8  dy
  real*8  ain(N, lot), a1(N, lot)
  integer j

  real*8 aaa
  real*8 bbb
  real*8 ccc

  aaa =  4.d0/3.d0
  bbb = -1.d0/12.d0
  ccc = -5.d0/2.d0

  if(N<4) then
     write(6,*)' problem, N too small for 1drv 5 points ',N
     stop
  endif

  a1 = ccc*ain

  do j=1, N-1
     a1(j,:) =  a1(j,:) + aaa*ain(j+1, :)
  enddo

  do j=1+1, N
     a1(j,:) =  a1(j,:) + aaa*ain(j-1, :)
  enddo

  do j=1, N-2
     a1(j,:) =  a1(j,:) + bbb* ain(j+2, :) 
  enddo

  do j=1+2, N
     a1(j,:) =  a1(j,:) + bbb* ain(j-2, :) 
  enddo

  a1 = a1/dy**2
end

subroutine c_second_drv_5_pnt_lot( ain, a1,  N, lot, dy)
  implicit none

  integer N, lot
  real*8  dy
  complex*16  ain(N, lot), a1(N, lot)
  integer j

  real*8 aaa
  real*8 bbb
  real*8 ccc

  aaa =  4.d0/3.d0
  bbb = -1.d0/12.d0
  ccc = -5.d0/2.d0

  if(N<4) then
     write(6,*)' problem, N too small for 1drv 5 points ',N
     stop
  endif

  a1 = ccc*ain

  do j=1, N-1
     a1(j,:) =  a1(j,:) + aaa*ain(j+1, :)
  enddo

  do j=1+1, N
     a1(j,:) =  a1(j,:) + aaa*ain(j-1, :)
  enddo

  do j=1, N-2
     a1(j,:) =  a1(j,:) + bbb* ain(j+2, :) 
  enddo

  do j=1+2, N
     a1(j,:) =  a1(j,:) + bbb* ain(j-2, :) 
  enddo

  a1 = a1/dy**2
end



!program    drive_r_2_drv
subroutine  drive_r_2_drv()

  implicit none
  integer, parameter :: N=10
  integer, parameter :: lot=6
  real*8,  parameter ::  dy = 0.24
  integer j, ilot
  real*8   ain(N, lot)
  real*8   a1(N, lot)
  real*8   xg(N)

  do j=1, N
     xg(j) = 0.3 + j*dy
  enddo
  
  ain(:, 1) = 1.5
  do ilot=1+1, lot
     ain(:, ilot) = xg**(ilot-1)
  enddo
  
  call r_second_drv_5_pnt_lot( ain, a1,  N, lot, dy)

  do ilot=1,lot
     write(8,*)' ilot, ',ilot, ' a1, analytical '
     do j=1, N
        write(8,*)  a1(j,ilot), ain(j,ilot)/xg(j)**2*(ilot-1)*(ilot-2)
     enddo
     write(8,*)
  enddo
  
end 



