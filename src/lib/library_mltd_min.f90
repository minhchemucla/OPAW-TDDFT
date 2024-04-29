! driver subroutine for min_mltd

subroutine min_mltd_driver
implicit none

integer, parameter :: nf=1, ny=2
real*8 vharm; external vharm
real*8 fixed(nf), y(ny), range

range = 7d0
fixed = 1.d0
y(:) =  0.2d0
call min_mltd(fixed, nf, y, ny, range, vharm)
write(6,*)' after ',y
end


! This f90 routine searches the minimum of a mutlidimensional function
!
! a very crude algortihm is employed

subroutine min_mltd(fixed, nf, y, ny, range, Vf)
  implicit none
  
  integer ny, nf, j, iter, isign, iy; external isign
  real*8  fixed(nf)
  real*8  y(ny)   ! initial y is used as a guess; final results reutned in y
  real*8  range, dy, Vf, vmid, vup, vdown, v_1, v_2, v_ratio
  
  real*8, allocatable, dimension(:) :: yup, ydown, vgrad_y
  
  allocate(yup(ny), ydown(ny), vgrad_y(ny),stat=j); if(j/=0)stop
  
  y=0.d0; dy = range *0.001
  
  do iter=1, 50
     
     vmid = Vf(fixed, nf, y, ny)
     
     vgrad_y = 0.d0
     do iy=1, ny
        yup = y  ;   yup(iy) = y(iy)+dy
        ydown = y; ydown(iy) = y(iy)-dy
        
        vgrad_y(iy) =   &
             (Vf(fixed, nf, yup, ny)-Vf(fixed,nf,ydown,ny))/  (2*dy)
        !! write(6,*)' yup, vf ',yup,Vf(fixed, nf, yup, ny)
        !! write(6,*)' ydn, vf ',yup,Vf(fixed, nf, ydown, ny)

     enddo

     if(sum(abs(vgrad_y))<1.d-10) goto 99

     vgrad_y = vgrad_y / dsqrt(sum(vgrad_y**2))   ! normalize gradient
     
     !
     !  now go until harmonic point reached.
     !

     yup = y + dy* vgrad_y;  vup   = Vf(fixed, nf, yup,   ny)
     ydown = y-dy* vgrad_y;  vdown = Vf(fixed, nf, ydown, ny)
     
     v_1 = (vup-vdown)/(2*dy)
     v_2 = (vup+vdown-2*vmid)/dy**2
     
     v_ratio=v_1/v_2;if(abs(v_ratio)>10*range)v_ratio=isign(v_ratio)*10*range

     y = y + vgrad_y *(-v_ratio)   
     
  enddo
99 continue

  vmid = Vf(fixed, nf, y, ny);  write(6,*)' vmid ', vmid

  deallocate(yup, ydown, vgrad_y)

end subroutine min_mltd


function isign(x)
  implicit none
  real*8 x
  integer isign
  if(abs(x)<1.d-15) then
     isign=0
  else if(x>0.d0) then
     isign=1
  else if(x<0.d0) then
     isign=-1
  else
     write(6,*)' x in isign ',x; stop
  endif
end function isign

function vharm(fixed, nf, y, ny)
  implicit none
  integer nf, ny
  real*8 fixed(nf), y(ny), vf, vharm
  if(nf/=1.or.ny/=2) stop
  
  !vf = 0.2*y(1)**2 + 0.2*y(2)**2
  vf = 0.3*(y(1)+2*y(2)-fixed(1))**2 + 0.2*(y(1)+y(2))**2
  vharm = vf
end function vharm
