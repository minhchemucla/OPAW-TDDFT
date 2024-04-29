!program c_acos_driver
subroutine c_acos_driver
  implicit none
  complex*16  zin, cout, c_acos
  real*8 x,y
  write(6,*)' give me x, y '
  read( 5,*) x, y
  zin = dcmplx(x, y)
  cout = c_acos(zin)
  write(6,*)' cout      ',cout
  write(6,*)' zin       ',zin
  write(6,*)' cos(cout) ',cos(cout)
end 

function c_acos(zin)
implicit none
  complex*16 zin, cz, diff, c_acos, c_acosh, u
  complex*16, parameter :: ci=(0.d0,1.d0)
  integer i

  if(real(zin)<-1) then
     u = -zin
     cz = c_acosh(u) / ci + dacos(-1.d0)
     goto 99
  else if(real(zin)>1) then
     cz = c_acosh(zin) / ci
     goto 99
  end if

  cz = 1-zin
  do i=1,20
     diff = cos(cz)-zin
     cz = cz+1/sin(cz)*diff
  enddo

99 continue
  diff = cos(cz)-zin
  if(abs(diff)>1.d-10) then
     write(6,*)' diff , cz, zin ',diff, cz, zin
     stop
  endif

  c_acos = cz
end

function c_acosh(zin)
implicit none
  complex*16 zin, cz, diff, c_acosh
  complex*16, parameter :: ci=(0.d0,1.d0)
  integer i

  if(real(zin)<=1) then
     write(6,*)' problem, real(zin)<= 1' ,zin
     stop
  endif

  cz = zin-1
  do i=1,20
     diff = c_cosh(cz)-zin
     cz = cz-1/c_sinh(cz)*diff
  enddo

  diff = c_cosh(cz)-zin
  if(abs(diff)>1.d-10) then
     write(6,*)' diff , cz, zin ',diff, cz, zin
     stop
  endif

  c_acosh = cz

contains
  function c_cosh(cz)
    implicit none
    complex*16 cz, c_cosh
    c_cosh = (exp(cz)+exp(-cz))/2
  end function c_cosh
  
  function c_sinh(cz)
    implicit none
    complex*16 cz, c_sinh
    c_sinh = (exp(cz)-exp(-cz))/2
  end function c_sinh

  
end function c_acosh
