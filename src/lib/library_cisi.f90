!program si_prog
!  implicit none
!  integer i
!  real*8 x, ci, sia
!  real*8, external :: si
!  do i=-1,5
!     if(i==-1) then; x=0d0
!     else;     x=2**i
!     endif
!     sia = si(x)
!     write(6,*) x, sia, ' x, si '
!  enddo
!end program si_prog

real*8 function si(x) ! sine integral
  implicit none
  real*8 x, ci, sia
  call cisi(x,ci,sia)
  si = sia
end function si

SUBROUTINE cisi(x,ci,si) ! from Numerical Recipes, tiny modification for double precision
  implicit none
  REAL*8 ci,si,x
  real*8, PARAMETER:: EPS=1d-24
  real*8, parameter :: EULER=0.57721566490153286d0
  integer, parameter :: MAXIT=200
  real*8, parameter :: fpmin = 1d-40, tmin=2d0
  real*8 pi, piby2
  INTEGER i,k  
  REAL*8 a,err,fact,sign,sum,sumc,sums,t,term
  COMPLEX*16 h,b,c,d,del  
  LOGICAL odd  

  pi = dacos(-1d0)
  piby2 = pi/2d0
  
  !write(6,*)' euler ', euler


  t=abs(x)  
  if(t.eq.0.d0)then  
     si=0.d0
     ci=-1.d0/FPMIN  
     return  
  endif
  if(t.gt.TMIN)then  
     b=dcmplx(1.d0,t)  
        c=1.d0/FPMIN  
        d=1.d0/b  
        h=d  
        do 11 i=2,MAXIT  
          a=-((dble(i-1))**2)
          b=b+2.d0
          d=1.d0/(a*d+b)  
          c=b+a/c  
          del=c*d  
          h=h*del  
          if(absc(del-1.d0).lt.EPS)goto 1  
11      continue  
        pause 'cf failed in cisi'  
1       continue  
        h=dcmplx(cos(t),-sin(t))*h  
        ci=-dble(h)  
        si=PIBY2+aimag(h)  
      else  
        if(t.lt.sqrt(FPMIN))then  
          sumc=0.d0  
          sums=t  
        else  
          sum=0.d0  
          sums=0.d0  
          sumc=0.d0  
          sign=1.d0  
          fact=1.d0  
          odd=.true.  
          do 12 k=1,MAXIT  
            fact=fact*t/k  
            term=fact/k  
            sum=sum+sign*term  
            err=term/abs(sum)  
            if(odd)then  
              sign=-sign  
              sums=sum  
              sum=sumc  
            else  
              sumc=sum  
              sum=sums  
            endif  
            if(err.lt.EPS)goto 2  
            odd=.not.odd  
12        continue  
          pause 'maxits exceeded in cisi'  
        endif  
2       si=sums  
        ci=sumc+log(t)+EULER  
      endif  
      if(x.lt.0.d0)si=-si  
      return  
      contains
        real*8 function absc(h)
          complex*16 h
          absc=abs(dble(h))+abs(aimag(h))  
        end function absc
     END
