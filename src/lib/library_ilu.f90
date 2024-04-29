
!
! ilu, for real matrices, converts a real matrix A to its LU form. Note that A is modified; L is the below diagonal of the output with 1 on the diagonal,
! and U is the diagonal and above.  See split_LU below
!
subroutine icholesky_r(A,L,P,n)  !  for symmetric matrices only
  implicit none
  integer n,i,j,k
  real*8  A(n,n), L(n,n), P(n,n),e,d
  
  write(6,*)' minval(p) ',minval(p),maxval(p)

  L = 0d0
  do i=1,n
     d= a(i,i)
     do k=1,i-1
        if(P(i,k)>0) d = d - L(i,k)**2
     end do
     if(d<1d-10) then
        write(6,*)' incomplete cheloesky failed, i, d ',i,d
        L(i,i)=1d0
     else
        L(i,i) = sqrt(d)
     end if
     do j=i+1,n
        if(P(j,i)>0d0) then
           L(j,i) = A(i,j)
           do k=1,i-1
              if(P(i,k)>0d0) L(j,i) = L(j,i)-L(i,k)*L(j,k)
           enddo
           L(j,i)=L(j,i)/L(i,i)
        end if
     end do
  end do

end subroutine icholesky_r


subroutine ilu_r(A,P,n)
  implicit none
  integer n,r,i,j
  real*8  A(n,n), P(n,n),e,d

  do r=1,n-1
     d= 1d0/A(r,r)
     do i=r+1,n
        if(P(i,r)>0) then
           e = d * A(i,r)
           A(i,r) = e
           do j=r+1,n
              if(P(i,j)>0.and.P(r,j)>0) A(i,j)=A(i,j)-e*A(r,j)
           enddo
        end if
     end do
  end do
end subroutine ilu_r


! given A, toll and n, gives L U and P where P is 1 where A,L,U are non-zero
subroutine split_LU_r(A,toll,n,L,U,P)
  implicit none
  integer i,n
  real*8  toll
  real*8  A(n,n), L(n,n), U(n,n), P(n,n)

  P = 0
  where (abs(A)>toll) P=1
  
  U=A
  L= 0d0
  call ilu_r(U,P,n)
  do i=1,n
     L(i,i)=1d0
     L(i,1:i-1)=U(i,1:i-1)
     U(i,1:i-1)=0d0
  end do
end subroutine split_LU_r

subroutine split_cholesky_r(A,toll,n,L,P)
  implicit none
  integer n
  real*8  toll
  real*8  A(n,n), L(n,n), P(n,n)

  P = 0
  where (abs(A)>toll) P=1
  
  call icholesky_r(A,L,P,n)
end subroutine split_cholesky_r
