subroutine inv_driver()   ! dummy routine to check r_inv, c_inv
implicit none
integer, parameter :: n=3
integer i, j
real*8      a(n,n), a_inv(n,n), unity(n,n), diff, ran_normal
complex*16  c(n,n), c_inv(n,n);  complex*16, parameter :: ci=(0.d0, 1.d0)

unity=0.d0; do i=1,n; unity(i,i) = 1.d0; enddo

do i=1, n
   do j=1, n
      a(i,j) = ran_normal()
      a(j,i) = a(i, j)
      c(i,j) = a(i, j) + ci*ran_normal()
      c(j,i) = c(i, j)
   enddo
enddo

write(6,*)' a ';       write(6,88)a ; call r_inv_mtrx(A,A_inv,N)
write(6,*)' a_after '; write(6,88)a
write(6,*)' a_inv   '; write(6,88)a_inv

88 format(' ',3f20.8)

diff = sum(abs(matmul(a, a_inv)-unity)); write(6,*)' diff ',diff


write(6,*)' c ';       write(6,99)c ; call c_inv_mtrx(C,C_inv,N)
write(6,*)' c_after '; write(6,99)c
write(6,*)' c_inv   '; write(6,99)c_inv

99 format(' ',6f13.5)

diff = sum(abs(matmul(c, c_inv)-unity)); write(6,*)' diff ',diff

end



