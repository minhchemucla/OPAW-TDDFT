subroutine r_mtrx_exp(r_A, r_exp_A, N) 
  implicit none
  integer N, j, matz, jp, i
  real*8  diff, delta ; external delta

  real*8 r_A(N,N), r_exp_A(N, N)

  real*8, dimension(:),   allocatable :: r_VL
  real*8, dimension(:,:), allocatable :: r_Vc, r_A_try, r_Vr_inv
  real*8    , dimension(:,:), allocatable :: unity

  integer isymm
  integer ifirst ; save ifirst; data ifirst /1 /

  ifirst = ifirst + 1
  matz = 1

  allocate(r_vc(N, N), r_vr_inv(N, N), r_VL(N), stat = j); if(j/=0) stop
  allocate(unity(N, N), r_A_try(N, N),          stat = j); if(j/=0) stop

  unity = 0.d0
  do i=1, N
     unity(i, i) = 1.d0
  enddo

  if(sum(abs(r_A - transpose(R_A)))<1.d-8) then
     isymm = 1
     call r_diag_normlz(r_A, r_Vc, r_VL, N)
     r_Vr_inv = transpose(r_Vc)
  else
     isymm = 0
!     call rg_diag(R_A, R_Vc, R_VL, N, matz)
!     call r_inv_mtrx(r_Vc, r_Vr_inv, N)   
      WRITE(6,*)' not ready for isymm ', isymm; stop
  endif
  
  diff =0
  do j=1, N
     diff = diff + sum(abs(matmul(r_A,r_Vc(:,j))-r_Vc(:,j)*r_vl(j)))
  enddo
  if(diff>1.d-8) then
     write(6,*)' diff (r_A*r_v-r_v*r_Vl) ', diff
     stop
  endif
  
  diff = sum(abs(matmul(r_vc,r_vr_inv)-unity))
  if(diff > 1.d-8) then 
     write(6,*)' overalp diff ',diff
     stop
  endif
  
  !
  ! check more  can erase
  ! 
  do j=1, N
     r_A_try(:, j) = r_vC(:, j)* r_Vl(j)
  enddo
  
  r_A_try = matmul( r_A_try, r_Vr_inv)
  diff = sum(abs(r_A_try - r_A))
  
  if(diff > 1.d-8) then
     write(6,*)' diff (r_A_try - R_A))',diff
     stop
  endif

  
  do j=1, N
     r_exp_A(:, j) = r_vC(:, j)* exp(r_Vl(j))
  enddo
  r_exp_A = matmul( r_exp_A, r_Vr_inv)
  
  deallocate(r_vc, r_vl, r_vr_inv, unity, r_A_try) 
  
end subroutine r_mtrx_exp


subroutine r_grahm_s_fast(ra, N, Nvec, Wght)  ! no checks, one loop only.
  implicit none                               ! when doubting, use r_grahm_s(
  integer N, Nvec, i, j
  real*8 ra(N, Nvec)
  real*8 Wght
  
  do i=1, Nvec
     do j=1, i-1
        ra(:, i) = ra(:, i) - sum(ra(:,i)*ra(:,j))*ra(:,j)*wght
     enddo
     ra(:, i) = ra(:, i)/ dsqrt(sum(ra(:,i)**2)*wght)
  enddo
end subroutine r_grahm_s_fast
  






