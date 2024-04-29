subroutine fft_prop_nsurface_1d(cfull, n, nstate, &
     mass, spacing, dt, c_exp_dtov2_pot)

  implicit none
  real*8   mass, dt, spacing
  integer  n,    nstate, ix, j
  complex*16        cfull(n, nstate)
  complex*16  c_exp_dtov2_pot(n, nstate, nstate)  ! no imag. pot. yet.

  do ix=1, n
     cfull(ix, :) = matmul(c_exp_dtov2_pot(ix, :, :), cfull(ix, :))
  enddo

  do j=1, nstate
  call fft_prop_1d_good(cfull(:, j), &
                   n, mass, spacing, dt) 
  enddo
  
  do ix=1, n
     cfull(ix, :) = matmul(c_exp_dtov2_pot(ix, :, :), cfull(ix, :))
  enddo

end subroutine fft_prop_nsurface_1d
  
