module cheby_coef_work
  implicit none; save
  integer                 :: N,m
  real*8,     allocatable :: d(:)
  complex*16, allocatable :: f(:)
end module cheby_coef_work

subroutine cheb_coeff_r_theta_mu(nnc, havg, DL, mu, Tp, ro)
! chebyshev coeff. for f(x)=theta_smooth(mu-x)=sum_k=0,1,2...nnc  co(k) T_n( (x-havg)/DL), for -DL+havg<x<DL+havg ! i think the havg is correct.

  use cheby_coef_work, only : n, m, f, d
  implicit none
  integer nnc,i
  real*8  havg, DL,  mu, Tp, ro(0:nnc) 
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  do i=0,size(f)-1; f(i) = erfc((D(i)-mu)/Tp)/2.d0;  enddo
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_theta_mu

subroutine cheb_coeff_r_theta_theta_mu(nnc, havg, DL, mu, Tp, ro)
  ! chebyshev coeff. for f(x)=(theta_smooth(mu-x))^2=sum_k=0,1,2...nnc  co(k) T_n( x/DL), for -DL<x<DL.
  use cheby_coef_work, only : n, m, f, d
  implicit none
  integer nnc,i
  real*8  havg, DL,  mu, Tp, ro(0:nnc) 
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  do i=0,size(f)-1; f(i) = erfc((D(i)-mu)/Tp)/2.d0;  enddo
  f = f*f  ! NOTE!
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_theta_theta_mu

subroutine cheb_coeff_r_theta_H_theta_mu(nnc, havg, DL, mu, Tp, ro)
  ! chebyshev coeff. for f(x)=(theta_smooth(mu-x))^2=sum_k=0,1,2...nnc  co(k) T_n( x/DL), for -DL<x<DL.
  use cheby_coef_work, only : n, m, f, d
  implicit none
  integer nnc,i
  real*8  havg, DL,  mu, Tp, ro(0:nnc) 
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  do i=0,size(f)-1; f(i) = erfc((D(i)-mu)/Tp)/2.d0;  enddo
  f = f*f*D  ! NOTE!
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_theta_H_theta_mu


subroutine cheb_coeff_r_theta_Ebottop(nnc, havg, DL, ebot, etop, Tp, ro)
  use cheby_coef_work, only : n, m, f, d
  implicit none
  integer nnc,i
  real*8  havg, DL,  ebot, etop, Tp, ro(0:nnc) 
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  do i=0,size(f)-1; f(i) = erfc((D(i)-etop)/Tp)*erfc((ebot-D(i))/Tp)/4.d0; enddo
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_theta_Ebottop

subroutine cheb_coeff_r_gaussian_E(nnc, havg, DL, E, Ewdth, ro)
  use cheby_coef_work, only : n, m,f, d
  implicit none
  integer nnc
  real*8  havg, DL,  E, Ewdth,  ro(0:nnc) 
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  f = exp(-((E-D)**2)/Ewdth**2/2.d0)
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_gaussian_E

subroutine cheb_coeff_r_gaussian_x_E(nnc, havg, DL, E, Ewdth, ro)
  use cheby_coef_work, only : n, m,f, d
  implicit none
  integer nnc
  real*8  havg, DL,  E, Ewdth,  ro(0:nnc) 
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  f = D*exp(-((E-D)**2)/Ewdth**2/2.d0)
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_gaussian_x_E

subroutine cheb_coeff_r_theta_Ebottop_gaussian_E(nnc, havg, DL, Ebot, Etop, Tp, E, Ewdth, ro)
  use cheby_coef_work, only : n, m,f, d
  implicit none
  integer nnc,i
  real*8  havg, DL,  Ebot, Etop, Tp, E, Ewdth,  ro(0:nnc)
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  do i=1, size(f); f(i) = exp(-(E-D(i))**2/Ewdth**2/2.d0)* erfc((D(i)-etop)/Tp)*erfc((ebot-D(i))/Tp)/4.d0;enddo
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_theta_Ebottop_gaussian_E

subroutine cheb_coeff_r_inv_sqrt(nnc, havg, DL, Tp, ro)
  use cheby_coef_work, only : n, m,f, d
  implicit none
  integer nnc
  real*8  havg, DL, Tp,ro(0:nnc)
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  where(D>0)
     f = 1/sqrt(D)
  elsewhere
     f=0d0
  end where
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_inv_sqrt


subroutine cheb_coeff_r_sqrt(nnc, havg, DL, Tp, ro)
  use cheby_coef_work, only : n, m,f, d
  implicit none
  integer nnc
  real*8  havg, DL, Tp,ro(0:nnc)
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  where(D>0)
     f = sqrt(D)
  elsewhere
     f=0d0
  end where
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_sqrt

subroutine cheb_coeff_sqrt(nnc, havg, DL, ro) ! like cheb_coeff_r_sqrt, without the unnencessary Tp
  use cheby_coef_work, only : n, m,f, d
  implicit none
  integer nnc
  real*8  havg, DL, ro(0:nnc)
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  where(D>0)
     f = sqrt(D)
  elsewhere
     f=0d0
  end where
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_sqrt

subroutine cheb_coeff_r_inv_power(nnc, havg, DL, gamma, power, ro)
  use cheby_coef_work, only : n, m,f, d
  implicit none
  integer nnc,i,power
  real*8  havg, DL, Tp,ro(0:nnc),gamma
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  f = 1d0/(gamma**power+D**power)**(1d0/dble(power))
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  call check_real(f,size(f)); ro   = f(0:nnc)
end subroutine cheb_coeff_r_inv_power


subroutine cheb_coeff_c_Phi(nnc, havg,dL, mu, Tp, t, co)
  use cheby_coef_work, only : n, m, f, d
  implicit none
  integer nnc, i
  real*8  havg, DL, mu, Tp, t
  complex*16              :: co(0:nnc) 
  complex*16, parameter   :: ci=(0d0,1d0)
  call prep_for_cheby_n_f_d(nnc, havg, dl); 
  if(    abs(t)<1d-8) then;  do i=0,size(f)-1; f(i) = exp(-ci*t*D(i)) * (0.5d0-erfc((D(i)-mu)/Tp)/2.d0 ); enddo
  elseif(    t< 0d0 ) then;  do i=0,size(f)-1; f(i) = exp(-ci*t*D(i)) * (     -erfc((D(i)-mu)/Tp)/2.d0 ); enddo
  else;                      do i=0,size(f)-1; f(i) = exp(-ci*t*D(i)) * (1d0  -erfc((D(i)-mu)/Tp)/2.d0 ); enddo;endif
     
  call fftsa(N,f,m);   f = 2.d0*f/dble(N);  f(0) = f(0)/2.d0;  co   = f(0:nnc)
  
  call possibly_print_coef

contains
  
  subroutine possibly_print_coef
    use mpi_lib_ours, only : rank
    implicit none
    integer, save :: cn=0
    if(cn<3) then
       if(rank==0) then
          write(6,*)' writing to file ',700000+cn; call flush(6)
          do i=0,nnc
             write(70000+cn,*)i,dble(f(i)),aimag(f(i))
          enddo
          call flush(70000+cn)
       end if
       cn=cn+1
    end if
    !call time_print(' fnsh poss prnt coef ')
  end subroutine possibly_print_coef
end subroutine cheb_coeff_c_phi

subroutine prep_for_cheby_n_f_d(nnc, havg, dl)
  use cheby_coef_work, only : n,m,f, d
  implicit none
  integer i, nnc
  real*8 pi,w, havg, dl
  pi = dacos(-1.d0)
  m = 3 + nint( log(dble(nnc))/log(2.d0));   N = 2**m;  
  if(allocated(f))  deallocate(f)
  if(allocated(d))  deallocate(d)
  allocate(f(0:n-1),d(0:n-1),stat=i); call check0(i,' prchbw ')
  do i=0,N-1; w = 2d0*pi/dble(N) * i; D(i) = DL * cos(w)+havg;   enddo
end subroutine prep_for_cheby_n_f_d
