     
!------
! this subroutine uses an iterative method to calculate
!  c_v_out = 1/(E-H)*c_v_in, using an c_hpsi routine (with name not passed)
!------

subroutine iterative_driver()
  implicit none
  integer, parameter :: n=3
  integer, parameter :: n_iter_max= 50000, n_en = 2
  real*8  energy(n_en); data energy / 2.d0, 0.1d0 /; save energy
  real*8 :: emin = -1,  emax=8
  real*8 :: vi_max = 0.2
  complex*16 c_v_in(n), c_v_out(n, n_en)
  complex*16,parameter :: ci= (0.d0, 1.d0)
  real*8 ran_ps_neg, evrangefactor
  integer i, ie
  external c_Hv_simple

  evrangefactor = 1.2
  do i=1, n
     c_v_in(i) = ran_ps_neg()
     c_v_in(i) = c_v_in(i) + ci*ran_ps_neg()
  enddo
  write(6,*)' c_v_in ',c_v_in

  call c_iter_inverse(c_v_in, c_v_out, n, emin, emax, vi_max, &
                          evrangefactor, &
                          n_iter_max, energy, n_en, c_Hv_simple)

  do ie=1,n_en
     write(96,*)' ie, c_v_out ',ie,c_v_out(:,ie)
  enddo

end    !  subroutine iterative_driver

!------------------------
! generalized and useful subroutine 
! to calculate 1/(energy-h)*c_v_in --> c_v_out.
! input: c_v_in(n), n, n_en, nergy(n_en),emin, emax, vi_max, n_itermax (play).
! uses a COMPLEX initial and final variables.
!------------------------
subroutine c_iter_inverse(c_v_in, c_v_out, n, emin, emax, vi_max, &
     EvRangeFactor,n_iter_max, energy, n_en, c_Hv)
  implicit none
  integer n_en, n, status, n_iter_max, j, ie, iter, is_it_power2
  real*8 emax, emin, energy(n_en), vi_max
  real*8      diff, diff_max
  real*8        EvRangeFactor

! note the complex on c_havg, dh, th, cs/snth

  complex*16, parameter :: ci=(0.d0,1.d0)
  complex*16 c_v_in(n), c_v_out(n, n_en), C_havg, C_dh
  complex*16 MaxEc, MinEc , Max_af, Min_af

  complex*16, dimension(:), allocatable :: c_v_0, c_v_1, c_v_2
  complex*16, dimension(:), allocatable :: csth, snth, expmith


  external c_Hv

  allocate (csth(n_en), snth(n_en), expmith(n_en), stat=j); if(j/=0)stop
  allocate (c_v_0(n), c_v_1(n), c_v_2(n),          stat=j); if(j/=0)stop


  write(6,*)' emin , emax ', emin, emax

  
  MaxEc = Emax + 0 
  MinEc = Emin - ci* Vi_max

  write(6,*)' MaxEc , MaxEc ', MinEc, MaxEc

  Max_af = (MaxEc*(1+EvRangeFactor)+MinEc*(1-EvRangeFactor))/2
  Min_af = (MinEc*(1+EvRangeFactor)+MaxEc*(1-EvRangeFactor))/2

  write(6,*)' Max_af , Max_af ', Min_af, Max_af

  c_havg = (Max_af + Min_af) / 2.0
  c_dh   = abs(Max_af - Min_af) / 2.0 
  
  write(6,*)' c_dh , c_havg ', c_dh, c_havg
 
  csth = (energy(:)-c_havg )/c_dh
  snth =   sqrt(1.d0-csth*csth)
  where(aimag(csth)<0.0d0)snth=-snth
  expmith = csth - ci*snth

  if(sum(abs(snth*snth+csth*csth-1.d0))>1.d-12) then 
     write(6,*)' problem '
     do ie=1, n_en
        write(6,*)' ie, snth, csth ', ie, snth(ie), csth(ie)
        write(6,*)' csth*csth+snth*snth ',csth*csth+snth*snth
     enddo
     stop
  endif
     
  do ie=1, n_en
     c_v_out(:,ie)= -ci/snth(ie)/c_dh *c_v_in(:)
  enddo

  c_v_0 = c_v_in(:)
  call c_scaled_ham_vec(c_v_0, c_v_1, n, c_havg, c_dh, c_Hv)

  do ie=1, n_en
     c_v_out(:, ie)= &
     c_v_out(:, ie) -2.d0*ci/snth(ie)/c_dh *expmith(ie)*c_v_1(:)
  enddo
  !
  !cheby
  !
  do iter=2, n_iter_max

     if(is_it_power2(iter) == 1) write(6,*)' iter ',iter

     call c_scaled_ham_vec(c_v_1, c_v_2, n, c_havg, c_dh, c_Hv)
     c_v_2 = c_v_2 *2 - c_v_0
     
     do ie=1, n_en
        c_v_out(:,ie) = &
        c_v_out(:,ie) - 2.d0*ci/snth(ie)/c_dh*((expmith(ie))**iter)*c_v_2(:)
     enddo

     c_v_0 = c_v_1; c_v_1 = c_v_2
  enddo

  !
  !  check
  !
  diff_max = 0.d0
  do ie=1, n_en
     c_v_1 = c_v_out(:, ie)
     call c_Hv(    c_v_1, c_v_2, n)  
     c_v_2 = energy(ie)* c_v_out(:, ie) -c_v_2
     write(96,*)' ie, (e-h)*result ', ie, c_v_2

     diff = sum(abs(c_v_2 - c_v_in))/sum(abs(c_v_in))
     write(96,*)' ie, diff ',ie, diff
     diff_max = max(diff_max, diff)
  enddo
  write(6,*)' diff_max ',diff_max

  deallocate (c_v_0, c_v_1, c_v_2, csth, snth, expmith)

end subroutine c_iter_inverse

subroutine c_scaled_ham_vec(c_v_0, c_v_1, n, c_havg, c_dh, c_Hv)
  implicit none
  integer n
  complex*16 c_dh
  complex*16 c_v_0(n), c_v_1(n)
  complex*16 c_havg  ! note !
  external c_Hv
  call c_Hv(c_v_0, c_v_1, n)
  c_v_1 = (c_v_1 - c_havg*c_v_0)/c_dh
end subroutine c_scaled_ham_vec
 
subroutine c_Hv_simple(c_v_0, c_v_1, n)
  implicit none
  integer n, i, j
  complex*16 c_v_0(n), c_v_1(n)
  complex*16, parameter :: ci = (0.d0, 1.d0)

!!  if(n /=1) then
!!     write(6,*)' problem, in n= ',n
!!     stop
!!  endif
!!  c_v_1 =  ( 2 - ci*0.1) *c_v_0

  complex*16 ch(3,3)
  save ch

  integer ifirst
  data ifirst / 1 /
  save ifirst

  if(n/=3) then; write(6,*)' n ',n;stop;endif

  ch = 0.d0
  do i=1, n
     do j=1, n
         ch(i,j)= i**2/50.d0+j**2/50.d0 
     enddo
     ch(i,i) =  ch(i, i) +  i*2.d0/3.d0 - ci*0.1d0/3.d0
  enddo

  if(ifirst < 2) then  
     write(6,*)' ch '
     do i=1,3
        write(6,*)(ch(i,j),j=1,3)
     enddo
     write(6,*)'  '
     ifirst = ifirst + 1
  endif

  if(sum(abs(transpose(ch)-ch))>1.d-8) stop

  c_v_1 = matmul(ch, c_v_0)

end subroutine c_Hv_simple











