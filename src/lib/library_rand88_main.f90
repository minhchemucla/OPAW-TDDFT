module rand_fast_module
  implicit none
  save
  logical :: called_prep_rand = .false.
  integer :: seed(3), seed0(3)
end module rand_fast_module

subroutine prep_rand
  use rand_fast_module, only : called_prep_rand
  implicit none
  called_prep_rand = .true.
  call read_seed0_fast
end subroutine prep_rand

subroutine read_seed0_fast
  use rand_fast_module, only : seed, seed0
  use simple_mpi, only : rank, bcast_i
  implicit none
  integer :: k, i, i15
  integer, save :: i1=1
  integer*8 :: m(3), i31
  character*10 ch

  if(i1/=1) return; i1=-1

  i31 = 2**16; i31 = i31*(2**15)
  i15 = 2**15
  !write(6,*)' i31 ',i31

  call check_le(i15,i15*i15,' i15, i30 ')
  if(.not.(i15*i15<i31))stop ' i15, i31 problems '

  call check(size(m),size(seed0),' szm, szseed0 ')

  if(rank==0) then
     open(4,file='random.inp',status='old')! dont erase: need counter.inp too.
     !do k=1,6
     !   read(4,*)
     !enddo
     read(4,*)ch,m(:)
     if(ch/='taus88')then
        write(6,*)' ERROR, stopping. 7th line in fort.4 should start with taus88 '
     end if
     
     do i=1,size(m)
        m(i) = mod(m(i),i31)
        seed0(i) = m(i)
     enddo
     write(6,*)' seed0 ',seed0
     close(4)
  end if
  
  call bcast_i(seed0, size(seed0), 0)
  seed = seed0
  write(6,*)' seed ',seed
end subroutine read_seed0_fast
  
subroutine assign_seed(j)
  use rand_fast_module, only : seed0, seed
  use rand_fast_module, only : called_prep_rand
  use Ecuyer_random_modified, only : s1, s2, s3
  use simple_mpi, only : rank, bcast_i
  implicit none
  integer st, j
  integer, external :: is_it_power2
  integer*8 m(3), i31, k1, k2, k3, k4, j8

  if(.not.called_prep_rand)call prep_rand

  i31 = 2**16; i31 = i31*(2**15)  ! 2^31 safely
  call check(size(m),size(seed0),' szm, szseed0 ')

  s1 = seed0(1)
  s2 = seed0(2)
  s3 = seed0(3)
  
  j8=j
  k1=17287
  k2=38478
  k3=93043
  k4=45588
  
  m(1) = s1 + k1*j8            ; seed(1)=mod(m(1),i31)
  m(2) = s2 + k2*j8 + k3*j8*j8 ; seed(2)=mod(m(2),i31)
  m(3) = s3 +         k4*j8*j8 ; seed(3)=mod(m(3),i31)

  write(6,*)' seed post assign ',seed
end subroutine assign_seed

subroutine rand_zv(z, m, wg)
  use rand_fast_module, only : called_prep_rand
  implicit none
  integer    m, st, i
  integer, allocatable :: i_a(:)
  real*8     a, pi, wg
  complex*16 ci, z(m)
  if(.not.called_prep_rand) stop ' prep_rand not called ! '

  allocate(i_a(m), stat=st)
  call check0(st,' i_apt ')
  !write(6,*)' m, size(i_a) ',m, size(i_a)

  a = 2d0**(-32d0)
  pi = dacos(-1d0)
  ci = dcmplx(0d0,1d0)
  call i_rand_fast_array(i_a, size(i_a)) ! rand # -2^31 to 2^31
  !write(6,*)' i_a ',i_a
  do i=1,m
     z(i) = exp( ci * 2d0* pi * (dble(i_a(i))* a + 0.5d0))/sqrt(wg)  ! r\in [0,1]
  enddo
  !write(6,*)' i_a, wg ',minval(dble(i_a)),maxval(dble(i_a)), wg
  !write(6,*)' z   ',maxval(abs(z))
  deallocate(i_a)
end subroutine rand_zv

subroutine  i_rand_fast_array(i_a, n) ! rand # list, -2^31 to 2^31
  use Ecuyer_random_modified, only : i_rand_fast_array_inmodule , s1, s2, s3, init_seeds
  use rand_fast_module, only : seed
  implicit none
  integer  n, i_a(n)
  call init_seeds(seed(1), seed(2), seed(3)) ! maybe needed? not sure - DN
  !s1 = seed(1); s2 = seed(2); s3 = seed(3)
  call i_rand_fast_array_inmodule(i_a, n)
  seed = (/ s1, s2, s3 /)
end subroutine i_rand_fast_array
