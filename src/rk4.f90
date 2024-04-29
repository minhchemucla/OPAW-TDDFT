subroutine rk4_prop_shs(p,pp,n,nb)
  !time independent H for now
  ! p = unperturbed psi
  ! pp = perturbed psi
  use main_mod,  only : dv,nocc,nx,ny,nz
  use tddft_mod, only : dt, state_map
  use tddft_mod, only : phi_bar, phi_bar_pert
  use mpi_lib_ours
  implicit none
  integer :: is, jb, n, nb
  real*8 ::  nrm0(nb), nrm1(nb)
  complex*16, dimension(n,nb) :: k1,k2,k4,y, p, pp, p2,pp2
  complex*16, parameter :: ci = (0d0,1d0)
  !real*8 :: wtime, wtime2
  p2 = p
  call norm0(p2)
  call rk4_step_shs(p2,k1,y,dt/2d0,.false.)
  call rk4_step_shs(y,k2,y,dt/2d0,.false.)
  call rk4_step_shs(y,k4,y,dt,.false.) !k4=k3
  k2 = k2 + k4 !k2+k3 
  call rk4_step_shs(y,k4,y,dt,.true.)
  p = p2 + dt/6d0*(k1+2d0*k2+k4)
  call norm1(p)

  pp2 = pp
  call norm0(pp2)
  call rk4_step_shs_pert(pp2,k1,y,dt/2d0,.false.)
  call rk4_step_shs_pert(y ,k2,y,dt/2d0,.false.)
  call rk4_step_shs_pert(y ,k4,y,dt,.false.)
  k2 = k2 + k4 !k2+k3 
  call rk4_step_shs_pert(y ,k4,y,dt,.true.)
  pp = pp2 + dt/6d0*(k1+2d0*k2+k4)
  call norm1(pp)
  contains 
    subroutine rk4_step_shs(pin,k,y,dt,last)
      implicit none
      real*8     :: dt
      logical    :: last
      complex*16 :: k(n,nb), y(n,nb), pin(n,nb)

      do is=1,nb
        jb=state_map(is)
        if(jb<=nocc) then
          call shs(pin(:,is),k(:,is),1)
          k(:,is) = -ci*k(:,is)
          if(.not. last) then
            y(:,is) = p2(:,is) + dt * k(:,is) !y1
            phi_bar(:,:,:,is,1) = reshape(y(:,is), (/nx,ny,nz/)) 
          endif
        endif
      enddo

      if(.not.last) then
        call ncpaw_make_hamiltonian_tddft_unpert
      endif
    end subroutine

    subroutine rk4_step_shs_pert(pin,k,y,dt,last)
      implicit none
      real*8     :: dt
      logical    :: last
      complex*16 :: k(n,nb), y(n,nb), pin(n,nb)

      do is=1,nb
        jb=state_map(is)
        if(jb<=nocc) then
          call shs_pert(pin(:,is),k(:,is),1)
          k(:,is) = -ci*k(:,is)
          if(.not. last) then
            y(:,is) = pp2(:,is) + dt * k(:,is) !y1
            phi_bar_pert(:,:,:,is,1) = reshape(y(:,is), (/nx,ny,nz/)) 
          endif
        endif
      enddo

      if(.not.last) then
        call ncpaw_make_hamiltonian_tddft_pert
      endif
    end subroutine

    subroutine norm0(p)
      implicit none
      complex*16 :: p(n,nb)
      do is=1,nb
        jb=state_map(is)
        if(jb<=nocc) then
          nrm0(is) = sqrt(dv*sum(abs(p(:,is))**2))
        endif
      enddo
    end subroutine

    subroutine norm1(p)
      implicit none
      complex*16 :: p(n,nb)
      do is=1,nb
        jb=state_map(is)
        if(jb<=nocc) then
          nrm1(is) = sqrt(dv*sum(abs(p(:,is))**2))
          p(:,is) = p(:,is) * nrm0(is)/nrm1(is)
        endif
      enddo
    end subroutine
end subroutine

subroutine rk4_prop_sh
  !time independent H for now
  ! p = unperturbed psi
  ! pp = perturbed psi
  use main_mod
  use tddft_mod 
  use mpi_lib_ours
  implicit none
  real*8  nrm0, nrm1
  complex*16 :: p(nx,ny,nz), pp(nx,ny,nz)
  complex*16, dimension(nx,ny,nz) :: k1,k2,k4,y
  complex*16, parameter :: ci = (0d0,1d0)
  real*8  :: wtime, wtime2
  

  !wtime = mpi_wtime()
  !wtime2 = mpi_wtime()
  nrm0 = sqrt(dv*sum(abs(p)**2))
  call sh(p,k1,1)
  k1 = -ci*k1
  y = p + dt/2d0 * k1 !y1
  !wtime = mpi_wtime() - wtime
  !write(6,*) 'rank, wtime one sh', rank, wtime

  call sh(y,k2,1)
  k2 = -ci*k2
  y = p + dt/2d0 * k2 !y2

  call sh(y,k4,1) !k4=k3
  k4 = -ci*k4
  y = p + dt * k4 !y3

  k2 = k2 + k4 !k2+k3 

  call sh(y,k4,1)
  k4 = -ci*k4

  p = p + dt/6d0*(k1+2d0*k2+k4)

  nrm1 = sqrt(dv*sum(abs(p)**2))
  p = p * nrm0/nrm1

  nrm0 = sqrt(dv*sum(abs(pp)**2))
  call sh_pert(pp,k1,1)
  k1 = -ci*k1
  y = pp + dt/2d0 * k1 !y1

  call sh_pert(y,k2,1)
  k2 = -ci*k2
  y = pp + dt/2d0 * k2 !y2

  call sh_pert(y,k4,1) !k4=k3
  k4 = -ci*k4
  y = pp + dt * k4 !y3

  k2 = k2 + k4 !k2+k3 

  call sh_pert(y,k4,1)
  k4 = -ci*k4

  pp = pp + dt/6d0*(k1+2d0*k2+k4)
  nrm1 = sqrt(dv*sum(abs(pp)**2))
  pp = pp * nrm0/nrm1

  !wtime2 = mpi_wtime() - wtime2
  !write(6,*) 'rank, wtime 8 sh', rank, wtime2
end subroutine
!subroutine rk4_prop_shs_debug(is, p)
!  !time independent H for now
!  ! p = unperturbed psi
!  ! pp = perturbed psi
!  use main_mod
!  use tddft_mod 
!  use mpi_lib_ours
!  implicit none
!  integer :: is
!  real*8  nrm0, nrm1
!  complex*16 :: p(nx,ny,nz)
!  complex*16, dimension(nx,ny,nz) :: k1,k2,k4,y
!  complex*16, parameter :: ci = (0d0,1d0)
!  !real*8 :: wtime, wtime2
!
!  !wtime = mpi_wtime()
!  !wtime2 = mpi_wtime()
!  nrm0 = sqrt(dv*sum(abs(p)**3d0))
!  call shs(p,k1,1)
!  write(6,*) 'shs 1', sum(abs(k1)**3d0)
!  k1 = -ci*k1
!  y = p + dt/2d0 * k1 !y1
!  !wtime = mpi_wtime() - wtime
!  !write(6,*) 'rank, wtime one shs', rank, wtime
!
!  call shs(y,k2,1)
!  write(6,*) 'shs 2', sum(abs(k2)**3d0)
!  k2 = -ci*k2
!  y = p + dt/2d0 * k2 !y2
!
!  call shs(y,k4,1) !k4=k3
!  write(6,*) 'shs 3', sum(abs(k4)**3d0)
!  k4 = -ci*k4
!  y = p + dt * k4 !y3
!
!  k2 = k2 + k4 !k2+k3 
!
!  call shs(y,k4,1)
!  write(6,*) 'shs 4', sum(abs(k4)**3d0)
!  k4 = -ci*k4
!
!  p = p + dt/6d0*(k1+2d0*k2+k4)
!  nrm1 = sqrt(dv*sum(abs(p)**3d0))
!  p = p * nrm0/nrm1
!  write(6,*) 'is, nrm0, nrm1', is, nrm0, nrm1
!end subroutine

subroutine rk4_prop_shs_debug(p, pp)
  !time independent H for now
  ! p = unperturbed psi
  ! pp = perturbed psi
  use main_mod
  use tddft_mod 
  use mpi_lib_ours
  implicit none
  real*8  nrm0, nrm1
  complex*16 :: p(nx,ny,nz), pp(nx,ny,nz)
  complex*16, dimension(nx,ny,nz) :: k1,k2,k4,y
  complex*16, parameter :: ci = (0d0,1d0)
  !real*8 :: wtime, wtime2

  !wtime = mpi_wtime()
  !wtime2 = mpi_wtime()
  nrm0 = sqrt(dv*sum(abs(p)**2))
  if(rank==0) write(6,*) 'pre k1', sum(abs(p))
  call shs(p,k1,1)
  if(rank==0) write(6,*) 'k1', sum(abs(k1))
  k1 = -ci*k1
  y = p + dt/2d0 * k1 !y1
  !wtime = mpi_wtime() - wtime
  !write(6,*) 'rank, wtime one shs', rank, wtime

  call shs(y,k2,1)
  if(rank==0) write(6,*) 'k2', sum(abs(k2))
  k2 = -ci*k2
  y = p + dt/2d0 * k2 !y2

  call shs(y,k4,1) !k4=k3
  if(rank==0) write(6,*) 'k3', sum(abs(k4))
  k4 = -ci*k4
  y = p + dt * k4 !y3

  k2 = k2 + k4 !k2+k3 

  call shs(y,k4,1)
  if(rank==0) write(6,*) 'k4', sum(abs(k4))
  k4 = -ci*k4

  p = p + dt/6d0*(k1+2d0*k2+k4)
  nrm1 = sqrt(dv*sum(abs(p)**2))
  p = p * nrm0/nrm1


  nrm0 = sqrt(dv*sum(abs(pp)**2))
  call shs_pert(pp,k1,1)
  k1 = -ci*k1
  y = pp + dt/2d0 * k1 !y1

  call shs_pert(y,k2,1)
  k2 = -ci*k2
  y = pp + dt/2d0 * k2 !y2

  call shs_pert(y,k4,1) !k4=k3
  k4 = -ci*k4
  y = pp + dt * k4 !y3

  k2 = k2 + k4 !k2+k3 

  call shs_pert(y,k4,1)
  k4 = -ci*k4

  pp = pp + dt/6d0*(k1+2d0*k2+k4)
  nrm1 = sqrt(dv*sum(abs(pp)**2))
  pp = pp * nrm0/nrm1
  !wtime2 = mpi_wtime() - wtime2
  !write(6,*) 'rank, wtime 8 shs', rank, wtime2
end subroutine

!subroutine rk4_prop_shs_debug(is, p)
!  !time independent H for now
!  ! p = unperturbed psi
!  ! pp = perturbed psi
!  use main_mod
!  use tddft_mod 
!  use mpi_lib_ours
!  implicit none
!  integer :: is
!  real*8  nrm0, nrm1
!  complex*16 :: p(nx,ny,nz)
!  complex*16, dimension(nx,ny,nz) :: k1,k2,k4,y
!  complex*16, parameter :: ci = (0d0,1d0)
!  !real*8 :: wtime, wtime2
!
!  !wtime = mpi_wtime()
!  !wtime2 = mpi_wtime()
!  nrm0 = sqrt(dv*sum(abs(p)**3d0))
!  call shs(p,k1,1)
!  write(6,*) 'shs 1', sum(abs(k1)**3d0)
!  k1 = -ci*k1
!  y = p + dt/2d0 * k1 !y1
!  !wtime = mpi_wtime() - wtime
!  !write(6,*) 'rank, wtime one shs', rank, wtime
!
!  call shs(y,k2,1)
!  write(6,*) 'shs 2', sum(abs(k2)**3d0)
!  k2 = -ci*k2
!  y = p + dt/2d0 * k2 !y2
!
!  call shs(y,k4,1) !k4=k3
!  write(6,*) 'shs 3', sum(abs(k4)**3d0)
!  k4 = -ci*k4
!  y = p + dt * k4 !y3
!
!  k2 = k2 + k4 !k2+k3 
!
!  call shs(y,k4,1)
!  write(6,*) 'shs 4', sum(abs(k4)**3d0)
!  k4 = -ci*k4
!
!  p = p + dt/6d0*(k1+2d0*k2+k4)
!  nrm1 = sqrt(dv*sum(abs(p)**3d0))
!  p = p * nrm0/nrm1
!  write(6,*) 'is, nrm0, nrm1', is, nrm0, nrm1
!end subroutine
