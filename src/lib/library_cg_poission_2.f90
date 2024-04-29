!
! calcPhi_Conjugate Gradient.  Solves for phi in:    -grad(eps grad phi) = rho
! Uses a gradient routine from the grad_div_curl module in the library
!
! IMPORTANT: BEFORE CALLING THIS ROUTINE CALL
! call  set_grad_div_curl(nd, nx, ny, nz, dx, dy, dz)! 
! ALSO: Phi should be set to something, zero or a guess, before the program is run
subroutine calcPhi_library_CG(phi, rho, eps_s,nx,ny,nz,dx,dy,dz,nd)
  use grad_div_curl !debug
  implicit none
  
  real*8 dx,dy,dz
  integer nx,ny,nz,nd
  real*8, intent(in) :: rho(Nx, Ny, Nz), &
                      eps_s(Nx, Ny, Nz)
  real*8             :: phi(Nx, Ny, Nz), &
                    phi_tmp(Nx, Ny, Nz), &
                          d(Nx, Ny, Nz), & ! search direction
                          r(Nx, Ny, Nz), & ! residual
                       rold(Nx, Ny, Nz), & ! residual
                         Ad(Nx, Ny, Nz), & ! operate A on d
                         dold(nx,ny,nz), &
                       phiold(nx,ny,nz), &
                        rtr_oldold,      &
                         rTr, &            ! inner product of r and r
                         rTr0, &           ! initial rTr
                         rTr_old, &        ! last rTr
                         lambda, gamma, &  ! coefficients
                         threshold = 1d-10  ! allowed error
  integer            :: iStepIn
  integer, parameter :: nStep = 200

  real*8             :: objFunc(Nx, Ny, Nz), &
                       tmp_grad(Nx, Ny, Nz, Nd), &
                          rho_d(Nx, Ny, Nz)
  integer ix, iy, iz, id

  call check(nd,3,' nd   3  ')

  r = 0d0
  iStepIn = 1
  call calcAxeps(r, phi,eps_s)
  
  r = rho-r
  d = r
  rTr = sum(r*r)
  rTr0 = rTr
  rtr_old = rtr

  do while(iStepIn==1 .or. iStepIn <= nStep .and. rTr > threshold**2 * rTr0) !sum(abs(objFunc))>1) 

     call calcAxeps(Ad, d, eps_s)
     lambda = rTr / sum(d*Ad) 
     phi = phi+lambda*d
     if(mod(iStepIn, 10)==0) then
        call calcAxeps(r, phi,eps_s)
        r = rho-r
     else
        r = r-lambda*Ad
     endif
     rTr_old = rTr
     rTr = sum(r*r)
     gamma = rTr / rTr_old
     write(6,*)' gamma,rtr,lambda ',real(gamma),real(rtr),real(lambda)
     d = r+gamma*d
     iStepIn = iStepIn+1
!     call calcObjFunc(objFunc, phi, rho)
  enddo

  write(*,*) istepin, 'steps'!, real(sum(abs(objFunc))), '|objFunc|'

contains

  ! Ax = -div(eps_s*grad(x))
  subroutine calcAxeps(Ax, x, eps_s)
    use grad_div_curl
    implicit none
    real*8  ::        x(Nx, Ny, Nz), &
                     Ax(Nx, Ny, Nz), &
               tmp_grad(Nx, Ny, Nz, Nd), &
                  eps_s(nx, ny, nz), eps
    integer :: id, ix, iy, iz
    
    call r_grad_gen(x, tmp_grad, 1)
    do id = 1, Nd
       do iz = 1, Nz
          do iy = 1, Ny
             do ix = 1, Nx
                eps = eps_s( ix,iy,iz)
                tmp_grad(ix, iy, iz, id) = eps*tmp_grad(ix, iy, iz, id)
             enddo
          enddo
       enddo
    enddo
    call  r_div_gen(tmp_grad, Ax, 1)

    Ax = -Ax  ! -grad(eps grad) : need for positive dentie
  end subroutine calcAxeps

end subroutine calcPhi_library_CG

subroutine driver_example_calc_phi_cg
  use grad_div_curl
  implicit none
  integer :: ifirst=1; save ifirst
  integer, parameter :: nx=64, ny=32, nz=32, nd=3
  real*8,  parameter :: dx=0.6d0, dy=dx, dz=dx
  real*8             :: phi(nx,ny,nz), &
                        rho(nx,ny,nz), &
                       phia(nx,ny,nz), &
                         xa(nx,ny,nz), &
                         ya(nx,ny,nz), &
                         za(nx,ny,nz), &
                         r2(nx,ny,nz), &
                        eps_s(nx,ny,nz)
  real*8            :: pi, w,eps0
  integer :: ix,iy,iz

  if(ifirst==1) then
     call set_grad_div_curl(nd, nx, ny, nz, dx, dy, dz)
     ifirst=-1
  end if

  pi = dacos(-1.d0)
  eps0 = 1.d0/(4.d0*pi)
  eps_s = eps0

  w = 1*dx
  do iz=1,nz
     do iy=1,ny
        do ix=1,nx
           xa(ix,iy,iz) = (ix-(nx-1)/2.d0)*dx
           ya(ix,iy,iz) = (iy-(ny-1)/2.d0)*dy
           za(ix,iy,iz) = (iz-(nz-1)/2.d0)*dz
        enddo
     enddo
  enddo

  r2 = xa**2 + ya**2 + za**2

  phia = exp(-r2/2d0/w**2)
!  rho  = -eps0*phia/w**4 * (r2-3d0*w**2)
  call calcAxeps(rho, phia,eps_s)

  call  calcPhi_library_CG(phi, rho, eps_s,nx,ny,nz,dx,dy,dz,nd)

  do ix=1,nx
     iy = ny/2
     iz = nz/2
     write(6,*)xa(ix,iy,iz), rho(ix,iy,iz), phi(ix,iy,iz),phia(ix,iy,iz), ' x, rho, phi, phia '
  enddo

contains
  subroutine calcAxeps(Ax, x, eps_s)
    use grad_div_curl
    implicit none
    real*8  ::        x(Nx, Ny, Nz), &
         Ax(Nx, Ny, Nz), &
         tmp_grad(Nx, Ny, Nz, Nd), &
         eps_s(nx, ny, nz), eps
    integer :: id, ix, iy, iz
    
    call r_grad_gen(x, tmp_grad, 1)
    do id = 1, Nd
       do iz = 1, Nz
          do iy = 1, Ny
             do ix = 1, Nx
                eps = eps_s( ix,iy,iz)
                tmp_grad(ix, iy, iz, id) = eps*tmp_grad(ix, iy, iz, id)
             enddo
          enddo
       enddo
    enddo
    call  r_div_gen(tmp_grad, Ax, 1)
    
    Ax = -Ax  ! -grad(eps grad) : need for positive dentie
  end subroutine calcAxeps

end subroutine driver_example_calc_phi_cg

! grad2(phi) = -rho/eps0
! phi = exp(-r^2/2w^2)
! gradphi = -rvec /w^2 phi
! divgradphi = - div(rvec)/w^2 * phi - rvec grad phi/w^2 = -3/w^2 phi + r^2/w^4 phi = phi/w^4 *(r^2-3w^2)

!program driver_example_calc_phi_cg_prog
!  call driver_example_calc_phi_cg
!end program driver_example_calc_phi_cg_prog

