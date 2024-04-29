!
! This file contains two main subroutines, with others contained in them.
!  First, feitfleck -- does propagation in one time step of
!   a kinteic hamiltonian in 3d (to use 1d, set ny=nz=1; for 2d nz=1)
!
! In addition, driver_feit_fleck -- a complete program (for compilation,
! I have a "program," which calls this sub.

subroutine feit_fleck(cpsi, vv, &
     nx, massx, dx, &
     ny, massy, dy, &
     nz, massz, dz , &
     dt)
  
  implicit none
  integer nx, ny, nz, ix, iy, iz
  real*8  massx, massy, massz, dx, dy, dz, dt
  
  complex*16 cpsi(nx, ny, nz)
  complex*16, parameter :: ci = (0.d0, 1.d0)
  
  real*8      vv(nx, ny, nz)

  cpsi = cpsi* exp(-ci* dt/2.d0* vv) 
  call prop_k_3d()

  cpsi = cpsi* exp(-ci* dt/2.d0* vv) 
  
contains
  subroutine prop_k_3d()
    implicit none
    complex*16 c1d(2000)
    
    call check_le(max(nx, ny, nz),  2000, 'mnxyz2   ')
    
    do iz=1,nz
       do iy=1,ny
          c1d(1:nx) = cpsi(1:nx, iy, iz)
          
          call fft_prop_1d_good(c1d, nx, massx, dx, dt) ! out could=in.
          cpsi(1:nx, iy, iz) = c1d(1:nx)
       enddo
    enddo
    
    do iy=1, ny
       do ix=1, nx
          c1d(1:nz) = cpsi(ix,iy,1:nz)
          
          call fft_prop_1d_good(c1d, nz, massz, dz, dt) 
          cpsi(ix, iy, 1:nz) = c1d(1:nz)
       enddo
    enddo
    
    do iz=1, nz
       do ix=1,nx
          c1d(1:ny) = cpsi(ix, 1:ny, iz) 
          
          call fft_prop_1d_good(c1d, ny, massy, dy, dt) 
          cpsi(ix, 1:ny, iz) = c1d(1:ny)
       enddo
    enddo
    
  end subroutine prop_k_3d
end subroutine feit_fleck
  
!program drive_feit_fleck_simple
!  call driver_feit_fleck_simple
!end program drive_feit_fleck_simple

subroutine driver_feit_fleck_simple
  implicit none
  integer nx, ny, nz, nt
  integer ix, iy, iz, it, st
  
  complex*16 cterm
  real*8  dt, norm
  real*8  massx, massy, massz, dx, dy, dz, x0, y0, z0, xmin, ymin, zmin
  real*8 x, y, z, dV
  real*8   wx, wy, wz
  real*8 ave_x, ave_y, ave_z
  
  real*8, allocatable     ::   vv(:,:,:)
  complex*16, allocatable :: cpsi(:,:,:), cpsi0(:,:,:)
  
  call read_param_psi
  allocate(vv(nx, ny, nz), cpsi(nx, ny, nz), stat=st);call check(st,0,'strp  ')
  allocate(               cpsi0(nx, ny, nz), stat=st);call check(st,0,'strp  ')
  call init_psi
  call prep_vv
  cpsi0 = cpsi
  
  do it=1,nt
     call feit_fleck(cpsi, vv, &
          nx, massx, dx, &
          ny, massy, dy, &
          nz, massz, dz , &
          dt)   
     
     norm = sum(abs(cpsi)**2)*dV
     cterm = sum(conjg(cpsi0)*cpsi)*dV
     write(2,* ) it*dt,      cterm
     write(21,*) it*dt, dble(cterm)
     
     call sub_ave_x
     call sub_ave_y
     call sub_ave_z


     write(31,*)it*dt, ave_x
     write(32,*)it*dt, ave_y
     write(33,*)it*dt, ave_z
  enddo
contains
  
  subroutine prep_vv
    implicit none
    vv = 0.d0
  end subroutine prep_vv



  subroutine read_param_psi
    implicit none
    
    read(400,*)massx, massy, massz
    read(400,*)nx, ny, nz
    read(400,*)xmin, ymin, zmin
    read(400,*)wx, wy, wz
    read(400,*)dx, dy, dz
    read(400,*)x0,y0, z0
    read(400,*)nt, dt
    
    dv = 1.d0
    if(nx>1) dv = dv*dx
    if(ny>1) dv = dv*dy
    if(nz>1) dv = dv*dz
    
    
  end subroutine read_param_psi
  
  subroutine init_psi
    implicit none
    
    do ix=1,nx
       do iy=1, ny
          do iz=1, nz
             x = xmin + dx*(ix-1)
             y = ymin + dy*(iy-1)
             z = zmin + dz*(iz-1)
             cpsi(ix, iy, iz) = &
                  exp(-(x-x0)**2/4.d0/wx**2) * &
                  exp(-(y-y0)**2/4.d0/wy**2) * &
                  exp(-(z-z0)**2/4.d0/wz**2) 
             
          enddo
       enddo
    enddo
    cpsi = cpsi /sqrt(dv*sum(abs(cpsi**2)))
    
  end  subroutine init_psi
  
  subroutine sub_ave_x()
    implicit none
    
    ave_x = 0.d0
    do ix=1, nx
       x = (ix-1)*dx + xmin
       ave_x = ave_x + sum(abs(cpsi(ix,:,:))**2)*dV*x
    enddo
    ave_x = ave_x/norm
  end subroutine sub_ave_x
  
  subroutine sub_ave_y()
    implicit none
    
    ave_y = 0.d0
    do iy=1, ny
       y = (iy-1)*dy + ymin
       ave_y = ave_y + sum(abs(cpsi(:,iy,:))**2)*dV*y
    enddo
    ave_y = ave_y/norm
  end subroutine sub_ave_y
  
  subroutine sub_ave_z()
    implicit none
    
    ave_z = 0.d0
    do iz=1, nz
       z = (iz-1)*dz + zmin
       ave_z = ave_z + sum(abs(cpsi(:,:,iz))**2)*dV*z
    enddo
    ave_z = ave_z/norm
  end subroutine sub_ave_z
  
end subroutine driver_feit_fleck_simple

