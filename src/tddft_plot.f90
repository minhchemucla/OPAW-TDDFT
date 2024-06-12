subroutine plot_dip
  use mpi_lib_ours
  use main_mod,  only : dx, dy, dz, nx, ny, nz, dv, nn!, dens
  use main_mod,  only : nb, nocc
  use main_mod,  only : dens, h_type
  use tddft_mod, only : dens_pert, ipol, state_map, n_restart
  use tddft_mod, only : it, dt, sm, it_restart!, dens_pert
  use tddft_mod, only : prop_bar, phi_bar, phi_bar_tot
  use tddft_mod, only : phi_bar_pert, phi_bar_tot_pert, nvirtual
  !use main_mod,  only : dens
  !use tddft_mod, only : dens_pert
  implicit none
  integer ix,iy,iz,i
  real*8 dip(3), tmp
  real*8 dip_pert(3)
  real*8 d_dip(3), tmp3d(3)
  real*8 rr(3), rr2(3)
  real*8 :: rho(nn), rho_pert(nn)
  complex*16 :: sp(nx,ny,nz,nb), sp_pert(nx,ny,nz,nb) !p = phi
  !complex*16 :: pb(nx,ny,nz,nb*nodes), pb_pert(nx,ny,nz,nb*nodes) !pb= psi_bar for calculating density
  complex*16 :: pb(nx,ny,nz,nb), pb_pert(nx,ny,nz,nb) !pb= psi_bar for calculating density
  logical :: restart_again, first_restart
  !real*8 dens(nn), dens_pert(nn)


  !call make_dens_bar
   
  !write(6,*) 'sum(dens), sum(dens_pert) dipole plot', sum(dens)*dv, sum(dens_pert)*dv

  !if (.not. prop_bar) call make_dens_bar
  call make_dens
  if(rank==0) then
    
    !if (.not. prop_bar) call make_dens_AE !all electron
    dip=0d0
    dip_pert=0d0
    d_dip = 0d0
    
    i=0
    do iz=1,nz
       do iy=1,ny
          do ix=1,nx
             i = i+1
             rr = (/ (ix-1-nx/2)*dx, (iy-1-ny/2)*dy, (iz-1-nz/2)*dz /)           
             !rr2 =0d0
             !rr2(ipol) = rr(ipol)
             !dip    = dip   + dv * rr(:) * dens(i)
             !dip_pert   = dip_pert  + dv * rr(:) * dens_pert(i)
             dip    = dip   + dv * rr(:) * rho(i)
             dip_pert   = dip_pert  + dv * rr(:) * rho_pert(i)
             !dip    = dip   + dv * rr2(:) * dens(i)
             !dip_pert   = dip_pert  + dv * rr2(:) * dens_pert(i)
          enddo
       enddo
    enddo
    
    d_dip = (dip_pert - dip)/sm
    !d_dip = (dip_pert - dip)

    
    
    if(n_restart>0 .and. it==it_restart+1) then
      inquire(file='./d_dip_cont.dat',exist=restart_again)
      if(restart_again) then
        open(unit=15, file='d_dip_cont.dat',position='append',status='old',action='write')
        open(unit=16, file='dip_cont.dat',position='append',status='old',action='write')
      else
        inquire(file='./d_dip.dat',exist=first_restart)
        if(first_restart) then
          open(unit=18, file='d_dip.dat',status='old') 
          rewind(18)
          open(unit=19, file='dip.dat',status='old') 
          rewind(19)
          open(unit=15, file='d_dip_cont.dat')
          rewind(15)
          open(unit=16, file='dip_cont.dat')
          rewind(16)
          do i=1,it_restart
            read(18,*) tmp, tmp3d
            write(15,888) tmp, tmp3d
            read(19,*) tmp, tmp3d
            write(16,*) tmp, tmp3d
          enddo
        else
          open(unit=15, file='d_dip.dat')
          rewind(15)
          open(unit=16, file='dip.dat')
          rewind(16)
        endif
      endif
    else
      if(it==1) then
        open(unit=15, file='d_dip.dat') 
        rewind(15)
        open(unit=16, file='dip.dat') 
        rewind(16)
      endif
    endif

    write( 6,*)it,(it-1)*dt,  d_dip(:),'it, t,     d_dip '; call flush(6)
    write( 6,*)it,(it-1)*dt, dip, dip_pert,'it, t,     dip,    dip_pert '; call flush(16)
    write(15,888)(it-1)*dt,  d_dip(:); call flush(15)
    write(16,*)(it-1)*dt, dip, dip_pert; call flush(16)

    888 format (' ',4e18.8)
  end if

  call sync_mpi

contains 
  subroutine make_dens
    implicit none
    integer :: jb
    real*8  :: rhok(nx,ny,nz),rhok_pert(nx,ny,nz)

    if (h_type .eq. 0 ) then
      do i=1,nb
         call sn_phi(phi_bar(:,:,:,i,1),pb(:,:,:,i),1,0.5d0)
         call sn_phi(phi_bar_pert(:,:,:,i,1),pb_pert(:,:,:,i),1,0.5d0)
      enddo
    else
      pb=phi_bar(:,:,:,:,1)
      pb_pert=phi_bar_pert(:,:,:,:,1)
    endif

    rhok=0d0
    rhok_pert=0d0
    !if(rank .le. nodes-2) then
    !  do i=1,nb
    !    rhok=rhok+abs(pb(:,:,:,i))**2
    !    rhok_pert=rhok_pert+abs(pb_pert(:,:,:,i))**2
    !  enddo
    !else if (rank==nodes-1) then
    !  do i=1,nb-nvirtual
    !    rhok=rhok+abs(pb(:,:,:,i))**2
    !    rhok_pert=rhok_pert+abs(pb_pert(:,:,:,i))**2
    !  enddo
    !  !write(6,*) 'nb-nvirtual', nb-nvirtual
    !  !write(6,*) 'rhok', 2d0*sum(rhok)*dv
    !endif
    do i=1,nb
      jb = state_map(i)
      if(jb<=nocc) then
        rhok=rhok+abs(pb(:,:,:,i))**2
        rhok_pert=rhok_pert+abs(pb_pert(:,:,:,i))**2
      endif
    enddo
    rho=2d0*reshape(rhok, (/nn/))
    rho_pert=2d0*reshape(rhok_pert, (/nn/))
    call allsum_r8(rho,size(rho))
    call allsum_r8(rho_pert,size(rho_pert))
    !if(rank==0)write(6,*) 'maxval(dens), maxval(dens_pert) dipole plot', maxval(rho), maxval(rho_pert)
    !if(rank==0)write(6,*) 'minval(dens), minval(dens_pert) dipole plot', minval(rho), minval(rho_pert)
  end subroutine
end subroutine plot_dip
  

