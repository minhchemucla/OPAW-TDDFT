subroutine prop_dt_sh_debug(flg, p, pp)
  use main_mod
  implicit none
  integer :: flg
  complex*16 :: p(nx,ny,nz), pp(nx,ny,nz)

  if (flg == 1) then
    call prop_dt_sh_debug_1(p)
    call prop_dt_sh_debug_1_pert(pp)
  else if (flg == 2) then
    call prop_dt_sh_debug_2(p)
    call prop_dt_sh_debug_2_pert(pp)
  else if (flg == 3) then
    call prop_dt_sh_debug_3(p)
    call prop_dt_sh_debug_3_pert(pp)
  else
    write(6,*) 'only debug prop 1,2 and 3 implemented'
    stop
  endif

end subroutine

subroutine prop_dt_sh_debug_1(p1)
  use tddft_mod
  use main_mod
  implicit none
  real*8  nrm0, nrm1
  real*8, parameter :: toll = 1d-14
  complex*16, parameter :: ci=(0d0,1d0)
  complex*16 :: p1(nx,ny,nz), p0(nx,ny,nz), p2(nx,ny,nz), p3(nx,ny,nz)


  p0 = p1
  nrm0 = sqrt(dv*sum(abs(p1)**2))

  if(nrm0< toll) then
    p1 = 0d0
  else
    call sh(p1,p3,1)
    !call h_phi(p1, p2, 1)
    !p3 = p2
    !if (prop_bar) call sn_phi(p2,p1,1,-1)
    !call sn_phi(p2,p3,1,-1)
   !write(6,*) 'test minh', p1(1,1,1)
    p0 = p0 - ci*dt*p3
    nrm1 = sqrt(dv*sum(abs(p0)**2))
    p1 = p0 * (nrm0/nrm1)
  endif

end subroutine

subroutine prop_dt_sh_debug_1_pert(p1)
  use tddft_mod
  use main_mod
  implicit none
  real*8  nrm0, nrm1
  real*8, parameter :: toll = 1d-14
  complex*16, parameter :: ci=(0d0,1d0)
  complex*16 :: p1(nx,ny,nz), p0(nx,ny,nz), p2(nx,ny,nz), p3(nx,ny,nz)


  p0 = p1
  nrm0 = sqrt(dv*sum(abs(p1)**2))

  if(nrm0< toll) then
    p1 = 0d0
  else
    call sh_pert(p1, p3, 1)
    !call h_phi_pert(p1, p2, 1)
    !p3 = p2
    !if (prop_bar) call sn_phi(p2,p1,1,-1)
    !call sn_phi(p2,p3,1,-1)
    p0 = p0 - ci*dt*p3
    nrm1 = sqrt(dv*sum(abs(p0)**2))
    p1 = p0 * (nrm0/nrm1)
  endif

end subroutine

subroutine prop_dt_sh_debug_2(pin)
  use tddft_mod
  use main_mod
  use atom_mod
  use atom_mod, only : p=> pawinfo, at=> atominfo
  implicit none
  integer :: ik
  real*8  nrm0, nrm1
  real*8, parameter :: toll = 1d-14
  complex*16, parameter :: ci=(0d0,1d0)
  complex*16 :: pin(nx,ny,nz), p0(nx,ny,nz), p1(nx,ny,nz), p2(nx,ny,nz), p3(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)

  ik = 1
  p0 = pin
  nrm0 = sqrt(dv*sum(abs(pin)**2))

  if(nrm0< toll) then
    pin = 0d0
  else
    call vd_phi(p0,p1,ik)
    call sn_phi(p1, p2, 1, -1d0)
    p0 = p0 - ci*dt/2d0*p2

    cin=p0
    call fft3d_forward(nx,ny,nz,cin,cout)
    cout=cout*ek(:,:,:,ik)
    call fft3d_backward(nx,ny,nz,cout,p1)
    call sn_phi(p1, p2, 1, -1d0)
    p0 = p0 -ci*dt*p2

    call vd_phi(p0,p1,ik)
    call sn_phi(p1, p2, 1, -1d0)
    p0 = p0 - ci*dt/2d0*p2

    nrm1 = sqrt(dv*sum(abs(p0)**2))
    pin = p0 * (nrm0/nrm1)
  endif
end subroutine

subroutine prop_dt_sh_debug_2_pert(pin)
  use tddft_mod
  use main_mod
  use atom_mod
  use atom_mod, only : p=> pawinfo, at=> atominfo
  implicit none
  integer :: ik
  real*8  nrm0, nrm1
  real*8, parameter :: toll = 1d-14
  complex*16, parameter :: ci=(0d0,1d0)
  complex*16 :: pin(nx,ny,nz), p0(nx,ny,nz), p1(nx,ny,nz), p2(nx,ny,nz), p3(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)

  ik = 1
  p0 = pin
  nrm0 = sqrt(dv*sum(abs(pin)**2))

  if(nrm0< toll) then
    pin = 0d0
  else
    call vd_pert_phi(p0,p1,ik)
    call sn_phi(p1, p2, 1, -1d0)
    p0 = p0 - ci*dt/2d0*p2

    cin=p0
    call fft3d_forward(nx,ny,nz,cin,cout)
    cout=cout*ek(:,:,:,ik)
    call fft3d_backward(nx,ny,nz,cout,p1)
    call sn_phi(p1, p2, 1, -1d0)
    p0 = p0 -ci*dt*p2

    call vd_pert_phi(p0,p1,ik)
    call sn_phi(p1, p2, 1, -1d0)
    p0 = p0 - ci*dt/2d0*p2

    nrm1 = sqrt(dv*sum(abs(p0)**2))
    pin = p0 * (nrm0/nrm1)
  endif
end subroutine

subroutine prop_dt_sh_debug_3(p1)
  use tddft_mod
  use main_mod
  implicit none
  real*8  nrm0, nrm1
  real*8, parameter :: toll = 1d-14
  complex*16, parameter :: ci=(0d0,1d0)
  complex*16 :: p1(nx,ny,nz), p0(nx,ny,nz), p2(nx,ny,nz), p3(nx,ny,nz)
  complex*16 :: tmp(nx,ny,nz)


  p0 = p1
  nrm0 = sqrt(dv*sum(abs(p1)**2))

  if(nrm0< toll) then
    p1 = 0d0
  else
    call vd_pert_phi(p0,tmp,1)
    call sn_phi(tmp,p3,1,-1d0)
    p2 = p0 - ci*dt*p3
    call prop_tbar_split_sh(p2,p3,1)
    call vd_pert_phi(p3,tmp,1)
    call sn_phi(tmp,p2,1,-1d0)
    p0 = p3 - ci*dt*p2

    nrm1 = sqrt(dv*sum(abs(p0)**2))
    p1 = p0 * (nrm0/nrm1)
  endif
  !write(6,*) 'nrm0,nrm1', nrm0, nrm1

end subroutine

subroutine prop_dt_sh_debug_3_pert(p1)
  use tddft_mod
  use main_mod
  implicit none
  real*8  nrm0, nrm1
  real*8, parameter :: toll = 1d-14
  complex*16, parameter :: ci=(0d0,1d0)
  complex*16 :: p1(nx,ny,nz), p0(nx,ny,nz), p2(nx,ny,nz), p3(nx,ny,nz)
  complex*16 :: tmp(nx,ny,nz)


  p0 = p1
  nrm0 = sqrt(dv*sum(abs(p1)**2))

  if(nrm0< toll) then
    p1 = 0d0
  else
    call vd_phi(p0,tmp,1)
    call sn_phi(tmp,p3,1,-1d0)
    p2 = p0 - ci*dt*p3
    call prop_tbar_split_sh(p2,p3,1)
    call vd_phi(p3,tmp,1)
    call sn_phi(tmp,p2,1,-1d0)
    p0 = p3 - ci*dt*p2

    nrm1 = sqrt(dv*sum(abs(p0)**2))
    p1 = p0 * (nrm0/nrm1)
  endif
  !write(6,*) 'nrm0,nrm1 pert', nrm0, nrm1

end subroutine



subroutine applySinv(pin,pout)
  use main_mod
  use atom_mod
  use atom_mod, only : p=> pawinfo, at=> atominfo
  implicit none
  complex*16 :: pin(nx,ny,nz), pout(nx,ny,nz)
  integer :: ib,jb,ia,is,js,it,ik
  integer :: ix,iy,iz,ix1,iy1,iz1,jx,jy,jz
  integer :: ixs,iys,izs,ixe,iye,ize
  complex*16, allocatable :: ca(:)

  ik = 1
  pout = pin
  do ia=1,natom
      it=atom_map(ia)
      allocate(ca(p(it)%mstates),stat=stat)
      if(stat/=0) stop 'ca alloc problem in sphi'
      ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
      izs=at(ia)%ir_start(3)
      ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
      ize=izs+p(it)%nrough(3)-1
      call proj1(ia,pin,ca,p(it)%mstates,ik)
       !call proj(ia,pout,ca,p(it)%mstates)

      if(.not.at(ia)%edge) then
          do is=1,p(it)%mstates
             pout(ixs:ixe,iys:iye,izs:ize)=&
              pout(ixs:ixe,iys:iye,izs:ize)+&
               at(ia)%local_p3d1_c(:,:,:,is,ik)*at(ia)%sinv(is,ik)*&
              ! at(ia)%local_p3d(:,:,:,is)*p(it)%ssqinvij(is)*&
              ca(is)
          enddo
      else
          do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
              jx=mod(ix+nx-1,nx)+1
              jy=mod(iy+ny-1,ny)+1
              jz=mod(iz+nz-1,nz)+1
              pout(jx,jy,jz)=pout(jx,jy,jz)+&
                 sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,ik)*&
                 at(ia)%sinv(:,ik)*ca)
          enddo;enddo;enddo
      endif
      deallocate(ca)
  enddo
end subroutine


!prop_bar

subroutine prop_bar_dt_debug_1(p1)
  use tddft_mod
  use main_mod
  implicit none
  real*8  nrm0, nrm1
  real*8, parameter :: toll = 1d-14
  complex*16, parameter :: ci=(0d0,1d0)
  complex*16 :: p1(nx,ny,nz), p0(nx,ny,nz), p2(nx,ny,nz), p3(nx,ny,nz)


  p0 = p1
  nrm0 = sqrt(dv*sum(abs(p1)**2))

  if(nrm0< toll) then
    p1 = 0d0
  else
    call shs(p1,p3,1)
    !call h_phi(p1, p2, 1)
    !p3 = p2
    !if (prop_bar) call sn_phi(p2,p1,1,-1)
    !call sn_phi(p2,p3,1,-1)
   !write(6,*) 'test minh', p1(1,1,1)
    p0 = p0 - ci*dt*p3
    nrm1 = sqrt(dv*sum(abs(p0)**2))
    p1 = p0 * (nrm0/nrm1)
  endif

end subroutine

subroutine prop_bar_dt_debug_1_pert(p1)
  use tddft_mod
  use main_mod
  implicit none
  real*8  nrm0, nrm1
  real*8, parameter :: toll = 1d-14
  complex*16, parameter :: ci=(0d0,1d0)
  complex*16 :: p1(nx,ny,nz), p0(nx,ny,nz), p2(nx,ny,nz), p3(nx,ny,nz)


  p0 = p1
  nrm0 = sqrt(dv*sum(abs(p1)**2))

  if(nrm0< toll) then
    p1 = 0d0
  else
    call shs_pert(p1, p3, 1)
    !call h_phi_pert(p1, p2, 1)
    !p3 = p2
    !if (prop_bar) call sn_phi(p2,p1,1,-1)
    !call sn_phi(p2,p3,1,-1)
    p0 = p0 - ci*dt*p3
    nrm1 = sqrt(dv*sum(abs(p0)**2))
    p1 = p0 * (nrm0/nrm1)
  endif

end subroutine

