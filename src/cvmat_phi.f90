subroutine cvmat_phi(pin,sp,ik)
    use main_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use mpi_lib_ours, only : rank
    integer :: ib,jb,ia,ix,iy,iz,is,js,it,ix1,iy1,iz1,ig,ms
    complex*16 :: pin(nx,ny,nz),sp(nx,ny,nz)
    complex*16, allocatable :: ca(:)

    sp=0d0
    do ia=1,natom
      it=atom_map(ia)
      ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
      izs=at(ia)%ir_start(3)
      ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
      ize=izs+p(it)%nrough(3)-1
      ms = p(it)%mstates

      allocate(ca(ms),stat=stat)
      if(stat/=0) stop 'ca alloc problem in prep_tbar'

      !call tbar_proj3(ia, pin, ca, ms)
      do is=1,ms
        ca(is) = sum(conjg(at(ia)%local_cv1_c(:,:,:,is))*pin)*dv
      enddo

      if(.not.at(ia)%edge) then
        do is=1,ms
          sp(ixs:ixe,iys:iye,izs:ize)=sp(ixs:ixe,iys:iye,izs:ize)+&
            at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cvmat(is,:)*ca(:))
        enddo
      else
        do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
          jx=mod(ix+nx-1,nx)+1
          jy=mod(iy+ny-1,ny)+1
          jz=mod(iz+nz-1,nz)+1
          do is=1,ms
            sp(jx,jy,jz)=sp(jx,jy,jz)+at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)&
                *sum(at(ia)%cvmat(is,:)*ca(:))
          enddo
        enddo;enddo;enddo
      endif
      deallocate(ca)
    enddo
end subroutine cvmat_phi

subroutine add_cvtaylor_phi(pin,sp,ik)
    use main_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use mpi_lib_ours, only : rank
    integer :: ib,jb,ia,ix,iy,iz,is,js,it,ix1,iy1,iz1,ig,ms
    complex*16 :: pin(nx,ny,nz),sp(nx,ny,nz)
    complex*16, allocatable :: ca(:)

    sp=pin
    do ia=1,natom
      it=atom_map(ia)
      ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
      izs=at(ia)%ir_start(3)
      ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
      ize=izs+p(it)%nrough(3)-1
      ms = p(it)%mstates

      allocate(ca(ms),stat=stat)
      if(stat/=0) stop 'ca alloc problem in prep_tbar'

      do is=1,ms
        ca(is) = sum(conjg(at(ia)%local_cv1_c(:,:,:,is))*pin)*dv
      enddo

      if(.not.at(ia)%edge) then
        do is=1,ms
          sp(ixs:ixe,iys:iye,izs:ize)=sp(ixs:ixe,iys:iye,izs:ize)+&
            at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cv_taylor(is,:)*ca(:))
        enddo
      else
        do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
          jx=mod(ix+nx-1,nx)+1
          jy=mod(iy+ny-1,ny)+1
          jz=mod(iz+nz-1,nz)+1
          do is=1,ms
            sp(jx,jy,jz)=sp(jx,jy,jz)+at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)&
                *sum(at(ia)%cv_taylor(is,:)*ca(:))
          enddo
        enddo;enddo;enddo
      endif
      deallocate(ca)
    enddo
end subroutine add_cvtaylor_phi

subroutine add_cvtaylor_phi1(pin,sp,ik)
    use main_mod
    use atom_mod
    use tddft_mod, only : tbar_taylor_power
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use mpi_lib_ours, only : rank
    integer :: ib,jb,ia,ix,iy,iz,is,js,it,ix1,iy1,iz1,ig,ms, i
    complex*16 :: pin(nx,ny,nz),sp(nx,ny,nz)
    complex*16 :: tmp(nx,ny,nz), tmp1(nx,ny,nz)
    complex*16, allocatable :: ca(:)

    sp=pin
    tmp=pin
    do i=1,tbar_taylor_power
      tmp1=0d0
      do ia=1,natom
        it=atom_map(ia)
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        ms = p(it)%mstates

        allocate(ca(ms),stat=stat)
        if(stat/=0) stop 'ca alloc problem in prep_tbar'

        do is=1,ms
          ca(is) = sum(conjg(at(ia)%local_cv1_c(:,:,:,is))*tmp)*dv/dble(i)
        enddo

        if(.not.at(ia)%edge) then
          do is=1,ms
            tmp1(ixs:ixe,iys:iye,izs:ize)=tmp1(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cvbarmat(is,:)*ca(:))
          enddo
        else
          do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
            jx=mod(ix+nx-1,nx)+1
            jy=mod(iy+ny-1,ny)+1
            jz=mod(iz+nz-1,nz)+1
            do is=1,ms
              tmp1(jx,jy,jz)=tmp1(jx,jy,jz)+at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)&
                  *sum(at(ia)%cvbarmat(is,:)*ca(:))
            enddo
          enddo;enddo;enddo
        endif
        deallocate(ca)
      enddo
      sp = sp + tmp1
      tmp = tmp1
    enddo
end subroutine add_cvtaylor_phi1

subroutine cv_phi(pin,sp,ik)
    use main_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use mpi_lib_ours, only : rank
    integer :: ib,jb,ia,ix,iy,iz,is,js,it,ix1,iy1,iz1,ig,ms
    complex*16 :: pin(nx,ny,nz),sp(nx,ny,nz)
    complex*16, allocatable :: ca(:)

    sp=0d0
    do ia=1,natom
      it=atom_map(ia)
      ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
      izs=at(ia)%ir_start(3)
      ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
      ize=izs+p(it)%nrough(3)-1
      ms = p(it)%mstates

      allocate(ca(ms),stat=stat)
      if(stat/=0) stop 'ca alloc problem in prep_tbar'

      !call tbar_proj2(ia, pin, ca, ms)
      do is=1,ms
        ca(is) = sum(conjg(at(ia)%local_cv_c(:,:,:,is))*pin)*dv
      enddo

      if(.not.at(ia)%edge) then
        do is=1,ms
          sp(ixs:ixe,iys:iye,izs:ize)=sp(ixs:ixe,iys:iye,izs:ize)+&
            at(ia)%local_p3d1_c(:,:,:,is,ik)*at(ia)%sinv(is,ik)*ca(is)
        enddo
      else
        do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
          jx=mod(ix+nx-1,nx)+1
          jy=mod(iy+ny-1,ny)+1
          jz=mod(iz+nz-1,nz)+1
          sp(jx,jy,jz)=sp(jx,jy,jz)+sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,ik)&
              *at(ia)%sinv(:,ik)*ca(:))
        enddo;enddo;enddo
      endif
      deallocate(ca)
    enddo
end subroutine cv_phi

subroutine cv_taylor_loop(pin,pout,ik)
  use main_mod
  use tddft_mod, only : tbar_taylor_power, dt
  use mpi_lib_ours
  use atom_mod, only : p=> pawinfo, at=> atominfo, natom, atom_map
  implicit none
  integer :: ik
  complex*16 :: pin(nx,ny,nz),pout(nx,ny,nz)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
  integer :: i
  complex*16 :: ci=(0d0,1d0)

  pout = pin
  cin = pin
  do i=1,tbar_taylor_power
    call cv_phi(cin,cout,ik) 
    cin = -ci*dt*0.5d0*cout/dble(i)
    pout = pout + cin
  enddo
end subroutine
