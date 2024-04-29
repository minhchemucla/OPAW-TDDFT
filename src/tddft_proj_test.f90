subroutine proj_phi(pin,sp,ik)
    use main_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use mpi_lib_ours, only : rank
    implicit none

    integer :: ib,jb,ia,ix,iy,iz,is,js,it,ix1,iy1,iz1,ig,ik
    complex*16 :: pin(nx,ny,nz),sp(nx,ny,nz)
    complex*16, allocatable :: ca(:)
    integer :: ixs,iys,izs,ixe,iye,ize,jx,jy,jz

    sp=0d0
    do ia=1,natom
      it=atom_map(ia)
      ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
      izs=at(ia)%ir_start(3)
      ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
      ize=izs+p(it)%nrough(3)-1
      allocate(ca(p(it)%mstates),stat=stat)
      if(stat/=0) stop 'ca alloc problem in sphi'
      call proj1(ia,pin,ca,p(it)%mstates,ik)
        if(.not.at(ia)%edge) then
          do is=1,p(it)%mstates
            sp(ixs:ixe,iys:iye,izs:ize)=&
             sp(ixs:ixe,iys:iye,izs:ize)+&
             at(ia)%local_p3d1_c(:,:,:,is,ik)*ca(is)
          enddo
        else
          do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
            jx=mod(ix+nx-1,nx)+1
            jy=mod(iy+ny-1,ny)+1
            jz=mod(iz+nz-1,nz)+1
            sp(jx,jy,jz)=sp(jx,jy,jz)+&
              sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,ik)*ca)
          enddo;enddo;enddo
        endif
      deallocate(ca)
    enddo
end subroutine proj_phi
