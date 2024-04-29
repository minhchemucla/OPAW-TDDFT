!subroutine tbar_proj(ia,pin,ca,ms)
!    use main_mod
!    use atom_mod
!    use atom_mod, only : p=> pawinfo, at=> atominfo
!    implicit none
!    
!    integer :: ia,ik
!    complex*16  :: ca(ms)
!    complex*16  :: pin(nx,ny,nz)
!    complex*16, allocatable :: local_pin(:,:,:)
!
!    integer :: ig,is,ms,ix,iy,iz,it,jx,jy,jz
!    integer :: ixs,iys,izs,ixe,iye,ize
!
!    it=atom_map(ia)
!    ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2);izs=at(ia)%ir_start(3)
!    ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1;ize=izs+p(it)%nrough(3)-1
!
!    if(.not.at(ia)%edge) then
!        do is=1,ms
!            ca(is)=sum(pin(ixs:ixe,iys:iye,izs:ize)&
!                *conjg(at(ia)%local_b_c(:,:,:,is)))*dv
!
!        enddo
!    else
!        allocate(local_pin(ixs:ixe,iys:iye,izs:ize),stat=stat)
!        if(stat/=0) stop 'local_pin alloc problem'
!        do iz=izs,ize
!            do iy=iys,iye
!                do ix=ixs,ixe
!                    jx=mod(ix+nx-1,nx)+1
!                    jy=mod(iy+ny-1,ny)+1
!                    jz=mod(iz+nz-1,nz)+1
!
!                    local_pin(ix,iy,iz)=pin(jx,jy,jz)
!                enddo
!            enddo
!         enddo
!        do is=1,ms
!            ca(is)=sum(local_pin(ixs:ixe,iys:iye,izs:ize)&
!                *conjg(at(ia)%local_b_c(:,:,:,is)))*dv
!        enddo
!        deallocate(local_pin)
!   endif
!end subroutine tbar_proj

subroutine tbar_proj2(ia,pin,ca,ms)
    use main_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    implicit none
    
    integer :: ia,ik
    complex*16  :: ca(ms)
    complex*16  :: pin(nx,ny,nz)
    complex*16, allocatable :: local_pin(:,:,:)

    integer :: ig,is,ms,ix,iy,iz,it,jx,jy,jz
    integer :: ixs,iys,izs,ixe,iye,ize

    it=atom_map(ia)
    ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2);izs=at(ia)%ir_start(3)
    ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1;ize=izs+p(it)%nrough(3)-1

    if(.not.at(ia)%edge) then
        do is=1,ms
            ca(is)=sum(pin(ixs:ixe,iys:iye,izs:ize)&
                *conjg(at(ia)%local_c_c(:,:,:,is)))*dv

        enddo
    else
        allocate(local_pin(ixs:ixe,iys:iye,izs:ize),stat=stat)
        if(stat/=0) stop 'local_pin alloc problem'
        do iz=izs,ize
            do iy=iys,iye
                do ix=ixs,ixe
                    jx=mod(ix+nx-1,nx)+1
                    jy=mod(iy+ny-1,ny)+1
                    jz=mod(iz+nz-1,nz)+1

                    local_pin(ix,iy,iz)=pin(jx,jy,jz)
                enddo
            enddo
         enddo
        do is=1,ms
            ca(is)=sum(local_pin(ixs:ixe,iys:iye,izs:ize)&
                *conjg(at(ia)%local_c_c(:,:,:,is)))*dv
        enddo
        deallocate(local_pin)
   endif
end subroutine tbar_proj2

subroutine tbar_proj3(ia,pin,ca,ms)
    use main_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    implicit none
    
    integer :: ia,ik
    complex*16  :: ca(ms)
    complex*16  :: pin(nx,ny,nz)
    complex*16, allocatable :: local_pin(:,:,:)

    integer :: ig,is,ms,ix,iy,iz,it,jx,jy,jz
    integer :: ixs,iys,izs,ixe,iye,ize

    it=atom_map(ia)
    ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2);izs=at(ia)%ir_start(3)
    ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1;ize=izs+p(it)%nrough(3)-1

    if(.not.at(ia)%edge) then
        do is=1,ms
            ca(is)=sum(pin(ixs:ixe,iys:iye,izs:ize)&
                *conjg(at(ia)%local_c1_c(:,:,:,is)))*dv

        enddo
    else
        allocate(local_pin(ixs:ixe,iys:iye,izs:ize),stat=stat)
        if(stat/=0) stop 'local_pin alloc problem'
        do iz=izs,ize
            do iy=iys,iye
                do ix=ixs,ixe
                    jx=mod(ix+nx-1,nx)+1
                    jy=mod(iy+ny-1,ny)+1
                    jz=mod(iz+nz-1,nz)+1

                    local_pin(ix,iy,iz)=pin(jx,jy,jz)
                enddo
            enddo
         enddo
        do is=1,ms
            ca(is)=sum(local_pin(ixs:ixe,iys:iye,izs:ize)&
                *conjg(at(ia)%local_c1_c(:,:,:,is)))*dv
        enddo
        deallocate(local_pin)
   endif
end subroutine tbar_proj3
