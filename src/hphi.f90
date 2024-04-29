subroutine h_phi(pin,hp,ik)
    use main_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use mpi_lib_ours
    implicit none

    integer :: ik
    complex*16 :: pin(nx,ny,nz)
    complex*16 :: hp(nx,ny,nz),tmp(nx,ny,nz)
    complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
    !real*8     :: wtime
    
    cin=pin
    !wtime = mpi_wtime() 
    call fft3d_forward(nx,ny,nz,cin,cout)
    cout=cout*ek(:,:,:,ik)
    call fft3d_backward(nx,ny,nz,cout,tmp)
    !wtime = mpi_wtime() - wtime
    !write(6,*) 'rank, wtime exp(T), hphi', rank, wtime
    hp=tmp+pin*vks
    call addvnl
contains
    subroutine addvnl
        implicit none

        integer :: ia,ig,ix,iy,iz,it,is,js,ix1,iy1,iz1
        integer :: ixs,iys,izs,ixe,iye,ize,jx,jy,jz
        complex*16, allocatable:: ca(:),arr(:)

        do ia=1,natom
            it=atom_map(ia)
            ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
            izs=at(ia)%ir_start(3)
            ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
            ize=izs+p(it)%nrough(3)-1
            allocate(ca(p(it)%mstates),arr(p(it)%mstates),stat=stat)
            if(stat/=0) stop 'ca alloc problem in addvnl'
!            call time_print('before proj')
            call proj(ia,pin,ca,p(it)%mstates,ik)
!            call time_print('after proj')
            if(.not.at(ia)%edge) then
                do is=1,p(it)%mstates
                 hp(ixs:ixe,iys:iye,izs:ize)=&
                     hp(ixs:ixe,iys:iye,izs:ize)+&
                     at(ia)%local_p3d_c(:,:,:,is,ik)*&
                     sum(at(ia)%dij(is,:)*ca(:))
                enddo
            else
                arr=matmul(at(ia)%dij,ca)
                do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
                    jx=mod(ix+nx-1,nx)+1
                    jy=mod(iy+ny-1,ny)+1
                    jz=mod(iz+nz-1,nz)+1
                      do is=1,p(it)%mstates
                        hp(jx,jy,jz)=hp(jx,jy,jz)+&
                        at(ia)%local_p3d_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,ik)*&
                        arr(is) 
                      enddo
                enddo;enddo;enddo
            endif
!            call time_print('after add dij')
            deallocate(ca,arr)
        enddo
    end subroutine
end subroutine h_phi
