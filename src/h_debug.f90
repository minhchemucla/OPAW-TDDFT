subroutine get_h_debug
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p=> pawinfo, at=> atominfo
    use mpi_lib_ours, only : rank
    implicit none
    integer :: ib,jb,ik
    complex*16  :: tmp(nx,ny,nz),hphi(nx,ny,nz)
    real*8  :: hm(nb,nb)
    complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz)
    real*8  :: vxc3d(nx,ny,nz)

    vxc3d=reshape(vxc,(/nx,ny,nz/))
    
    do ik=1,nk_loc
        write(14000+rank+1,*) 'ek'
        do ib=1,nb
            cin=phit(:,:,:,ib,ik)
            call fft3d_forward(nx,ny,nz,cin,cout)
            cout=cout*ek(:,:,:,ik)
            call fft3d_backward(nx,ny,nz,cout,cin)
!            write(*,*) 'max,min imag kin',maxval(dimag(cin)),minval(dimag(cin))
            hphi=cin
            do jb=1,nb
                hm(ib,jb)=sum(hphi*conjg(phit(:,:,:,jb,ik)))*dv
            enddo
            write(14000+rank+1,*) hm(ib,ib)
        enddo       

        write(14000+rank+1,*) 'vion'
        do ib=1,nb
            hphi=phit(:,:,:,ib,ik)*vloc_tot
            do jb=1,nb
                hm(ib,jb)=sum(hphi*conjg(phit(:,:,:,jb,ik)))*dv
            enddo
            write(14000+rank+1,*) hm(ib,ib)
        enddo        

        write(14000+rank+1,*) 'vxc'
        do ib=1,nb
            hphi=phit(:,:,:,ib,ik)*vxc3d
            do jb=1,nb
                hm(ib,jb)=sum(hphi*conjg(phit(:,:,:,jb,ik)))*dv
            enddo
            write(14000+rank+1,*) hm(ib,ib)
        enddo        

        write(14000+rank+1,*) 'vh'
        
        do ib=1,nb
            hphi=phit(:,:,:,ib,ik)*vh
            do jb=1,nb
                hm(ib,jb)=sum(hphi*conjg(phit(:,:,:,jb,ik)))*dv
            enddo
            write(14000+rank+1,*) hm(ib,ib)
        enddo        

        write(14000+rank+1,*) 'vloc'
        do ib=1,nb
            hphi=phit(:,:,:,ib,ik)*vks
            do jb=1,nb
                hm(ib,jb)=sum(hphi*conjg(phit(:,:,:,jb,ik)))*dv
            enddo
            write(14000+rank+1,*) hm(ib,ib)
        enddo        

        write(14000+rank+1,*) 'vnl'
        do ib=1,nb
            hphi=0d0
            call addvnl
            do jb=1,nb
                hm(ib,jb)=sum(hphi*conjg(phit(:,:,:,jb,ik)))*dv
            enddo
            write(14000+rank+1,*) hm(ib,ib)
        enddo        

        write(14000+rank+1,*) 'eorb'
        write(15000+rank+1,*) 'hamiltonian'
        do ib=1,nb
            cin=phit(:,:,:,ib,ik)
            call fft3d_forward(nx,ny,nz,cin,cout)
            cout=cout*ek(:,:,:,ik)
            call fft3d_backward(nx,ny,nz,cout,cin)
!        write(*,*) 'max,min imag kin',maxval(dimag(cin)),minval(dimag(cin))
            tmp=cin
            hphi=tmp+phit(:,:,:,ib,ik)*vks
            call addvnl
            do jb=1,nb
                hm(ib,jb)=sum(hphi*conjg(phit(:,:,:,jb,ik)))*dv
            enddo
            write(14000+rank+1,*) hm(ib,ib)
            write(15000+rank+1,'(100(1x,f9.5))') hm(ib,:)
        enddo        
    enddo
contains
    subroutine addvnl
        implicit none

        integer :: ia,ig,ix,iy,iz,it,is,js,ix1,iy1,iz1
        integer :: ixs,iys,izs,ixe,iye,ize,jx,jy,jz
        complex*16, allocatable:: ca(:)

        do ia=1,natom
            it=atom_map(ia)
            ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
            izs=at(ia)%ir_start(3)
            ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
            ize=izs+p(it)%nrough(3)-1
            allocate(ca(p(it)%mstates),stat=stat)
            if(stat/=0) stop 'ca alloc problem in addvnl'
!            call time_print('before proj')
            call proj(ia,phit(:,:,:,ib,ik),ca,p(it)%mstates,ik)
!            call time_print('after proj')
            if(.not.at(ia)%edge) then
                do is=1,p(it)%mstates
                 hphi(ixs:ixe,iys:iye,izs:ize)=&
                     hphi(ixs:ixe,iys:iye,izs:ize)+&
                     at(ia)%local_p3d_c(:,:,:,is,ik)*&
                     sum(at(ia)%dij(is,:)*ca(:))
                enddo
            else
                do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
                    jx=mod(ix+nx-1,nx)+1
                    jy=mod(iy+ny-1,ny)+1
                    jz=mod(iz+nz-1,nz)+1
                    do is=1,p(it)%mstates
                        hphi(jx,jy,jz)=hphi(jx,jy,jz)+&
                        at(ia)%local_p3d_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,ik)*&
                        sum(at(ia)%dij(is,:)*ca(:))
                    enddo
                enddo;enddo;enddo
            endif
!            call time_print('after add dij')
            deallocate(ca)
        enddo
    end subroutine
end subroutine get_h_debug
