!get all-electron partial wavefunctions on rough grid
subroutine get_local_aep3d
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use mpi_lib_ours, only : rank, nodes, bcast_r8
    implicit none

    integer :: ia, ig, it, ik, jk, is
    integer :: ix,iy,iz
    integer ip
    real*8  :: x,y,z,xc,yc,zc,r
    real*8  :: xx,yy,zz
    real*8, allocatable :: tmp3(:,:,:,:)
    !tmp3 -> nfine, nfine, nfine

    !if(rank==0) ip=6
    !if(rank>0)  ip = rank+86000

    !if(rank==0) write(*,*) 'sum of projectors*dv'
    do ia=1,natom
        it=atom_map(ia)
        call alloc_aep3d
        if(rank==0) then

            !tmp3 is phi on the fine grid

            allocate(tmp3(p(it)%nfine(1),p(it)%nfine(2),&
                p(it)%nfine(3),p(it)%mstates),stat=stat)
            if(stat/=0) stop 'tmp3 alloc problem'

            !obtain projector on fine grid (tmp3)
            do iz=1,p(it)%nfine(3)
                do iy=1,p(it)%nfine(2)
                    do ix=1,p(it)%nfine(1)
                        call get_xyz 
                        call get_radialpart !interpolate from radial to fine grid
                        call get_angularpart
                    enddo
                enddo
            enddo

            !spline the projector onto rough grid(local_aep3d)
            call get_aep3d_r !r:rough
            deallocate(tmp3)
        endif
        call bcast_r8(at(ia)%local_aep3d,size(at(ia)%local_aep3d),0)
        call aep3d_c !c:complex
    enddo       
    call check_ae_vs_pp

contains
    subroutine alloc_aep3d
        implicit none
        allocate(at(ia)%local_aep3d(p(it)%nrough(1),&
            p(it)%nrough(2),p(it)%nrough(3),p(it)%mstates),stat=stat)
        if(stat/=0) stop 'problem alloc local_aep3d'

    end subroutine alloc_aep3d

    subroutine aep3d_c
        implicit none

        integer :: nrough(3),irs(3),ms
        real*8  :: phase

        nrough=p(it)%nrough
        irs   =at(ia)%ir_start
        ms    =p(it)%mstates
        allocate(at(ia)%local_aep3d_c(nrough(1),nrough(2),nrough(3),ms,nk_loc)&
            ,stat=stat)
        if(stat/=0) stop 'local_aep3d_c alloc problem'
        do iz=1,nrough(3)
            do iy=1,nrough(2)
                do ix=1,nrough(1)
                    x=dble(irs(1)+ix-2)*dx-xmax
                    y=dble(irs(2)+iy-2)*dy-ymax
                    z=dble(irs(3)+iz-2)*dz-zmax
                    do ik=1,nk_loc
                        jk=(ik-1)*nodes+rank+1
                        phase=x*kpt(1,jk)+y*kpt(2,jk)+&
                            z*kpt(3,jk)
                        at(ia)%local_aep3d_c(ix,iy,iz,:,ik)=at(ia)%local_aep3d&
                            (ix,iy,iz,:)*cmplx(cos(phase),-sin(phase))
                    enddo
                enddo
            enddo
        enddo
    end subroutine aep3d_c

    subroutine get_aep3d_r
        implicit none

        real*8,  allocatable :: tmp(:,:,:),tmp1(:,:,:)

        integer :: nrough(3),nfine(3)
        integer :: i,j,is,ms

        nrough=p(it)%nrough
        nfine =p(it)%nfine
        ms    =p(it)%mstates

        allocate(tmp(nrough(1),nfine(2),nfine(3)),&
            tmp1(nrough(1),nrough(2),nfine(3)),stat=stat)
        if(stat/=0) stop 'tmp alloc problem aep3d_r'

        do is=1,ms
            do i=1,nfine(2)
                do j=1,nfine(3)
                    tmp(:,i,j)=matmul(p(it)%bx,tmp3(:,i,j,is))
                enddo
            enddo

            do i=1,nrough(1)
                do j=1,nfine(3)
                    tmp1(i,:,j)=matmul(p(it)%by,tmp(i,:,j))
                enddo
            enddo

            do i=1,nrough(1)
                do j=1,nrough(2)
                    at(ia)%local_aep3d(i,j,:,is)=matmul(p(it)%bz,tmp1(i,j,:))
                enddo
            enddo
        enddo

        deallocate(tmp1,tmp)
    end subroutine

    subroutine get_xyz
        implicit none
        ! calculates the x, y, z coordinates of the fine grid point from the origin of grid
        x=-xmax+dble(at(ia)%ir_start(1)-1)*dx+(ix-1)*p(it)%dxf(1)
        y=-ymax+dble(at(ia)%ir_start(2)-1)*dy+(iy-1)*p(it)%dxf(2)
        z=-zmax+dble(at(ia)%ir_start(3)-1)*dz+(iz-1)*p(it)%dxf(3)
        
        ! calculates the distance from atom to fine grid hile accounting for periodicity of system
        ! minbyabs calculates which of the three quantities is the smallest
        xc=minbyabs(x-at(ia)%coord(1),x+2d0*xmax-at(ia)%coord(1),x-2d0*xmax-at(ia)%coord(1))
        yc=minbyabs(y-at(ia)%coord(2),y+2d0*ymax-at(ia)%coord(2),y-2d0*ymax-at(ia)%coord(2))
        zc=minbyabs(z-at(ia)%coord(3),z+2d0*zmax-at(ia)%coord(3),z-2d0*zmax-at(ia)%coord(3))
        
!        write(160+ia,*) ix,iy,iz,x,y,z,xc,yc,zc

        r=max(sqrt(xc**2+yc**2+zc**2),1d-8)
!        write(160+ia,*) r
        ! used for the angular portion
        xx=xc/r
        yy=yc/r
        zz=zc/r
    end subroutine get_xyz

    subroutine get_radialpart
        implicit none
        integer :: istate,ir,jr,jm,ii !rr(ir)<r<rr(jr),jr=ir+1
        integer :: ms, ls
        real*8  :: wr1,wr2   !use linear interpolation, weight of ir,ir+1

        !simple linear fit
        call interpolate(p(it)%rr,p(it)%nr,r,ir,jr,wr1,wr2)

        ms=p(it)%mstates

        do istate=1,ms
            ls=p(it)%ms_ls(istate)
            tmp3(ix,iy,iz,istate)=p(it)%phi(ls,ir)*wr1+p(it)%phi(ls,jr)*wr2
        enddo
    end subroutine get_radialpart

    subroutine get_angularpart
        implicit none
        
        integer :: l,m,istate,ms,ls
        real*8  :: ylm
        
        ms=p(it)%mstates

        do istate=1,ms
            ls=p(it)%ms_ls (istate)
            m =p(it)%mstate(istate)
            l =p(it)%lstate(ls)
!            write(*,*) m,l
            
            call sphericalharmonics(l,m,xx,yy,zz,ylm)
!            write(18,*) 'xx,yy,zz,ylm',xx,yy,zz,ylm
            
            tmp3(ix,iy,iz,istate)= &
            tmp3(ix,iy,iz,istate)*ylm !!!note!!
        enddo
    end subroutine get_angularpart

    real*8 function minbyabs(x,y,z)
        implicit none
    
        real*8,intent(in)  :: x,y,z

        minbyabs=x    
        if (abs(y)<abs(minbyabs)) then !i and j are the index of the two closest radial grid  to the fine grid


            minbyabs=y
        endif
        if (abs(z)<abs(minbyabs)) then
            minbyabs=z
        endif
    end function        

    subroutine check_ae_vs_pp
      implicit none
      integer :: ia,ig,ix,iy,iz,it,is,js,ix1,iy1,iz1
      integer :: ixs,iys,izs,ixe,iye,ize,jx,jy,jz
      complex*16 :: test_wf(nx,ny,nz)

      test_wf = 0d0 
      
      do ia=1,natom
        it=atom_map(ia)
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        if(.not.at(ia)%edge) then
          do is=1,p(it)%mstates
            test_wf(ixs:ixe,iys:iye,izs:ize)=test_wf(ixs:ixe,iys:iye,izs:ize)&
              +at(ia)%local_aep3d_c(ixs:ixe,iys:iye,izs:ize,is,1)&
              -at(ia)%local_pp3d_c(ixs:ixe,iys:iye,izs:ize,is,1)
              !+at(ia)%local_pp3d_c(ixs:ixe,iys:iye,izs:ize,is,1)
          enddo
        else
          do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
            jx=mod(ix+nx-1,nx)+1
            jy=mod(iy+ny-1,ny)+1
            jz=mod(iz+nz-1,nz)+1
            do is=1,p(it)%mstates
              test_wf(jx,jy,jz)=test_wf(jx,jy,jz)&
              +at(ia)%local_aep3d_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)&
              -at(ia)%local_pp3d_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)
              !+at(ia)%local_pp3d_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)
            enddo
          enddo;enddo;enddo
        endif
      enddo

      if(rank==0) write(340,*) real(test_wf)

    end subroutine

end subroutine get_local_aep3d
