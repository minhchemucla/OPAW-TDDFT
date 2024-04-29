!projector and shape function
subroutine get_local_p3d
    ! calculates the local 3d projector functions around each atom
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
    real*8, allocatable :: tmp3(:,:,:,:),tmp3_1(:,:,:,:)
    !tmp3 -> nfine, nfine, nfine

    if(rank==0) ip=6
    if(rank>0)  ip = rank+86000

    if(rank==0) write(*,*) 'sum of projectors*dv'
    do ia=1,natom
        it=atom_map(ia)
        call alloc_p3d
        if(rank==0) then

            !tmp3 and tmp3_1 are ptilde and ptilde1 on the fine grid

            allocate(tmp3(p(it)%nfine(1),p(it)%nfine(2),&
                p(it)%nfine(3),p(it)%mstates),stat=stat)
            if(stat/=0) stop 'tmp3 alloc problem'

            allocate(tmp3_1(p(it)%nfine(1),p(it)%nfine(2),&
                p(it)%nfine(3),p(it)%mstates),stat=stat)
            if(stat/=0) stop 'tmp3_1 alloc problem'

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

            !spline the projector onto rough grid(local_p3d)
            call get_p3d_r !r:rough
            deallocate(tmp3,tmp3_1)
!            do is=1,p(it)%mstates
!                write(*,*) 'ia,is',sum(at(ia)%local_p3d(:,:,:,is))*dv
!            enddo 
            ! transform 
            call transform_p3d(ia)
        endif
        call bcast_r8(at(ia)%local_p3d,size(at(ia)%local_p3d),0)
        call bcast_r8(at(ia)%local_p3d1,size(at(ia)%local_p3d1),0)
        call p3d_c !c:complex
    enddo       

contains
    subroutine alloc_p3d
        implicit none
        allocate(at(ia)%local_p3d(p(it)%nrough(1),&
            p(it)%nrough(2),p(it)%nrough(3),p(it)%mstates),stat=stat)
        if(stat/=0) stop 'problem alloc localp3d'
        allocate(at(ia)%local_p3d1(p(it)%nrough(1),&
            p(it)%nrough(2),p(it)%nrough(3),p(it)%mstates),stat=stat)
        if(stat/=0) stop 'problem alloc localp3d1'

    end subroutine alloc_p3d

    subroutine p3d_c
        implicit none

        integer :: nrough(3),irs(3),ms
        real*8  :: phase

        nrough=p(it)%nrough
        irs   =at(ia)%ir_start
        ms    =p(it)%mstates
        allocate(at(ia)%local_p3d_c(nrough(1),nrough(2),nrough(3),ms,nk_loc)&
            ,stat=stat)
        if(stat/=0) stop 'local_p3d_c alloc problem'
        allocate(at(ia)%local_p3d1_c(nrough(1),nrough(2),nrough(3),ms,nk_loc)&
            ,stat=stat)
        if(stat/=0) stop 'local_p3d1_c alloc problem'
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
                        at(ia)%local_p3d_c(ix,iy,iz,:,ik)=at(ia)%local_p3d&
                            (ix,iy,iz,:)*cmplx(cos(phase),-sin(phase))
                        at(ia)%local_p3d1_c(ix,iy,iz,:,ik)=at(ia)%local_p3d1&
                            (ix,iy,iz,:)*cmplx(cos(phase),-sin(phase))
                    enddo
                enddo
            enddo
        enddo
    end subroutine p3d_c

    subroutine get_p3d_r
        implicit none

        real*8,  allocatable :: tmp(:,:,:),tmp1(:,:,:)
        real*8,  allocatable :: tmp_1(:,:,:),tmp1_1(:,:,:)
!        real*8,  allocatable :: xr1(:),y2a1(:),b1(:,:)
!        real*8,  allocatable :: xr2(:),y2a2(:),b2(:,:)
!        real*8,  allocatable :: xr3(:),y2a3(:),b3(:,:)
!        real*8,  allocatable :: yr1(:),yr2(:),yf3(:)

        integer :: nrough(3),nfine(3)
        integer :: i,j,is,ms

        nrough=p(it)%nrough
        nfine =p(it)%nfine
        ms    =p(it)%mstates

        allocate(tmp(nrough(1),nfine(2),nfine(3)),&
            tmp1(nrough(1),nrough(2),nfine(3)),stat=stat)
        if(stat/=0) stop 'tmp alloc problem p3d_r'
        allocate(tmp_1(nrough(1),nfine(2),nfine(3)),&
            tmp1_1(nrough(1),nrough(2),nfine(3)),stat=stat)
        if(stat/=0) stop 'tmp alloc problem p3d_r'

        do is=1,ms
            do i=1,nfine(2)
                do j=1,nfine(3)
                    tmp(:,i,j)=matmul(p(it)%bx,tmp3(:,i,j,is))
                    tmp_1(:,i,j)=matmul(p(it)%bx,tmp3_1(:,i,j,is))
                enddo
            enddo

            do i=1,nrough(1)
                do j=1,nfine(3)
                    tmp1(i,:,j)=matmul(p(it)%by,tmp(i,:,j))
                    tmp1_1(i,:,j)=matmul(p(it)%by,tmp_1(i,:,j))
                enddo
            enddo

            do i=1,nrough(1)
                do j=1,nrough(2)
                    at(ia)%local_p3d(i,j,:,is)=matmul(p(it)%bz,tmp1(i,j,:))
                    at(ia)%local_p3d1(i,j,:,is)=matmul(p(it)%bz,tmp1_1(i,j,:))
                enddo
            enddo
        enddo

        deallocate(tmp1,tmp1_1,tmp,tmp_1)
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
!            call splint_dn(p(it)%rr,p(it)%ptilde(ls,:),&
!                p(it)%y2a_p(:,ls),p(it)%nr,r,tmp3(ix,iy,iz,istate))
!            call splint_dn(p(it)%rr,p(it)%ptilde1(ls,:),&
!                p(it)%y2a_p1(:,ls),p(it)%nr,r,tmp3_1(ix,iy,iz,istate))
            tmp3(ix,iy,iz,istate)=p(it)%ptilde(ls,ir)*wr1 &
                +p(it)%ptilde(ls,jr)*wr2
            tmp3_1(ix,iy,iz,istate)=p(it)%ptilde1(ls,ir)*wr1 &
                +p(it)%ptilde1(ls,jr)*wr2
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
            tmp3_1(ix,iy,iz,istate)= &
            tmp3_1(ix,iy,iz,istate)*ylm !!!note!!
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

end subroutine get_local_p3d
