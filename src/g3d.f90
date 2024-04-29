subroutine get_localg3d
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use mpi_lib_ours, only : rank, nodes, bcast_r8, bcast_i
    implicit none

    integer :: ia,ig,it,ix,iy,iz
    real*8  :: x,y,z,xc,yc,zc,r,xx,yy,zz
    integer :: qmmax

    do ia=1,natom
        it=atom_map(ia)
        if(rank==0) then
            do ig=1,ngrid(ia)
!               if(is_it_power2(ig)==1) call writef(' grid_ia ')
                call get_xyz
                call get_radialpart
                call get_angularpart
            enddo
            qmmax=(2*p(it)%nl-1)**2
            do ig=1,qmmax
              write(420,*) 'ia,ig,sumlocalg',ia,ig,sum(at(ia)%local_g3d(:,ig))
            enddo
        else
            qmmax=(2*p(it)%nl-1)**2
            allocate(at(ia)%local_grid(ngrid(ia),3),stat=stat)
            if(stat/=0) stop 'problem alloc localgrids'
            allocate(at(ia)%local_g3d(ngrid(ia),qmmax),stat=stat)
            if(stat/=0) stop 'problem alloc localg3d'
        endif
        !write(*,*) 'sizeg3d',rank,ia,size(at(ia)%local_g3d)
        call bcast_r8(at(ia)%local_g3d,size(at(ia)%local_g3d),0)
        call bcast_i (at(ia)%local_grid,size(at(ia)%local_grid),0)
    enddo
    if(rank==0) then
        write(*,*) 'finished mapping projectors on local grids'
        write(*,*) '=========================================='
        write(*,*)
    endif
contains        
    subroutine get_angularpart
        implicit none

        integer :: l,m,istate,ms,ls
        real*8  :: ylm
        
        do istate=1,(2*p(it)%nl-1)**2
            l=floor(sqrt(dble(istate-1)))
            m=istate-l**2-l-1
            call sphericalharmonics(l,m,xx,yy,zz,ylm)
!            write(6000+istate*10,*) ig,at(ia)%local_g3d(ig,istate),ylm
            at(ia)%local_g3d(ig,istate)=at(ia)%local_g3d(ig,istate)*ylm
        enddo
    end subroutine get_angularpart

    subroutine get_radialpart
        implicit none

        integer :: istate,ii,ir,jr!,jm,ii !rr(ir)<r<rr(jr),jr=ir+1
        integer :: l
        real*8  :: wr1,wr2   !use linear interpolation, weight of ir,ir+1

        if(r>p(it)%rcomp) then
            at(ia)%local_g3d(ig,:)=0d0
        else
            call interpolate(p(it)%rr,p(it)%nr,r,ir,jr,wr1,wr2)
            do istate=1,(2*p(it)%nl-1)**2
                l=floor(sqrt(dble(istate-1)))
                at(ia)%local_g3d(ig,istate)=&
                    p(it)%gl(ir,l+1)*wr1+p(it)%gl(jr,l+1)*wr2
            enddo
        endif

    end subroutine get_radialpart

    subroutine get_xyz
        implicit none

        ix=at(ia)%local_grid(ig,1)
        iy=at(ia)%local_grid(ig,2)
        iz=at(ia)%local_grid(ig,3)

        x=-xmax+dble(ix-1)*dx
        y=-ymax+dble(iy-1)*dy
        z=-zmax+dble(iz-1)*dz

!        xc=x-at(ia)%coord(1)
!        yc=y-at(ia)%coord(2)
!        zc=z-at(ia)%coord(3)
        xc=minbyabs(x-at(ia)%coord(1),x+2d0*xmax-at(ia)%coord(1),x-2d0*xmax-at(ia)%coord(1))
        yc=minbyabs(y-at(ia)%coord(2),y+2d0*ymax-at(ia)%coord(2),y-2d0*ymax-at(ia)%coord(2))
        zc=minbyabs(z-at(ia)%coord(3),z+2d0*zmax-at(ia)%coord(3),z-2d0*zmax-at(ia)%coord(3))
        
!        write(112,*) 'ix,iy,iz,xc,yc,zc',ix,iy,iz,xc,yc,zc

        r=max(sqrt(xc**2+yc**2+zc**2),1d-8)
!        write(160+ia,*) r
        xx=xc/r
        yy=yc/r
        zz=zc/r
    end subroutine get_xyz

    real*8 function minbyabs(x,y,z)
        implicit none
    
        real*8,intent(in)  :: x,y,z

        minbyabs=x    
        if (abs(y)<abs(minbyabs)) then
            minbyabs=y
        endif
        if (abs(z)<abs(minbyabs)) then
            minbyabs=z
        endif
    end function        

end subroutine      
