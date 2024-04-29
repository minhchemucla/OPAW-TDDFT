! Note : see
! https://www.abinit.org/sites/default/files/infos/7.8/tutorial/lesson_paw1.html
! for overlap of paw spheres
subroutine map_localgrid
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use mpi_lib_ours 
    
    implicit none

    integer :: ix,iy,iz,ia,it,ig
    real*8  :: x,y,z,r,xc,yc,zc
    logical :: found

    ngrid=0

    do iz=1,nz;do iy=1,ny;do ix=1,nx
        call get_xyz
        do ia=1,natom
            call get_r

            it=atom_map(ia)
            if(r<p(it)%rcut) then
                ngrid(ia)=ngrid(ia)+1
            endif
        enddo
    enddo;enddo;enddo
    
    if(rank==0) write(*,*) '#. local grid points belonging to atoms: ',ngrid

    call alloc_localgrids

    do ia=1,natom
        
        ig=1
        it=atom_map(ia)

!        write(111,*) 'local grids for atom: ',ia

        do iz=1,nz;do iy=1,ny;do ix=1,nx
            call get_xyz
            call get_r
            if(ig>ngrid(ia)+1) then
                write(*,*) 'number of local grids near atom: ',ia,' wrong'
                stop
            endif
            if(r<p(it)%rcut) then
                at(ia)%local_grid(ig,:)=(/ix,iy,iz/)
!                write(110+ia,*) ig,ix,iy,iz
                ig=ig+1
            endif
        enddo;enddo;enddo
    enddo        

    if(rank==0) then
        write(*,*) 'finished mapping local grids'
        write(*,*) '============================'
        write(*,*)
    endif
contains

    subroutine alloc_localgrids
        implicit none
        integer :: it,ms,ng,qmmax

        do ia=1,natom
            it=atom_map(ia)
            ms=p(it)%mstates
            ng=ngrid(ia)
            qmmax=(2*p(it)%nl-1)**2
            
            allocate(at(ia)%local_grid(ng,3),stat=stat)
            if(stat/=0) stop 'problem alloc localgrids'
            allocate(at(ia)%local_g3d(ng,qmmax),stat=stat)
            if(stat/=0) stop 'problem alloc localg3d'
!            allocate(at(ia)%local_p3d(p(it)%nrough(1),&
!                p(it)%nrough(2),p(it)%nrough(3),ms),stat=stat)
!            if(stat/=0) stop 'problem alloc localp3d'
!            allocate(at(ia)%local_p3d1(p(it)%nrough(1),&
!                p(it)%nrough(2),p(it)%nrough(3),ms),stat=stat)
!            if(stat/=0) stop 'problem alloc localp3d1'
        enddo
    end subroutine alloc_localgrids

    subroutine get_r
        implicit none

        xc=min(abs(x-at(ia)%coord(1)),abs(x+xmax*2d0-at(ia)%coord(1)),abs(x-xmax*2d0-at(ia)%coord(1)))
        yc=min(abs(y-at(ia)%coord(2)),abs(y+ymax*2d0-at(ia)%coord(2)),abs(y-ymax*2d0-at(ia)%coord(2)))
        zc=min(abs(z-at(ia)%coord(3)),abs(z+zmax*2d0-at(ia)%coord(3)),abs(z-zmax*2d0-at(ia)%coord(3)))
        
        r=sqrt(xc**2+yc**2+zc**2)
    end subroutine get_r

    subroutine get_xyz
        implicit none

        x=-xmax+dble(ix-1)*dx
        y=-ymax+dble(iy-1)*dy
        z=-zmax+dble(iz-1)*dz
    end subroutine get_xyz

    subroutine check_2s
        implicit none

        integer :: istate,ir,jr,jm
        integer :: ms,ls
        real*8  :: wr1,wr2,ou

        call interpolate(p(it)%rr,p(it)%nr,r,ir,jr,wr1,wr2)
        ou=p(it)%ptilde(1,ir)*wr1+p(it)%ptilde(1,jr)*wr2
!        write(15,*) r,ou
    end subroutine check_2s

end subroutine map_localgrid
