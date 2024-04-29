subroutine get_nhat
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    implicit none

    real*8 :: a,b
    integer :: ia,ig,it,is,js,igl,ix,iy,iz

    nhat=0d0

    do ia=1,natom
        it=atom_map(ia)
        do ig=1,ngrid(ia)
            ix=at(ia)%local_grid(ig,1)
            iy=at(ia)%local_grid(ig,2)
            iz=at(ia)%local_grid(ig,3)
            do is=1,p(it)%mstates
                do js=1,p(it)%mstates
                    do igl=1,(2*p(it)%nl-1)**2
                        nhat(ix,iy,iz)=nhat(ix,iy,iz)+at(ia)%rhoij(is,js) &
                            *p(it)%qijlm(igl,is,js)*at(ia)%local_g3d(ig,igl)
                    enddo
                enddo
            enddo
        enddo
    enddo

    !write(*,*) 'nhat',sum(nhat**3)
    if(iscf>0) then
        a=sum(nhat)*dv
        b=sum(dens)*dv
        !write(*,*) 'sum dens',b
        !write(*,*) 'sum dens+nhat',a+b
        nhat=nhat*(rnel)/(a+b)
        dens=dens*(rnel)/(a+b)
        a=sum(nhat)*dv
        b=sum(dens)*dv
        !write(*,*) 'adjusted sum dens+nhat',a+b
    endif

end subroutine get_nhat      

subroutine get_nhat_elem(psi_i,psi_j,nhatij)
    !nhat_ij(r) = sum_LM sum_a sum_is,ij Q_ij^{a,LM}(r)<phit_i|p_is><p_ij|phit_j>
    use main_mod
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    implicit none

    integer :: i,j,ia,ig,it,is,js,igl,ix,iy,iz,ms
    !real*8 :: nhatij(nx,ny,nz)
    complex*16 :: nhatij(nx,ny,nz)
    complex*16 :: psi_i(nx,ny,nz), psi_j(nx,ny,nz)
    complex*16, allocatable :: bra(:), ket(:)

    if(any(shape(nhatij)/=(/nx,ny,nz/))) then
      write(6,*) 'shape(nhatij)', shape(nhatij)
      write(6,*) 'not matching nx,ny,nz:', nx,ny,nz
    endif


    nhatij=0d0
    !write(6,*) 'ia, i, j', ia, i,j
    do ia=1,natom
      it=atom_map(ia)
      ms=p(it)%mstates
      allocate(bra(ms),ket(ms),stat=stat)
      if(stat/=0) stop 'bra alloc problem in nhat'
      !only programmed for gamma point and nspin=1
      !bra = <p_is|
      call proj(ia,psi_i,bra,ms,1) 
      call proj(ia,psi_j,ket,ms,1) 
      !if(i.ne.j) then
      !else
      !  ket = bra
      !endif
      do ig=1,ngrid(ia)
        ix=at(ia)%local_grid(ig,1)
        iy=at(ia)%local_grid(ig,2)
        iz=at(ia)%local_grid(ig,3)
        do is=1,ms
          do js=1,ms
            do igl=1,(2*p(it)%nl-1)**2
              nhatij(ix,iy,iz)=nhatij(ix,iy,iz)+conjg(bra(is))*ket(js) &
                *p(it)%qijlm(igl,is,js)*at(ia)%local_g3d(ig,igl)
            enddo
          enddo
        enddo
      enddo
      deallocate(bra,ket)
    enddo
end subroutine get_nhat_elem
