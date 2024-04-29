subroutine tddft_rhoij
  use main_mod 
  !use tddft_mod, only : phi_bar_tot, phi_bar_tot_pert
  use tddft_mod, only : phi_bar, phi_bar_pert, nvirtual, state_map
  use paw_mod
  use atom_mod, only : at => atominfo, p => pawinfo, atom_map, ngrid, natom
  use mpi_lib_ours
  implicit none

  integer :: ia,ig,ix,iy,iz,is,ms,it,js,ib,ik,jk,jb
  complex*16, allocatable :: ca(:),tmp(:,:)
  complex*16, allocatable :: ca_pert(:),tmp_pert(:,:)
  !ca = <p_i | psi >'s

  ik=1
  !do ik=1,nk_loc
      !jk=(ik-1)*nodes+rank+1
      !do ib=1,nocc
      do ia=1,natom
          it=atom_map(ia)
          ms=p(it)%mstates
          allocate(ca(ms),tmp(ms,ms),stat=stat)
          if(stat/=0) stop 'ca alloc problem in coeff'
          allocate(ca_pert(ms),tmp_pert(ms,ms),stat=stat)
          if(stat/=0) stop 'ca_pert alloc problem in coeff'

          at(ia)%rhoij=0d0
          at(ia)%rhoij_pert=0d0
          !if (rank .le. nodes-2) then
          do ib=1,nb
            jb=state_map(ib)
            if(jb <= nocc) then
              call proj(ia,phi_bar(:,:,:,ib,ik),ca,ms,ik) 
              call proj(ia,phi_bar_pert(:,:,:,ib,ik),ca_pert,ms,ik) 
          
              tmp=0d0
              tmp_pert=0d0
              do is=1,ms
                  do js=1,ms
                      tmp(is,js)=tmp(is,js)+conjg(ca(is))*ca(js)*2d0 !Note
                      tmp_pert(is,js)=tmp_pert(is,js)+conjg(ca_pert(is))*ca_pert(js)*2d0 !Note
                  enddo
              enddo
              at(ia)%rhoij=at(ia)%rhoij+dble(tmp)!*wk(jk)
              at(ia)%rhoij_pert=at(ia)%rhoij_pert+dble(tmp_pert)!*wk(jk)
            endif
          enddo
          !else
          !  do ib=1,nb-nvirtual
          !      call proj(ia,phi_bar(:,:,:,ib,ik),ca,ms,ik) 
          !      call proj(ia,phi_bar_pert(:,:,:,ib,ik),ca_pert,ms,ik) 
          !  
          !      tmp=0d0
          !      tmp_pert=0d0
          !      do is=1,ms
          !          do js=1,ms
          !              tmp(is,js)=tmp(is,js)+conjg(ca(is))*ca(js)*2d0 !Note
          !              tmp_pert(is,js)=tmp_pert(is,js)+conjg(ca_pert(is))*ca_pert(js)*2d0 !Note
          !          enddo
          !      enddo
          !      at(ia)%rhoij=at(ia)%rhoij+dble(tmp)!*wk(jk)
          !      at(ia)%rhoij_pert=at(ia)%rhoij_pert+dble(tmp_pert)!*wk(jk)
          !  enddo
          !endif
          deallocate(ca,tmp)
          deallocate(ca_pert,tmp_pert)
          call allsum_r8(at(ia)%rhoij,size(at(ia)%rhoij))
          call allsum_r8(at(ia)%rhoij_pert,size(at(ia)%rhoij_pert))
          !at(ia)%rhoij=at(ia)%rhoij+at(ia)%rho0ij
          !at(ia)%rhoij_pert=at(ia)%rhoij_pert+at(ia)%rho0ij
      enddo
  !enddo
end subroutine

subroutine tddft_rhoij_2(pert)
  use main_mod , only : nocc, nb
  use tddft_mod, only : phi_bar, phi_bar_pert, nvirtual, state_map
  use paw_mod
  use atom_mod, only : at => atominfo, p => pawinfo, atom_map, ngrid, natom
  use mpi_lib_ours
  implicit none
  integer :: ia,ig,ix,iy,iz,is,ms,it,js,ib,ik,jk,jb
  complex*16, allocatable :: ca(:),tmp(:,:)
  logical :: pert
  !ca = <p_i | psi >'s

  ik=1
  do ia=1,natom
    it=atom_map(ia)
    ms=p(it)%mstates
    allocate(ca(ms),tmp(ms,ms),stat=stat)
    if(stat/=0) stop 'ca alloc problem in coeff'

    if(.not.pert) at(ia)%rhoij=0d0
    if(pert)      at(ia)%rhoij_pert=0d0
    do ib=1,nb
      jb=state_map(ib)
      if(jb <= nocc) then
        if(.not.pert) call proj(ia,phi_bar(:,:,:,ib,ik),ca,ms,ik) 
        if(pert)      call proj(ia,phi_bar_pert(:,:,:,ib,ik),ca,ms,ik) 
        tmp=0d0
        do is=1,ms
          do js=1,ms
              tmp(is,js)=tmp(is,js)+conjg(ca(is))*ca(js)*2d0 !Note
          enddo
        enddo
        if(.not.pert) at(ia)%rhoij=at(ia)%rhoij+dble(tmp)!*wk(jk)
        if(pert)      at(ia)%rhoij_pert=at(ia)%rhoij_pert+dble(tmp)!*wk(jk)
      endif
    enddo
    deallocate(ca,tmp)
    if(.not.pert) call allsum_r8(at(ia)%rhoij,size(at(ia)%rhoij))
    if(pert)      call allsum_r8(at(ia)%rhoij_pert,size(at(ia)%rhoij_pert))
  enddo
end subroutine
