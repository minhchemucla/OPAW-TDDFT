subroutine interpolate_vks_dij
  use main_mod
  use tddft_mod
  use atom_mod, only : at => atominfo, p => pawinfo,&
          natom, atom_map!, at_old => atominfo_old
  use mpi_lib_ours
  implicit none
  integer :: ia,it2,ms
  real*8  :: vks_tmp(nx,ny,nz)
  real*8, allocatable :: dij_tmp(:,:)

  if(.not.allocated(vks_old)) allocate(vks_old(nx,ny,nz))
  if(.not.allocated(vks_pert_old)) allocate(vks_pert_old(nx,ny,nz))

  if(it.eq.1) then
    vks_old = vks
    vks_pert_old = vks_pert
    do ia=1,natom
      at(ia)%dij_old(:,:) = at(ia)%dij(:,:)
      at(ia)%dij_pert_old(:,:) = at(ia)%dij_pert(:,:)
    enddo
  else 
     vks_tmp = vks
     vks = vks*1.5d0-vks_old*0.5d0
     vks_old = vks_tmp

     vks_tmp = vks_pert
     vks_pert = vks_pert*1.5d0-vks_pert_old*0.5d0
     vks_pert_old = vks_tmp
     do ia=1,natom
       it2=atom_map(ia)
       ms=p(it2)%mstates

       allocate(dij_tmp(ms,ms))

       dij_tmp(:,:) = at(ia)%dij(:,:)
       at(ia)%dij(:,:) = at(ia)%dij(:,:)*1.5d0-at(ia)%dij_old(:,:)*0.5d0
       at(ia)%dij_old(:,:) = dij_tmp(:,:)

       dij_tmp(:,:) = at(ia)%dij_pert(:,:)
       at(ia)%dij_pert(:,:) = at(ia)%dij_pert(:,:)*1.5d0-at(ia)%dij_pert_old(:,:)*0.5d0
       at(ia)%dij_pert_old(:,:) = dij_tmp(:,:)

       deallocate(dij_tmp)
     enddo
  endif

end subroutine
