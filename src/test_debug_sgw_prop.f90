subroutine test_debug_sgw_prop
  use tddft_mod, only : dt
  use main_mod, only : nx,ny,nz,nn,vks,dens
  use atom_mod, only : at=>atominfo
  use mpi_lib_ours, only : rank, sync_mpi
  implicit none
  complex*16, allocatable :: pt(:,:)
  integer :: natom,ia,is,ns,st
  

  if(rank==0) then
    open(unit=302,file='debug_prop.bin',form='unformatted',status='old')
    read(302) dt; write(6,*) 'dt', dt
    read(302) ns; write(6,*) 'ns', ns
    call flush(6)
    allocate(pt(nn,ns),stat=st); if(st/=0) stop 'pt debug_prop'
    do is=1,ns
      read(302) pt(:,is)
    enddo
    read(302) vks
    read(302) dens
    read(302) natom
    do ia=1,natom
      read(302) at(ia)%dij
    enddo
    close(302)

    do is=1,ns
      call rk4_prop_shs_debug(is, pt(:,is))
    enddo

    call flush(6)
  endif

  call sync_mpi

  stop "developing test_debug_sgw_prop"

end subroutine
