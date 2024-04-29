subroutine read_phi_bar
  use tddft_mod
  use main_mod
  use mpi_lib_ours
  implicit none
  integer :: nds, nds_read, ip,ip2,i,j,is,st
  character(100) :: file_phi_bar, file_phi_bar_pert,x

  close(34)
  open(unit=34,file='./tddft_save/info',status='old')
  read(34,*) nb, nds_read, it_restart

  nds = max(nodes,1)
  if(nds_read.ne.nds) stop 'nds_read/=nds'

  call alloc_tddft

  allocate(state_map(nb),stat=st); if(st/=0) stop 'state_map'
  do is=1,nb
    state_map(is) = is + rank*nb 
  enddo

  file_phi_bar = './tddft_save/phi_bar.bin'
  file_phi_bar_pert = './tddft_save/phi_bar_pert.bin'
  !file_phi_bar = './tddft_save/phi_bar.txt'
  !file_phi_bar_pert = './tddft_save/phi_bar_pert.txt'

  ip = 234
  ip2 = 341

  close(ip )
  close(ip2)
  open(ip ,file=file_phi_bar     ,status='old',form='unformatted'); rewind(ip )
  open(ip2,file=file_phi_bar_pert,status='old',form='unformatted'); rewind(ip2)
  !open(ip ,file=file_phi_bar     ,status='old'); rewind(ip )
  !open(ip2,file=file_phi_bar_pert,status='old'); rewind(ip2)

  call sync_mpi
  do i=0,nds-1
    do j=1,nb
      if(rank==i) then
        read(ip ) phi_bar(:,:,:,j,1)
        read(ip2) phi_bar_pert(:,:,:,j,1)
      else
        read(ip ) 
        read(ip2) 
      endif
    enddo
  enddo

  call sync_mpi
end subroutine read_phi_bar
