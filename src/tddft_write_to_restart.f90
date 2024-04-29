subroutine write_to_restart
  !writes phi and phi_bar of each rank to the respective wavefunction
  use tddft_mod
  use main_mod
  use mpi_lib_ours
  implicit none
  integer :: nds, ip,ip2,i,j
  complex*16 :: tmp(nx,ny,nz)
  character(100) :: file_phi_bar, file_phi_bar_pert,x
  logical fex
  

  nds = max(nodes,1)
  ip  = 3000+rank
  ip2  = 4000+rank


  close(305); rewind(305)
  if(rank==0) then
    inquire(file='./tddft_save',exist=fex)
    if(.not.fex) then
      call system("mkdir tddft_save") 
      open(unit=305,file='./tddft_save/info',status='new')
    else
      open(unit=305,file='./tddft_save/info',status='replace')
    end if

    write(305,*) nb, nds, it
    close(305)

    !file_phi_bar = './tddft_save/phi_bar.txt'
    !file_phi_bar_pert = './tddft_save/phi_bar_pert.txt'
    file_phi_bar = './tddft_save/phi_bar.bin'
    file_phi_bar_pert = './tddft_save/phi_bar_pert.bin'

    !write(6,*) 'rank, file_phi', rank, file_phi_bar
    !write(6,*) 'rank, file_phi_bar', rank, file_phi_bar_pert

    close(ip); rewind(ip)
    close(ip2); rewind(ip2)
    open(ip,file=file_phi_bar,status='replace',form='unformatted')
    open(ip2,file=file_phi_bar_pert,status='replace',form='unformatted')
    !open(ip,file=file_phi_bar,status='unknown')
    !open(ip2,file=file_phi_bar_pert,status='unknown')
  endif

  if (rank==0) then
    do i=1,nb
      write(ip ) phi_bar(:,:,:,i,:)
      write(ip2) phi_bar_pert(:,:,:,i,:)
    enddo
  endif
  do i=1,nds-1
    do j=1,nb
      tmp = phi_bar(:,:,:,j,1)
      call sr_mpi_c16(tmp, size(tmp), i, 0)
      if(rank==0) write(ip ) tmp

      call sync_mpi

      tmp = phi_bar_pert(:,:,:,j,1)
      call sr_mpi_c16(tmp, size(tmp), i, 0)
      if(rank==0) write(ip2) tmp
    enddo
  enddo

  !stop 'debug write to restart'
end subroutine write_to_restart
