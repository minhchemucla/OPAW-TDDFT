subroutine vk_prep
  use main_mod
  use mpi_lib_ours
  implicit none
  integer nxb, nyb, nzb, stat

  nxb = scale_vh * nx
  nyb = scale_vh * ny
  nzb = scale_vh * nz
    
  if(rank==0) then
      call vk_prep_b(nxb,nyb,nzb)
  else
      if(allocated(vk)) deallocate(vk)
      allocate(vk(nxb*nyb*nzb),stat=stat)
      if(stat/=0) stop 'vk alloc problem'
  endif
  call sync_mpi
  call bcast_r8(vk,size(vk),0)

end subroutine vk_prep

subroutine vk_prep_b(nxb,nyb,nzb)
  use main_mod
  implicit none
  integer nxb, nyb, nzb, ngb,nd,ix,iy,iz,i
  real*8, allocatable:: grdb(:,:,:,:)

  nd = 3
  ngb = nxb*nyb*nzb; if(abs(ngb-dble(nxb)*dble(nyb)*dble(nzb))>1d-4) stop ' ngb not proper ; maybe need integer*8 '

  if(allocated(vk)) deallocate(vk);  allocate(vk(ngb), stat=i); if(i/=0) stop ' vk_prep_b '

  allocate(grdb(nxb,nyb,nzb,3), stat=i); if(i/=0) stop ' grdb '

  do iz=1,nzb
     do iy=1,nyb
        do ix=1,nxb
           grdb(ix,iy,iz,:) = (/ (ix-1-nxb/2)*dx, (iy-1-nyb/2)*dy, (iz-1-nzb/2)*dz /); 
        enddo
     enddo
  enddo
  
  !call prep_vk_TuMa_2(ngb, nd, nxb, nyb, nzb, dx,dy,dz, grdb, vk)
  call  prep_vk( nxb,nyb,nzb,dx,dy,dz,vk,periodic)
  vk= vk/dble(nxb)/dble(nyb)/dble(nzb)
  deallocate(grdb)
end subroutine vk_prep_b

