subroutine vh_sub(dens,vh,scale_vh)
  use main_mod, only : nx,ny,nz,vk,nn
  implicit none
  integer  nxb,nyb,nzb,i,n
  real*8 dens(nn)
  real*8   vh(nn)
  real*8, allocatable :: densb(:)
  real*8, allocatable ::   vhb(:)
  integer  :: scale_vh

  nxb = scale_vh*nx
  nyb = scale_vh*ny
  nzb = scale_vh*nz

  n=nn

  call alloc_densb_vhb
  call dens_to_densb(dens,densb)
  call get_vcon(densb,nxb,nyb,nzb,vk,vhb,0)
  call vhb_to_vh(vhb,vh)
  call dealloc_densb_vhb
contains
  subroutine alloc_densb_vhb
    implicit none
    integer nb
    nb = n*(scale_vh)**3
    allocate(densb(nb), stat=i); if(i/=0) stop ' densb '
    allocate(vhb(nb),   stat=i); if(i/=0) stop ' vhb '
  end subroutine alloc_densb_vhb

  subroutine dealloc_densb_vhb
    implicit none
    deallocate(densb, vhb)
  end subroutine dealloc_densb_vhb

  subroutine dens_to_densb(dens,densb)
    implicit none
    real*8 dens( nx, ny, nz)
    real*8 densb(nxb,nyb,nzb)

    integer mx,my,mz,px,py,pz
    densb = 0d0
    
    mx= 1 + nxb/2 - nx/2
    px=     mx+nx-1

    my= 1 + nyb/2 - ny/2
    py=     my+ny-1

    mz= 1 + nzb/2 - nz/2
    pz=     mz+nz-1

    densb(mx:px,my:py,mz:pz) = dens
  end subroutine dens_to_densb

  subroutine vhb_to_vh(vhb,vh)
    implicit none
    real*8 vh( nx, ny, nz)
    real*8 vhb(nxb,nyb,nzb)
    integer mx,my,mz,px,py,pz

    mx= 1 + nxb/2 - nx/2
    px=     mx+nx-1

    my= 1 + nyb/2 - ny/2
    py=     my+ny-1

    mz= 1 + nzb/2 - nz/2
    pz=     mz+nz-1

    vh = vhb(mx:px,my:py,mz:pz) 

  end subroutine vhb_to_vh

end subroutine vh_sub

