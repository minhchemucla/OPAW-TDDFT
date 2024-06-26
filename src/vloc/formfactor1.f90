subroutine add_formfactor(ma,na,nx,ny,nz)  ! nx, etc -- refer to "big".
  use main_mod,         only : dx,dy,dz, periodic
  use atom_mod,        only : atom_map, atom_Z, atominfo
  use form,         only : cvkb, cvn, crho_n
  use mpi_lib_ours, only : rank, nodes, gather_c16
  implicit none
  integer ma,na,nx,ny,nz
  integer ix,iy,iz,ia,islc,slcwdth,izb,izt,i
  integer, external :: is_it_power2
  real*8 dkx,dky,dkz
  real*8 kx,ky,kz
  real*8 pi
  real*8 rstart(3), zc
  complex*16, allocatable :: cx(:), cy(:), cform(:,:), cform_slc(:,:,:)
  complex*16,   parameter :: ci = (0d0,1d0)

  pi = dacos(-1d0)
  zc = atom_Z(ma)
  
  rstart = -0.5d0*(/ dx*nx,dy*ny,dz*nz /)

  dkx = 2.d0*pi/(nx*dx)
  dky = 2.d0*pi/(ny*dy)
  dkz = 2.d0*pi/(nz*dz)

  allocate(cx(nx),cy(ny),cform(nx,ny),stat=i); if(i/=0) stop ' cx,cy '

  if(rank==0)then
     allocate(cform_slc(nx,ny,max(nodes,1)),stat=i); if(i/=0) stop ' cform_slc '
  else
     allocate(cform_slc(1,1,1),stat=i); if(i/=0) stop ' cform_slc_b '
  end if

  if(rank==0.or.rank==nodes-1) then
     !write(6,*)' cnt ', cnt(:,:)
     write(6,*)' rstart ',rstart 
  end if

  slc : do islc = 1, (nz-1)/max(nodes,1)+1
     izb = (islc-1)*max(nodes,1) + 1
     izt = (islc-1)*max(nodes,1) + max(nodes,1)
     if(izb>nz) stop ' izb '
     if(izt>nz) izt=nz
     slcwdth = izt-izb+1
     iz = izb + rank
     if(rank==0.or.rank==nodes-1)then;&
          write(6,*)'rank, islc, iz, slcwdth ',rank,islc,iz,slcwdth ;call flush(6)
     endif
     iznz : if(iz.le.nz) then
        kz = (iz-1)* dkz
        if(kz>pi/dz) kz = kz-2*pi/dz
        !if(iz==nz/2+1) kz=0.d0   ! remove??
        cform(:,:) = 0d0
        ialoop : do ia=1, na
           maif : if( ma == atom_map(ia) ) then
              do ix=1,nx
                 kx = (ix-1)* dkx
                 if(kx>pi/dx) kx = kx-2*pi/dx
                 !if(ix==nx/2+1) kx=0.d0 ! remove??
                 cx(ix) = exp(-ci* (atominfo(ia)%coord(1)-rstart(1))*kx )
                 if(ix==nx/2+1) cx(ix)=0d0 ! remove??
              enddo

              do iy=1,ny
                 ky = (iy-1)* dky
                 if(ky>pi/dy) ky = ky-2*pi/dy
                 !if(iy==ny/2+1) ky=0.d0  ! remove, see above
                 cy(iy) = exp(-ci* ( (atominfo(ia)%coord(2)-rstart(2))*ky + (atominfo(ia)%coord(3)-rstart(3))*kz  ) )
                 if(iy==ny/2+1.or.iz==nz/2+1) cy(iy)=0d0 !removed
              enddo

              do iy=1,ny
                 do ix=1,nx
                    cform(ix,iy) = cform(ix,iy) + cx(ix)*cy(iy)
                 enddo
              enddo
           end if maif
        end do ialoop
     end if iznz
     
     call  gather_c16(cform, cform_slc, nx*ny, 0)
     if(rank==0) then
        cvn(:,:,izb:izt) = cvn(:,:,izb:izt) + cform_slc(:,:,1:slcwdth) * cvkb(:,:,izb:izt)
     end if
     
     if(periodic.and.rank==0) then
        crho_n(:,:,izb:izt) =  crho_n(:,:,izb:izt) + cform_slc(:,:,1:slcwdth) * zc 
     end if
  end do slc
  if(izt.lt.nz) stop ' izt<nz '

  deallocate(cx,cy,cform)
  if(allocated(cform_slc))deallocate(cform_slc)
end subroutine add_formfactor
