subroutine vloc_tuma_prep
  use main_mod,         only : nn, nx, ny, nz, periodic
  use main_mod,         only : vloc_tot
  use mpi_lib_ours,     only : rank, sync_mpi, bcast_r8
  implicit none
  integer :: n
  n=nn
  call vloc_prep
!  call vloc_cnst
!  vloc_tot = vloc_tot + vloc_c
!  call add_bare_charges_qq(vloc_tot, nx, ny, nz)
  call bcast_r8(vloc_tot, n, 0)
end subroutine vloc_tuma_prep

subroutine vloc_prep
  use main_mod
  use atom_mod,        only : atom_Z, natom, ntype, pawinfo
  use form,         only : cvkb, cvn, crho_n 
  use mpi_lib_ours, only : rank
  implicit none
  integer                  :: ch,ma,nrr,ia,ig,iz,iy,ix,nxb,nyb,nzb,st
  integer, external        :: is_it_power2
  complex*16, parameter    :: ci = (0.d0,1.d0)

  integer :: n ,it
  n=nn

  nxb = nx * scale_vh
  nyb = ny * scale_vh
  nzb = nz * scale_vh
  
  if(rank==0) then
     allocate(cvkb(nxb,nyb,nzb),cvn(nxb,nyb,nzb), stat=st)
     if(st/=0) stop ' cvkb cvn  allocation problem '
     cvn = 0d0
     if(periodic) then
        allocate(crho_n(nx,ny,nz), stat=st)
        call check0(st,' crho_n ')
        crho_n = 0d0
     endif
  endif

  atomtypeloop : do ma=1,ntype
     ch = atom_Z(ma)
     if(rank==0) then
        write(6,*)' kb pseudpt formfactor. atomtype(ma) =',ma,' charge= ',ch,' nrr ',pawinfo(ma)%nr
        call flush(6)
     endif
     
     !nrr = pawinfo(ma)%nr                                    !minh
     !if(nrr<50.or.nrr>3000) stop ' nr values problem '       !minh

     if(rank==0) then
        nrr = pawinfo(ma)%nr                                    !minh
        if(nrr<50.or.nrr>3000) stop ' nr values problem '       !minh
        call make_cvkb(ma, ch, nxb, nyb, nzb, nrr, nrr, pawinfo(ma)%rr,cvkb)
        call check_real(cvkb, size(cvkb))
        call plot_cvkb
     end if
     call add_formfactor(ma,natom,nxb,nyb,nzb)  ! cvn = cvn + form*wg_big*vk_big, for rank=0
  enddo atomtypeloop
  
  if(allocated(cvkb))deallocate(cvkb)
!  allocate(vloc_tot(n),stat=st)
!  if(st/=0) stop' vloc_tot alloc '

  vloc_tot = 0d0
  if(rank==0) then
     call fft3d_general(  cvn,nxb,nyb,nzb,1)
     cvn = cvn/(dble(nxb)*dble(nyb)*dble(nzb)*dv)
     call check_c_real(   cvn,size(cvn),' cvn ')     ! after fft v should be real
     call c_to_v_compress(cvn, nxb,nyb,nzb,vloc_tot,nx,ny,nz)
     write(6,*)' sum_vloc_tot(post kb) ',sum(abs(vloc_tot))
  end if
  
  if(allocated(cvn))deallocate(cvn)
contains
  subroutine plot_cvkb
    implicit none
    integer ikx,iky,ikz
    real*8  dkx,dky,dkz,kx,ky,kz,pi,k2
    pi=dacos(-1d0)
    write(6,*)' maxval, minval cvkb ',maxval(dble(cvkb)),minval(dble(cvkb))
    do ikz=1,nzb
       do iky=1,nyb
          do ikx=1,nxb
             dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(nxb)*dx)
             dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(nyb)*dy)
             dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(nzb)*dz)
           
             kx = dble(ikx-1)* dkx
             ky = dble(iky-1)* dky
             kz = dble(ikz-1)* dkz
             
             if(kx>pi/dx) kx = kx-2d0*pi/dx  
             if(ky>pi/dy) ky = ky-2d0*pi/dy
             if(kz>pi/dz) kz = kz-2d0*pi/dz
             
             k2 = kx**2+ky**2+kz**2
          enddo
       enddo
    enddo
  end subroutine plot_cvkb

end subroutine vloc_prep

subroutine  c_to_v_compress(cv, nxb, nyb, nzb, vv, nx, ny, nz)
  implicit none
  integer nxb, nyb, nzb
  integer nx, ny, nz
  integer mx, my, mz
  integer px, py, pz
  complex*16 cv(nxb, nyb, nzb)
  real*8     vv(nx,    ny,    nz)
  
  mx= 1 + nxb/2 - nx/2
  px=     nxb/2 + nx/2
  
  my= 1 + nyb/2 - ny/2
  py=     nyb/2 + ny/2    ! dont remove,mistake
  
  mz= 1 + nzb/2 - nz/2
  pz=     nzb/2 + nz/2

  vv = cv(mx:px,my:py,mz:pz)
end subroutine c_to_v_compress

