subroutine make_cvkb(ma,ch,nxb,nyb,nzb,nnr,nr,rr,cvkb)
  use main_mod,   only : dx,dy,dz,dv,periodic
  use atom_mod, only : pawinfo

  implicit none
  integer nxb,nyb,nzb,ch,nnr,st,nr,ix,iy,iz,ir,ma !, l
  real*8  rr(*), rv(3), vrlg, rscalar, rlg
  real*8  aa_here, charge_here, dr1
  real*8, allocatable, dimension(:):: vr, vr2, rrlog
  complex*16 cvkb( nxb, nyb, nzb )
  
  if(periodic) then
     aa_here = max(7d0/min(nxb*dx,nyb*dy,nzb*dz), 5d0/rr(nr))
     !write(6,*)' aa_here ',aa_here
  end if

  allocate(vr(nnr), stat=st); call check0(st,' vrrr ')
  call vloc_prep_atom(ma, ch, pawinfo(ma)%vloc, rr, vr, nr, nnr)

  call check_le(nr, nnr,' nr, nnr ')
  !
  ! OK, now we have vr on a long grid.  Now we need to put in on a big 3-d grid.
  !  First, fit to spline. Let's fit log(r), attention to r=0; dr/10. And i removed the log.
  !
  allocate(vr2(nr), rrlog(nr), stat=st); call check0(st,' vr2 ')

  dr1 = rr(2)-rr(1)
  do ir=1,nr
     rrlog(ir) = (max(rr(ir),0.1*dr1))
  enddo
  call SPLINE_dn(rrlog,vr,nr,vr2)
  charge_here = vr(nr)*rr(nr)

  do iz=1,nzb
     do iy=1,nyb
        do ix=1,nxb
           rv = (/ (ix-1-nxb/2)*dx, (iy-1-nyb/2)*dy, (iz-1-nzb/2)*dz /)
           rscalar = max(1d-10,sqrt(sum(rv**2)))
           rlg = (max(rscalar,0.1*dr1)) ! no log now -- DN

           if(rscalar<rr(1)) then;      vrlg = vr(1)
           elseif(rscalar>rr(nr)) then; vrlg = charge_here/rscalar
           else
              call SPLINT_dn(rrlog,vr,vr2,nr,rlg,vrlg) 
           end if
            
           if(periodic) vrlg = vrlg - charge_here * erf(aa_here*rscalar)/rscalar 
           cvkb(ix,iy,iz)    = vrlg
        end do
     end do
  end do

  ! now we have the 3d potential in r space. convert to k:
  call fft3d_forward_many( nxb,nyb, nzb,1,cvkb)
  cvkb = cvkb*dv
  call shift_cvkb
  if(periodic) call add_long_range_periodic
  call check_real(cvkb, size(cvkb))
  deallocate(vr,vr2,rrlog)
contains
  subroutine interp_lin(rrlog,vr,nr,rlg,vrlg)
    use main_mod, only : dv
    implicit none
    integer nr, i, j,k
    real*8 rrlog(nr), vr(nr), rlg, vrlg
    i=1
    j=nr
    do while(j>i+1)
       k=(i+j)/2
       if(rlg<rrlog(k)) then
          j=k
       else
          i=k
       endif
    end do
    vrlg = (vr(i)*(rrlog(j)-rlg)+vr(j)*(rlg-rrlog(i)))/(rrlog(j)-rrlog(i)) 
  end subroutine interp_lin

  subroutine shift_cvkb
    use main_mod, only : dx, dy, dz
    implicit none
    integer ikx,iky,ikz
    real*8  dkx,dky,dkz,kx,ky,kz,pi,r0(3)
    complex*16, parameter :: ci = (0d0,1d0)
    pi=dacos(-1d0)
    r0 = (/ -nxb/2d0*dx, -nyb/2d0*dy, -nzb/2d0*dz /)
    do ikz=1,nzb
       do iky=1,nyb
          do ikx=1,nxb
             dkx = 0d0 ; if(dx>1e-8) dkx = 2.d0*pi/(dble(Nxb)*dx)
             dky = 0d0 ; if(dy>1e-8) dky = 2.d0*pi/(dble(Nyb)*dy)
             dkz = 0d0 ; if(dz>1e-8) dkz = 2.d0*pi/(dble(Nzb)*dz)
             
             kx = dble(ikx-1)* dkx
             ky = dble(iky-1)* dky
             kz = dble(ikz-1)* dkz
             
             if(kx>pi/dx) kx = kx-2d0*pi/dx  
             if(ky>pi/dy) ky = ky-2d0*pi/dy
             if(kz>pi/dz) kz = kz-2d0*pi/dz
             
             cvkb(ikx,iky,ikz) = &
             cvkb(ikx,iky,ikz) * exp(-ci*( r0(1)*kx+r0(2)*ky+r0(3)*kz))
          enddo
       enddo
    end do
  end subroutine shift_cvkb

  subroutine add_long_range_periodic
    use main_mod, only : dx, dy, dz
    implicit none
    integer ikx,iky,ikz
    real*8  dkx,dky,dkz,kx,ky,kz,pi,r0(3),k2
    complex*16, parameter :: ci = (0d0,1d0)
    pi=dacos(-1d0)
    r0 = (/ -nxb/2d0*dx, -nyb/2d0*dy, -nzb/2d0*dz /)
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
             if(k2>1d-4) then
                cvkb(ikx,iky,ikz) = &
                cvkb(ikx,iky,ikz) + &
                                  charge_here * 4d0*pi/k2*exp(-k2/4d0/aa_here**2)
             else
                cvkb(ikx,iky,ikz) = 0d0
             endif
          enddo
       enddo
    enddo


  end subroutine add_long_range_periodic
end subroutine make_cvkb

