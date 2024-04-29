subroutine cheby_Tn_rp0_rH_ket( bra, ket, cb_tmax, ket_f, &
                                dV, ng, havg, dh, rH_sub, Tn, ktop, &
                                polynom_type)

  implicit none

  integer     polynom_type
  integer     ng, ktop, k
  real*8      havg, dh,  dV
  real*8      bra(ng), ket(ng), ket_f(ng), Tn(0:ktop)
  real*8      pa(ng), pb(ng), pc(ng)
  complex*16  cb_tmax(0:ktop)
  external rH_sub

  write(6,*)' ktop in faber ',ktop

  ket_f =0

  pa = ket
  call rH_sub(pa, pb, ng)     
  pb = pb/dh-havg/dh*pa
  
  do k=0, ktop

     Tn(k) = sum(bra*pa)*dV; write(669,*)k,Tn(k); call flush(669)
     
     call rH_sub(pb, pc, ng); pc = pc/dh-havg/dh*pb

     ket_f = ket_f + cb_tmax(k) * pa

     select case(polynom_type)
        
     case(1)
        pc = 2.d0*pc - pa
     case(2)
        pc = 2.d0*pc + pa

     case default
        write(6,*)' polynom_type =  ',polynom_type
        stop
     end select

     pa = pb
     pb = pc
  end do

end subroutine cheby_Tn_rp0_rH_ket

subroutine cheby_Tn_rp0_rH(bra, ket, dV, ng, havg, dh, rH_sub, Tn, ktop, &
                           polynom_type)
  implicit none

  integer     polynom_type
  integer     ng, ktop, k
  real*8      havg, dh,  dV
  real*8      bra(ng), ket(ng), Tn(0:ktop)
  real*8      pa(ng), pb(ng), pc(ng)
  external rH_sub

  pa = ket
  call rH_sub(pa, pb, ng)     
  pb = pb/dh-havg/dh*pa
  
  do k=0, ktop

     Tn(k) = sum(bra*pa)*dV; write(669,*)k,Tn(k)
     
     call rH_sub(pb, pc, ng); pc = pc/dh-havg/dh*pb

     select case(polynom_type)
        
     case(1)
        pc = 2.d0*pc - pa
     case(2)
        pc = 2.d0*pc + pa
     case default
        write(6,*)' polynom_type =  ',polynom_type
        stop
     end select

     pa = pb
     pb = pc
  end do

end subroutine cheby_Tn_rp0_rH

subroutine c_bessel(cn, ktop, c_x)
  implicit none
  integer ktop, k
  complex*16 cn(0:ktop), c_x
  complex*16, parameter :: ci=(0.d0,1.d0)
  call c_bessel_fft_sub(cn, ktop, c_x) 
  
  cn(1:ktop) = cn(1:ktop)/2.d0
  do k=0,ktop
     cn(k) = cn(k)/(ci**k)
  enddo
end subroutine c_bessel

subroutine c_bessel_fft_sub(cn, ktop, c_x)  ! possibly multiply by (-1)**k
                                            ! or simply conjugate
  implicit none
  
  integer ktop, mx, nx, ix,st
  complex*16 c_x, cn(0:ktop)
  complex*16, allocatable :: ca(:), cra(:), cia(:)
  real*8 pi
  real*8,     allocatable :: xa(:), fa(:)
  complex*16  :: ci=(0.d0,1.d0)

  pi = dacos(-1.d0)

  
  mx = dlog(dble(ktop))/dlog(2.d0) + 3
  nx = 2**mx
  
  allocate(ca(0:nx-1), xa(0:nx-1), fa(0:nx-1), &
       cra(0:nx-1), cia(0:nx-1), &
       stat=st); call check(st,0,' caxa   ')
  
  do ix=0,nx-1
     fa(ix) = ix*2.d0*pi/dble(nx)
  end do
  xa = cos(fa)
  
  ca = exp(-ci*c_x*xa)
  cra = dble(ca)
  cia = aimag(ca)
  
  call fftsa(nx, cra, mx)
  call fftsa(nx, cia, mx)
  
  cra =  dble(cra)
  cia = -dble(cia)
  
  ca = (cra + ci*cia)*2.d0/nx
  
  ca(0) = ca(0)/2.d0
  
  cn = ca(0:ktop)
  
  deallocate(ca, xa, fa, cra, cia)
  
end subroutine c_bessel_fft_sub

subroutine get_faber_ce(havg, dh, ea, ne, aa, k, ce, polynom_type)
  implicit none

  integer    ne, k, polynom_type
  real*8     havg
  real*8     dh
  real*8     ea(ne)
  real*8     aa
  complex*16 ce(ne), ca(ne), cb(ne), cq(ne)
  complex*16, parameter :: ci=(0.d0,1.d0)

  ca(:)  = (ea(:)+ci*aa-Havg)/dh
  cb(:)  = sqrt(1.d0-ca(:)**2)
  cq(:)  = ca(:)-ci*cb(:)

  select case(polynom_type)
     
  case(1) ! chebyshev
     if(k>0) then
        ce(:)   = 2*cq(:)**dble(k)*(-ci/cb(:)/dh)
     else
        ce(:)   =                  (-ci/cb(:)/dh)
     end if

  case(2)
     if(k>0) then
        ce(:)   = 2*cq(:)**dble(k)*(-ci/cb(:)/dh)* (ci)**dble(k)
     else
        ce(:)   =                  (-ci/cb(:)/dh)
     end if

  case default
     write(6,*)' polynom_type is incorrect ',polynom_type
     stop
  end select

  write(7654,*)k,abs(ce(1)); call flush(7654)

end subroutine get_faber_ce

subroutine faber_from_rTn(  havg, dh, ea, ne, aa, Tn, ktop, cge, &
                            polynom_type)
  implicit none
  
  integer ne, ktop, k, polynom_type
  real*8     havg, dh
  real*8     ea(ne), aa
  real*8     Tn(0:ktop)
  complex*16 cge(ne), ce(ne)

  cge = 0.d0
  do k=0,ktop
     call get_faber_ce(havg, dh, ea, ne, aa, k, ce, polynom_type)
     cge(:) = cge(:) + Tn(k)*ce(:)
  end do
end

!program cheby_drvr_program
!call cheby_Tn_rp0_drvr
!end

subroutine cheby_Tn_rp0_drvr
  implicit none

  integer            :: ie
  integer, parameter :: ng=1, ne=10, ktop=5000
  external rh_simple_sub
  complex*16, parameter :: ci=(0.d0,1.d0)

  integer  k
  real*8   dh, havg
  real*8   ea(ne), aa, dV
  real*8   ket(ng), bra(ng), Tn(ktop)
  complex*16  cge(ne)
  
  ea(1)=-0.9d0; do ie=1,ne; ea(ie)=ea(1)+(ie-1.d0)*0.2d0; enddo
  
  aa = 0.1

  ket=1.d0
  bra = 1

  havg = 0.d0
  dh = 1d0
  dV = 1.d0

  call cheby_Tn_rp0_rH(bra, ket, dV, ng, havg, dh , &
       rH_simple_sub, Tn, ktop, 1)
  call faber_from_rTn(  havg, dh, ea, ne, aa, Tn, ktop, cge, 1)

  do  ie=1, ne
     write(6,*)' ie, diff ',real(ea(ie)), &
                        1.d0/(ea(ie)+ci*aa-0.9d0)-cge(ie) 
     write(6,*)' ie,  cge(ie)         ',real(ea(ie)), cge(ie)
     write(6,*)
  end do
   
end subroutine cheby_Tn_rp0_drvr

subroutine rh_simple_sub(pa, pb, ng)
  implicit none

  integer ng
  real*8 pa(ng), pb(ng)
  
  pb = 0.9d0*pa
end subroutine rh_simple_sub


subroutine faber_general_c( bra, ket, cb_tmax, ket_f, &
                                dV, ng, havg, dh, cH_sub, Tn, ktop)

  implicit none

  real*8, parameter :: xm = 0.d0, xd = 0.25d0

  integer     ng, ktop, k
  real*8      havg, dh,  dV
  complex*16  term
  complex*16  bra(ng), ket(ng), ket_f(ng), Tn(0:ktop)
  complex*16  pa(ng), pb(ng), pc(ng)
  complex*16  cb_tmax(0:ktop)
  external cH_sub

  write(6,*)' ktop in faber ',ktop
  ket_f = 0.d0
  pa    = ket
  call cH_sub(pa, pb, ng)     
  pb = pb/dh-havg/dh*pa
  pb = pb   -xm*     pa
  
  do k=0, ktop
     Tn(k) = sum(bra*pa)*dV; write(669,*)k,Tn(k); call flush(669)
     call cH_sub(pb, pc, ng); pc = pc/dh-havg/dh*pb

     ket_f = ket_f + cb_tmax(k) * pa
     term = xd
     if(k>0) term = 2*xd
     pc    = pc - xm * pb  - term* pa
     pa = pb
     pb = pc
  end do

end subroutine faber_general_c

subroutine faber_driver
  implicit none

  integer, parameter :: ng = 20, ktop=300
  real*8,  parameter :: dV = 1.d0

  real*8     :: havg = 0.d0, dh = 2.5d0
  complex*16 :: bra(ng), ket(ng), cb_tmax(0:ktop), ket_f(ng), Tn(0:ktop)
  external ch_simple_sub

  ket = 1.d0
  bra = 1.d0
  call  faber_general_c( bra, ket, cb_tmax, ket_f, &
                                dV, ng, havg, dh, cH_simple_sub, Tn, ktop)
end subroutine faber_driver

subroutine ch_simple_sub(ca, cb, ng)
  implicit none
  integer ig, ng
  complex*16 ca(ng)
  complex*16 cb(ng)
  real*8     xg
  complex*16 :: ci = (0.d0,1.d0)

  do ig=1,ng
     xg = dble(ig-ng/2)/(ng/2)
     cb(ig) = ca(ig)*(xg+ci*xg**3)
  end do
  
end subroutine ch_simple_sub
