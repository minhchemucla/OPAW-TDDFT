! first formalism: for Chebsyehv (T_n(H), H hermitian or close to it.)
!
! Next formalism: for T_tilde(L), where L=iH (i.e., e.g., L=d/dt), L close
! to antihermitian, and T_titlde_n+2 = 2L T_tilde_n+1 + T_tilde_n (i.e., + 
! sign on the last term rather than - sign)

! First FOrmalism is:
! Formalism: f(H) ~ sum_n a_n T_n(W) 
! where
!  W= H_scaled = (H-I*Hbar)/dH
!  T_n(W) = cos(n acos(W))
!  Define
!  Z=acos(W)
! Then
!  W=cos(Z), T_n(W) = cos(nz)
! Now a_n = 1/(pi)*(2-delta_n0) integral cos(n z)f(h(z)) dz  (0=<z<pi)
!(since integral(cos(nz)cos(mz)dz)= 1/2 *integral(cos(n-m)z+cos(n+m)z)dz=
!(if #z's bigger than 2N) delta_nm * pi /(2-delta_n0)
! so 
!
! f(H) ~ sum_n (2-delta_n0)   b_n
!
! where
!       b_n = integral cos(nz) f(h(z)) dz/pi
! Let's calculate therefore expectation values. 
!
! <f(H)> = sum_n (2-delta_n0) b_n c_n
! where c_n = <psi0|T_n(W)|psi0> or an ensemble average.
!
! Then, the question is: do we need to calculate a_n for each function f
! (e.g., if f_E(H) ~g(H-E) 
! Fortunately, not.
! Instead, calculate:
! <f> = sum_n (2-delta_n0)/2 (d_n + q_n) c_n
!
! where d_n = integral exp( inz) f(h(z)) dz/pi
!       q_n = integral exp(-inz) f(h(z)) dz/pi
!
! if f is an inherently real function (i.e., real for real arguments) then
! we can simply write
! <f_E> = Re(s(E)) 
! where s(E) = sum_n (2-delta_n0)/2  * d_n c_n
!
! Le'ts plug
!  s(E)= integral v(z) f(h(z)-e) dz/pi] 
!   where 
!   v(z) = sum_n (2-delta_n0)/2 c_n exp(inz) 
!        
! also, since dz=constant=pi/nz, S(E)=sum(v*f)/nz  
!
!  Example: f(h(z)) = delta_smooth(h(z)-e)
!
!  Pick 1000 values of E. 
!
! Note one more trick --
!  reduce the values of c_n exponentially (or through a Gaussian).
!
!   One more thing: v(z).  
!   The summation is from 0 to pi-dz. 
!   So let's change 0 to 2*pi-dz.  This is OK, since w is then just doubled;
!   so we have exp(inz)+exp(in(-z)) = cos(nz)/2.d0   

!
! Input file should be in chb_res.txt
! To call: just write 
!  call chebyshev_spec.  
!  end
!
!
! Next formalism: Note that H==-iL is close to hermitian.
!  So we can look at T_n(H). Now T_n = 2H T_n-1 - T_n-2.
! T_n = -2i L T_n-1 - T_n-2.  Define: T_tilde_n(L) = (i)^n T_n(H)
! so T_tilde_n(L) = (i)^n T_n = (i)^n (-2i T_n-1 - T_n-2)
! = (2i^(n-1) T_n-1 + (i^(n-2)) T_n-2 (as -i*i^n = i^(n-1), -i^n = i^(n-2))
! = 2 T_tilde_n-1 + T_tilde_n
!
! So if we are given T_tilde, we should just multiply them by
! (-i)^n to get T_n(H); it wil be simpler if the user prepares it.


subroutine chebyshev_spec
  implicit none
  integer, parameter :: nemax    = 100000
  integer  nmax, nch, nz, ie, ne
  real*8   egrd(nemax), ee, hbar, dh, ww_res
  complex*16 sE
  complex*16, allocatable :: c_chb(:), v_chb(:), zgrd(:), ff(:)
  
  call get_param  
  nz = nmax

  call allocate_c_v_z
  call get_chb
  call make_zgrid
  call convert_chb
  call make_egrid

  do ie=1,ne
     ee = egrd(ie)
     call make_func(ww_res)
     SE = sum(ff*v_chb)/nz
     write(1000,*)egrd(ie),dble(se)
  enddo

contains
  subroutine get_param
    implicit none
    integer ii,m, nch_read
    real*8 tempx
    open(900,file='chb_res.txt',status='old')
    
    read(900,*)hbar,dh
    write(6,*)' hbar, dh ',hbar,dh

    write(6,*)' give me ww_res , nch_read (give 0 if unknown) '
    read( 5,*)  ww_res, nch_read
    write(6,*)' = ',  ww_res,nch_read
    nmax=0
    do 
       read(900,*,end=99)ii,tempx 
    end do
    
99  continue
    rewind(900)
    if(nch_read<1.or.nch_read>ii+1) then
       nmax=ii+1
       nch =ii+1
    else
       nmax = nch_read
       nch  = nch_read
       ii   = nch_read-1
    end if
    m=int(dlog(dble(nmax)+0.001)/dlog(2.d0))
    write(6,*)' m ',m
    if(nmax/=2**m) then
       nmax=2**(m+2)
    else
       nmax=2**(m+1)
    endif
    if(nmax>4*ii.or.nmax<2*ii) then
       write(6,*)' nmax, ii ',nmax,ii
       stop
    endif
    nmax = nmax/2  ! erase

  end subroutine get_param

  subroutine allocate_c_v_z
    implicit none
    integer st
    allocate (c_chb(nmax), &
         v_chb(nz),   stat=st); call check(st,0,' st0 ')
    allocate (zgrd(nz), ff(nz),   stat=st); call check(st,0,' st0 ')
  end subroutine allocate_c_v_z

  subroutine get_chb
    implicit none
    integer i,im
    real*8  tempx
    read(900,*)
    c_chb = 0.d0
    do i=1,nch
       read(900,*)im,tempx 
       call check(im,i-1,' im,i-1 ')
       c_chb(i) = tempx
    end do

  end subroutine get_chb

  subroutine convert_chb
    implicit none
    v_chb = 0.d0
    v_chb(1:nch) = c_chb(1:nch)
    v_chb(2:nch) = v_chb(2:nch) *2.d0
    call fft_good(v_chb,nmax)
    v_chb = conjg(v_chb)  

  end subroutine convert_chb

  subroutine make_egrid
    implicit none
    real*8 emax,emin,de
    integer ie
    write(6,*)' give me emin, emax,ne '
    read( 5,*)  emin,emax,ne
    write(6,*)' = ',emin,emax,ne
    if(emin<hbar-dh) then
       write(6,*)' rescaling emin from ',emin,' to hbar-dh '
       emin = hbar-dh
       write(6,*)' new emin = ',emin
    end if
    if(emax>hbar+dh) then
       write(6,*)' rescaling emax from ',emax,' to hbar+dh '
       emax = hbar+dh
       write(6,*)' new emax = ',emax
    endif
    
    if(emin.ge.emax) then
       write(6,*)' reorder emin,emax ; now they are ',emin,emax
       stop
    endif

    call check_le(ne,nemax,' ne, nemax    ')

    de = (emax-emin)/max((ne-1),1)
    write(6,*)' de ',de
    do ie=1, ne
       egrd(ie) = emin+(ie-1)*de  
    enddo
  end subroutine make_egrid

  subroutine make_zgrid
    implicit none
    integer i
    real*8 pi
    pi = dacos(-1.d0)
    do i=1, nz
       zgrd(i) = 2*pi/nz*(i-1)
    end do
  end subroutine make_zgrid

  subroutine make_func(ww_res)
    implicit none
    integer i
    real*8 ww, h_of_z(nz), ww_res, pi
    pi = acos(-1.d0)
    h_of_z = hbar + dh* cos(zgrd)
    ww = ww_res* dh
    ff = exp(-(h_of_z-ee)**2/2.d0/ww**2.d0) / sqrt(2.d0*pi) / ww
  end subroutine make_func
end subroutine chebyshev_spec












