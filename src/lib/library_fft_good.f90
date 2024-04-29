!program drive_fft_easy
!  call drive_fft_easy_mul()
!end

!
! This subroutine reads from a file data in a specfieid format (see below). 
! Note that if tgivenornot==2, dt is not used.
! The output is given in output files.

subroutine drv_fft_plot
  implicit none
  integer realorcomplex,tgivenornot,ifile,ifltr, ifout, ntin
  real*8  dt,gamma,strech

  write(6,*)' versatile fft routine '
  write(6,*)' Has 9 input parameters, as follows.  All values, if given 0, revert to default '
  write(6,*)' iin,iout -- # of input/output file ifileout default -- ifilin+10000 '
  write(6,*)' damping (data multiplied by exp(-(damping*(t-t0))**2/2d0) ) '
  write(6,*)' strech -- by how much do you want to strech the file. default -- 4 '
  write(6,*)' realorcomplex -- real(defualt)  unless realorcomplex=2 (complex) or 3 (t,x,y) ' 
  write(6,*)' tgivenornot -- 1 for given, 2 not; default 1 '
  write(6,*)' ifltr -- sin filter if 1; default 0 .  nt, dt -- if 0, extracted from input '
  write(6,*)
  write(6,*)' Example input: 600, 6001, 1e-3, 6*0        (i.e., strech=4,real,tgiven,nofilter )'
  write(6,*)' so give me iin, iout, damping,strech, realorcomplex,tgivenornot, ifltr,nt,dt '
  read(5,*)  ifile, ifout, gamma, strech,realorcomplex, tgivenornot, ifltr,ntin,dt
  call fft_plot(realorcomplex,tgivenornot,ifile,ifout, ifltr, ntin, dt, strech, gamma)
end subroutine drv_fft_plot


subroutine fft_plot(realorcomplex,tgivenornot,ifile,ifout, ifltr,ntin, dt, strech, gamma)
  implicit none

  integer i, ntot, realorcomplex, tgivenornot, m, ifile, nfull,ifltr,ifout, ntin
  real*8  dt , xterm,time,t0, kmin, dk, gamma, pi, strech, yterm
  complex*16, allocatable :: ca(:), caout(:)

  pi = dacos(-1.d0)

  if(strech<0.01) strech=4

  do i=1,1000000
     read(ifile,*,end=99)xterm
  enddo
  write(6,*)' too many lines '
  stop


99 ntot=i-1
  write(6,*)' ntot = ',ntot
  rewind(ifile)
  
  !  nfull  >= 2*ntot

  if(ntin/=0) then
     if(ntin<ntot) ntot=ntin
  endif

  nfull = ntot*nint(strech)
  do m=0,int(dlog(dble(nfull))/dlog(2.d0)) + 4
     if(2**m.le.nfull*2.and.2**(m+1)>nfull*2) goto 66
  enddo
  66 continue
  nfull = 2**m

  allocate(ca(0:ntot-1),   stat=i); call check(i,0,' st_i ')
  allocate(caout(0:nfull-1),stat=i); call check(i,0,' st_i ')
  ca = 0.d0
  caout = 0.d0

  do i=0,ntot-1
     select case(tgivenornot)
     case(0,1)
        select case(realorcomplex)
        case(0,1)
           read(ifile,*)time, xterm; ca(i) = xterm
        case(2)
           read(ifile,*)time, ca(i)
        case(3)
           read(ifile,*)time, xterm, yterm; ca(i)=dcmplx(xterm, yterm)
        case default
           write(6,*)' realorcomplex '
           write(6,*)  realorcomplex 
           stop
        end select

        if(i==0) t0=time
        if(i==1) dt=time-t0

        if(abs(time-t0-i*dt)>1.e-4*dt) then
           write(6,*)' problem in time : time, t0,i,dt ',time,t0,i,dt
           stop
        endif
     case(2)
        t0 = 0.d0
        select case(realorcomplex)
        case(1)
           read(ifile,*)xterm; ca(i) = xterm
        case(2)
           read(ifile,*)ca(i)
        case default
           write(6,*)' realorcomplex '
           write(6,*)  realorcomplex 
           stop
        end select

     case default
        write(6,*)' tgivenornot '
        stop
     end select
  end do

  write(6,*)' ntot,nfull ',ntot,nfull

  do i=0,ntot-1
     time = dt*i
     ca(i) = ca(i)* exp(-(gamma*dt*i)**2/2d0)
     if(ifltr>0) then
        ca(i)= ca(i)* sin(time*pi/(dt*ntot))**ifltr
     end if
  enddo

  call fft_easy_multiple(ca, caout, ntot, nfull, 1, dt, t0,kmin, dk)

  if(ifout==0) then
     ifout = 10000+ifile
     open(10000+ifile,file='fft_out.txt',status='unknown')
  end if
  if(ifout==ifile) then
     write(6,*)' adding 1000 to ifout so it doesnt equal ifile '
     ifout = ifile+1000
  end if

  do i=0,nfull-1
     write(ifout,88) i*dk+kmin, abs(caout(i)),dble(caout(i)),aimag(caout(i))
88   format(' ',4e19.9)
  enddo
end subroutine fft_plot
subroutine drive_fft_easy_mul()
  implicit none
  integer, parameter ::  nin= 57, nout=128, lot=2
  integer            ::  j, l, ilot
  real*8, parameter  ::  dx = 0.95, xmin = 0.1234
  real*8             ::  pi 
  real*8             ::  x
  real*8             ::  kmin
  real*8             ::  dk, kl
  complex*16,parameter :: ci = (0.d0, 1.d0)
  complex*16         :: caout(nout, lot), ca(nin, lot)
  complex*16         :: cterm

! random checking value
  do j=1, nin
     do ilot=1, lot
        ca(j, ilot) = j**2*0.3+ilot*j*0.5+ ci*j**2*ilot*0.9
     enddo
  enddo
  
  call fft_easy_multiple(ca, caout, nin, nout, lot, dx, xmin,kmin, dk)

  write(6,*)' kmin, dk ',kmin,dk
  
  do ilot=1, lot
     do l=1, nout
        kl = kmin+(l-1)*dk
        cterm= 0
        do j=1, nin
           x = xmin+(j-1)*dx
           cterm = cterm+ dx*exp(-ci*kl*x)*ca(j,ilot)
        enddo
        if(abs(caout(l,ilot)-cterm)>1.d-8) then
           write(6,*)
           write(6,*)l,ilot,' j, ilot, fft, ft'
           write(6,*)caout(l,ilot)
           write(6,*)cterm
           stop
        endif
     enddo
  enddo

end subroutine drive_fft_easy_mul
  
!
! input: ca, nin, xmin, nout (nout should be 2**power), dx, lot.  
! Output: caout, kmin, dk.  Orders correctly (negative k, then positive k)
!
! this routine returns:
!  caout(l,:) = dx*sum_j exp(-ci*xj*kl) ca(j, :)
!  where xj = xmin+j*dx (or (j-1)*dx if we count from 1)
!  and   kl = kmin+l*dk (or (l-1)*dk if we count from 1
!
!  Note that kmin and dk are output, no need to prespecify them;
!  also, typically, nout=nin;
!  Finally, if you fft a single vector, then give lot=1.

subroutine fft_easy_multiple(ca, caout, nin, nout, lot, dx, xmin,kmin, dk)
  implicit none
  
  integer                          ::nin, nout, lot, j, ilot, nbeg, nend,l
  integer is_it_power2
  external is_it_power2


  complex*16 , parameter           :: ci=(0.d0,1.d0) 
  complex*16                       :: ca(0:nin-1, lot), caout(0:nout-1, lot)
  complex*16, allocatable, dimension(:,:) :: ctemp
  real*8                           :: dx, kmin, dk, pi, kl,xmin
  
  allocate(ctemp(0:nout-1,lot), stat=j); call check(j,0,' ctempj ')

  pi = dacos(-1.d0)
  dk = 2*pi/nout/dx
  
  nbeg = -(nout)/2+1
  nend =  (nout)/2

!
! first we need to verify that nout is a product of 2.
!  

  if(is_it_power2(nout) /=1) then
     write(6,*)' problem with nout = ',nout,' not power of 2 '
  end if

  call check_le(nin, nout,' ninout ')
  caout = 0
  caout(0:nin-1, :) = ca
  
  call fft_multiple(caout, nout, lot)

  caout = caout*dx


!
! now we need to shift the array.
!

!  
 
  ctemp = caout

  caout(nend-1:nout-1, :) = ctemp(0:nend,   :)
  caout(0:nend-2 , :) = ctemp(nend+1:nout-1,:)
       
  kmin = -(nend-1)*dk

  do l=0,nout-1
     kl = kmin+dk*l
     caout(l, :) = caout(l, :) * exp(-ci*xmin*kl)
  enddo

  deallocate(ctemp)

end subroutine fft_easy_multiple

subroutine fft3(ca, Nx, Ny, Nz, Nv, sign)
  implicit none
  integer sign, Nx, Ny, Nz, Nv, iv
  complex*16    ca(Nx, Ny, Nz, Nv)

  do iv=1, Nv
     call fft_regular_3d(ca(1,1,1,iv),Nx,Ny,Nz,sign)
  enddo
end subroutine fft3

subroutine fft_regular_3d(ca,Nx,Ny,Nz,sign)
  integer sign, Nx, Ny, Nz, istop, st, iy, iz
  complex*16    ca(Nx, Ny, Nz)
  complex*16,   allocatable, dimension(:,:,:) :: cw1, cw2

  istop = 0
  if(abs(sign)/=1) istop =1
  if(Nx<1)         istop = 1
  if(Ny<1)         istop = 1
  if(Nz<1)         istop = 1
  if(Nx*Ny*Nz > 1e8) istop = 1

  if(sign == 1) ca = conjg(ca)

  if(istop == 1) then
     write(6,*)sign,Nx,Ny,Nz,' sign, Nx, Ny, Nz '
     stop
  endif

  if(Nx>1) call fft_multiple(ca,Nx,Ny*Nz)

  if(Ny>1) then
     allocate(cw1(Ny, Nx, Nz), stat=st); call check(st,0,'cw1st ')
     do iy=1, Ny
        cw1(iy,:,:) = ca(:,iy,:)
     enddo
     
     call fft_multiple(cw1,Ny,Nx*Nz)
     
     do iy=1, Ny
        ca(:,iy,:) = cw1(iy,:,:) 
     enddo
     
     deallocate(cw1)
  endif
  
  if(Nz>1) then
     allocate(cw2(Nz,Nx,Ny), stat=st); call check(st,0,' cw2 ')
     
     do iz=1, Nz
        cw2(iz,:,:) = ca(:,:,iz)
     enddo

     call fft_multiple(cw2,Nz,Nx*Ny)
     
     do iz=1, Nz
        ca(:, :, iz) = cw2(iz, :, :)
     enddo

     deallocate(cw2)
  endif

  if(sign == 1) ca = conjg(ca)

end subroutine


subroutine fft_multiple(ca, n, LOT)
  implicit none
  integer n, m, lot
  complex*16 ca(n, lot)

  m = nint(dlog(n*1.d0)/dlog(2.d0))
  call check(2**m, n, ' 2**m,n ')
  call check_le(1, lot, ' 1_lot ')
  call check_le(lot,101202,' lot_big ')

  call fftsa_lot(n, ca, m,  lot)
end subroutine fft_multiple

subroutine fft_good(ca, n)
  implicit none
  integer n, m
  complex*16 ca(n)

  m = nint(dlog(n*1.d0)/dlog(2.d0))
  call check(2**m, n, ' 2**m,n ')

  call fftsa(n, ca, m)  ! in this case
  
end subroutine fft_good

!***********************************************************************
!
!                      SUBROUTINE FFTSA
!
!     THIS IS A SIMPLE FFT ROUTINE FOR COMPLEX VECTORS OF 
!     LENGTH N=2**M.
!
!     ARGUEMENTS:  
!        N         LENGTH OF THE VECTOR
!        A         COMPLEX VECTOR OF LENGTH N. ON INPUT, A(k)=f[r(k)]
!                  ON OUTPUT, A(k)=sum(j=1,N) f[r(j)] * exp[i*r(j)*p(k)]  ! actually, exp(-i*r(j)*p(k))
!        M         N=2**M
!
!***********************************************************************
      SUBROUTINE FFTSA(N,A,M)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 A(*),U,W,T,ZI
      PIPI = DACOS(-1.D0)
      N11=2**M
      IF(N.NE.N11) THEN
         WRITE(6,*) ' *** FATAL ERROR IN ROUTINE FFTSA '
         WRITE(6,*) ' *** N NOT EQUAL 2**M: (N,M) - ', N, M
         STOP
      ENDIF
      NM1=N11-1
      NV2=N11/2
      J=1
      DO 7 I=1,NM1
      IF (I.GE.J) GOTO 5
      T=A(J)
      A(J) = A(I)
      A(I) = T
5     K=NV2
6     IF (K.GE.J) GOTO 7
      J =J- K
      K=K/2
      GOTO 6
7     J=J+K
      DO 20 L=1, M
      LE = 2** L
      LE1 = LE/2
      U=(1.D0,0.D0)
      ANG=PIPI / LE1
      ZI=(0.D0,1.D0)
      W= DCOS(ANG)+ZI* DSIN(ANG)
      DO 20 J=1, LE1
      DO 10 I=J,N11, LE
      IP=I + LE1
      T= A(IP)*CONJG(U)
      A(IP) = A(I) - T
10    A(I)  = A(I) + T
20    U=U*W
      RETURN
      END

!
!     Like fftsa for many vectors.  Did not test yet.
!

      SUBROUTINE FFTSA_lot(N,A,M,lot)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 A(n,lot)


      complex*16  w, zi
      complex*16, allocatable, dimension(:) :: U, TA
      allocate(U(lot),TA(lot), stat=ist)
      call check(ist,0,' istff')

      PIPI = DACOS(-1.D0)
      N11=2**M
      IF(N.NE.N11) THEN
         WRITE(6,*) ' *** FATAL ERROR IN ROUTINE FFTSA '
         WRITE(6,*) ' *** N NOT EQUAL 2**M: (N,M) - ', N, M
         STOP
      ENDIF
      NM1=N11-1
      NV2=N11/2
      J=1
      DO 7 I=1,NM1
      IF (I.GE.J) GOTO 5
      TA=A(J,:)
      A(J,:) = A(I,:)
      A(I,:) = TA
5     K=NV2
6     IF (K.GE.J) GOTO 7
      J =J- K
      K=K/2
      GOTO 6
7     J=J+K
      DO 20 L=1, M
      LE = 2** L
      LE1 = LE/2
      U=(1.D0,0.D0)
      ANG=PIPI / LE1
      ZI=(0.D0,1.D0)
      W= DCOS(ANG)+ZI* DSIN(ANG)   ! really -ZI*sin since used in conjg(U)
      DO 20 J=1, LE1
      DO 10 I=J,N11, LE
      IP=I + LE1
      TA= A(IP,:)*CONJG(U)
      A(IP,:) = A(I,:) - TA
10    A(I,:)  = A(I,:) + TA
20    U=U*W
      deallocate(U,TA)
      RETURN
      END
