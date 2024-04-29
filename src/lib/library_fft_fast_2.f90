!
! efficient fft for multiple vectors.  a(nl,n,nr) is fft along the middle
! dimension (n) only. Nl and or nr can be any value, in particulalr 1 (each
! or both; method
! is efficent only when nl is big (say bigger than 20)
!
subroutine fftk(NL,N,nr,A)
  implicit none
  real*8 pipi
  integer n,m,n11,nm1,nv2,i,j,NL,k,le,le1,l,ip,ig,ir,nr
  COMPLEX*16 A(NL,n,nr), U, ta(NL,nr)
  complex*16  w, zi
  real*8 ang
  
  PIPI = DACOS(-1.D0)
  m = nint(dlog(dble(n))/dlog(2.d0))
  N11=2**M
  IF(N.NE.N11) THEN
     WRITE(6,*) ' *** FATAL ERROR IN ROUTINE fftk '
     WRITE(6,*) ' *** N NOT EQUAL 2**M: (N,M) - ', N, M
     STOP
  ENDIF
  NM1=N11-1
  NV2=N11/2

  J=1
  DO I=1,NM1
     IF (I.GE.J) GOTO 5
     do ir=1, nr
        do ig=1,NL
           TA(ig,ir)=A(ig,J,ir)
           A(ig,J,ir) = A(ig,I,ir)
           A(ig,I,ir) = TA(ig,ir)
        enddo
     enddo

5    K=NV2
6    IF (K.GE.J) GOTO 7
     J =J- K
     K=K/2
     GOTO 6
7    J=J+K
  enddo

  DO L=1, M
     LE = 2** L
     LE1 = LE/2
     U=(1.D0,0.D0)
     ANG=PIPI / LE1
     ZI=(0.D0,1.D0)
     W= DCOS(ANG)-ZI* DSIN(ANG)
     DO J=1, LE1
        DO I=J,N11, LE
           IP=I + LE1
           do ir=1,nr
              do ig=1, NL
                 TA(ig,  ir) = A(ig,IP,ir)*U
                 A(ig,IP,ir) = A(ig,I, ir) - TA(ig,ir)
                 A(ig,I, ir) = A(ig,I, ir) + TA(ig,ir)
              enddo
           enddo
        enddo
        U=U*W
     enddo
  enddo

end SUBROUTINE Fftk

SUBROUTINE Fftka(NL,N,nr,A)
  implicit none
  real*8 pipi
  integer n,m,n11,nm1,nv2,i,j,NL,k,le,le1,l,ip,ir,nr
  COMPLEX*16 A(NL,n,nr), U, ta(NL,nr)
  complex*16  w, zi
  real*8 ang
  
  PIPI = DACOS(-1.D0)
  m = nint(dlog(dble(n))/dlog(2.d0))
  N11=2**M
  IF(N.NE.N11) THEN
     WRITE(6,*) ' *** FATAL ERROR IN ROUTINE fftka '
     WRITE(6,*) ' *** N NOT EQUAL 2**M: (N,M) - ', N, M
     STOP
  ENDIF
  NM1=N11-1
  NV2=N11/2

  do ir=1, nr
     J=1
     DO I=1,NM1
        IF (I.GE.J) GOTO 5
        TA(:,ir)=A(:,J,ir)
        A(:,J,ir) = A(:,I,ir)
        A(:,I,ir) = TA(:,ir)
5       K=NV2
6      IF (K.GE.J) GOTO 7
       J =J- K
       K=K/2
       GOTO 6
7      J=J+K
    enddo

    DO L=1, M
       LE = 2** L
       LE1 = LE/2
       U=(1.D0,0.D0)
       ANG=PIPI / LE1
       ZI=(0.D0,1.D0)
       W= DCOS(ANG)-ZI* DSIN(ANG)
       DO J=1, LE1
          DO I=J,N11, LE
             IP=I + LE1
             TA(:,  ir) = A(:,IP,ir)*U
             A(:,IP,ir) = A(:,I, ir) - TA(:,ir)
             A(:,I, ir) = A(:,I, ir) + TA(:,ir)
          enddo
          U=U*W
       enddo
    enddo
 enddo
end SUBROUTINE Fftka

!program fft3_check_top_prog
!  call fft3_check_top
!end program fft3_check_top_prog

subroutine fftk_3d(ca,lot, Nx,Ny,Nz,sign)
  integer lot, sign, Nx, Ny, Nz, istop, st, iy, iz
  complex*16    ca(lot, Nx, Ny, Nz)

  istop = 0
  if(abs(sign)/=1) istop =1
  if(lot<1)        istop = 1
  if(Nx<1)         istop = 1
  if(Ny<1)         istop = 1
  if(Nz<1)         istop = 1
  if(Nx*Ny*Nz > 1e8) istop = 1

  if(istop == 1) then
     write(6,*)lot,sign,Nx,Ny,Nz,' lot,sign, Nx, Ny, Nz '
     stop
  endif
  
  if(sign == 1) ca = conjg(ca)

  if(Nx>1) call fftk(lot,Nx,Ny*Nz,ca)

  if(Ny>1) call fftk(lot*Nx,ny,Nz,ca)

  if(Nz>1) call fftk(lot*Nx*Ny,nz,1,ca)

  if(sign == 1) ca = conjg(ca)

end subroutine fftk_3d

subroutine fftk_3d_g1st(ca, Nx,Ny,Nz,lot,sign)
  integer lot, sign, Nx, Ny, Nz, istop, st, iy, iz
  complex*16    ca(Nx, Ny, Nz,lot)

  istop = 0
  if(abs(sign)/=1) istop =1
  if(lot<1)        istop = 1
  if(Nx<1)         istop = 1
  if(Ny<1)         istop = 1
  if(Nz<1)         istop = 1
  if(Nx*Ny*Nz > 1e8) istop = 1

  if(istop == 1) then
     write(6,*)lot,sign,Nx,Ny,Nz,' lot,sign, Nx, Ny, Nz '
     stop
  endif
  
  if(sign == 1) ca = conjg(ca)

  if(Nx>1) call fftk(1,Nx,Ny*Nz*lot,ca)

  if(Ny>1) call fftk(Nx,ny,Nz*lot,ca)

  if(Nz>1) call fftk(Nx*Ny,nz,lot,ca)

  if(sign == 1) ca = conjg(ca)

end subroutine fftk_3d_g1st

