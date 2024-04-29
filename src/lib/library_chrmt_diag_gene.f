!      call driver_c_diag_gen()  ! dummy subroutine;checks cdiag_gen
!      end

      subroutine ch_generalized(F, aC, S, eps, M, L, tol_s)
      implicit none

      integer, parameter :: flag_check = 1
      integer, parameter :: method     = 1
      integer M, i, L, ism
      integer j                         
      complex*16 F(M, M), aC(M, M), S(M, M)
      real*8     tol_s, diff, eps(M), tol_m

      real*8,  allocatable   ::s_evl(:) , E_diag(:,:)
      complex*16, allocatable,dimension(:,:)::X,F_rot,S_rot,
     x                                          Y,aCsm
      integer,    allocatable,dimension(:)  :: imap


      allocate(X(M, M),s_evl(M), stat=j);if(j/=0) stop
      allocate(imap( M),          stat=j);if(j/=0) stop

      if(M<1.or.M>4000) stop  ! change

      call c_check_hrmt(F, M)
      call c_check_hrmt(S, M)

      call ch_diag_normlz(S, X, s_evl, M)

      write(6,*)'s_evl ',s_evl

      tol_m = tol_s * maxval(abs(s_evl))
      L = 0
      do j=1, M
         if(s_evl(j) > tol_m ) L=L+1
      enddo

      if(L == 0) return

      allocate (Y(M, L), stat=j); if(j /=0) stop
      
      ism = 0
      do i=1, M
         if((s_evl(i)) > tol_m) then
            ism = ism + 1
            Y(:, ism) = X(:, i)/  sqrt(s_evl(i))
         endif
      enddo

      if(ism /= L) then
         write(6,*)' ism ',ism;stop
      endif

      deallocate(X)

      allocate(F_rot(L, L), aCsm(L, L), stat=j);if(j/=0) stop
      F_rot = matmul(matmul(conjg(transpose(Y)), F), Y)  
      F_rot = (F_rot + conjg(transpose(F_rot)) )/2.d0

      eps = 0.d0
      call ch_diag_normlz(F_rot, aCsm, eps, L)

      aC = 0.d0; do j=1,L;aC(:, j) = matmul(Y, aCsm(:, j));enddo
      deallocate(F_rot,  aCsm, Y)

c------
c now check - skip if wanted
c-----
      check_if_2 : if(flag_check == 1) then
         allocate(E_diag(M, M), stat=j);if(j/=0) stop
      
         E_diag = 0.d0 ; do i=1, L; E_diag(i, i)=eps(i); enddo !note: to L
         diff = sum(abs(matmul(F, aC) - matmul(matmul(S,aC),E_diag)))
         
         diff = diff / sum(abs(aC))/sum(abs(F))*M**2
         write(6,*)' scaled diff in gene_diag ',diff
         
ccc   call rcheck_small(diff, ' r_diag_gen ')
         
         deallocate(E_diag)
      endif check_if_2

 88   format(1x,4f15.6)

      deallocate(s_evl, imap)

      end



      subroutine driver_chrmt_diag_gen()  ! dummy subroutine;checks cdiag_gene
      implicit none
      integer M; parameter(M=4)
      integer i, j, L
      complex*16 A(M,M), Vc(M,M), U(M,M), A_dum(M, M), S(M, M)
      complex*16 B(M, M)
      complex*16, parameter  :: c_i=(0.d0, 1.d0)
      real*8      deltanew, w(M)
      real*8      tol_s, ID(M,M)

      S = 0
      A = 0
      ID = 0.d0
      do i=1,M
         ID(i,i) = 1.d0
      enddo

      S(1,1) = 1.d-2
      S(2,2) = 1.d-2
      S(3,3) = 3
      S(4,4) = 3

      A(1,1) = 0.5
      A(2,2) = 0.5
      A(3,3) = 0.5
      A(4,4) = 0.8


      do i=1, M
         do j=1, M
            deltanew=0
            if(i==j)deltanew=1
            B(i, j) = 
     x       i*5+j*5 + (i-j)**2+i*j+(i**2+j**2+4.0)**2
     x       + c_i* (i+j)
!            A(i, j) = i + j + 2*i**3 + 2*j**3 + (i+j+1+i*j)**2 
!     x                + c_i *i*j/10.d0   +0.1*deltanew*i
!            S(i, j) = deltanew + c_i* (i-j)/300.d0 
         enddo
      enddo

      do i=1,M
         do j=1,i-1
            B(:,i) = B(:,i)-B(:,j)*dot_product(B(:,j),B(:,i))
         enddo
         B(:,i) = B(:,i)/sqrt(sum(abs(B(:,i))**2))
      enddo

      write(6,*)' deviation of B from unitarity '
      write(6,*)sum(abs( matmul(conjg(transpose(B)),B)-id))
      write(6,*)
      

!      B = B/sum(abs(B))
!      A = A+conjg(transpose(A))
!      A = A / sum(abs(A))*M**2
      

      S = matmul(B, matmul(S, conjg(transpose(B))))  ! note
      A = matmul(B, matmul(A, conjg(transpose(B))))  ! note



      write(6,*)' A '
      write(6,88)A

      write(6,*)' S '
      write(6,88)S

      tol_s = 1.d-8
      write(6,*)' tol_s ',tol_s
      call ch_generalized(A, Vc, S, w, M, L, tol_s)
      write(6,*)' after '
      write(6,*)' L ', L

      write(6,*)' A '
      write(6,88) A
      write(6,*)' S '
      write(6,88) S
      write(6,*)' w '
      write(6,88) w
      write(6,*)' Vc '
      do i=1, M
      write(6,88) Vc(i,:) 
      enddo

 88   format(1x,4d15.6)
      
      end

      subroutine c_check_hrmt(S, M)
      implicit none
      integer M
      complex*16 S(M,M)
      real*8 dev
      
      dev = sum(abs(S-conjg(transpose(S))))/sum(abs(S))
      if(dev>1.d-8) then
         write(6,*)' deviation from hermiticity is ',dev
         stop
      endif
      end

      subroutine ch_diag_normlz(F, C, eps, L)
      implicit none
      
      integer L,matz,ierr
      
      complex*16 F(L,L),C(L,L)
      complex*16, parameter ::ci=(0.d0,1.d0)
      real*8  eps(L)

      real*8, allocatable :: zr(:,:),zi(:,:),fv1(:),fv2(:),fm1(:,:)
      real*8, allocatable :: ar(:,:),ai(:,:)

      allocate(ar(L,L),ai(L,L),zr(L,L),zi(L,L),fv1(L),fv2(L),
     x         fm1(2,L),stat=ierr); call check(ierr,0,' ierr0     ')
      
      matz=1
      ar = F
      ai = 0.5d0/ci*(F-conjg(F))
      call ch(L,L,ar,ai,eps,matz,zr,zi,fv1,fv2,fm1,ierr)
      if(ierr/=0) then
         write(6,*)' ierr '
         stop
      endif
      C = zr+ci*zi

      end

