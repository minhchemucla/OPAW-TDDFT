!      call driver_c_diag_gen()  ! dummy subroutine;checks cdiag_gen
!      end

c-----
c ! solves F aC = S_i aC E_diag. Wasteful on space.
c----
      subroutine c_diag_gnrlz_Shrmt_Fnot(F, aC, S, eps, M, L, tol_s)
      implicit none

      integer, parameter :: flag_check = 0
      integer, parameter :: method     = 1
      integer M, i, L, ism
      integer j                         
      complex*16 F(M, M), aC(M, M), S(M, M), eps(M)
      real*8     tol_s, diff, tol_m


      complex*16, allocatable,dimension(:,:)::X,F_rot,S_rot,E_diag,
     x                                          Y,aCsm
      integer,    allocatable,dimension(:)  :: imap
      real*8,     allocatable,dimension(:)  :: abs_s, srl_evl

      allocate(X(M, M),srl_evl(M), stat=j);if(j/=0) stop
      allocate(imap( M),          stat=j);if(j/=0) stop
      allocate(abs_s(M),          stat=j);if(j/=0) stop

      if(M<1.or.M>1000) stop  ! change

!!!!!!!!      call c_check_symm(F, M)  ! F not necc. symm!

      call ch_diag_normlz(S, X, srl_evl, M)
      tol_m = tol_s *maxval(abs(srl_evl))
      write(6,*)' srl_evl, tol_m ',srl_evl, tol_m
      L = 0
      do j=1, M
         if(abs(srl_evl(j)) > tol_m) L=L+1
      enddo

      if(L == 0) return

      allocate (Y(M, L), stat=j); if(j /=0) stop
      
      ism = 0
      do i=1, M
         if(abs(srl_evl(i)) > tol_m) then
            ism = ism + 1
            write(66,*)' ism, i ',ism,i
            Y(:, ism) = X(:, i)/  sqrt(abs(srl_evl(i)))
         endif
      enddo

      if(ism /= L) then
         write(6,*)' ism ',ism;stop
      endif

      deallocate(X)
      allocate(F_rot(L, L), aCsm(L, L), S_rot(L, L), stat=j)
      if(j/=0) stop
      F_rot = matmul(matmul(conjg(transpose(Y)), F), Y)  
      S_rot = matmul(matmul(conjg(transpose(Y)), S), Y)  

      eps = 0.d0
      if(sum(abs(F_rot - transpose(F_rot)))<1.d-8) then
         F_rot = (F_rot+transpose(F_rot))/2.d0
         write(66,*)' s_rot ',s_rot
           write(66,*)' F_rot,',F_rot
         call cs_diag_normlz(F_rot, aCsm, eps, L, 1, method)
           write(6,*)' F_rot,',F_rot
           write(6,*)' acsm  ',acsm
           write(6,*)' L,eps   ',L,eps

      else

       
         call cg_diag(F_rot, acSm, eps, L, 1, 0)
           write(6,*)' F_rot,',F_rot
           write(6,*)' acsm  ',acsm
           write(6,*)' L,eps   ',L,eps
         
      endif
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

      deallocate(srl_evl, imap, abs_s)

      end
c-----
c ! solves F aC = S_i aC E_diag. Wasteful on space.
c----
      subroutine c_diag_generalized(F, aC, S, eps, M, L, tol_s)
      implicit none

      integer, parameter :: flag_check = 0
      integer, parameter :: method     = 1
      integer M, i, L, ism
      integer j                         
      complex*16 F(M, M), aC(M, M), S(M, M), eps(M)
      real*8     tol_s, diff, tol_m

      complex*16, allocatable,dimension(:)   ::s_evl 
      complex*16, allocatable,dimension(:,:)::X,F_rot,S_rot,E_diag,
     x                                          Y,aCsm

      allocate(X(M, M),s_evl(M), stat=j);if(j/=0) stop

      if(M<1.or.M>1000) stop  ! change

!!!!!!!!      call c_check_symm(F, M)  ! F not necc. symm!
      call c_check_symm(S, M)


      call cs_diag_normlz(S, X, s_evl, M, 1, method)

      tol_m = tol_s *maxval(abs(s_evl))
      L = 0
      do j=1, M
         if(abs(s_evl(j)) > tol_m ) L=L+1
      enddo

      if(L == 0) return

      allocate (Y(M, L), stat=j); if(j /=0) stop
      
      ism = 0
      do i=1, M
         if(abs(s_evl(i)) > tol_m) then
            ism = ism + 1
            Y(:, ism) = X(:, i)/  sqrt(s_evl(i))
         endif
      enddo

      if(ism /= L) then
         write(6,*)' ism ',ism;stop
      endif

      deallocate(X)
c----
c erase
c----

      check_if : if(flag_check == 1) then
         allocate(S_rot(L, L),  stat=j);if(j/=0) stop
         S_rot = matmul(matmul(transpose(Y), S), Y)
         S_rot = (S_rot + transpose(S_rot))/2.d0
         
         do j=1,L; S_rot(j, j) = S_rot(j, j) - 1.d0; enddo
         diff = sum(abs(S_rot))
            
         if(diff>1.d-4) then    ! erase
            write(6,*)' diff in s_rot ', diff
            write(6,*) ' S_rot '
            write(6,88) S_rot
            write(6,*)' Y '
            write(6,88) Y
            stop
         endif
         deallocate(S_rot)

      endif check_if
c---
c end of erase
c---
      allocate(F_rot(L, L), aCsm(L, L), stat=j);if(j/=0) stop
      F_rot = matmul(matmul(transpose(Y), F), Y)  

      eps = 0.d0
      if(sum(abs(F_rot - transpose(F_rot)))<1.d-8) then
         F_rot = (F_rot+transpose(F_rot))/2.d0

         call cs_diag_normlz(F_rot, aCsm, eps, L, 1, method)
      else
         
         call cg_diag(F_rot, acSm, eps, L, 1, 0)
         
      endif
      
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

      deallocate(s_evl)

      end


      subroutine driver_c_diag_gen()  ! dummy subroutine;checks cdiag_gene
      implicit none
      integer M; parameter(M=4)
      integer i, j, L
      complex*16 A(M,M), Vc(M,M), w(M), U(M,M), A_dum(M, M), S(M, M)
      complex*16 B(M, M)
      complex*16, parameter  :: c_i=(0.d0, 1.d0)
      real*8      delta; external delta
      real*8      tol_s

      do i=1, M
         do j=1, M
            B(i, j) = 
     x       i*5+j*5 + (i-j)**2+i*j+(i**2+j**2+4.0)**2
     x       + c_i* (i+j)
            A(i, j) = i + j + 2*i**3 + 2*j**3 + (i+j+1+i*j)**2 
     x                + c_i *i*j/10.d0   +0.1*delta(i, j)* i
            S(i, j) = delta(i, j) + c_i* (i*j)/300.d0 
         enddo
      enddo

      B = B/sum(abs(B))
      A = A+transpose(A)
      A = A / sum(abs(A))*M**2

!      S = matmul(B, matmul(S, transpose(B)))  ! note
!      A = matmul(B, matmul(A, transpose(B)))  ! note

      write(6,*)' A '
      write(6,88)A

      write(6,*)' S '
      write(6,88)S

      tol_s = 1.d-8
      write(6,*)' tol_s ',tol_s
      call cs_diag_normlz(A, Vc, w, M, 1, 1)
      write(6,*)' w_pre '
      write(6,88) w
      call c_diag_generalized(A, Vc, S, w, M, L, tol_s)
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



