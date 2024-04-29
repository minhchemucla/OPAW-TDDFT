      subroutine driver_diag_gen(M) ! dummy subroutine;checks rdiag_gene
      implicit none
      integer M
      integer i, j, L
      real*8 A(M,M), Vc(M,M), w(M), U(M,M), A_dum(M, M), S(M, M),delta
      real*8 B(M, M)

      do i=1, M
         do j=1, M
            B(i, j) = 
     x       i*5+j*5 + (i-j)**2+i*j+(i**2+j**2+4.0)**2
            A(i, j) = i + j + 2*i**3 + 2*j**3 + (i+j+1+i*j)**2
            S(i, j) = delta(i, j) 
         enddo
      enddo

      B = B/sum(abs(B))
      A = A+transpose(A)
      A = A / sum(abs(A))*M**2

      S = matmul(B, matmul(S, transpose(B)))  ! note
      A = matmul(B, matmul(A, transpose(B)))  ! note

!      write(6,*)' A '
!      write(6,88)A

!      write(6,*)' S '
!      write(6,88)S

!      write(6,*)' A '
!      write(6,88)A

cc    call r_diag_normlz(A, Vc, w, M)
      call r_diag_generalized(A, Vc, S, w, M, L) 
      write(6,*)' after '

!      write(6,*)' A '
!      write(6,88) A
!      write(6,*)' S '
!      write(6,88) S
      write(6,*)' w '
      write(6,88) w
!      write(6,*)' Vc '
!      do i=1, M
!      write(6,88) Vc(i,:) 
!      enddo

 88   format(1x,4f15.6)
      
      end
c-----
c ! solves F aC = S_i aC E_diag. Wasteful on space.
c----

      subroutine r_diag_generalized(F, aC, S, eps, M, L)
      implicit none
      real*8, parameter :: tol_s = 1d-4
      integer M, i, iflg, L, ism
      integer j                         
      real*8 F(M, M), aC(M, M), S(M, M), eps(M), diff

      real*8,allocatable,dimension(:)   ::s_evl 
      real*8,allocatable,dimension(:,:)::X,F_rot,S_rot,E_diag,
     x                                          Y,aCsm
      allocate(X(M, M),s_evl(M), stat=j);if(j/=0) stop


      if(M<1.or.M>1000) stop  ! change

      call r_check_symm(F, M)
      call r_check_symm(S, M)
      call r_diag_normlz(S, X, s_evl, M)

      if(minval(s_evl) <0 . and. 
     x       abs(minval(s_evl))> tol_s*maxval(s_evl)) then
         write(6,*)' minval(s_evl), maxval(s_evl) '
         write(6,*)  minval(s_evl), maxval(s_evl)
         stop
      endif

      L = 0
      do i=1,M;if(s_evl(i)/maxval(s_evl)>tol_s)L=L+1;enddo
      if(L == 0)write(6,*) 'problem, L=0, tol_s_max ',maxval(s_evl)
      write(6,*)' reduced matrix size in can. diag.  ',L,' vs. ',M
      allocate (Y(M, L), stat=j); if(j /=0) stop
      
      ism = 0
      do i=1, M
         if(s_evl(i)/maxval(s_evl)>tol_s) then
            ism = ism + 1
            Y(:, ism) = X(:, i)/dsqrt(s_evl(i))
         endif
      enddo

      if(ism /= L) then
         write(6,*)' ism ',ism;stop
      endif

      deallocate(X)

      write(6,*)' s_evl ',s_evl
c----
c erase
c----
      allocate(S_rot(L, L),  stat=j);if(j/=0) stop
      S_rot = matmul(matmul(transpose(Y), S), Y)
      S_rot = (S_rot + transpose(S_rot))/2.d0

      do j=1,L; S_rot(j, j) = S_rot(j, j) - 1.d0; enddo
      diff = sum(abs(S_rot))

      if(diff>1.d-4) then  ! erase
         write(6,*)' diff in s_rot ', diff
         write(6,*) ' S_rot '
         write(6,88) S_rot
         write(6,*)' Y '
         write(6,88) Y
         stop
      endif

c---
c end of erase
c---
      deallocate(S_rot)
      allocate(F_rot(L, L), aCsm(L, L), stat=j);if(j/=0) stop
      F_rot = matmul(matmul(transpose(Y), F), Y)  
      F_rot = (F_rot + transpose(F_rot) )/2.d0

      eps = 0.d0
      call r_diag_normlz(F_rot, aCsm, eps, L)
      aC = 0.d0; do j=1,L;aC(:, j) = matmul(Y, aCsm(:, j));enddo
      deallocate(F_rot,  aCsm, Y)

c------
c now check - skip if wanted
c-----
      allocate(E_diag(M, M), stat=j);if(j/=0) stop
      
      iflg = 1
      if(iflg==1) then 

         E_diag = 0.d0 ; do i=1, L; E_diag(i, i)=eps(i); enddo !note: to L
         diff = sum(abs(matmul(F, aC) - matmul(matmul(S,aC),E_diag)))

         diff = diff / sum(abs(aC))/sum(abs(F))*M**2
         write(6,*)' scaled diff in gene_diag ',diff

ccc         call rcheck_small(diff, ' r_diag_gen ')

      endif
      deallocate(E_diag, s_evl)
     

 88   format(1x,4f15.6)

      end


