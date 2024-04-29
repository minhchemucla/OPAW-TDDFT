      subroutine equate_cr(cca,    rra,       nn)
      implicit none
      integer nn, i
      real*8     rra(2,nn)
      complex*16 cca(nn)


      do  i=1, nn
         rra(1, i) = cca(i)
         rra(2, i) = aimag(cca(i))
      enddo

      write(6,*)' stage 33334 '

      end subroutine

c real*8 matrix plot

      subroutine r_plot_mt(x_mt, a, n1, n2, iout)
      implicit none
      integer n1, n2, i, j, iout, m1, m2
      character(*) a
      real*8 x_mt(n1, n2)
      write(6,*)' plottig matrix ',a, ' n1*n2 = ',n1,'*',n2
      
      if(n1*n2>14*14) then
         write(6,*)' plotting only up to min(n1,14), min(n2,14) '
         m1 = min(n1, 14)
         m2 = min(n2, 14)
      else
         m1 = n1
         m2 = n2
      endif
         
      if(n1<1.or.n2<1) then
         write(6,*)' problem in x_mt with ',a,' n1,n2=',n1,n2
         stop
      endif

      write(iout,*)

      do i=1, m1
         do j=1, m2
            write(iout,88)i,j,x_mt(i, j)
 88         format(' ',2i8,4x,f22.9)

         enddo
         write(iout,*)
      enddo
      end


c

c complex*16 matrix plot

      subroutine c_plot_mt(c_mt, a, n1, n2, iout)
      implicit none
      integer n1, n2, i, j, iout, m1, m2
      character(*) a
      complex*16 c_mt(n1, n2)
      write(6,*)' plottig matrix ',a, ' n1*n2 = ',n1,'*',n2
      
      if(n1*n2>14*14) then
         write(6,*)' plotting only up to min(n1,14), min(n2,14) '
         m1 = min(n1, 14)
         m2 = min(n2, 14)
      else
         m1 = n1
         m2 = n2
      endif
         
      if(n1<1.or.n2<1) then
         write(6,*)' problem in c_mt with ',a,' n1,n2=',n1,n2
         stop
      endif

      write(iout,*)

      do i=1, m1
         do j=1, m2
            write(iout,88)i,j,dble(c_mt(i,j)),aimag(c_mt(i,j))
 88         format(' ',2i8,4x,f22.9,' +ci*',f22.9)

         enddo
         write(iout,*)
      enddo
      end


c complex grahm schmidt. No weight. NO CONJUGATION.
      subroutine c_Grahm_Schmidt_wght_1( Cv, N, Nvec)   ! no conjguation
      implicit none                                     ! thus no dot.prdct
      integer N, Nvec, i, j
      complex*16 Cv(N, Nvec)
      
      do i=1, Nvec
         do j=1, i-1
            Cv(:, i) = Cv(:, i)- sum(Cv(:,i)*Cv(:, j)) * Cv(:,j)
            call  rcheck_small(abs( sum(cv(:, i)*cv(:, j))),'cGs_ij ')
         enddo
         Cv(:, i) = Cv(:, i) /   sqrt(sum(Cv(:, i)*Cv(:, i)))
         call  rcheck_small(abs( sum(cv(:, i)*cv(:, i)) - 1.d0),
     x                                           'cGs_norm   ')
      enddo
      end

 
      subroutine r_Grahm_Schmidt_wght( pr, wg, N, Nvec) 
      implicit none                                     ! thus no dot.prdct
      integer N, Nvec, i, j
      real*8 pr(n, nvec), wg(n)
      
      do i=1, Nvec
         do j=1, i-1
            Pr(:, i) = Pr(:, i)- sum((Pr(:,j))*Pr(:, i)*wg) * Pr(:,j)
             call rcheck_small(
     x          abs( sum((pr(:, i))*pr(:, j)*wg)), ' cGs_ij  ')
         enddo
         Pr(:, i) = Pr(:, i) /   sqrt(sum((Pr(:, i))*Pr(:, i)*wg))
         call 
     x   rcheck_small(abs( sum((pr(:, i))*pr(:, i)*wg) - 1.d0),
     x                                                   ' cGs_norm')
      enddo
      end

c complex grahm schmidt. No weight. 
      subroutine c_Grahm_Schmidt_wght_conj( Cv, N, Nvec) 
      implicit none                                     ! thus no dot.prdct
      integer N, Nvec, i, j
      complex*16 Cv(N, Nvec)
      
      do i=1, Nvec
         do j=1, i-1
            Cv(:, i) = Cv(:, i)- sum(conjg(Cv(:,j))*Cv(:, i)) * Cv(:,j)
             call rcheck_small(
     x          abs( sum(conjg(cv(:, i))*cv(:, j))), ' cGs_ij  ')
         enddo
         Cv(:, i) = Cv(:, i) /   sqrt(sum(conjg(Cv(:, i))*Cv(:, i)))
         call 
     x   rcheck_small(abs( sum(conjg(cv(:, i))*cv(:, i)) - 1.d0),
     x                                                   ' cGs_norm')
      enddo
      end


      subroutine r_Grahm_Schmidt( Rv, N, Nvec, dv) 
      implicit none             
      integer N, Nvec, i, j
      real*8  Rv(N, Nvec)
      real*8     dv
      
      do i=1, Nvec
         do j=1, i-1
            Rv(:, i) = 
     &       Rv(:, i) -sum(Rv(:,j)*Rv(:, i))*dv * Rv(:,j)
         enddo
         Rv(:, i) = Rv(:, i) /   sqrt(sum(Rv(:, i)*Rv(:, i))*dv)
      enddo
      end

      subroutine c_Grahm_Schmidt( Cv, N, Nvec, dv) 
      implicit none             
      integer N, Nvec, i, j
      complex*16 Cv(N, Nvec)
      real*8     dv
      
      do i=1, Nvec
         do j=1, i-1
            Cv(:, i) = 
     &       Cv(:, i) -sum(conjg(Cv(:,j))*Cv(:, i))*dv * Cv(:,j)
         enddo
         Cv(:, i) = Cv(:, i) /   sqrt(sum(conjg(Cv(:, i))*Cv(:, i))*dv)
      enddo
      end

c checks if a number is small; if yes, print A_char and stops

      subroutine rcheck_small(r, A_char)
      implicit none
      real*8 r
      character(*) A_char
      if(abs(r)>1.d-7) then
         write(6,*)' in r_check_small_problem. r= ',r
         write(6,*)' it is ', A_char
         stop
      endif
      end
     

c like rcheck, but for complex.

      subroutine c_check_symm(Ca, N)
      implicit none
      integer N, i, j
      complex*16 Ca(N,N)
      real*8 maxc
      maxc = maxval(abs(ca))
      do i=1, N
         do j=1, N
            if( abs(Ca(i, j)-Ca(j, i))>1.d-8*maxC) then
               write(6,*)' in check_symm, N = ',N
               write(6,*)' i, j, Ca(i,j ), Ca(j, i)',i, j
               write(6,*)        Ca(i, j), Ca(j, i)
               stop
            endif
         enddo
      enddo
      end

c checks if A is symmetrc.
      subroutine r_check_symm(A, N)
      implicit none
      integer N, i, j
      real*8 A(N,N),maxa
      maxa = maxval(abs(a))
      do i=1, N
         do j=1, N
            if( abs(A(i, j)-A(j, i))>1.d-8) then
               write(6,*)' in check_symm, N = ',N
               write(6,*)' i, j, A(i,j ), A(j, i)', i, j
               write(6,*)        A(i, j), A(j, i)
               stop
            endif
         enddo
      enddo
      end

c trace
      function r_trace(rm,N)
      implicit none
      integer N, i
      real*8 rm(N,N),r_trace

      r_trace = 0.d0
      do i=1, N
         r_trace = r_trace + rm(i, i)
      enddo
      end
c-----
c     
c----- 

c      subroutine cinv_vec_hmE(cA, c_src, chi, E, Md)!  chi=1/(ca-E)*c_src
c      implicit none
c      integer Md, i
c      real*8 E
c      complex*16 cA(Md, Md), c_src(Md), chi(Md)
c
c      do i=1,Md; cA(i,i) = cA(i,i)-E ; enddo
c      call cinv_vec(cA, c_src, chi, Md)  ! calcs. chi = ca_inv*c_src
c      do i=1,Md; cA(i,i) = cA(i,i)+E ; enddo
c      end
c
c      subroutine cinv_vec(cA, c_src, chi, Md)  ! calcs. chi = ca_inv*c_src
c      implicit none
c      integer Md, st, i, j
c      complex*16 cA(Md, Md), c_src(Md), chi(Md)
c      complex*16, dimension(:,:), allocatable :: cA_inv
c      allocate (cA_inv(Md,Md), stat=st) ; if(st.ne.0) stop
c      cA_inv= cA
c
c      call c_check_symm(cA, Md)
c
c      call zinvmtrx(cA_inv, Md, Md)
c      
c      chi = matmul(cA_inv,c_src)
c
c      deallocate(cA_inv)
c      end

c returns 1 if i, j, k, l, are different

      function i_4alldiff( i, j, k, l)
      implicit none
      integer i, j, k, l, i_4alldiff
      
      if(i.eq.j.or.i.eq.k.or.i.eq.l.or.j.eq.k.or.j.eq.l.or.k.eq.l)then
         i_4alldiff=0
      else
         i_4alldiff=1
      endif
      end

c retruns 1 if i and j and k are all different

      function i_3alldiff( i, j, k)
      implicit none
      integer i, j, k,  i_3alldiff
      
      if(i.eq.j.or.i.eq.k.or.j.eq.k) then
         i_3alldiff=0
      else
         i_3alldiff=1
      endif
      end

c  returns 1 i or j or k or l are full.

      function i_4anyfull( i, j, k, l, iO)
      implicit none
      integer i, j, k, l, iO(*), i_4anyfull
      
      if(iO(i).ne.0.or.iO(j).ne.0.or.iO(k).ne.0.or.iO(l).ne.0)then
         i_4anyfull = 1
      else
         i_4anyfull = 0  ! all empty
      endif
      end

c similar to previous one

      function i_3anyfull( i, j, k, iO)
      implicit none
      integer i, j, k,  iO(*), i_3anyfull
      
      if(iO(i).ne.0.or.iO(j).ne.0.or.iO(k).ne.0)then
         i_3anyfull = 1
      else
         i_3anyfull = 0  ! all empty
      endif
      end

c kronecker delta; reutrns REAL*8 results.
      function delta(i,j)
      implicit none
      integer i,j
      real*8 delta
      if(i.eq.j) then
         delta = 1.d0
      else
         delta = 0.d0
      endif
      end

c if i not equal to j, stop and print achar.

      subroutine check(i,j,achar)
      implicit none
      integer i,j
      character(*) achar
      if(i.ne.j) then 
         write(6,*)' stopping since ',i,' n.e. ', j,' in ', achar 
         stop
      endif
      end

      subroutine check0(i,achar)
      implicit none
      integer i
      character(*) achar
      if(i.ne.0) then 
         write(6,*)' stopping since,',i,' n.e.0 in ', achar
         stop
      endif
      end


c if x not close to y, stop and print achar.  Depends on tolerance.
      subroutine check_r(x,y,achar)
      implicit none
      real*8 x, y
      character(*) achar
      if(abs(x-y)/max(abs(x),abs(y),1.d-20)>1.d-12) then 
         write(6,*)x, y, achar 
         stop
      endif
      end


c if i not less or eqaul to j, stop and print achar.



      subroutine check_LE(i,j,achar)
      implicit none
      integer i,j
      character(*) achar
      if(i > j) then 
         write(6,*)' problem in check_le ',i, j, achar 
         stop
      endif
      end


      subroutine check_lele(i,j,k,achar)
      implicit none
      integer i,j,k
      character(*) achar
      if(i > j .or. j>k) then 
         write(6,*)' problem in check_lele ',i, j, k,achar 
         stop
      endif
      end

      subroutine check_lelele(i,j,k,l,achar)
      implicit none
      integer i,j,k,l
      character(*) achar
      if(i > j .or. j>k .or. k>l) then 
         write(6,*)' problem in check_lelele ',i, j, k,l,achar 
         stop
      endif
      end

      subroutine check_lelelele(i,j,k,l,m,achar)
      implicit none
      integer i,j,k,l,m
      character(*) achar
      if(i > j .or. j>k .or. k>l .or. l>m) then 
         write(6,*)' problem in check_lelelele ',i, j, k,l,m,achar 
         stop
      endif
      end


      subroutine check_r_LE(x,y,achar)
      implicit none
      real*8 x,y
      character(*) achar
      if(x > y) then 
         write(6,*)x, y, achar 
         stop
      endif
      end

      subroutine check_r_lele(x,y,z,achar)
      implicit none
      real*8 x,y,z
      character(*) achar
      if(x > y.or.y>z) then 
         write(6,*)x, y, z, achar 
         stop
      endif
      end


      subroutine check_Lt(i,j,achar)
      implicit none
      integer i,j
      character(*) achar
      if(i .ge.  j) then 
         write(6,*)i, j, achar 
         stop
      endif
      end

c if i not less or eqaul to j, stop and print achar.

      subroutine check_real_LE(x,y,achar)
      implicit none
      real*8 x,y
      character(*) achar
      if(x > y) then 
         write(6,*)x, y, achar 
         stop
      endif
      end

      subroutine check_real_Lt(x,y,achar)
      implicit none
      real*8 x,y
      character(*) achar
      if(x .ge.  y) then 
         write(6,*)x, y, achar 
         stop
      endif
      end



     
c------
c cu = cu+E*Unity
c------
      subroutine cpls( cu_rd, E, Md)
      implicit none
      integer Md, i
      real*8 E
      complex*16 cu_rd(Md, Md)
      do i=1, Md
         cu_rd(i, i) = cu_rd(i, i)+ E
      enddo
      end

c------
c cu = cu-E*Unity
c------
      subroutine cmns( cu_rd, E, Md)
      implicit none
      integer Md, i
      real*8 E
      complex*16 cu_rd(Md, Md)
      do i=1, Md
         cu_rd(i, i) = cu_rd(i, i)- E
      enddo
      end
c-----------------------------------
c ! checks if idet=1,2,4,8,16,etc. (1 then; -1 otherwise)
c-----------------------------------
      function is_it_power2(i_in) 
      implicit none
      integer is_it_power2,i,i_in
      
      i = abs(i_in)
      if(i==0.or.i==1) then
         is_it_power2 = 1
         return
      endif

      is_it_power2= -1
      if(  2**nint((dlog(dble(i))/dlog(2.d0)))==i) 
     &                                     is_it_power2 = 1
      end function is_it_power2

      function is_i_near_power_x(ii,x)
      implicit none
      integer ii, is_i_near_power_x, j, ans
      real*8 x,y
      
      j = int(log(dble(ii))/log(x))
      y = x**dble(j)
      
      if(y.le.ii.and.y+1.gt.ii) then
         ans = 1
      else
         ans = 0
      end if
      
      is_i_near_power_x = ans
      end function is_i_near_power_x

      subroutine fetch_r(x,aa,iplt) ! routine to read real*8 variable
                                ! Also read character string and
                                ! checks that it equals to input string.
      use mpi_lib_ours, only : rank
      implicit none
      integer iplt
      real*8 x
      character (*) aa
      character*20  bb
      
      if(rank==0) then;
         write(6,*)' now reading ',aa,' from ',iplt; call flush(6)
      endif
      read(iplt, *) bb, x
      if(bb/=aa) then
         write(6,*)' problem; supposed read ',aa,' instead reads ',bb; 
         call flush(6)
         stop
      endif
      if(rank==0)write(6,*)' = ',x
      end subroutine fetch_r

      subroutine fetch_i(ii,aa,iplt) ! same for integer.
      use mpi_lib_ours, only : rank
      implicit none
      integer ii, iplt
      character (*) aa
      character*20 bb
      
      if(rank==0)then;
         write(6,*)' now reading ',aa,' from ',iplt; call flush(6)
      endif
      read(iplt, *) bb, ii
      if(bb/=aa) then
         write(6,*)' problem; supposed read ',aa,' instead reads ',bb
         call flush(6)
         stop
      endif
      if(rank==0)write(6,*)' = ',ii
      end subroutine fetch_i

      subroutine fetch_L(L,aa,iplt) 
      use mpi_lib_ours, only : rank
      implicit none
      logical L
      integer iplt
      character(*) aa
      character*20 bb
      
      if(rank==0)then
         write(6,*)' now reading ',aa,' from ',iplt
         call flush(6)
      endif
      read(iplt, *) bb, L
      if(bb/=aa) then
         write(6,*)' problem; supposed read ',aa,' instead reads ',bb
         call flush(6)
         stop
      endif
      if(rank==0)write(6,*)' = ',L
      end subroutine fetch_L

      subroutine get_I( N, Achar, ifile)
      use mpi_lib_ours, only : rank
      implicit none
      integer N, ifile
      character(*) Achar
      if(rank==0) write(6,*   )' give me ',Achar
      read(ifile,*) N
      if(rank==0)write(6,*)' = ', N
      end

      subroutine get_r_mtrx( xmat,     n1, n2, achar ,ifile)
      implicit none
      integer  n1, n2, ifile
      real*8   xmat(n1, n2), term
      character(*) Achar

      integer i1, i1dum
      integer i2, i2dum
      integer iread

      write(6,*)' give me ',Achar
      write(6,*)' which is a matrix of size',n1,' times ',n2
      write(6,*)' use the format i1, i2, xmat(i1,i2)  (i1 vfaster) '

      iread = 0
      do i2=1,n2
         do i1=1,n1
            read(ifile, *,end=99)i1dum, i2dum, term
            call check(i1, i1dum, ' i1xmt ')
            call check(i2, i2dum, ' i2xmt ')
            xmat(i1,i2) = term
            iread = 1

         enddo
      enddo

      goto 100

 99   continue
      write(6,*)' if reached here there was problem '
      write(6,*)' i1, i2, i1dum, i2dum, term '
      write(6,*)  i1, i2, i1dum, i2dum, term
      stop
      
 100  continue
      
      if(iread==0) then
         write(6,*)' xmat not read '; stop
      endif

      write(6,*)' finished reading the matrix.  Printing some info: '
      
      write(6,*)' minval(xmat) ',minval(xmat)
      write(6,*)' maxval(xmat) ',maxval(xmat)
      write(6,*)' minabs(xmat) ',minval(abs(xmat))
      end

      subroutine get_r_vctr( xvec,     n1, achar ,ifile)
      implicit none
      integer  n1, ifile
      integer iread
      real*8   xvec(n1), term
      character(*) Achar

      integer i1, i1dum

     
      write(6,*)' give me ',Achar
      write(6,*)' which is a vector of size',n1


      iread = 0
      do i1=1,n1
         read(ifile, *,end=99)i1dum,  term
         call check(i1, i1dum, ' i1vmt ')
         xvec(i1) = term
         iread    = 1
      enddo
      goto 100

 99   continue
      write(6,*)' if reached here there was problem '
      write(6,*)' i1,  i1dum, term '
      write(6,*)  i1,  i1dum, term
      stop
      
 100  continue

      if(iread==0) then
         write(6,*)' xvec not read '; stop
      endif

      write(6,*)' finished reading the vector.  Printing some info: '
      
      write(6,*)' minval(xvec) ',minval(xvec)
      write(6,*)' maxval(xvec) ',maxval(xvec)
      write(6,*)' minabs(xvec) ',minval(abs(xvec))
      end


      subroutine get_i_vctr( ivec,     n1, achar ,ifile)
      implicit none
      integer  n1, ifile
      integer iread
      integer   ivec(n1), term
      character(*) Achar

      integer i1, i1dum

     
      write(6,*)' give me ',Achar
      write(6,*)' which is a vector of size',n1


      iread = 0
      do i1=1,n1
         read(ifile, *,end=99)i1dum,  term
         call check(i1, i1dum, ' i1imt ')
         ivec(i1) = term
         iread    = 1
      enddo
      goto 100

 99   continue
      write(6,*)' if reached here there was problem '
      write(6,*)' i1,  i1dum, term '
      write(6,*)  i1,  i1dum, term
      stop
      
 100  continue

      if(iread==0) then
         write(6,*)' ivec not read '; stop
      endif

      write(6,*)' finished reading the vector.  Printing some info: '
      
      write(6,*)' minval(ivec) ',minval(ivec)
      write(6,*)' maxval(ivec) ',maxval(ivec)
      write(6,*)' minabs(ivec) ',minval(abs(ivec))
      end

      subroutine get_r( X, Achar, ifile)
      implicit none
      integer ifile
      real*8  X
      character(*) Achar
      write(6,*   )' give me ',Achar
      read(ifile,*) X
      write(6,*)' = ', X
      end
    
      subroutine get_A( A, len, Achar, ifile)
      implicit none
      integer ifile, len
      character   A(len)
      character(*) Achar
      write(6,*   )' give me ',Achar
      read(ifile,*) A
      write(6,*)' = ', A
      end

      subroutine timerd( A)
      implicit none
      character(*) A
      real*8 :: time = 0.d0, tnew
      integer :: ifirst = 1
      save time
      save ifirst

      if(ifirst==1) then
         ifirst=-1
         call cpu_time(time)
         write(6,*)' first time cal to timer, measured from here '
         return
      else
         call cpu_time(tnew)
         write(6,*)'time for',A,' = ',tnew-time
         time = tnew
      endif
      
      end subroutine

      subroutine check_real(cvec, n)
      implicit none
      integer n
      complex*16 cvec(n)
      real*8, parameter :: thresh_abs = 1.d-10
      real*8, parameter :: thresh_rel = 1.d-8
      real*8 normi, norma
      
      norma= sum(abs(cvec))
      normi= sum(abs(aimag(cvec)))
      
      if(normi>norma*thresh_rel.and.normi>thresh_abs) then
         write(6,*)' in check_real, normi = ',normi,' norma = ',norma
         stop
      endif
      end subroutine check_real
  
      function dble_interp(n, va,ra,r)   ! inear interpolation to give v(r), given two arrays ra, va, from va(ra).
      implicit none
      integer n, k,j,jm
      real*8 va(n), ra(n), r, dble_interp, dr,ddr
      if(r>ra(n)) then
         dble_interp = va(n) ; return
      else if (r<ra(1)) then
         dble_interp = va(1) ; return
      else 
         j=1
         k=n
         searchi : do
            if(k-j>1) then
               jm = (j+k)/2
               if(r<ra(jm)) then
                  k=jm
               else
                  j=jm
               end if
            else
               exit searchi
            end if
         end do searchi

! determined the point below r
         call check(k-j,1,'ikmj  ')
         
! now linear interpolation
         ddr = r      -ra(j)
         dr  = ra(k)  -ra(j)
         dble_interp      = va(j) * (dr-ddr)/dr + va(k) * ddr/dr
      endif
      end
      
      subroutine stopsub(a)
      implicit none
      character(*) a
      write(6,*)' problem, stopping :   ',a
      call flush(6)
      stop
      end 

      subroutine stamp(a)
      use mpi_lib_ours, only : rank
      implicit none
      character(*) a
      if(rank==0)then
         write(6,*)' stage-stamp : ',a
         call flush(6)
      end if
      end subroutine stamp

      real*8 function d2_3d(r,rp)
      implicit none
      real*8 r(3), rp(3)
      d2_3d = sum((r-rp)**2)
      end
      
