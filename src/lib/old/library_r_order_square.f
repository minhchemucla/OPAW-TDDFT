      subroutine order_square_driver  ! a dummy subroutine to drive order .

      real*8 w(8), w_dum(8)
      integer i_old_ofnew(8), M

      data  w / 0.8d0, 0.8d0, 0.3d0, 0.5d0, 0.8d0, 0.3d0, 0.9d0,0.9d0/
      save w
      M = 8

      call order_square_r_index(w,i_old_ofnew,8)

      do i=1, 8
         write(6,*)' w_old ',w(i)
      enddo

      w_dum = w; do i=1,M ; w(i) = w_dum(i_old_ofnew(i)) ;enddo

      do j=1,8
         write(6,*)' i_old_ofnew,w_new ',i_old_ofnew(j),
     x                                 w(j)
      enddo

      end

      subroutine order(w,m)
      implicit none
      integer m,st,i
      integer, allocatable :: i_old_ofnew(:)
      real*8 w(m)
      real*8, allocatable :: w_dum(:)
      
      allocate(i_old_ofnew(m),stat=st); call check0(st,'i_oldnew ')
      allocate(       w_dum(m),stat=st); call check0(st,'w_dum    ')
      call order_square_r_index(w,i_old_ofnew,m)
      w_dum=w
      do i=1,m
         w(i)=w_dum(i_old_ofnew(i))
      enddo
      deallocate(i_old_ofnew)
      deallocate(w_dum)
      end
      
      subroutine order_square_r_index(w,i_old_ofnew,M)
c----------------------
c A general ordering subroutine. Takes a vector w(M), and produces
c an index vector i_old_ofnew(M). Note: w is left intact. if you want
c to order it, define w_dum(m), and say: 
c   w_dum = w; do i=1,M ; w(i) = w_dum(i_old_ofnew(i)) ;enddo.
c---------------------
      implicit none
      integer   M
      integer  i_old_ofnew(M)
      real*8   w(M)

      integer ia, ib, i, j,st 
      real*8 wmin
      real*8, save :: big = 1e30

      integer, dimension(:), allocatable :: i_occ

      if(M>1e6.or.M<1) then
         write(6,*)' M ',M
         stop
      endif

      if(maxval(abs(w))>big) then
         write(6,*)' maxval(abs(w)) ',maxval(abs(w))
         stop
      endif

      allocate ( i_occ(M),stat=st ); call check(st,0,' i_occ ')
      
      i_occ = 0
      do i=1, M
         !if(mod(i,100)==1) write(6,*)' debug, i= ',i; 
         !if(mod(i,100)==1) call flush(6)

         wmin = big   ! note
         i_old_ofnew(i) = 0
         do j=1, M
            if(w(j).lt.wmin. and. i_occ(j).eq.0) then
               wmin = w(j)
               i_old_ofnew(i) = j
            endif
         enddo
         
         if(i_old_ofnew(i).eq.0) then
            write(6,*)' problem, ordering i ',i
            stop
         endif

         i_occ(i_old_ofnew(i)) = 1
      enddo
           
      if( maxval(i_occ).ne.1. or. minval(i_occ).ne.1) then
         write(6,*)'  i_occ ',i_occ
         stop
      endif

      if( maxval(i_old_ofnew ) > M. or. minval (i_old_ofnew) <1)then
         write(6,*)' problem. i_old_ofnew ',i_old_ofnew
         stop
      endif

      do i=1, M-1
         if(w(i_old_ofnew(i+1)).lt.w(i_old_ofnew(i))) then
            write(6,*)' problem . i, i+1 ',i, i+1
            write(6,*)' i_of_ofnew ',i_old_ofnew(i:i+1)
            write(6,*)' w_of_ofnew(i, i+1 ',w(i_old_ofnew(i)),
     x                                      w(i_old_ofnew(i+1))
            stop
         endif
      enddo

      do ia=1, M-1
         do ib=1, M-1
            if(ia.ne.ib.and.i_old_ofnew(ia).eq.i_old_ofnew(ib)) then
               write(6,*)' ia, ib, i_old_ofnew ',ia, ib, i_old_ofnew 
               stop
            endif
         enddo
      enddo


      deallocate( i_occ) 

      end


