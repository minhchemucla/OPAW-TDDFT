subroutine r_diag_subparts(A, Vc, eps, n)
!
! this subouritne checks for isolated subparts
!
  implicit none
  integer n

  real*8 A(n, n), Vc(n, n), eps(n), delta
  real*8, dimension(:, :), allocatable :: Apart, Vpart
  real*8, parameter :: tol = 1.d-10

  integer i, j, l, block_present, block_include, block_remove
  integer mn, mx, iblock, max_block, isymm, ipntr, nb, flag_check

  integer, allocatable, dimension(:) ::  block, map_shuffle, num_in_block
  integer, allocatable, dimension(:) ::  min_of_block, max_of_block, &
                                         inv_map_shuffle, i_old_ofnew
  allocate(          block(n), stat=i); if(i/=0) stop
  allocate(    map_shuffle(n), stat=i); if(i/=0) stop
  allocate(num_in_block(n),    stat=i); if(i/=0) stop
  allocate(min_of_block(n),    stat=i); if(i/=0) stop
  allocate(max_of_block(n),    stat=i); if(i/=0) stop
  allocate(inv_map_shuffle(n), stat=i); if(i/=0) stop
  allocate(i_old_ofnew(n), stat=i); if(i/=0) stop

  block_present = 0;   block(:) = 0

  isymm=0;  if(maxval(abs(A-transpose(A))) < tol) isymm = 1

  do i=1, n

     if(block(i) == 0) then
        block_present = block_present + 1
        block(i)      = block_present
     endif
        
     do j=1, n
        if(j/=i) then
           if(abs(A(i, j)) > tol .or. abs(A(j, i)) > tol ) then 
              if(block(j) == 0) then
                 block(j) = block(i)
              else if(block(j) /= block(i)) then
                 block_remove  = max(block(i), block(j))
                 block_include = min(block(i), block(j))
                 block_present =  block_present - 1
                 do l=1, n
                    if(block(l) == block_remove) block(l) = block_include
                 enddo
              endif
           endif
        endif
     enddo
  enddo

  max_block = block_present

!
! now that you divided into blocks, divide to several parts. 
!  

  map_shuffle(:) = 0
  inv_map_shuffle(:) = 0

  ipntr = 0
  do iblock = 1, max_block

     num_in_block(iblock) = 0
     min_of_block(iblock) = 100000
     max_of_block(iblock) = 0

     do i=1, n
        if(block(i) == iblock) then
           ipntr = ipntr + 1
           map_shuffle(i) = ipntr
           inv_map_shuffle(ipntr) = i

           num_in_block(iblock) = num_in_block(iblock) + 1
           min_of_block(iblock) = min(ipntr, min_of_block(iblock) )
           max_of_block(iblock) = max(ipntr, max_of_block(iblock) )
        endif
     enddo
  enddo

  do iblock=1, max_block
     if(max_of_block(iblock)-min_of_block(iblock) /= num_in_block(iblock)-1) &
          stop
  enddo
  do iblock=2, max_block-1
     if(max_of_block(iblock)   /=min_of_block(iblock+1)-1 .or. &
        max_of_block(iblock-1) /=min_of_block(iblock)-1 ) stop
  end do
  if(min_of_block(1) /=1) stop
  if(max_of_block(max_block) /= n) stop

  if (maxval(map_shuffle) /=n .or.minval(map_shuffle) /= 1) stop
  if (maxval(inv_map_shuffle) /= n .or. minval(inv_map_shuffle) /= 1) stop
  if( sum(inv_map_shuffle) /= ((n+1)*n)/2 ) stop
  if( sum(    map_shuffle) /= ((n+1)*n)/2 ) stop

  !
  ! check: can erase
  ! 

  do i=1, n
     do j=1, n
        if(block(i) /= block(j) .and. (abs(A(i,j))>tol.or.abs(A(j,i))>tol))then
           write(6,*)' i, block(i), j, block(j) ',i,block(i),j, block(j)
           write(6,*)' A ',A(i, j), A(j, i)
           stop
        endif
     enddo
  enddo

  A(:, :) = A(inv_map_shuffle(:), inv_map_shuffle(:) )  ! check !

  !
  ! end of check
  !
  Vc = 0.d0
  ipntr = 0
  do iblock = 1, max_block
     
     nb = num_in_block(iblock)
     mn = min_of_block(iblock)
     mx = max_of_block(iblock)
     ipntr = ipntr + nb

     if(mx-mn+1 /= nb) stop
     if(ipntr   /= mx) stop

     allocate(Apart(nb, nb), stat=i); if(i/=0) stop
     allocate(Vpart(nb, nb), stat=i); if(i/=0) stop

     Apart = A(mn:mx, mn:mx)
     if(isymm == 1 .and. maxval(abs(transpose(Apart)-Apart))>tol) stop
     
     select case (isymm)
        case(1) 
        call r_diag_normlz(Apart, Vpart, eps(mn) , nb)

        case(0)
        write(6,*)' put here a matrix diagonalization for non-symm '
        stop
        
     end select

     Vc(mn:mx, mn:mx) = Vpart

     deallocate(Apart)
     deallocate(Vpart)

  enddo

  A (:, :) =  A(map_shuffle(:), map_shuffle(:) )  ! check 
  Vc(:, :) = Vc(map_shuffle(:), :)

  call order_r_index(eps,i_old_ofnew,n)
  eps =       eps( i_old_ofnew(:))
  Vc(:, :) = Vc(:, i_old_ofnew(:))

  !
  ! check diagonalization - can be erased
  !
  
  flag_check = 1
  if(flag_check == 1) then
     do j=1, n
        if(maxval(abs(matmul(A, Vc(:, j))-eps(j)*Vc(:, j)))>tol) then
           write(6,*)' problem in diag ',j
           stop
        endif
     enddo
  endif

  if(flag_check == 1) then
     do j=1, n
        do l=1, n
           delta = 0.d0
           if(j==l) delta = 1.d0

           if(abs(sum(Vc(:,j)*Vc(:,l))-delta)>1.d-6) then
              write(6,*)' problem in orthogonality ',j,l, &
                          sum(Vc(:,j)*Vc(:,l))
              stop
           endif
        enddo
     enddo
  end if
              
           

  deallocate( block, map_shuffle, num_in_block, min_of_block, max_of_block, &
              inv_map_shuffle, i_old_ofnew)
end subroutine r_diag_subparts

  subroutine r_diag_subparts_driver()
!!program    r_diag_subparts_driver

  implicit none
  integer, parameter :: n=7
  real*8   Vc(n, n), A(n, n), eps(n)
  
  integer i, j

  A=0.d0
  do j=1, 7, 3
     do i=1, 7, 3
        A(j, i) = 2
        if(j/=i) A(j, i) = 1
     enddo
  enddo

  do j=2, 5, 3
     do i=2, 5, 3
        A(j, i) = 2
        if(j/=i) A(j, i) = 1
     enddo
  enddo

  do j=3, 6, 3
     do i=3, 6, 3
        A(j, i) = 2
        if(j/=i) A(j, i) = 1
     enddo
  enddo

  write(6,*)' A '
  write(6,888)transpose(A)
  write(6,*)' '

  call r_diag_subparts(A, Vc, eps, n) 
 
  write(6,*)' eps '
  write(6,888)eps
  write(6,*)' '

  write(6,*)' Vc '
  write(6,888)transpose(Vc)
  write(6,*)' '

888 format(' ',7f11.2)

end
