module diis1_modu
  use main_mod, only : mix_di => mix_diis1
  implicit none
  save
  integer :: stage_diis1=0
  integer  kb1,k1
  integer, parameter :: ndi1_start=0 
  real*8, parameter   :: tiny1 = 1e-9  ! not used
  real*8, allocatable :: bank1(:,:), err1(:,:)
end module diis1_modu

subroutine diis1_general(vec, nn)
  use param_mod, only : L=>ndi1
  use diis1_modu,    only : mix_di,  bank1, err1,k1, kb1, stage_diis1
  use main_mod,         only : dv,iscf
  use mat_module,   only : mat_diag
  use mpi_lib_ours, only : rank
  implicit none
  integer st, i,j,m,nn
  real*8 vec(nn)
  real*8, allocatable :: B(:,:), Bv(:,:), Bw(:), C(:)

  call alloc_bank1_err(nn)

  if(stage_diis1==0) then
     stage_diis1=1
     !if(rddns) then;  
     !   if(rank==0) then
     !      write(6,*)' not ready for read_diis -- ignoring it '
     !   endif
     !   !call read_diis_all ! will do it only if not read earlier.
     !endif
     
     K1=0
     Kb1=1
     write(6,*)' size(bank1(:,kb1) ',size(bank1(:,kb1)) 
     write(6,*)' size(vec)       ',size(vec)
     bank1(:,Kb1)=vec
     return
  end if

!  if(rank==0) then
!     write(6,*)' L, mix_di, Kb1, K, iscf ', L, mix_di,Kb1,K,iscf
!  endif

  call time_print(' post L      ')

  call stackerr;                              call time_print(' poststack   ')
  call minimize;                              call time_print(' postminim   ')
  call predict_output;                        call time_print(' postpred    ')
  call bankstack;                             call time_print(' postbank    ')
  !call print_diis;                            call time_print(' postprin    ')
  
  call time_print(' post diis ')
  deallocate(B,Bv,Bw,C)
contains
  
  subroutine stackerr
    implicit none
    write(6,*)' K1, L,KB1 ',K1,L,KB1;    call flush(6)
    write(6,*)' allocated(err1,bank1) ',allocated(err1),allocated(bank1); call flush(6)
    write(6,*)' size(vec) ;',size(vec); call flush(6)
    write(6,*)' size(err1) ;',size(err1,1),size(err1,2); call flush(6)
    write(6,*)' size(bnk) ;',size(bank1,1),size(bank1,2); call flush(6)
    write(6,*)' |vec| ',sum(abs(vec)); call flush(6)
    write(6,*)' |bank1_1:Kb1| ',sum(abs(bank1(:,1:KB1))); call flush(6)

    if(K1==L) err1(:,1:K1-1)=err1(:,2:K1)
    K1=min(K1+1,L)
    err1(:,K1)= vec-bank1(:,Kb1)
  end subroutine stackerr

  subroutine minimize 
    use mat_module, only : mat_inv
    implicit none
    !
    ! Make diis matrix 
    !
    M=K1+1
    allocate (B(M,M), Bv(M,M), Bw(M), C(M), stat=st)
    call check0(st,'  BBw ')
    
    do i=1,K1
       do j=1,K1
          B(i,j)=sum(err1(:,i)*err1(:,j))*dv
       enddo
    enddo
    
    B(1:K1,M)= 1d0
    B(M,1:K1)= 1d0
    B(M, M) = 0d0
 
    C       = 0d0
    C(M)    = 1d0
    
    write(6,*)' M= ',size(B,1),M,' B= '
    do i=1,size(B,1)
       write(6,*)B(i,:)
    enddo
    write(6,*); call flush(6)
    call mat_inv(B,Bv)
    C  = matmul( Bv,C)

    if(rank==0) then
       write(6,*) ' sum(C(1:K1) ',sum(C(1:K1))
       write(6,*) '     C(1:K1) ',C(1:K1)
       write(6,*) ' iscf,norm(error1) ',iscf,B(K1,K1); write(6,*) ; call flush(6)
    end if

    C(1:K1)= C(1:K1)/sum(C(1:K1))
  end subroutine minimize

  subroutine predict_output
    implicit none
    vec = matmul(bank1(:,1:K1)+mix_di*err1(:,1:K1),C(1:K1))
  end subroutine predict_output

  subroutine bankstack
    implicit none
    if(Kb1==L) bank1(:,1:L-1)=bank1(:,2:L)
    call check(Kb1,K1,' Kb1 K1 ')
    Kb1=min(Kb1+1,L)
    bank1(:,Kb1) = vec
  end subroutine bankstack

  subroutine print_diis
    use main_mod, only : dens
    implicit none
    integer, parameter :: print_d = 1
    integer :: ik, ikb1

    if(print_d/=1) return
    if(rank/=0) return

    open(32,file='diis.txt',status='replace')
    rewind(32)
    
    write(32,*)Kb1,K1, ' Kb1, K1 '; call flush(32)
    write(32,*)nn,' nn '; call flush(32)
    write(32,*)' dens '; call flush(32)
    write(32,*)  dens; call flush(32)
    write(32,*)' bank1(:, 1:kb1) of v '; call flush(32)
    do ikb1=1,kb1
       write(32,*)  bank1(:,ikb1); call flush(32)
    enddo
    write(32,*)' err1(:,1:kb1)  of v'
    do ik=1,k1
       write(32,*)  err1(:,ik); call flush(32)
    enddo
    call flush(32)
    close(32)
  end subroutine print_diis
end subroutine diis1_general

subroutine alloc_bank1_err(nn)
  use param_mod, only : L=>ndi1
  use diis1_modu,    only : bank1, err1
  use mpi_lib_ours, only : rank
  implicit none
  integer st,nn
  
  if(.not.allocated(bank1)) then
     allocate(bank1(nn, L),stat=st)
     call check0(st,' bank1  ')
  end if
  
  if(.not.allocated(err1)) then
     allocate(err1(nn, L),stat=st)
     call check0(st,' err1   ')
  end if
end subroutine alloc_bank1_err


