module diis_modu
  use main_mod, only : mix_di => mix_diis
  implicit none
  save
  integer :: stage_diis=0
  integer  kb,k
  integer, parameter :: ndi_start=0 
  real*8, parameter   :: tiny = 1e-9  ! not used
  real*8, allocatable :: bank(:,:), err(:,:)
end module diis_modu

subroutine diis_general(vec, nn)
  use param_mod, only : L=>ndi
  use diis_modu,    only : mix_di,  bank, err,k, kb, stage_diis
  use main_mod,         only : dv,iscf
  use mat_module,   only : mat_diag
  use mpi_lib_ours, only : rank
  implicit none
  integer st, i,j,m,nn
  real*8 vec(nn)
  real*8, allocatable :: B(:,:), Bv(:,:), Bw(:), C(:)

  call alloc_bank_err(nn)

  if(stage_diis==0) then
     stage_diis=1
     !if(rddns) then;  
     !   if(rank==0) then
     !      write(6,*)' not ready for read_diis -- ignoring it '
     !   endif
     !   !call read_diis_all ! will do it only if not read earlier.
     !endif
     
     K=0
     Kb=1
     write(6,*)' size(bank(:,kb) ',size(bank(:,kb)) 
     write(6,*)' size(vec)       ',size(vec)
     bank(:,Kb)=vec
     return
  end if

!  if(rank==0) then
!     write(6,*)' L, mix_di, Kb, K, iscf ', L, mix_di,Kb,K,iscf
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
    write(6,*)' K, L,KB ',K,L,KB;    call flush(6)
    write(6,*)' allocated(err,bank) ',allocated(err),allocated(bank); call flush(6)
    write(6,*)' size(vec) ;',size(vec); call flush(6)
    write(6,*)' size(err) ;',size(err,1),size(err,2); call flush(6)
    write(6,*)' size(bnk) ;',size(bank,1),size(bank,2); call flush(6)
    write(6,*)' |vec| ',sum(abs(vec)); call flush(6)
    write(6,*)' |bank_1:Kb| ',sum(abs(bank(:,1:KB))); call flush(6)

    if(K==L) err(:,1:K-1)=err(:,2:K)
    K=min(K+1,L)
    err(:,K)= vec-bank(:,Kb)
  end subroutine stackerr

  subroutine minimize 
    use mat_module, only : mat_inv
    implicit none
    !
    ! Make diis matrix 
    !
    M=K+1
    allocate (B(M,M), Bv(M,M), Bw(M), C(M), stat=st)
    call check0(st,'  BBw ')
    
    do i=1,K
       do j=1,K
          B(i,j)=sum(err(:,i)*err(:,j))*dv
       enddo
    enddo
    
    B(1:K,M)= 1d0
    B(M,1:K)= 1d0
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
       write(6,*) ' sum(C(1:K) ',sum(C(1:K))
       write(6,*) '     C(1:K) ',C(1:K)
       write(6,*) ' iscf,norm(error0) ',iscf,B(K,K); write(6,*) ; call flush(6)
    end if

    C(1:K)= C(1:K)/sum(C(1:K))
  end subroutine minimize

  subroutine predict_output
    implicit none
    vec = matmul(bank(:,1:K)+mix_di*err(:,1:K),C(1:K))
  end subroutine predict_output

  subroutine bankstack
    implicit none
    if(Kb==L) bank(:,1:L-1)=bank(:,2:L)
    call check(Kb,K,' Kb K ')
    Kb=min(Kb+1,L)
    bank(:,Kb) = vec
  end subroutine bankstack

  subroutine print_diis
    use main_mod, only : dens
    implicit none
    integer, parameter :: print_d = 1
    integer :: ik, ikb

    if(print_d/=1) return
    if(rank/=0) return

    open(32,file='diis.txt',status='replace')
    rewind(32)
    
    write(32,*)Kb,K, ' Kb, K '; call flush(32)
    write(32,*)nn,' nn '; call flush(32)
    write(32,*)' dens '; call flush(32)
    write(32,*)  dens; call flush(32)
    write(32,*)' bank(:, 1:kb) of v '; call flush(32)
    do ikb=1,kb
       write(32,*)  bank(:,ikb); call flush(32)
    enddo
    write(32,*)' err(:,1:kb)  of v'
    do ik=1,k
       write(32,*)  err(:,ik); call flush(32)
    enddo
    call flush(32)
    close(32)
  end subroutine print_diis
end subroutine diis_general

subroutine alloc_bank_err(nn)
  use param_mod, only : L=>ndi
  use diis_modu,    only : bank, err
  use mpi_lib_ours, only : rank
  implicit none
  integer st,nn
  
  if(.not.allocated(bank)) then
     allocate(bank(nn, L),stat=st)
     call check0(st,' bank  ')
  end if
  
  if(.not.allocated(err)) then
     allocate(err(nn, L),stat=st)
     call check0(st,' err   ')
  end if
end subroutine alloc_bank_err

