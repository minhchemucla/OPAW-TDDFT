subroutine set_hminmax
  use main_mod, only : hmin, hmax, hmin_0, hmax_0, read_hminmax
  use main_mod, only : iscf
  use param_mod, only : iscf_hminmax
  use mpi_lib_ours, only : rank
  implicit none
  
  if(iscf_hminmax<0) then
     hmin(:) = hmin_0
     hmax(:) = hmax_0
  else
     if (iscf<iscf_hminmax) call get_hminmax
  end if
  if(rank==0)write(6,*)' iscf: rank, hmin, hmax, iscf_hminmax ',iscf,rank, hmin,hmax, iscf_hminmax
  !write(6,*)' iscf: rank, hmin, hmax, iscf_hminmax ',iscf,rank, hmin,hmax, iscf_hminmax
end subroutine set_hminmax

subroutine get_hminmax
  use param_mod, only : niter,fctr,k_flg
  use main_mod,     only : n=>nn, vks, hmin, hmax, dv,nk_loc, h_type 
  use mpi_lib_ours
  implicit none
  integer i, st, ik
  real*8             :: e1(nk_loc), e2(nk_loc), havg,dh
  complex*16, allocatable :: z(:), hz(:)
  real*8, allocatable :: tmp(:)

  call debug_dhhavg

  allocate( z(n), stat=st)
  if(st/=0) stop 'z alloc problem'
  allocate(hz(n), stat=st)
  if(st/=0) stop 'hz alloc problem'
  allocate(tmp(n),stat=st)
  if(st/=0) stop 'tmp alloc in hminmax'
  
  do ik=1,nk_loc
      call rand_r(tmp,n,dv)
      if(.not.k_flg) call bcast_r8(tmp,n,0)
      z=tmp
    !  if(rank==0) then
         do i=1,niter
            z =z/sqrt(sum(abs(z)**2))
    !        call model_h(z,hz)
            if(h_type .eq. 1) then
              call shs(z,hz,ik)
            else
              call sh(z,hz,ik)
            endif
            !call hc(z, hz)
            !e1=sum(conjg(z)*hz)
            e1(ik)=sum(conjg(z)*hz)
            if(rank==0) then
                write(17,*) ' i, e1 ' ,i,e1(ik); call flush(17)
            endif
            z=hz
         enddo
    !  endif
    
      call rand_r(tmp,n,dv)
      if(.not.k_flg) call bcast_r8(tmp,n,0)
      z=tmp
    !  if(rank==0) then
         do i=1,niter
            z =z/sqrt(sum(abs(z)**2))
    !        call model_h(z,hz)
            if(h_type .eq. 1) then
              call shs(z,hz,ik)
            else
              call sh(z,hz,ik)
            endif
            !call shs(z,hz,ik)
    !        call hc(z, hz)
            hz = hz-e1(ik)*z
            e2(ik)=sum(conjg(z)*hz)
            !e2=sum(conjg(z)*hz)
            if(rank==0) then
                write(17,*) ' i, e2 ' , i,e2(ik); call flush(17)
            endif
            z=hz
         enddo
         
         e2(ik) = e2(ik) + e1(ik)
         hmax(ik)=max(e1(ik),e2(ik))
         hmin(ik)=min(e1(ik),e2(ik))
         havg=(hmin(ik)+hmax(ik))/2d0
         dh  =(hmax(ik)-hmin(ik))*fctr
         hmax(ik)=havg+dh/2d0
         hmin(ik)=havg-dh/2d0

     enddo
!  end if

!  call bcast_scalar_r8(hmin)
!  call bcast_scalar_r8(hmax)

  deallocate(z,hz,tmp)

contains
  subroutine debug_dhhavg
    use main_mod, only : vks
    implicit none
    if(rank==0) then
       write(6,*)' vks_minmax ',minval(vks),maxval(vks)
       call flush(6)
    endif
  end subroutine debug_dhhavg
end subroutine get_hminmax 

