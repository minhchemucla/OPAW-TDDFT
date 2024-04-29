subroutine chebyr_m(hmax,hbot)
  use param_mod, only : nch
  use main_mod, only : phit,nb,nx,ny,nz,dv,nk_loc, h_type
  use mpi_lib_ours
  implicit none
  integer :: st, k, ib, ik, ip
  real*8 hmax(nk_loc), hbot(nk_loc), havdum, dhdum
  complex*16, allocatable, dimension(:,:,:,:) :: ca, cb, cc
  

  allocate(ca(nx,ny,nz,nb), cb(nx,ny,nz,nb),&
          cc(nx,ny,nz,nb), stat=st)
  call check0(st,' ca,cb,cc ')


  !ip = 2000+rank
  !write(ip,*) 
  !write(*,*) 'rank, nb', rank, nb

  do ik=1,nk_loc
      ca = phit(:,:,:,:,ik)
      cb = phit(:,:,:,:,ik)
      havdum= (hmax(ik)+hbot(ik))/2d0
      dhdum = (hmax(ik)-hbot(ik))/2d0
      do k=1,nch
         do ib=1,nb
           if(h_type .eq. 0) then
                call sh(cb(:,:,:,ib),cc(:,:,:,ib),ik)
           else if (h_type .eq. 1) then
                call shs(cb(:,:,:,ib),cc(:,:,:,ib),ik)
           else
             write(*,*) 'h_type only sh (=0) or shs (=1)'
             stop
           endif
         enddo

         if(k>1) cc= 2d0/dhdum * (cc - havdum*cb) - ca
         if(k==1)cc= 1d0/dhdum * (cc - havdum*cb)
     
         ca=cb
         cb=cc
         call time_print(' post_scaling ')
      enddo
      phit(:,:,:,:,ik) = cc
  enddo
  !write(6,*) 'rank, phit_cheby test minh:', rank, phit(1,1,1,1,1)
  deallocate(ca,cb,cc)
end subroutine chebyr_m

subroutine chebyr_m_tddft_init(hmax,hbot)
  use param_mod, only : nch
  use main_mod, only : phit,nb,nx,ny,nz,dv,nk_loc, h_type
  use tddft_mod, only : nstates, nc_init, nb_init
  use mpi_lib_ours
  implicit none
  integer :: st, k, ib, ik
  real*8 hmax(nk_loc), hbot(nk_loc), havdum, dhdum
  complex*16, allocatable, dimension(:,:,:,:) :: ca, cb, cc
  


  allocate(ca(nx,ny,nz,nb_init), cb(nx,ny,nz,nb_init),&
          cc(nx,ny,nz,nb_init), stat=st)
  call check0(st,' ca,cb,cc ')

  if (rank<nc_init) then
    do ik=1,nk_loc
      ca = phit(:,:,:,:,ik)
      cb = phit(:,:,:,:,ik)
      havdum= (hmax(ik)+hbot(ik))/2d0
      dhdum = (hmax(ik)-hbot(ik))/2d0
      do k=1,nch
        do ib=1,nb_init
          if(h_type .eq. 0) then
            call sh(cb(:,:,:,ib),cc(:,:,:,ib),ik)
          else if (h_type .eq. 1) then
            call shs(cb(:,:,:,ib),cc(:,:,:,ib),ik)
          else
            write(*,*) 'h_type only sh (=0) or shs (=1)'
            stop
          endif
        enddo

        if(k>1) cc= 2d0/dhdum * (cc - havdum*cb) - ca
        if(k==1)cc= 1d0/dhdum * (cc - havdum*cb)
     
        ca=cb
        cb=cc
        !call time_print(' post_scaling ')
      enddo
      phit(:,:,:,:,ik) = cc
    enddo
    deallocate(ca,cb,cc)
  endif

  call sync_mpi

end subroutine chebyr_m_tddft_init

subroutine chebyr_m_tddft_init_old(hmax,hbot)
  use param_mod, only : nch
  use main_mod, only : phit_tot,nb,nx,ny,nz,dv,nk_loc, h_type
  use tddft_mod, only : nstates
  use mpi_lib_ours
  implicit none
  integer :: st, k, ib, ik
  real*8 hmax(nk_loc), hbot(nk_loc), havdum, dhdum
  complex*16, allocatable, dimension(:,:,:,:) :: ca, cb, cc
  


  allocate(ca(nx,ny,nz,nstates), cb(nx,ny,nz,nstates),&
          cc(nx,ny,nz,nstates), stat=st)
  call check0(st,' ca,cb,cc ')

  do ik=1,nk_loc
      ca = phit_tot(:,:,:,:,ik)
      cb = phit_tot(:,:,:,:,ik)
      havdum= (hmax(ik)+hbot(ik))/2d0
      dhdum = (hmax(ik)-hbot(ik))/2d0
      do k=1,nch
         do ib=1,nstates
           if(h_type .eq. 0) then
                call sh(cb(:,:,:,ib),cc(:,:,:,ib),ik)
           else if (h_type .eq. 1) then
                call shs(cb(:,:,:,ib),cc(:,:,:,ib),ik)
           else
             write(*,*) 'h_type only sh (=0) or shs (=1)'
             stop
           endif
         enddo

         if(k>1) cc= 2d0/dhdum * (cc - havdum*cb) - ca
         if(k==1)cc= 1d0/dhdum * (cc - havdum*cb)
     
         ca=cb
         cb=cc
         call time_print(' post_scaling ')
      enddo
      phit_tot(:,:,:,:,ik) = cc
  enddo
  !write(6,*) 'rank, phit_cheby test minh:', rank, phit(1,1,1,1,1)
  deallocate(ca,cb,cc)

end subroutine chebyr_m_tddft_init_old

