subroutine orthog_phi
      use param_mod, only : k_flg
      use mat_module, only : mat_diag
      use main_mod
      use mpi_lib_ours
      implicit none
      
      integer :: i,j,ik,stat
      real*8  :: norm
      real*8,allocatable  :: eig(:)
      complex*16,allocatable :: hm(:,:),vec(:,:)
      complex*16,allocatable :: tmp(:,:,:,:),tmp1(:,:,:)
      complex*16,allocatable :: sp1(:,:),sptot(:,:),p1(:,:),pt1(:,:)
      
      if(stat/=0) stop 'tmp1 allc problem orthogphi'
      if(k_flg) then
          allocate(tmp(nx,ny,nz,nb),hm(nb,nb),vec(nb,nb),eig(nb),stat=stat)
          allocate(tmp1(nx,ny,nz),stat=stat)
      else
          allocate(hm(nb*nodes,nb*nodes),&
                  vec(nb*nodes,nb*nodes),eig(nb*nodes),stat=stat)
          if(stat/=0) stop 'tmp111 allc problem orthogphi'
          allocate(sp1(nn/nodes,nb*nodes),sptot(nn,nb*nodes),&
                  p1(nn/nodes,nb*nodes),pt1(nn/nodes,nb*nodes),stat=stat) 
      endif
      if(stat/=0) stop 'tmp allc problem orthogphi'


      do ik=1,nk_loc
          call sync_mpi
          if(k_flg) then
            do i=1,nb
              if (h_type .eq. 0) then
                call sn_phi(phit(:,:,:,i,ik),tmp1,ik,1d0) !minh 
              else if (h_type .eq. 1) then
                tmp1 =  phit(:,:,:,i,ik) !minh
              endif
              do j=1,nb
                hm(i,j)=sum(conjg(tmp1)*phit(:,:,:,j,ik))*dv
              enddo
            enddo
            call mat_diag(hm,vec,eig)
          else
            call gather_c16(phit,phit_tot,size(phit),0)
            if(rank==0) then
              do i=1,nb*nodes
                if (h_type .eq. 0) then
                  call sn_phi(phit_tot(:,:,:,i,ik),sptot(:,i),ik,1d0) !minh 
                else if (h_type .eq. 1) then
                  sptot(:,i) = reshape(phit_tot(:,:,:,i,ik), (/nn/)) !minh
                endif
!                do j=1,nb*nodes
!                  hm(i,j)=sum(conjg(tmp1)*phit_tot(:,:,:,j,ik))*dv
!                enddo
              enddo
!              call mat_diag(hm,vec,eig)
            endif

            do i=1,nb*nodes
              call scatter_c16(sptot(:,i),sp1(:,i),size(sp1(:,i)),0)
              call scatter_c16(phit_tot(:,:,:,i,ik),p1(:,i),nn/nodes,0)
            enddo
!            write(*,*) 'check scatter'
!            write(*,*) 'sptot1-5',sptot(1:5,2),sptot(nn/2+1:nn/2+5,2)
!            write(*,*) 'sp11-5',sp1(1:5,2)

!            if(rank==0) then
!              allocate(tmp2(nn))
!              tmp2=reshape(phit_tot(:,:,:,1,1),(/nn/))
!            if(rank==0) write(*,*) 'phit_tot1-5',&
!                   phit_tot(1:5,1,1,25,1),phit_tot(1:5,1,nz/2+1,25,1)
!            endif
!            write(*,*) 'p11-5',p1(1:5,25)
!            if(rank==0) deallocate(tmp2)

            do i=1,nb*nodes
              do j=1,nb*nodes
                  hm(i,j)=sum(conjg(sp1(:,i))*p1(:,j))*dv
              enddo
              !if(rank==0)write(6,*) hm(i,:)
            enddo
!            call sync_mpi
!
!            do i=1,nb*nodes
!              write(7000+rank,*) hm(i,:)
!            enddo
            call allsum_c16(hm,size(hm)) 
!            do i=1,nb*nodes
!              write(7100+rank,*) hm(i,:)
!            enddo
            if(rank==0) call mat_diag(hm,vec,eig)
                
!            if(rank==0) write(*,*) 'eig',eig
!            call sync_mpi
!            stop

            call bcast_c16(vec,size(vec),0)
            call bcast_r8(eig,size(eig),0)

!            pt1=p1
!            p1=0d0
!            do i=1,nb*nodes
!              do j=1,nb*nodes
!                p1(:,j)=p1(:,j)+pt1(:,i)*vec(i,j)
!              enddo
!            enddo
            p1=matmul(p1,vec)
            do i=1,nb*nodes
              p1(:,i)=p1(:,i)/sqrt(eig(i))
              call gather_c16(p1(:,i),phit_tot(:,:,:,i,ik),nn/nodes,0)
            enddo

          endif

          if(k_flg) then
            tmp=phit(:,:,:,:,ik)
            phit(:,:,:,:,ik)=0d0
            do i=1,nb
              do j=1,nb
                 phit(:,:,:,j,ik)=phit(:,:,:,j,ik)+tmp(:,:,:,i)*vec(i,j)
              enddo
            enddo

            do i=1,nb
              if (h_type .eq. 0) then  !minh
                call sn_phi(phit(:,:,:,i,ik),tmp1,ik,1d0) !minh
              else if (h_type .eq. 1) then   !minh
                tmp1 =  phit(:,:,:,i,ik)  !minh
              endif
              norm=sum(conjg(phit(:,:,:,i,ik))*tmp1)*dv
              phit(:,:,:,i,ik)=phit(:,:,:,i,ik)/sqrt(norm)
            enddo
!          else
!            if(rank==0) then
!              allocate(tmp1(nx,ny,nz))
!              tmp=phit_tot(:,:,:,:,ik)
!              phit_tot(:,:,:,:,ik)=0d0
!              do i=1,nb*nodes
!                do j=1,nb*nodes
!                   phit_tot(:,:,:,j,ik)=phit_tot(:,:,:,j,ik)+tmp(:,:,:,i)*vec(i,j)
!                enddo
!              enddo

!              do i=1,nb*nodes
!                call sn_phi(phit_tot(:,:,:,i,ik),tmp1,ik,1d0)
!                do j=1,nb*nodes
!                  norm=sum(conjg(phit_tot(:,:,:,j,ik))*tmp1)*dv
!                  write(720,*) 'i,j,norm',i,j,norm
!                enddo
!!                phit_tot(:,:,:,i,ik)=phit_tot(:,:,:,i,ik)/sqrt(norm)
!              enddo
!
!              deallocate(tmp1)
!            endif
!            call sync_mpi
!            stop
          endif      
      enddo

      !MINH TESTING START PLEASE REMOVE LATER
      !  ik = 1
      !  if(rank==0) then
      !    do i=1,nb*nodes
      !      if (h_type .eq. 0) then
      !        call sn_phi(phit_tot(:,:,:,i,ik),sptot(:,i),ik,1d0) !minh 
      !      else if (h_type .eq. 1) then
      !        sptot(:,i) = reshape(phit_tot(:,:,:,i,ik), (/nn/)) !minh
      !      endif
      !    enddo
      !  endif

      !  do i=1,nb*nodes
      !    call scatter_c16(sptot(:,i),sp1(:,i),size(sp1(:,i)),0)
      !    call scatter_c16(phit_tot(:,:,:,i,ik),p1(:,i),nn/nodes,0)
      !  enddo

      !  do i=1,nb*nodes
      !    do j=1,nb*nodes
      !        hm(i,j)=sum(conjg(sp1(:,i))*p1(:,j))*dv
      !    enddo
      !  enddo
      !  call allsum_c16(hm,size(hm)) 
      !  if(rank==0)write(6,*) 'check orthogonality of <phi|S|phi>'
      !  do i=1,nb*nodes
      !    if(rank==0)write(6,*) hm(i,:)
      !  enddo
      !MINH TESTING END PLEASE REMOVE
      if(k_flg) then
           deallocate(tmp1,tmp,hm,vec,eig,stat=stat)
      else
           deallocate(hm,vec,eig,sp1,sptot,p1,stat=stat)
      endif
      if(stat/=0) stop 'tmp deallc problem orthogphi'
end subroutine      

subroutine orthog_phi_tddft_init
  use param_mod, only : k_flg
  use mat_module, only : mat_diag
  use main_mod
  use tddft_mod, only: nstates, nc_init, nb_init
  use mpi_lib_ours
  implicit none
  
  integer :: i,j,ik,stat
  real*8  :: norm
  real*8,allocatable  :: eig(:)
  complex*16,allocatable :: hm(:,:),vec(:,:)
  complex*16,allocatable :: tmp(:,:,:,:),tmp1(:,:,:)
  complex*16,allocatable :: sp1(:,:),sptot(:,:),p1(:,:)
  complex*16, allocatable :: tmp2(:), ptot(:,:)
  
  if(stat/=0) stop 'tmp1 allc problem orthogphi'
  if(k_flg) then
      allocate(tmp(nx,ny,nz,nb_init),hm(nstates,nb_init),vec(nstates,nb_init),eig(nb_init),stat=stat)
      allocate(tmp1(nx,ny,nz),stat=stat)
  else
      allocate(hm(nstates,nstates),&
              vec(nstates,nstates),eig(nstates),stat=stat)
      if(stat/=0) stop 'tmp111 allc problem orthogphi'
      allocate(sp1(nn/nc_init,nstates),sptot(nn,nstates),&
              p1(nn/nc_init,nstates),ptot(nn,nstates), stat=stat) 
  endif
  if(stat/=0) stop 'tmp allc problem orthogphi'

  ik=1
  call sync_mpi
  if(rank==0) then
    do i=1,nb_init
      phit_tot(:,:,:,i,1)=phit(:,:,:,i,1)
    enddo
  endif
  do i=1,nc_init-1
    do j=1,nb_init
      call send_receive_c16_array_in_mpi(phit(:,:,:,j,ik), size(phit(:,:,:,j,ik)), i, 0, i)
      if(rank==0)then
        phit_tot(:,:,:,i*nb_init+j,ik)=phit(:,:,:,j,ik)
      endif
    enddo
  enddo


  if(rank==0) then
    do i=1,nstates
      if (h_type .eq. 0) then
        call sn_phi(phit_tot(:,:,:,i,ik),sptot(:,i),ik,1d0) !minh 
      else if (h_type .eq. 1) then
        sptot(:,i) = reshape(phit_tot(:,:,:,i,ik), (/nn/)) !minh
      endif
      ptot(:,i) = reshape(phit_tot(:,:,:,i,ik), (/nn/)) !minh
    enddo
  endif

  call sync_mpi

  do i=1,nstates
    do j=1,nc_init-1
      if (rank==0) then
        sp1(:,i) = sptot(j*nn/nc_init+1:(j+1)*nn/nc_init,i)
        p1(:,i) = ptot(j*nn/nc_init+1:(j+1)*nn/nc_init,i)
      endif
      call send_receive_c16_array_in_mpi(sp1(:,i), size(sp1(:,i)), 0, j, j)
      call send_receive_c16_array_in_mpi(p1(:,i), size(p1(:,i)), 0, j, j)
    enddo
  enddo
  if(rank .gt. nc_init-1) then
    sp1 = 0
    p1 = 0
  endif
!stop 'test 1'
  if (rank==0) then
    do i=1,nstates
      sp1(:,i) = sptot(1:nn/nc_init,i)
      p1(:,i) = ptot(1:nn/nc_init,i)
    enddo
  endif

  call sync_mpi

  do i=1,nstates
    do j=1,nstates
      hm(i,j)=sum(conjg(sp1(:,i))*p1(:,j))*dv
    enddo
  enddo

  call allsum_c16(hm,size(hm)) 
  if(rank==0) call mat_diag(hm,vec,eig)
  if(rank==0) write(6,*) eig
  if (rank==0 .and. any(eig < 0)) stop 'negative eigenvalues orthog_phi'
  call bcast_c16(vec,size(vec),0)
  call bcast_r8(eig,size(eig),0)

  p1=matmul(p1,vec)
  do i=1,nstates
    p1(:,i)=p1(:,i)/sqrt(eig(i))
  enddo

  if (rank==0) then
    do i=1,nstates
      ptot(1:nn/nc_init,i) = p1(:,i)
    enddo
  endif
  do i=1,nstates
    do j=1,nc_init-1
      call send_receive_c16_array_in_mpi(p1(:,i), size(p1(:,i)), j, 0, j)
      if (rank==0) then
        ptot(j*nn/nc_init+1:(j+1)*nn/nc_init,i) = p1(:,i)
      endif
    enddo
  enddo

  do i=1,nstates
    if (rank==0) phit_tot(:,:,:,i,ik) = reshape(ptot(:,i), (/nx,ny,nz/))
  enddo


  if(rank==0) then
    do i=1,nstates
      if (h_type .eq. 0) then
        call sn_phi(phit_tot(:,:,:,i,ik),sptot(:,i),ik,1d0) !minh 
      else if (h_type .eq. 1) then
        sptot(:,i) = reshape(phit_tot(:,:,:,i,ik), (/nn/)) !minh
      endif
      ptot(:,i) = reshape(phit_tot(:,:,:,i,ik), (/nn/)) !minh
    enddo
  endif



  !if (rank==0) then
  !  do i=1,nstates
  !    write(6,*) i, sum(conjg(ptot(:,i))*sptot(:,i))*dv, sum(ptot(:,i))*dv
  !  enddo
  !endif


  deallocate(hm,vec,eig,sp1,sptot,p1,stat=stat)
  if(stat/=0) stop 'tmp deallc problem orthogphi'
end subroutine      

subroutine orthog_phi_tddft_init_old
      use param_mod, only : k_flg
      use mat_module, only : mat_diag
      use main_mod
      use tddft_mod, only: nstates
      use mpi_lib_ours
      implicit none
      
      integer :: i,j,ik,stat
      real*8  :: norm
      real*8,allocatable  :: eig(:)
      complex*16,allocatable :: hm(:,:),vec(:,:)
      complex*16,allocatable :: tmp(:,:,:,:),tmp1(:,:,:)
      complex*16,allocatable :: sp1(:,:),sptot(:,:),p1(:,:)
      complex*16, allocatable :: tmp2(:)
      
      if(stat/=0) stop 'tmp1 allc problem orthogphi'
      if(k_flg) then
          allocate(tmp(nx,ny,nz,nstates),hm(nstates,nstates),vec(nstates,nstates),eig(nstates),stat=stat)
          allocate(tmp1(nx,ny,nz),stat=stat)
      else
          allocate(hm(nstates,nstates),&
                  vec(nstates,nstates),eig(nstates),stat=stat)
          if(stat/=0) stop 'tmp111 allc problem orthogphi'
          allocate(sp1(nn,nstates),sptot(nn,nstates),&
                  p1(nn,nstates),stat=stat) 
      endif
      if(stat/=0) stop 'tmp allc problem orthogphi'


      do ik=1,nk_loc
          !call sync_mpi
          if(rank==0) then
            do i=1,nstates
              if (h_type .eq. 0) then
                call sn_phi(phit_tot(:,:,:,i,ik),sptot(:,i),ik,1d0) !minh 
              else if (h_type .eq. 1) then
                sptot(:,i) = reshape(phit_tot(:,:,:,i,ik), (/nn/)) !minh
              endif
            enddo
!            call mat_diag(hm,vec,eig)
          endif



          do i=1,nstates
            do j=1,nstates
                p1(:,j) = reshape(phit_tot(:,:,:,j,ik), (/nn/))
                hm(i,j)=sum(conjg(sptot(:,i))*p1(:,j))*dv
            enddo
            !if(rank==0)write(6,*) hm(i,:)
          enddo

          if(rank==0) call mat_diag(hm,vec,eig)
              

          p1=matmul(p1,vec)
          do i=1,nstates
            p1(:,i)=p1(:,i)/sqrt(eig(i))
            phit_tot(:,:,:,i,ik) = reshape(p1(:,i), (/nx,ny,nz/))
          enddo
      enddo

      deallocate(hm,vec,eig,sp1,sptot,p1,stat=stat)
      if(stat/=0) stop 'tmp deallc problem orthogphi'
end subroutine      
