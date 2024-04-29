subroutine diagh
      use param_mod, only : k_flg
      use mat_module, only : mat_diag
      use main_mod
      use mpi_lib_ours
      implicit none

      complex*16,allocatable :: hm(:,:),vec(:,:)
      complex*16, allocatable :: hp(:,:,:),tmp(:,:,:,:)
      complex*16, allocatable :: hp1(:,:),hptot(:,:),p1(:,:)
      real*8, allocatable :: eig_tot(:),eig(:)
      integer,allocatable :: ord(:)
      integer :: i,j,ik,jk,stat

      if(k_flg) then
        allocate(hp(nx,ny,nz),tmp(nx,ny,nz,nb),stat=stat)
        if(stat/=0) stop 'hp tmp alloc problem diag'
        allocate(hm(nb,nb),vec(nb,nb),eig(nb*nkpt),stat=stat)
        if(stat/=0) stop 'hp tmp1 alloc problem diag'
      else
        allocate(hp1(nn/nodes,nb*nodes),hptot(nn,nb*nodes),&
                p1(nn/nodes,nb*nodes),stat=stat)
        if(stat/=0) stop 'hp tmp alloc problem diag'
        allocate(hm(nb*nodes,nb*nodes),vec(nb*nodes,nb*nodes),&
                eig(nb*nodes),stat=stat)
        if(stat/=0) stop 'hp tmp1 alloc problem diag'
      endif

      if(k_flg) then
        if(rank==0) then
          allocate(eig_tot(nb*nkpt),ord(nb*nkpt),stat=i)
          if(i/=0) stop 'eig_tot alloc problem'
          if(rank==0)write(*,*) 'diagh'
        else
          allocate(eig_tot(1),stat=i)
          if(i/=0) stop 'eig_tot alloc problem'
        endif
        if(rank==0)write(*,*) '1. the hamiltonian matrix'
      endif

      do ik=1,nk_loc
        if(k_flg) then
          do i=1,nb
                !call h_phi(phit(:,:,:,i,ik),hp,ik) !minh comment out
              if(h_type .eq. 0) then !minh 
                call h_phi(phit(:,:,:,i,ik),hp,ik) !minh 
                !call sh(phit(:,:,:,i,ik),hp,ik) !minh 
              else if (h_type .eq. 1) then !minh 
                call shs(phit(:,:,:,i,ik),hp,ik) !minh 
              endif
            do j=1,nb
                hm(i,j)=sum(conjg(phit(:,:,:,j,ik))*hp)*dv
            enddo
          enddo

          if(rank==0)write(*,*) 'check hermiticity'
          if(rank==0)write(*,*) maxval(dble(hm-conjg(transpose(hm)))),&
                maxval(dimag(hm-conjg(transpose(hm))))
          call mat_diag(hm,vec,eig((ik-1)*nb+1:ik*nb))

          tmp=phit(:,:,:,:,ik)
          phit(:,:,:,:,ik)=0d0
          do i=1,nb
            do j=1,nb
              phit(:,:,:,j,ik)=phit(:,:,:,j,ik)+tmp(:,:,:,i)*vec(i,j)
            enddo
          enddo
          mu(ik)=eig((ik-1)*nb+nocc+1)+dmua
          if(rank==0)write(*,*) 'updated mu',rank,ik,mu(ik)
        else
          if(rank==0) then
            do i=1,nb*nodes
              !call h_phi(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh comment out
              if(h_type .eq. 0) then !minh 
                call h_phi(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh 
                !call sh(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh 
              else if (h_type .eq. 1) then !minh 
                call shs(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh 
              endif
!              do j=1,nb*nodes
!                hm(i,j)=sum(conjg(phit_tot(:,:,:,j,ik))*hp)*dv
!              enddo
            enddo
          endif
        
          do i=1,nb*nodes
            call scatter_c16(hptot(:,i),hp1(:,i),size(hp1(:,i)),0)
            call scatter_c16(phit_tot(:,:,:,i,ik),p1(:,i),nn/nodes,0)
          enddo

          do i=1,nb*nodes
            do j=1,nb*nodes
              hm(i,j)=sum(conjg(hp1(:,i))*p1(:,j))*dv
            enddo
          enddo

          call allsum_c16(hm,size(hm))

          if(rank==0) then
            write(*,*) 'flag0'
            write(*,*) 'flag1'
            write(*,*) 'check hermiticity'
            write(*,*) maxval(dble(hm-conjg(transpose(hm)))),&
                maxval(dimag(hm-conjg(transpose(hm))))

            call mat_diag(hm,vec,eig)
            write(*,*) 'eig',eig
            mu(ik)=eig(nocc+1)+dmua
            if(h_type .eq. 0) then
              open(unit=1,file='eig')
            else
              open(unit=1,file='eig_bar')
            endif  
            !open(unit=1,file='eig')
            write(1,*) eig
            close(1)
          endif

          call bcast_c16(vec,size(vec),0)
          call bcast_r8(eig,size(eig),0)
          call bcast_r8(mu,size(mu),0)
          if(rank==0)write(*,*) 'updated mu',rank,ik,mu(ik)

          p1=matmul(p1,vec)
          do i=1,nb*nodes
            call gather_c16(p1(:,i),phit_tot(:,:,:,i,ik),nn/nodes,0)
          enddo

!          call writewf
!          tmp=phit_tot(:,:,:,:,ik)
!          phit_tot(:,:,:,:,ik)=0d0
!          do i=1,nb*nodes
!            do j=1,nb*nodes
!              phit_tot(:,:,:,j,ik)=phit_tot(:,:,:,j,ik)+tmp(:,:,:,i)*vec(i,j)
!            enddo
!          enddo
!        endif
        endif
      enddo 
      
      if(rank==0)write(*,*) 'eigenvalues'

      if(k_flg) then
        call gatherv_r8(eig,eig_tot,size(eig),size(eig_tot),0)
        if(rank==0) then
          call order_r_index(eig_tot,ord,nkpt*nb)
          write(*,*) 'homo,lumo,gap,in Ha',eig_tot(ord(nkpt*nocc))&
              ,eig_tot(ord(nkpt*nocc+1)),eig_tot(ord(nkpt*nocc+1))-eig_tot(ord(nkpt*nocc))
          write(*,*) 'homo,lumo,gap,in eV',eig_tot(ord(nkpt*nocc))*27.211396641308&
              ,eig_tot(ord(nkpt*nocc+1))*27.211396641308,&
              (eig_tot(ord(nkpt*nocc+1))-eig_tot(ord(nkpt*nocc)))*27.211396641308
          deallocate(ord)
        endif
      else
        if(rank==0) then
          write(*,*) eig
          write(*,*) 'homo,lumo,gap,in Ha',eig(nocc),eig(nocc+1),eig(nocc+1)-eig(nocc)
          write(*,*) 'homo,lumo,gap,in eV',eig(nocc)*27.211396641308,eig(nocc+1)*27.211396641308,&
                  (eig(nocc+1)-eig(nocc))*27.211396641308
        endif
        call scatter_c16(phit_tot,phit,size(phit),0)
      endif
!      write(*,*) '3. <psi|shs|psi>'
!      do i=1,nb
!        call h_phi(phit(:,:,:,i),hp)
!!        call model_h(phi(:,:,:,i),hp)
!        do j=1,nb
!            hm(i,j)=sum(phit(:,:,:,j)*hp)*dv
!        enddo
!!        write(*,'(100(1x,f9.5))') hm(i,:)
!        write(*,*) hm(i,:)
!      enddo
      if(k_flg) then
          deallocate(ord,eig_tot,hp,tmp,hm,vec,eig,stat=stat)
      else
          deallocate(hp1,hptot,p1,hm,vec,eig,stat=stat)
      endif
      if(stat/=0) stop 'hp tmp dealloc problem diag'

end subroutine diagh

subroutine diagh_tddft_init
      use param_mod, only : k_flg
      use mat_module, only : mat_diag
      use main_mod
      use tddft_mod, only : nstates, nc_init, nb_init
      use mpi_lib_ours
      implicit none

      complex*16,allocatable :: hm(:,:),vec(:,:)
      complex*16, allocatable :: hp(:,:,:),tmp(:,:,:,:)
      complex*16, allocatable :: hp1(:,:),hptot(:,:),p1(:,:)
      complex*16, allocatable :: ptot(:,:)
      real*8, allocatable :: eig_tot(:),eig(:)
      integer,allocatable :: ord(:)
      integer :: i,j,ik,jk,stat

      allocate(hp1(nn/nc_init,nstates),hptot(nn,nstates),&
              p1(nn/nc_init,nstates),ptot(nn,nstates), stat=stat)
      if(stat/=0) stop 'hp tmp alloc problem diag'
      allocate(hm(nstates,nstates),vec(nstates,nstates),&
              eig(nstates),stat=stat)
      if(stat/=0) stop 'hp tmp1 alloc problem diag'


      do ik=1,nk_loc
          if(rank==0) then
            do i=1,nstates
              if(h_type .eq. 0) then !minh 
                call h_phi(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh 
              else if (h_type .eq. 1) then !minh 
                call shs(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh 
              endif
              ptot(:,i) = reshape(phit_tot(:,:,:,i,ik), (/nn/)) !minh
            enddo
          endif
          
          call sync_mpi
          hp1 = 0
          p1 = 0
          do i=1,nstates
            do j=1,nc_init-1
              if (rank==0) then
                hp1(:,i) = hptot(j*nn/nc_init+1:(j+1)*nn/nc_init,i)
                p1(:,i) = ptot(j*nn/nc_init+1:(j+1)*nn/nc_init,i)
              endif
              call send_receive_c16_array_in_mpi(hp1(:,i), size(hp1(:,i)), 0, j, j)
              call send_receive_c16_array_in_mpi(p1(:,i), size(p1(:,i)), 0, j, j)
            enddo
          enddo
          if (rank==0) then
            do i=1,nstates
              hp1(:,i) = hptot(1:nn/nc_init,i)
              p1(:,i) = ptot(1:nn/nc_init,i)
            enddo
          endif

          do i=1,nstates
            do j=1,nstates
              hm(i,j)=sum(conjg(hp1(:,i))*p1(:,j))*dv
            enddo
          enddo

          call allsum_c16(hm,size(hm))
          if(rank==0) then
            write(*,*) 'flag0'
            write(*,*) 'flag1'
            write(*,*) 'check hermiticity'
            write(*,*) maxval(dble(hm-conjg(transpose(hm)))),&
                maxval(dimag(hm-conjg(transpose(hm))))

            !write(6,*) hm
            call mat_diag(hm,vec,eig)
            !write(*,*) 'eig',eig
            mu(ik)=eig(nocc+1)+dmua
          endif

          !if(rank==0)write(*,*) 'updated mu',rank,ik,mu(ik)
          call bcast_c16(vec,size(vec),0)
          call bcast_r8(eig,size(eig),0)
          call bcast_r8(mu,size(mu),0)
          if(rank==0)write(*,*) 'updated mu',rank,ik,mu(ik)

          p1=matmul(p1,vec)
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
      enddo 
      
      if(rank==0)write(*,*) 'eigenvalues'

      if(rank==0) then
        write(*,*) eig
        write(*,*) 'homo,lumo,gap,in Ha',eig(nocc),eig(nocc+1),eig(nocc+1)-eig(nocc)
        write(*,*) 'homo,lumo,gap,in eV',eig(nocc)*27.211396641308,eig(nocc+1)*27.211396641308,&
                (eig(nocc+1)-eig(nocc))*27.211396641308
      endif

      deallocate(hp1,hptot,p1,hm,vec,eig,stat=stat)
      if(stat/=0) stop 'hp tmp dealloc problem diag'
end subroutine diagh_tddft_init

subroutine diagh_tddft_init_old
      use param_mod, only : k_flg
      use mat_module, only : mat_diag
      use main_mod
      use tddft_mod, only : nstates
      use mpi_lib_ours
      implicit none

      complex*16,allocatable :: hm(:,:),vec(:,:)
      complex*16, allocatable :: hp(:,:,:),tmp(:,:,:,:)
      complex*16, allocatable :: hp1(:,:),hptot(:,:),p1(:,:)
      real*8, allocatable :: eig_tot(:),eig(:)
      integer,allocatable :: ord(:)
      integer :: i,j,ik,jk,stat

      allocate(hp1(nn,nstates),hptot(nn,nstates),&
              p1(nn,nstates),stat=stat)
      if(stat/=0) stop 'hp tmp alloc problem diag'
      allocate(hm(nstates,nstates),vec(nstates,nstates),&
              eig(nstates),stat=stat)
      if(stat/=0) stop 'hp tmp1 alloc problem diag'


      do ik=1,nk_loc
          !if(rank==0) then
            do i=1,nstates
              !call h_phi(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh comment out
              if(h_type .eq. 0) then !minh 
                call h_phi(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh 
                !call sh(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh 
              else if (h_type .eq. 1) then !minh 
                call shs(phit_tot(:,:,:,i,ik),hptot(:,i),ik) !minh 
              endif
!              enddo
            enddo
          !endif

          do i=1,nstates
            do j=1,nstates
              p1(:,j) = reshape(phit_tot(:,:,:,j,ik), (/nn/))
              hm(i,j)=sum(conjg(hptot(:,i))*p1(:,j))*dv
            enddo
          enddo

          !if(rank==0) then
            write(*,*) 'flag0'
            write(*,*) 'flag1'
            write(*,*) 'check hermiticity'
            write(*,*) maxval(dble(hm-conjg(transpose(hm)))),&
                maxval(dimag(hm-conjg(transpose(hm))))

            call mat_diag(hm,vec,eig)
            write(*,*) 'eig',eig
            mu(ik)=eig(nocc+1)+dmua
            !if(h_type .eq. 0) then
            !  open(unit=1,file='eig')
            !else
            !  open(unit=1,file='eig_bar')
            !endif  
            !open(unit=1,file='eig')
            !write(1,*) eig
            !close(1)
          !endif

          !if(rank==0)write(*,*) 'updated mu',rank,ik,mu(ik)
          write(*,*) 'updated mu',rank,ik,mu(ik)

          p1=matmul(p1,vec)
          do i=1,nstates
            phit_tot(:,:,:,i,ik) = reshape(p1(:,i), (/nx,ny,nz/))
          enddo
      enddo 
      
      !if(rank==0)write(*,*) 'eigenvalues'
      write(*,*) 'eigenvalues'

      !if(rank==0) then
        write(*,*) eig
        write(*,*) 'homo,lumo,gap,in Ha',eig(nocc),eig(nocc+1),eig(nocc+1)-eig(nocc)
        write(*,*) 'homo,lumo,gap,in eV',eig(nocc)*27.211396641308,eig(nocc+1)*27.211396641308,&
                (eig(nocc+1)-eig(nocc))*27.211396641308
      !endif

      deallocate(hp1,hptot,p1,hm,vec,eig,stat=stat)
      if(stat/=0) stop 'hp tmp dealloc problem diag'
end subroutine diagh_tddft_init_old

