subroutine update_dens
      use param_mod, only : k_flg
      use main_mod
      use atom_mod, only : natom, at => atominfo, atom_map
      use mpi_lib_ours
!                at_old => atominfo_old
      implicit none

      integer :: i,ia,it
      real*8  :: sp(nn)
      real*8  :: n1,n2,n3

!      call get_sij
      if(k_flg) then
         call acc_coeff !calculates rho_ij for each atom. See Eq. after Eq. (8) 
         call acc_dens(dens,nx,ny,nz) !calculates valence density from wf. Eq. (7)
      else
         if(rank==0)call acc_coeff1
         if(rank==0)call acc_dens1(dens,nx,ny,nz)
      endif

      if(rank==0) write(*,*) 'before allsum,rank,dens',rank,sum(dens)*dv
      if(k_flg) then
        call allsum_r8(dens,size(dens))
        do ia=1,natom
          call allsum_r8(at(ia)%rhoij,size(at(ia)%rhoij))
        enddo
!      else
!        call bcast_r8(dens,size(dens),0)
!        do ia=1,natom
!          call bcast_r8(at(ia)%rhoij,size(at(ia)%rhoij),0)
!        enddo
      endif

      !write(1800+rank,*) 'rhoij'
      !do ia=1,natom
      !    do i=1,size(at(ia)%rhoij,1)
      !        write(1800+rank,*) at(ia)%rhoij(i,i)
      !    enddo
      !enddo
!      dens=dens*0.1d0+dens_p*0.9d0
!      dens_p=dens
!
!      do ia=1,natom
!        at(ia)%rhoij=at(ia)%rhoij*0.1+at_old(ia)%rhoij*0.9
!        at_old(ia)%rhoij=at(ia)%rhoij
!      enddo

!      if(k_flg) then
        if(rank==0) then 
          call get_nhat
        endif
        call bcast_r8(nhat,size(nhat),0) 
!      else
!        call get_nhat
!      endif
      !if(rank==0) write(*,*) 'max,min(dens)',rank,maxval(dens), minval(dens)
      !if(rank==0) write(*,*) 'max,min(nhat)',rank,maxval(nhat), minval(nhat)
contains
  subroutine acc_dens1(rho,nx,ny,nz)
    implicit none
    integer :: nx,ny,nz,ik,jk,iib
    real*8  :: rho(nx,ny,nz),rhok(nx,ny,nz)

    rho=0d0
    do ik=1,nk_loc
        jk=(ik-1)*nodes+rank+1
        rhok=0d0
!        do i=1,nb
!            iib=rank*nb+i
!            if(iib.le.nocc) rhok=rhok+abs(phit(:,:,:,i,ik))**2
!        enddo
        do i=1,nocc
           rhok=rhok+abs(phit_tot(:,:,:,i,ik))**2
        enddo
        rhok=rhok*wk(jk)
        rho=rho+rhok
    enddo
    rho=rho*2d0
  end subroutine acc_dens1

  subroutine acc_dens(rho,nx,ny,nz)
      implicit none
      integer :: nx,ny,nz,ik,jk
      real*8  :: rho(nx,ny,nz),rhok(nx,ny,nz)

      rho=0d0
      do ik=1,nk_loc
          jk=(ik-1)*nodes+rank+1
          rhok=0d0
          do i=1,nocc
              rhok=rhok+abs(phit(:,:,:,i,ik))**2
          enddo
          write(*,*) 'jk,rhok',jk,sum(rhok)*dv
          rhok=rhok*wk(jk)
          rho=rho+rhok
      enddo
      rho=rho*2d0
  end subroutine acc_dens
end subroutine update_dens      
