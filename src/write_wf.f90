subroutine writewf
   use main_mod
   use param_mod
   use mpi_lib_ours

   implicit none
  
   integer :: i,ix,iy,iz
   real*8 :: eig(nb*nodes)
   real*8 :: tmp(nx,ny,nz)

   if(k_flg) stop 'not printing wf for kpoint'
  
  
   !call update_dens
   !call get_pot1(dens)
   !if (rank==0) then
   !  write(*,*) 'ncoret',sum(ncoret**3)
   !  write(*,*) 'vks',sum(vks**3)
   !  write(*,*) 'vxc',sum(vxc**3)
   !endif
!   write(*,*) 'sum dens', sum(abs(phit_tot(:,:,:,1:nocc,1)**2))*dv

   if(flg_bin) then
     if(rank==0) then

       if(h_type .eq. 0) then
         open(unit=1,file='eig')
       else
         open(unit=1,file='eig_bar')
       endif  
       !open(unit=1,file='eig')
       read(1,*) eig
       close(1)
       
       if(h_type .eq. 0) then
         open(unit=1,file='wf.bin',form='unformatted')
       else
         open(unit=1,file='wf_bar.bin',form='unformatted')
       endif  
       write(1) 'nx      ',nx
       write(1) 'ny      ',ny
       write(1) 'nz      ',nz
       write(1) 'dx      ',dx
       write(1) 'dy      ',dy
       write(1) 'dz      ',dz
       write(1) 'nsp     ',1
       write(1) 'nstates ',nb*nodes
       write(1) 'orbmn_wf',1
       write(1) 'ihomo   ',nocc
       write(1) 'evls    '
       write(1) eig
       write(1) 'orbitals'
       do i=1,nb*nodes
         write(1) i,1
         do iz=1,nz
           do iy=1,ny
             do ix=1,nx
               if (dble(phit_tot(ix,iy,iz,i,1))>0d0) then
                 tmp(ix,iy,iz)=abs(phit_tot(ix,iy,iz,i,1))
               else
                 tmp(ix,iy,iz)=-1d0*abs(phit_tot(ix,iy,iz,i,1))
               endif
             enddo
           enddo
         enddo
         write(1) tmp 
       enddo
       close(1)
     endif
   else
     if(rank==0) then

       !open(unit=1,file='eig')
       if(h_type .eq. 0) then
         open(unit=1,file='eig')
       else
         open(unit=1,file='eig_bar')
       endif  
       read(1,*) eig
       close(1)

       if(h_type .eq. 0) then
         open(unit=1,file='wf.txt')
       else
         open(unit=1,file='wf_bar.txt')
       endif  
       write(1,*) 'nx      ',nx
       write(1,*) 'ny      ',ny
       write(1,*) 'nz      ',nz
       write(1,*) 'dx      ',dx
       write(1,*) 'dy      ',dy
       write(1,*) 'dz      ',dz
       write(1,*) 'nsp     ',1
       write(1,*) 'nstates ',nb*nodes
       write(1,*) 'orbmn_wf',1
       write(1,*) 'ihomo   ',nocc
       write(1,*) 'evls    ' 
       write(1,*) eig
       write(1,*) 'orbitals'
       do i=1,nb*nodes
         write(1,*) i,1
         do iz=1,nz
           do iy=1,ny
             do ix=1,nx
               if (dble(phit_tot(ix,iy,iz,i,1))>0d0) then
                 tmp(ix,iy,iz)=abs(phit_tot(ix,iy,iz,i,1))
               else
                 tmp(ix,iy,iz)=-1d0*abs(phit_tot(ix,iy,iz,i,1))
               endif
             enddo
           enddo
         enddo
         write(1,*) tmp 
         !if(rank==0 ) write(6,*) 'is, norm', i, sum(abs(tmp)**2d0)*dv
       enddo
       close(1)
     endif
   endif
end subroutine writewf      

subroutine write_wf_for_sgw_debug(eig,ns)
   use main_mod
   use param_mod
   use mpi_lib_ours

   implicit none
  
   integer :: i,ix,iy,iz,ns
   real*8 :: eig(ns)
   real*8 :: tmp(nx,ny,nz)

   if(k_flg) stop 'not printing wf for kpoint'
  
  
   !call update_dens
   !call get_pot1(dens)
   !if (rank==0) then
   !  write(*,*) 'ncoret',sum(ncoret**3)
   !  write(*,*) 'vks',sum(vks**3)
   !  write(*,*) 'vxc',sum(vxc**3)
   !endif
!   write(*,*) 'sum dens', sum(abs(phit_tot(:,:,:,1:nocc,1)**2))*dv

   if(flg_bin) then
     if(rank==0) then

       if(h_type .eq. 0) then
         open(unit=1,file='wf_sgw.bin',form='unformatted')
       else
         open(unit=1,file='wf_bar_sgw.bin',form='unformatted')
       endif  
       write(1) 'nx      ',nx
       write(1) 'ny      ',ny
       write(1) 'nz      ',nz
       write(1) 'dx      ',dx
       write(1) 'dy      ',dy
       write(1) 'dz      ',dz
       write(1) 'nsp     ',1
       write(1) 'nstates ',ns
       write(1) 'evls    '
       write(1) eig
       write(1) 'orbitals'
       do i=1,ns
         write(1) i,1
         do iz=1,nz
           do iy=1,ny
             do ix=1,nx
               if (dble(phit_tot(ix,iy,iz,i,1))>0d0) then
                 tmp(ix,iy,iz)=abs(phit_tot(ix,iy,iz,i,1))
               else
                 tmp(ix,iy,iz)=-1d0*abs(phit_tot(ix,iy,iz,i,1))
               endif
             enddo
           enddo
         enddo
         write(1) tmp 
       enddo
       close(1)
     endif
   else
     if(rank==0) then

       if(h_type .eq. 0) then
         open(unit=1,file='wf_sgw.txt')
       else
         open(unit=1,file='wf_bar_sgw.txt')
       endif  
       write(1,*) 'nx      ',nx
       write(1,*) 'ny      ',ny
       write(1,*) 'nz      ',nz
       write(1,*) 'dx      ',dx
       write(1,*) 'dy      ',dy
       write(1,*) 'dz      ',dz
       write(1,*) 'nsp     ',1
       write(1,*) 'nstates ',ns
       write(1,*) 'evls    '
       write(1,*) eig
       write(1,*) 'orbitals'
       do i=1,ns
         write(1,*) i,1
         do iz=1,nz
           do iy=1,ny
             do ix=1,nx
               if (dble(phit_tot(ix,iy,iz,i,1))>0d0) then
                 tmp(ix,iy,iz)=abs(phit_tot(ix,iy,iz,i,1))
               else
                 tmp(ix,iy,iz)=-1d0*abs(phit_tot(ix,iy,iz,i,1))
               endif
             enddo
           enddo
         enddo
         write(1,*) tmp 
         !if(rank==0 ) write(6,*) 'is, norm', i, sum(abs(tmp)**2d0)*dv
       enddo
       close(1)
     endif
   endif
end subroutine write_wf_for_sgw_debug
