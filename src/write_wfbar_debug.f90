subroutine writewfbar_debug
   use main_mod
   use param_mod
   use mpi_lib_ours

   implicit none
  
   integer :: i,ix,iy,iz
   real*8 :: eig(nb*nodes)
   real*8 :: tmp(nx,ny,nz)
   complex*16 :: sp(nx,ny,nz,nb)

   if(k_flg) stop 'not printing wf for kpoint'

   call scatter_c16(phit_tot,phit,size(phit),0)
   do i=1,nb
      call sn_phi(phit(:,:,:,i,1),sp(:,:,:,i),1,0.5d0)
   enddo
   call gather_c16(sp,phit_tot,size(phit),0)

!   write(*,*) 'sum dens', sum(abs(phit_tot(:,:,:,1:nocc,1)**2))*dv

   if(flg_bin) then
     if(rank==0) then

       open(unit=1,file='eig')
       read(1,*) eig
       close(1)

       open(unit=1,file='wfbar_debug.bin')
       write(1) 'nx      ',nx
       write(1) 'ny      ',ny
       write(1) 'nz      ',nz
       write(1) 'dx      ',dx
       write(1) 'dy      ',dy
       write(1) 'dz      ',dz
       write(1) 'nsp     ',1
       write(1) 'nstates ',nb*nodes
       write(1) 'evls'
       write(1) eig
       write(1) 'orbitals'
       do i=1,nb*nodes
         write(1) i,'1'
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

       open(unit=1,file='eig')
       read(1,*) eig
       close(1)

       open(unit=1,file='wfbar_debug.txt')
       write(1,*) 'nx      ',nx
       write(1,*) 'ny      ',ny
       write(1,*) 'nz      ',nz
       write(1,*) 'dx      ',dx
       write(1,*) 'dy      ',dy
       write(1,*) 'dz      ',dz
       write(1,*) 'nsp     ',1
       write(1,*) 'nstates ',nb*nodes
       write(1,*) 'evls'
       write(1,*) eig
       write(1,*) 'orbitals'
       do i=1,nb*nodes
         write(1,*) i,'1'
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
       enddo
       close(1)
     endif
   endif
end subroutine writewfbar_debug
