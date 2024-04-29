subroutine test_eig
  use main_mod
  use mpi_lib_ours
  use tddft_mod
  implicit none
  complex*16 :: tmp(nx,ny,nz), tmp2(nx,ny,nz), tmp3(nx,ny,nz)
  real*8  :: mean, std, tmp4(nx,ny,nz)
  integer :: is, ns
  real*8, allocatable :: eig(:)

  if (tddft_flag > -1) then
    ns = nstates
  else
    ns = nb*nodes
  endif

  allocate(eig(ns))

  call ncpaw_make_hamiltonian

  if(h_type .eq. 0) then
    open(unit=1,file='eig')
  else
    open(unit=1,file='eig_bar')
  endif  
  !open(unit=1,file='eig')
  read(1,*) eig
  close(1)

  write(120,*) 'Hpsi/psi'
  write(121,*) 'S^(-1)Hpsi/psi'
  write(122,*) 'Hpsibar/psibar    psibar = S^{1/2}psi'
  write(123,*) 'S^(-1/2)HS^(-1/2)psibar/psibar    psibar = S^{1/2}psi'
  write(124,*) 'S^(-1/2)HS^(-1/2)psi/psi'  
  if(rank==0) then
    do is=1,ns
      tmp = phit_tot(:,:,:,is,1)
      call h_phi(tmp,tmp2,1)
      tmp4 = real(tmp2/tmp)
      mean = 1d0/nn*sum(tmp4)
      std = sqrt(1d0/nn*sum((tmp4-mean)**2))
      write(120,*) real(tmp2/tmp)
      write(120,*) mean, std, eig(is)
      call flush(120)

      !call sn_phi(tmp2,tmp3,1,-1d0)
      tmp = phit_tot(:,:,:,is,1)
      call sh(tmp,tmp3,1)
      tmp4 = real(tmp3/tmp)
      mean = 1d0/nn*sum(tmp4)
      std = sqrt(1d0/nn*sum((tmp4-mean)**2))
      write(121,*) real(tmp3/tmp)
      write(121,*) mean, std, eig(is)
      call flush(121)

      tmp = phit_tot(:,:,:,is,1)
      call sn_phi(tmp,tmp2,1,0.5d0)
      call h_phi(tmp2,tmp3,1)
      tmp4 = real(tmp3/tmp2)
      mean = 1d0/nn*sum(tmp4)
      std = sqrt(1d0/nn*sum((tmp4-mean)**2))
      write(122,*) real(tmp3/tmp2)
      write(122,*) mean, std, eig(is)
      call flush(122)

      tmp = phit_tot(:,:,:,is,1)
      call sn_phi(tmp,tmp2,1,0.5d0)
      call shs(tmp2,tmp3,1)
      tmp4 = real(tmp3/tmp2)
      mean = 1d0/nn*sum(tmp4)
      std = sqrt(1d0/nn*sum((tmp4-mean)**2))
      write(123,*) real(tmp3/tmp2)
      write(123,*) mean, std, eig(is)
      call flush(123)

      tmp = phit_tot(:,:,:,is,1)
      call shs(tmp,tmp2,1)
      tmp4 = real(tmp2/tmp)
      mean = 1d0/nn*sum(tmp4)
      std = sqrt(1d0/nn*sum((tmp4-mean)**2))
      write(124,*) real(tmp2/tmp)
      write(124,*) mean, std, eig(is)
      call flush(124)
    enddo
  endif

  deallocate(eig)

end subroutine
