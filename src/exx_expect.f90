! Label the psuedo wavefunctions by |psi_i>
!The paw exact exchange operator is given by
  !Check Carsten Rostgaard's master thesis: Exact Exchange in Density functional calculations
  !section 6.7 
! v_xx|psi_i> = sum_j^Nocc [v_xx]ij |psi_j> +
!               - sum_a sum_i1,i2 |p_i1>[X_i1,i2+sum_i3,i4 e_i1,i3,i2,i4*P_i3,i4]<p_i2|psi_i>
!               + sum_j^Nocc sum_a sum_i1,i2 |p_i1>Dijfockhat(i1,i2,j,i)<p_i2|psi_j>
!
! where [v_xx]_ij = int (r-r')^(-1) (n_ij+sum_a nhat_ij)) dr
subroutine exx_expect(psi_i,psi_j,psi_n,nx,ny,nz,nocc,exx)
  use main_mod, only : scale_vh,dv,iscf, &
                       nscf
  use atom_mod
  use atom_mod, only : p=>pawinfo, at=> atominfo
  use libpaw_mod 
  use mpi_lib_ours
  implicit none
  integer :: nn,nx,ny,nz,nocc
  complex*16  :: psi_i(nx,ny,nz), psi_j(nx,ny,nz), psi_n(nx,ny,nz,nocc)
  integer :: i,j,k,l,it,is,js,ia,jlm,il,jl
  integer :: ms,ind,st
  integer :: idij,irhoij,ispden,istate
  integer :: klmn_kl,klmn_ij,klmn_il,klmn_kj,lmn2_size,klmn,klmn1,klmn2
  integer :: ilmn_k,jlmn_l,ilmn_i,jlmn_j
  real*8  :: ro(1)
  real*8  :: exx !expectation value of PAW exact exchange operator
  real*8  :: exx_soft_vv
  real*8  :: exx_soft_vv_1
  real*8  :: exx_soft_vv_2
  real*8  :: exx_vv_cv
  real*8  :: exx_vv
  real*8  :: exx_cv
  real*8  :: exx_hat
  real*8, allocatable  :: mat(:,:),mat2(:,:),tensor4(:,:,:,:)
  real*8, allocatable  :: dij_vv(:),dij_cv(:),tmp(:)
  real*8, allocatable  :: dij_vv2(:,:),dij_cv2(:,:)
  real*8, allocatable  :: dijfock(:,:),dijfock_hat(:,:)
  real*8, allocatable  :: nij(:,:,:),pot(:)
  complex*16, allocatable :: tmpc_i(:),tmpc_j(:),tmp2c(:),nhatij(:,:,:)

  if(rank==0) then
    nn=nx*ny*nz
    allocate(nhatij(nx,ny,nz),nij(nx,ny,nz),tmp(nn),stat=st); if(st/=0) stop 'nhatij'
    allocate(pot(nn),tmpc_i(nn),tmpc_j(nn),tmp2c(nn),stat=st); if(st/=0) stop 'pot'
    exx = 0d0
    exx_soft_vv   = 0d0
    exx_soft_vv_1 = 0d0
    exx_soft_vv_2 = 0d0
    exx_vv_cv     = 0d0
    exx_hat       = 0d0
    exx_cv        = 0d0
    exx_vv        = 0d0
    exx_vv_cv     = 0d0
    iscf = nscf+1
    istate = 1 !can be any number since dijfockhat(is,js,j,istate)

    call allocate_dijfock
    call set_rhoij
    call prep_dijfock
    call prep_dijfockhat_and_exx !calculates sum_j <psi_i|[v_x^F(r)]_ji|psi_j> as well 
    call calc_exx
  endif
  contains
    subroutine prep_dijfock
      implicit none
      real*8 :: r
      real*8, allocatable :: eijkl(:,:)

      idij = 1
      ispden = 1
      do ia=1,natom
        it=atom_map(ia)
        ms=p(it)%mstates
        lmn2_size = pawtab(it)%lmn2_size
        allocate(dij_vv(lmn2_size),dij_cv(lmn2_size),stat=st); if(st/=0) stop 'dij_vv,dij_cv'
        allocate(eijkl(lmn2_size,lmn2_size),stat=st); if(st/=0) stop 'eijkl'
        dij_vv=0d0
        dij_cv=0d0
        eijkl = pawtab(it)%eijkl
        do irhoij=1,pawrhoij(ia)%nrhoijsel
          klmn_kl=pawrhoij(ia)%rhoijselect(irhoij)
          ro(1)=pawrhoij(ia)%rhoijp(irhoij,ispden)*pawtab(it)%dltij(klmn_kl)
          ilmn_k=pawtab(it)%indklmn(7,klmn_kl)
          jlmn_l=pawtab(it)%indklmn(8,klmn_kl)
!* Fock ontribution to the element (k,l) of dijfock
          dij_vv(klmn_kl)=dij_vv(klmn_kl)-ro(1)*eijkl(klmn_kl,klmn_kl)
!* Fock ontribution to the element (i,j) of dijfock with (i,j) < (k,l)
!* We reind that i<j and k<l by construction
          do klmn_ij=1,klmn_kl-1
            ilmn_i=pawtab(it)%indklmn(7,klmn_ij)
            jlmn_j=pawtab(it)%indklmn(8,klmn_ij)
!* In ths case, i < l
            klmn_il=jlmn_l*(jlmn_l-1)/2+ilmn_i
!* If k j, one must consider the index of the symmetric element (j,k) ; otherwise, the index of the element (k,j) is calculated.
            if (ilmn_k>jlmn_j) then
              klmn_kj=ilmn_k*(ilmn_k-1)/2+jlmn_j
            else
              klmn_kj=jlmn_j*(jlmn_j-1)/2+ilmn_k
            end if
!* In ths case, (i,l) >= (k,j)
            dij_vv(klmn_ij)=dij_vv(klmn_ij)-ro(1)*eijkl(klmn_il,klmn_kj)
          end do
!* Fock contribution to the element (i,j) of dijfock with (i,j) > (k,l)
!* We reind that i<j and k<l by construction
          do klmn_ij=klmn_kl+1,lmn2_size
            ilmn_i=pawtab(it)%indklmn(7,klmn_ij)
            jlmn_j=pawtab(it)%indklmn(8,klmn_ij)
!* In ths case, k < j
            klmn_kj=jlmn_j*(jlmn_j-1)/2+ilmn_k
!* If i l, one must consider the index of the symmetric element (l,i) ; otherwise, the index of the element (i,l) is calculated.
            if (ilmn_i>jlmn_l) then
              klmn_il=ilmn_i*(ilmn_i-1)/2+jlmn_l
            else
              klmn_il=jlmn_l*(jlmn_l-1)/2+ilmn_i
            end if
!* In ths case, (k,j) >= (i,l)
            dij_vv(klmn_ij)=dij_vv(klmn_ij)-ro(1)*eijkl(klmn_kj,klmn_il)
          end do
        end do
        !write(6,*) 'ia, sum(dij_vv), maxval(dij_vv), minval(dij_vv)', ia, sum(dij_vv), maxval(dij_vv), minval(dij_vv)
! Add the core-valence contribution
        do klmn_ij=1,lmn2_size
          dij_cv(klmn_ij)=dij_cv(klmn_ij)+pawtab(it)%ex_cvij(klmn_ij)
        end do

        !paw_ij(ia)%dijfock(:,1) = dij_vv(:) + dij_cv(:)
        at(ia)%dijfock = 0d0
        at(ia)%dijfock_vv = 0d0
        at(ia)%dijfock_cv = 0d0
        ind=0d0
        do j=1,ms
            do i=1,j
                ind=ind+1
                !at(ia)%dijfock(i,j)=2d0*dij_vv(ind)+dij_cv(ind)
                !at(ia)%dijfock(j,i)=2d0*dij_vv(ind)+dij_cv(ind)
                at(ia)%dijfock_vv(i,j)=dij_vv(ind) !the factor of 2 is already included
                at(ia)%dijfock_vv(j,i)=dij_vv(ind)
                at(ia)%dijfock_cv(i,j)=dij_cv(ind)
                at(ia)%dijfock_cv(j,i)=dij_cv(ind)
            enddo
        enddo
        at(ia)%dijfock = at(ia)%dijfock_vv + at(ia)%dijfock_cv
        deallocate(dij_vv,dij_cv,eijkl)
      enddo

    end subroutine

    subroutine prep_dijfockhat_and_exx
      implicit none
      integer :: ia,ig,it,is,js,igl,ix,iy,iz
      real*8  :: vv1,vv2

      tmpc_i  = reshape(psi_i,(/nn/))
      tmpc_j  = reshape(psi_j,(/nn/))
      do j=1,nocc
        !write(6,*) 'j, norm(phit_tot)', j, sum(abs(phit_tot(:,:,:,j,1))**2d0)*dv
        !write(6,*) 'i,j dijfockhat', istate,j
        !call flush(6)
        tmp2c = reshape(psi_n(:,:,:,j),(/nn/))
        call get_nhat_elem(psi_j,psi_n(:,:,:,j),nhatij)

        nij(:,:,:) = conjg(psi_j)*psi_n(:,:,:,j)
        tmp = reshape(nij,(/nn/))
        call vh_sub(tmp,pot,scale_vh)
        exx_soft_vv_1 = exx_soft_vv_1 - sum(conjg(tmpc_i)*pot*tmp2c)*dv
        vv1 = - sum(conjg(tmpc_i)*pot*tmp2c)*dv

        tmp = reshape(nhatij,(/nn/))
        call vh_sub(tmp,pot,scale_vh)
        exx_soft_vv_2 = exx_soft_vv_2 - sum(conjg(tmpc_i)*pot*tmp2c)*dv
        vv2 = - sum(conjg(tmpc_i)*pot*tmp2c)*dv
        write(6,*) 'difjockhat prep nhatij', sum(nhatij), sum(psi_j), sum(psi_n(:,:,:,j))

        write(6,*) 'j, vv_1, vv_2', j,vv1,vv2

        nij(:,:,:) = conjg(psi_j)*psi_n(:,:,:,j)+ nhatij
        tmp = reshape(nij,(/nn/))
        call vh_sub(tmp,pot,scale_vh)

        exx = exx - sum(conjg(tmpc_i)*pot*tmp2c)*dv
        exx_soft_vv = exx_soft_vv - sum(conjg(tmpc_i)*pot*tmp2c)*dv


        do ia=1,natom
          !at(ia)%dijfockhat = 0d0
          at(ia)%dijfockhat(:,:,j,istate) = 0d0
          it = atom_map(ia)
          ms = p(it)%mstates
          do is=1,ms
            do js=1,ms
              nhatij = 0d0
              do igl=1,(2*p(it)%nl-1)**2
                do ig=1,ngrid(ia)
                  ix=at(ia)%local_grid(ig,1)
                  iy=at(ia)%local_grid(ig,2)
                  iz=at(ia)%local_grid(ig,3)
                  nhatij(ix,iy,iz)=nhatij(ix,iy,iz)+p(it)%qijlm(igl,is,js)*at(ia)%local_g3d(ig,igl)
                enddo
              enddo
              tmp = reshape(nhatij,(/nn/))
              at(ia)%dijfockhat(is,js,j,istate) = at(ia)%dijfockhat(is,js,j,istate) - sum(pot*tmp)*dv
            enddo
          enddo
        enddo

      enddo

      write(6,*) 'exx_soft_vv: ', exx_soft_vv
      write(6,*) 'exx_soft_vv_1: ', exx_soft_vv_1
      write(6,*) 'exx_soft_vv_2: ', exx_soft_vv_2
    end subroutine

    subroutine calc_exx 
      implicit none
      complex*16, allocatable :: ca_i(:), ca_j(:), ca2(:)

      do ia=1,natom
        it = atom_map(ia)
        ms = p(it)%mstates

        allocate(ca_i(ms),ca_j(ms),ca2(ms),stat=st); if(st/=0) stop 'ca,ca2'
        
        call proj(ia,psi_i,ca_i,ms,1) 
        call proj(ia,psi_j,ca_j,ms,1) 

        do is=1,ms
          do js=1,ms
            exx = exx + conjg(ca_i(is))*at(ia)%dijfock(is,js)*ca_j(js)
            exx_vv_cv = exx_vv_cv + conjg(ca_i(is))*at(ia)%dijfock(is,js)*ca_j(js)
            exx_vv = exx_vv + conjg(ca_i(is))*at(ia)%dijfock_vv(is,js)*ca_j(js)
            exx_cv = exx_cv + conjg(ca_i(is))*at(ia)%dijfock_cv(is,js)*ca_j(js)
          enddo
        enddo

        do j=1,nocc
          call proj(ia,psi_n(:,:,:,j),ca2,ms,1) 
          do is=1,ms
            do js=1,ms
              exx = exx + conjg(ca_i(is))*at(ia)%dijfockhat(is,js,j,istate)*ca2(js)
              exx_hat = exx_hat + conjg(ca_i(is))*at(ia)%dijfockhat(is,js,j,istate)*ca2(js)
            enddo
          enddo
        enddo

        deallocate(ca_i,ca_j,ca2)
      enddo

      write(6,*) 'exx_vv_cv: ', exx_vv_cv
      write(6,*) 'exx_vv: ', exx_vv
      write(6,*) 'exx_cv: ', exx_cv
      write(6,*) 'exx_hat: ', exx_hat
      write(6,*) 'exx: ', exx
    end subroutine
end subroutine
