subroutine prep_tbar_sh
  use main_mod
  use tddft_mod, only : tbar_taylor_power, dt, phi_bar_tot
  use mat_module
  use mpi_lib_ours
  use paw_mod
  use atom_mod
  use atom_mod, only : p=> pawinfo, at=> atominfo
  implicit none
  integer :: it, ia, i,j, is, js, k
  integer :: ix, iy, iz, ik, ms, ns, nr
  integer :: ixs,iys,izs,ixe,iye,ize,jx,jy,jz
  integer :: ifs(3), ife(3), nxyzf(3)
  complex*16, allocatable :: tmp_wf3d(:,:,:,:), tmp_wf3d_fine(:,:,:,:)
  complex*16, allocatable :: ttmp_wf3d(:,:,:,:)
  complex*16, allocatable :: mat_chk(:,:), mat_chk2(:,:)
  complex*16 :: cin(nx,ny,nz),cout(nx,ny,nz), tmp(nx,ny,nz)
  complex*16 :: ttmp(nx,ny,nz)
  complex*16 :: tmp_wf(nn)
  complex*16 :: ci = (0d0,1d0)
  real*8  :: x,y,z,xc,yc,zc,r
  real*8  :: xx,yy,zz
  real*8  :: nxyzf_r(3)
  complex*16  :: a, b, c, aa, bb, cc, dd
  logical :: edge

  !below: 1st and 2nd deriv, tmp ptilde, t applied on ptilde
  real*8, allocatable :: dptilde(:), ddptilde(:), tmpp(:)
  complex*16, allocatable :: tmp3(:,:,:,:),tmp3t(:,:,:,:)
  integer :: counter
  integer, allocatable :: index_state(:)
  real*8,  allocatable :: mat(:,:),vec(:,:),eig(:)
  real*8,  allocatable :: mat1(:,:),mat2(:,:),mat3(:,:)
  real*8,  allocatable :: v(:,:), ymat(:,:), l(:,:), bmat(:,:)

  ik=1
  do ia=1,natom
    it=atom_map(ia)
    ms = p(it)%mstates
    
    call alloc_c3d

    if(rank==0) then
      allocate(tmp_wf3d(nx,ny,nz,ms),stat=stat)
      if(stat/=0) stop 'tmp_wf3d alloc problem in prep_tbar'
      allocate(ttmp_wf3d(nx,ny,nz,ms),stat=stat)
      if(stat/=0) stop 'ttmp_wf3d alloc problem in prep_tbar'

      allocate(mat_chk(ms,ms),stat=stat)
      if(stat/=0) stop 'mat_chk alloc problem in prep_tbar'
      allocate(mat_chk2(ms,ms),stat=stat)
      if(stat/=0) stop 'mat_chk2 alloc problem in prep_tbar'

      ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
      izs=at(ia)%ir_start(3)
      ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
      ize=izs+p(it)%nrough(3)-1

      tmp_wf3d=0d0
      if(.not.at(ia)%edge) then
        tmp_wf3d(ixs:ixe,iys:iye,izs:ize,:)=at(ia)%local_p3d1_c(:,:,:,:,ik)
        !b = sum(conjg(at(ia)%local_p3d1_c(:,:,:,1,ik))*ttmp(ixs:ixe,iys:iye,izs:ize))*dv
      else
        do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
          jx=mod(ix+nx-1,nx)+1
          jy=mod(iy+ny-1,ny)+1
          jz=mod(iz+nz-1,nz)+1
          tmp_wf3d(jx,jy,jz,:)=at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,ik)
        enddo;enddo;enddo
      endif

      do i=1,ms
        cin=tmp_wf3d(:,:,:,i)
        call fft3d_forward(nx,ny,nz,cin,cout)
        cout=cout*ek(:,:,:,ik)
        call fft3d_backward(nx,ny,nz,cout,ttmp_wf3d(:,:,:,i))
      enddo
      at(ia)%local_c_c(:,:,:,:) = ttmp_wf3d

      do i=1,ms
        do j=1,ms
          at(ia)%wmat(i,j) = sum(conjg(tmp_wf3d(:,:,:,i))*at(ia)%local_c_c(:,:,:,j))*dv
        enddo
      enddo

      call mat_inv(at(ia)%wmat, at(ia)%amat)

      at(ia)%local_c1_c = 0d0
      do i=1,ms; do j=1,ms
        at(ia)%local_c1_c(:,:,:,i) = at(ia)%local_c1_c(:,:,:,i) + at(ia)%local_c_c(:,:,:,j)*at(ia)%amat(j,i)
      enddo;  enddo

      write(6,*) '  '
      write(6,*) 'overlap matrix of p3d1 and c1 for ia: ', ia
      do i=1,ms
        do j=1,ms
          mat_chk(i,j) = sum(conjg(tmp_wf3d(:,:,:,i))*&
                            at(ia)%local_c1_c(:,:,:,j))*dv
        enddo
        write(6,*) abs(mat_chk(i,:))
      enddo
      write(6,*) '  '

      at(ia)%local_c=0d0
      write(6,*) 'checking if c = w c1 matches c_c for ia: ', ia
      do i=1,ms
        do j=1,ms
          at(ia)%local_c(:,:,:,i) = at(ia)%local_c(:,:,:,i) + at(ia)%wmat(j,i)*at(ia)%local_c1_c(:,:,:,j)
        enddo
        mat_chk(i,:) = at(ia)%sinv(i,1)*conjg(at(ia)%wmat(:,i))
      enddo
      if (sum(abs(at(ia)%local_c - at(ia)%local_c_c))*dv > 1e-11) then
        write(6,*) 'diff local_c = wc1 and local_c_c greater than 1e-12:', &
              sum(abs(at(ia)%local_c - at(ia)%local_c_c))*dv
        stop
      else
        write(6,*) 'check passed'
      endif
      write(6,*) '  '

      at(ia)%c_taylor = 0d0
      mat_chk2 = -ci*dt*0.5d0*mat_chk
      at(ia)%cbarmat = mat_chk2
      at(ia)%cmat = mat_chk
      do i=1,tbar_taylor_power
        at(ia)%c_taylor = at(ia)%c_taylor + mat_chk2
        at(ia)%c_taylor_terms(:,:,i) = mat_chk2
        if(i.ne.tbar_taylor_power) mat_chk2 = 1/dble(i+1)*matmul(at(ia)%cbarmat,mat_chk2)
      enddo
      deallocate(tmp_wf3d, ttmp_wf3d, mat_chk, mat_chk2)
    endif
    call bcast_c16(at(ia)%local_c_c,size(at(ia)%local_c_c),0)
    call bcast_c16(at(ia)%local_c1_c,size(at(ia)%local_c1_c),0)
    call bcast_c16(at(ia)%c_taylor_terms,size(at(ia)%c_taylor_terms),0)
    call bcast_c16(at(ia)%c_taylor,size(at(ia)%c_taylor),0)
    call bcast_c16(at(ia)%cbarmat,size(at(ia)%cbarmat),0)
    call bcast_c16(at(ia)%cmat,size(at(ia)%cmat),0)
    call bcast_c16(at(ia)%wmat,size(at(ia)%wmat),0)
    call bcast_c16(at(ia)%amat,size(at(ia)%amat),0)

  enddo

  !if(rank==0)call debug_tbar_sh
  call debug_tbar_sh_multi
contains
  subroutine alloc_c3d
    implicit none

    allocate(at(ia)%local_c(nx,ny,nz,ms))
    if(stat/=0) stop 'local_c alloc problem in prep_tbar'

    allocate(at(ia)%local_c_c(nx,ny,nz,ms))
    if(stat/=0) stop 'local_c_c alloc problem in prep_tbar'

    allocate(at(ia)%local_c1_c(nx,ny,nz,ms))
    if(stat/=0) stop 'local_c1_c alloc problem in prep_tbar'

    allocate(at(ia)%wmat(ms,ms), at(ia)%amat(ms,ms))  
    if(stat/=0) stop 'wmat, amat'

    allocate(at(ia)%c_taylor(ms, ms),stat=stat)
    if(stat/=0) stop 'c_taylor alloc problem in prep_tbar'

    allocate(at(ia)%c_taylor_terms(ms, ms, tbar_taylor_power),stat=stat)
    if(stat/=0) stop 'c_taylor_terms alloc problem in prep_tbar'

    allocate(at(ia)%cmat(ms, ms),stat=stat)
    if(stat/=0) stop 'cmat alloc problem in prep_tbar'

    allocate(at(ia)%cbarmat(ms, ms),stat=stat)
    if(stat/=0) stop 'cbarmat alloc problem in prep_tbar'
  end subroutine alloc_c3d

  subroutine debug_tbar_sh
    use tddft_mod, only : phi_bar_tot
    implicit none
    real*8, allocatable :: overlap(:,:)
    complex*16 :: tmp(nx,ny,nz), tmp1(nx,ny,nz), tmp2(nx,ny,nz)
    complex*16 :: tmp3(nx,ny,nz), tmp4(nx,ny,nz), tmp5(nx,ny,nz)
    complex*16 :: tmp6(nx,ny,nz), tmp7(nx,ny,nz), tmp8(nx,ny,nz)
    complex*16 :: tmp9(nx,ny,nz), tmp10(nx,ny,nz), tmp11(nx,ny,nz)
    complex*16 :: tmp12(nx,ny,nz), tmp13(nx,ny,nz), tmp14(nx,ny,nz)
    complex*16 :: tmp15(nx,ny,nz), tmp16(nx,ny,nz), tmp17(nx,ny,nz)
    complex*16 :: tmp18(nx,ny,nz), tmp19(nx,ny,nz), tmp20(nx,ny,nz)
    complex*16 :: tmp21(nx,ny,nz), tmp22(nx,ny,nz), tmp23(nx,ny,nz)
    complex*16 :: wf_term(nx,ny,nz,tbar_taylor_power)
    complex*16 :: wf_term1(nx,ny,nz,tbar_taylor_power)
    complex*16, allocatable :: ca(:), sn(:), y(:)
    complex*16, allocatable :: ca2(:), ca3(:), mat_tmp(:,:)
    complex*16, allocatable :: ca4(:), ca5(:)
    complex*16 :: wf1(nn), wf2(nn), wf3(nn), wf4(nn), wf5(nn)
    complex*16 :: wf6(nn), wf7(nn), wf8(nn), wf9(nn), wf10(nn)
    complex*16 :: aa, bb, cc, dd, ee, ff, gg, hh
    complex*16 :: aaa(tbar_taylor_power)
    
    real*8  :: nrm0, nrm1, nrm2, nrm3, nrm4
    integer :: js, ks, iis, jjs, kks, ja, jt

    write(6,*) 'debug_tbar_sh'
    write(6,*) 'now checking s^-1t vs tbar=t+sum_i |p_i>y_i<c_i| various methods'

    ! wf1 = s^-1t psi_0
    do js=1,nocc
      !tmp = phi_bar_tot(:,:,:,2,ik)
      tmp = phi_bar_tot(:,:,:,js,ik)
      call fft3d_forward(nx,ny,nz,tmp,tmp1)
      tmp1 = tmp1*ek(:,:,:,ik)
      call fft3d_backward(nx,ny,nz,tmp1,tmp2)
      call sn_phi(tmp2, tmp1, 1,-1d0)
      !tmp1=tmp1-tmp2
      wf1 = reshape(tmp1, (/nn/))

      !tmp = phi_bar_tot(:,:,:,1,ik)
      !call fft3d_forward(nx,ny,nz,tmp,tmp1)
      !tmp1 = tmp1*ek(:,:,:,ik)
      !call fft3d_backward(nx,ny,nz,tmp1,tmp2)
      !tmp1 = tmp
      tmp3 = 0d0
      tmp4 = 0d0
      call cmat_phi(tmp,tmp5,ik)
      call c_phi(tmp,tmp3,ik)
      !wf5 = reshape(tmp2+tmp5, (/nn/))
      !wf2 = reshape(tmp3, (/nn/))
      !wf5 = reshape(tmp5, (/nn/))
      tmp5 = 0d0
      do ia=1,natom
        it=atom_map(ia)
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        ms = p(it)%mstates

        !allocate(ca(ms))
        allocate(ca(ms),ca2(ms), ca3(ms))
        allocate(y(ms),stat=stat)
        if(stat/=0) stop 'y alloc problem in prep_tbar'


        !call tbar_proj2(ia, tmp, ca, ms)
        !call tbar_proj3(ia, tmp, ca3, ms)
        do is=1,ms
          ca(is) = sum(conjg(at(ia)%local_c_c(:,:,:,is))*tmp)*dv
          ca3(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv
        enddo
        call proj1(ia, tmp2, ca2, ms, 1)
        y=at(ia)%sinv(:,ik)

        !write(6,*) 'ia, ca=<zeta T|psi>', ia, ca(1:2)
        !write(6,*) 'is, ca2=<zeta|T psi>', ia, ca2(1:2)
        !write(6,*) 'is, maxval(abs(ca2-ca1))', ia, maxval(abs(ca-ca2))
        !write(6,*) 'is, ca2-ca1', ia, ca-ca2
  !      write(6,*) 'test1', y
        if(.not.at(ia)%edge) then
          do is=1,ms
           ! tmp3(ixs:ixe,iys:iye,izs:ize)=tmp3(ixs:ixe,iys:iye,izs:ize)+&
           ! at(ia)%local_p3d1_c(:,:,:,is,1)*y(is)*ca(is)
            tmp4(ixs:ixe,iys:iye,izs:ize)=tmp4(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*y(is)*sum(conjg(at(ia)%wmat(:,is))*ca3(:))
            tmp5(ixs:ixe,iys:iye,izs:ize)=tmp5(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cmat(is,:)*ca3(:))
          enddo
        else
          do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
            jx=mod(ix+nx-1,nx)+1
            jy=mod(iy+ny-1,ny)+1
            jz=mod(iz+nz-1,nz)+1
            !tmp3(jx,jy,jz)=tmp3(jx,jy,jz)+sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,1)*y(:)*ca(:))
            do ks=1,ms
              tmp4(jx,jy,jz)=tmp4(jx,jy,jz)+sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,1)*y(:)*&
                conjg(at(ia)%wmat(ks,:)))*ca3(ks)
              tmp5(jx,jy,jz)=tmp5(jx,jy,jz)+sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,1)*&
                at(ia)%cmat(:,ks))*ca3(ks)
            enddo
          enddo;enddo;enddo
        endif

        !deallocate(ca, y)
        deallocate(ca, ca2, ca3, y)
      enddo

      wf2 = reshape(tmp2+tmp3, (/nn/))
      wf3 = reshape(tmp2+tmp4, (/nn/))
      wf4 = reshape(tmp2+tmp5, (/nn/))
      wf5 = reshape(tmp2+tmp5, (/nn/))
      !wf3 = reshape(tmp4, (/nn/))
      !wf4 = reshape(tmp5, (/nn/))
      aa = sum(abs(wf1-wf2))*dv
      bb = sum(abs(wf1-wf3))*dv
      cc = sum(abs(wf1-wf4))*dv
      dd = sum(abs(wf1-wf5))*dv
      nrm0 = sum(conjg(wf1)*wf1)*dv
      nrm1 = sum(conjg(wf2)*wf2)*dv
      nrm2 = sum(conjg(wf3)*wf3)*dv
      nrm3 = sum(conjg(wf4)*wf4)*dv
      nrm4 = sum(conjg(wf5)*wf5)*dv
      write(6,*) 'js, aa, bb, cc, dd', js, real(aa), real(bb), real(cc), real(dd)
      write(6,*) 'nrm0,1,2,3,4 ', nrm0, nrm1, nrm2, nrm3, nrm4
      !close(310+js)
      !open(unit=310, file='t_bar_2_debug')
      !do i=1,nn
      !  write(310+js,*) i, real(wf1(i)), real(wf2(i)), real(wf1(i)-wf2(i))!, real(wf3(i))
      !  !write(310,*) i, real(wf1(i)-wf2(i))
      !enddo
      !close(310+js)
    enddo
    !close(310)
    !open(unit=310, file='t_bar_2_debug')
    !do i=1,nn
    !  write(310,*) i, real(wf1(i)), real(wf2(i)), real(wf1(i)-wf2(i))!, real(wf3(i))
    !  !write(310,*) i, real(wf1(i)-wf2(i))
    !enddo
    !close(310)
  !------------------------------------
    !write(6,*) 'check orthogonality on grid for different atoms'
    !do ia=1,natom
    !  do ja=1,natom
    !    write(6,*) 'ia, ja', ia, ja
    !    jt=atom_map(ja)
    !    ixs=at(ja)%ir_start(1);iys=at(ja)%ir_start(2)
    !    izs=at(ja)%ir_start(3)
    !    ixe=ixs+p(jt)%nrough(1)-1;iye=iys+p(jt)%nrough(2)-1
    !    ize=izs+p(jt)%nrough(3)-1
    !    it=atom_map(ia)
    !    !ms = p(it)%mstates

    !    allocate(mat_tmp(p(it)%mstates,p(jt)%mstates),stat=stat)
    !    if(stat/=0) stop 'mat_tmp alloc problem in prep_tbar'

    !    do is=1,p(it)%mstates
    !      do js=1,p(jt)%mstates
    !        mat_tmp(is,js) = sum(conjg(at(ia)%local_c1_c(ixs:ixe,iys:iye,izs:ize,is))*&
    !                at(ja)%local_p3d1_c(:,:,:,js,1))*dv
    !      enddo
    !      write(6,*) real(mat_tmp(is,:))
    !    enddo
    !    deallocate(mat_tmp)
    !    write(6,*) '  '
    !  enddo
    !enddo

  !------------------------------------
    write(6,*) 'now checking e^sum vs e^sum approx 1 + sum'

    do js=1,nocc
      wf_term = 0d0
      wf_term1 = 0d0
      tmp = phi_bar_tot(:,:,:,js,ik)
      nrm0 = sum(conjg(tmp)*tmp)*dv

      call add_ctaylor_phi(tmp,tmp5,ik)
      call add_ctaylor_phi1(tmp,tmp16,ik)
      !call c_taylor_loop(tmp,tmp21,ik)
      !tmp1 = tmp
      tmp1 = 0d0
      tmp17 = 0d0
      tmp18 = 0d0
      tmp19 = 0d0
      do ia=1,natom
        it=atom_map(ia)
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        ms = p(it)%mstates

        allocate(ca(ms),stat=stat)
        if(stat/=0) stop 'ca alloc problem in prep_tbar'
        allocate(ca2(ms),stat=stat)
        if(stat/=0) stop 'ca2 alloc problem in prep_tbar'

        do is=1,ms; ca(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv; enddo
        !write(6,*) 'ia, ca', ia, ca
        !write(6,*) 'ia', ia
        !write(6,*) 'ixs, ixe', ixs, ixe
        !write(6,*) 'iys, iye', iys, iye
        !write(6,*) 'izs, ize', izs, ize

        if(.not.at(ia)%edge) then
          tmp20 = 0d0
          do iis=1,ms
            tmp20(ixs:ixe,iys:iye,izs:ize) = tmp20(ixs:ixe,iys:iye,izs:ize)+&
               at(ia)%local_p3d1_c(:,:,:,iis,ik)*sum(at(ia)%cbarmat(iis,:)*ca(:))
          enddo
          do is=1,ms; ca2(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp20)*dv; enddo
          !write(6,*) 'ia, ca2',ia, ca2
          do is=1,ms
            tmp1(ixs:ixe,iys:iye,izs:ize)=tmp1(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%c_taylor_terms(is,:,1)*ca(:))
            tmp19(ixs:ixe,iys:iye,izs:ize)=tmp19(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,ik)*sum(at(ia)%cbarmat(is,:)*ca2(:))/2d0
            tmp17(ixs:ixe,iys:iye,izs:ize)=tmp17(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cbarmat(is,:)*ca(:))
          enddo
        else
          do is=1,ms
            do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
              jx=mod(ix+nx-1,nx)+1
              jy=mod(iy+ny-1,ny)+1
              jz=mod(iz+nz-1,nz)+1
              tmp1(jx,jy,jz)=tmp1(jx,jy,jz)+at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)&
                  *sum(at(ia)%c_taylor(is,:)*ca(:))
            enddo;enddo;enddo
          enddo
        endif

        deallocate(ca)
        deallocate(ca2)
      enddo
      do ia=1,natom
        it=atom_map(ia)
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        ms = p(it)%mstates

        allocate(ca(ms),stat=stat)
        if(stat/=0) stop 'ca, y alloc problem in prep_tbar'
        do is=1,ms; ca(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp17)*dv; enddo
        !write(6,*) 'ia, ca_18', ia, ca
        if(.not.at(ia)%edge) then
          do is=1,ms
            tmp18(ixs:ixe,iys:iye,izs:ize)=tmp18(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cbarmat(is,:)*ca(:))
          enddo
        endif
        deallocate(ca)
      enddo
      write(6,*) 'i, diff 1 17', sum(abs(tmp17-tmp1))*dv
      write(6,*) 'i, diff 1 19', sum(abs(tmp19-tmp1))*dv
      write(6,*) 'i, diff 1 18', sum(abs(tmp18-tmp1))*dv
      write(6,*) 'i, diff 18 19', sum(abs(tmp19-tmp18))*dv
      wf_term(:,:,:,1) = tmp17
      wf_term(:,:,:,1) = tmp18

      tmp3 = 0d0
      tmp4 = tmp
      tmp6 = tmp
      tmp7 = tmp
      tmp8 = 0d0
      tmp9 = tmp
      tmp10 = 0d0
      tmp11 = tmp
      tmp13 = tmp
      tmp14 = tmp
      tmp15 = tmp

      do i=1,tbar_taylor_power
        tmp8=0d0
        cin=tmp6
        call fft3d_forward(nx,ny,nz,cin,cout)
        cout=cout*ek(:,:,:,ik)
        call fft3d_backward(nx,ny,nz,cout,tmp6)
        call sn_phi(tmp6,tmp8,ik,-1d0)
        tmp6=tmp8-tmp6
        !cin=tmp6
        !call sum_inv_T_phi(cin,tmp6,ik)
        tmp6=-ci*dt*0.5d0/dble(i)*tmp6
        !wf_term(:,:,:,i) = tmp6
        !wf_term1(:,:,:,i) = tmp6
        !write(6,*) 'i tmp8', sum(conjg(tmp8)*tmp8)*dv
        !write(6,*) 'i tmp6', i, sum(conjg(tmp6)*tmp6)*dv
        tmp7=tmp7+tmp6

        if(i .eq. 1) then 
          cin = tmp
          call fft3d_forward(nx,ny,nz,cin,cout)
          cout=cout*ek(:,:,:,ik)
          call fft3d_backward(nx,ny,nz,cout,tmp10)
        else
          cin = tmp10
          call fft3d_forward(nx,ny,nz,cin,cout)
          cout=cout*ek(:,:,:,ik)
          call fft3d_backward(nx,ny,nz,cout,tmp10)
        endif
        !wf_term(:,:,:,i) = tmp10
        !write(6,*) 'i, wf_diff_term', i, sum(abs(wf_term(:,:,:,i)-wf_term1(:,:,:,i)))*dv
        tmp12=0d0
        tmp8=0d0
        tmp3=0d0
        !tmp16=0d0
        do ia=1,natom
          it=atom_map(ia)
          ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
          izs=at(ia)%ir_start(3)
          ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
          ize=izs+p(it)%nrough(3)-1
          ms = p(it)%mstates

          allocate(ca(ms),ca2(ms), ca3(ms), ca4(ms))
          if(stat/=0) stop 'ca alloc problem in prep_tbar'
          allocate(y(ms),stat=stat)
          if(stat/=0) stop 'y alloc problem in prep_tbar'
          y=at(ia)%sinv(:,ik)

          !write(6,*) 'ia', ia, ca
          !write(6,*) 'ia2', ia, ca2
    !      write(6,*) 'test1', y
          if(.not.at(ia)%edge) then
            !write(6,*) 'i', i
            do is=1,ms
              !ca(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv
              ca2(is) = sum(conjg(at(ia)%local_c_c(:,:,:,is))*tmp13)*dv
              ca3(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp14)*dv
              ca4(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv
            enddo
            call proj1(ia, tmp10, ca, ms, ik)
            !write(6,*) 'ia,i ca3', ia, i, ca3
            do is=1,ms
              tmp3(ixs:ixe,iys:iye,izs:ize)=tmp3(ixs:ixe,iys:iye,izs:ize)+&
                at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cbarmat(is,:)*ca3(:))/dble(i)
              tmp8(ixs:ixe,iys:iye,izs:ize)=tmp8(ixs:ixe,iys:iye,izs:ize)+&
                (-ci*dt*0.5d0)/dble(i)*at(ia)%local_p3d1_c(:,:,:,is,1)*y(is)*ca2(is)
              tmp12(ixs:ixe,iys:iye,izs:ize)=tmp12(ixs:ixe,iys:iye,izs:ize)+&
                (-ci*dt*0.5d0)/dble(i)*at(ia)%local_p3d1_c(:,:,:,is,1)*y(is)*ca(is)
              !tmp16(ixs:ixe,iys:iye,izs:ize)=tmp16(ixs:ixe,iys:iye,izs:ize)+&
              !  at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%c_taylor_terms(is,:,i)*ca4(:))
            enddo
            !tmp13(ixs:ixe,iys:iye,izs:ize)=tmp12(ixs:ixe,iys:iye,izs:ize)
            !write(6,*) 'ia, i, tmp3', ia, i, sum(conjg(tmp3)*tmp3)*dv
            !write(6,*) 'ia, i, ca-c2', ia, i, maxval(real(ca-ca2)), maxval(aimag(ca-ca2))
            !write(6,*) 'ia, i, ca', ia, i, ca(1:2)
            !aaa(i) = aaa(i) + sum(conjg(tmp3)*tmp3)*dv
          else
            if(i .eq. 1) then 
              do is=1,ms
                ca3(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv
              enddo
              !call proj1(ia,tmp,ca,ms,ik)
            else
              do is=1,ms
                ca3(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp3)*dv/dble(i)
              enddo

              !call proj1(ia,-ci*dt/2d0/dble(i)*tmp6,ca,ms,ik)
            endif
            do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
              do is=1,ms
                tmp3(jx,jy,jz)=tmp3(jx,jy,jz)+at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)&
                    *sum(at(ia)%cbarmat(is,:)*ca3(:))
              enddo
            enddo;enddo;enddo
          endif
          deallocate(ca, ca2, ca3, ca4, y)
        enddo
        !write(6,*) 'i, tmp8-tmp12 ', i, maxval(real(tmp8-tmp12)), maxval(aimag(tmp8-tmp12))
        tmp14=tmp3
        tmp4=tmp4+tmp14
        tmp13=tmp8
        tmp9=tmp9+tmp13
        tmp10=tmp12
        tmp11=tmp11+tmp12
        !tmp15=tmp15+tmp16
        wf_term1(:,:,:,i) = tmp3
      enddo

      !do i=1,tbar_taylor_power
      !  write(6,*) 'i, test c^2',i,  sum(abs(wf_term(:,:,:,i)-wf_term1(:,:,:,i)))*dv
      !enddo
      !write(6,*) 'i, test c^2',i,  sum(abs(wf_term(:,:,:,1)-wf_term1(:,:,:,1)))*dv
      wf1 = reshape(tmp1, (/nn/)) ! same as add_C_taylor
      wf2 = reshape(tmp, (/nn/))  !intput
      wf3 = reshape(tmp4, (/nn/)) !e^sum = 1 + sum + ...
      wf4 = reshape(tmp5, (/nn/)) ! add_C_taylor
      wf5 = reshape(tmp7, (/nn/)) !e^-idt/2S^-1T=1-idt/2S^-2T...
      wf6 = reshape(tmp9, (/nn/)) !
      wf7 = reshape(tmp11, (/nn/)) !
      wf8 = reshape(tmp16, (/nn/)) !
      wf9 = reshape(tmp21, (/nn/)) !
      
      !nrm1 = sum(conjg(wf1)*wf1)*dv
      !nrm2 = sum(conjg(wf4)*wf4)*dv
      !nrm3 = sum(conjg(wf3)*wf3)*dv
      !nrm4 = sum(conjg(wf5)*wf5)*dv
      aa = sum(abs(wf5-wf1))*dv !exact taylor - add_C_taylor
      bb = sum(abs(wf5-wf4))*dv
      cc = sum(abs(wf5-wf3))*dv
      dd = sum(abs(wf5-wf6))*dv
      ee = sum(abs(wf5-wf7))*dv
      ff = sum(abs(wf5-wf8))*dv
      !gg = sum(abs(wf5-wf9))*dv
      !write(6,*) 'js, wf1', js, maxval(real(wf1)), maxval(aimag(wf1))
      write(6,*) 'js, aa, bb, cc, :', js, real(aa), real(bb), real(cc)
      write(6,*) 'js, dd, ee, ff:', js, real(dd), real(ee), real(ff)
      !write(6,*) 'js, gg:', js, real(gg)
      !write(6,*) 'js, nrm0,1,2,3,4:', js, nrm1, nrm2, nrm3, nrm4
      !write(6,*) 'js, wf3', js, maxval(real(wf3)), maxval(aimag(wf3))
      !write(6,*) 'js, diff e', js, aa
    enddo

    write(6,*) 'now checking e^s^-1t vs split step'

    do js=1,nocc
      tmp = phi_bar_tot(:,:,:,js,ik)
      tmp1 = tmp
      tmp2 = tmp
      do i=1,tbar_taylor_power
        call st_phi(-ci*dt/dble(i)*tmp1, tmp3, ik)
        !call t_phi(-ci*dt/dble(i)*tmp1, tmp3, ik)
        !call sum_inv_T_phi(-ci*dt/dble(i)*tmp1, tmp3, ik)
        !aa = aa/dble(i)
        tmp2 = tmp2 +tmp3
        tmp1 = tmp3
      enddo
      call prop_tbar_split_sh(tmp, tmp4, ik)  
     
      wf1 = reshape(tmp2, (/nn/))
      wf3 = reshape(tmp4, (/nn/))
      nrm0 = sum(conjg(wf1)*wf1)*dv
      nrm1 = sum(conjg(wf3)*wf3)*dv
      aa = sum(abs(wf1-wf3))*dv
      write(6,*) 'js, aa, nrm0, nrm1:', js, real(aa), nrm0, nrm1
    enddo
  end subroutine debug_tbar_sh

  subroutine debug_tbar_sh_multi
    use tddft_mod, only : phi_bar_tot, phi_bar
    implicit none
    real*8, allocatable :: overlap(:,:)
    complex*16 :: tmp(nx,ny,nz), tmp1(nx,ny,nz), tmp2(nx,ny,nz)
    complex*16 :: tmp3(nx,ny,nz), tmp4(nx,ny,nz), tmp5(nx,ny,nz)
    complex*16 :: tmp6(nx,ny,nz), tmp7(nx,ny,nz), tmp8(nx,ny,nz)
    complex*16 :: tmp9(nx,ny,nz), tmp10(nx,ny,nz), tmp11(nx,ny,nz)
    complex*16 :: tmp12(nx,ny,nz), tmp13(nx,ny,nz), tmp14(nx,ny,nz)
    complex*16 :: tmp15(nx,ny,nz), tmp16(nx,ny,nz), tmp17(nx,ny,nz)
    complex*16 :: tmp18(nx,ny,nz), tmp19(nx,ny,nz), tmp20(nx,ny,nz)
    complex*16 :: tmp21(nx,ny,nz), tmp22(nx,ny,nz), tmp23(nx,ny,nz)
    complex*16 :: wf_term(nx,ny,nz,tbar_taylor_power)
    complex*16 :: wf_term1(nx,ny,nz,tbar_taylor_power)
    complex*16, allocatable :: ca(:), sn(:), y(:)
    complex*16, allocatable :: ca2(:), ca3(:), mat_tmp(:,:)
    complex*16, allocatable :: ca4(:), ca5(:)
    complex*16 :: wf1(nn), wf2(nn), wf3(nn), wf4(nn), wf5(nn)
    complex*16 :: wf6(nn), wf7(nn), wf8(nn), wf9(nn), wf10(nn)
    complex*16 :: aa, bb, cc, dd, ee, ff, gg, hh
    complex*16 :: aaa(tbar_taylor_power)
    real*8  :: nrm0, nrm1, nrm2, nrm3, nrm4
    real*8  :: w1time, w2time
    integer :: js, ks, iis, jjs, kks, ja, jt

    write(6,*) 'debug_tbar_sh_mult'
    write(6,*) 'now checking s^-1t vs tbar=t+sum_i |p_i>y_i<c_i| various methods'

    call scatter_c16(phi_bar_tot,phi_bar,size(phi_bar),0)
    ! wf1 = s^-1t psi_0
    do js=1,nb
      !tmp = phi_bar_tot(:,:,:,2,ik)
      tmp = phi_bar(:,:,:,js,ik)
      call fft3d_forward(nx,ny,nz,tmp,tmp1)
      tmp1 = tmp1*ek(:,:,:,ik)
      call fft3d_backward(nx,ny,nz,tmp1,tmp2)
      call sn_phi(tmp2, tmp1, 1,-1d0)
      !tmp1=tmp1-tmp2
      wf1 = reshape(tmp1, (/nn/))

      !tmp = phi_bar_tot(:,:,:,1,ik)
      !call fft3d_forward(nx,ny,nz,tmp,tmp1)
      !tmp1 = tmp1*ek(:,:,:,ik)
      !call fft3d_backward(nx,ny,nz,tmp1,tmp2)
      !tmp1 = tmp
      tmp3 = 0d0
      tmp4 = 0d0
      call cmat_phi(tmp,tmp5,ik)
      call c_phi(tmp,tmp3,ik)
      !wf5 = reshape(tmp2+tmp5, (/nn/))
      !wf2 = reshape(tmp3, (/nn/))
      !wf5 = reshape(tmp5, (/nn/))
      tmp5 = 0d0
      do ia=1,natom
        it=atom_map(ia)
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        ms = p(it)%mstates

        !allocate(ca(ms))
        allocate(ca(ms),ca2(ms), ca3(ms))
        allocate(y(ms),stat=stat)
        if(stat/=0) stop 'y alloc problem in prep_tbar'


        !call tbar_proj2(ia, tmp, ca, ms)
        !call tbar_proj3(ia, tmp, ca3, ms)
        do is=1,ms
          ca(is) = sum(conjg(at(ia)%local_c_c(:,:,:,is))*tmp)*dv
          ca3(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv
        enddo
        call proj1(ia, tmp2, ca2, ms, 1)
        y=at(ia)%sinv(:,ik)

        !write(6,*) 'ia, ca=<zeta T|psi>', ia, ca(1:2)
        !write(6,*) 'is, ca2=<zeta|T psi>', ia, ca2(1:2)
        !write(6,*) 'is, maxval(abs(ca2-ca1))', ia, maxval(abs(ca-ca2))
        !write(6,*) 'is, ca2-ca1', ia, ca-ca2
  !      write(6,*) 'test1', y
        if(.not.at(ia)%edge) then
          do is=1,ms
           ! tmp3(ixs:ixe,iys:iye,izs:ize)=tmp3(ixs:ixe,iys:iye,izs:ize)+&
           ! at(ia)%local_p3d1_c(:,:,:,is,1)*y(is)*ca(is)
            tmp4(ixs:ixe,iys:iye,izs:ize)=tmp4(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*y(is)*sum(conjg(at(ia)%wmat(:,is))*ca3(:))
            tmp5(ixs:ixe,iys:iye,izs:ize)=tmp5(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cmat(is,:)*ca3(:))
          enddo
        else
          do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
            jx=mod(ix+nx-1,nx)+1
            jy=mod(iy+ny-1,ny)+1
            jz=mod(iz+nz-1,nz)+1
            !tmp3(jx,jy,jz)=tmp3(jx,jy,jz)+sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,1)*y(:)*ca(:))
            do ks=1,ms
              tmp4(jx,jy,jz)=tmp4(jx,jy,jz)+sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,1)*y(:)*&
                conjg(at(ia)%wmat(ks,:)))*ca3(ks)
              tmp5(jx,jy,jz)=tmp5(jx,jy,jz)+sum(at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,:,1)*&
                at(ia)%cmat(:,ks))*ca3(ks)
            enddo
          enddo;enddo;enddo
        endif

        !deallocate(ca, y)
        deallocate(ca, ca2, ca3, y)
      enddo

      wf2 = reshape(tmp2+tmp3, (/nn/))
      wf3 = reshape(tmp2+tmp4, (/nn/))
      wf4 = reshape(tmp2+tmp5, (/nn/))
      wf5 = reshape(tmp2+tmp5, (/nn/))
      !wf3 = reshape(tmp4, (/nn/))
      !wf4 = reshape(tmp5, (/nn/))
      aa = sum(abs(wf1-wf2))*dv
      bb = sum(abs(wf1-wf3))*dv
      cc = sum(abs(wf1-wf4))*dv
      dd = sum(abs(wf1-wf5))*dv
      nrm0 = sum(conjg(wf1)*wf1)*dv
      nrm1 = sum(conjg(wf2)*wf2)*dv
      nrm2 = sum(conjg(wf3)*wf3)*dv
      nrm3 = sum(conjg(wf4)*wf4)*dv
      nrm4 = sum(conjg(wf5)*wf5)*dv
      write(6,*) 'rank, js, aa, bb, cc, dd',rank,  js, real(aa), real(bb), real(cc), real(dd)
      write(6,*) 'rank, nrm0,1,2,3,4 ',rank,  nrm0, nrm1, nrm2, nrm3, nrm4
      !close(310+js)
      !open(unit=310, file='t_bar_2_debug')
      !do i=1,nn
      !  write(310+js,*) i, real(wf1(i)), real(wf2(i)), real(wf1(i)-wf2(i))!, real(wf3(i))
      !  !write(310,*) i, real(wf1(i)-wf2(i))
      !enddo
      !close(310+js)
    enddo
    !close(310)
    !open(unit=310, file='t_bar_2_debug')
    !do i=1,nn
    !  write(310,*) i, real(wf1(i)), real(wf2(i)), real(wf1(i)-wf2(i))!, real(wf3(i))
    !  !write(310,*) i, real(wf1(i)-wf2(i))
    !enddo
    !close(310)
  !------------------------------------
    !write(6,*) 'check orthogonality on grid for different atoms'
    !do ia=1,natom
    !  do ja=1,natom
    !    write(6,*) 'ia, ja', ia, ja
    !    jt=atom_map(ja)
    !    ixs=at(ja)%ir_start(1);iys=at(ja)%ir_start(2)
    !    izs=at(ja)%ir_start(3)
    !    ixe=ixs+p(jt)%nrough(1)-1;iye=iys+p(jt)%nrough(2)-1
    !    ize=izs+p(jt)%nrough(3)-1
    !    it=atom_map(ia)
    !    !ms = p(it)%mstates

    !    allocate(mat_tmp(p(it)%mstates,p(jt)%mstates),stat=stat)
    !    if(stat/=0) stop 'mat_tmp alloc problem in prep_tbar'

    !    do is=1,p(it)%mstates
    !      do js=1,p(jt)%mstates
    !        mat_tmp(is,js) = sum(conjg(at(ia)%local_c1_c(ixs:ixe,iys:iye,izs:ize,is))*&
    !                at(ja)%local_p3d1_c(:,:,:,js,1))*dv
    !      enddo
    !      write(6,*) real(mat_tmp(is,:))
    !    enddo
    !    deallocate(mat_tmp)
    !    write(6,*) '  '
    !  enddo
    !enddo

  !------------------------------------
    write(6,*) 'now checking e^sum vs e^sum approx 1 + sum'

    do js=1,nb
      wf_term = 0d0
      wf_term1 = 0d0
      tmp = phi_bar(:,:,:,js,ik)
      nrm0 = sum(conjg(tmp)*tmp)*dv

      call add_ctaylor_phi(tmp,tmp5,ik)
      call add_ctaylor_phi1(tmp,tmp16,ik)
      !call c_taylor_loop(tmp,tmp21,ik)
      !tmp1 = tmp
      tmp1 = 0d0
      tmp17 = 0d0
      tmp18 = 0d0
      tmp19 = 0d0
      do ia=1,natom
        it=atom_map(ia)
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        ms = p(it)%mstates

        allocate(ca(ms),stat=stat)
        if(stat/=0) stop 'ca alloc problem in prep_tbar'
        allocate(ca2(ms),stat=stat)
        if(stat/=0) stop 'ca2 alloc problem in prep_tbar'

        do is=1,ms; ca(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv; enddo
        !write(6,*) 'ia, ca', ia, ca
        !write(6,*) 'ia', ia
        !write(6,*) 'ixs, ixe', ixs, ixe
        !write(6,*) 'iys, iye', iys, iye
        !write(6,*) 'izs, ize', izs, ize

        if(.not.at(ia)%edge) then
          tmp20 = 0d0
          do iis=1,ms
            tmp20(ixs:ixe,iys:iye,izs:ize) = tmp20(ixs:ixe,iys:iye,izs:ize)+&
               at(ia)%local_p3d1_c(:,:,:,iis,ik)*sum(at(ia)%cbarmat(iis,:)*ca(:))
          enddo
          do is=1,ms; ca2(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp20)*dv; enddo
          !write(6,*) 'ia, ca2',ia, ca2
          do is=1,ms
            tmp1(ixs:ixe,iys:iye,izs:ize)=tmp1(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%c_taylor_terms(is,:,1)*ca(:))
            tmp19(ixs:ixe,iys:iye,izs:ize)=tmp19(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,ik)*sum(at(ia)%cbarmat(is,:)*ca2(:))/2d0
            tmp17(ixs:ixe,iys:iye,izs:ize)=tmp17(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cbarmat(is,:)*ca(:))
          enddo
        else
          do is=1,ms
            do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
              jx=mod(ix+nx-1,nx)+1
              jy=mod(iy+ny-1,ny)+1
              jz=mod(iz+nz-1,nz)+1
              tmp1(jx,jy,jz)=tmp1(jx,jy,jz)+at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)&
                  *sum(at(ia)%c_taylor(is,:)*ca(:))
            enddo;enddo;enddo
          enddo
        endif

        deallocate(ca)
        deallocate(ca2)
      enddo
      do ia=1,natom
        it=atom_map(ia)
        ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
        izs=at(ia)%ir_start(3)
        ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
        ize=izs+p(it)%nrough(3)-1
        ms = p(it)%mstates

        allocate(ca(ms),stat=stat)
        if(stat/=0) stop 'ca, y alloc problem in prep_tbar'
        do is=1,ms; ca(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp17)*dv; enddo
        !write(6,*) 'ia, ca_18', ia, ca
        if(.not.at(ia)%edge) then
          do is=1,ms
            tmp18(ixs:ixe,iys:iye,izs:ize)=tmp18(ixs:ixe,iys:iye,izs:ize)+&
              at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cbarmat(is,:)*ca(:))
          enddo
        endif
        deallocate(ca)
      enddo
      !write(6,*) 'i, diff 1 17', sum(abs(tmp17-tmp1))*dv
      !write(6,*) 'i, diff 1 19', sum(abs(tmp19-tmp1))*dv
      !write(6,*) 'i, diff 1 18', sum(abs(tmp18-tmp1))*dv
      !write(6,*) 'i, diff 18 19', sum(abs(tmp19-tmp18))*dv
      wf_term(:,:,:,1) = tmp17
      wf_term(:,:,:,1) = tmp18

      tmp3 = 0d0
      tmp4 = tmp
      tmp6 = tmp
      tmp7 = tmp
      tmp8 = 0d0
      tmp9 = tmp
      tmp10 = 0d0
      tmp11 = tmp
      tmp13 = tmp
      tmp14 = tmp
      tmp15 = tmp

      do i=1,tbar_taylor_power
        tmp8=0d0
        cin=tmp6
        call fft3d_forward(nx,ny,nz,cin,cout)
        cout=cout*ek(:,:,:,ik)
        call fft3d_backward(nx,ny,nz,cout,tmp6)
        call sn_phi(tmp6,tmp8,ik,-1d0)
        tmp6=tmp8-tmp6
        !cin=tmp6
        !call sum_inv_T_phi(cin,tmp6,ik)
        tmp6=-ci*dt*0.5d0/dble(i)*tmp6
        !wf_term(:,:,:,i) = tmp6
        !wf_term1(:,:,:,i) = tmp6
        !write(6,*) 'i tmp8', sum(conjg(tmp8)*tmp8)*dv
        !write(6,*) 'i tmp6', i, sum(conjg(tmp6)*tmp6)*dv
        tmp7=tmp7+tmp6

        if(i .eq. 1) then 
          cin = tmp
          call fft3d_forward(nx,ny,nz,cin,cout)
          cout=cout*ek(:,:,:,ik)
          call fft3d_backward(nx,ny,nz,cout,tmp10)
        else
          cin = tmp10
          call fft3d_forward(nx,ny,nz,cin,cout)
          cout=cout*ek(:,:,:,ik)
          call fft3d_backward(nx,ny,nz,cout,tmp10)
        endif
        !wf_term(:,:,:,i) = tmp10
        !write(6,*) 'i, wf_diff_term', i, sum(abs(wf_term(:,:,:,i)-wf_term1(:,:,:,i)))*dv
        tmp12=0d0
        tmp8=0d0
        tmp3=0d0
        !tmp16=0d0
        do ia=1,natom
          it=atom_map(ia)
          ixs=at(ia)%ir_start(1);iys=at(ia)%ir_start(2)
          izs=at(ia)%ir_start(3)
          ixe=ixs+p(it)%nrough(1)-1;iye=iys+p(it)%nrough(2)-1
          ize=izs+p(it)%nrough(3)-1
          ms = p(it)%mstates

          allocate(ca(ms),ca2(ms), ca3(ms), ca4(ms))
          if(stat/=0) stop 'ca alloc problem in prep_tbar'
          allocate(y(ms),stat=stat)
          if(stat/=0) stop 'y alloc problem in prep_tbar'
          y=at(ia)%sinv(:,ik)

          !write(6,*) 'ia', ia, ca
          !write(6,*) 'ia2', ia, ca2
    !      write(6,*) 'test1', y
          if(.not.at(ia)%edge) then
            !write(6,*) 'i', i
            do is=1,ms
              !ca(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv
              ca2(is) = sum(conjg(at(ia)%local_c_c(:,:,:,is))*tmp13)*dv
              ca3(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp14)*dv
              ca4(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv
            enddo
            call proj1(ia, tmp10, ca, ms, ik)
            !write(6,*) 'ia,i ca3', ia, i, ca3
            do is=1,ms
              tmp3(ixs:ixe,iys:iye,izs:ize)=tmp3(ixs:ixe,iys:iye,izs:ize)+&
                at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%cbarmat(is,:)*ca3(:))/dble(i)
              tmp8(ixs:ixe,iys:iye,izs:ize)=tmp8(ixs:ixe,iys:iye,izs:ize)+&
                (-ci*dt*0.5d0)/dble(i)*at(ia)%local_p3d1_c(:,:,:,is,1)*y(is)*ca2(is)
              tmp12(ixs:ixe,iys:iye,izs:ize)=tmp12(ixs:ixe,iys:iye,izs:ize)+&
                (-ci*dt*0.5d0)/dble(i)*at(ia)%local_p3d1_c(:,:,:,is,1)*y(is)*ca(is)
              !tmp16(ixs:ixe,iys:iye,izs:ize)=tmp16(ixs:ixe,iys:iye,izs:ize)+&
              !  at(ia)%local_p3d1_c(:,:,:,is,1)*sum(at(ia)%c_taylor_terms(is,:,i)*ca4(:))
            enddo
            !tmp13(ixs:ixe,iys:iye,izs:ize)=tmp12(ixs:ixe,iys:iye,izs:ize)
            !write(6,*) 'ia, i, tmp3', ia, i, sum(conjg(tmp3)*tmp3)*dv
            !write(6,*) 'ia, i, ca-c2', ia, i, maxval(real(ca-ca2)), maxval(aimag(ca-ca2))
            !write(6,*) 'ia, i, ca', ia, i, ca(1:2)
            !aaa(i) = aaa(i) + sum(conjg(tmp3)*tmp3)*dv
          else
            if(i .eq. 1) then 
              do is=1,ms
                ca3(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp)*dv
              enddo
              !call proj1(ia,tmp,ca,ms,ik)
            else
              do is=1,ms
                ca3(is) = sum(conjg(at(ia)%local_c1_c(:,:,:,is))*tmp3)*dv/dble(i)
              enddo

              !call proj1(ia,-ci*dt/2d0/dble(i)*tmp6,ca,ms,ik)
            endif
            do iz=izs,ize;do iy=iys,iye;do ix=ixs,ixe
              do is=1,ms
                tmp3(jx,jy,jz)=tmp3(jx,jy,jz)+at(ia)%local_p3d1_c(ix-ixs+1,iy-iys+1,iz-izs+1,is,1)&
                    *sum(at(ia)%cbarmat(is,:)*ca3(:))
              enddo
            enddo;enddo;enddo
          endif
          deallocate(ca, ca2, ca3, ca4, y)
        enddo
        !write(6,*) 'i, tmp8-tmp12 ', i, maxval(real(tmp8-tmp12)), maxval(aimag(tmp8-tmp12))
        tmp14=tmp3
        tmp4=tmp4+tmp14
        tmp13=tmp8
        tmp9=tmp9+tmp13
        tmp10=tmp12
        tmp11=tmp11+tmp12
        !tmp15=tmp15+tmp16
        wf_term1(:,:,:,i) = tmp3
      enddo

      !do i=1,tbar_taylor_power
      !  write(6,*) 'i, test c^2',i,  sum(abs(wf_term(:,:,:,i)-wf_term1(:,:,:,i)))*dv
      !enddo
      !write(6,*) 'i, test c^2',i,  sum(abs(wf_term(:,:,:,1)-wf_term1(:,:,:,1)))*dv
      wf1 = reshape(tmp1, (/nn/)) ! same as add_C_taylor
      wf2 = reshape(tmp, (/nn/))  !intput
      wf3 = reshape(tmp4, (/nn/)) !e^sum = 1 + sum + ...
      wf4 = reshape(tmp5, (/nn/)) ! add_C_taylor
      wf5 = reshape(tmp7, (/nn/)) !e^-idt/2S^-1T=1-idt/2S^-2T...
      wf6 = reshape(tmp9, (/nn/)) !
      wf7 = reshape(tmp11, (/nn/)) !
      wf8 = reshape(tmp16, (/nn/)) !
      wf9 = reshape(tmp21, (/nn/)) !
      
      !nrm1 = sum(conjg(wf1)*wf1)*dv
      !nrm2 = sum(conjg(wf4)*wf4)*dv
      !nrm3 = sum(conjg(wf3)*wf3)*dv
      !nrm4 = sum(conjg(wf5)*wf5)*dv
      aa = sum(abs(wf5-wf1))*dv !exact taylor - add_C_taylor
      bb = sum(abs(wf5-wf4))*dv
      cc = sum(abs(wf5-wf3))*dv
      dd = sum(abs(wf5-wf6))*dv
      ee = sum(abs(wf5-wf7))*dv
      ff = sum(abs(wf5-wf8))*dv
      !gg = sum(abs(wf5-wf9))*dv
      !write(6,*) 'js, wf1', js, maxval(real(wf1)), maxval(aimag(wf1))
      write(6,*) 'rank, js, aa, bb, cc, :', rank, js, real(aa), real(bb), real(cc)
      write(6,*) 'rank, js, dd, ee, ff:', rank, js, real(dd), real(ee), real(ff)
      !write(6,*) 'js, gg:', js, real(gg)
      !write(6,*) 'js, nrm0,1,2,3,4:', js, nrm1, nrm2, nrm3, nrm4
      !write(6,*) 'js, wf3', js, maxval(real(wf3)), maxval(aimag(wf3))
      !write(6,*) 'js, diff e', js, aa
    enddo

    write(6,*) 'now checking e^s^-1t vs split step'

    w1time = 0d0
    do js=1,nb
      tmp = phi_bar(:,:,:,js,ik)
      tmp1 = tmp
      tmp2 = tmp
      do i=1,tbar_taylor_power
        call st_phi(-ci*dt/dble(i)*tmp1, tmp3, ik)
        !call t_phi(-ci*dt/dble(i)*tmp1, tmp3, ik)
        !call sum_inv_T_phi(-ci*dt/dble(i)*tmp1, tmp3, ik)
        !aa = aa/dble(i)
        tmp2 = tmp2 +tmp3
        tmp1 = tmp3
      enddo
      w2time = mpi_wtime()
      call prop_tbar_split_sh(tmp, tmp4, ik)  
      w2time = mpi_wtime() - w2time
      w1time = w1time + w2time
     
      wf1 = reshape(tmp2, (/nn/))
      wf3 = reshape(tmp4, (/nn/))
      nrm0 = sum(conjg(wf1)*wf1)*dv
      nrm1 = sum(conjg(wf3)*wf3)*dv
      aa = sum(abs(wf1-wf3))*dv
      write(6,*) 'rank, js, aa, nrm0, nrm1:', rank, js, real(aa), nrm0, nrm1
    enddo
    call sync_mpi
    write(6,*) 'rank, walltime apply tbar', rank, w1time
  end subroutine debug_tbar_sh_multi
end subroutine prep_tbar_sh
