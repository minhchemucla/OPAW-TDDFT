program main
    ! Any referenced equations are from "Implementation of the projector 
        !augmented-wave method in the ABINIT code: Application to the 
        !study of iron under pressure"
        ! or
        !"Real Space Orthogonal Projector-Augmented-Wave Method"
        !Use context to figure out which paper. - Minh
    use param_mod, only : iscf_start,iscf_hminmax,k_flg,flg_bin
    use libpaw_mod
    use main_mod
    use tddft_mod 
    use mpi_lib_ours
    use paw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    implicit none

    integer :: i,ib,ik,ia
    real*8, allocatable :: dens3d(:,:,:)
    real*8 :: tstart, tend, tstart_tddft, tend_tddft
    real*8 :: ts, te
    real*8 :: wtime

    iscf=0d0
    call cpu_time(tstart)
    wtime = mpi_wtime()
    call prepare_mpi_lib
    call read_input

    call ncpaw_libpaw_prepare

    call get_rnel

    if (tddft_flag < 1) then
      allocate(phit(nx,ny,nz,nb,nk_loc))
      !phit are the pseudo wfs
      if(.not.k_flg) then
         if(rank==0) then
             allocate(phit_tot(nx,ny,nz,nb*nodes,nk_loc))
             !for parallelization
         else
              !phit_tot is stored only on rank=3
             allocate(phit_tot(0,0,0,0,0))
         endif
      endif
      call prep_phi0 !fills phit or phit_tot with random numbers

  !    if(k_flg) then
      call get_nhat
      if(rank==0)write(*,*) 'sum nhat0',sum(nhat)*dv
      allocate(dens3d(nx,ny,nz))
      dens3d=reshape(dens,(/nx,ny,nz/))
      dens3d=dens3d-nhat
      dens=reshape(dens3d,(/nn/))
  !    else
  !        call get_nhat
  !        write(*,*) 'sum nhat0',sum(nhat)*dv
  !        dens3d=reshape(dens,(/nx,ny,nz/))
  !        dens3d=dens3d-nhat
  !        dens=reshape(dens3d,(/nn/))
  !    endif
  !    call debug_phi_dens
  !    call get_sij
  !    call test_sss


      do iscf=1,nscf
!call cpu_time(ts)
          call ncpaw_make_hamiltonian
!call cpu_time(te)
!if(rank==0) write(*,*) 'iscf, ncpaw_make_hamiltonian', iscf, te-ts
!call cpu_time(ts)
          call set_hminmax
!call cpu_time(te)
!if(rank==0) write(*,*) 'iscf, sethminmax            ', iscf, te-ts
           
  !        call get_h_debug 
          !write(*,*) 'rank,hmax,mu',rank,hmax,mu
!call cpu_time(ts)
          call chebyr_m(hmax,mu)
!call cpu_time(te)
!if(rank==0) write(*,*) 'iscf, chebyr_m              ', iscf, te-ts
!call cpu_time(ts)
          call orthog_phi
!call cpu_time(te)
!if(rank==0) write(*,*) 'iscf, orthog_phi            ', iscf, te-ts
!call cpu_time(ts)
          call diagh
!call cpu_time(te)
!if(rank==0) write(*,*) 'iscf, diagh                 ', iscf, te-ts
  !        call get_h_debug 
  !      if(iscf==(nscf-1)) then
  !        call writewf_debug
          !call writewfbar_debug
        !endif
      enddo
      !call test_eig !minh
      !if(h_type.eq.0) call writewfbar

      call writewf
      !call writewfbar
    endif

    !call calc_exx_homo
     
    call cpu_time(tend_tddft)
    call sync_mpi
    write(6,*) 'rank, cpu total time for DFT portion: ', rank, tend_tddft-tstart
    !if (tddft_flag > -1) call opaw_tddft_old
    if (tddft_flag > -1) then
      tstart_tddft = tend_tddft
      call opaw_tddft
      call cpu_time(tend_tddft)
      write(6,*) 'rank, cpu total time for TDDFT portion: ', rank, tend_tddft-tstart_tddft
    endif
  
    call cpu_time(tend)
    wtime= mpi_wtime()-wtime
    call sync_mpi
    write(6,*) 'rank, cpu total time: ', rank, tend-tstart
    call sync_mpi
    if(rank==0) write(6,*) 'approximate total wall time for rank==0: ', wtime
    call finalize_mpi_lib
contains  
    subroutine check_sij
        implicit none

        integer :: ia,it
        write(*,*) 'check sij'
        do ia=1,natom
            write(*,*) 'atom',ia
            it=atom_map(ia)
            write(*,*) at(ia)%s
            write(*,*) p(it)%sij
        enddo
    end subroutine

    subroutine debug_rhoij
        implicit none
        integer :: ia,is,it

        if(nodes>1) stop 'debug rhoij, nodes should be 1'

        open(unit=1,file='rhoij')
        do ia=1,natom
            read(1,*) at(ia)%rhoij
        enddo
        close(1)
    end subroutine 

    subroutine debug_phi_dens
        implicit none

        real*8 :: tmp(nx,ny,nz)
        integer :: is,it

        if(nodes/=nkpt) then
            write(*,*) 'for debugging, let nodes = nkpt'
            stop
        endif

        phit=0d0
        nb=0
        do while(.true.)
            read(90,*,end=12)
            nb=nb+1
        enddo
12      rewind(90)
!
        do ib=1,nb
        read(90,*) tmp
            phit(:,:,:,ib,1)=tmp
        enddo

        call acc_dens(dens,nx,ny,nz)
        call acc_coeff

        call allsum_r8(dens,size(dens))
        do ia=1,natom
            call allsum_r8(at(ia)%rhoij,size(at(ia)%rhoij))
        enddo

        if(rank==0) then
            do ia=1,natom
                write(117,*) at(ia)%rhoij
            enddo
            write(117,*)
        endif
!        call debug_rhoij
        if(rank==0) then
            call get_nhat
            open(unit=9,file='nhat')
            write(9,*) nhat
            close(9)
        endif

        if(rank==0) then
            write(*,*) 'sum dens',sum(dens)*dv
            write(*,*) 'sum dens+nhat',sum(dens)*dv+sum(nhat)*dv
        endif

        call bcast_r8(nhat,size(nhat),0)
!        write(101,*) dens
        call get_pot(dens)
        call get_dij

        call get_h_debug 
        call get_sij
        call finalize_mpi_lib
        stop
    end subroutine debug_phi_dens

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


    subroutine check_pij
        implicit none
        integer :: ia,ms,is,it,js
        write(*,*) 'checking <pi|pj>'
        do ia=1,natom
            write(*,*) 'atom',ia
            it=atom_map(ia)
            ms=p(it)%mstates
            do is=1,ms
                do js=1,ms
                    write(*,'(f9.5)',advance='no') sum(at(ia)%local_p3d(:,:,:,is)*&
                        at(ia)%local_p3d(:,:,:,js))*dv

                enddo
                write(*,*)
            enddo
        enddo
    end subroutine
    
    subroutine read_input
        implicit none
        open(unit=13,file='dim_paw.inp')

        if (rank==0) then
            ! number of grid points
            call fetch_i (nx,'nx',13)
            call fetch_i (ny,'ny',13)
            call fetch_i (nz,'nz',13)

            !length of the box that encloses system
            call fetch_r (box_x,'box_x',13)
            call fetch_r (box_y,'box_y',13)
            call fetch_r (box_z,'box_z',13)

            call fetch_L (periodic,'periodic',13)

            call fetch_r (mu0    , 'mu',13)
            !dmua is a small shift of the chemical potential to 
            !ensure Filter includes all occupied states
            call fetch_r (dmua   , 'dmua',13)

            call fetch_r (mix_diis , 'mix_diis' ,13)
            call fetch_r (mix_diis1, 'mix_diis1',13)

            if(mix_diis1<0d0) diis_dij=.false.

            call fetch_i (nscf, 'nscf',13)
            !Functional 0=lda, 1=pbe
            call fetch_i (funct,'funct',13)

            call fetch_i (iscf_hminmax,'iscf_hminmax',13)
            call fetch_i (iscf_start,  'iscf_start'  ,13)
            call fetch_r (hmax_0,'hmax',13)
            call fetch_r (hmin_0,'hmin',13)

            !finegrid spacing
            call fetch_r (p_fg, 'finegrid',13)

            call fetch_r (ekcut, 'ekcut',13)
            !
            call fetch_i(h_type, 'h_type',13)
            call fetch_l(flg_bin   , 'flg_bin',13)
            ! -1 = no tddft, 0 tddft with no dft,tddft>0= tddft and dft
            call fetch_i(tddft_flag, 'tddft',13)

            if (tddft_flag < -1) stop 'tddft: -1=dft, 0=dft+tddft, >0 tddft'
            if (tddft_flag > -1) then
              call fetch_i (nt, 'nt', 13)
              call fetch_r (dt, 'dt', 13)
              call fetch_i(ipol, 'exct_pol',13)
              call fetch_r(sm, 'strength ',13)
              call fetch_i(prop_type, 'prop_type',13)
              !call fetch_l(prop_bar, 'prop_bar',13)
              call fetch_i(debug_prop_flag, 'debug_prop',13)
              call fetch_i(dyn, 'theory',13)
              call fetch_i(n_restart, 'n_restart',13)
              !call fetch_i(nvirtual_read, 'nvirtual',13)
              !call fetch_r(pfrozen, 'percent_frozen',13)

              !if (pfrozen .ge. 1d0) stop 'please make percent_frozen less than 1'
              !if (pfrozen .lt. 0d0) stop 'please make percent_frozen a positive number'
              !pfrozen = 0d0 !frozen orbital propagation not making for now
            endif

            if(h_type/=0 .and. h_type/=1) then
                write(*,*) 'h_type should be 0 or 1'
                stop
            endif

            if(funct/=0 .and. funct/=1) then
                write(*,*) 'funct should be 0 or 1'
                stop
            endif
            mu0=mu0+dmua
   
            dx=box_x/nx 
            dy=box_y/ny 
            dz=box_z/nz 
            write(*,*) 'dx,dy,dz',dx,dy,dz
            nn=nx*ny*nz
            dv=dx*dy*dz
           
            !if((.not. k_flg) .and. mod(nn,nodes)/=0) then
            if((.not. k_flg) .and. mod(nn,nodes)/=0 .and. tddft_flag<1) then
                write(*,*) 'parallelize over bands, nn should be multiple of #. nodes'
                write(*,*) 'nodes,nn',nodes,nn
                stop
            endif

            
            nfovnr=max(ceiling(dx/p_fg),1)
            nfovnr=max(nfovnr,ceiling(dy/p_fg))
            nfovnr=max(nfovnr,ceiling(dz/p_fg))
            write(*,*) 'nfovnr',nfovnr

            kkx=6.283185307180/(nx*dx)
            kky=6.283185307180/(ny*dy)
            kkz=6.283185307180/(nz*dz)
    
            if(periodic) scale_vh=1
            if(.not.periodic) scale_vh=2
            xmax=dble(nx)*dx/2d0
            ymax=dble(ny)*dy/2d0
            zmax=dble(nz)*dz/2d0
            write(*,*) 'xmax,ymax,zmax',xmax,ymax,zmax

            if(ekcut > 0d0) ekread=.true.

            close(13)
        endif

        call bcast_scalar_i(nx)
        call bcast_scalar_i(ny)
        call bcast_scalar_i(nz)
        call bcast_scalar_i(scale_vh)
        call bcast_scalar_i(nscf)
        call bcast_scalar_i(nn)
        call bcast_scalar_i(funct)
        call bcast_scalar_i(nfovnr)
        call bcast_scalar_i(iscf_hminmax)
        call bcast_scalar_i(iscf_start)

        call bcast_scalar_l(periodic)
        call bcast_scalar_l(diis_dij)

        call bcast_scalar_r8(p_fg)
        call bcast_scalar_r8(dx)
        call bcast_scalar_r8(dy)
        call bcast_scalar_r8(dz)
        call bcast_scalar_r8(mu0)
        call bcast_scalar_r8(dmua)
        call bcast_scalar_r8(mix_diis1)
        call bcast_scalar_r8(mix_diis)
        call bcast_scalar_r8(dv)
        call bcast_scalar_r8(kkx)
        call bcast_scalar_r8(kky)
        call bcast_scalar_r8(kkz)
        call bcast_scalar_r8(xmax)
        call bcast_scalar_r8(ymax)
        call bcast_scalar_r8(zmax)

        call bcast_scalar_r8(hmin_0)
        call bcast_scalar_r8(hmax_0)
        call bcast_scalar_r8(ekcut)
        call bcast_scalar_i(h_type)
        call bcast_scalar_i(tddft_flag)
        call bcast_scalar_i(nt)
        call bcast_scalar_r8(dt)
        call bcast_scalar_i(ipol)
        call bcast_scalar_r8(sm)
        !call bcast_scalar_l(prop_bar)
        call bcast_scalar_i(prop_type)
        call bcast_scalar_i(debug_prop_flag)
        call bcast_scalar_i(dyn)
        call bcast_scalar_l(flg_bin)
        call bcast_scalar_i(n_restart)
        !call bcast_scalar_i(nvirtual_read)
    end subroutine
end program main      
