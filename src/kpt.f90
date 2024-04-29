subroutine prep_kpt
    use mpi_lib_ours
    use main_mod
    use param_mod, only : k_flg
    implicit none

    if(k_flg) then
      open(unit=19,file='kpoints')
      call count_kpt
      call alloc_kpt
      call read_kpt
      close(19)
    else
      ! uses kpt = (0,0,0) in Brilluoin zone
      nk_loc=1
      nkpt=nodes
      allocate(wk(nkpt),kpt(3,nkpt))
      allocate(mu(nk_loc),hmin(nk_loc),hmax(nk_loc))
      wk=1d0
      kpt=0d0
    endif
contains
    subroutine read_kpt
        implicit none

        integer :: i
        if(rank==0) then
            do i=1,nkpt
                read(19,*) kpt(:,i),wk(i)
            enddo
            wk=wk/2d0
            kpt(1,:)=kpt(1,:)*kkx
            kpt(2,:)=kpt(2,:)*kky
            kpt(3,:)=kpt(3,:)*kkz
            write(*,*) 'kpoints'
            do i=1,nkpt
                write(*,*) 'i,kpt',i,kpt(:,i),wk(i)
            enddo
        endif

        call bcast_r8(kpt,size(kpt),0)
        call bcast_r8(wk,size(wk),0)

    end subroutine read_kpt

    subroutine alloc_kpt
        implicit none

        integer :: stat

        allocate(kpt(3,nkpt),wk(nkpt),stat=stat)
        if(stat/=0) stop 'kpt alloc problem'
        allocate(mu(nk_loc),hmin(nk_loc),hmax(nk_loc),stat=stat)
        if(stat/=0) stop 'mu alloc problem'
        mu=mu0
    end subroutine alloc_kpt

    subroutine count_kpt
        implicit none
        
        real*8 :: inp
        integer :: i,ik

        if(rank==0) then
            i=0
            do while (.true.)
                read(19,*,end=11) inp
                i=i+1
            enddo
11          rewind(19)
        
            nkpt=i

            if(nkpt<nodes) then
                write(*,*) 'nkpt<nodes',nkpt,nodes
                stop
            endif
        endif

        call bcast_scalar_i(nkpt)
        nk_loc=0
        do ik=rank+1,nkpt,nodes
            nk_loc=nk_loc+1
        enddo

    end subroutine    
end subroutine prep_kpt      
