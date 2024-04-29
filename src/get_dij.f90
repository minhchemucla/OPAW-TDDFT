subroutine get_dij
    use param_mod, only : param_dij_max ,k_flg
    use mpi_lib_ours
    use main_mod 
    use libpaw_mod
    use atom_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    use paw_mod
    implicit none

    integer :: ia   
!    if(k_flg) then 
      if(rank==0) then
        call set_rhoij
    !    call rd_input

        call pawdenpot(compch_sph,epaw,epawdc,0,ixc,&
            natom,natom,nspden,ntypat,nucdipmom,0,0,paw_an,paw_an,&
            paw_ij,pawang,pawprtvol,pawrad,pawrhoij,0,pawtab,xcdev,1.0d0,xclevel,1d-14,ucvol,znucl)
    
    !    write(461,*) 'rfgd'
    !    write(461,*) pawfgrtab(1)%rfgd
    !    call pawgylm(pawfgrtab(1)%gylm,pawfgrtab(1)%gylmgr,pawfgrtab(1)%gylmgr2,&
    !        paw_an(1)%lm_size,pawfgrtab(1)%nfgd,1,1,1,pawtab(1),pawfgrtab(1)%rfgd)
    
        call pawdij(cplex,enunit,gprimd,ipert,natom,natom,nfft,nfftot,&
    &     nspden,ntypat,paw_an,paw_ij,pawang,pawfgrtab,pawprtvol,&
    &     pawrad,pawrhoij,pawspnorb,pawtab,xcdev,qphon,spnorbscl,&
    &     ucvol,charge,vks,vxc,xred)
    
        call transfer_dij
      endif

      do ia=1,natom
        call bcast_r8(at(ia)%dij,size(at(ia)%dij),0)
      enddo
!    else
!        call set_rhoij
!
!        call pawdenpot(compch_sph,epaw,epawdc,0,ixc,&
!            natom,natom,nspden,ntypat,nucdipmom,0,0,paw_an,paw_an,&
!            paw_ij,pawang,pawprtvol,pawrad,pawrhoij,0,pawtab,xcdev,1.0d0,xclevel,1d-14,ucvol,znucl)
!    
!        call pawdij(cplex,enunit,gprimd,ipert,natom,natom,nfft,nfftot,&
!    &     nspden,ntypat,paw_an,paw_ij,pawang,pawfgrtab,pawprtvol,&
!    &     pawrad,pawrhoij,pawspnorb,pawtab,xcdev,qphon,spnorbscl,&
!    &     ucvol,charge,vks,vxc,xred)
!    
!        call transfer_dij
!    endif

contains
    subroutine set_rhoij
        implicit none

        integer :: ndim,ndim1,i,j,ind,nselect

        do ia=1,natom
            pawrhoij(ia)%rhoijp=0d0
            pawrhoij(ia)%rhoijselect=0d0
            pawrhoij(ia)%nrhoijsel=0
            ndim=size(at(ia)%rhoij,1)
            ndim1=size(pawrhoij(ia)%rhoijp,1)
            if (ndim*(ndim+1)/2/=ndim1) then
                write(*,*) 'dim of rhoij not match,ia,ndim,ndim1',ia,ndim,ndim1
                stop
            endif
            
            ind=0

            do j=1,ndim
                do i=1,j
                    ind=ind+1
                    pawrhoij(ia)%rhoijp(ind,1)=at(ia)%rhoij(i,j)
                enddo
            enddo

            nselect=0
            do i=1,ndim1
                if(abs(pawrhoij(ia)%rhoijp(i,1))>1d-10) then
                    nselect=nselect+1
                    pawrhoij(ia)%rhoijselect(nselect)=i
                    pawrhoij(ia)%rhoijp(nselect,1)=pawrhoij(ia)%rhoijp(i,1)
                endif
            enddo

            pawrhoij(ia)%nrhoijsel=nselect
!        write(600,*) pawrhoij(ia)%rhoijp
        enddo

    end subroutine

!    subroutine rd_input
!        implicit none
!
!        read(20,*) vxc_r
!        read(21,*) vtrial_r
!    end subroutine        
    
    subroutine transfer_dij
        implicit none

        integer :: i,j,ind,ndim,ndim1
        real*8 dd ! DN

        siz_dijall=0
        do ia=1,natom
            ndim=size(at(ia)%rhoij,1)
            ndim1=size(paw_ij(ia)%dij,1)
            if (ndim*(ndim+1)/2/=ndim1) then
                write(*,*) 'dim of dij not match,ia,ndim,ndim1',ia,ndim,ndim1
                stop
            endif
            siz_dijall=siz_dijall+ndim1
        enddo

        if(.not.allocated(dijall)) then
            allocate(dijall(siz_dijall),stat=i)
            if(i/=0) stop 'dijall alloc problem'
        endif
       
        !write(*,*) 'siz_dijall',siz_dijall 

        siz_dijall=1
        do ia=1,natom
!            at(ia)%dij(i,j)=min(paw_ij(ia)%dij(ind,1),100d0)
!            at(ia)%dij(j,i)=min(paw_ij(ia)%dij(ind,1),100d0)
            ndim1=size(paw_ij(ia)%dij,1)
            dijall(siz_dijall:siz_dijall+ndim1-1)=paw_ij(ia)%dij(:,1)
            siz_dijall=siz_dijall+ndim1
        enddo

        if (siz_dijall/=size(dijall)+1) stop 'siz dij problem' 
        if(diis_dij .and. k_flg) call diis1_general(dijall,size(dijall))

        siz_dijall=1
        do ia=1,natom
            ndim=size(at(ia)%rhoij,1)
            ndim1=size(paw_ij(ia)%dij,1)
            do j=1,ndim
                do i=j,ndim
                   ind=i*(i-1)/2+j
                   dd = dijall(siz_dijall+ind-1) !DN
                   if(abs(dd)>param_dij_max) dd = dd /abs(dd)*param_dij_max !DN 
                   at(ia)%dij(i,j)=dd ! DN 
                   at(ia)%dij(j,i)=dd ! DN
                   !if(i==ndim.and.j==ndim) then
                   !   write(6,*)' for atom ',ia,' dij_18_18 debug ',dd
                   !endif
                enddo
            enddo
            siz_dijall=siz_dijall+ndim1
        enddo

    end subroutine
end subroutine get_dij      
