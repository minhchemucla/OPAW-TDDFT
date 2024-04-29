subroutine get_rnel
      use param_mod, only : rnel_1, rnel_2, rnel_3, k_flg
      use main_mod
      use atom_mod
      use mpi_lib_ours
      implicit none

      integer :: ia,it,nds

      nds = max(nodes,1)

      if(rank==0) then
          rnel=0d0
          do ia=1,natom
            it=atom_map(ia)
            rnel=rnel+pawinfo(it)%val
          enddo
          write(*,*) 'total valence charge is',rnel
          nocc=int(rnel/2d0+1d-8)
          write(*,*) 'number of occupied states is',nocc
          nb=max(ceiling((dble(nocc))*rnel_1)+rnel_2,rnel_3)
          write(*,*) 'total number of states is',nb
          if(.not.k_flg) nb=(nb-1)/nds+1
          write(*,*) 'adjusted total number of states per node is',nb
      endif

      call bcast_scalar_i(nb)
      call bcast_scalar_i(nocc)
      call bcast_scalar_r8(rnel)
      call sync_mpi
end subroutine get_rnel      

