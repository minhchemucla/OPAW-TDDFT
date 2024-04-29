subroutine tddft_get_nb
  use main_mod
  use tddft_mod
  use mpi_lib_ours
  implicit none
  integer :: ns,nds,st,is,ip

  ip = 76000+rank
  nds = max(nodes,1)

  !some cores could be propagating virtual states. I could proceed as before by restricting the number of 
  !cores so that this doesn't happen but this would trade off the potential extra parallelization of
  !slicing the grid points. It's easier to code it this way as well so I'll just deal with it for now.
  !- Minh

  if(rank==0) then
    !if (mod(nocc+nvirtual,nodes) .ne. 0) stop 'please adjust nodes so mod(nocc+nvirtual,nodes)=0'
    if (mod(nocc,nds) .ne. 0 .and. nds < nocc) then
      ns = nocc/nds !should always round down
      do nvirtual=1,nstates-ns*nds+1
        if(mod(nocc+nvirtual,nds) .eq. 0) goto 10
      enddo
      write(6,*) 'nvirtual too big. nvirtual :', nvirtual
      write(6,*) 'nstates', nstates
      stop
    else if (mod(nocc,nds) .ne. 0 .and. nds >= nocc) then
      nvirtual = nds - nocc 
    else
      ns = 0
      nvirtual = 0
    endif
10  write(6,*) 'nocc, nvirtual, nodes, ns, nstates', nocc, nvirtual, nds, ns, nstates
    if(nocc+nvirtual>nstates) stop 'nocc+nvirtual > nstates'
    nb=(nocc+nvirtual)/nds
    !if (nvirtual .ge. nb) stop 'nvirtual > states per nodes i.e one or more nodes are propagating &
    !   virtual states unessecarily. Please reduce number of cores'
    write(6,*) 'adjusted total number of states to be propagated per node for tddft is',nb
  endif

  call bcast_scalar_i(nb)
  call bcast_scalar_i(nvirtual)
  allocate(state_map(nb),stat=st); if(st/=0) stop 'state_map'

  do is=1,nb
    state_map(is) = is + rank*nb 
  enddo
  !call bcast_scalar_i(nfrozen)
  !call bcast_scalar_i(nprop)
  write(ip,*) 'state_map, nb',nb
  do is=1,nb
    write(ip,*) is, state_map(is)
  enddo
  call flush(ip)

  !stop
end subroutine tddft_get_nb

