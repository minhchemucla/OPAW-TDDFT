subroutine allocate_dijfock
  use main_mod, only : nocc
  use atom_mod
  use atom_mod, only : p=>pawinfo, at=> atominfo
  use libpaw_mod 
  implicit none
  integer :: ia,st,it,ms
  
  do ia=1,natom
    it = atom_map(ia)
    ms = p(it)%mstates
    if(.not. allocated(at(ia)%dijfockhat)) then
      allocate(at(ia)%dijfockhat(ms,ms,nocc,nocc),stat=st)
      if(st/=0) then
        write(6,*) 'allocate error at(ia)%dijhat'
        write(6,*) 'ia:', ia
        stop
      endif
    endif

    if(.not. allocated(at(ia)%dijfock)) then
      allocate(at(ia)%dijfock(ms,ms),stat=st)
      if(st/=0) then
        write(6,*) 'allocate error at(ia)%dijfock'
        write(6,*) 'ia:', ia
        stop
      endif
      allocate(at(ia)%dijfock_vv(ms,ms),stat=st); if(st/=0) stop 'vv'
      allocate(at(ia)%dijfock_cv(ms,ms),stat=st); if(st/=0) stop 'cv'
    endif
  enddo
end subroutine

