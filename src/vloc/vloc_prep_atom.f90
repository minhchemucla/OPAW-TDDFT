!subroutine  vloc_prep_atom(ch, l, rr, vr, nr, nnr)
subroutine  vloc_prep_atom(ma, ch, l, rr, vr, nr, nnr) !minh - I added ma 
!  use main, only : pp_multiple
  implicit none
  integer nnr,nr,lp,ir,jr,l,ch, ma
  real*8  vr(nnr),r,p,v,rr(nnr)

!  if(pp_multiple) then
     call vloc_prep_atom_upf
!  else
!     call vloc_prep_atom_fort5000
!  endif

contains
  subroutine vloc_prep_atom_upf
    use atom_mod
    implicit none
    !integer ma
    logical fit

    fit=.false.
    !do ma=1,ntype !minh comment out
       if(atom_Z(ma)==ch) then
          fit=.true.
          !write(6,*) 'MINH PAWINFO', ma, pawinfo(ma)%nr, nr, atom_Z(ma), ch
          call check(pawinfo(ma)%nr,nr,        ' nrpp_ma, nr     ')
          call check(size(pawinfo(ma)%vloc), nnr,' sz_vploc_1, nnr ')
          call check_le(nr,   nnr,       ' nr, nnr         ')
          vr(1:nr) = pawinfo(ma)%vloc
       endif
    !enddo !minh comment out
    if(.not.fit) stop ' error: didnt read vpploc '
  end subroutine vloc_prep_atom_upf

  subroutine vloc_prep_atom_fort5000
    implicit none
    rewind(ch+5000)
    read(ch+5000,*) 
    do lp=0,l
       read(ch+5000,*)nr; call check_le(nr,nnr,' nr nnr   ')
       do ir=1,nr
          read(ch+5000,*)jr,r,p,v; call check(jr,ir,' jr vs.ir');
          if(lp==l) then
             rr(ir)=r; vr(ir) = v 
             !write(798,*)rr(ir),vr(ir) ! removed
          end if
       enddo
    enddo
  end subroutine vloc_prep_atom_fort5000

end subroutine vloc_prep_atom
