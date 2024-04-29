subroutine findtype(z,chz)
    implicit none

    character*2      :: chz
    integer          :: z,i
    character*2      :: atoms(118)
    logical          :: atom_found

    atoms( 1: 2)=&
      (/'H ',                                                                                'He'/)
    atoms( 3:10)=& 
      (/'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne'/)
    atoms(11:18)=&
      (/'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar'/)
    atoms(19:36)=&
      (/'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr'/)
    atoms(37:54)=&
      (/'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe'/)
    atoms(55:56)=&
      (/'Cs','Ba'/)
    atoms(72:86)=&
                     (/'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'/)
    atoms(87:88)=&
      (/'Fr','Ra'/)
    atoms(104:118)=     (/'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/)

    atoms(57:71)=       (/'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'/)
    atoms(89:103)=      (/'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/)

    call to_upper(chz(1:1))
    call to_lower(chz(2:2))

    atom_found=.false.

    do i=1,118
        if(chz==atoms(i)) then
            z=i
            atom_found=.true.
        endif
    enddo

    if(.not.atom_found) then
        write(*,*) 'Error : atom type not found for atom with name: ',chz
        stop
    endif
end subroutine findtype

subroutine to_lower(str)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character :: str
     integer :: i,j

     do i = 1, len(str)
          j = iachar(str(i:i))
          if (j>= iachar("A") .and. j<=iachar("Z") ) then
               str(i:i) = achar(iachar(str(i:i))+32)
          else
               str(i:i) = str(i:i)
          end if
     end do

end subroutine to_lower

subroutine to_upper(str)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page

     implicit none

     character :: str
     integer :: i,j

     do i = 1, len(str)
          j = iachar(str(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               str(i:i) = achar(iachar(str(i:i))-32)
          else
               str(i:i) = str(i:i)
          end if
     end do

end subroutine to_upper
