!for an increasing array, find the position of target
!and find the weight associated with the two endpoints
!linear interpolation
subroutine interpolate(array,n,targ,i,j,wi,wj)
    implicit none

    integer :: i,j,n,jm !i and j are the index of the two closest radial grid  to the fine grid
    real*8  :: array(n)
    real*8  :: wi,wj,targ
    
    i=1
    j=n

    do while (j-i>1)
        jm=(i+j)/2
        if(targ.lt.array(jm)) then
            j=jm
        else
            i=jm
        endif
    enddo
    
    wi=(array(j)-targ)/(array(j)-array(i))
    wj=(targ-array(i))/(array(j)-array(i))

end subroutine interpolate      
