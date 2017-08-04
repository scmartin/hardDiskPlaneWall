!*****************************************************************
!
!  This module is used in the creation of the initial coordinates
!  for the simulation. It is used by the planar_init modules.
!  It detects particle overlaps when attempting to insert
!  particles; if this function returns "1", the attempted insertion
!  is rejected.
!
!*****************************************************************
module overlap

contains

integer function ovrlap(coord,x,y,lx,npart,pcount)

IMPLICIT NONE

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15,307)
integer(dp) :: npart,pcount
real(dp) ::  coord(npart,2),r2,lx,delx,x,y
integer(dp) ::j

ovrlap = 0
do j = 1,pcount-1
    delx = coord(j,1) - x
    if(delx.gt.0.5d0*lx) then
        delx = delx - lx
    else  
        if(delx.lt.-0.5d0*lx) delx = delx + lx
    endif 
    r2 = delx**2 + (coord(j,2) - y)**2
    if(r2 <= 1.0d0) then
        ovrlap = 1
        return
    endif
enddo
return

end function ovrlap

end module overlap
