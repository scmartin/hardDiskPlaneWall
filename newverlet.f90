!************************************************************************!
!  This module has two subroutines. verlet() builds the neighbor list    !
!  and check is used to determine if an update is necessary.             !
!  The neighbor lists are constructed using two arrays. vlist holds      !
!  the indices of all neighbors of all particles. The point array holds  !
!  index of the first neighbor of the particle whose coordinate index    !
!  is the same as the point index. For example, if I want to know the    !
!  neighbors of particle 4 (coords(4,:)) I would look at point(4); the   !
!  number stored at point(4) tells me the index of the first neighbor of !
!  particle 4 in the vlist array. I then look up the number stored at    !
!  point(5). This tells me the index of the first neighbor of particle 5 !
!                                                                        !
!  TL;DR The neighbors of particle i are stored in                       !
!        vlist(point(i):point(i+1)-1)                                    !
!************************************************************************!

module verletmod

implicit none

CONTAINS

subroutine verlet(update,rv,npart,coords,ocoords,vlist,point,lx,listmax)
! creates new verlet neighbor lists

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15,307)
real(dp) :: r2,disx
integer(dp) :: listcount,i,j,npart,listmax
real(dp), dimension(npart,2) :: coords,ocoords
integer(dp), dimension(listmax) :: vlist
integer(dp), dimension(npart+1) :: point
real(dp) :: lx,rv
logical :: update

listcount = 0  ! listcount was used to verify that the length of the
               !   neighbor list array (vlist) was long enough to hold the full
               !   list of neighbors for all particles 
point = 0      ! point is an array which holds the indices for the first
               !   particle in each neighbor list
vlist = 0      ! vlist holds the indices for all neighbors
update = .false. ! update is used by main program to determine whether to call
                 !   this subroutine. update=.true. calls verlet()
ocoords = coords ! saves the coords at the point that the neighbor list is
                 !   created. used by check() subroutine to see how far
                 !   particles have moved since last neighbor list was generated

! The main loop to generate the neighbor lists. selects each particle and loops
! over all other particles to find those that are within the verlet cutoff
! radius

do i = 1,npart
  point(i) = listcount + 1
!   if (listcount == listmax) then
!     write (*,*) "listmax = ",listmax
!     write (*,*) "listcount =",listcount
!     write(6,'(1000020f10.5)') coords(:,1)
!     write(6,'(10000F10.5)') coords(:,1)
!     stop 'verlet list too short'
!   endif
  do j = 1,npart
    if (j/=i) then
      ! test for inclusion of particle j in verlet list of i
      disx = (coords(i,1)-coords(j,1))
      if (abs(disx) > lx*0.5d0) disx = disx - sign(lx,disx)
      r2 = disx**2+(coords(i,2)-coords(j,2))**2
      if (r2<rv**2) then
        listcount = listcount + 1
        vlist(listcount) = j
      endif
      if (r2<1) stop " a particle impact was missed!" !this if statement was
                           ! included to ensure that particle impacts are not
                           ! missed. Missed impacts indicate that the verlet
                           ! cutoff is too small relative to the step size
    endif
  enddo
  point(npart+1) = listcount + 1
enddo
!   print *, "coords "
!   write(6,'(20f10.5)') coords(:,1)
!   write(6,'(20F10.5)') coords(:,2)
!   print *, "point ",point
!   print *, "vlist ",vlist
!   print *, "liscount",listcount
end subroutine verlet

! check() subroutine calculates how far each particle has moved since the
! last neighbor list was generated. If any particle has moved more than half the
! minimum distance to impact a particle not previously considered a neighbor,
! check sets update=.true. and verlet() is called

subroutine check(npart,coords,ocoords,rv,update,lx2)

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15,307)
integer(dp) :: i,npart,n
real(dp), dimension(npart,2) :: coords,ocoords
real(dp), dimension(npart) :: dispx,dispy
real(dp) :: dismax,rv,lx2
logical :: update

dismax = 0.0d0
dispx = coords(:,1) - ocoords(:,1)
do i = 1,npart
  if (abs(dispx(i)) > lx2) then
    dispx(i) = dispx(i) - sign(lx2*2.0d0,dispx(i))
  endif
enddo
dispy = coords(:,2) - ocoords(:,2)
! print *, "old coords "
! write(6,'(20F10.5)') ocoords(1,:)
! write(6,'(20F10.5)') ocoords(2,:)
! print *, "current coords "
! write(6,'(20f10.5)') coords(1,:)
! write(6,'(20F10.5)') coords(2,:)
! print *, "movement since last update "
! write(6,'(20f10.5)') sqrt(disp)
! dismax = max(maxval(dispx),maxval(dispy))
! print *, "max displacement ",dismax
do i = 1,npart
  dismax = max(dispx(i)**2+dispy(i)**2,dismax)
enddo
update = (SQRT(dismax)>=(rv-1.5d0)/2.0d0)
! print *, "hit enter to continue"
! read *
return
end subroutine check

end module
