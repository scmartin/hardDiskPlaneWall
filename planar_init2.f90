!*********************************************************************************
!
! This module creates the initial coordinates for the hard disk simulation. It
! starts by randomly trying to insert particles in the simulation box. After a 
! critical ratio of successes/attempts is reached, particles are shifted towards
! the hard walls and new insertion attempts are limited to the middle region of
! the simulation box. This process is repeated until all particles are inserted.
! This process allows number densities of at least 0.9 to be used without a 
! crystal lattice input file. For a system of 10000 particles at a density of
! 0.9, this process takes less than five minutes. Some optimization has been
! done, but further optimization could yield faster initialization.
!
!
!*********************************************************************************

module ellipse_init

use overlap  ! overlap module checks 

CONTAINS

subroutine coordbuild(coords,lx,ly,rho,npart)

IMPLICIT NONE
INTEGER, PARAMETER  ::  sp = SELECTED_REAL_KIND(6,37)
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15,307)
real (dp) :: coords(npart,2)
real (dp) :: r(3),r2
real (dp) :: rho,lx,ly,x,y,xx,yy,delx
integer(dp) :: i,j,k,pcount,npart,attempts=10,accepts=10,resets=0,moves=0, &
               stuff,part,misses=0,start,endd

!   
!   place first particle
!   n = 1 
10   call random_number(r)
coords(1,1) = lx*(r(1)-0.5d0)         ! insert first particle
coords(1,2) = (ly-1.0d0)*r(2) + 0.5d0
!open(100,file="fileshit.txt")
do pcount = 2,npart
20  call random_number(r)
!    print *, "attempts",attempts
!    print *, "accepts",accepts
!    print *, "resets",resets
    if (resets == 0) then  ! resets counts the # of times the accept/attempt
                           ! ratio reached the critical value.

      x = lx*(r(1)-0.5d0)         ! until the reset ratio is reached the first time
      y = (ly-1.0d0)*r(2) + 0.5d0 ! the particles are inserted anywhere within the box
    ! after the first reset, particles are inserted closer to the center of the
    ! box as previous particles are moved towards the walls    
    else
      x = lx*(r(1)-0.5d0)
      y = ((r(2)-0.5d0)*(ly-2.0d0)/dble(resets+1)) + (r(3)-0.5d0)*1.0d0 + ly*0.5d0
    endif
    attempts = attempts + 1
!    print *, "accept/attempt",dble(accepts)/dble(attempts)
    if(ovrlap(coords,x,y,lx,npart,pcount).eq.0) then
!  overlp() iterates over previous particles, tests
!  for overlap based on center and radii
        coords(pcount,1) = x 
        coords(pcount,2) = y
        accepts = accepts+1
!	print *, 'accept'
    else
      ! if the acceptance ratio for particle insertions is too low, this loop
      ! pushes previous particles towards the walls, checking for particle
      ! overlaps at each step. This is the most computationally expensive part
      ! of the subroutine
      if(dble(accepts)/dble(attempts) < 0.5) then
        misses = 0
        do i = 1,resets
          ! for the first half iterations, all particles are moved
          if (i < resets/2) then
            start = 1
          ! the second half of the loop only moves particles that have been
          ! added more recently, as earlier particles are likely more closely
          ! packed near the wall and attempting to move them is likely to fail
          else
            start = misses
          endif
          if (misses >= pcount) misses = int(misses/2)
          do j = start,pcount-1
            call random_number(r)
!            j = int(r(3)*(pcount))
            xx = coords(j,1)
            yy = coords(j,2)
            ! particles added earlier in the simulation are moved smaller
            ! distances, as they are more likely to have already been pushed
            ! towards the wall
            coords(j,1) = coords(j,1)+(r(1)-0.5d0)*dble(i)/dble(resets)
            ! particles are moved away from the center of the box in the y
            ! direction
            if (coords(j,2) >= ly/2.0d0) then
              coords(j,2) = coords(j,2)+(r(2)-0.1d0)*dble(i)/dble(resets)
            else
              coords(j,2) = coords(j,2)-(r(2)-0.1d0)*dble(i)/dble(resets)
            endif
            ! don't allow moves which wrap particles in the x direction or
            ! impact the hard walls
            if (abs(coords(j,1))>=lx/2.0d0.or.coords(j,2)>ly-0.5d0.or.coords(j,2)<0.5d0) then
              coords(j,1) = xx
              coords(j,2) = yy
              misses = misses+1
            else
              ! check for particle overlap with each move
              do k = 1,pcount-1
                if (k /= j) then
                  delx = coords(j,1)-coords(k,1)
                  if (abs(delx) >= lx*0.5d0) delx = delx - sign(lx,delx)
                  r2 = (delx)**2 + (coords(j,2)-coords(k,2))**2
                  if (r2 <= 1.0d0) then
                    coords(j,1) = xx
                    coords(j,2) = yy
                    misses = misses+1
                    exit
                  endif
                endif
                if (k == pcount) then
                  moves = moves+1
                endif
              enddo
            endif
          enddo
        enddo
      attempts = 10
      resets = resets + 1
!      print *, "resets",resets
!      print *, "accept/attempt",dble(accepts)/dble(attempts)
!      write(100,"(5000f18.7)") coords(:pcount-1,1) 
!      write(100,"(5000f18.7)") coords(:pcount-1,2)
!      write(100,*) "\n"
      endif
!      print *, 'reject'
      go to 20
    endif
enddo
!close(100)
!do pcount = 1,npart
!    print *, coords(pcount,1),coords(pcount,2)
!enddo
end subroutine coordbuild

end module ellipse_init
