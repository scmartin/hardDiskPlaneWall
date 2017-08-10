!*******************************************************
!
!  ToDo: Add ability to take in input file
!
!*******************************************************

program diskMC

use overlap
use ellipse_init
use verletmod

implicit none

!****************************************************************************
!
!  
!
!
!
!****************************************************************************

INTEGER, PARAMETER  ::  sp = SELECTED_REAL_KIND(6,37)
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15,307)
real(dp), dimension(:,:), allocatable :: coords,ocoords
real(dp), dimension(3) :: r
integer(dp), dimension(:), allocatable :: ndens,instdens,chunkdens,vlist,point
real(dp), dimension(:), allocatable :: rhox,rhoblock
integer(dp) :: icycle,i,j,nbin,ibin,nsamp,ibinmin,ibinmax,listmax,chkcount=0, &
           vlcount=0
real(dp) :: movemax,lx,ly,dx,dy,xold,yold,rv,vxsblock=0.0
real(dp) :: bin,rho,lx2,binvol,rhoave,tstart,tfin,timer,tottime=0.0d0, &
            vstart,vfin,vtimer,vtime=0.0d0,cstart,cfin,ctimer,ctime=0.0d0, &
            vtimemax=0.0d0
integer(dp) :: n,iaccept,nupdate,npart,col,   &
             ncycle,neq,neqstep,nsteps,istep
character(24) :: filename,ndprof,incoords,coordfile,errfile, &
                 infile,formstring,input,densfile,volfile
character(10) :: time
real(dp) :: v,vxs,nface
logical :: update

integer(dp) :: values(8)

! give a random seed using the time

call date_and_time(time=time,values=values)
call random_seed(put=values(7:8))

!*****************************************************************************
!
!  Asks user for input file. If no input file is used starting coordinates
!  are generated and saved to a file. The starting coordinates generated are
!  random and therefore this method does not work well at high densities.
!
!*****************************************************************************
CALL get_environment_variable("PBS_JOBID",filename)
write(*,*) "Do you have an input file? y or n"
read(*,*) input

!  not yet able to read in an input file

if (input == "y") then
  write(*,*) "input file name?"
  read(*,*) infile
  write(*,*) infile
  open(20,file=infile)
  read(20,"(1I10)") ncycle ! # of production cycles
  read(20,"(1I10)") neq    ! # of equillibrium cycles
  read(20,"(1I10)") nbin   ! # of density histogram bins
  read(20,"(1I10)") npart
  read(20,"(F25.0)") ly   ! length of the hard wall sides of the sim box
  read(20,"(f25.0)") lx   ! length of the open sides of the sim box
  read(20,"(f25.0)") movemax  ! max MC step size
  read(20,"(1I10)") ibinmin  ! first histogram bin of the bulk region
  read(20,"(1F15.0)") rv  ! verlet neighbor list radius
  write(6,"(A15,1I10)") 'prod_cycles',ncycle ! # of production cycles
  write(6,"(A15,1I10)") 'eq_cycles',neq    ! # of equillibrium cycles
  write(6,"(A15,1I10)") 'bins',nbin   ! # of density histogram bins
  write(6,"(A15,1I10)") 'particles',npart
  write(6,"(A15,F10.2)") 'y',ly   ! length of the hard wall sides of the sim box
  write(6,"(A15,F10.2)") 'x',lx   ! length of the open sides of the sim box
  write(6,"(A15,F10.2)") 'move_max',movemax  ! max MC step size
  write(6,"(A15,1I10)") 'bulkbin',ibinmin  ! first histogram bin of the bulk region
  write(6,"(A15,1F10.2)") 'verlet_radius',rv  ! verlet neighbor list radius
  update = .true.     ! update==true means that the neighbor lists need updated
                      ! and triggers a call to the appropriate subroutine
!  npart = nint(rho*ly*lx)  ! # of particles in the simulation cell
  rho = dble(npart)/(lx*ly)
  ibinmax = nbin-ibinmin+1  ! mirror of ibinmin
  listmax = int(rho*npart*4*rv**2)  ! sets the maximum size of the neighbor list
  allocate(coords(npart,2))      ! current particle coordinates
  allocate(ocoords(npart,2))     ! particle coordinates at the time of the last
                                 !   verlet list udate
  allocate(ndens(nbin))          ! counts the number of particles counted in a 
                                 !   a given histogram bin over the whole
                                 !   simulation
  allocate(instdens(nbin))       ! instantaneous number of particles in a bin
                                 !   as calculated at the end of each cycle.
  allocate(rhoblock(nbin))       ! average density in each bin over a small
                                 !   number of cycles, used for block averaging
  allocate(chunkdens(nbin))      ! holds the particle count over ten-cycle
                                 !   blocks, which are used for the density
                                 !   which is output for block averaging
  allocate(vlist(listmax))       ! verlet list
  allocate(point(npart+1))       ! point list holds the beginning of each
                                 !   particles neighbor list so that the
                                 !   neighbor lists can be in a single array
  vlist = 0
  point = 0
  ! write(6,*) 'npart = ',npart
  lx2 = lx*0.5d0
  bin = ly/dfloat(nbin)   ! height of bin (y-direction)
  binvol = lx*bin         ! volume of each bin
  nsteps = ncycle*npart   ! # of production steps
  
  
  nupdate = 0     ! tracks the number of times the verlet list is updated
!  neqstep = neq*npart  deprecated; number of equillibrium particle moves 
  iaccept = 0  ! number of accepted particle moves, used in optimizting
               !acceptance ratio
  nsamp = 0
  
!  print *, "nparts = ",npart
  ndens = 0
  instdens = 0
  chunkdens = 0
!  write(*,*) "initial coordinate file?"
  read(20,"(A24)") incoords
  close(20)
  open(60,file=incoords)
  do i = 1,npart
    read(60,100) coords(i,1),coords(i,2)
  enddo
  close(60)
!  stop
else
  if (input == "n") then
    write(*,*) "building coordinates and writing to file"
  else
    write(*,*) "not a valid option. building coordinates"
  endif
  write(6,*) "enter: prod. cycles, eq. cycles, nbin, rho, ly, lx, movemax, ibinmin, verlet radius"
  read(5,*) ncycle,neq,nbin,rho,ly,lx,movemax,ibinmin,rv
  
 
  update = .true.     ! update==true means that the neighbor lists need updated
                      ! and triggers a call to the appropriate subroutine
  npart = nint(rho*ly*lx)  ! # of particles in the simulation cell
  ibinmax = nbin-ibinmin+1  ! mirror of ibinmin
  listmax = int(rho*npart*4*rv**2)  ! sets the maximum size of the neighbor list

  allocate(coords(npart,2))      ! current particle coordinates
  allocate(ocoords(npart,2))     ! particle coordinates at the time of the last
                                 !   verlet list udate
  allocate(ndens(nbin))          ! counts the number of particles counted in a 
                                 !   a given histogram bin over the whole
                                 !   simulation
  allocate(instdens(nbin))       ! instantaneous number of particles in a bin
                                 !   as calculated at the end of each cycle.
  allocate(rhoblock(nbin))       ! average density in each bin over a small
                                 !   number of cycles, used for block averaging
  allocate(chunkdens(nbin))      ! holds the particle count over ten-cycle
                                 !   blocks, which are used for the density
                                 !   which is output for block averaging
  allocate(vlist(listmax))       ! verlet list
  allocate(point(npart+1))       ! point list holds the beginning of each
                                 !   particles neighbor list so that the
                                 !   neighbor lists can be in a single array
  vlist = 0
  point = 0
  ! write(6,*) 'npart = ',npart
  lx2 = lx*0.5d0
  bin = ly/dfloat(nbin)   ! height of bin (y-direction)
  binvol = lx*bin         ! volume of each bin
  nsteps = ncycle*npart   ! # of production steps
  
  write(6,"(A15,1I10)") 'prod_cycles',ncycle ! # of production cycles
  write(6,"(A15,1I10)") 'eq_cycles',neq    ! # of equillibrium cycles
  write(6,"(A15,1I10)") 'bins',nbin   ! # of density histogram bins
  write(6,"(A15,1I10)") 'particles',npart
  write(6,"(A15,F10.2)") 'y',ly   ! length of the hard wall sides of the sim box
  write(6,"(A15,F10.2)") 'x',lx   ! length of the open sides of the sim box
  write(6,"(A15,F10.2)") 'move_max',movemax  ! max MC step size
  write(6,"(A15,1I10)") 'bulkbin',ibinmin  ! first histogram bin of the bulk region
  write(6,"(A15,1F10.2)") 'verlet_radius',rv  ! verlet neighbor list radius
  
  nupdate = 0     ! tracks the number of times the verlet list is updated
!  neqstep = neq*npart  deprecated; number of equillibrium particle moves 
  iaccept = 0  ! number of accepted particle moves, used in optimizting
               !acceptance ratio
  nsamp = 0
  
  ndens = 0
  instdens = 0
  chunkdens = 0
  
  !   
  !  Generate starting coords
  !
  call coordbuild(coords,lx,ly,rho,npart)
  
  !  Print starting coords to file
  
  incoords = trim(filename)//'in.txt'
  open(20, file=incoords)
  do col = 1,npart
    write(20,100) coords(col,:)
  enddo
  close(20)
  
  ! write(6,'(50F8.3)') coords(:,1)
  ! write(6,'(50F8.3)') coords(:,2)
endif
!
!  begin monte carlo moves
!  neqsteps of equilibration are followed by nsteps of production to
!  collect averages
!  number of steps per cycle = npart
!
densfile = trim(filename)//"dens.dat"
volfile = trim(filename)//"vol.dat"
open(40, file=densfile)
open(80, file=volfile)
update = .true.
call verlet(update,rv,npart,coords,ocoords,vlist,point,lx,listmax)

print *, 'starting MC steps'
do icycle = 1,ncycle + neq
  call cpu_time(cstart)
  call check(npart,coords,ocoords,rv,update,lx2)
  call cpu_time(cfin)
  ctimer = cfin-cstart
  ctime = ctime + ctimer
  chkcount = chkcount + 1
  call cpu_time(vstart)
  if (update) then
    vlcount = vlcount + 1
    call verlet(update,rv,npart,coords,ocoords,vlist,point,lx,listmax)
  endif
  call cpu_time(vfin)
  vtimer = vfin-vstart
  vtimemax = max(vtimemax,vtimer)
  vtime = vtime + vtimer
  do istep = 1,npart
  !  if(mod(istep,(nsteps+neqstep)/10).eq.0) then
  !    write(6,*) (istep-neqstep)/npart
  !    write (*,*) "listmax = ",listmax
  !    write(6,'(1000020f10.5)') coords(:,1)
  !    write(6,'(10000F10.5)') coords(:,2)
  !  endif
  ! ^ writes the number of production cycles ran
  !
    call random_number(r)
  ! 
  !       Choose random particle 
  !
    n = int(npart*r(1)) + 1  ! n=index of particle to move
    xold = coords(n,1)
    yold = coords(n,2)
  !   print *, "particle ",n,"has moved"!*********************************************
  !   print *, "old coords are ",xold,yold !*********************************************
  !
  !       move random particle by random displacement
  !
    dx = movemax*(r(2) - 0.5d0)
    dy = movemax*(r(3) - 0.5d0)
    coords(n,1) = coords(n,1) + dx
    coords(n,2) = coords(n,2) + dy
    if(abs(coords(n,1)).gt.lx2) then
      coords(n,1) = coords(n,1) - sign(lx,coords(n,1))                   ! minimum image convention
  !    print *, "old coords are ",xold,yold
  !    print *, "new coords are ",coords(n,1),coords(n,2)
  !    read *
    endif
  !   print *, "new coords are ",coords(n,1),coords(n,2)
  !   print *, "step # ",istep
  !
  !     reject move if new particle position overlaps with wall or another 
  !     particle, otherwise, accept
  !
  !   print *, "preparing to test for impact. hit enter"
  !   print *, "coords are"
  !   write(6,'(50F8.3)') coords(:,1)
  !   write(6,'(50F8.3)') coords(:,2)
  !   read *
    call cpu_time(tstart)
    if((coords(n,2)<0.5d0).or.(coords(n,2)>(ly-0.5d0))  &
       .or.(vimpact(n) == .true.)) then
  !     write(6,*) " REJECTED! particle ",n
      coords(n,1) = xold
      coords(n,2) = yold
    else
      if(icycle>neq) iaccept = iaccept + 1
    endif
    call cpu_time(tfin)
    timer = tfin-tstart
    tottime = tottime + timer
  !
  ! update density profile: ndens is a running total of the number of particles in
  ! each bin over all the production cycles. The average occurs in the loop 
  ! below "output density profile"
  !
  enddo
  if(icycle.gt.neq) then  ! start after equillibration 
    do i = 1,npart
      ibin = int(coords(i,2)/bin)+ 1          ! sets which bin particle i is in
      ndens(ibin) = ndens(ibin)+1             ! this loops creates running total
      instdens(ibin) = instdens(ibin)+1       ! of particles in each bin
      chunkdens(ibin) = chunkdens(ibin)+1
    enddo
    vxs = (lx*2.0d0*bin*dble(ibinmin-1) &
          - lx*bin*dble(ibinmax-ibinmin+1)*(sum(instdens(:ibinmin-1))+sum(instdens(ibinmax+1:)))/sum(instdens(ibinmin:ibinmax))) &
          /(lx*2.0d0)
    vxsblock = vxsblock + vxs
    vxs = 0.0d0
    nsamp = nsamp + 1  ! nsamp is number of production cycles so far
!    write(*,*) nsamp
    if(mod(icycle,10)==0) then  ! this block outputs averages over 10 cycles
      do ibin = 1,nbin          ! which can then be block averaged
        rhoblock(ibin) = (dfloat(chunkdens(ibin))/10.0d0)/binvol 
      enddo
      write(40,300) rhoblock
      write(80,"(F18.12)") vxsblock/10.0d0
      chunkdens = 0
      vxsblock = 0.0d0
      rhoblock = 0.0d0
    endif
    instdens = 0
  endif
enddo
close(40)
close(80)
deallocate(chunkdens)
deallocate(rhoblock)
deallocate(instdens)
deallocate(vlist)
deallocate(point)
deallocate(ocoords)
!  output density profile

allocate(rhox(nbin))
ndprof = trim(filename)//'dens.txt'
open(22, file=ndprof)
do ibin  = 1,nbin
  rhox(ibin) = (dble(ndens(ibin))/dble(nsamp))/binvol  ! rhox(ibin)=average density of each bin
  write(22,100) (ibin-0.5d0)*bin,rhox(ibin)
!  write(4,*) (ibin-0.5d0)*bin,ndens(ibin)
enddo
! print *, sum(rhox*binvol)
close(22)
!  calulate bulk density by averaging over the middle of the simulation box

rhoave = 0.0d0
do ibin = ibinmin,ibinmax
  rhoave = rhoave + rhox(ibin)
enddo
rhoave = rhoave/dfloat(ibinmax-ibinmin+1)
v = binvol*(dble(ibinmin-1)*2.0d0)
nface = 0.0d0

! do ibin = 1,nbin
!   if(ibin<ibinmin.or.ibin>ibinmax) then
! !    write(*,*) ibin,ndens(ibin),dble(ndens(ibin))/dble(nsamp)
!     nface = nface + dble(ndens(ibin))/dble(nsamp)
!   endif
! enddo

! nface below is the average density of the interfacial region
nface = (sum(rhox(:ibinmin-1))+sum(rhox(ibinmax+1:)))/(dble(ibinmin-1)*2.0d0)

!vxs = (v - nface/rhoave)/(2.0d0*lx)

! vxs is the excess volume calculated using the interface and bulk densities
! averaged over the entire simulation
vxs = dble(ibinmin-1)*bin - (nface*dble(ibinmin-1)*bin/rhoave)

coordfile = trim(filename)//'out.txt'
open(21, file=coordfile)

do col = 1,npart
  write(21,100) coords(col,:)
enddo

close(21)

write(6,*) '# bulk density = ',rhoave
write(6,*) '# v_excess = ',vxs
write(6,*) "# impact time = ",tottime
write(6,*) "# max vlist time = ", vtimemax
write(6,*) "# list update time = ",vtime
write(6,*) "# check time = ",ctime
write(6,*) "# of checks = ",chkcount
write(6,*) "# of vlist updates = ",vlcount

!
!     output acceptance ratio
!
write(6,*) "# acceptance ratio = ",dble(iaccept)/dble(nsteps)
100  format(3F18.7)
200  format(3A18)
300  format(10000F20.16)
deallocate(coords)
deallocate(ndens)
deallocate(rhox)

CONTAINS
!**********************************************************************
!  old function for testing an MC move for rejection. This function
!  was used before the verlet list was implemented
!
!**********************************************************************

!     integer function impact(coords,n,npart,point,vlist)  ! n=index of particle that was moved 
!       
!     INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15,307)
!     integer :: npart,n
!     real(dp) :: coords(npart,2),r2,delx
!     integer :: i
! 
!     impact = 0
! !   write(6,*) "got to overlap ",n,npart
!     do i = 1,npart
!       if(i.eq.n) go to 10
!       delx = coords(i,1) - coords(n,1)
!       if(delx.gt.0.5d0*lx) then
!         delx = delx - lx
!       endif 
!       r2 = delx**2 + (coords(n,2)-coords(i,2))**2
!       if(r2.lt.1.0d0) then
!         impact = 1
! !       write(6,*) "overlap ",r2,n,i
!         return
!       endif
!  10   continue
!     enddo
!     return
!     end function impact

!*************************************************************************
!
! This function takes in a particle just moved and checks it against its
! neighbor list to look for impacts. If an impact is found, vimpact=.true.
! is returned and the move is rejected
!
!
!*************************************************************************

  logical function vimpact(n)  ! n=index of particle that was moved 
    
  INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15,307)
  integer(dp),intent(in) :: n
  real(dp) :: r2,delx
  integer(dp) :: i,j,start,fin
!  logical(dp) :: test


  vimpact = .false.
  start = point(n)
  fin = point(n+1)-1
  if (start <= fin) then
    do i = start,fin
      j = vlist(i)
      delx = coords(j,1) - coords(n,1)
!      test = .false.
      if(abs(delx).gt.0.5d0*lx) then
!        test = .true.
        delx = delx - sign(lx,delx)
      endif 
      r2 = delx**2 + (coords(j,2)-coords(n,2))**2
!      if (test) then
!        print *, "r^2 = ",r2
!        read *
!      endif
      if(r2.lt.1.0d0) then
        vimpact = .true.
!        if (test) then
!          write(6,*) "particle ",n," hit part ",j
!          read *
!        endif
        return
      endif
    enddo
  else
!     print *, "no neighbors"
  endif
  return
  end function vimpact


end program diskmc
