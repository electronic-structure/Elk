
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genstress
use modmain
use modmpi
use modstore
implicit none
! local variables
integer i
real(8) a(3,3),et0,t1
! store original parameters
avec0(:,:)=avec(:,:)
tforce0=tforce
tforce=.false.
! restore original symmetries
call symmetry
! generate the strain tensors
call genstrain
! run the ground-state calculation
call gndstate
! store the total energy
et0=engytot
! subsequent calculations will read STATE.OUT
trdstate=.true.
! loop over strain tensors
do i=1,nstrain
  if (mp_mpi) then
    write(*,'("Info(genstress): strain tensor ",I1," of ",I1)') i,nstrain
  end if
! displace lattice vectors
  call r3mm(strain(:,:,i),avec,a)
  avec(:,:)=avec0(:,:)+deltast*a(:,:)
! run the ground-state calculation
  call gndstate
! compute the stress tensor component
  stress(i)=(engytot-et0)/deltast
end do
! compute the maximum stress magnitude over all lattice vectors
stressmax=0.d0
do i=1,nstrain
  t1=abs(stress(i))
  if (t1.gt.stressmax) stressmax=t1
end do
! restore original parameters
avec(:,:)=avec0(:,:)
tforce=tforce0
return
end subroutine

