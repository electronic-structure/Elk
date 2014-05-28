
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rotzflm
! !INTERFACE:
subroutine rotzflm(rot,lmin,lmax,n,ld,zflm1,zflm2)
! !INPUT/OUTPUT PARAMETERS:
!   rot   : rotation matrix (in,real(3,3))
!   lmin  : minimum angular momentum (in,integer)
!   lmax  : maximum angular momentum (in,integer)
!   n     : number of functions to rotate (in,integer)
!   ld    : leading dimension (in,integer)
!   zflm1 : coefficients of the complex spherical harmonic expansion for each
!           function (in,complex(ld,n))
!   zflm2 : coefficients of rotated functions (out,complex(ld,n))
! !DESCRIPTION:
!   Rotates a set of complex functions
!   $$ f_i({\bf r})=\sum_{lm}f_{lm}^iY_{lm}(\hat{\bf r}) $$
!   for all $i$, given the coefficients $f_{lm}^i$ and a rotation matrix $R$.
!   This is done by first the computing the Euler angles $(\alpha,\beta,\gamma)$
!   of $R^{-1}$ (see routine {\tt roteuler}) and then applying the spherical
!   harmonic rotation matrix generated by the routine {\tt ylmrot}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Modified, December 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rot(3,3)
integer, intent(in) :: lmin,lmax
integer, intent(in) :: n
integer, intent(in) :: ld
complex(8), intent(in) :: zflm1(ld,n)
complex(8), intent(out) :: zflm2(ld,n)
! local variables
integer lmmax,l,lm1,lm2,nm,p
real(8) det,rotp(3,3)
real(8) ang(3),angi(3)
complex(8), parameter :: zzero=(0.d0,0.d0),zone=(1.d0,0.d0)
! allocatable arrays
complex(8), allocatable :: d(:,:)
if (lmin.lt.0) then
  write(*,*)
  write(*,'("Error(rotzflm): lmin < 0 : ",I8)') lmin
  write(*,*)
  stop
end if
if (lmin.gt.lmax) then
  write(*,*)
  write(*,'("Error(rotzflm): lmin > lmax : ",2I8)') lmin,lmax
  write(*,*)
  stop
end if
if (n.eq.0) return
if (n.lt.0) then
  write(*,*)
  write(*,'("Error(rotzflm): n < 0 : ",I8)') n
  write(*,*)
  stop
end if
lmmax=(lmax+1)**2
allocate(d(lmmax,lmmax))
! find the determinant
det=rot(1,2)*rot(2,3)*rot(3,1)-rot(1,3)*rot(2,2)*rot(3,1) &
   +rot(1,3)*rot(2,1)*rot(3,2)-rot(1,1)*rot(2,3)*rot(3,2) &
   +rot(1,1)*rot(2,2)*rot(3,3)-rot(1,2)*rot(2,1)*rot(3,3)
! make the rotation proper
p=1
if (det.lt.0.d0) p=-1
rotp(:,:)=dble(p)*rot(:,:)
! compute the Euler angles of the rotation matrix
call roteuler(rotp,ang)
! inverse rotation: the function is to be rotated, not the spherical harmonics
angi(1)=-ang(3)
angi(2)=-ang(2)
angi(3)=-ang(1)
! determine the rotation matrix for complex spherical harmonics
call ylmrot(p,angi,lmax,lmmax,d)
! apply rotation matrix (d and zflm may have different starting indices)
lm1=lmin**2+1
lm2=1
do l=lmin,lmax
  nm=2*l+1
  call zgemm('N','N',nm,n,nm,zone,d(lm1,lm1),lmmax,zflm1(lm2,1),ld,zzero, &
   zflm2(lm2,1),ld)
  lm1=lm1+nm
  lm2=lm2+nm
end do
deallocate(d)
return
end subroutine
!EOC

