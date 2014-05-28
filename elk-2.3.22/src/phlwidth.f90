
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phlwidth
use modmain
use modphonon
implicit none
! local variables
integer i,j,iq,iv
integer lwork,info
real(8) gmin,gmax
! allocatable arrays
real(8), allocatable :: wq(:),gq(:,:),gp(:,:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: dynq(:,:,:)
complex(8), allocatable :: ev(:,:),b(:,:)
complex(8), allocatable :: gmq(:,:,:),gmr(:,:,:)
complex(8), allocatable :: gmp(:,:),work(:)
! initialise universal variables
call init0
call init2
allocate(wq(nbph),gq(nbph,nqpt),gp(nbph,npp1d))
allocate(rwork(3*nbph))
allocate(dynq(nbph,nbph,nqpt))
allocate(ev(nbph,nbph),b(nbph,nbph))
allocate(gmq(nbph,nbph,nqpt))
allocate(gmr(nbph,nbph,nqptnr))
allocate(gmp(nbph,nbph))
lwork=2*nbph
allocate(work(lwork))
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! read in the phonon linewidths for each q-point
call readgamma(gq)
! loop over phonon q-points
do iq=1,nqpt
! diagonalise the dynamical matrix
  call dynev(dynq(:,:,iq),wq,ev)
! construct a complex matrix from the phonon eigenvectors such that its
! eigenvalues squared are the phonon linewidths
  do i=1,nbph
    do j=1,nbph
      b(i,j)=sqrt(abs(gq(i,iq)))*conjg(ev(j,i))
    end do
  end do
  call zgemm('N','N',nbph,nbph,nbph,zone,ev,nbph,b,nbph,zzero,gmq(:,:,iq),nbph)
end do
! Fourier transform the gamma matrices to real-space
call dynqtor(gmq,gmr)
! generate a set of q-point vectors along a path in the Brillouin zone
call connect(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
gmin=1.d8
gmax=0.d0
! compute the linewidths along the path
do iq=1,npp1d
! compute the gamma matrix at this particular q-point
  call dynrtoq(vplp1d(:,iq),gmr,gmp)
! diagonalise the gamma matrix
  call zheev('N','U',nbph,gmp,nbph,gp(:,iq),work,lwork,rwork,info)
! square the eigenvalues to recover the linewidths
  gp(:,iq)=gp(:,iq)**2
  gmin=min(gmin,gp(1,iq))
  gmax=max(gmax,gp(nbph,iq))
end do
gmax=gmax+(gmax-gmin)*0.5d0
gmin=gmin-(gmax-gmin)*0.5d0
! output the vertex location lines
open(50,file='PHLWLINES.OUT',action='WRITE',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),gmin
  write(50,'(2G18.10)') dvp1d(iv),gmax
  write(50,'("     ")')
end do
close(50)
! output the phonon linewidth dispersion
open(50,file='PHLWIDTH.OUT',action='WRITE',form='FORMATTED')
do i=1,nbph
  do iq=1,npp1d
    write(50,'(2G18.10)') dpp1d(iq),gp(i,iq)
  end do
  write(50,'("     ")')
end do
close(50)
write(*,*)
write(*,'("Info(phlwidth):")')
write(*,'(" phonon linewidth dispersion written to PHLWIDTH.OUT")')
write(*,'(" vertex location lines written to PHLWLINES.OUT")')
deallocate(wq,gq,gp,rwork,dynq)
deallocate(ev,b,gmq,gmr,gmp,work)
return
end subroutine

