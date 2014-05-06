subroutine sirius_ground_state
use modmain
implicit none
#ifdef SIRIUS
integer is,ia,rank
real(8) rms,eold,de

call sirius_platform_initialize(.false.)

call sirius_clear
call sirius_set_lattice_vectors(avec(:, 1), avec(:, 2), avec(:, 3))
call sirius_set_pw_cutoff(gmaxvr)
call sirius_set_aw_cutoff(rgkmax)
call sirius_set_auto_rmt(1)

do is=1,nspecies
  call sirius_add_atom_type(trim(spfname(is)), trim(spfname(is)))
enddo

do is=1,nspecies
  do ia=1,natoms(is)
    call sirius_add_atom(trim(spfname(is)), atposl(1,ia,is), bfcmt0(1,ia,is))
  enddo
enddo

call init0
call init1

call sirius_platform_mpi_rank(rank)

! init ground state class
call sirius_ground_state_initialize(kset_id)

! generate initial density
call sirius_initial_density
call sirius_density_mixer_initialize

! generate initial effective poential
call sirius_generate_effective_potential
call sirius_potential_mixer_initialize

eold = 0.d0
! main DFT loop
do iscl=1,maxscl
  ! generate potential
  call sirius_generate_effective_potential

  ! symmetrize potential
  call symrf(1,vsmt,vsir)
  if (spinpol) then
    call symrvf(1,bxcmt,bxcir)
  end if

  ! mix potential if possible
  call sirius_mix_potential
  
  ! compute new eigen states
  call sirius_find_eigen_states(kset_id, 1)

  ! find band occupancies 
  call sirius_find_band_occupancies(kset_id)

  ! generate new density
  call sirius_generate_density(kset_id)

  ! symmetrize density
  call symrf(1,rhomt,rhoir)
  if (spinpol) then
    call symrvf(1,magmt,magir)
  end if
  
  ! mix density
  call sirius_mix_density(rms)

  ! get total energy
  call sirius_get_total_energy(engytot)

  ! write some info
  call sirius_ground_state_print_info

  de = abs(engytot - eold)
  ! use epspot for rho RMS
  if (rms.lt.epspot.and.de.lt.epsengy) exit

  if (rank.eq.0) then
    write(*,'("iteration : ", I3, ", density RMS : ", G18.10, ", energy difference : ", G18.10)')iscl,rms,de
  end if

  eold = engytot
end do

call sirius_create_storage_file
call sirius_save_density
call sirius_save_potential
call sirius_print_timers
call sirius_write_json_output
call sirius_platform_barrier

call sirius_clear
return
#endif
end subroutine

