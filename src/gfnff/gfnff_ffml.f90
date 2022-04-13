! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Module  for FF calculation with ML correction
module xtb_gfnff_ffml
  use xtb_mctc_accuracy, only : wp
  use xtb_gfnff_topology, only : TGFFTopology 
  use xtb_type_molecule                        
  use xtb_type_environment
  use xtb_type_restart, only : TRestart
  use xtb_io_reader, only : readMolecule
  use xtb_io_writer, only : writeMolecule
  use mctc_io_filetype, only : filetype, getFileType => get_filetype

  implicit none
  private

  public :: Tffml, calc_ML_correction, set_ffml

  ! holds info for ML correction of GFN-FF calculation
  type :: Tffml
    logical  :: fixMD
    !logical  :: runML
! contains  ! procedures...
  end type Tffml

contains

! setup type for FF calculation with ML correction
subroutine set_ffml(self)
  class(Tffml), intent(out) :: self
  ! run deterministic MD simulation (0K temperature and fixed shifts)
  self%fixMD = .true.
!  self%runML = .true.
end subroutine set_ffml

! the actual routine for the ML correction
  !@thomas check if all input is needed
subroutine calc_ML_correction(ffml,fname,topo, mol)
        type(Tffml), intent(in)        :: ffml
  character(len=*), intent(in)    :: fname
  type(TGFFTopology), intent(in)  :: topo
  type(TMolecule), intent(in)     :: mol

  character(len=:), allocatable   :: cmd  !@thomas 
  integer :: ich, sdf_ftype, maxoptcycle, scciter
  type(TEnvironment) :: envtmp
  type(TMolecule)    :: moltmp
  type(TRestart) :: chktmp
  real(wp) :: egap, etemp, etot, stmp(3,3)
  real(wp), allocatable :: gtmp(:,:)
  ! trivial approach: just calling the gen_ref_struc.sh script
if (.false.) then !<delete !@thomas 
  cmd="~/bin/gen_ref_struc_orig.sh "//fname
!  write(*,*) 'cmd |>',cmd,'<|'
  call execute_command_line(cmd)
else !<delete
  ! hopefully more suffisitcated approach

  !(1)! convert 3D to SMILE 
  write(*,*) 'Coordinate file: ',fname
  write(*,*) 'Generating SMILE'
  cmd="obabel "//fname//" -O out.smi -b"
  !write(*,*) 'cmd |>',cmd,'<|'
  call execute_command_line(cmd)

  !(2)! remove cis/trans info from SMILE (represented as / \ )
  !@thomas TODO not good to use external script here?!
  cmd="~/bin/rmCisTrans.sh"
  call execute_command_line(cmd)

  !(3)! convert SMILE to 2D   geo_2D.sdf
  ! file names are hardcoded from here on
  write(*,*) 'Convert SMILE to 2D'
  cmd="obabel out.smi -O geo_2D.sdf --gen2d -h"
  call execute_command_line(cmd)

  !(4)! load geo_2D.sdf geometry
  ! need env, ftype for sdf is 6 (see comment below), 
  !sdf_ftype = getFileType('geo_2D.sdf')
  !write(*,*) 'sdf_ftype =',sdf_ftype  ! gives sdf_ftype = 6
  sdf_ftype = 6
  call init(envtmp)
  call open_file(ich, 'geo_2D.sdf', 'r')
  call readMolecule(envtmp, moltmp, ich, sdf_ftype)
  call close_file(ich)

  !(5)! run structure converter (2D to 3D) with deterministic (fixed) MD
  ! retrieved Data: maxscciter=300.0_wp, maxoptcycle=0,
  ! ml_struc_converter also optimizes structure:
  !(6)! optimize reference structure
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! IMPORTANT NOTE !!! IF OPTIMIZATION FAILS: gfnff_topo file might be wrong !
  !  maybe reference topo is too different from original topo
  egap=0.0_wp
  etemp=0.0_wp
  scciter=250
  maxoptcycle=0
  etot=0.0_wp
  stmp=0.0_wp
  allocate(gtmp(3,moltmp%n), source=0.0_wp)
  call ml_struc_convert(envtmp, .false., moltmp, chktmp, egap, &
          etemp, scciter, maxoptcycle, etot, gtmp, stmp)

  !(7)! convert to original file format  ! use mctc file reader and writer
  call open_file(ich, 'ref_struc.xyz', 'w')
  call writeMolecule(moltmp, ich, format=fileType%xyz, energy=etot, &
   & gnorm=NORM2(gtmp))                                                   
  call close_file(ich)
  write(*,*) 'Optimized reference structue.'
  write(*,*) 'Geometry is written to ref_struc.xyz'
  write(*,*)
  write(*,*) '__etot= ', etot
  write(*,*) '__gnorm= ', NORM2(gtmp)
endif !<delete


  write(*,*)
end subroutine calc_ML_correction


subroutine ml_struc_convert( &
         & env,restart,mol,chk,egap,et,maxiter,maxcycle,&
         & etot,g,sigma)
  use xtb_mctc_accuracy, only : wp
  use xtb_gfnff_param
  use xtb_gfnff_setup
  use xtb_disp_dftd3param
  use xtb_type_environment
  use xtb_type_molecule
  use xtb_type_restart
  use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
  use xtb_type_data
  use xtb_restart
  use xtb_setmod
  use xtb_setparam
  use xtb_dynamic
  use xtb_geoopt
  use xtb_readin, only : xfind
  implicit none
! Dummy -----------------------------------------------------------------------
  type(TEnvironment),intent(inout)            :: env
  type(TMolecule),intent(inout)               :: mol
  type(TRestart),intent(inout)                :: chk
  integer,intent(in)                          :: maxiter
  integer,intent(in)                          :: maxcycle
  real(wp),intent(inout)                      :: etot
  real(wp),intent(in)                         :: et
  real(wp),intent(inout)                      :: egap
  real(wp),intent(inout)                      :: g(3,mol%n)
  real(wp),intent(inout)                      :: sigma(3,3)
  logical,intent(in)                          :: restart
  character(len=:),allocatable                :: fnv
! Stack -----------------------------------------------------------------------
  type(TEnvironment)                          :: env2
  type(TGFFCalculator)                        :: calc, calc2
  integer                                     :: ich
  integer                                     :: idum
  integer                                     :: mode_input
  real(wp)                                    :: time_in
  real(wp)                                    :: temp_in
  real(wp)                                    :: step_in
  real(wp)                                    :: dump_in
  real(wp)                                    :: hmass_in
  logical                                     :: exist
  logical                                     :: fail
  integer, allocatable :: opt_in
  character(len=*),parameter                  :: p_fname_param_gfnff = '.param_gfnff.xtb'
! loop geoopt -----------------------------------------------------------------
  integer                                     :: i, j
  integer, parameter                          :: num_shift_runs = 3
  real(wp), allocatable                       :: shift(:)
  real(wp), allocatable                       :: mol_xyz_arr(:,:,:)
  real(wp)                                    :: etot_arr(num_shift_runs)
  real(wp)                                    :: g_arr(num_shift_runs)
  real(wp)                                    :: sign_threshold
  type(TMolecule)                             :: mol_shifted
  allocate(shift(mol%n), source=0.0_wp)
  allocate(mol_xyz_arr(3, mol%n, num_shift_runs), source=0.0_wp)
  etot_arr = 0.0_wp
!------------------------------------------------------------------------------
! set up force field
  call struc_convert_header(env%unit)
  if (allocated(set%opt_engine)) then
    opt_in = set%opt_engine
  end if
  set%opt_engine = p_engine_rf
  mode_input = set%mode_extrun
  set%mode_extrun = p_ext_gfnff
  if (.not.allocated(fnv)) fnv=xfind(p_fname_param_gfnff)
  call newGFFCalculator(env, mol, calc, fnv, restart, gffVersion%harmonic2020)

!===============================
! Set Block
  time_in  = set%time_md
  set%time_md  = 5.0_wp              ! short 5 ps MD to move it in 3D
  temp_in  = set%temp_md
  set%temp_md  = 298.0_wp            ! md temperature 298 K
  step_in  = set%tstep_md
  set%tstep_md = 2.5_wp              ! md time step 2.5 fs
  dump_in  = set%dump_md2
  set%dump_md2 = 100.0_wp            ! md dump 100 fs
  hmass_in = set%md_hmass
  set%md_hmass = 4.0_wp              ! md hydrogen mass
  set%nvt_md = .true.                ! md thermostat
!===============================
  if (allocated(set%opt_logfile)) then
    fnv = set%opt_logfile
  else
    deallocate(fnv)
  endif
  set%opt_logfile = 'convert.log'
!------------------------------------------------------------------------------
! force field geometry optimization
  ! loop runs 3 geoopt with different shifts in the new 3rd coordinate
  ! and then keeps the mol%xyz with lowest etot for md
  mol_shifted = mol
  do i=1, num_shift_runs
    mol_shifted = mol
    ! create array with deterministic alternating shifts between -1 and 1 using sin()
    do j=1, size(shift)
      shift(j) = ((-1)**j)*sin(0.1*real(i)*real(j))**2
    enddo
    ! set temperature to 0K
    set%temp_md  = 0.0_wp    ! md temperature 0 K
    mol_shifted%xyz(3,:) = mol_shifted%xyz(3,:) + shift  ! apply shifts
    
    call geometry_optimization &
        &     (env,mol_shifted,chk,calc,   &
        &      egap,set%etemp,maxiter,maxcycle,etot,g,sigma,p_olev_crude,.false.,.true.,fail)
    mol_xyz_arr(:,:,i) = mol_shifted%xyz  ! store optimized xyz
    etot_arr(i) = etot                    ! store energy etot
    g_arr(i) = NORM2(g)                   ! store gradient norm
  enddo

  ! always take geometry from first optimization unless unreasonable !@thomas TODO maybe take lowest e again??
  mol%xyz(:,:) = mol_xyz_arr(:,:,minloc(etot_arr, DIM=1))  ! keep xyz with lowest etot
  ! TODO TODO TODO
  ! if(g_arr(1).lt.50.0_wp) then
  !   mol%xyz=mol_xyz_arr(:,:,1)
  ! elseif(g_arr(2).lt.50.0_wp) then
  !   mol%xyz=mol_xyz_arr(:,:,2)
  ! elseif(g_arr(3).lt.50.0_wp) then
  !   mol%xyz=mol_xyz_arr(:,:,3)
  ! else
  !   call env%error('Structure converter could not optimize geometry properly.')
  ! endif

    if (allocated(fnv)) then
      set%opt_logfile = fnv
    else
      deallocate(set%opt_logfile)
    endif
    write(*,*)
!------------------------------------------------------------------------------
! force field md simulation
  idum = 0
  call md                &
      &   (env,mol,chk,calc, &
      &    egap,set%etemp,maxiter,etot,g,sigma,0,set%temp_md,idum)
!------------------------------------------------------------------------------
! set all back to input
  set%time_md  = time_in
  set%temp_md  = temp_in
  set%tstep_md = step_in
  set%dump_md2 = dump_in
  set%md_hmass = hmass_in
  if (allocated(opt_in)) then
    set%opt_engine = opt_in
  else
    deallocate(set%opt_engine)
  end if
  ! optimize the final geometry
  etot=0.0_wp
  g=0.0_wp
  sigma=0.0_wp
  if (.not.allocated(fnv)) fnv=xfind(p_fname_param_gfnff)
  if (allocated(fnv)) write(*,*) 'test fnv allocation delete this'
  call newGFFCalculator(env, mol, calc2, fnv, .false.)
  call geometry_optimization &
      &     (env,mol,chk,calc2,   &
      &      egap,set%etemp,maxiter,maxcycle,etot,g,sigma,p_olev_crude,.false.,.true.,fail)
!------------------------------------------------------------------------------
  write(*,*)
  write(*,'(10x," ------------------------------------------------- ")')
  write(*,'(10x,"|           2D => 3D conversion done!             |")')
  write(*,'(10x," ------------------------------------------------- ")')
  write(*,*)
  set%mode_extrun = mode_input
  mol%info%two_dimensional=.false.
  call gfnff_param_dealloc(calc%topo)

end subroutine ml_struc_convert

end module xtb_gfnff_ffml
