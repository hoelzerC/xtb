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
  use xtb_gfnff_topology, only : TGFFTopology, Tffml 
  use xtb_type_molecule                        
  use xtb_type_environment
  use xtb_type_restart, only : TRestart
  use xtb_io_reader, only : readMolecule
  use xtb_io_writer, only : writeMolecule
  use mctc_io_filetype, only : filetype, getFileType => get_filetype
  use forpy_mod

  implicit none
  private

  public :: calc_ML_correction

  !@thomas_ffml
  ! the Tffml type is now in topology due to circular dependencies problem

contains

! the actual routine for the ML correction
  !@thomas check if all input is needed
subroutine calc_ML_correction(env,ffml,fname,topo, mol)
  type(TEnvironment),intent(inout)            :: env
  type(Tffml), intent(in)        :: ffml
  character(len=*), intent(in)    :: fname
  type(TGFFTopology), intent(in)  :: topo
  type(TMolecule), intent(in)     :: mol

  type(TGFFTopology)              :: refTopo
  character(len=:), allocatable   :: cmd  !@thomas 
  integer :: ich, sdf_ftype, maxoptcycle, scciter, i
  type(TEnvironment) :: envtmp
  type(TMolecule)    :: moltmp
  type(TRestart) :: chktmp
  real(wp) :: egap, etemp, etot, stmp(3,3)
  real(wp), allocatable :: gtmp(:,:)
  integer :: ierror
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
          etemp, scciter, maxoptcycle, etot, gtmp, stmp, refTopo)

  !(7)! convert to original file format  ! use mctc file reader and writer
  call open_file(ich, 'ref_struc.xyz', 'w')
  call writeMolecule(moltmp, ich, format=fileType%xyz, energy=etot, &
   & gnorm=NORM2(gtmp))                                                   
  call close_file(ich)
  write(*,*) 'Optimized reference structue.'
  write(*,*) 'Geometry is written to ref_struc.xyz'
  write(*,*)
  write(*,*) '__ total energy= ', etot
  write(*,*) '__ gradient norm= ', NORM2(gtmp)
endif !<delete


  write(*,*)

  !@thomas testing interaction with python functions !@thomas_mark01
  call send_input_to_python_ML_receive_eg(env, mol%n, mol, topo, ffml)
  ! finalize forpy 
  call forpy_finalize

end subroutine calc_ML_correction


subroutine ml_struc_convert( &
         & env,restart,mol,chk,egap,et,maxiter,maxcycle,&
         & etot,g,sigma,refTopo)
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
  type(TGFFTopology), intent(out)             :: refTopo
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

  mol%xyz(:,:) = mol_xyz_arr(:,:,minloc(etot_arr, DIM=1))  ! keep xyz with lowest etot

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
  ! save reference topology
  refTopo=calc2%topo

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


subroutine send_input_to_python_ML_receive_eg(env,n, mol, topo, ffml)
  use forpy_mod
  use iso_fortran_env, only: real64
  type(TEnvironment),intent(inout)            :: env
  integer, intent(in)    :: n  ! number atoms
  type(TMolecule), intent(in)     :: mol
  type(TGFFTopology), intent(in)  :: topo
  type(Tffml), intent(in)        :: ffml
  !real(wp), intent(in)   :: xyz(3,n) ! xyz coordinates
  real(wp)   :: xyz_local(3,n)  ! xyz coordinates copy
  character(len=*), parameter :: source = "gfnff_ffml"
  integer, dimension(2) :: shape_grad
  real(wp) :: ml_energy
  real(wp) :: ml_gradient(3,n), ml_gradientT(n,3)
  real(wp) :: calc_grad(3,n) !@thomas delete debug
  integer :: i,j !@thomas delete debug
  real(wp), pointer ,dimension(:,:) :: tmp_point
  ! forpy types and needed variables
  type(module_py)  :: ml_module    ! module with all the python functions etc for ML part
  type(list)       :: paths      ! for adding path to python module
  type(object)     :: receive_obj  ! python object received from function call or similar
  type(object)     :: obj1, obj2
  type(dict)       :: kwargs, receive_dict
  type(tuple)      :: args  ! python tuple to feed into ml_function
  type(ndarray)    :: xyz_arr, tmp1, eatom_arr
  integer          :: ierror  ! return value for forpy methods (error value)

  ! initialize forpy
  ierror = forpy_initialize()
  write(*,'(a40,i5)') 'Initiallized forpy with ierror=',ierror
  
  !@thomas_mark01 goto call

  ! add path of python module file to paths
  !@thomas TODO relativer pfad!?! Geht nicht wenn jmd die xtb binary verwendet -> error handling
  ! ist halt nicht gegeben das der user das file überhaupt hat. Also unnötig das 
  ! irgendwie weiter zu suchen, ggf kann man gucken obs ne möglichkeit gibt das 
  ! zu "installieren" oder sowas
  ierror = get_sys_path(paths)
  ierror = paths%append(".") ! the module ml_mod.py should be in cwd
  write(*,'(a40,i5)') 'Added path with ierror=',ierror
  if(ierror.eq.-1)then ! tell user to get module if not found
    !@thomas TODO adjust name of ml_mod.py
    call env%error("Could not import ml_mod.py. Please copy the file into your current working &
      &directory. It is available at https://github.com/grimme-lab/xtb/tree/main/src/gfnff",source)
  endif

  ! import python ML module
  ierror = import_py(ml_module, "ml_mod") ! omit the .py
  write(*,'(a40,i5)') 'Imported ml_module with ierror=',ierror

  ! create tuple (args) containing arguments for ML python function
  ierror = tuple_create(args, 0)  !@thomas empty tuple needed for call_py
  ! create dict containing kwargs
  ierror = dict_create(kwargs)
  ! add xyz coordinates
  xyz_local=mol%xyz 
  ierror = ndarray_create(xyz_arr, xyz_local)
  ierror = kwargs%setitem("xyz", xyz_arr)
  ! add atom wise energy
  ierror = ndarray_create(eatom_arr, ffml%eatoms)
  ierror = kwargs%setitem("eatoms", eatom_arr)

  ! call ML function from the python module
  ierror = call_py(receive_obj, ml_module, "receive_ml_input_send_output", args, kwargs)
  write(*,'(a40,i5)') 'Called receive_send_fct with ierror=',ierror

  ! unpack received object
  ierror = dict_create(receive_dict)
  ierror = cast(receive_dict, receive_obj)
  ! retrieve energy
  ierror = receive_dict%getitem(obj1, "energy")
  ierror = cast(ml_energy, obj1)
  ! retrieve gradient
  shape_grad(1) = 3  ! define shape 3 dimensions (x y z)
  shape_grad(2) = n  ! n atoms
  ierror = ndarray_create_empty(tmp1, shape_grad)
  ierror = receive_dict%getitem(obj2, "gradient")  ! write gradient into object format
  ierror = cast(tmp1, obj2)                        ! cast object to ndarray
  ierror = tmp1%get_data(tmp_point,'A')                ! get pointer to ndarray data

  ml_gradient = tmp_point                          ! get ml_gradient through pointer

  ! calculate gradient directly in fortran
  do i=1, 3
    do j=1, n
      calc_grad(i,j) = xyz_local(i,j)*ffml%eatoms(j)
    enddo
  enddo

  ! warning for future me !@thomas
  write(*,*) ''
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*) '  when using actual ML python module  '
  write(*,*) '    check what format gradient has    '
  write(*,*) '           (3,n) vs (n,3)             '
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!@thomas test for data transfer
  write(*,*) ''
  write(*,*) '-- -- -- Testing of data transfer -- -- --'
  write(*,'(a,f25.15)') 'Fortran: Receivd ml_energy=',ml_energy
  write(*,'(a,f25.15)') 'Fortran: calculated energy=',SUM(abs(xyz_local))
  write(*,'(a)') 'Fortran: difference between calculated and received gradient'
  write(*,'(a,f20.15)') 'Sum should be zero:', SUM(abs(calc_grad - ml_gradient))
  write(*,*) '-- -- -- -- -- -- -- -- -- -- -- -- -- -- '
  call args%destroy
  call paths%destroy
  call receive_obj%destroy
  call obj1%destroy
  call obj2%destroy
  call receive_dict%destroy
  call tmp1%destroy
  call xyz_arr%destroy
  call ml_module%destroy

end subroutine

end module xtb_gfnff_ffml
