!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--------MAIN PROGRAM TO CALCULATE ATMOSPHERIC STRUCTURE----------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!purpose: calculate atmospheric structure and store everything in
!         an external hdf5-file
!
!v0: including model 0-5:
!     0: plane parallel test model
!     1: read in snapshot from Nicos 2D WR models
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
program main
!
use prog_type
use params_input, only: input_mod
!
implicit none
!
! ... local scalars
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
!
!---------------------read input parameters-----------------------------
!
call read_input
!
!---------------------print input parameters----------------------------
!
call print_input
!
!-------------------calculating model atmosphere------------------------
!
write(*,*) '-------------------------calculate model atmosphere----------------------------'
write(*,*)
!
select case(input_mod)
!
   case(0)
      write(*,*) '------calculating standard model 0--------'
      write(*,*)
      call calc_model3d_pp
      call output_mod3d
!
   case(1)
      write(*,*) '------snapshot from Nicos 2d WR sims------'
      write(*,*)
      call calc_model3d_nicowr2d
      call output_mod3d
!
   case(2)
      write(*,*) '------snapshot from Nicos 3d WR sims------'
      write(*,*)
      call calc_model3d_nicowr3d
      call output_mod3d      
!
!
    case default
       stop 'error: unvalid input-model specified, check input_mod'
!
end select
!
!
end program main
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!---------------------SUBROUTINES AND FUNCTIONS-------------------------
!
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!-------------------read in all input parameters------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
use prog_type
use mod_directories, only: indat_file, model_dir
use fund_const, only: rsu, yr, xmsu  , pi
use params_input, only: input_mod, unit_velocity, unit_density, unit_length, unit_temperature
!
implicit none
!
! ... local scalars
integer(i4b) :: opt_alo_cont, opt_angint_method, opt_grey, opt_opac, opt_epsc, opt_gridxyz
logical :: opt_ng_cont, opt_ait_cont
!
! ... local characters
character(len=100) :: opal_dir
!
! ... local logicals
logical :: verbose
!
! ... local functions
!
! ... namelist
namelist / input_options / model_dir, verbose
namelist / input_model / input_mod
namelist / input_units / unit_velocity, unit_density, unit_length, unit_temperature
!
!-----------------------------------------------------------------------
!
write(*,*) '----------------------------read input-----------------------------------------'
write(*,*) 'input file name (*.nml) to define model-atmosphere'
read(*,*) indat_file
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(indat_file), status='old', form='formatted')
!
!read options
   rewind 1
   read(1, nml=input_options)
   
!read model parameters
   rewind 1
   read(1, nml=input_model)
!
!-----------------------read model units--------------------------------
!
   rewind 1
   read(1, nml=input_units)
!
close(1)
!
!
end subroutine read_input

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_input
!
use prog_type
use fund_const
use params_input
!
IMPLICIT NONE
!
write(*,*) '-------------------------summary of input parameter----------------------------'
write(*,*)
!
write(*,'(a20, i20)') 'input model', input_mod
write(*,*)
write(*,*) 'units'
write(*,'(a20, es20.8)') 'unit_length [r_sun]', unit_length
write(*,'(a20, es20.8)') 'unit_density [g/cm^3]', unit_density
write(*,'(a20, es20.8)') 'unit_temperature [K]', unit_temperature
write(*,'(a20, es20.8)') 'unit_velocity [cm/s]', unit_velocity
write(*,*)
!
!
!
end subroutine print_input
