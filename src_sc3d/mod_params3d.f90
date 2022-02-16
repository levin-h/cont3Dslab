!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_directories
!
!---------------------module for directories----------------------------
!
implicit none
!
character(len=100) :: input_file
character(len=100) :: output_file
character(len=11), parameter :: output_dir='outputFILES'
!input_file:      namelist-file with all input-data
!output_file:     name of the final output file (has to end with .h5)
!output_dir:      directory of final output
!output_dir_temp: directory where solution after each iteration step is stored
!all strings are only the directory-name (without '/')
!
!character(len=10), parameter :: model1d_file='model1d.h5'
!character(len=10), parameter :: model2d_file='model2d.h5'
!character(len=10), parameter :: model3d_file='model3d.h5'
!modelxd_file:   file where model-atmosphere is stored
!
!
end module mod_directories
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module options
!
use prog_type  
!
implicit none
!
integer(i4b) :: opt_alo_cont   !(affects only 3d solution scheme)
!opt_alo_cont = 0 if classical lambda iteration shall be used for continuum transfer
!opt_alo_cont = 1 if diagonal alo is used for ALI scheme for continuum transfer
!opt_alo_cont = 2 if direct neighbour alo is used for ALI scheme for continuum transfer (6 neighbours)
!opt_alo_cont = 3 if nearest neighbour alo is used for ALI scheme for continuum transfer (26 neighbours)
!
logical :: opt_ng_cont, opt_ait_cont
!opt_ng_cont = true if ng-extrapolation shall be performed for continuum transfer
!opt_ait_cont = true if aitken extrapolation shall be performed for continuum transfer

integer(i4b) :: opt_method
!
integer(i4b) :: opt_opac
!opt_opac = 0   some analytic opacity laws
!opt_opac = 1   OPAL opacities

integer(i4b) :: opt_epsc
!opt_epsc = 0   eps_cont3d = eps_cont = constant
!opt_opac = 1   eps_cont3d = (chi_tot-chi_thomson)/chi_tot
!
!

end module options
!
!-----------------------------------------------------------------------
!------------------------input (from *.nml file)------------------------
!-----------------------------------------------------------------------
!
module params_input
!
use prog_type
!
implicit none
!
!for model atmosphere and calculation of opacities
real(dp) :: yhe, hei
!
!for modelling the continuum
real(dp) :: kcont, eps_cont
!
!units
real(dp) :: unit_length, unit_density, unit_temperature, unit_velocity, &
     unit_length_cgs
!
!logical to decide for verbose output
logical :: verbose = .true.
!
end module params_input
!
!-----------------------------------------------------------------------
!--------------------------for debugging--------------------------------
!-----------------------------------------------------------------------
!
module mod_debug
!
use prog_type
!
implicit none
!
integer(i4b) :: iindx, kindx
!
end module mod_debug
!
!-----------------------------------------------------------------------
!-------------------------for timing------------------------------------
!-----------------------------------------------------------------------
!
module mod_timing
!
use prog_type
!
implicit none
!
real(dp) :: ts_tot, te_tot, ttot_tot
real(dp) :: ts_formal, te_formal, ttot_formal
real(dp) :: ts_snew, te_snew, ttot_snew
!
!*_tot: total computation time
!*_formal: for each formal solution
!*_snew: for calculating a new source function
!
end module mod_timing
!
!-----------------------------------------------------------------------
!----------------------------for benchmarks-----------------------------
!-----------------------------------------------------------------------
!
module mod_benchmark
  !
  use prog_type
  use fund_const
  
  implicit none
  !
  integer(i4b) :: benchmark_mod
  !
  !for searchlight beam test
  real(dp) :: n_y, n_z, nn_x, nn_y, nn_z, thetab, phib
  real(dp), dimension(:,:,:), allocatable :: int3d_theo
  !
  !for plane-parallel diffusion
  real(dp), parameter :: tau_min=1.d-3, tau_max=1.d6

  end module mod_benchmark
