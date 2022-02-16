!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_directories
!
implicit none
!
character(len=100) :: indat_file
!
character(len=300) :: model_dir
!
character(len=10), parameter :: model1d_file='model1d.h5'
character(len=10), parameter :: model2d_file='model2d.h5'
character(len=10), parameter :: model3d_file='model3d.h5'
!modelxd_file:   file where model-atmosphere is stored
!
end module mod_directories
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_modext
!
!--------------------dimensions of external atmoshpere------------------
!
use prog_type
!
implicit none
!
integer(i4b) :: nx, ny, nz
!
!
end module dime_modext
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module model1d
!
!----------------------1d model atmpsphere------------------------------
!
use prog_type
!
implicit none
!
real(dp), dimension(:), allocatable :: z_modext1d, velz_modext1d, &
                                       rho_modext1d, t_modext1d
!
!
end module model1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module model2d
!
!----------------------2d model atmpsphere------------------------------
!
use prog_type
!
implicit none
!
real(dp), dimension(:), allocatable :: x_modext2d, z_modext2d
real(dp), dimension(:,:), allocatable :: velx_modext2d, vely_modext2d, velz_modext2d, &
                                         rho_modext2d, t_modext2d
!
end module model2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module model3d
!
!----------------------3d model atmpsphere------------------------------
!
use prog_type
!
implicit none
!
real(dp), dimension(:), allocatable :: x_modext3d, y_modext3d, z_modext3d
real(dp), dimension(:,:,:), allocatable :: velx_modext3d, vely_modext3d, velz_modext3d, &
                                           rho_modext3d, tgas_modext3d, trad_modext3d
!
end module model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_input
!
!--------------------------input parameters-----------------------------
!
use prog_type
!
implicit none
!
integer(i4b) :: input_mod
!
!spatial domain
real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
!
!units
real(dp) :: unit_length, unit_density, unit_velocity, unit_temperature
!
!stellar surface and atmospheric properties
real(dp) :: teff, trad, tmin, rstar, yhe, hei
real(dp) :: rstar_cgs
!
!
!for beta velocity law
real(dp) :: vmin, vmax, beta, xmloss
real(dp) :: vmin_cgs, vmax_cgs, xmloss_cgs
!
!
!
end module params_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_stellar
!
!--------------------------stellar parameters---------------------------
!
use prog_type
!
implicit none
!
real(dp) :: sr
!
end module params_stellar
