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
end module options
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module freq
!
!   frequencies
!
use prog_type
!
implicit none
!
! ... scalars
integer(i4b) :: nnue
real(dp) :: xnue0
!
! ... arrays
real(dp) :: xic1_nue
!
!
end module freq
!
!-----------------------------------------------------------------------
!----------------------------1d - grids---------------------------------
!-----------------------------------------------------------------------
!
module dime1d
!
use prog_type
use omp_lib
!
implicit none
!
integer(i4b) :: nz
!
real(dp) :: zmin, zmax
!
real(dp), dimension(:), allocatable :: z
integer, dimension(:), allocatable :: imask1d
!
real(dp), dimension(:), allocatable :: b1d, rho1d, t1d, scont1d, scont1d_diff, scont1d_old, &
                                       mint1d, opac1d, int1d, velz1d, normalization1d, bnue1d, tau1d
!
!for all alo terms
real(dp), dimension(:), allocatable :: alo_diag, alo_subdiag, alo_supdiag
!
!for paralleliztion
!
end module dime1d
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
character(len=100) :: input_file
!
!for model atmosphere
real(dp) :: teff, trad, tmin, rstar, xmloss, vmin, vmax, beta, yhe, hei
!
!for modelling the continuum
real(dp) :: kcont, eps_cont
!
real(dp) :: rstar_cgs, xmloss_cgs, vmin_cgs, vmax_cgs
!
real(dp), parameter :: tau_min=1.d-3, tau_max=1.d6
!
end module params_input
!
!-----------------------------------------------------------------------
!-----------------------iteration-variables-----------------------------
!-----------------------------------------------------------------------
!
module iter
!
use prog_type
implicit none
!
integer(i4b), parameter :: itmaxc=350
integer(i4b) :: it_start_cont
real(dp), parameter :: devmaxc=1.d-5
real(dp), dimension(itmaxc) :: epsmaxc_arr
!
real(dp), dimension(:), allocatable :: solm0, solm1, solm2, solm3
!devmax(l): maximum required percentage-deviation (continuum, line)
!eps1d: error for central-ray mean intensities
!
end module iter
!
!-----------------------------------------------------------------------
!----------------------------ng-extrapolation---------------------------
!-----------------------------------------------------------------------
!
module ng_extra
!
use prog_type
!
implicit none
!
integer(i4b), parameter :: ng_const=6
!
!ng_const: ng (or aitken)-extrapolation is performed at each ng_const iterate
!
end module ng_extra
