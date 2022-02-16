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
!-----------------------diffusion approximation-------------------------
!-----------------------------------------------------------------------
!
module bcondition
!
use prog_type
!
implicit none
!
real(dp) :: xic1
!
end module bcondition
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module angles
!
! angular integration parameters and arrays
!
use prog_type
use fund_const
!
implicit none
!
integer(i4b) :: ntheta
integer(i4b) :: dim_mu, dim_omega
!
real(dp), dimension(:), allocatable :: n_x, n_z, nodes_mu, weight_mu, weight_omega
!
integer, dimension(:,:), allocatable :: q_alo
!
!
!dim_mu: total number of mu-integration-nodes
!nodes_mu: mu-integration-nodes
!weight_mu: mu-integration-weight
!
!n_theta: number of integration nodes for theta=[0,pi/2]
!
!q_alo: indices for alo-calculations for each direction
!       to store nearest neighbours correctly
!
end module angles
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
real(dp), dimension(:), allocatable :: nodes_nue, weight_nue, xic1_nue
!
!
end module freq
!
!-----------------------------------------------------------------------
!----------------------------2d - grids---------------------------------
!-----------------------------------------------------------------------
!
module dime2d
!
use prog_type
use omp_lib
!
implicit none
!
integer(i4b) :: nx, nz
!
real(dp) :: zmin, zmax, xmin, xmax
!
real(dp), dimension(:), allocatable :: x, z
integer, dimension(:,:), allocatable :: imask2d, imaskb2d
!imaskb2d = 0 ghost-zones at iz=1 and iz=nz layer
!imaskb2d = 1 ghost-zones at iz=2 layer
!imaskb2d = 2 ghost-zones at iz=nz-1 layer
!imaskb2d = 3 ghost-zones at ix=1 layer (left left boundary)
!imaskb2d = 4 ghost-zones at ix=2 layer (left boundary)
!imaskb2d = 5 zones at ix=3 layer (first layer in computational domain)
!imaskb2d = 6 zones at ix=nx-2 layer (last layer in computational domain)
!imaskb2d = 7 ghost-zones at ix=nx-1 layer (right boundary)
!imaskb2d = 8 ghost-zones at ix=nx layer (right right boundary)
!imaskb2d = 9 zones within computational domain
!
real(dp), dimension(:,:), allocatable :: b2d, rho2d, t2d, scont2d, mint2d, opac2d, int2d, velz2d, normalization2d, bnue2d, tau2d
!
real(dp), dimension(:,:,:), allocatable :: intbound2d
!
!for all alo terms
real(dp), dimension(:,:,:), allocatable :: alocont_nn2d, alocont_o_nn2d
real(dp), dimension(:), allocatable :: alocont_data, alocont_data_diag, aloline_data, aloline_data_diag
integer(i4b), dimension(:), allocatable :: alocont_colindx, alocont_rowindx, aloline_colindx, aloline_rowindx
!
!for paralleliztion
real(dp), dimension(:,:,:), allocatable :: alocont_nn2d_tmp
real(dp), dimension(:,:), allocatable :: mint2d_tmp, normalization2d_tmp
!$omp threadprivate(int2d, alocont_o_nn2d, alocont_nn2d_tmp, mint2d_tmp)
!
end module dime2d
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
integer(i4b), parameter :: itmaxc=1350
integer(i4b) :: it_start_cont, nconvc
real(dp), parameter :: devmaxc=1.d-5
real(dp), dimension(itmaxc) :: epsmaxc_arr


!
!devmax(l): maximum required percentage-deviation (continuum, line)
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
!
!-----------------------------------------------------------------------
!----------------------------for benchmarks-----------------------------
!-----------------------------------------------------------------------
!
module mod_benchmark
!
use prog_type

implicit none
!
integer(i4b) :: benchmark_mod
!
!for searchlight beam test
real(dp) :: n_z, nn_z, nn_x
real(dp), dimension(:,:), allocatable :: int2d_theo
integer(i4b), parameter :: itmaxi=1000
integer(i4b) :: nconvi
real(dp), parameter :: devmaxi=1.d-5
real(dp), dimension(itmaxi) :: epsmaxi_arr
!
!for plane-parallel diffusion
real(dp), parameter :: tau_min=1.d-3, tau_max=1.d6

end module mod_benchmark
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
