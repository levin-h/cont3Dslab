!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module options_surfb
!
use prog_type
!  
implicit none
!
character(len=100) :: indat_file, input_file1, input_file2, output_file
!
integer(i4b) :: input_mod
!input_mod=0   1d slab
!input_mod=1   2d slab
!input_mod=2   3d slab
!
integer(i4b) :: opt_opac
!opt_opac=0    continuum opacity = 0
!opt_opac=1    OPAL tables
!opt_opac=2    continuum opacity = thomson-opacity

integer(i4b) :: opt_scont
!opt_scont=0    continuum source function = 0
!opt_scont=1    continuum source function = planck function(trad)
!opt_scont=2    continuum source function = optically thin dilution
!opt_scont=3    N/A
!opt_scont=4    N/A
!opt_scont=5    continuum source function from sc3d.eo (not yet implemented)

integer(i4b) :: opt_sline
!opt_sline=0    line source function = 0
!opt_sline=1    line source function = planck function (trad)
!opt_sline=2    line source function from optically thin approach
!opt_sline=3    line source function from nlte-departure coefficients
!opt_sline=4    line source function from from Sobolev (not yet implemented)
!opt_sline=5    line source function from sc3d.eo (not yet implemented)
!
integer(i4b) :: nxobs_surfb
!nxobs_surfb    number of xobs-points for which surface brightness shall be calculated
!
!
end module options_surfb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_surfb
!
use prog_type
!
implicit none
!
!--------------------dimension of x, y, z grids-------------------------
!
integer(i4b) :: nxobs_fs
integer(i4b) :: nxp, nyp
!
!nxobs_fs: number of used frequency points to calculate surface-brightness
!nxp: number of 'p-rays' in x-dimension
!nyp: number of 'p-rays' in y-dimension
!
end module dime_surfb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module dime_model3d
!
!--------------------------3d-model atmosphere--------------------------
!
use prog_type
!
implicit none
!
integer(i4b) :: nx, ny, nz         !cartesian coordinates
!
real(dp), dimension(:,:,:), allocatable :: rho3d, opac3d, opalbar3d, tgas3d, &
                                           trad3d, velx3d, vely3d, velz3d, &
                                           sline3d, scont3d
!
real(dp), dimension(:), allocatable :: x, y, z
!
end module dime_model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module mod_surfb
!
use prog_type
use fund_const
!
implicit none
!
integer(i4b) :: nz_ray
real(dp), dimension(:), allocatable :: z_ray, opalbar_ray, opac_ray, &
                                       sline_ray, scont_ray, &
                                       profile_ray, velz_ray, vth_ray, &
                                       temp_ray
!$omp threadprivate(nz_ray, z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, velz_ray, &
!$omp               vth_ray, temp_ray)
!nz_ray:      number of data points along an arbitrary ray
!z_ray:       line of sight coordinates of an arbitrary ray
!opalbar_ray: frequency integrated line opacity along ray
!sline_ray:   line-source-function along ray
!profile_ray: profile function along ray
!velz_ray:    line of sight velocity along ray
!vth_ray:     thermal velocity along ray
!temp_ray:    temperature along ray
!
real(dp), parameter :: del_xobs=1.d0/3.d0
!del_xobs:  maximum allowed frequency steps 
!           (shifted frequency from line center in fiducial doppler widths)
!
real(dp), dimension(:), allocatable :: xobs, xnue, xic_nue, xicc_nue, flux_tot, flux_cont, &
                                       flux_emi, flux_abs, femi_x, fabs_x, &
                                       ftot_x, fcont_x, ftot_errx, fcont_errx, &
                                       ftot_err, fcont_err, normt_x, normt, &
                                       acoeff_xic, bcoeff_xic, ccoeff_xic, dcoeff_xic, &
                                       flux_tot_tmp, flux_cont_tmp, normt_tmp, ftot_err_tmp, &
                                       fcont_err_tmp, flux_emi_tmp, flux_abs_tmp
!
!$omp threadprivate(flux_tot_tmp, flux_cont_tmp, normt_tmp, ftot_err_tmp, fcont_err_tmp, &
!$omp               ftot_x, fcont_x, ftot_errx, fcont_errx, normt_x, flux_emi_tmp, flux_abs_tmp, &
!$omp               femi_x, fabs_x)
!
!xobs:      grid of frequency points
!flux_tot:  total emergent flux as function of frequency
!flux_cont: continuum emergent flux as function of frequency
!flux_emi:  only emission part of profile
!flux_abs:  only absorption part of profile
!normt:     test normalization of x-y integration
!ftot_errx:  error of total flux in p-integration (for each xobs and each phi)
!fcont_errx: error of continuum flux in p-integration (for each xobs and each phi)
!ftot_err:   error of total flux in zeta-integration (for each xobs)
!fcont_err:  error of continuum flux in zeta-integration (for each xobs)
!*_tmp:      temporary arrays for integration in parallelized code
!
real(dp) :: relerr_contx, relerr_totx, relerr_cont, relerr_tot, hmax
!$omp threadprivate(relerr_contx, relerr_totx)
!
!relerr_contx: maximum error of continuum flux in x-integration
!relerr_totx:  maximum error of total flux in x-integration
!relerr_cont:  maximumn error of continuum flux in y-integration
!relerr_tot:   maximum error of total flux in y-integration
!hmax:         maximum step size of p-grid
!
real(dp), parameter :: del_vel=1.d0/3.d0
real(dp) :: del_vel2
!del_vel: maximum allowed velocity steps along ray in thermal velocities (in order to resolve resonance zones)
!del_vel2: maximum allowed velocity steps along ray in fiducial thermal velocities
!
real(dp) :: iin, iin_c, iem, iem_c, iemi, iabs
!$omp threadprivate(iin, iin_c, iem, iem_c, iemi, iabs)
!iin:   core intensity (read in planck-function and set iin to that value
!                       or use photospheric profile)
!iin_c: core intensity in the absence of line
!iem:   emergent intensity at outer boundaries (from line and continuum)
!iem_c: emergent intensity at outer boundaries (from continuum alone)
!iemi:  only emission part of profile
!iabs:  only absorption part of profile
!
real(dp) :: phinorm, vphot_proj
!$omp threadprivate(phinorm, vphot_proj)
!phinorm:    normalized profile
!
logical :: lcore
!$omp threadprivate(lcore)
!lcore: logical to describe core and non-core rays
!
real(dp), dimension(:), allocatable :: xp, xpw, xpw1, xpw_err
real(dp), dimension(:), allocatable :: yp, ypw, ypw1, ypw_err
!xp, yp:    x,y 'p-ray' grid
!xpw, ypw:  corresponding integration weights
!
!for the surface brightness
real(dp), dimension(:), allocatable :: xobs_surface, xnue_surface
real(dp), dimension(:,:), allocatable :: iem_surface, iemi_surface, iabs_surface, icont_surface
!iem_surface: emergent intensity on a x-y-surface
!             calculated at xobs_surface
!iemi_surface: only emission part
!iabs_surface: only absorption part
!
end module mod_surfb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
module params_surfb
!
use prog_type
!
implicit none
!
!
!input parameter
real(dp) :: trad, tmin, yhe, hei
real(dp) :: xic1
real(dp) :: vmin, vmax, vmicro, vth_fiducial, vth_min
real(dp) :: eps_cont, kcont, eps_line, kline
real(dp) :: xnue0
!transition frequency of considered line-transition
!
integer(i4b) :: na
!atomic mass number of considered atom (needed for correct calculation of
!   thermal velocity)
!
real(dp) :: gl, gu, flu
!
integer(i4b) :: iline
!
real(dp), parameter :: vth_lim=5.d5
!vth_lim is minimum adopted thermal velocity if vmicro is small and temperature is small (for computational reasons)
!
real(dp), parameter :: xlim=3.d0
!xlim is the maximum allowed shift of (nue_obs-nue_0)/vmax
!
real(dp) :: unit_length, unit_density, unit_velocity, unit_temperature
!units of input model file
!
end module params_surfb
