!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!------------------------------program----------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!v1: formal solution of intensity along z-direction for 3d slab
!    calculating emergent fluxes
!    calculating surface brightness
!
!notes:   including a minimum allowed thermal velocity vth_lim=5 km/s
!         to increase numerical stability
!
!-----------------------------------------------------------------------
!-----------------------start-------------------------------------------
!-----------------------------------------------------------------------
!
program surfb3d
!
use prog_type
use fund_const, only: pi, cgs_clight
use options_surfb, only: input_mod, nxobs_surfb
use omp_lib
!
implicit none
!
! ... local scalars
integer(i4b) :: indx_x, indx_y, indx_xobs, indx_zray
integer(i4b) :: i, err
!
! ... local arrays
!
! ... local logicals
!
! ... local functions
real(dp) :: vthermal
!
!
!
call read_input
!
!
!------------------------read in model atmosphere-----------------------
!
select case(input_mod)
!
!---------------------------3d model------------------------------------
!
   case(2) 
      call read_model3d
      call print_model3d
!
   case default
      stop 'input model not specified'
!
end select
!
!------------------------setting up xobs grid---------------------------
!
call grid_xobs
!
!------------------------setting up p-ray grid--------------------------
!
call grid_xyray(0)
!
!--------------------setting up photospheric profile--------------------
!
call calc_photprof
call output_photprof
!
!--------------check weights for analytical functions-------------------
!
call check_weights
!
!-----------------------calculate surface-brightness--------------------
!
do i=1, nxobs_surfb
   call calc_iem_surface(i)
   call output_surface(i)
enddo
!
!--------------------calculate emergent flux profile--------------------
!
call calc_fluxem
call output_fluxem(1)
call print_flux
call print_ferr
!
!
!-----------------------------------------------------------------------
!
write(*,*) 'program done'
!
end program surfb3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_xobs
!
!-----------------calculates frequency integration nodes----------------
!
use prog_type
use fund_const, only: cgs_clight, cgs_kb, cgs_mp
use dime_surfb, only: nxobs_fs
use mod_surfb, only: del_xobs, xobs, xnue, xic_nue, xicc_nue, xobs_surface, xnue_surface
use params_surfb, only: xnue0, vmax, vmin, vth_fiducial, vth_min, xlim
use options_surfb, only: nxobs_surfb
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: err
real(dp) :: xmax, xmin
real(dp) :: deldop_fiducial, delta, dum_dxobs
!
!calculate ratio of actual (minimum) thermal velocity and scaling used in code
delta=vth_min/vth_fiducial
!
!calculate xmax (corresponding to vmax/vth_min + xlim, but in units of vth_fiducial)
xmin = vmin/vth_fiducial - xlim*delta
xmax = vmax/vth_fiducial + xlim*delta
!
!calculate fiducial doppler width
deldop_fiducial = xnue0 * vth_fiducial/cgs_clight
!
!calculate number of needed frequency points 
!   note: (nue-nue0)/deldop_min = del_xobs
!         nxobs_fs needs to be odd in order to include xobs=0
dum_dxobs = del_xobs*delta
nxobs_fs=2*nint(abs(xmin)/dum_dxobs/2) + 2*nint(xmax/dum_dxobs/2) + 1
!
!write(*,*) nxobs_fs
!write(*,*) vmin, vmax
!write(*,*) xmin, xmax
!write(*,*) xlim
!write(*,*) delta, vth_min, vth_fiducial, 2*nint(xmin/dum_dxobs/2)
!write(*,*) 
!
!allocate xobs-array, xnue-array, xic-nue-array
allocate(xobs(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xobs'
allocate(xnue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xnue'
allocate(xic_nue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xic_nue'
allocate(xicc_nue(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error grid_xobs: xicc_nue'
!
!calculate observers frame frequency grid
!
do i=1, nxobs_fs
   xobs(i) = xmin + (i-1)*(xmax-xmin)/(nxobs_fs-1)
   xnue(i) = xobs(i)*deldop_fiducial + xnue0
enddo
!
!--------------------create xobs-grid for surface brightness----------
!
allocate(xobs_surface(nxobs_surfb), stat=err)
allocate(xnue_surface(nxobs_surfb), stat=err)
!
!
do i=1, nxobs_surfb
   xobs_surface(i) = xmin + (i-1)*(xmax-xmin)/(nxobs_surfb-1)
   xnue_surface(i) = xobs_surface(i)*deldop_fiducial + xnue0   
!   write(*,*) xobs_surface(i), xnue_surface(i)
enddo
!
!
!
end subroutine grid_xobs
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine grid_xyray(opt_grid)
!
!-----------------sets up grid of impact parameter----------------------
!-------------------------for x-y plane---------------------------------
!
use prog_type
use dime_surfb, only: nxp, nyp
use dime_model3d, only: nx, ny, x, y
use mod_surfb, only: xp, xpw, xpw1, xpw_err, yp, ypw, ypw1, ypw_err
use mod_integ1d, only: precalc_weight_trapez_err
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: opt_grid
!
! ... local scalars
integer(i4b) :: i,j,k, ni
integer(i4b) :: err, nx_dum, ny_dum
real(dp) :: del, sum, sum1, sum_err, sumex
!
! ... local arrays
real(dp), dimension(:), allocatable :: x_dum, y_dum
!
!
select case(opt_grid)
   case(0)
!calculate xp-yp directly from input model-grid
      nx_dum=nx
      ny_dum=ny
      allocate(x_dum(nx_dum),stat=err)
      allocate(y_dum(ny_dum),stat=err)
      x_dum=x
      y_dum=y
   case default
      stop 'error in grid_xyray: opt_grid not properly specified'
end select
!
!-------------------include equidistant subintervals--------------------
!------------in order to calculate error estimation weights-------------
!
!add ni equidistant points between each interval
ni=2
nxp = ni*nx_dum-ni+1
nyp = ni*ny_dum-ni+1
allocate(xp(nxp), stat=err)
allocate(xpw(nxp), stat=err)
allocate(xpw1(nxp), stat=err)
allocate(xpw_err(nxp), stat=err)
!
allocate(yp(nyp), stat=err)
allocate(ypw(nyp), stat=err)
allocate(ypw1(nyp), stat=err)
allocate(ypw_err(nyp), stat=err)
!
!for x dimension
k=1
do i=1, nx_dum-1
!
   del = (x_dum(i+1)-x_dum(i))/ni
   do j=1, ni
      xp(k+j-1) = x_dum(i) + float(j-1)*del
   enddo
   k=k+ni
enddo
xp(nxp) = x_dum(nx_dum)
!
!
!
!for y dimension
k=1
do i=1, ny_dum-1
!
   del = (y_dum(i+1)-y_dum(i))/ni
   do j=1, ni
      yp(k+j-1) = y_dum(i) + float(j-1)*del
   enddo
   k=k+ni
enddo
yp(nyp) = y_dum(ny_dum)
!
!-------------------------print out p-grid------------------------------
!
write(*,*) '-----------------------------calculating xp-yp-grid----------------------------'
write(*,*)
!
write(*,*) 'xp-grid:'
write(*,'(8es20.8)') xp
write(*,*)
write(*,*) 'yp-grid:'
write(*,'(8es20.8)') yp
write(*,*)
!
!-----------------------------------------------------------------------
!
!calculating weights for trapezoidal rule including error weights
!  (need to have at least 2 subsequent equidistant subintervals)
call precalc_weight_trapez_err(xp, nxp, ni, xpw, xpw1, xpw_err)
call precalc_weight_trapez_err(yp, nyp, ni, ypw, ypw1, ypw_err)
!
!calculating weights for simpson's rule including error weights
!  (need to have at least 4 subsequent equidistant subintervals)
!call precalc_weight_simps_err(xp, nxp, ni, xpw, xpw1, xpw_err)
!call precalc_weight_simps_err(yp, nyp, ni, ypw, ypw1, ypw_err)
!
!calculating weights for trapezoidal rule
!   (no need of equidistant subintervals since no error weights are calculated)
!call precalc_weight_trapez(xp, nxp, xpw)
!call precalc_weight_trapez(yp, nyp, ypw)
!
!------------------------test integration-------------------------------
!
!integrate f(x,y)=1
sumex = (xp(nxp)-xp(1))*(yp(nyp)-yp(1))
sum=0.d0
sum1=0.d0
sum_err=0.d0
do i=1, nxp
   do j=1, nyp
      sum = sum + xpw(i) * ypw(j)
      sum1 = sum1 + xpw1(i) * ypw1(j)
      sum_err = sum_err + xpw_err(i)*ypw_err(j)
   enddo
enddo
!
!
if(abs(sumex-sum)/sumex.gt.1.d-10) then
   write(*,*) '--------------error in xp-yp-integration---------'
   write(*,'(6a30)') 'integral(num)', 'integral(num), half step', 'integral(exact)', 'abs error_1(num)', 'abs error_2(num', 'error(exact)'
   write(*,'(6(e30.8))') sum, sum1, sumex, sum_err, (sum-sum1)/6.d0, (sumex-sum)
!error_2(num) is (sum-sum1)/6.d0  for trapezoidal rule and
!                (sum-sum1)/15.d0 for simpson's rule
   stop
endif
!
!
!
end subroutine grid_xyray
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_photprof
!
!----------set photospheric profile to constant value: xic1-------------
!----------set photospheric continuum to constant value: xic1-----------
!
use prog_type
use fund_const, only: cgs_clight
use params_surfb, only: xic1
use mod_surfb, only: xic_nue, xicc_nue
!
implicit none
!
! ... arguments
!
! ... local scalars
!
! ... local functions
!

write(*,*) '----------------------no photospheric profile is included----------------------'
write(*,*)
!
xic_nue=xic1
xicc_nue=xic1
!
end subroutine calc_photprof
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_ray3d(xp, yp, lcore)
!
!------sets up a ray along x for the given xp,yp coordinate-------------
!
!to calculate:
!   1. z-ray         from z_ray = z of 3d input coordinate system
!   3. opalbar_ray   trilinear interpolation from given 3d-grid
!   4. sline_ray     trilinear interpolation from given 3d-grid
!   5. velx(y,z)_ray trilinear interpolation from given 3d-grid
!   5. velz_ray      projected velocity along the ray
!
!input: xp, yp  
!
!output: lcore    logical to describe if core or non-core ray
!        all physical quantities along ray
!
use prog_type
!use dime_surfb, only:  n
use dime_model3d, only: nx, ny, nz, x, y, z, &
                        opac3d, scont3d, opalbar3d, sline3d, velx3d, vely3d, velz3d, tgas3d
use mod_surfb, only: nz_ray, z_ray, opalbar_ray, opac_ray, scont_ray, sline_ray, &
                     vth_ray, temp_ray,  profile_ray, velz_ray, &
                     del_vel2, vphot_proj
use params_surfb, only: vth_fiducial, vth_min, vth_lim, xic1, vmin, vmax, tmin, na, vmicro
use mod_interp1d, only: interpol_yp_spline, interpol_yp, find_index
use mod_interp3d, only: coeff3d_8p_lin
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xp, yp
logical, intent(out) :: lcore
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: iz, iz_dum, nadd
integer(i4b) :: iim2, iim1, ii, iip1, &
                jjm2, jjm1, jj, jjp1, &
                kkm2, kkm1, kk, kkp1
integer(i4b), parameter :: nz_max = 20000
real(dp) :: acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
real(dp) :: velx, vely, velz, dum_vel2, dum_vel1, dum_delv
!!
!! ... local arrays
real(dp), dimension(:), allocatable :: zdum_ray, veldum_ray, vthdum_ray, &
                               opacdum_ray, opalbardum_ray, scontdum_ray, slinedum_ray, tempdum_ray
real(dp), dimension(:), allocatable :: zdum2_ray, veldum2_ray, vthdum2_ray, opacdum2_ray, opalbardum2_ray, scontdum2_ray, &
                               slinedum2_ray, tempdum2_ray
!!
!! ... local logicals
!logical :: expol, linfo_phot, linfo_max, llogr, llogt, llogp, llogf, lr2, ldum
!!
!! ... local functions
real(dp) :: vthermal
!real(dp) :: calc_vmicro
!logical :: boundary
!!

allocate(zdum_ray(nz_max), veldum_ray(nz_max), vthdum_ray(nz_max), &
     opacdum_ray(nz_max), opalbardum_ray(nz_max), &
     scontdum_ray(nz_max), slinedum_ray(nz_max), tempdum_ray(nz_max))
allocate(zdum2_ray(nz_max), veldum2_ray(nz_max), vthdum2_ray(nz_max), &
     opacdum2_ray(nz_max), opalbardum2_ray(nz_max), scontdum2_ray(nz_max), &
     slinedum2_ray(nz_max), tempdum2_ray(nz_max))

!!
iz=0
!
!
do i=nz, 1, -1
   zdum_ray(i) = z(nz-i+1)
   iz=iz+1
enddo
!
!check if core ray (for now, always core ray)
lcore=.true.
!
!for now, no photospheric rotation
vphot_proj = 0.d0
!
!find indices for given xp, yp
call find_index(xp, x, nx, iim2, iim1, ii, iip1)
call find_index(yp, y, ny, jjm2, jjm1, jj, jjp1)
!
!-----------------------------------------------------------------------
!
do i=1, iz
!
!calculate z_ray in carthesian coordinates
!
!search for indices of a the surrounding grid-cell (or neighbouring grid-cells for extrapolation)
   call find_index(zdum_ray(i), z, nz, kkm2, kkm1, kk, kkp1)

   call coeff3d_8p_lin(x(iim1), x(ii), y(jjm1), y(jj), z(kkm1), z(kk), xp, yp, zdum_ray(i), &
                       acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff)
!
!
!----------------------interpolation of velocity components-------------
!
   velx = acoeff*velx3d(iim1,jjm1,kkm1) + bcoeff*velx3d(ii,jjm1,kkm1) + &
          ccoeff*velx3d(iim1,jj,kkm1) + dcoeff*velx3d(ii,jj,kkm1) + &
          ecoeff*velx3d(iim1,jjm1,kk) + fcoeff*velx3d(ii,jjm1,kk) + &
          gcoeff*velx3d(iim1,jj,kk) + hcoeff*velx3d(ii,jj,kk)
   
   vely = acoeff*vely3d(iim1,jjm1,kkm1) + bcoeff*vely3d(ii,jjm1,kkm1) + &
          ccoeff*vely3d(iim1,jj,kkm1) + dcoeff*vely3d(ii,jj,kkm1) + &
          ecoeff*vely3d(iim1,jjm1,kk) + fcoeff*vely3d(ii,jjm1,kk) + &
          gcoeff*vely3d(iim1,jj,kk) + hcoeff*vely3d(ii,jj,kk)   

   velz = acoeff*velz3d(iim1,jjm1,kkm1) + bcoeff*velz3d(ii,jjm1,kkm1) + &
          ccoeff*velz3d(iim1,jj,kkm1) + dcoeff*velz3d(ii,jj,kkm1) + &
          ecoeff*velz3d(iim1,jjm1,kk) + fcoeff*velz3d(ii,jjm1,kk) + &
          gcoeff*velz3d(iim1,jj,kk) + hcoeff*velz3d(ii,jj,kk)
!
!calculation of velocity projected onto ray
   veldum_ray(i) = velz
!
!----------------------interpolation of line opacity--------------------
!
!!line opacity in units of 1/unit_length
   opalbardum_ray(i) = acoeff*opalbar3d(iim1,jjm1,kkm1) + bcoeff*opalbar3d(ii,jjm1,kkm1) + &
                       ccoeff*opalbar3d(iim1,jj,kkm1) + dcoeff*opalbar3d(ii,jj,kkm1) + &
                       ecoeff*opalbar3d(iim1,jjm1,kk) + fcoeff*opalbar3d(ii,jjm1,kk) + &
                       gcoeff*opalbar3d(iim1,jj,kk) + hcoeff*opalbar3d(ii,jj,kk)
!
!------------------------line source function---------------------------
!
!get sline on cube vertices
   slinedum_ray(i) = acoeff*sline3d(iim1,jjm1,kkm1) + bcoeff*sline3d(ii,jjm1,kkm1) + &
                     ccoeff*sline3d(iim1,jj,kkm1) + dcoeff*sline3d(ii,jj,kkm1) + &
                     ecoeff*sline3d(iim1,jjm1,kk) + fcoeff*sline3d(ii,jjm1,kk) + &
                     gcoeff*sline3d(iim1,jj,kk) + hcoeff*sline3d(ii,jj,kk)
!
!----------------------interpolation of continuum opacity---------------
!
   opacdum_ray(i) = acoeff*opac3d(iim1,jjm1,kkm1) + bcoeff*opac3d(ii,jjm1,kkm1) + &
                    ccoeff*opac3d(iim1,jj,kkm1) + dcoeff*opac3d(ii,jj,kkm1) + &
                    ecoeff*opac3d(iim1,jjm1,kk) + fcoeff*opac3d(ii,jjm1,kk) + &
                    gcoeff*opac3d(iim1,jj,kk) + hcoeff*opac3d(ii,jj,kk)
!
!---------------------continuum source function-------------------------
!
   scontdum_ray(i) = acoeff*scont3d(iim1,jjm1,kkm1) + bcoeff*scont3d(ii,jjm1,kkm1) + &
                     ccoeff*scont3d(iim1,jj,kkm1) + dcoeff*scont3d(ii,jj,kkm1) + &
                     ecoeff*scont3d(iim1,jjm1,kk) + fcoeff*scont3d(ii,jjm1,kk) + &
                     gcoeff*scont3d(iim1,jj,kk) + hcoeff*scont3d(ii,jj,kk)
!
!---------------------temperature---------------------------------------
!
   tempdum_ray(i) = acoeff*tgas3d(iim1,jjm1,kkm1) + bcoeff*tgas3d(ii,jjm1,kkm1) + &
                     ccoeff*tgas3d(iim1,jj,kkm1) + dcoeff*tgas3d(ii,jj,kkm1) + &
                     ecoeff*tgas3d(iim1,jjm1,kk) + fcoeff*tgas3d(ii,jjm1,kk) + &
                     gcoeff*tgas3d(iim1,jj,kk) + hcoeff*tgas3d(ii,jj,kk)
!   tempdum_ray(i)=0.d0
!   
!-----------------------------------------------------------------------   
!
!calculate corresponding thermal velocity (for a fixed tmin)
!      vthdum_ray(i) = vthermal(vmicro, tmin, na)   
   vthdum_ray(i) = max(vth_lim,vthermal(vmicro, tempdum_ray(i), na))
!   write(*,*) vthdum_ray(i)/1.d5
!
!-----------------------------------------------------------------------
!
!
enddo
!
!-----------------------------------------------------------------------
!
!----------------check if resonance zones are resolved------------------
!-----------if not resolved: add grid points and interpolate------------
!
!velocity at beginning point
dum_vel1=veldum_ray(1)
!
veldum2_ray(1)=veldum_ray(1)
zdum2_ray(1)=zdum_ray(1)
opacdum2_ray(1)=opacdum_ray(1)
opalbardum2_ray(1)=opalbardum_ray(1)
scontdum2_ray(1)=scontdum_ray(1)
slinedum2_ray(1)=slinedum_ray(1)
tempdum2_ray(1)=tempdum_ray(1)
vthdum2_ray(1)=vthdum_ray(1)
!
k=0
iz_dum=iz
!
!
do i=2, iz
!
   dum_vel2=veldum_ray(i)
!
!calculation of velocity steps
   dum_delv=dum_vel2-dum_vel1
!
!  if(dum_vel1*dum_vel2.lt.0.d0) then
!     write(*,'(a96)') 'error in setup_ray_scp: velocity jump along ray from positive velocities <-> negative velocities'
!     write(*,'(a82)') '   => 1. trust the interpolation scheme for shear velocites, etc (not recommended)'
!     write(*,'(a68)') '         at least, you should check against an analytic velocity law'
!     write(*,'(a82)') '   => 2. use an analytic description of the velocity law, instead of interpolation'
!     stop
!  endif
!
   if(abs(dum_delv).gt.del_vel2) then
!always round up note: ceiling is only f90
      nadd = ceiling(abs(dum_delv)/del_vel2)

      do j=1, nadd-1
!add nadd-1 additional points in velocity space
         veldum2_ray(i+k)=veldum_ray(i-1) + j*(veldum_ray(i)-veldum_ray(i-1))/nadd
!interpolation of all variables
         zdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), zdum_ray(i-1), zdum_ray(i), veldum2_ray(i+k))
         opacdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opacdum_ray(i-1), opacdum_ray(i), veldum2_ray(i+k))
         opalbardum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), opalbardum_ray(i-1), opalbardum_ray(i), veldum2_ray(i+k))
         scontdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), scontdum_ray(i-1), scontdum_ray(i), veldum2_ray(i+k))
         slinedum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), slinedum_ray(i-1), slinedum_ray(i), veldum2_ray(i+k))
         tempdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), tempdum_ray(i-1), tempdum_ray(i), veldum2_ray(i+k))
         vthdum2_ray(i+k)=interpol_yp(veldum_ray(i-1), veldum_ray(i), vthdum_ray(i-1), vthdum_ray(i), veldum2_ray(i+k))
!
         k=k+1
         iz_dum=iz_dum+1
      enddo
   endif

   veldum2_ray(i+k)=veldum_ray(i)
   zdum2_ray(i+k)=zdum_ray(i)
   opacdum2_ray(i+k)=opacdum_ray(i)
   opalbardum2_ray(i+k)=opalbardum_ray(i)
   scontdum2_ray(i+k)=scontdum_ray(i)
   slinedum2_ray(i+k)=slinedum_ray(i)
   tempdum2_ray(i+k)=tempdum_ray(i)
   vthdum2_ray(i+k)=vthdum_ray(i)

   dum_vel1=dum_vel2
! 
enddo

!stop 'go on here'
!
!------------------store everything in global arrays--------------------
!
nz_ray=iz_dum
!
call allocate_fs1d
!
do i=1, nz_ray
   z_ray(i) = zdum2_ray(nz_ray+1-i)
   opalbar_ray(i) = opalbardum2_ray(nz_ray+1-i)
   opac_ray(i) = opacdum2_ray(nz_ray+1-i)
   scont_ray(i) = scontdum2_ray(nz_ray+1-i)
   sline_ray(i) = slinedum2_ray(nz_ray+1-i)
   vth_ray(i) = vthdum2_ray(nz_ray+1-i)
   temp_ray(i) = tempdum2_ray(nz_ray+1-i)
   velz_ray(i) = veldum2_ray(nz_ray+1-i)
enddo
!
!
!
end subroutine setup_ray3d
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine get_iin(xobsi, lcore, iin, iin_c)
!
use prog_type
use fund_const, only: cgs_clight
use dime_surfb, only: nxobs_fs
use params_surfb, only: xnue0, xic1
use mod_surfb, only: xic_nue, xicc_nue, xnue, vphot_proj, acoeff_xic, bcoeff_xic, ccoeff_xic, dcoeff_xic, xobs
use omp_lib
use mod_interp1d, only: interpol_yp, interpol_yp_spline, find_index
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xobsi
logical, intent(in) :: lcore
real(dp), intent(out) :: iin, iin_c
!
! ... local scalars
integer(i4b) ::  iim2, iim1, ii, iip1
real(dp) :: nue_cmf, xcmf, delta
!
! ... local functions


if(lcore) then
!note: need to interpolate from photospheric profile, since this is only given in cmf, 
!      whereas i am calculating in observers frame (shift of profile due to rotational velocities)
!calculate comoving frame frequency
   xcmf = xobsi - vphot_proj
!
!constant continuum
   iin_c = xic1
!
!find index within xobs-grid, for which xcmf matches and interpolate
   call find_index(xcmf, xobs, nxobs_fs, iim2, iim1, ii, iip1)
!
   iin = interpol_yp(xobs(iim1), xobs(ii), xic_nue(iim1), xic_nue(ii), xcmf)
   iin_c = interpol_yp(xobs(iim1), xobs(ii), xicc_nue(iim1), xicc_nue(ii), xcmf) !for frequency dependent continuum
!
else
!set non-core ray
   iin=0.d0
   iin_c=0.d0
endif

!
end subroutine get_iin
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine print_ray(xobs, zeta, p, iin, iin_c, fname)
!!
!use prog_type
!use mod_surfb, only: nhat, nz_ray, z_ray, opac_ray, opalbar_ray, velz_ray, &
!                        scont_ray, sline_ray, temp_ray, profile_ray
!!
!implicit none
!!
!! ... arguments
!real(dp), intent(in) :: xobs, p, zeta, iin, iin_c
!character(len=*), intent(in) :: fname
!!
!! ... local scalars
!integer(i4b) :: i
!!
!write(*,*) '-----------atmospheric structure along ray-------------'
!!
!open(1,file=trim(fname))
!   write(1,'(a15, 3e20.8)') 'for nhat:    ', nhat
!   write(1,'(a15, 2e20.8)') 'for p, zeta: ', p ,zeta
!   write(1,'(a15, e20.8)')  'for xobs:    ', xobs
!   write(1,'(a15, e20.8)')  'boundary iin ', iin
!   write(1,'(a15, e20.8)')  'boundary iin_c ', iin_c
!   write(1,*)
!   write(1,'(a4, 8(a20))') '#', 'z [r_star]', 'opac [1/cm]', 'opalbar [1/cm]', 'vel_z [cm/s]', 's_cont', 's_line', 'temp', 'profile'
!   do i=1, nz_ray
!      write(1,'(i4, 8(e20.8))')  i, z_ray(i), opac_ray(i), opalbar_ray(i), velz_ray(i), scont_ray(i), sline_ray(i), temp_ray(i), profile_ray(i)
!   enddo
!   write(1,*)
!close(1)
!!
!write(*,'(a15, 3e20.8)') 'for nhat:    ', nhat
!write(*,'(a15, 2e20.8)') 'for p, zeta: ', p ,zeta
!write(*,'(a15, e20.8)')  'for xobs:    ', xobs
!write(*,'(a15, e20.8)')  'boundary iin ', iin
!write(*,'(a15, e20.8)')  'boundary iin_c ', iin_c
!write(*,*)
!write(*,'(a4, 8(a20))') '#', 'z [r_star]', 'opac [1/cm]', 'opalbar [1/cm]', 'vel_z [cm/s]', 's_cont', 's_line', 'temp', 'profile'
!do i=1, nz_ray
!   write(*,'(i4, 8(e20.8))')  i, z_ray(i), opac_ray(i), opalbar_ray(i), velz_ray(i), scont_ray(i), sline_ray(i), temp_ray(i), profile_ray(i)
!enddo
!write(*,*)
!
!end subroutine print_ray
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine print_flux
!
use prog_type
use fund_const, only: pi
use dime_surfb, only: nxobs_fs
use mod_surfb, only: flux_tot, flux_cont, xobs, normt
use params_surfb, only: vmax, vth_fiducial
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i
!
!xobs in units of vmax, insted of vthfiducial
!
write(*,*) '-------------------------emergent flux profile---------------------------------'
!
write(*,'(a4, 5(a20))') '#', 'xobs', 'flux_tot', 'flux_cont', 'f_tot/f_cont', 'normt'
do i=1, nxobs_fs
   write(*,'(i4, 5(e20.8))')  i, xobs(i)*vth_fiducial/vmax, flux_tot(i), flux_cont(i), &
                              flux_tot(i)/flux_cont(i), normt(i)
enddo
write(*,*)
!
end subroutine print_flux
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_ferr
!
use prog_type
use dime_surfb, only: nxp, nyp
use mod_surfb, only: relerr_cont, relerr_tot, relerr_contx, relerr_totx
!
implicit none
!
! ... local scalars
!

write(*,*) '---------------------error estimation of flux integral-------------------------'
write(*,'(a50, 2es20.8)') 'max rel error x-integ f_tot, f_cont', relerr_totx, relerr_contx
write(*,'(a50, 2es20.8)') 'max rel error y-integ f_tot, f_cont', relerr_tot, relerr_cont
write(*,'(a50, i10, i10, es20.8)') '# xp/yp-grid-points', nxp, nyp
write(*,*)
!
!
end subroutine print_ferr
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_err(nd, sol1, sol2, abserr1, abserr2, relerr1, relerr2)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nd
real(dp), dimension(nd), intent(in) :: sol1, sol2, abserr1, abserr2
real(dp) :: relerr1, relerr2
!
! ... local scalars
integer(i4b) :: i
real(dp) :: dum_err1, dum_err2
!
relerr1=0.d0
relerr2=0.d0
!
do i=1, nd
!
   if(sol1(i).eq.0.d0) then
      dum_err1=0.d0
   else
      dum_err1=abs(abserr1(i)/sol1(i))
   endif
!
   if(dum_err1.gt.relerr1) then
      relerr1=dum_err1
   endif
!
   if(sol2(i).eq.0.d0) then
      dum_err2=0.d0
   else
      dum_err2=abs(abserr2(i)/sol2(i))
   endif
!
   if(dum_err2.gt.relerr2) then
      relerr2=dum_err2
   endif
!
enddo
!
!
!
end subroutine calc_err
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
subroutine allocate_fluxes
!
use prog_type
use dime_surfb, only: nxobs_fs, nyp
use mod_surfb, only: flux_tot_tmp, flux_cont_tmp, ftot_x, fcont_x, ftot_errx, fcont_errx, &
                     normt_x, normt_tmp, ftot_err_tmp, fcont_err_tmp, flux_emi_tmp, flux_abs_tmp, &
                     femi_x, fabs_x
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
if(allocated(flux_tot_tmp)) deallocate(flux_tot_tmp)
if(allocated(flux_cont_tmp)) deallocate(flux_cont_tmp)
if(allocated(normt_tmp)) deallocate(normt_tmp)
if(allocated(ftot_err_tmp)) deallocate(ftot_err_tmp)
if(allocated(fcont_err_tmp)) deallocate(ftot_err_tmp)
if(allocated(ftot_x)) deallocate(ftot_x)
if(allocated(fcont_x)) deallocate(fcont_x)
if(allocated(ftot_errx)) deallocate(ftot_errx)
if(allocated(fcont_errx)) deallocate(fcont_errx)
if(allocated(normt_x)) deallocate(normt_x)
if(allocated(flux_emi_tmp)) deallocate(flux_emi_tmp)
if(allocated(flux_abs_tmp)) deallocate(flux_abs_tmp)
if(allocated(femi_x)) deallocate(femi_x)
if(allocated(fabs_x)) deallocate(fabs_x)
!
!
allocate(flux_tot_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_tot_tmp'
!
allocate(flux_cont_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_cont_x'
!
allocate(normt_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: normt_tmp'
!
allocate(ftot_err_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: ftot_err_tmp'
!
allocate(fcont_err_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: fcont_err_tmp'
!
allocate(ftot_x(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: ftot_x'
!
allocate(fcont_x(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: fcont_x'
!
allocate(ftot_errx(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: ftot_errx'
!
allocate(fcont_errx(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: fcont_errx'
!
allocate(normt_x(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: normt_x'
!
allocate(flux_emi_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_emi_tmp'
!
allocate(flux_abs_tmp(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: flux_abs_tmp'
!
allocate(femi_x(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: femi_x'
!
allocate(fabs_x(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error allocate_fluxes: fabs_x'
!
!
!
end subroutine allocate_fluxes
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine deallocate_fluxes
!
use prog_type
use dime_surfb, only: nxobs_fs, nyp
use mod_surfb, only: flux_tot_tmp, flux_cont_tmp, ftot_x, fcont_x, ftot_errx, fcont_errx, &
                        normt_x, normt_tmp, ftot_err_tmp, fcont_err_tmp, flux_emi_tmp, flux_abs_tmp, &
                        femi_x, fabs_x
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
deallocate(flux_tot_tmp)
deallocate(flux_cont_tmp)
deallocate(normt_tmp)
deallocate(ftot_err_tmp)
deallocate(fcont_err_tmp)
deallocate(ftot_x)
deallocate(fcont_x)
deallocate(ftot_errx)
deallocate(fcont_errx)
deallocate(normt_x)
deallocate(flux_emi_tmp)
deallocate(flux_abs_tmp)
deallocate(femi_x)
deallocate(fabs_x)
!
!
!
end subroutine deallocate_fluxes
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_weights
!
!-------------------testing weights and error weights-------------------
!----------------------------for xp-yp-grids----------------------------
!
use prog_type
use fund_const, only: pi
use dime_surfb, only: nxp, nyp
use mod_surfb, only: xp, xpw, xpw1, xpw_err, yp, ypw, ypw1, ypw_err
!
implicit none
!
! ... local scalars
integer(i4b) :: err
real(dp) :: sum1, sum2, sumerr, sumex
!
! ... local arrays
real(dp), dimension(:), allocatable :: yvalue
!
!
write(*,*) '--------------------checking weights and error weights-------------------------'
write(*,*)
!
allocate(yvalue(nxp), stat=err)
!
write(*,*) 'integration in x'
write(*,'(6a20)') 'function', 'half step', 'full step', 'exact', 'error(num)', 'error(exact)'
!
!linear function in x
yvalue=xp
sumex=0.5d0*(xp(nxp)**2.d0 - xp(1)**2.d0)
sum1=sum(xpw*yvalue)
sum2=sum(xpw1*yvalue)
sumerr=sum(xpw_err*yvalue)
write(*,'(a20, 5es20.8)') 'linear', sum1, sum2, sumex, sumerr, sumex-sum1
!
!quadratic function
yvalue = xp*xp
sumex=1.d0/3.d0 * (xp(nxp)**3.d0 - xp(1)**3.d0)
sum1=sum(xpw*yvalue)
sum2=sum(xpw1*yvalue)
sumerr=sum(xpw_err*yvalue)
write(*,'(a20, 5es20.8)') 'quadratic', sum1, sum2, sumex, sumerr, sumex-sum1
!
!cubic function
yvalue=xp*xp*xp
sumex=1.d0/4.d0 * (xp(nxp)**4.d0 - xp(1)**4.d0)
sum1=sum(xpw*yvalue)
sum2=sum(xpw1*yvalue)
sumerr=sum(xpw_err*yvalue)
write(*,'(a20, 5es20.8)') 'cubic', sum1, sum2, sumex, sumerr, sumex-sum1
!
!to the four
yvalue=xp*xp*xp*xp / xp(nxp)
sumex=1.d0/5.0/xp(nxp) * (xp(nxp)**5.d0 - xp(1)**5.d0)
sum1=sum(xpw*yvalue)
sum2=sum(xpw1*yvalue)
sumerr=sum(xpw_err*yvalue)
write(*,'(a20, 5es20.8)') 'to the four', sum1, sum2, sumex, sumerr, sumex-sum1
!
!exponential
yvalue=exp(xp/xp(nxp)*5.d0)
sumex=(exp(5.d0)-exp(xp(1)/xp(nxp)*5.d0))*xp(nxp)/5.d0
sum1=sum(xpw*yvalue)
sum2=sum(xpw1*yvalue)
sumerr=sum(xpw_err*yvalue)
write(*,'(a20, 5es20.8)') 'exponential', sum1, sum2, sumex, sumerr, sumex-sum1
!
!cosine function
yvalue=cos(xp/xp(nxp) * 2.d0*pi)
sumex=(sin(2.d0*pi)-sin(xp(1)/xp(nxp) * 2.d0*pi))*xp(nxp)/(2.d0*pi)
sum1=sum(xpw*yvalue)
sum2=sum(xpw1*yvalue)
sumerr=sum(xpw_err*yvalue)
write(*,'(a20, 5es20.8)') 'cos', sum1, sum2, sumex, sumerr, sumex-sum1
write(*,*)
!
!------------------------------------------------------------------------
!
deallocate(yvalue)
allocate(yvalue(nyp), stat=err)
!
write(*,*) 'integration in y'
write(*,'(6a20)') 'function', 'half step', 'full step', 'exact', 'error(num)', 'error(exact)'
!
!linear function in y
yvalue=yp
sumex=0.5d0*(yp(nyp)**2.d0 - yp(1)**2.d0)
sum1=sum(ypw*yvalue)
sum2=sum(ypw1*yvalue)
sumerr=sum(ypw_err*yvalue)
write(*,'(a20, 5es20.8)') 'linear', sum1, sum2, sumex, sumerr, sumex-sum1
!
!quadratic function
yvalue = yp*yp
sumex=1.d0/3.d0 * (yp(nyp)**3.d0 - yp(1)**3.d0)
sum1=sum(ypw*yvalue)
sum2=sum(ypw1*yvalue)
sumerr=sum(ypw_err*yvalue)
write(*,'(a20, 5es20.8)') 'quadratic', sum1, sum2, sumex, sumerr, sumex-sum1
!
!cubic function
yvalue=yp*yp*yp
sumex=1.d0/4.d0 * (yp(nyp)**4.d0 - yp(1)**4.d0)
sum1=sum(ypw*yvalue)
sum2=sum(ypw1*yvalue)
sumerr=sum(ypw_err*yvalue)
write(*,'(a20, 5es20.8)') 'cubic', sum1, sum2, sumex, sumerr, sumex-sum1
!
!to the four
yvalue=yp*yp*yp*yp / yp(nyp)
sumex=1.d0/5.0/yp(nyp) * (yp(nyp)**5.d0 - yp(1)**5.d0)
sum1=sum(ypw*yvalue)
sum2=sum(ypw1*yvalue)
sumerr=sum(ypw_err*yvalue)
write(*,'(a20, 5es20.8)') 'to the four', sum1, sum2, sumex, sumerr, sumex-sum1
!
!exponential
yvalue=exp(yp/yp(nyp)*5.d0)
sumex=(exp(5.d0)-exp(yp(1)/yp(nyp)*5.d0))*yp(nyp)/5.d0
sum1=sum(ypw*yvalue)
sum2=sum(ypw1*yvalue)
sumerr=sum(ypw_err*yvalue)
write(*,'(a20, 5es20.8)') 'exponential', sum1, sum2, sumex, sumerr, sumex-sum1
!
!cosine function
yvalue=cos(yp/yp(nyp) * 2.d0*pi)
sumex=(sin(2.d0*pi)-sin(yp(1)/yp(nyp) * 2.d0*pi))*yp(nyp)/(2.d0*pi)
sum1=sum(ypw*yvalue)
sum2=sum(ypw1*yvalue)
sumerr=sum(ypw_err*yvalue)
write(*,'(a20, 5es20.8)') 'cos', sum1, sum2, sumex, sumerr, sumex-sum1
write(*,*)
!
!
end subroutine check_weights
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_fluxem
!
use prog_type
use fund_const
use options_surfb, only: input_mod
use dime_surfb, only: nxobs_fs, nxp, nyp
use mod_surfb, only: flux_tot, flux_cont, normt, ftot_err, fcont_err, flux_emi, flux_abs, &
     fabs_x, fcont_err_tmp, fcont_errx, fcont_x, femi_x, flux_abs_tmp, flux_cont_tmp, &
     flux_emi_tmp, flux_tot_tmp, ftot_err_tmp, ftot_errx, ftot_x, &
     normt_tmp, normt_x, lcore, xp, yp, xpw, ypw, xpw_err, ypw_err, xobs, &
     iabs, iem, iem_c, iemi, iin, iin_c, &
     nz_ray, opac_ray, opalbar_ray, phinorm, scont_ray, sline_ray, z_ray, velz_ray, profile_ray, vth_ray, &
     relerr_contx, relerr_totx, relerr_cont, relerr_tot
use params_surfb, only: vth_fiducial
use mod_math, only: calc_phinorm
!
implicit none
!
! ... local scalars
integer(i4b) :: indx_xp, indx_yp, indx_xobs, indx_zray
integer(i4b) :: err, i
real(dp) :: err1, err2

!--------------------allocate fluxes (not threadprivate)----------------
!
allocate(flux_tot(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: flux_tot'
allocate(flux_cont(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: flux_cont'
allocate(normt(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: normt'
allocate(ftot_err(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: ftot_err'
allocate(fcont_err(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: fcont_err'
allocate(flux_emi(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: flux_emi'
allocate(flux_abs(nxobs_fs), stat=err)
   if(err.ne.0) stop 'allocation error: flux_abs'
!   
!
write(*,*) '----------------------calculating emergent flux profile------------------------'
write(*,*)
!
!----------------set up a single xp,yp ray from the grid----------------
!
flux_tot=zero
flux_cont=zero
normt=zero
ftot_err=zero
fcont_err=zero
flux_emi=zero
flux_abs=zero
!
!$omp parallel &
!$omp private(indx_xp, indx_yp, indx_xobs, indx_zray, err1, err2)
!
!allocation of global (threadprivate) arrays
call allocate_fluxes
!
flux_tot_tmp=zero 
flux_cont_tmp=zero
flux_emi_tmp=zero
flux_abs_tmp=zero
normt_tmp=zero
ftot_err_tmp=zero
fcont_err_tmp=zero

!$omp do schedule(dynamic)
do indx_yp=1, nyp
!
   write(*,'(a30, i5, a2, i5)') 'calculating indx_y (nyp)', indx_yp, '/', nyp
   ftot_x=zero
   fcont_x=zero
   ftot_errx=zero
   fcont_errx=zero
   normt_x=zero
   femi_x=zero
   fabs_x=zero
!
   do indx_xp=1, nxp
!
      select case(input_mod)            
!         case(0)
!            call setup_ray1d(zeta(indx_zeta), p(indx_p), lcore)
!         case(1)
!             call setup_ray3d_rot(zeta(indx_zeta), p(indx_p), lcore)
!            call setup_ray3d(zeta(indx_zeta), p(indx_p), lcore)
         case(2)
            call setup_ray3d(xp(indx_xp),yp(indx_yp),lcore)
         case default
            stop 'error in spec main: input_mod not specified'
      end select
!
      do indx_xobs=1, nxobs_fs

!          call print_ray(xobs(indx_xobs), zeta(indx_zeta), p(indx_p), iin, iin_c, 'TRASH/testa.dat')
!
         !calculation of profile function for given xobs along ray
         do indx_zray=1, nz_ray
            call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs(indx_xobs), phinorm)
            profile_ray(indx_zray)=phinorm
!             write(*,'(10es20.8)') velz_ray(indx_zray), profile_ray(indx_zray), vth_ray(indx_zray), vth_fiducial
         enddo
!          stop 'go on bla'
!
!         !setting intensity from core
         call get_iin(xobs(indx_xobs), lcore, iin, iin_c)
!          write(*,*) iin, iin_c, xic1
!          if(omp_get_thread_num().eq.1.and.indx_xobs.eq.1) write(*,'(i5,l5,5es20.8)') indx_xobs, lcore, iin, iin_c, p(indx_p)
!
         !formal solution along ray
         call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem, iem_c, iemi, iabs)
!         if(iem.gt.iem_c) then
!            do i=1, nz_ray
!               write(*,*) z_ray(i), opalbar_ray(i), opac_ray(i), sline_ray(i), scont_ray(i), profile_ray(i)
!            enddo
!            write(*,*) iem_c, iem, iin, iin_c
!            stop
!         endif
!
!x-integration
         ftot_x(indx_xobs) = ftot_x(indx_xobs) + iem * xpw(indx_xp)
         fcont_x(indx_xobs) = fcont_x(indx_xobs) + iem_c * xpw(indx_xp)
         normt_x(indx_xobs) = normt_x(indx_xobs) + xpw(indx_xp)
         femi_x(indx_xobs) = femi_x(indx_xobs) + iemi * xpw(indx_xp)
         fabs_x(indx_xobs) = fabs_x(indx_xobs) + iabs * xpw(indx_xp)
!corresponding errors
         ftot_errx(indx_xobs) = ftot_errx(indx_xobs) + iem * xpw_err(indx_xp)
         fcont_errx(indx_xobs) = fcont_errx(indx_xobs) + iem_c * xpw_err(indx_xp)
!
!debug open
!               sum=sum+iem_c * p(indx_p) * pw(indx_p)
!               write(*,*) ftot_p(indx_xobs), sum, p(indx_p)
!debug close
!
       !given p-ray is calculated for all xobs
       enddo
   !all xp-points are calculated for all xobs
   enddo
!
!output error p-integration
   call calc_err(nxobs_fs, ftot_x, fcont_x, ftot_errx, fcont_errx, err1, err2)
   write(*,'(a60, 3es20.8)') 'max rel error in xp-integration (tot, cont) for given yp', yp(indx_yp), err1, err2
   relerr_totx=max(err1, relerr_totx)
   relerr_contx=max(err2, relerr_contx)
!
!yp-integration
    flux_tot_tmp=flux_tot_tmp + ypw(indx_yp)*ftot_x
    flux_cont_tmp=flux_cont_tmp + ypw(indx_yp)*fcont_x
    normt_tmp=normt_tmp + ypw(indx_yp) * normt_x
    ftot_err_tmp=ftot_err_tmp + ypw_err(indx_yp) * ftot_errx
    fcont_err_tmp=fcont_err_tmp + ypw_err(indx_yp) * fcont_errx
    flux_emi_tmp=flux_emi_tmp + ypw(indx_yp)*femi_x
    flux_abs_tmp=flux_abs_tmp + ypw(indx_yp)*fabs_x
!
   !all yp-points are calculated for all xp points and all xobs
   enddo
!
!!now: adding up all thread-private arrays
!$omp critical
flux_tot=flux_tot + flux_tot_tmp
flux_cont=flux_cont + flux_cont_tmp
normt=normt + normt_tmp
flux_emi=flux_emi + flux_emi_tmp
flux_abs=flux_abs + flux_abs_tmp
!!corresponding errors
ftot_err=ftot_err + ftot_err_tmp
fcont_err=fcont_err + fcont_err_tmp
!$omp end critical
!
!deallocate all arrays that have been allocated in parallel region
call deallocate_fluxes
call deallocate_fs1d

!$omp end parallel
!
!output error yp-integration
call calc_err(nxobs_fs, flux_tot, flux_cont, ftot_err, fcont_err, relerr_tot, relerr_cont)
write(*,*)
write(*,'(a50, 2es20.8)') 'max rel error in yp-integration (tot, cont)', relerr_tot, relerr_cont
write(*,*)
!
!
end subroutine calc_fluxem
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_iem_surface(ixobs)
!
use prog_type
use options_surfb, only: input_mod
use dime_surfb, only: nxp, nyp
use mod_surfb, only: nz_ray, vth_ray, velz_ray, z_ray, temp_ray, &
                     profile_ray, nz_ray, opalbar_ray, opac_ray, sline_ray, scont_ray
use params_surfb, only: vth_fiducial
use mod_surfb, only: xp, yp, iem_surface, iemi_surface, iabs_surface, xobs_surface, phinorm, icont_surface
use mod_math, only: calc_phinorm
!!
implicit none
!!
!! ... arguments
integer(i4b), intent(in) :: ixobs
!!
!! ... local scalars
integer(i4b) :: ixp, iyp, indx_zray, err
real(dp) :: xobs
real(dp) :: iin, iin_c, iem, iem_c, iemi, iabs
!!
!! ... local arrays
!!
!! ... local logicals
logical :: lcore
!!
!! ... local functions
!real(dp) :: interpol_yp
!!
!!-----------------------------------------------------------------------
!!
!
xobs = xobs_surface(ixobs)
!
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------calculating emergent intensity on surface-------------------'
write(*,*) 'for xobs', xobs
write(*,*)
!
!-----------------------------------------------------------------------
!
if(allocated(iem_surface)) deallocate(iem_surface)
if(allocated(iemi_surface)) deallocate(iemi_surface)
if(allocated(iabs_surface)) deallocate(iabs_surface)
if(allocated(icont_surface)) deallocate(icont_surface)
!
allocate(iem_surface(nxp, nyp), stat=err)
allocate(iemi_surface(nxp, nyp), stat=err)
allocate(iabs_surface(nxp, nyp), stat=err)
allocate(icont_surface(nxp, nyp), stat=err)
!
!-----------------------------------------------------------------------
!
do ixp=1,nxp
   do iyp=1, nyp
!
      select case(input_mod)
!         case(0)
!            call setup_ray1d(zeta(indx_zeta), p(indx_p), lcore)
!         case(1)
!!            call setup_ray3d(zeta(indx_zeta), p(indx_p), lcore)
         !            call setup_ray3d_rot(zeta(indx_zeta), p(indx_p), lcore)
         case(2)
            call setup_ray3d(xp(ixp), yp(iyp), lcore)
         case default
            stop 'error in calc_iem_surface: input_mod not surfbified'
      end select
!
!------------------------calculation of indx_xobs1----------------------
!
!calculation of profile function for given xobs along ray
      do indx_zray=1, nz_ray
         call calc_phinorm(velz_ray(indx_zray), vth_ray(indx_zray), vth_fiducial, xobs, phinorm)
         profile_ray(indx_zray)=phinorm
      enddo
!
!setting intensity from core
      call get_iin(xobs, lcore, iin, iin_c)
!
!!formal solution along ray
      call formal_ray(z_ray, opalbar_ray, opac_ray, sline_ray, scont_ray, profile_ray, nz_ray, iin, iin_c, iem, iem_c, iemi, iabs)
!      write(*,*) z_ray
!      write(*,'(es20.8, l5, 3es20.8)') p(indx_p), lcore, iin, iin_c, iem1
!
!------------------interpolation of both xobs-solutions-----------------
!
      iem_surface(ixp,iyp) = iem
      iemi_surface(ixp,iyp) = iemi
      iabs_surface(ixp,iyp) = iabs
      icont_surface(ixp,iyp) = iem_c

!      write(*,*) iem/iem_c, iemi, iabs, iem_c
!
   enddo
enddo
!!
!!
!!
end subroutine calc_iem_surface
!
