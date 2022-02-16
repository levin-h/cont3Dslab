function bvel(rad, vinf, bconst, beta)
!
! calculates velocities from beta velocity law in cgs
!           v(rad) = vinf*(1.-bconst/rad)^beta
!
!INPUT: vinf in cgs
!       rad  in r_star
!       beta, bconst
!OUTPUT: v(rad) in cgs
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rad, vinf, bconst, beta
real(dp) :: bvel
!
bvel=vinf*(1.d0-bconst/rad)**beta
!
end function bvel
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine bvel3d(vmin, vinf, beta, x, y, z, velx, vely, velz, gradv)
!
!calculates velocity components of beta velocity law
!for test (spherical symmetric) purposes
!
!INPUT: vinf, vmin in cgs
!       x, y, z coordinates in r_star
!       rad  in r_star
!       beta
!OUTPUT: velx, vely, velz in cgs
!        gradv in cgs
!
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vmin, vinf, beta, x, y, z
real(dp), intent(out) :: velx, vely, velz, gradv
!
! ... local scalars
real(dp) :: rad, velr, bconst
!
bconst=1.d0-(vmin/vinf)**(1.d0/beta)
rad=sqrt(x**2 + y**2 + z**2)
!
velr=vinf*(1.d0-bconst/rad)**beta
!
velx=velr*x/rad
vely=velr*y/rad
velz=velr*z/rad
!
gradv = vinf*beta*(1.d0-bconst/rad)**(beta-1)*bconst/rad**2
!
end subroutine bvel3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine vinf_bc(vmin, vrot, vesc, vcrit, beta, zeta, gamma, theta0, vinf, b)
!
!   calculates terminal velocity and b-factors for
!         bjorkman&casinelli (1993) model
!
!INPUT: vmin:   minimum velocity at stellar surface
!       vrot:   rotational velocity
!       vesc:   escape velocity
!       vcrit:  break up velocity
!       beta:   beta parameter (of beta velocity law)
!       zeta:   zeta parameter
!       gamma:  gamma parameter
!       theta0: co-latitude of plane
! 
!OUTPUT: vinf:   terminal velocity at theta0
!        b:      b-coefficient of radial velocity law at theta0
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vmin, vrot, vesc, vcrit, &
                        beta, zeta, gamma, theta0
real(dp), intent(out) :: vinf, b
!
! ... local scalars
!
vinf = zeta*vesc*(1.d0-sin(theta0)*vrot/vcrit)**gamma
b = 1.d0-(vmin/vinf)**(1.d0/beta)
!
end subroutine vinf_bc
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function vthermal(vmicro, temp, matom)
!
!   calculates thermal velocity
!
!INPUT: vmicro: micro-turbulent velocity in cgs
!       temp: temperature in Kelvin
!       matom: atomic number
!

use prog_type
use fund_const, ONLY: cgs_kb, cgs_mp
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vmicro, temp
integer(i4b), intent(in) :: matom
real(dp) :: vthermal
!
vthermal = sqrt(2.d0*cgs_kb*temp/(matom*cgs_mp) + vmicro**2)
!
end function vthermal
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function deldop(xnue0, vth)
!
!   calculates thermal doppler width
!
!INPUT: xnue0: transition frequency in 1/s
!       vth: thermal velocity in cm/s
!
use prog_type
use fund_const, ONLY: cgs_clight
!
implicit none
!
! ... arguments
real(dp), intent(in) :: xnue0, vth
real(dp) :: deldop
!
deldop = xnue0 * vth / cgs_clight
!
end function deldop
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function calc_vmicro(vmin, vmax, vel)
!
!   calculates microt-turbulent velocity as a linear
!      function of absolute value of velocity
!
!INPUT:
!   vmin   minimum abs(velocity) = vmicro_min in cgs
!   vmax   maximum abs(velocity) = 10.*vmicro_max in cgs
!   vel    abs(velocity) at current position in cgs
!
USE prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: vmin, vmax, vel
real(dp) :: calc_vmicro
!
! ... local scalars
real(dp) :: vmicro_min, vmicro_max
!
!set minimum and maximum microturbulent velocities
vmicro_min=vmin
vmicro_max=0.1d0*vmax
!
calc_vmicro = vmicro_min + (vmicro_min-vmicro_max)*(vel-vmin)/(vmin-vmax)
!
!
end function calc_vmicro
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opac_thomson(yhe, hei, rho, kcont)
!
use prog_type
use fund_const, only: cgs_mp, sigmae
!
implicit none
!
! ... arguments
real(dp), intent(in) :: yhe, hei, rho, kcont
real(dp) :: opac_thomson
!
! ... local scalars
real(dp) :: c1, c2, ne
!
c1=(1.d0+4.d0*yhe)*cgs_mp
c2=(1.d0+hei*yhe)/c1
ne=c2*rho
opac_thomson=sigmae*ne*kcont
!
end function opac_thomson
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opac_thomson2(sr, vinf, mdot, yhe, hei, rho, kcont, alpha, rad)
!
!same as opac_thomson, however, opacity multiplied by additional factor (v(r)/vinf)**alpha
!
!  INPUT: sr:      stellar radius in cgs
!         vinf:    terminal velocity in cgs
!         mdot:    mass-loss rate in cgs
!         rad:     radial coordinate in cgs
!         rho:     density in cgs
!  OUTPUT: frequency integrated line opacity in frequency space (in 1/cm)
!
use prog_type
use fund_const, only: pi, cgs_mp, sigmae
!
implicit none
!
! ... arguments
real(dp), intent(in) :: sr, vinf, mdot, yhe, hei, rho, kcont, alpha, rad
real(dp) :: opac_thomson2
!
! ... local scalars
real(dp) :: vr, c1, c2, ne
!
! ... local arrays
!
!-----------------------------------------------------------------------
!
vr = mdot/4.d0/pi/rad/rad/rho
!
c1=(1.d0+4.d0*yhe)*cgs_mp
c2=(1.d0+hei*yhe)/c1
ne=c2*rho
opac_thomson2=sigmae*ne*kcont*(vr/vinf)**alpha


!open(1, file='outputFILES_DEBUG/data_test.dat', access='append')
!   write(1,*), rad/sr, vr/vinf, opac_thomson2
!close(1)
!
end function opac_thomson2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opac_opal(kcont, yhe, hei, rho, temp, nrho, ntemp, rho_opal, temp_opal, kappa_opal)
!
!all temperatures and densities in log-space
!kappa-table in log-space
!
use prog_type
use mod_interp1d, only: find_index  
use mod_interp2d, only: interpol2d_4p_lin
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nrho, ntemp
real(dp), intent(in) :: kcont, rho, temp, yhe, hei
real(dp), dimension(nrho), intent(in) :: rho_opal
real(dp), dimension(ntemp), intent(in) :: temp_opal
real(dp), dimension(ntemp,nrho) :: kappa_opal
real(dp) :: opac_opal
!
! ... local scalars
integer(i4b) :: i, iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1
real(dp) :: kappa
real(dp) :: rho_local, t_local
!
! ... local functions
real(dp) :: opac_thomson
!
!
t_local=10.**temp
rho_local = 10.d0**rho/(t_local*1.d-6)**3
rho_local = log10(rho_local)
t_local=temp
!
!
!use thomson opacity if temperatures too large
if(t_local.gt.maxval(temp_opal)) then
   opac_opal = opac_thomson(yhe, hei, 10.d0**rho, kcont)
   write(*,*) opac_opal/10.d0**rho
else
!
!adapt rho-local values if outside of opacity table
   do i=1, 1000
      if(rho_local.lt.minval(rho_opal)) then
         rho_local = rho_local + 0.5
      elseif(rho_local.gt.maxval(rho_opal)) then
         rho_local = rho_local - 0.5
      else
         exit
      endif
   enddo
!
!adapt t_local values if outside of opacity table
   do i=1, 1000
      if(t_local.lt.minval(temp_opal)) then
         t_local = t_local + 0.05
!   elseif(t_local.gt.maxval(temp_opal)) then
!      t_local = t_local - 0.05
      else
         exit
      endif
   enddo
!
!
!   
   call find_index(temp, temp_opal, ntemp, iim2, iim1, ii, iip1)
   call find_index(rho_local, rho_opal, nrho, jjm2, jjm1, jj, jjp1)

   kappa = interpol2d_4p_lin(kappa_opal(iim1, jjm1), kappa_opal(ii,jjm1), &
                             kappa_opal(iim1, jj),   kappa_opal(ii,jj), &
                             temp_opal(iim1), temp_opal(ii), &
                             rho_opal(jjm1), rho_opal(jj), temp, rho_local)
!
!
!write(*,*) rho, temp, rho_local, rho_local_orig, t_local, kappa, 10.d0**kappa
!
   opac_opal = kcont * 10.d0**kappa * 10.d0**rho
!
endif  
!
end function opac_opal
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbar_model_kline(yhe, hei, rho, kline)
!
!  calculates line opacity according to line-strength parameter
!          (see definition in puls & springmann 2000)
!
!  input: yhe   helium abundance by number
!         hei   number of free electrons per helium atom
!         rho   density in cgs
!         kline  line-strength parameter
!  output: frequency integrated line opacity in frequency-space (in 1/cm)
!
use prog_type
use fund_const, only: sigmae, cgs_mp
!
implicit none
!
! ... arguments
real(dp), intent(in) :: yhe, hei, rho, kline
real(dp) :: opalbar_model_kline
!
! ... local scalars
real(dp) :: c1, c2, ne, chi_thomson
!
! ... local arrays
!
!calculate thomson-opacity
c1=(1.d0+4.d0*yhe)*cgs_mp
c2=(1.d0+hei*yhe)/c1
ne=c2*rho
chi_thomson=sigmae*ne
!
!calculate line opacity
opalbar_model_kline = kline*chi_thomson
!
end function opalbar_model_kline
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbar_model_kline2(yhe, hei, rho, kline)
!
!  opalbar = kline*rho^2
!
!  input: yhe   helium abundance by number
!         hei   number of free electrons per helium atom
!         rho   density in cgs
!         kline  line-strength parameter
!  output: frequency integrated line opacity in frequency-space (in 1/cm)
!
use prog_type
use fund_const, only: sigmae, cgs_mp
!
implicit none
!
! ... arguments
real(dp), intent(in) :: yhe, hei, rho, kline
real(dp) :: opalbar_model_kline2
!
! ... local scalars
real(dp) :: c1, c2, ne, chi_thomson
!
! ... local arrays
!
!calculate thomson-opacity
c1=(1.d0+4.d0*yhe)*cgs_mp
c2=(1.d0+hei*yhe)/c1
ne=c2*rho
chi_thomson=sigmae*ne
!
!calculate line opacity
opalbar_model_kline2 = kline*rho**2
!
end function opalbar_model_kline2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbar_model_hamann(sr, vinf, mdot, kappa0, alpha, vth_fid, rad, rho)
!
!  calculates line opacity according to definition by hamann 1981
!(see also Puls,Owocki,Fullerton 1993, Sundqvist Diss (Appendix))
!
!  INPUT: sr:      stellar radius in cgs
!         vinf:    terminal velocity in cgs
!         mdot:    mass-loss rate in cgs
!         kappa0:  line-strength parameter (dimensionless) as defined by Hamann
!         alpha:   see Hamann
!         vth_fid: fiducial thermal velocity in cgs
!         rad:     radial coordinate in cgs
!         rho:     density in cgs
!  OUTPUT: frequency integrated line opacity in frequency space (in 1/cm)
!
use prog_type
use fund_const, only: pi
!
implicit none
!
! ... arguments
real(dp), intent(in) :: sr, vinf, mdot, kappa0, vth_fid, rho, rad, alpha
real(dp) :: opalbar_model_hamann
!
! ... local scalars
!
! ... local arrays
!
!-----------------------------------------------------------------------
!
!old version: without alpha
!in own units
!opalbar_model_hamann = 4.d0*pi*sr*vinf*vinf*kappa0*rho/mdot/vth_fid
!
!new version: including alpha: opalbar=kappa0*v(r)^alpha / r^2 / v(r) (with v in vmax, r in rstar)
!in own units
opalbar_model_hamann = 4.d0*pi*sr*vinf*vinf*kappa0*rho/mdot/vth_fid * (mdot/4.d0/pi/rad/rad/rho/vinf)**alpha
!
!
end function opalbar_model_hamann

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbar_halpha(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
!
!assuming complete H-ionization and recombination line  
!all input in cgs
!
!-----------------see personal notes on opacities-----------------------
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rho, temp, b2, b3, yhe, hei, sr, vth_fiducial, xnue0
real(dp) :: opalbar_halpha
!
! ... local scalars
real(dp), parameter :: energy_leveln=cgs_clight*cgs_planck*109678.77d0, &
                       energy_level2=cgs_clight*cgs_planck*82259.158d0, &
                       energy_level3=cgs_clight*cgs_planck*97492.305d0, &
                       c1 = half*(cgs_planck**2/two/pi/cgs_me/cgs_kb)**(three/two), &
                       c2 = pi*cgs_e**2 / cgs_me/cgs_clight /cgs_mp**2, &
                       gf=5.1286d0
real(dp) ::deldop_fiducial, dum1
!
! ... local arrays
!
dum1=(one+yhe*hei)/(one+four*yhe)**2
dum1=dum1*rho*rho*temp**(-three/two)
dum1=dum1*c1*c2*gf
dum1=dum1*(b2*exp((energy_leveln-energy_level2)/cgs_kb/temp) - &
           b3*exp((energy_leveln-energy_level3)/cgs_kb/temp))
!
!write(*,*) rho, temp, yhe, hei, dum1
!so far in cgs, now in own units
!
!finally need to divide this by deldop_fiducial
deldop_fiducial=xnue0*vth_fiducial/cgs_clight
opalbar_halpha=dum1*sr/deldop_fiducial
!
!
end function opalbar_halpha
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbar_hbeta(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b4, rho)
!
!assuming complete H-ioniazation and recombination line  
!all input in cgs
!
!-----------------see personal notes on opacities-----------------------
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rho, temp, b2, b4, yhe, hei, sr, vth_fiducial, xnue0
real(dp) :: opalbar_hbeta
!
! ... local scalars
real(dp), parameter :: energy_leveln=cgs_clight*cgs_planck*109678.77d0, &
                       energy_level2=cgs_clight*cgs_planck*82259.158d0, &
                       energy_level4=cgs_clight*cgs_planck*102832.904d0, &
                       c1 = half*(cgs_planck**2/two/pi/cgs_me/cgs_kb)**(three/two), &
                       c2 = pi*cgs_e**2 / cgs_me/cgs_clight /cgs_mp**2, &
                       gf = 0.95508055
real(dp) ::deldop_fiducial, dum1
!
! ... local arrays
!
dum1=(one+yhe*hei)/(one+four*yhe)**2
dum1=dum1*rho*rho*temp**(-three/two)
dum1=dum1*c1*c2*gf
dum1=dum1*(b2*exp((energy_leveln-energy_level2)/cgs_kb/temp) - &
           b4*exp((energy_leveln-energy_level4)/cgs_kb/temp))
!
!write(*,*) rho, temp, yhe, hei, dum1
!so far in cgs, now in own units
!
!finally need to divide this by deldop_fiducial
deldop_fiducial=xnue0*vth_fiducial/cgs_clight
opalbar_hbeta=dum1*sr/deldop_fiducial
!
!
end function opalbar_hbeta
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbar_halpha2(sr, yhe, temp, vth_fiducial, xnue0, b2, b3, rho)
!
!assuming neutral H, and all elements in ground state  
!all input in cgs
!
!-----------------see personal notes on opacities-----------------------
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rho, temp, b2, b3, yhe, sr, vth_fiducial, xnue0
real(dp) :: opalbar_halpha2
!
! ... local scalars
real(dp), parameter :: energy_level2=cgs_clight*cgs_planck*82259.158d0, &
                       energy_level3=cgs_clight*cgs_planck*97492.305d0, &
                       f23 = 6.4108d-1, &
                       g1 = 2.d0, g2 = 8.d0, &
                       c1 = pi*cgs_e**2 / cgs_me/cgs_clight /cgs_mp
real(dp) ::deldop_fiducial, dum1
!
! ... local arrays
!
dum1 = c1 * f23 * g2/g1 * rho/(one+four*yhe)
dum1=dum1*(b2*exp(-energy_level2/cgs_kb/temp) - &
           b3*exp(-energy_level3/cgs_kb/temp))
!write(*,*) rho, temp, yhe, hei, dum1
!so far in cgs, now in own units
!
!
!finally need to divide this by deldop_fiducial
deldop_fiducial=xnue0*vth_fiducial/cgs_clight
opalbar_halpha2=dum1*sr/deldop_fiducial
!
!
end function opalbar_halpha2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbar_halpha3(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
!
!calculate hydrogen ionization balance from LTE
!calculate level populations with departure coefficients
!all input in cgs
!
!-----------------see personal notes on opacities-----------------------
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rho, temp, b2, b3, yhe, hei, sr, vth_fiducial, xnue0
real(dp) :: opalbar_halpha3
!
! ... atomic data (from NIST)
integer(i4b), parameter :: nlist = 40
real(dp), dimension(nlist), parameter :: gvalue= (/ 2.d0, 8.d0, 18.d0, 32.d0, 50.d0, 72.d0, &
     98.d0, 128.d0, 162.d0, 200.d0, 242.d0, 288.d0, 338.d0, 392.d0, 450.d0, 512.d0, 578.d0, &
     648.d0, 722.d0, 800.d0, 882.d0, 968.d0, 1058.d0, 1152.d0, 1250.d0, 1352.d0, 1458.d0, &
     1568.d0, 1682.d0, 1800.d0, 1922.d0, 2048.d0, 2178.d0, 2312.d0, 2450.d0, 2592.d0, 2738.d0, &
     2888.d0, 3042.d0, 3200d0 /)
real(dp), dimension(nlist), parameter :: energy= (/ 0.d0, 82259.158d0, 97492.304d0, 102823.904d0, &
     105291.657d0, 106632.1681d0, 107440.4508d0, 107965.0568d0, 108324.7253d0, 108581.9945d0, &
     108772.3445d0, 108917.1209d0, 109029.7913d0, 109119.1917d0, 109191.3154d0, 109250.3433d0, &
     109299.2642d0, 109340.2604d0, 109374.9555d0, 109404.5776d0, 109430.0696d0, 109452.1650d0, &
     109471.4416d0, 109488.3592d0, 109503.2875d0, 109516.5267d0, 109528.3223d0, 109538.8768d0, &
     109548.3584d0, 109556.9077d0, 109564.6431d0, 109571.6647d0, 109578.0577d0, 109583.8949d0, &
     109589.2390d0, 109594.1439d0, 109598.6566d0, 109602.8177d0, 109606.6628d0, 109610.2232d0 /) * cgs_clight*cgs_planck
real(dp), parameter :: energy_ion = 109678.77d0 * cgs_clight*cgs_planck
real(dp), parameter :: f23 = 6.4108d-1, &
                       c1 = half*(cgs_planck**2/two/pi/cgs_me/cgs_kb)**(three/two), &
                       c2 = pi*cgs_e**2 / cgs_me/cgs_clight

! ... local scalars
integer(i4b) :: i
real(dp) :: energy_level2, energy_level3, g2, g3
real(dp) :: deldop_fiducial, fdum1, fdum2, dum1
real(dp) :: nh, nhi, nhii, n2, n3, nel, zhi, zhii, q, disc
!
! ... local arrays
!
!hydrogen density
nh = rho/(one+four*yhe)/cgs_mp
!
!partition functions
zhi = zero
zhii = one
do i=1, nlist
   zhi = zhi + gvalue(i) * exp(-energy(i)/cgs_kb/temp)
enddo
!
!calculate q value (see personal notes)
q = temp**(3./2.)/c1 * zhii/zhi * exp(-energy_ion/cgs_kb/temp)
!
!calculate electron density
fdum1 = q - hei*yhe*nh
fdum2 = -q*nh*(one+hei*yhe)
disc = fdum1**2 - four*fdum2
if(disc.lt.zero) stop 'error in opalbar_halpha3: cannot calculate electron density with discriminant < zero'
nel = (-fdum1 + sqrt(disc))/two

!write(*,*) temp, nh, q, disc, nel, zhii, zhi
!
!calculate hydrogen ionization stages
nhii = nel - hei*yhe*nh
nhi = nh - nhii
!
!calculate level populations of level two and three
energy_level2 = energy(2)
energy_level3 = energy(3)
g2 = gvalue(2)
g3 = gvalue(3)
n2 = b2 * g2 * nhi/zhi * exp(-energy_level2/cgs_kb/temp)
n3 = b3 * g3 * nhi/zhi * exp(-energy_level3/cgs_kb/temp)

dum1 = c2*f23*(n2 - g2/g3 * n3)


!write(*,'(10es20.8)') temp, nel, nh, nhi, nhii, n2, n3, dum1
!write(*,*) rho, temp, yhe, hei, dum1
!so far in cgs, now in own units
!
!
!finally need to divide this by deldop_fiducial
deldop_fiducial=xnue0*vth_fiducial/cgs_clight
opalbar_halpha3=dum1*sr/deldop_fiducial
!
if(opalbar_halpha3.lt.zero) stop 'error in opalbar_halpha3: opacity < zero not allowed'
!
end function opalbar_halpha3
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE COMPONENT_W(rad, mu, rho, temp, vel, vel_r, vel_theta, &
                       teff, v_esc, v_inf, rhoc_star)
!
!           CALCULATES DENSITY, TEMPERATURE AND VELOCITY
!                  FOR WIND UPFLOW COMPONENT
!
!INPUT: rad:   radius where values shall be calculated
!       mu:    co-latitude (w.r.t magnetic pole axis) where values shall be calculated
!       teff, v_esc, v_inf, rhoc_star: stellar parameter
!
!OUTPUT: density, temperature, velocity and velocity components
!        at those points
!
!NOTE: all quantities are defined in a carthesian coordinate system, where the z-axis
!      is aligned with the magnetic-pole-axis
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: rad, mu, teff, v_esc, v_inf, rhoc_star
REAL(DP) :: rho, temp, vel, vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: dum1, dum2, dum3, dum4
!
!----------------------------density------------------------------------
!
dum1 = sqrt(rad-1.d0+mu*mu)*sqrt(1.d0+3.d0*mu*mu)
dum2 = (1.d0/rad)**(3.d0/2.d0)
dum3 = rad-1.d0
dum4 = 4.d0*rad-3.d0+3.d0*mu*mu
!
rho = 2.d0*(rhoc_star*v_esc/v_inf)*dum1*dum2/dum3/dum4
!
!--------------------------temperature----------------------------------
!
temp = TEFF
!
!----------------------------velocity-----------------------------------
!
vel = v_inf*(1.d0-1.d0/rad)
!
if(mu.ge.0.d0) then
   dum1 = 2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
else
   dum1 = -2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = -1.d0*sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
endif
!
!
vel_r = vel*dum1
vel_theta = vel*dum2
!
!
!
END SUBROUTINE COMPONENT_W
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE COMPONENT_C(rad, mu, r_apex, rho, temp, vel, vel_r, vel_theta, &
                       teff, v_esc, rhoc_star, delta)
!
!           CALCULATES DENSITY, TEMPERATURE AND VELOCITY
!                  FOR COOL DOWNFLOW COMPONENT
!
!INPUT: rad:    radius where values shall be calculated
!       mu:     co-latitude where values shall be calculated
!       r_apex: radius of loop apex
!
!OUTPUT: rho:    density at the given point
!        temp:   temperature at the given pont
!        vel, vel_r, vel_theta: velocity components and absoulute velocity
!
!NOTE: all quantities are defined in carthensian coordinate system, 
!      where z-axis is aligned with magnetic pole axis
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: rad, mu, r_apex, teff, v_esc, rhoc_star, delta
REAL(DP) :: rho, temp, vel, vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: dum1, dum2, dum3, dum4, mu_dum
!
!----------------------------density------------------------------------
!
mu_dum=mu
if(mu_dum.eq.0.d0) then
!avoid division by zero
   mu_dum=1.d-13
endif
!
dum1 = sqrt(rad-1.d0+mu_dum*mu_dum)*sqrt(1.d0+3.d0*mu_dum*mu_dum)
dum2 = 1.d0/rad/rad
dum3 = sqrt(mu_dum*mu_dum+delta*delta/rad/rad)
dum4 = 4.d0*rad-3.d0+3.d0*mu_dum*mu_dum
rho = 2.d0*rhoc_star*dum1*dum2/dum3/dum4
!
!--------------------------temperature----------------------------------
!
temp = TEFF
!
!----------------------------velocity-----------------------------------
!
!write(*,*) rad, r_apex
if(abs(rad-r_apex).lt.1.d-14) then
   vel = 0.d0
else
   vel = v_esc * sqrt(1.d0/rad - 1.d0/r_apex)
endif
!
if(mu.ge.0.d0) then
   dum1 = -2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = -1.d0*sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
else
   dum1 = 2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
endif
!
vel_r = vel*dum1
vel_theta = vel*dum2
!
!
!
END SUBROUTINE COMPONENT_C
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE COMPONENT_S(rad, mu, mu_star, mu_shock, rho, temp, vel, vel_r, vel_theta, &
                       teff, t_inf, v_esc, v_inf, rhoc_star, delta)
!
!           CALCULATES DENSITY, TEMPERATURE AND VELOCITY
!                FOR POST-SHOCK COMPONENT
!
!INPUT: rad: radius where values shall be calculated
!       mu:  co-latitude where values shall be calculated
!       mu_star: intersection point of photosphere with closed loop (attached to given point)
!       mu_shock: corresponding to the shock-radius on the closed loop
!
!OUTPUT: rho: density
!        temp: temperature
!        vel, vel_r, vel_theta: velocity-components and absolute velocity
!
!NOTE: all quantities are defined in carthensian coordinate system, 
!      where z-axis is aligned with magnetic pole axis
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: rad, mu, mu_star, mu_shock, teff, t_inf, v_esc, v_inf, rhoc_star, delta
REAL(DP) :: rho, temp, vel, vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: dum0, dum1, dum2, dum3, dum4, dum_temp, rho_w
!
!--------------------------temperature----------------------------------
!
dum1 = 1.d0 - (1.d0-mu_star*mu_star)/(1.d0-mu_shock*mu_shock)
dum1 = dum1*dum1
dum2 = abs(mu - mu**3.d0 + 3.d0*mu**5.d0/5.d0 - mu**7.d0/7.d0)
dum3 = abs(mu_shock - mu_shock**3.d0 + 3.d0*mu_shock*5.d0/5.d0 - mu_shock**7.d0/7.d0)
dum4 = (dum2/dum3)**(1.d0/3.d0)
!
dum_temp= T_INF*dum1
temp = dum_temp*dum4
temp = max(temp, TEFF)
!
!----------------------------density------------------------------------
!
!calculate wind component at r_shock, mu_shock
dum0 = (1.d0-mu_shock*mu_shock)/(1.d0-mu_star*mu_star)
dum1 = sqrt(dum0-1.d0+mu_shock*mu_shock)*sqrt(1.d0+3.d0*mu_shock*mu_shock)
dum2 = (1.d0/dum0)**(3.d0/2.d0)
dum3 = dum0-1.d0
dum4 = 4.d0*dum0-3.d0+3.d0*mu*mu
rho_w = (rhoc_star*v_esc/v_inf)*dum1*dum2/dum3/dum4
!
!calculate density in post-shock region
rho = 4.d0 * rho_w * dum_temp / temp
!
!----------------------------velocity-----------------------------------
!
dum1 = (1.d0-(1.d0-mu_star*mu_star)/(1.d0-mu_shock*mu_shock))*V_INF/4.d0
dum2 = temp/dum_temp
dum3 = (1.d0-mu_shock*mu_shock)/(1.d0-mu*mu)
dum3 = dum3**3.d0
dum4 = sqrt(1.d0+3.d0*mu*mu)/sqrt(1.d0+3.d0*mu_shock*mu_shock)
vel = dum1*dum2*dum3*dum4
!
if(mu.ge.0.d0) then
   dum1 = 2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
else
   dum1 = -2.d0*mu/sqrt(1.d0+3.d0*mu*mu)
   dum2 = -1.d0*sqrt(1.d0-mu*mu)/sqrt(1.d0+3.d0*mu*mu)
endif
!
vel_r = vel*dum1
vel_theta = vel*dum2
!
!
END SUBROUTINE COMPONENT_S
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE VEL_ADM_WCS(rad, mu, v_inf, v_esc, ralfven, chi_inf, t_inf, teff, vel_r, vel_theta)
!
!--------------------CALCULATE ADM VELOCITY COMPONENTS------------------
!
!INPUT: rad      current radial position in r_star mu current latitudinal position
!                (measured w.r.t magnetic pole axis) 
!       v_inf    terminal velocity in cm/s
!       v_esc    escape velocity in cm/s
!       ralfven  alfven radius in r_star
!       chi_inf  cooling parameter
!       t_inf, teff
!
!OUTPUT: vel_r      velocity component in radial direction
!        vel_theta  velocity component in latitudinal direction
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: v_inf, v_esc, rad, mu, ralfven, chi_inf, teff, t_inf
REAL(DP) :: vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: mmu, mu_star, mu_lower, mu_upper, mu_shock, r_apex, r_shock
REAL(DP) :: dum1, dum2, vel, temp, temp_shock
!
!-----------------------------------------------------------------------
!
mmu=mu
if(mmu.eq.1.d0) mmu=0.9999999d0
if(mmu.eq.-1.d0) mmu=-0.9999999d0
!
!--------------------------wind upflow----------------------------------
!
vel = v_inf*(1.d0-1.d0/rad)
!
if(mmu.ge.0.d0) then
   dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
   dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
else
   dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
   dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
endif
!
vel_r = vel*dum1
vel_theta = vel*dum2
!
!------------------------cooled downflow--------------------------------
!
mu_star=sqrt(1.d0-(1.d0-mmu*mmu)/rad)
mu_star=min(mu_star, 0.999999999d0)               
r_apex=1.d0/(1.d0-mu_star*mu_star)
!
if(r_apex.lt.ralfven) then
     
   if(abs(rad-r_apex).lt.1.d-14) then
      vel = 0.d0
   else
      vel = v_esc * sqrt(1.d0/rad - 1.d0/r_apex)
   endif
!
   if(mmu.ge.0.d0) then
      dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   else
      dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   endif

   vel_r = vel*dum1
   vel_theta = vel*dum2
!
!-------------------------shock region----------------------------------
!
   mu_lower=0.d0
   mu_upper=mu_star
   CALL GET_MU_SHOCK(mu_star, chi_inf, mu_lower, mu_upper, mu_shock)
   r_shock=r_apex*(1.d0-mu_shock*mu_shock)

   if(rad.gt.r_shock) then

!temperature
      temp_shock = t_inf * (1.d0-1.d0/r_shock)**2
      temp = (abs(mmu - mmu**3 + 3.d0*mmu**5/5.d0 - mmu**7/7.d0) / &
              abs(mu_shock - mu_shock**3 + 3.d0*mu_shock*5/5.d0 - mu_shock**7/7.d0))**(1.d0/3.d0)
      temp = temp*temp_shock
      temp = max(temp, TEFF)

!velocity
      dum1 = (1.d0-1.d0/r_shock)*v_inf/4.d0 * temp/temp_shock
      vel = dum1 * (r_shock/rad)**3 * sqrt((1.d0+3.d0*mmu**2)/(1.d0+3.d0*mu_shock**2))
!
      if(mmu.ge.0.d0) then
         dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
         dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
      else
         dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
         dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
      endif
!
!      vel_r = vel*dum1
!      vel_theta = vel*dum2
   endif
!
endif


END SUBROUTINE VEL_ADM_WCS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE VEL_ADM_WCS_BVEL(rad, mu, v_inf, v_esc, ralfven, chi_inf, t_inf, teff, vel_r, vel_theta)
!
!--------------------CALCULATE ADM VELOCITY COMPONENTS------------------
!---NOTE: IN WIND UPFLOW REGION, BETA-VELOCITY LAW
!
!INPUT: rad      current radial position in r_star mu current latitudinal position
!                (measured w.r.t magnetic pole axis) 
!       v_inf    terminal velocity in cm/s
!       v_esc    escape velocity in cm/s
!       ralfven  alfven radius in r_star
!       chi_inf  cooling parameter
!       t_inf, teff
!
!OUTPUT: vel_r      velocity component in radial direction
!        vel_theta  velocity component in latitudinal direction
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: v_inf, v_esc, rad, mu, ralfven, chi_inf, teff, t_inf
REAL(DP) :: vel_r, vel_theta
!
! ... local scalars
REAL(DP) :: mmu, mu_star, mu_lower, mu_upper, mu_shock, r_apex, r_shock
REAL(DP) :: dum1, dum2, vel, temp, temp_shock
!
!-----------------------------------------------------------------------
!
mmu=mu
if(mmu.eq.1.d0) mmu=0.9999999d0
if(mmu.eq.-1.d0) mmu=-0.9999999d0
!
!--------------------------wind upflow----------------------------------
!
vel = v_inf*(1.d0-1.d0/rad)
!
vel_r = vel
vel_theta = 0.d0
!
!------------------------cooled downflow--------------------------------
!
mu_star=sqrt(1.d0-(1.d0-mmu*mmu)/rad)
mu_star=min(mu_star, 0.999999999d0)               
r_apex=1.d0/(1.d0-mu_star*mu_star)
!
if(r_apex.lt.ralfven) then
     
   if(abs(rad-r_apex).lt.1.d-14) then
      vel = 0.d0
   else
      vel = v_esc * sqrt(1.d0/rad - 1.d0/r_apex)
   endif
!
   if(mmu.ge.0.d0) then
      dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   else
      dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   endif

   vel_r = vel*dum1
   vel_theta = vel*dum2
!
!-------------------------shock region----------------------------------
!
   mu_lower=0.d0
   mu_upper=mu_star
   CALL GET_MU_SHOCK(mu_star, chi_inf, mu_lower, mu_upper, mu_shock)
   r_shock=r_apex*(1.d0-mu_shock*mu_shock)

   if(rad.gt.r_shock) then

!temperature
      temp_shock = t_inf * (1.d0-1.d0/r_shock)**2
      temp = (abs(mmu - mmu**3 + 3.d0*mmu**5/5.d0 - mmu**7/7.d0) / &
              abs(mu_shock - mu_shock**3 + 3.d0*mu_shock*5/5.d0 - mu_shock**7/7.d0))**(1.d0/3.d0)
      temp = temp*temp_shock
      temp = max(temp, TEFF)

!velocity
      dum1 = (1.d0-1.d0/r_shock)*v_inf/4.d0 * temp/temp_shock
      vel = dum1 * (r_shock/rad)**3 * sqrt((1.d0+3.d0*mmu**2)/(1.d0+3.d0*mu_shock**2))
!
      if(mmu.ge.0.d0) then
         dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
         dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
      else
         dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
         dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
      endif
!
!      vel_r = vel*dum1
!      vel_theta = vel*dum2
   endif
!
endif


END SUBROUTINE VEL_ADM_WCS_BVEL
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE RHOVEL_ADM_C(rad, mu, v_esc, ralfven, rhoc_star, delta, vel_r, vel_theta, rho)
!
!--------------------CALCULATE ADM VELOCITY COMPONENTS------------------
!
!INPUT: rad      current radial position in r_star mu current latitudinal position
!                (measured w.r.t magnetic pole axis) 
!       v_esc    escape velocity in cm/s
!       ralfven  alfven radius in r_star
!
!OUTPUT: vel_r      velocity component in radial direction
!        vel_theta  velocity component in latitudinal direction
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) ::  v_esc, rad, mu, ralfven, rhoc_star, delta
REAL(DP) :: vel_r, vel_theta, rho
!
! ... local scalars
REAL(DP) :: mmu, mu_star, r_apex
REAL(DP) :: dum1, dum2, dum3, dum4, vel
!
!------------------------cooled downflow--------------------------------
!
mmu=mu
if(mmu.eq.1.d0) mmu=0.9999999d0
if(mmu.eq.-1.d0) mmu=-0.9999999d0
!
mu_star=sqrt(1.d0-(1.d0-mmu*mmu)/rad)
mu_star=min(mu_star, 0.999999999d0)               
r_apex=1.d0/(1.d0-mu_star*mu_star)
!
!-----------------------------------------------------------------------
!
if(r_apex.ge.ralfven) then
   vel_r=0.d0
   vel_theta=0.d0
   rho=0.d0
   return
else
!
!----------------------velocity-----------------------------------------
!
   if(abs(rad-r_apex).lt.1.d-14) then
      vel = 0.d0
   else
      vel = v_esc * sqrt(1.d0/rad - 1.d0/r_apex)
   endif
!
   if(mmu.ge.0.d0) then
      dum1 = -2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = -1.d0*sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   else
      dum1 = 2.d0*mmu/sqrt(1.d0+3.d0*mmu*mmu)
      dum2 = sqrt(1.d0-mmu*mmu)/sqrt(1.d0+3.d0*mmu*mmu)
   endif
!
   vel_r = vel*dum1
   vel_theta = vel*dum2
!
!----------------------density------------------------------------------
!
   dum1 = sqrt(rad-1.d0+mmu**2)*sqrt(1.d0+3.d0*mmu**2)
   dum2 = 1.d0/rad**2
   dum3 = sqrt(mmu**2+(delta/rad)**2)
   dum4 = 4.d0*rad-3.d0+3.d0*mmu**2
   rho = 2.d0*rhoc_star*dum1*dum2/dum3/dum4
!
endif


END SUBROUTINE RHOVEL_ADM_C
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
SUBROUTINE GET_MU_SHOCK(mu_star, chi_inf, xl, xu, mu_shock)
!
!-----------------------------------------------------------------------
!
!   CALCULATE MU_SHOCK (CORRESPONDING TO SHOCK-RADIUS) FOR GIVEN 
!                       r_apex, chi_inf
!
!   A TRANSCENDENTAL FUNCTION (SEE Owocki et al 2016, eq. 19)
!         NEEDS TO BE SOLVED
!
!   THIS IS DONE BY APPLYING THE REGULA-FALSI METHOD
!
!ON INPUT: mu_star: INTERSECTION POINT OF DIPOLE-LOOP WITH PHOTOSPHERE
!          chi_inf: COOLING PARAMETER
!          xl: LOWER GUESS OF MU_SHOCK (WILL BE DESTROYED)
!          xu: UPPER GUESS OF MU_SHOCK (WILL BE DESTROYED)
!
!ON OUTPUT: mu_shock: SOLUTION OF THE PROBLEM
!
!-----------------------------------------------------------------------
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: mu_star, chi_inf
REAL(DP) :: xl, xu, mu_shock
!
!
! ... local scalars
INTEGER(I4B), PARAMETER :: nit=100
INTEGER(I4B) :: i
REAL(DP), PARAMETER :: eps=1.d-14
REAL(DP) :: x_lower, x_upper, y_lower, y_upper, x_new, y_new
REAL(DP) :: swap
!
! ... local functions
REAL(DP) :: g_fct
!
!-----------------------------------------------------------------------
!
!swap lower and upper values, if necessary
if(xl.gt.xu) then
   swap=xl
   xl=xu
   xu=swap
endif
!
x_lower=xl
x_upper=xu
y_lower = g_fct(x_lower, mu_star, chi_inf)
y_upper = g_fct(x_upper, mu_star, chi_inf)
!
if(y_lower*y_upper .gt. 0.d0) then
!error if no interval is found
   STOP 'ERROR IN GET_MU_SHOCK: NO INTERVAL FOR WHICH f(xl) < f(xu) IS FOUND'
endif
!
!--now, when start values give unequal signs, make regula falsi method--
!
do i=1, nit
!
   x_new = x_lower - y_lower * (x_upper-x_lower)/(y_upper-y_lower)
   y_new = g_fct(x_new, mu_star, chi_inf)
!
   if(abs(y_new).le.eps) exit
!
   if(y_new*y_lower.lt.0.d0 ) then
      x_upper=x_new
      y_upper=y_new
   else
      x_lower=x_new
      y_lower=y_new
   endif
!
enddo
!
mu_shock=x_new
!
!
!
END SUBROUTINE GET_MU_SHOCK
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
FUNCTION G_FCT(x, mu_star, chi_inf)
!
USE prog_type
!
IMPLICIT NONE
!
! ... arguments
REAL(DP), INTENT(IN) :: x, mu_star, chi_inf
REAL(DP) :: g_fct
!
! ... local scalars
REAL(DP) :: term1, term2, term3, term4, dum1, dum2
!
dum1=1.d0-x*x
dum2=1.d0-mu_star*mu_star
!
term1 = x - x**3.d0 + 3.d0/5.d0 * x**5.d0 - x**7.d0 / 7.d0
term2 = chi_inf*(1.d0+3.d0*mu_star**2.d0)/6.d0/mu_star/dum2**2.d0
term3 = (1.d0-dum2/dum1)**4.d0
term4 = dum1**6.d0 / (1.d0+3.d0*x**2.d0)
!
g_fct = abs(term1) - term2*term3*term4
!
!
!
END FUNCTION G_FCT

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine depcoeff_petrenz(velr, b2, b3)
!
!--calculate departure coefficients according to petrenz & puls 1995----
!-----------------from puls et al 1995, eq 45---------------------------
!
!input:  v_r in units of v_inf
!output: departure coefficients for h-alpha: b2, b3
!
use prog_type
!
implicit none
!
! ... arguments
real(dp), intent(in) :: velr
real(dp) :: b2, b3
!
! ... local scalars
real(dp), parameter :: b2_in=1.5d0, b3_in=1.2d0, &
                       b2_min=1.2d0, b3_min=1.1d0, &
                       b2_inf=1.3d0, b3_inf=1.1d0
!
! ... local functions
!
if(velr.lt.0.01d0) then
   b2 = 1.d0 + (b2_in - 1.d0)*velr/0.01d0
   b3 = 0.9d0 + (b3_in - .9d0)*velr/0.01d0
else if (velr.lt.0.1d0) then
   b2 = b2_in + (b2_min-b2_in)*(velr-0.01d0)/0.09d0
   b3 = b3_in + (b3_min-b3_in)*(velr-0.01d0)/0.09d0
else if (velr.le.1.d0) then
   b2 = b2_min + (b2_inf-b2_min)*(velr-0.1d0)/(0.9d0)
   b3 = b3_min + (b3_inf-b3_min)*(velr-0.1d0)/(0.9d0)
else
   stop 'error in depcoeff_petrenz: wrong velocity range'
endif
!
!
end subroutine depcoeff_petrenz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function opalbar_petrenz(sr, yhe, hei, temp, vth_fiducial, xnue0, b2, b3, rho)
!
!-----------------according to petrenz & puls 1995----------------------
!
use prog_type
use fund_const, only: cgs_clight, pi, cgs_me, cgs_mp, cgs_e
!
implicit none
!
! ... arguments
real(dp), intent(in) :: rho, temp, b2, b3, yhe, hei, sr, vth_fiducial, xnue0
real(dp) :: opalbar_petrenz
!
! ... local scalars
real(dp), parameter :: c1=2.07d-16, gf=5.1256d0
real(dp) ::deldop_fiducial, dum1
!
! ... local arrays
!
dum1=(1.d0+yhe*hei)/(1.d0+4.d0*yhe)**2
dum1=dum1*rho*rho*temp**(-1.5d0)
dum1=dum1*c1*gf*pi*cgs_e*cgs_e/cgs_me/cgs_clight/cgs_mp/cgs_mp
dum1=dum1*(b2*exp(3.954d4/temp) - b3*exp(1.753d4/temp))
!
!so far in cgs, now in own units
!
!finally need to divide this by deldop_fiducial
deldop_fiducial=xnue0*vth_fiducial/cgs_clight
opalbar_petrenz=dum1*sr/deldop_fiducial
!
!
end function opalbar_petrenz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function sline_petrenz(temp, b2, b3)
!
!-----------------according to petrenz & puls 1995----------------------
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: b2, b3, temp
real(dp) :: sline_petrenz
!
! ... local scalars
real(dp) :: c1, c2
!
!
!
!2*h*nu^3/c^2 (for h_alpha)
c1=1.404d-3
!h*nu/k (for h_alpha)
c2=2.192d4
!
sline_petrenz=c1 / (b2*exp(c2/temp)/b3 - 1.d0)
!
end function sline_petrenz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function calc_req(r_pole, m_star, v_rot)
!
!this function calculates the equatorial radius of a rotating star from
! which is approximated by an ellipsoid
!
!inputs:
!   r_pole: polar radius in r_sun
!   m_star: stellar mass in m_sun (without eddington factor!!!)
!   v_rot:  rotational velocity at equator in km/s
!
!outputs:
!   r_eq:   equatorial radius in units of r_pole
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: rsu, xmsu, cgs_grav
!
implicit none
!
! ... arguments
real(dp), intent(in) :: r_pole, m_star, v_rot
real(dp) :: calc_req
!
! ... local scalars
real(dp) :: r_pole_cgs, m_star_cgs, v_rot_cgs, r_eq_cgs, w0
!
r_pole_cgs = r_pole*rsu
m_star_cgs = m_star*xmsu
v_rot_cgs = v_rot*1.d5
!
w0 = v_rot_cgs**2*r_pole_cgs/2.d0/cgs_grav/m_star_cgs
!
r_eq_cgs = r_pole_cgs/(1.d0-w0)
calc_req = r_eq_cgs/r_pole_cgs
!
return
!
end function calc_req
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine surface_rot1(theta, r_pole, m_star, v_rot, r_eq, v_crit, omega, r_surface)
!
!       This function calculates the surface radius of a rotating star from
!       a roche model (e.g. Cranmer/Owocki 1995) for a given co-latitude
!
!input:
!       r_pole: polar radius in r_sun
!       m_star: stellar mass in m_sun (without eddington factor!!!)
!       v_rot:  rotational velocity at equator in km/s
!       theta:       co-latitude (for which radius is calculated)
!
!output:
!       r_eq:   equatorial radius in units of r_pole
!       v_crit: critical velocity in km/s
!       omega:  omega = omega_eq / omega_crit   (omega_eq, omega_crit the
!               angular velocities at equator and critical)
!       r_surface:   radius(theta) of stellar surface in units of r_pole
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: rsu, xmsu, cgs_grav, pi
!
implicit none
!
! ... arguments
real(dp), intent(in) :: theta, r_pole, m_star, v_rot
real(dp), intent(out) :: r_eq, v_crit, omega, r_surface
!
! ... local scalars
real(dp) :: r_pole_cgs, m_star_cgs, v_rot_cgs, r_eq_cgs, w0, v_crit_cgs
!
r_pole_cgs = r_pole*rsu
m_star_cgs = m_star*xmsu
v_rot_cgs = v_rot*1.d5
!
w0 = v_rot_cgs**2*r_pole_cgs/2.d0/cgs_grav/m_star_cgs
!
r_eq_cgs = r_pole_cgs/(1.d0-w0)
r_eq = r_eq_cgs/r_pole_cgs
!
v_crit_cgs = sqrt(2.d0*cgs_grav*m_star_cgs/3.d0/r_pole_cgs)
v_crit = v_crit_cgs/1.d5
!
omega = v_rot_cgs*1.5d0*r_pole_cgs/r_eq_cgs/v_crit_cgs
!
!calculate surface (in units of r_pole)
r_surface=3.d0/omega/sin(theta)*cos((pi+acos(omega*sin(theta)))/3.d0)
!
end subroutine surface_rot1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function sline_depcoeff(xnue0, temp, b2, b3)
!
!calculated line source function for given departure coefficients
!of any given transition with frequency xnue0
!
!all input in cgs
!
!-----------------see personal notes on opacities-----------------------
!
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: b2, b3, temp, xnue0
real(dp) :: sline_depcoeff
!
! ... local scalars
real(dp) :: c1, c2
!
!2*h*nu^3/c^2
c1=two*cgs_planck*xnue0**3 / cgs_clight**2
!h*nu/k
c2=cgs_planck*xnue0/cgs_kb
!
sline_depcoeff = c1 / (b2*exp(c2/temp)/b3 - 1.d0)
!
end function sline_depcoeff
