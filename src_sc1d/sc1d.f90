!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
program prog
!
use prog_type
!
implicit none
!
! ... local scalars
integer(i4b), parameter :: nz=32041, niter=10000
integer(i4b) :: i, nconv
real(dp), parameter :: dev_max=1.d-10, mu=1.d0/sqrt(3.d0)
real(dp) :: eps, zmin, zmax, chi_min, chi_max, err, tau_max, tau_min
real(dp) :: s1, s2, s3, s4, s5
integer(i4b) :: opt_alo, opt_ng, opt_method
!
! ... local arrays
real(dp), dimension(nz) :: z, bnue1d, scont1d, &
                           mint1d, scont1d_old, scont1d_diff, err1d, &
                           tau1d, alo_diag, alo_subdiag, alo_supdiag, &
                           solm0, solm1, solm2, solm3
real(dp), dimension(niter) :: scont_maxcorr
!
! ... local characters
!
! ... local logicals
!
! ... local functions
!
!-----------------------------------------------------------------------
!
!read input
call read_input
!
!create z-grid
call gridz
!
!create model
call setup_mod1d
!
call setup_opacities1
!
call check_grid1d
!
!create angular grid
call calcnodes_mu
!
!create frequency grid
call calcnodes_nue
!
!boundary condition from diffusion approximation
call calc_bcondition
!
!perform radiative transfer
call conttrans_sc1d
!
!
call output
!
!
end program prog
!
!-----------------------------------------------------------------------
!-------------------read in all input parameters------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
use prog_type
use options, only: opt_ng_cont, opt_ait_cont, opt_alo_cont, opt_method
use fund_const, only: rsu, yr, xmsu  
use params_input, only: teff, trad, tmin, xmloss, vmin, vmax, beta, &
                        rstar, input_file, eps_cont, kcont, yhe, hei, &
                        rstar_cgs, xmloss_cgs, vmin_cgs, vmax_cgs
use freq, only: xnue0, nnue
use dime1d, only: nz, zmin, zmax
use angles, only: ntheta
!
implicit none
!
! ... local scalars
real(dp) :: mstar
integer(i4b) :: opt_angint_method
!
! ... local characters
!
! ... local functions
real(dp) :: calc_req
!
! ... namelist
namelist / input_model / teff, trad, tmin, xmloss, vmin, vmax, beta, rstar, xnue0, yhe, hei
namelist / input_cont / eps_cont, kcont
namelist / dimensions_freq / nnue
namelist / dimensions_angles / ntheta
namelist / dimensions_3dz / nz, zmin, zmax
namelist / input_options / opt_ng_cont, opt_ait_cont, opt_alo_cont, opt_angint_method
namelist / input_sc1d / opt_method
!
!-----------------------------------------------------------------------
!
write(*,*) '----------------------------read input-----------------------------------------'
write(*,*) 'input file name (*.nml) to define model-atmosphere'
read(*,*) input_file
write(*,*) 'reading input from file: ', trim(input_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(input_file), status='old', form='formatted')
!
   rewind 1
   read(1, nml=input_options)
!   
!read 1d model parameters
   rewind 1
   read(1, nml=input_model)
   rstar_cgs=rstar*rsu
   vmin_cgs=vmin*1.d5
   vmax_cgs=vmax*1.d5
   xmloss_cgs=xmloss*xmsu/yr
   tmin=teff*tmin
!
!--------------------read dimensions for 1d grid------------------------
!
!
   rewind 1
   read(1, nml=dimensions_3dz)
!
!---------read dimensions/grid-parameters for frequency grid------------
!
   rewind 1
   read(1, nml=dimensions_freq)
!
!-----------read dimensions/grid-parameters for angle grid--------------
!
   rewind 1
   read(1, nml=dimensions_angles)
!
!--------read model parameters for modelling continuum and line---------
!
   rewind 1
   read(1, nml=input_cont)
!
!--------read model parameters for method to be used--------------------
!
   rewind 1
   read(1, nml=input_sc1d)
!
close(1)
!
!
if(opt_method.ne.4.and.opt_method.ne.5) opt_method=4
!
write(*,*) 'input opt_method: ', opt_method
write(*,*) '   0 - '
write(*,*) '   1 - '
write(*,*) '   2 - '
write(*,*) '   3 - '
write(*,*) '   4 - 1st order sc'
write(*,*) '   5 - 2nd order sc'
write(*,*)
!
write(*,*) 'input opt_alo: ', opt_alo_cont
write(*,*) '   0 - classical lambda iteration (alo=0)'
write(*,*) '   1 - diagonal alo'
write(*,*) '   2 - tridiagonal alo'
write(*,*)
!
write(*,*) 'input opt_ng: ', opt_ng_cont
write(*,*) '   0 - NG extrapolation off'
write(*,*) '   1 - NG extrapolation on'
write(*,*)
!
!
end subroutine read_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine gridz
!
use prog_type
use dime1d, only: nz, z, zmin, zmax, imask1d
!
implicit none
!
integer(i4b) :: i, err
real(dp) :: zmin_dum, zmax_dum, zshift, del
!
allocate(z(nz), stat=err)
   if(err.ne.0) stop 'error in gridz: allocation'
!
allocate(imask1d(nz), stat=err)
   if(err.ne.0) stop 'error in gridz: allocation'
!
!----------------------equidistant grid---------------------------------
!
if(nz.le.5) stop 'error in gridxz: nz needs to be larger than 5'   
do i=3, nz-2
   z(i) = zmin + (i-3)*(zmax-zmin)/(nz-5)
enddo
!ghost zones
z(2)=2.*z(3)-z(4)
z(1)=2.*z(2)-z(3)
z(nz-1)=2*z(nz-2)-z(nz-3)
z(nz)=2*z(nz-1)-z(nz-2)
!
imask1d=1
imask1d(1)=0
imask1d(2)=0
imask1d(nz-1)=0
imask1d(nz)=0
!
end subroutine gridz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_mod1d
!
!-------calculate model atmosphere from beta velocity law--------------
!
use prog_type
use fund_const, only: cgs_mp, pi, xmsu, sigmae, one, four
use dime1d, only: z, rho1d, velz1d, t1d, nz
use params_input, only: vmin_cgs, vmax_cgs, beta, xmloss_cgs, yhe, hei, beta, rstar_cgs, tmin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: bconst, zshift
real(dp) :: velr, rad
!
! ... local characters
!
! ... local functions
real(dp) :: bvel
!
!
!-----------------------allocate arrays---------------------------------
!
allocate(rho1d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_mod1d'
allocate(velz1d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_mod1d'
allocate(t1d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_mod1d'
!
!------------------calculate/define required constants------------------
!
!b-factor for beta-velocity-law
bconst=1.d0-(vmin_cgs/vmax_cgs)**(one/beta)
!
if(z(3).lt.0.) stop 'radial beta-velocity law not meaniningful for such z-range'
!
zshift=0.d0
if(z(3).lt.1.) then
   write(*,*) 'shifting z-values to get a radial beta-velocity law starting at r=1.'
   zshift=one-z(3)
endif
!
do i=3, nz-2 
   rad = z(i)+zshift
   velr = bvel(rad, vmax_cgs, bconst, beta)
   velz1d(i) = velr
   rho1d(i) = xmloss_cgs/four/pi/(rstar_cgs*rad)**2/velr
   t1d(i) = tmin
enddo
!ghost zones
velz1d(2)=velz1d(3)
rho1d(2)=rho1d(3)
t1d(2)=t1d(3)
velz1d(1)=velz1d(2)
rho1d(1)=rho1d(2)
t1d(1)=t1d(2)

velz1d(nz-1)=velz1d(nz-2)
rho1d(nz-1)=rho1d(nz-2)
t1d(nz-1)=t1d(nz-2)
velz1d(nz)=velz1d(nz-1)
rho1d(nz)=rho1d(nz-1)
t1d(nz)=t1d(nz-1)
!
!
end subroutine setup_mod1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_opacities1
!
use prog_type
use params_input, only: kcont, yhe, hei
use params_input, only: rstar_cgs, tau_max, tau_min
use dime1d, only: nz, z, rho1d, opac1d, zmax, zmin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: rho, dtau, grad, a, b
!
! ... local logicals
!
! ... local functions
real(dp) :: opac_thomson
!
allocate(opac1d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_opacities1'
!
do i=1, nz
   rho=rho1d(i)
!in units 1/rstar
   opac1d(i) = opac_thomson(yhe, hei, rho, kcont)*rstar_cgs
!   write(*,*) opac1d(i)
enddo
!
!
!use constant opacity (k=1 corresponds to tau=1) for test cases
dtau=1.d0
opac1d=kcont*dtau/(zmax-zmin)
!
!
!
!exponentially decreasing opacity
b=log(tau_max/tau_min)/(zmax-zmin)
a=b*tau_max*exp(b*zmin)
opac1d=a*exp(-b*z)

!
!
end subroutine setup_opacities1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_mu
!
use prog_type
use fund_const, only: pi
use angles, only: ntheta, dim_mu, nodes_mu, weight_mu, q_alo
use mod_integ1d, only: precalc_oweight_trapez
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: rho, dtau
real(dp) :: theta_min, theta_max
!
! ... local logicals
!
! ... local functions
real(dp) :: opac_thomson
!
dim_mu=2*ntheta 
allocate(nodes_mu(dim_mu), stat=err)
  if(err.ne.0) stop 'error: allocation in calcnodes_mu'
allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'error: allocation in calcnodes_mu'
!
!
if(dim_mu.eq.2) then
!to be consistent with diffusion equation (two-stream approximation)
   theta_min=acos(1.d0/sqrt(3.d0))
   theta_max=acos(-1.d0/sqrt(3.d0))
else
   theta_min=1.d-6
   theta_max=pi-1.d-6
endif
!
!
!
do i=1, dim_mu
   nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
   nodes_mu(i) = cos(nodes_mu(i))
   if(abs(nodes_mu(i)).lt.1.d-15) nodes_mu(i)=0.d0
enddo
!
call precalc_oweight_trapez(nodes_mu, dim_mu, -1.d0, 1.d0, weight_mu)
!
!normalization
weight_mu=weight_mu/2.d0
!
!
!set indices for nearest neighbour alo-calculations (direction dependent)
allocate(q_alo(dim_mu,3), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_mu'
q_alo=0
do i=1, dim_mu
   if(nodes_mu(i).gt.0.d0) then
      q_alo(i,:) = (/ 1, 2, 3 /)
   elseif(nodes_mu(i).lt.0.d0) then
      q_alo(i,:) = (/ 3, 2, 1 /)
   else
      stop 'error in calcnodes_mu: mu=0 not allowed'
   endif
enddo
!
end subroutine calcnodes_mu
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_bcondition
!
!calculation of boundary condition: diffusion approximation with
!
!   xic1=bnue(xnue, trad)
!   xic2=dB/dtau(xnue, teff)
!
use prog_type
use params_input, only: trad
use freq, only: nnue, xic1_nue, nodes_nue
use mod_math, only: bnue
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
!
! ... for testing
!
! ... local arrays
!
! ... local functions
!
allocate(xic1_nue(nnue), stat=err)
  if(err.ne.0) stop 'error: allocation in calc_bcondition'
!
do i=1, nnue
   xic1_nue(i)=bnue(nodes_nue(i),trad)
enddo
!
end subroutine calc_bcondition
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calcnodes_nue
!
use prog_type
use fund_const, only: pi
use freq, only: nnue, nodes_nue, weight_nue, xnue0
use dime1d, only: t1d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: rho, dtau
real(dp) :: theta_min, theta_max
!
! ... local logicals
!
! ... local functions
real(dp) :: opac_thomson
!
allocate(nodes_nue(nnue), stat=err)
  if(err.ne.0) stop 'error: allocation in calcnodes_nue'
allocate(weight_nue(nnue), stat=err)
   if(err.ne.0) stop 'error: allocation in calcnodes_nue'
!
if(nnue.eq.1) then
   nodes_nue=xnue0
   weight_nue=1.d0
else
   stop 'todo: calcnodes_nue'
endif
!
end subroutine calcnodes_nue
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine conttrans_sc1d
!
!-------------------------continuum transport---------------------------
!
!---------calculating intensities for different angles mu---------------
!-------------------in plane parallel symmetry-------------------------------
!--------------performing angle integration on z-axis-------------------
!
use prog_type
use iter, only: itmaxc, devmaxc, epsmaxc_arr, nconv
use dime1d, only: nz, z, scont1d, t1d, mint1d, normalization1d, alocont_nn1d, imask1d, bnue1d
use params_input, only: eps_cont
use ng_extra, only: ng_const
use freq, only: nodes_nue
use options, only: opt_ng_cont, opt_ait_cont
use mod_interp2d, only: wp_integ1d
use mod_math, only: bnue, calc_dev1d
use mod_ng_extrapol, only: ng_expol1d, ait_expol1d
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k, l, err, iz
integer(i4b) :: s1, s2, s3, s4
integer(i4b) :: iim2, iim1, ii, iip1
real(dp) :: dummy1, dummy2, rad
real(dp) :: eps_max
!
! ... local arrays
real(dp), dimension(:), allocatable :: eps1d
real(dp), dimension(:,:), allocatable :: scont1d_ng
!
!... local characters
!
! ... local functions
!
! ... local logicals
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------continuum transport (1d)------------------------------'
write(*,*)
!
allocate(bnue1d(nz), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc1d'
allocate(eps1d(nz), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc1d'
allocate(scont1d_ng(4,nz), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc1d'
allocate(mint1d(nz), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc1d'
allocate(scont1d(nz), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc1d'
allocate(normalization1d(nz), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc1d'
allocate(alocont_nn1d(nz,3), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc1d'
!
!-----------------------------------------------------------------------
!
!initialisation of iteration step at which the old solution has to be 
!stored (for ng-extrapolation/aitken extrapolation)
!
s1=1
s2=2
s3=3
s4=4
!
!-----------------------------------------------------------------------
!
!calculating start values of continuum source function
call calc_startval_cont
!
!calculating log(r3d)
!do i=1, ndxmax
!   do j=1, ndymax
!      do k=1, ndzmax
!         log3d(i,j,k)=log10(sqrt(x(i)**2 + y(j)**2 + z(k)**2))
!      enddo
!   enddo
!enddo
!
!calculating planck function along z direction
do i=1, nz
   bnue1d(i) = bnue(nodes_nue(1), t1d(i))
enddo
!
normalization1d=0.d0
alocont_nn1d=0.d0
!
eps1d=0.d0
!
!************************start iteration scheme*************************
!
do i=1, itmaxc
!
!****************************output*************************************
!
   write(*,fmt='(a5, 6(a20))') '#', 'z', 'j', 's_c', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=1, nz
      write(*,fmt='(i5, 10(e20.8))') j, z(j), mint1d(j), scont1d(j), eps1d(j), normalization1d(j), alocont_nn1d(j,1), alocont_nn1d(j,2), alocont_nn1d(j,3)
   end do
   write(*,*)
!
!*************************output end************************************
!
   eps1d=mint1d
!
   call output_itstep_cont(i)
!
!----------------calculating mean intensities on central ray------------
!-------------------------------step i----------------------------------
!
   write(*,*) '--------------calculating angle integrated mean intensities in 1d--------------'
   write(*,*) 'step', i
   write(*,*)
!
!   wp_integ1d=wp_integ1d*0.
   call mint_sc1d
!   wp_integ1d=wp_integ1d*0.95
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)
   call calc_dev1d(eps1d, mint1d, imask1d, nz, eps_max, iz)
   epsmaxc_arr(i)=eps_max
   write(*,'(a30, i4, f8.4, es18.8)') 'max (dev) at grid-point:', iz, z(iz) , eps_max
   write(*,*) 'wp_integ1d', wp_integ1d
   write(*,*)
!
   nconv=i
   if(abs(eps_max).lt.devmaxc) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
!      call cpu_time(te_it)
!      ttot_it = ttot_it + te_it - ts_it
      exit
   else if(i.eq.itmaxc) then
      write(*,*) "no convergence after iteration no. ", i
   end if
!
!-------calculating alo-corrected source-functions on central ray-------
!
   call scont_new1d   
!
!------------extrapolation of old subsequent source functions-----------
!
   if(opt_ng_cont.or.opt_ait_cont) then
!
!----------storing old source-functions for ng-extrapolation------------
!
      if(i.eq.s1) then
         write(*,*) '----------------------storing source fct. at step n-3--------------------------'
         write(*,*)
         scont1d_ng(1,:)=scont1d
         s1=s1+ng_const
      elseif(i.eq.s2) then
         write(*,*) '----------------------storing source fct. at step n-2--------------------------'
         write(*,*)
         scont1d_ng(2,:)=scont1d
         s2=s2+ng_const
      elseif(i.eq.s3) then
         write(*,*) '----------------------storing source fct. at step n-1--------------------------'
         write(*,*)
         scont1d_ng(3,:)=scont1d
         s3=s3+ng_const
      elseif(i.eq.s4) then
         write(*,*) '----------------------storing source fct. at step n----------------------------'
         write(*,*)
         scont1d_ng(4,:)=scont1d
         s4=s4+ng_const
         if(opt_ng_cont) call ng_expol1d(scont1d_ng,nz, verbose=.true.)
         if(opt_ait_cont) call ait_expol1d(scont1d_ng,nz, verbose=.true.)
         scont1d=scont1d_ng(1,:)
      endif

   endif
   !
!  if(i.gt.100) stop 'go on in conttrans'
enddo
!
!
!
end subroutine conttrans_sc1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_startval_cont
!
!-----------------------------------------------------------------------
!---------calculates start value of continuum source function-----------
!-----------------------------------------------------------------------
!
use prog_type
use dime1d, only: nz, z, scont1d, t1d, mint1d, opac1d
use params_input, only: eps_cont
use freq, only: nodes_nue
use mod_math, only: bnue
!
implicit none
!
! ... local scalars
integer(i4b) :: i,j,k
integer(i4b) :: indx_therm3d, indx_therm1d
real(dp) :: opac1, opac2, r1, r2, dtau, tau, thdepth, rad
real(dp) :: r_therm, t_therm, i_core, dum_bnue, dilfac, zshift
!
! ... local arrays
!
! ... local functions
!
!-----------------------calculate thermalization depth------------------
!
if(eps_cont.eq.0.d0) then
   !set thermalization depth to arbitrary value
   thdepth = 10.d0
else
   thdepth = 1.d0/sqrt(eps_cont)
endif
!
!---------------------calculate spherical symmetric tau-----------------
!
!start at outer boundary with tau=0
tau=0.d0
indx_therm1d=1
!
do i=nz-1, 2, -1
   opac1=opac1d(i+1)
   opac2=opac1d(i)
   r1=z(i+1)
   r2=z(i)
   dtau=(opac2+opac1)*(r1-r2)/2.d0
   tau=tau+dtau
   indx_therm1d=indx_therm1d+1
!
   write(*,'(8es20.8,3l5)') z(i), tau, opac1, opac2, r1, r2, dtau, thdepth
   if(tau.gt.thdepth) exit
enddo
!
zshift=1.d0-z(2)
!
r_therm=z(nz+1-indx_therm1d)+zshift
t_therm=t1d(nz+1-indx_therm1d)
i_core=bnue(nodes_nue(1), t_therm)
!
!-----------------------------------------------------------------------
!
scont1d=0.d0
do i=2, nz
   if(z(i).le.r_therm) then
!when thermalized, source function = planck function
      scont1d(i)=bnue(nodes_nue(1), t1d(i))
   else
!when not thermalized, calculate source function from
!mean intensity (i_core * dilution factor)
      rad=z(i)+zshift
      dilfac = r_therm**2 / rad**2
      dilfac = 0.5d0*(1.d0 - sqrt(1.d0 - dilfac))
      dum_bnue=bnue(nodes_nue(1), t1d(i))
      scont1d(i)=(1.d0-eps_cont)*dilfac*i_core + eps_cont*dum_bnue
   endif

enddo
scont1d(2)=bnue(nodes_nue(1),t1d(2))
scont1d(1)=bnue(nodes_nue(1),t1d(1))
!
mint1d=0.d0
!
scont1d=0.d0
scont1d(2)=bnue(nodes_nue(1),t1d(2))
scont1d(1)=bnue(nodes_nue(1),t1d(1))
!
end subroutine calc_startval_cont
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_itstep_cont(itnr)
!
!------------------output after each iteration step---------------------
!input: itnr: current iteration step
!
use prog_type
use fund_const
use dime1d, only: scont1d, mint1d
use iter, only: epsmaxc_arr, itmaxc, devmaxc
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: itnr
character(len=21), parameter :: output_dir_temp='outputFILES_TEMP/sc1d'
!
! ... local scalars
!
! ... local arrays
!
!-----------------------------------------------------------------------
!
write(*,*) '-----------------------storing current iterate in------------------------------'
write(*,*) output_dir_temp//'/iterparam_cont.dat'
write(*,*) output_dir_temp//'/iterepsmax_cont.dat'
write(*,*) output_dir_temp//'/solution_scont1d.dat'
write(*,*) output_dir_temp//'/solution_mint1d.dat'
write(*,*)
!
!----------------------output iteration parameter-----------------------
!
open(1, file=trim(output_dir_temp)//'/iterparam_cont.dat', form='formatted')
   write(1,'(3(a20))') 'itmaxc', 'current it', 'devmax'
   write(1,'(2i20, e20.8)') itmaxc, itnr, devmaxc
close(1)

open(1, file=trim(output_dir_temp)//'/iterepsmax_cont.dat', form='unformatted')
   write(1) epsmaxc_arr
close(1)
!
!-----------------output 3-d solution grids-----------------------------
!
open(1, file=trim(output_dir_temp)//'/solution_scont1d.dat', form='unformatted')
   write(1) scont1d
close(1)
!
open(1, file=trim(output_dir_temp)//'/solution_mint1d.dat', form='unformatted')
   write(1) mint1d
close(1)
!
!
!
end subroutine output_itstep_cont
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output

use prog_type
use dime1d, only: nz, z, scont1d, tau1d, opac1d, bnue1d
use iter, only: nconv, epsmaxc_arr
use params_input, only: tau_min
!
implicit none
!
integer(i4b) :: i, err
real(dp) :: dtau
!
allocate(tau1d(nz), stat=err)
tau1d(nz-2)=tau_min
i=nz-1
dtau=0.5d0*(opac1d(i)+opac1d(i-1))*(z(i)-z(i-1))
tau1d(i)=tau1d(i-1)-dtau
i=nz
dtau=0.5d0*(opac1d(i)+opac1d(i-1))*(z(i)-z(i-1))
tau1d(i)=tau1d(i-1)-dtau
do i=nz-3, 1, -1
   dtau=0.5d0*(opac1d(i)+opac1d(i+1))*(z(i+1)-z(i))
   tau1d(i)=tau1d(i+1)+dtau
enddo
!
open(1, file='outputFILES/sc1d/output.dat', form='formatted')
   do i=1, nz
      write(1,'(i5, 5es20.8)') i, z(i), opac1d(i), tau1d(i), scont1d(i), bnue1d(i)
   enddo
close(1)
!
open(1, file='outputFILES/sc1d/convergence.dat', form='formatted')
   do i=1, nconv
      write(1,'(i5, es20.8)') i, epsmaxc_arr(i)
   enddo
close(1)


end subroutine output
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_grid1d
!
!test the upwind and downwind delta-tau steps
!such that alo-coefficients are not 'extreme'  
!
use prog_type
use dime1d, only: nz, z, opac1d
use options, only: opt_method
use mod_integ1d, only: integ1d_tau_ud
use mod_math, only: calc_dtau_crit
!
implicit none
!
! ... local scalars
integer(i4b) :: k
real(dp) :: delt_u, delt_d, dels_u, dels_d
real(dp) :: nn_z
real(dp) :: opac_u, opac_p, opac_d
real(dp) :: delt_crita, delt_critb, delt_critc
!
!linear interpolations give always stable alo-coefficients
if(opt_method.eq.4) return
!
do k=2, nz-1
   opac_p=opac1d(k)
   opac_u=opac1d(k-1)
   opac_d=opac1d(k+1)

!check n_z=1 direction
   nn_z=1.d0
   dels_u=(z(k)-z(k-1))/nn_z
   dels_d=(z(k+1)-z(k))/nn_z
   call integ1d_tau_ud(dels_u, dels_d, opac_u, opac_p, opac_d, delt_u, delt_d)
   call calc_dtau_crit(delt_u, delt_crita, delt_critb, delt_critc)
 !  write(*,*) k, delt_u, delt_d, delt_crita, delt_critb, delt_critc
   if(delt_d.gt.delt_crita) stop 'error in check_grid1d: delt_d > delt_crita => increase resolution'
   if(delt_d.lt.delt_critb) stop 'error in check_grid1d: delt_d < delt_critb => decrease resolution'
   if(delt_d.lt.delt_critc) stop 'error in check_grid1d: delt_d < delt_critc => decrease resolution'
!
!check n_z=0.707 direction
   nn_z=0.707d0
   dels_u=(z(k)-z(k-1))/nn_z
   dels_d=(z(k+1)-z(k))/nn_z
   call integ1d_tau_ud(dels_u, dels_d, opac_u, opac_p, opac_d, delt_u, delt_d)
   call calc_dtau_crit(delt_u, delt_crita, delt_critb, delt_critc)
!   write(*,*) k, delt_u, delt_d, delt_crita, delt_critb, delt_critc
   if(delt_d.gt.delt_crita) stop 'error in check_grid1d: delt_d > delt_crita => increase resolution'
   if(delt_d.lt.delt_critb) stop 'error in check_grid1d: delt_d < delt_critb => decrease resolution'
   if(delt_d.lt.delt_critc) stop 'error in check_grid1d: delt_d < delt_critc => decrease resolution'
!
!check n_z=1/sqrt(3) direction
   nn_z=1.d0/sqrt(3.d0)
   dels_u=(z(k)-z(k-1))/nn_z
   dels_d=(z(k+1)-z(k))/nn_z
   call integ1d_tau_ud(dels_u, dels_d, opac_u, opac_p, opac_d, delt_u, delt_d)
   call calc_dtau_crit(delt_u, delt_crita, delt_critb, delt_critc)
!   write(*,*) k, delt_u, delt_d, delt_crita, delt_critb, delt_critc
   if(delt_d.gt.delt_crita) stop 'error in check_grid1d: delt_d > delt_crita => increase resolution'
   if(delt_d.lt.delt_critb) stop 'error in check_grid1d: delt_d < delt_critb => decrease resolution'
   if(delt_d.lt.delt_critc) stop 'error in check_grid1d: delt_d < delt_critc => decrease resolution'
   
enddo



end subroutine check_grid1d
