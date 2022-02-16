!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
program prog
!
use prog_type
use options, only: opt_method
use fund_const, only: one
!
implicit none
!
! ... local scalars
!
! ... local arrays
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
!create x,z-grid
call gridxz
!
!create model
call setup_mod2d
!
call setup_opacities1
!
call check_grid2d
!
!create angular grid
call calcnodes_mu
call check_nodes_omega
!
!create frequency grid
call calcnodes_nue
!
!boundary condition from diffusion approximation
call calc_bcondition
!
call allocate_global2d
!
!
call make_benchmark
!
!check if alo coefficients are correct
!call check_alo2d
!call testlambs_direct
!
!perform radiative transfer
call conttrans_sc2d
!
call output
!
!
end program prog
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine allocate_global2d
!
use prog_type
use fund_const
use dime2d, only: nx, nz, bnue2d, mint2d, scont2d, normalization2d, alocont_nn2d, intbound2d
use angles, only: dim_omega
! 
implicit none
integer(i4b) :: err
!
write(*,*) '-------------------------allocating all global 2d arrays-----------------------'
write(*,*)
!
allocate(bnue2d(nx,nz), stat=err)
   if(err.ne.0) stop 'allocation error in allocate_global2d'
allocate(mint2d(nx,nz), stat=err)
   if(err.ne.0) stop 'allocation error in allocate_global2d'
allocate(scont2d(nx,nz), stat=err)
   if(err.ne.0) stop 'allocation error in allocate_global2d'
allocate(normalization2d(nx,nz), stat=err)
   if(err.ne.0) stop 'allocation error in allocate_global2d'
allocate(alocont_nn2d(nx,nz,27), stat=err)
   if(err.ne.0) stop 'allocation error in allocate_global2d'
allocate(intbound2d(2,nz,dim_omega), stat=err)
if(err.ne.0) stop 'allocation error in allocate_global2d'
intbound2d=zero

end subroutine allocate_global2d
!
!-----------------------------------------------------------------------
!-------------------read in all input parameters------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
use prog_type
use options, only: opt_ng_cont, opt_ait_cont, opt_alo_cont, opt_method
use fund_const, only: rsu, yr, xmsu, pi
use params_input, only: teff, trad, tmin, xmloss, vmin, vmax, beta, &
                        rstar, input_file, eps_cont, kcont, yhe, hei, &
                        rstar_cgs, xmloss_cgs, vmin_cgs, vmax_cgs
use freq, only: xnue0, nnue
use dime2d, only: nx, nz, xmin, xmax, zmin, zmax
use angles, only: ntheta
use mod_benchmark, only: n_z, nn_z, benchmark_mod
!
implicit none
!
! ... local scalars
real(dp) :: mstar, theta, phi
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
namelist / dimensions_3dx / nx, xmin, xmax
namelist / dimensions_3dz / nz, zmin, zmax
namelist / input_options / opt_ng_cont, opt_ait_cont, opt_alo_cont, opt_angint_method
namelist / input_sc2d / opt_method
namelist / benchmark / benchmark_mod, theta, phi
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
!read 2d model parameters
   rewind 1
   read(1, nml=input_model)
   rstar_cgs=rstar*rsu
   vmin_cgs=vmin*1.d5
   vmax_cgs=vmax*1.d5
   xmloss_cgs=xmloss*xmsu/yr
   tmin=teff*tmin
!
!--------------------read dimensions for 2d grid------------------------
!
   rewind 1
   read(1, nml=dimensions_3dx)
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
   read(1, nml=input_sc2d)
!
!--------read benchmark models------------------------------------------
!
   rewind 1
   read(1, nml=benchmark)
   n_z=cos(theta*pi/180.d0)
   nn_z=cos(theta*pi/180.d0)
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
subroutine gridxz
!
use prog_type
use fund_const  
use dime2d, only: nx, nz, x, z, xmin, xmax, zmin, zmax, imask2d, imaskb2d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, indx_1d, indx_x, indx_z, nxp, nxm
real(dp) :: zmin_dum, zmax_dum, zshift, del, dxp, dxm, xmin_dum, xmax_dum
!
! ... local arrays
real(dp), dimension(:), allocatable :: x_dump, x_dumm
!
write(*,*) '-------------------------creating x and z grid---------------------------------'
write(*,*)

allocate(x(nx), stat=err)
   if(err.ne.0) stop 'error in gridxz: allocation'
allocate(z(nz), stat=err)
   if(err.ne.0) stop 'error in gridxz: allocation'
!
allocate(imask2d(nx,nz), stat=err)
   if(err.ne.0) stop 'error in gridz: allocation'
allocate(imaskb2d(nx,nz), stat=err)
   if(err.ne.0) stop 'error in gridz: allocation'
!
!----------------------equidistant grid in x----------------------------
!
if(nx.le.5) stop 'error in gridxz: nx needs to be larger than 5'   
do i=3, nx-2
   x(i) = xmin + (i-3)*(xmax-xmin)/(nx-5)
enddo
!ghost zones
x(2)=two*x(3)-x(4)
x(1)=two*x(2)-x(3)
x(nx-1)=two*x(nx-2)-x(nx-3)
x(nx)=two*x(nx-1)-x(nx-2)
!
!
!----------------------logarithmic grid in x----------------------------
!
!if(nx.le.5) stop 'error in gridxz: nx needs to be larger than 5'   
!if(xmin.ge.zero) stop 'error in gridxz: xmin needs to be less than 0 here'
!!
!x=zero
!!
!dxp=xmax-zero
!dxm=zero-xmin
!
!nxp=(nx-5)*dxp/dxm/(one+dxp/dxm)
!nxm=nx-5-nxp
!
!xmin_dum=zero+one
!xmax_dum=xmax+five+one
!allocate(x_dump(nxp))
!call grid_log(xmin_dum, xmax_dum, nxp, x_dump)
!x_dump = (x_dump-one)/(xmax_dum-one) * xmax
!x_dump(1) = x_dump(2)/two
!!
!!
!!
!xmax_dum=abs(xmin-five-one)
!xmin_dum=one
!allocate(x_dumm(nxm))
!call grid_log(xmin_dum,xmax_dum,nxm, x_dumm)
!x_dumm = (x_dumm-one)/(xmax_dum-one) * xmin
!x_dumm(1) = x_dumm(2)/two
!!
!do i=3, nxm+2
!   j=nxm+3-i
!   x(i) = x_dumm(j)
!enddo
!x(nxm+3)=zero
!do i=nxm+4, nx-2
!   j=i-nxm-3
!   x(i) = x_dump(j)
!enddo
!
!ghost zones
!x(2)=two*x(3)-x(4)
!x(1)=two*x(2)-x(3)
!x(nx-1)=two*x(nx-2)-x(nx-3)
!x(nx)=two*x(nx-1)-x(nx-2)
!
!
!----------------------equidistant grid in z----------------------------
!
if(nz.le.5) stop 'error in gridxz: nz needs to be larger than 5'   
do i=3, nz-2
   z(i) = zmin + (i-3)*(zmax-zmin)/(nz-5)
enddo
!ghost zones
z(2)=two*z(3)-z(4)
z(1)=two*z(2)-z(3)
z(nz-1)=two*z(nz-2)-z(nz-3)
z(nz)=two*z(nz-1)-z(nz-2)
!
!
imask2d=1
do i=1, nx
   imask2d(i,1)=0
   imask2d(i,2)=0
   imask2d(i,nz-1)=0
   imask2d(i,nz)=0
enddo
!do i=1, nz
!   imask2d(1,i)=0
!   imask2d(2,i)=0
!   imask2d(nx-1,i)=0
!   imask2d(nx,i)=0
!enddo
!
!standard points within computational domain are labelled 9
imaskb2d=9
do i=1, nx
   imaskb2d(i,1)=0     !lower most boundary is labelled 0
   imaskb2d(i,nz)=0    !upper most boundary is labelled 0
   imaskb2d(i,2)=1     !next lower most boundary is labelled 1
   imaskb2d(i,nz-1)=2  !next upper most boundary is labelled 2
enddo
do i=3, nz-2
   imaskb2d(1,i)=3    !left left boundary is labelled 3
   imaskb2d(2,i)=4    !left boundary is labelled 4
   imaskb2d(3,i)=5    !left main point is labelled 5
   imaskb2d(nx-2,i)=6 !right main point is labelled 6
   imaskb2d(nx-1,i)=7 !right boundary is labelled 7
   imaskb2d(nx,i)=8   !right right boundary is labelled 8
enddo
!
end subroutine gridxz
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_mod2d
!
!-------calculate model atmosphere from beta velocity law--------------
!
use prog_type
use fund_const, only: cgs_mp, pi, xmsu, sigmae, one, four, zero
use dime2d, only: x, z, rho2d, velz2d, t2d, nx, nz
use params_input, only: vmin_cgs, vmax_cgs, beta, xmloss_cgs, yhe, hei, beta, rstar_cgs, tmin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, k, err
real(dp) :: bconst, zshift
real(dp) :: velr, rad
!
! ... local characters
!
! ... local functions
real(dp) :: bvel
!
write(*,*) '-----------------------creating 2d model atmosphere----------------------------'
write(*,*)
!
!-----------------------allocate arrays---------------------------------
!
allocate(rho2d(nx,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_mod2d'
allocate(velz2d(nx,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_mod2d'
allocate(t2d(nx,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_mod2d'
!
!------------------calculate/define required constants------------------
!
!b-factor for beta-velocity-law
bconst=one-(vmin_cgs/vmax_cgs)**(one/beta)
!
if(z(3).lt.0.) stop 'radial beta-velocity law not meaniningful for such z-range'
!
zshift=zero
if(z(3).lt.1.) then
   write(*,*) 'shifting z-values to get a radial beta-velocity law starting at r=1.'
   write(*,*)
   zshift=one-z(3)
endif
!
do i=1, nx
!
   do k=3, nz-2 
      rad = z(k)+zshift
      velr = bvel(rad, vmax_cgs, bconst, beta)
      velz2d(i,k) = velr
      rho2d(i,k) = xmloss_cgs/four/pi/(rstar_cgs*rad)**2/velr
      t2d(i,k) = tmin
   enddo
!
!ghost zones
   velz2d(i,2)=velz2d(i,3)
   rho2d(i,2)=rho2d(i,3)
   t2d(i,2)=t2d(i,3)
   velz2d(i,1)=velz2d(i,2)
   rho2d(i,1)=rho2d(i,2)
   t2d(i,1)=t2d(i,2)

   velz2d(i,nz-1)=velz2d(i,nz-2)
   rho2d(i,nz-1)=rho2d(i,nz-2)
   t2d(i,nz-1)=t2d(i,nz-2)
   velz2d(i,nz)=velz2d(i,nz-1)
   rho2d(i,nz)=rho2d(i,nz-1)
   t2d(i,nz)=t2d(i,nz-1)

enddo
!
!
end subroutine setup_mod2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine setup_opacities1
!
use prog_type
use fund_const  
use params_input, only: kcont, yhe, hei
use params_input, only: rstar_cgs
use dime2d, only: nx, nz, x, z, rho2d, opac2d, xmin, xmax, zmax, zmin
use mod_benchmark, only: tau_min, tau_max
!
implicit none
!
! ... local scalars
integer(i4b) :: i, k, err
real(dp) :: rho, dtau, grad, a, b
!
! ... local logicals
!
! ... local functions
real(dp) :: opac_thomson
!
write(*,*) '-------------------------setting up opacities----------------------------------'
write(*,*)
!
allocate(opac2d(nx, nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_opacities1'
!
do i=1, nx
   do k=1, nz
      rho=rho2d(i,k)
!in units 1/rstar
      opac2d(i,k) = opac_thomson(yhe, hei, rho, kcont)*rstar_cgs
!   write(*,*) opac1d(i)
   enddo
enddo
!
!
!use constant opacity (k=1 corresponds to tau=1) for test cases
dtau=one
opac2d=kcont*dtau/(zmax-zmin)
!
!exponentially decreasing opacity
b=log(tau_max/tau_min)/(zmax-zmin)
a=b*tau_max*exp(b*zmin)
do i=1, nx
   opac2d(i,:)=a*exp(-b*z)
enddo
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
use fund_const
use angles, only: ntheta, dim_mu, dim_omega, nodes_mu, weight_mu, q_alo, n_x, n_z, weight_omega
use mod_integ1d, only: precalc_oweight_trapez
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err, indx
real(dp) :: rho, dtau
real(dp) :: theta_min, theta_max
!
! ... local logicals
!
! ... local functions
real(dp) :: opac_thomson

write(*,*) '---------------------------creating angular grid-------------------------------'
write(*,*)
!
dim_mu=2*ntheta
dim_omega=2*dim_mu
allocate(nodes_mu(dim_mu), stat=err)
  if(err.ne.0) stop 'error: allocation in calcnodes_mu'
allocate(weight_mu(dim_mu), stat=err)
   if(err.ne.0) stop 'error: allocation in calcnodes_mu'
allocate(n_x(dim_omega), stat=err)
   if(err.ne.0) stop 'error: allocation in calcnodes_mu'
allocate(n_z(dim_omega), stat=err)
   if(err.ne.0) stop 'error: allocation in calcnodes_mu'
allocate(weight_omega(dim_omega), stat=err)
   if(err.ne.0) stop 'error: allocation in calcnodes_mu'
!
!
if(dim_mu.eq.2) then
!to be consistent with diffusion equation (two-stream approximation)
   theta_min=acos(one/sqrt(three))
   theta_max=acos(-one/sqrt(three))
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
   if(abs(nodes_mu(i)).lt.small_number) nodes_mu(i)=zero
enddo
call precalc_oweight_trapez(nodes_mu, dim_mu, -one, one, weight_mu)
weight_mu=weight_mu/two
!
!write(*,*) weight_mu
!write(*,*) nodes_mu
!stop

indx=1
do i=1, dim_mu
   n_z(indx) = nodes_mu(i)
   n_z(indx+1) = nodes_mu(i)
   n_x(indx) = sqrt(one-nodes_mu(i)**2)
   n_x(indx+1) = -sqrt(one-nodes_mu(i)**2)
   weight_omega(indx)=weight_mu(i)/two
   weight_omega(indx+1)=weight_mu(i)/two
   indx=indx+2
enddo
!
!
!normalization
!
!
!set indices for nearest neighbour alo-calculations (direction dependent)
allocate(q_alo(dim_omega,27), stat=err)
   if(err.ne.0) stop 'allocation error in calcnodes_mu'
q_alo=0
do i=1, dim_omega
   if(n_x(i).gt.0..and.n_z(i).gt.0.) then
      q_alo(i,:) = (/ 1, 1, 1, 4, 5, 6, 1, 1, 1, &
                      1, 1, 1, 13, 14, 15, 1, 1, 1, &
                      1, 1, 1, 22, 23, 24, 1, 1, 1 /)
   elseif(n_x(i).lt.0..and.n_z(i).gt.0.) then
      q_alo(i,:) = (/ 1, 1, 1, 6, 5, 4, 1, 1, 1, &
                      1, 1, 1, 15, 14, 13, 1, 1, 1, &
                      1, 1, 1, 24, 23, 22, 1, 1, 1 /)
   elseif(n_x(i).gt.0..and.n_z(i).lt.0.) then
      q_alo(i,:) = (/ 1, 1, 1, 22, 23, 24, 1, 1, 1, &
                      1, 1, 1, 13, 14, 15, 1, 1, 1, &
                      1, 1, 1, 4, 5, 6, 1, 1, 1 /)
   elseif(n_x(i).lt.0..and.n_z(i).lt.0.) then
      q_alo(i,:) = (/ 1, 1, 1, 24, 23, 22, 1, 1, 1, &
                      1, 1, 1, 15, 14, 13, 1, 1, 1, &
                      1, 1, 1, 6, 5, 4, 1, 1, 1 /)
   else
      stop 'error in calcnodes_mu: mu=0 not allowed'
   endif
enddo
!
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
write(*,*) '----------------------setting the boundary condition---------------------------'
write(*,*)
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
use fund_const
use freq, only: nnue, nodes_nue, weight_nue, xnue0
use dime2d, only: t2d
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
write(*,*) '---------------------------creating frequency grid-----------------------------'
write(*,*)
!
allocate(nodes_nue(nnue), stat=err)
  if(err.ne.0) stop 'error: allocation in calcnodes_nue'
allocate(weight_nue(nnue), stat=err)
   if(err.ne.0) stop 'error: allocation in calcnodes_nue'
!
if(nnue.eq.1) then
   nodes_nue=xnue0
   weight_nue=one
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
subroutine conttrans_sc2d
!
!-------------------------continuum transport---------------------------
!
!---------calculating intensities for different angles mu---------------
!-------------------in plane parallel symmetry-------------------------------
!--------------performing angle integration on z-axis-------------------
!
use prog_type
use fund_const  
use iter, only: itmaxc, devmaxc, epsmaxc_arr, nconvc
use dime2d, only: nx, nz, x, z, scont2d, t2d, mint2d, normalization2d, alocont_nn2d, imask2d, imaskb2d, bnue2d
use params_input, only: eps_cont
use ng_extra, only: ng_const
use freq, only: nodes_nue
use options, only: opt_ng_cont, opt_ait_cont, opt_method
use mod_interp2d, only: wp_integ1d, wpa_integ1d, wpb_integ1d, wp_interp1d, wpa_interp1d, wpb_interp1d, wp_interp2d, &
     wpa_interp2d, wpb_interp2d, interpolation_threshold
use mod_math, only: calc_dev2d, bnue
use mod_ng_extrapol, only: store_ng2d, ng_expol2d, ait_expol2d
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k, l, err, ix, iz, indx_threshold
integer(i4b) :: s1, s2, s3, s4, s4b
integer(i4b) :: iim2, iim1, ii, iip1
real(dp) :: dummy1, dummy2, rad
real(dp) :: eps_max
!
! ... local arrays
real(dp), dimension(:,:), allocatable :: eps2d
real(dp), dimension(:,:), allocatable :: scont2d_ng
!
!... local characters
!
! ... local functions
!
! ... local logicals
!
!-----------------------------------------------------------------------
!
write(*,*) '-------------------------continuum transport (2d sc)---------------------------'
write(*,*)
!
allocate(eps2d(nx,nz), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc2d'
allocate(scont2d_ng(4,nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error in conttrans_sc2d'
!-----------------------------------------------------------------------
!
!initialisation of iteration step at which the old solution has to be 
!stored (for ng-extrapolation/aitken extrapolation)
!
s1=1
s2=2
s3=3
s4=4
s4b=4
!
!setting index for threshold
indx_threshold=5
!
!-----------------------------------------------------------------------
!
!calculating start values of continuum source function
call calc_startval_cont
!
!calculating planck function along z direction
do i=1, nx
   do j=1, nz
      bnue2d(i,j) = bnue(nodes_nue(1), t2d(i,j))
   enddo
enddo
!
normalization2d=zero
alocont_nn2d=zero
!
eps2d=zero
scont2d_ng=zero
!
!************************start iteration scheme*************************
!
do i=1, itmaxc
!
!****************************output*************************************
!
   write(*,fmt='(a5, 6(a20))') '#', 'z', 'j', 's_c', 'deviation(old-new)', 'normalization', 'alo_diag'
   do j=1, nz
      write(*,fmt='(i5, 10(e20.8))') j, z(j), mint2d(nx/2+1,j), scont2d(nx/2+1,j), eps2d(nx/2+1,j), normalization2d(nx/2+1,j), alocont_nn2d(nx/2+1,j,14)
   end do
   write(*,*)
!   do j=1, nz
!      write(*,fmt='(i5, 10(e20.8))') j, alocont_nn2d(nx/2+1-6,j,5), alocont_nn2d(nx/2+1,j,14), alocont_nn2d(nx/2+1,j,23)
!   end do
!   write(*,*)
!
!*************************output end************************************
!
   eps2d=mint2d
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
!   if(i.le.s4b) then
!use linear interpolations for the first s4b iteration steps
!in order to get a first guess of a 'smooth' solution.
!otherwise: convergence behaviour might be oscillatory
!      call mint_sc2d(4)
!   else
      call mint_sc2d(opt_method)
!   endif
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)
   call calc_dev2d(eps2d, mint2d, imask2d, nx, nz, eps_max, ix, iz)
   epsmaxc_arr(i)=eps_max
   write(*,'(a30, 2i4, 2f14.4, es18.8)') 'max (dev) at grid-point:', ix, iz, x(ix), z(iz) , eps_max
   write(*,*)
!
   nconvc=i
!
   if(abs(eps_max).lt.devmaxc) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
      write(*,*)
!for thin continua, ensure that higher order interpolation scheme
!   has been used, and not only linear approach
      if(i.gt.s4b) exit
   else if(i.eq.itmaxc) then
      write(*,*) "no convergence after iteration no. ", i
      write(*,*)
   end if
   
   if(i.gt.indx_threshold) then
      if(abs(epsmaxc_arr(i)).ge.abs(epsmaxc_arr(i-1)).and. &
         abs(epsmaxc_arr(i-1)).le.abs(epsmaxc_arr(i-2)).and. &
         abs(epsmaxc_arr(i-2)).ge.abs(epsmaxc_arr(i-3)).and. &
         abs(epsmaxc_arr(i-3)).le.abs(epsmaxc_arr(i-4)).and. &
         abs(epsmaxc_arr(i-4)).ge.abs(epsmaxc_arr(i-5))) then
!
         write(*,*) 'error in iteration scheme: oscillations!!!'
         write(*,*) 'max deviation at iteration i-5', epsmaxc_arr(i-5)
         write(*,*) 'max deviation at iteration i-4', epsmaxc_arr(i-4)
         write(*,*) 'max deviation at iteration i-3', epsmaxc_arr(i-3)
         write(*,*) 'max deviation at iteration i-2', epsmaxc_arr(i-2)
         write(*,*) 'max deviation at iteration i-1', epsmaxc_arr(i-1)
         write(*,*) 'max deviation at iteration i  ', epsmaxc_arr(i)
         write(*,*) 'possible solutions: '
         write(*,*) '   1. use linear interpolations for upwind/downwind source function'
         write(*,*) '      and for upwind intensities to have same solution procedure in'
         write(*,*) '      each iteration step (independent of the source function and intensities'
         write(*,*) '      themselves.'
         write(*,*) '   2. avoid monotonicity constraint in integration of source contribution'
         write(*,*) '      (try linear approximation of source function along ray)'
         write(*,*) '   3. increase the spatial grid resolution, in order to avoid'
         write(*,*) '      monotonicity constraints in quadratic interpolation procedures'
         write(*,*) '   4. increasing interpolation threshold in quadratic interpolation procedure'
         interpolation_threshold=min(interpolation_threshold+0.1d0,one)
         write(*,*) 'setting interpolation threshold to', interpolation_threshold
         wpa_interp2d=wpa_interp2d+one
         wpb_interp2d=wpb_interp2d+one
         wp_interp2d=wpa_interp2d/wpb_interp2d
         wpa_interp1d=wpa_interp1d+one
         wpb_interp1d=wpb_interp1d+one
         wp_interp1d=wpa_interp1d/wpb_interp1d
         wpa_integ1d=wpa_integ1d+one
         wpb_integ1d=wpb_integ1d+one
         wp_integ1d=wpa_integ1d/wpb_integ1d
         write(*,*) 'setting derivative weights for 2d bezier interpolation from, to', (wpa_interp2d-one)/(wpb_interp2d-one), wp_interp2d
         write(*,*) 'setting derivative weights for 1d bezier interpolation from, to', (wpa_interp1d-one)/(wpb_interp1d-one), wp_interp1d
         write(*,*) 'setting derivative weights for 1d bezier integration   from, to', (wpa_integ1d-one)/(wpb_integ1d-one), wp_integ1d

!start ng-extrapolation from beginning
         s1=i+1
         s2=i+2
         s3=i+3
         s4=i+4
         indx_threshold=i+5
!
      endif
   endif
!
!-------calculating alo-corrected source-functions on central ray-------
!
   call scont_new2d   
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
         call store_ng2d(1,scont2d_ng,nx,nz,scont2d)
         s1=s1+ng_const
      elseif(i.eq.s2) then
         write(*,*) '----------------------storing source fct. at step n-2--------------------------'
         write(*,*)
         call store_ng2d(2,scont2d_ng,nx,nz,scont2d)
         s2=s2+ng_const
      elseif(i.eq.s3) then
         write(*,*) '----------------------storing source fct. at step n-1--------------------------'
         write(*,*)
         call store_ng2d(3,scont2d_ng,nx,nz,scont2d)
         s3=s3+ng_const
      elseif(i.eq.s4) then
         write(*,*) '----------------------storing source fct. at step n----------------------------'
         write(*,*)
         call store_ng2d(4,scont2d_ng,nx,nz,scont2d)
         s4=s4+ng_const

         
         if(opt_ng_cont) call ng_expol2d(scont2d_ng,nx,nz,scont2d)
         if(opt_ait_cont) call ait_expol2d(scont2d_ng,nx,nz,scont2d)
      endif

   endif
!
!   if(i.ge.1) stop 'go on in conttrans'
enddo
!
write(*,*)
write(*,*) 'finally used derivative weights for 2d bezier interpolation', wp_interp2d
write(*,*) 'finally used derivative weights for 1d bezier interpolation', wp_interp1d
write(*,*) 'finally used derivative weights for 1d bezier integration  ', wp_integ1d
write(*,*)
!
!
end subroutine conttrans_sc2d
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
use fund_const  
use dime2d, only: nx, nz, x, z, scont2d, t2d, mint2d, opac2d
use params_input, only: eps_cont
use freq, only: nodes_nue
use mod_math, only: bnue
!
implicit none
!
! ... local scalars
integer(i4b) :: i,k
integer(i4b) :: indx_therm3d, indx_therm1d, indx_x
real(dp) :: opac1, opac2, r1, r2, dtau, tau, thdepth, rad
real(dp) :: r_therm, t_therm, i_core, dum_bnue, dilfac, zshift
!
! ... local arrays
!
! ... local functions
!
!-----------------------calculate thermalization depth------------------
!
if(eps_cont.eq.zero) then
   !set thermalization depth to arbitrary value
   thdepth = 10.d0
else
   thdepth = one/sqrt(eps_cont)
endif
!
!---------------------calculate spherical symmetric tau-----------------
!
!start at outer boundary with tau=0
tau=zero
indx_therm1d=1
indx_x=nx/2+1
!
do i=nz-1, 2, -1
   opac1=opac2d(indx_x,i+1)
   opac2=opac2d(indx_x,i)
   r1=z(i+1)
   r2=z(i)
   dtau=(opac2+opac1)*(r1-r2)/two
   tau=tau+dtau
   indx_therm1d=indx_therm1d+1
!
!   write(*,'(8es20.8,3l5)') z(i), tau, opac1, opac2, r1, r2, dtau, thdepth
   if(tau.gt.thdepth) exit
enddo
!
zshift=one-z(3)
!
r_therm=z(nz+1-indx_therm1d)+zshift
t_therm=t2d(indx_x,nz+1-indx_therm1d)
i_core=bnue(nodes_nue(1), t_therm)
!
!-----------------------------------------------------------------------
!
scont2d=zero
do i=1, nx
   do k=3, nz
      if(z(k).le.r_therm) then
!when thermalized, source function = planck function
         scont2d(i,k)=bnue(nodes_nue(1), t2d(i,k))
      else
!when not thermalized, calculate source function from
!mean intensity (i_core * dilution factor)
         rad=z(k)+zshift
         dilfac = r_therm**2 / rad**2
         dilfac = half*(one - sqrt(one - dilfac))
         dum_bnue=bnue(nodes_nue(1), t2d(i,k))
         scont2d(i,k)=(one-eps_cont)*dilfac*i_core + eps_cont*dum_bnue
      endif
   enddo
!
   scont2d(i,2)=bnue(nodes_nue(1),t2d(i,2))
   scont2d(i,1)=bnue(nodes_nue(1),t2d(i,1))

enddo
!
mint2d=zero
!
!
!or just zero initial condition
scont2d=zero
do i=1, nx
   scont2d(i,2)=bnue(nodes_nue(1),t2d(i,2))
   scont2d(i,1)=bnue(nodes_nue(1),t2d(i,1))
enddo


!scont2d=zero
!do i=1, nx
!   scont2d(i,2)=one
!   scont2d(i,1)=one
!enddo
!!scont2d=scont2d(1,1)
!scont2d=one
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
use dime2d, only: scont2d, mint2d
use iter, only: epsmaxc_arr, itmaxc, devmaxc
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: itnr
character(len=21), parameter :: output_dir_temp='outputFILES_TEMP/sc2d'
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
open(1, file=trim(output_dir_temp)//'/solution_scont2d.dat', form='unformatted')
   write(1) scont2d
close(1)
!
open(1, file=trim(output_dir_temp)//'/solution_mint2d.dat', form='unformatted')
   write(1) mint2d
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
subroutine make_benchmark
!
use prog_type
use mod_benchmark, only: benchmark_mod
use options, only: opt_method
!
implicit none
!
! ... local scalars
integer(i4b) :: i
!
!calculate delx, dely, delz grids if not done yet
!
select case(benchmark_mod)
   case(1)
      write(*,*) '-----------performing benchmark model 1: searchligh beam test 2d---------------'
      write(*,*)
!calculating solution
      call benchmark01_solution
!output to file
      call output_benchmark01

!
   case default
      write(*,*) '----------------------no benchmark is being performed--------------------------'
      write(*,*)
      return
end select
!
!stop program, because opacities and source functions have been overwritten
if(benchmark_mod.ne.0) then
   write(*,*)
   stop '----------------------benchmark and main program done--------------------------'
endif
!
!
end subroutine make_benchmark
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine benchmark01_solution
!
use prog_type
use fund_const
use dime2d, only: z, int2d, scont2d, opac2d, alocont_o_nn2d, alocont_nn2d_tmp, &
                  normalization2d_tmp, mint2d_tmp, nx, nz, t2d, bnue2d, intbound2d, zmin, zmax
use angles, only: dim_omega, nodes_mu, n_z, n_x
use freq, only: nodes_nue
use ng_extra, only: ng_const
use options, only: opt_ng_cont, opt_ait_cont, opt_method
use mod_benchmark, only: int2d_theo, itmaxi, nconvi, devmaxi, epsmaxi_arr, nn_z, nn_x
use params_input, only: kcont
use mod_math, only: bnue, calc_dev2d
use mod_ng_extrapol, only: store_ng2d, ng_expol2d, ait_expol2d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err, ix, iz, oindx
integer(i4b) :: startz, endz, gamma
integer(i4b) :: s1, s2, s3, s4, s4b
real(dp) :: mu, eps_max, dtau, ts, te
!
! ... local arrays
real(dp), dimension(:,:), allocatable :: eps2d
real(dp), dimension(:,:), allocatable :: intbound2d_ng, intbound2d_single
integer(i4b), dimension(:,:), allocatable :: imaskbound2d
!
! ... local functions
!
! ... local characters
!
call cpu_time(ts)
!
!use constant opacity (k=1 corresponds to tau=1) for test cases
dtau=one
opac2d=kcont*dtau/(zmax-zmin)
!
nodes_mu=nn_z
n_z=nn_z
n_x=-sqrt(one-n_z**2)
oindx=1

if(.not.allocated(int2d)) allocate(int2d(nx,nz))
if(.not.allocated(alocont_o_nn2d)) allocate(alocont_o_nn2d(nx,nz,27))
if(.not.allocated(alocont_nn2d_tmp)) allocate(alocont_nn2d_tmp(nx,nz,27))
if(.not.allocated(normalization2d_tmp)) allocate(normalization2d_tmp(nx,nz))
if(.not.allocated(mint2d_tmp)) allocate(mint2d_tmp(nx,nz))
allocate(eps2d(2,nz), stat=err)
allocate(imaskbound2d(2,nz), stat=err)
allocate(intbound2d_single(2,nz), stat=err)
allocate(intbound2d_ng(4,2*nz), stat=err)
eps2d=zero
intbound2d_ng=zero
imaskbound2d=1
!
scont2d=zero
!calculating planck function
do i=1, nx
   do j=1, nz
      bnue2d(i,j) = bnue(nodes_nue(1), t2d(i,j))
   enddo
enddo
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
s4b=4

epsmaxi_arr=zero
!
!-----------------------------------------------------------------------
!
!************************start iteration scheme*************************
!
do i=1, itmaxi
!
   if(.not.allocated(alocont_nn2d_tmp)) allocate(alocont_nn2d_tmp(nx,nz,27))
   if(.not.allocated(normalization2d_tmp)) allocate(normalization2d_tmp(nx,nz))
   if(.not.allocated(mint2d_tmp)) allocate(mint2d_tmp(nx,nz))
!
   eps2d=intbound2d(:,:,oindx)
!
   write(*,*) '--------------calculating intensity for a given ray----------------------------'
   write(*,*) 'step', i
!   write(*,*)
!   write(*,fmt='(a5, 4(a20))') '#', 'z', 'intboundary(1)', 'intboundary(2)', 'int_theo (chi=const)'
!   do j=1, nz
!      write(*,fmt='(i5, 10(e20.8))') j, z(j), intbound2d(1,j,oindx), intbound2d(2,j,oindx), int2d(nx/2+1,2)*exp(-opac2d(nx/2+1,2)*(z(j)-z(2))/n_z(oindx))
!   end do
!   write(*,*)
!
   if(opt_method.eq.4) then
      call fsc_cont2d_lin(oindx,1)
   elseif(opt_method.eq.5) then
      call fsc_cont2d(oindx,1)
   endif
   intbound2d_single=intbound2d(:,:,oindx)
!
!-------------calculating percentage-error of mean intensities----------
!
   write(*,*) '-----------------------------calculating deviation-----------------------------'
   write(*,*)
   call calc_dev2d(eps2d, intbound2d_single, imaskbound2d, 2, nz, eps_max, ix, iz)
   epsmaxi_arr(i)=eps_max
   write(*,'(a30, 2i4, f8.4, es18.8)') 'max (dev) at grid-point:', ix, iz, z(iz) , eps_max
   write(*,*)
!
   nconvi=i
!
   if(abs(eps_max).lt.devmaxi) then
      write(*,*) "convergence after iteration no. ", i
      write(*,*) "max (dev): ", eps_max
      write(*,*)
      exit
   else if(i.eq.itmaxi) then
      write(*,*) "no convergence after iteration no. ", i
      write(*,*)
   end if
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
         call store_ng2d(1,intbound2d_ng,2,nz,intbound2d_single)
         s1=s1+ng_const
      elseif(i.eq.s2) then
         write(*,*) '----------------------storing source fct. at step n-2--------------------------'
         write(*,*)
         call store_ng2d(2,intbound2d_ng,2,nz,intbound2d_single)
         s2=s2+ng_const
      elseif(i.eq.s3) then
         write(*,*) '----------------------storing source fct. at step n-1--------------------------'
         write(*,*)
         call store_ng2d(3,intbound2d_ng,2,nz,intbound2d_single)
         s3=s3+ng_const
      elseif(i.eq.s4) then
         write(*,*) '----------------------storing source fct. at step n----------------------------'
         write(*,*)
         call store_ng2d(4,intbound2d_ng,2,nz,intbound2d_single)
         s4=s4+ng_const

         
         if(opt_ng_cont) call ng_expol2d(intbound2d_ng,2,nz,intbound2d_single)
         if(opt_ait_cont) call ait_expol2d(intbound2d_ng,2,nz,intbound2d_single)
         intbound2d(:,:,oindx)=intbound2d_single
      endif
   endif
   
enddo
!
deallocate(alocont_nn2d_tmp)
deallocate(normalization2d_tmp)
deallocate(mint2d_tmp)
!
call cpu_time(te)
write(*,*)
write(*,*) 'total computation time', te-ts
write(*,*)
!
!-----------------set theoretical intensities (for plane-parallel atmosphere)--------------
!
if(n_z(oindx).gt.zero) then
   startz = 3
   endz = nz-2
   gamma =  1
elseif(n_z(oindx).lt.zero) then
   startz = nz-2
   endz = 3
endif

allocate(int2d_theo(nx,nz),stat=err)
int2d_theo=zero
int2d_theo(:,startz-gamma)=int2d(:,startz-gamma)
int2d_theo(:,startz-2*gamma)=int2d(:,startz-2*gamma)
!
i=nx/2+1
do k=startz, endz, gamma
   dtau = half*(opac2d(i,k)+opac2d(i,k-gamma))*(z(k)-z(k-gamma))/n_z(oindx)
   int2d_theo(:,k)=int2d_theo(:,k-gamma)*exp(-dtau)
enddo
!
end subroutine benchmark01_solution
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_alo2d
!
use prog_type
use fund_const
use dime2d, only: int2d, alocont_o_nn2d, alocont_nn2d, mint2d, normalization2d, &
                  alocont_nn2d_tmp, mint2d_tmp, normalization2d_tmp, bnue2d, t2d, &
                  nx, nz, imask2d, scont2d, imaskb2d
use angles, only: dim_omega, n_z, n_x
use freq, only: nodes_nue
use mod_debug, only: iindx, kindx
use options, only: opt_method
!
implicit none
!
integer(i4b) :: i, k, oindx
!
alocont_nn2d=zero
normalization2d=zero
mint2d=zero
!bnue2d=zero
!
if(allocated(int2d)) deallocate(int2d)
if(allocated(alocont_o_nn2d)) deallocate(alocont_o_nn2d)
if(allocated(alocont_nn2d_tmp)) deallocate(alocont_nn2d_tmp)
if(allocated(normalization2d_tmp)) deallocate(normalization2d_tmp)
if(allocated(mint2d_tmp)) deallocate(mint2d_tmp)
!
allocate(int2d(nx,nz))
allocate(alocont_o_nn2d(nx,nz,27))
allocate(alocont_nn2d_tmp(nx,nz,27))
allocate(normalization2d_tmp(nx,nz))
allocate(mint2d_tmp(nx,nz))

alocont_nn2d_tmp=zero
normalization2d_tmp=zero
mint2d_tmp=zero

do i=nx-1, nx-1!nx
   do k=1, nz
!do i=1,1
!   do k=19,19!nz!20,20
      scont2d=zero
      alocont_nn2d_tmp=zero
      mint2d_tmp=zero
      normalization2d_tmp=zero
      iindx=i
      kindx=k
      select case(imaskb2d(i,k))
         case(5,6,9)
            scont2d(i,k)=one
         case(3,7)
            scont2d(1,k)=one
            scont2d(nx-1,k)=one
         case(4,8)
            scont2d(2,k)=one
            scont2d(nx,k)=one
         case(1,2)
            scont2d(i,k)=one
      end select

      do oindx=dim_omega, 1, -1
         if(opt_method.eq.4) then
            call fsc_cont2d_lin(oindx,1)
         elseif(opt_method.eq.5) then
            call fsc_cont2d(oindx,1)
         endif
         select case(imaskb2d(i,k))
!            case(0) !at k=0 or k=nz
!               write(*,*) 14, i, k, oindx, int2d(i,k), alocont_o_nn2d(i,k,14), alocont_o_nn2d(i,k,14)-int2d(i,k)
!               write(*,*)
!            case(1,2) !at k=1 or k=nz-1
!               write(*,*) 5, i, k, oindx, int2d(i,k-1), alocont_o_nn2d(i,k,5), alocont_o_nn2d(i,k,5)-int2d(i,k-1)
!               write(*,*) 14, i, k, oindx, int2d(i,k), alocont_o_nn2d(i,k,14), alocont_o_nn2d(i,k,14)-int2d(i,k)
!               write(*,*) 23, i, k, oindx, int2d(i,k+1), alocont_o_nn2d(i,k,23), alocont_o_nn2d(i,k,23)-int2d(i,k+1)
!               write(*,*)
            case(3) !at i=1
               write(*,*) 4, i, k, oindx, int2d(nx-2,k-1), alocont_o_nn2d(nx-1,k,4), alocont_o_nn2d(nx-1,k,4)-int2d(nx-2,k-1)
               write(*,*) 5, i, k, oindx, int2d(nx-1,k-1), alocont_o_nn2d(nx-1,k,5), alocont_o_nn2d(nx-1,k,5)-int2d(nx-1,k-1)
               write(*,*) 6, i, k, oindx, int2d(i+1,k-1), alocont_o_nn2d(i,k,6), alocont_o_nn2d(i,k,6)-int2d(i+1,k-1)
               write(*,*) 13, i, k, oindx, int2d(nx-2,k), alocont_o_nn2d(nx-1,k,13), alocont_o_nn2d(nx-1,k,13)-int2d(nx-2,k)               
               write(*,*) 14, i, k, oindx, int2d(i,k), alocont_o_nn2d(i,k,14), alocont_o_nn2d(i,k,14)-int2d(i,k), n_x(oindx), n_z(oindx)
               write(*,*) 15, i, k, oindx, int2d(i+1,k), alocont_o_nn2d(i,k,15), alocont_o_nn2d(i,k,15)-int2d(i+1,k)
               write(*,*) 22, i, k, oindx, int2d(nx-2,k+1), alocont_o_nn2d(nx-1,k,22), alocont_o_nn2d(nx-1,k,22)-int2d(nx-2,k+1)
               write(*,*) 23, i, k, oindx, int2d(i,k+1), alocont_o_nn2d(i,k,23), alocont_o_nn2d(i,k,23)-int2d(i,k+1)
               write(*,*) 24, i, k, oindx, int2d(i+1,k+1), alocont_o_nn2d(i,k,24), alocont_o_nn2d(i,k,24)-int2d(i+1,k+1)
               write(*,*)
!            case(8) !at i=nx
!               write(*,*) 4, i, k, oindx, int2d(i-1,k-1), alocont_o_nn2d(i,k,4), alocont_o_nn2d(i,k,4)-int2d(i-1,k-1)
!               write(*,*) 5, i, k, oindx, int2d(i,k-1), alocont_o_nn2d(i,k,5), alocont_o_nn2d(i,k,5)-int2d(i,k-1)
!               write(*,*) 13, i, k, oindx, int2d(i-1,k), alocont_o_nn2d(i,k,13), alocont_o_nn2d(i,k,13)-int2d(i-1,k)
!               write(*,*) 14, i, k, oindx, int2d(i,k), alocont_o_nn2d(i,k,14), alocont_o_nn2d(i,k,14)-int2d(i,k)
!               write(*,*) 22, i, k, oindx, int2d(i-1,k+1), alocont_o_nn2d(i-1,k,22), alocont_o_nn2d(i,k,22)-int2d(i-1,k+1)               
!               write(*,*) 23, i, k, oindx, int2d(i,k+1), alocont_o_nn2d(i,k,23), alocont_o_nn2d(i,k,23)-int2d(i,k+1)
!               write(*,*)
!            case(4,5,6,7,9) !standard alo
            case(7) !standard alo
               write(*,*) 4, i, k, oindx, int2d(i-1,k-1), alocont_o_nn2d(i,k,4), alocont_o_nn2d(i,k,4)-int2d(i-1,k-1)
               write(*,*) 5, i, k, oindx, int2d(i,k-1), alocont_o_nn2d(i,k,5), alocont_o_nn2d(i,k,5)-int2d(i,k-1)
               write(*,*) 6, i, k, oindx, int2d(i+1,k-1), alocont_o_nn2d(i,k,6), alocont_o_nn2d(i,k,6)-int2d(i+1,k-1)
               write(*,*) 13, i, k, oindx, int2d(i-1,k), alocont_o_nn2d(i,k,13), alocont_o_nn2d(i,k,13)-int2d(i-1,k)
               write(*,*) 14, i, k, oindx, int2d(i,k), alocont_o_nn2d(i,k,14), alocont_o_nn2d(i,k,14)-int2d(i,k), n_x(oindx), n_z(oindx)
               write(*,*) 15, i, k, oindx, int2d(i+1,k), alocont_o_nn2d(i,k,15), alocont_o_nn2d(i,k,15)-int2d(i+1,k)
               write(*,*) 22, i, k, oindx, int2d(i-1,k+1), alocont_o_nn2d(i,k,22), alocont_o_nn2d(i,k,22)-int2d(i-1,k+1)
               write(*,*) 23, i, k, oindx, int2d(i,k+1), alocont_o_nn2d(i,k,23), alocont_o_nn2d(i,k,23)-int2d(i,k+1)
               write(*,*) 24, i, k, oindx, int2d(i+1,k+1), alocont_o_nn2d(i,k,24), alocont_o_nn2d(i,k,24)-int2d(i+1,k+1)
               write(*,*)
               
!               if(abs(int2d(i,k)-alocont_o_nn2d(i,k,14)).gt.small_number) stop 'error in check alo2d: 14'
!               if(abs(int2d(i-1,k)-alocont_o_nn2d(i,k,13)).gt.small_number) stop 'error in check alo2d: 13'
!               if(abs(int2d(i+1,k)-alocont_o_nn2d(i,k,15)).gt.small_number) stop 'error in check alo2d: 15'               
!               if(abs(int2d(i,k-1)-alocont_o_nn2d(i,k,5)).gt.small_number) stop 'error in check alo2d: 5'
!               if(abs(int2d(i,k+1)-alocont_o_nn2d(i,k,23)).gt.small_number) stop 'error in check alo2d: 23'               
!               if(abs(int2d(i-1,k-1)-alocont_o_nn2d(i,k,4)).gt.small_number) stop 'error in check alo2d: 4'
!               if(abs(int2d(i+1,k+1)-alocont_o_nn2d(i,k,24)).gt.small_number) stop 'error in check alo2d: 24'
!               if(abs(int2d(i-1,k+1)-alocont_o_nn2d(i,k,22)).gt.small_number) stop 'error in check alo2d: 22'
!               if(abs(int2d(i+1,k-1)-alocont_o_nn2d(i,k,6)).gt.small_number) stop 'error in check alo2d: 6'
            case default
         end select
      enddo
   enddo
enddo


stop 'go on in check_alo2d'


end subroutine check_alo2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine testlambs_direct
!
use prog_type
use fund_const
use dime2d, only: int2d, alocont_o_nn2d, alocont_nn2d, mint2d, normalization2d, &
                  alocont_nn2d_tmp, mint2d_tmp, normalization2d_tmp, bnue2d, t2d, &
                  nx, nz, imask2d, scont2d, imaskb2d, z
use angles, only: dim_omega, n_z
use freq, only: nodes_nue
use mod_debug, only: iindx, kindx
use options, only: opt_method
use params_input, only: eps_cont
use iter, only: itmaxc, devmaxc, epsmaxc_arr, nconvc
use mod_math, only: conv_indx_2d_to_1d, conv_indx_1d_to_2d, calc_dev, bnue

!
implicit none
!
!local scalars
integer(i4b) :: nit, j, i, ii, k, kk, oindx, indx_1d_col, indx_1d_row, err
real(dp) :: eps_max, fdum1, fdum2, fdum3, scont_new, wor, sum, sum2
!
!local arrays
real(dp), dimension(:,:), allocatable :: lambda_mat
real(dp), dimension(:), allocatable :: phi_boundary
real(dp), dimension(:,:), allocatable :: a_mat, unity_mat, eps_mat, zeta_mat
real(dp), dimension(:), allocatable :: bnue_vec
real(dp), dimension(:), allocatable :: scont_vec, b_vec, mint_vec, eps_vec
!
!local functions

!
allocate(lambda_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
allocate(a_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
allocate(zeta_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
allocate(unity_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
allocate(eps_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
allocate(bnue_vec(nx*nz-2*nz), stat=err)
allocate(scont_vec(nx*nz-2*nz), stat=err)
allocate(b_vec(nx*nz-2*nz), stat=err)
allocate(phi_boundary(nx*nz-2*nz), stat=err)
allocate(mint_vec(nx*nz-2*nz), stat=err)
allocate(eps_vec(nx*nz-2*nz), stat=err)
!
!
!prepare bnue array
do i=1, nx
   do k=1, nz
!      bnue2d(i,k) = one!bnue(nodes_nue(1), t2d(i,k))
      bnue2d(i,k) = bnue(nodes_nue(1), t2d(i,k))
   enddo
enddo
!
alocont_nn2d=zero
normalization2d=zero
mint2d=zero
!
if(allocated(int2d)) deallocate(int2d)
if(allocated(alocont_o_nn2d)) deallocate(alocont_o_nn2d)
if(allocated(alocont_nn2d_tmp)) deallocate(alocont_nn2d_tmp)
if(allocated(normalization2d_tmp)) deallocate(normalization2d_tmp)
if(allocated(mint2d_tmp)) deallocate(mint2d_tmp)
!
allocate(int2d(nx,nz))
allocate(alocont_o_nn2d(nx,nz,27))
allocate(alocont_nn2d_tmp(nx,nz,27))
allocate(normalization2d_tmp(nx,nz))
allocate(mint2d_tmp(nx,nz))

alocont_nn2d_tmp=zero
normalization2d_tmp=zero
mint2d_tmp=zero

goto 40
!
!--------------------calculate lambda matrix-----------------------
!
10 continue
!note: set boundary to zero
do i=1, nx-2
   do k=1, nz

     alocont_nn2d_tmp=zero
     normalization2d_tmp=zero
     mint2d_tmp=zero
     scont2d=zero
!     scont2d(i,k)=one
     call conv_indx_2d_to_1d(i,k, nx-2, indx_1d_col)      
!
     select case(imaskb2d(i,k))
         case(5,6,9)
            scont2d(i,k)=one
         case(3,7)
            scont2d(1,k)=one
            scont2d(nx-1,k)=one
!            scont2d(i,k)=one            
         case(4,8)
            scont2d(2,k)=one
            scont2d(nx,k)=one
!            scont2d(i,k)=one            
         case(1,2)
            scont2d(i,k)=one
            if(i.eq.1) scont2d(nx-1,k)=one
            if(i.eq.2) scont2d(nx,k)=one
      end select

      do oindx=1, dim_omega
         if(opt_method.eq.4) then
            call fsc_cont2d_lin(oindx,1)
         elseif(opt_method.eq.5) then
            call fsc_cont2d(oindx,1)
         endif
      enddo

!      write(*,*) i, k, mint2d_tmp(i,k), alocont_nn2d_tmp(i,k,14), normalization2d_tmp(i,k)
!
!store corresponding matrix elements      
      do ii=1, nx-2
         do kk=1, nz
            call conv_indx_2d_to_1d(ii,kk, nx-2, indx_1d_row)
            lambda_mat(indx_1d_row,indx_1d_col)=mint2d_tmp(ii,kk)
         enddo
      enddo
   enddo
enddo

open(1,file='TRASH/lambda_mat', form='unformatted')
   write(1) lambda_mat
close(1)
!
stop 'lambda matrix calculated, returning'
!
!--------------------calculate boundary contribution---------------
!
20 continue
!note: set source function to zero, and only boundary intensities
scont2d=zero

do oindx=1, dim_omega
   if(opt_method.eq.4) then
      call fsc_cont2d_lin(oindx,1)
   elseif(opt_method.eq.5) then
      call fsc_cont2d(oindx,1)
   endif
enddo
      
do ii=1, nx-2
   do kk=1, nz
      call conv_indx_2d_to_1d(ii,kk, nx-2, indx_1d_row)
      phi_boundary(indx_1d_row)=mint2d_tmp(ii,kk)
      bnue_vec(indx_1d_row)=bnue2d(ii,kk)
      if(kk.eq.nz) bnue_vec(indx_1d_row)=zero
      if(kk.eq.nz-1) bnue_vec(indx_1d_row)=zero
      if(kk.eq.1) phi_boundary(indx_1d_row)=bnue2d(ii,kk)
      if(kk.eq.2) phi_boundary(indx_1d_row)=bnue2d(ii,kk)
!      write(*,*) ii, kk, bnue_vec(indx_1d_row), bnue2d(ii,kk)
   enddo
enddo
!
open(1,file='TRASH/phi_boundary', form='unformatted')
   write(1) phi_boundary
close(1)
!
open(1,file='TRASH/bnue_vec', form='unformatted')
   write(1) bnue_vec
close(1)
!
stop 'boundary contribution calculated, returning'
!
!------------------set up matrix system----------------------------
!
30 continue
!
!
!
open(1,file='TRASH/lambda_mat', form='unformatted')
   read(1) lambda_mat
close(1)
open(1,file='TRASH/phi_boundary', form='unformatted')
   read(1) phi_boundary
close(1)
open(1,file='TRASH/bnue_vec', form='unformatted')
   read(1) bnue_vec
close(1)
!
!
!
eps_mat=zero
unity_mat=zero
do i=1, nx*nz-2*nz
   eps_mat(i,i)=eps_cont
   unity_mat(i,i)=one
enddo
zeta_mat=zero
zeta_mat=unity_mat-eps_mat
!
a_mat=matmul(zeta_mat,lambda_mat)
a_mat=unity_mat-a_mat
!
bnue_vec=matmul(eps_mat,bnue_vec)
phi_boundary=matmul(zeta_mat,phi_boundary)
b_vec=phi_boundary+bnue_vec
!
!note: bnue_vec here serves only as dummy argument
call dgesv(nx*nz-2*nz, 1, a_mat, nx*nz-2*nz, bnue_vec, b_vec, nx*nz-2*nz, err)
scont_vec=b_vec
!if(err.ne.0) stop 'error in testlambs_direct: no inversion possible'
!
!transform solution to 2d array
do i=1, nx*nz-2*nz
   call conv_indx_1d_to_2d (i, nx-2, ii, kk)
   scont2d(ii,kk) = scont_vec(i)
enddo
!
do k=1, nz
   write(*,*) z(k), scont2d(1,k), scont2d(nx/2+1,k), scont2d(nx,k), scont2d(nx/2+1,k)/bnue2d(1,1)
enddo
!
!dummy output
epsmaxc_arr=one
nconvc=10

write(*,*) bnue2d(1,:)
call output

stop 'testlambs_direct done, returning'
!
!------------------solve iteratively-------------------------------
!
40 continue
!
open(1,file='TRASH/lambda_mat', form='unformatted')
   read(1) lambda_mat
close(1)
open(1,file='TRASH/phi_boundary', form='unformatted')
   read(1) phi_boundary
close(1)
open(1,file='TRASH/bnue_vec', form='unformatted')
   read(1) bnue_vec
close(1)
!
!check the lambda-matrix
do i=1, nx*nz-2*nz
   sum=zero
   do j=1, nx*nz-2*nz
      sum=sum+lambda_mat(i,j)
   enddo
   call conv_indx_1d_to_2d (i, nx-2, ii, kk)
   write(*,*) i, ii, kk, sum, lambda_mat(i,i)
enddo
!stop
!
!start value
scont_vec=zero
do i=1, nx-2
   do k=1, nz
      call conv_indx_2d_to_1d(i,k,nx-2,indx_1d_row)
      scont_vec(indx_1d_row)=bnue2d(i,k)
   enddo
enddo
!

!stop
ii=1
kk=17
call conv_indx_2d_to_1d(ii,kk, nx-2, indx_1d_row)
sum=zero
sum2=zero
do j=1, nx*nz-2*nz
   call conv_indx_1d_to_2d (j, nx-2, ii, kk)
   write(*,*) ii, kk, lambda_mat(indx_1d_row,j), scont_vec(j), lambda_mat(indx_1d_row,j)*scont_vec(j), phi_boundary(j)
   sum=sum+lambda_mat(indx_1d_row,j)
   sum2=sum2+lambda_mat(indx_1d_row,j)*scont_vec(j)
enddo
write(*,*) 'sum', sum, sum2, bnue2d(ii,kk), phi_boundary(indx_1d_row), sum*bnue2d(ii,kk)
!stop

!
!
mint_vec=zero
eps_vec=zero
!
!for underestimation of diagonal
wor=one!.99d0
!
do nit=1, itmaxc
   mint_vec=matmul(lambda_mat,scont_vec)+phi_boundary

   call calc_dev(eps_vec, mint_vec, nx*nz-2*nz, eps_max)
   eps_vec=mint_vec
   write(*,*) 'maximum deviation', nit, eps_max
   epsmaxc_arr(nit)=eps_max
   nconvc=nit
   if(abs(eps_max).lt.devmaxc) exit
   
   do i=1, nx*nz-2*nz
      fdum1=(one-eps_cont)
      fdum2=one-fdum1*lambda_mat(i,i)
      scont_new=fdum1/fdum2*mint_vec(i)-fdum1/fdum2*lambda_mat(i,i)*scont_vec(i) + eps_cont/fdum2*bnue_vec(i)

      fdum3=one/(one/fdum1-lambda_mat(i,i))
      scont_new=fdum3*mint_vec(i)-fdum3*lambda_mat(i,i)*scont_vec(i) + eps_cont/fdum2*bnue_vec(i)
      call conv_indx_1d_to_2d (i, nx-2, ii, kk)

!      write(*,*) ii, kk, fdum1/fdum2, fdum3, lambda_mat(i,i), mint_vec(i), phi_boundary(i), scont_new/bnue2d(ii,kk), one-lambda_mat(i,i)
!      if(scont_new/bnue2d(ii,kk).gt.one) write(*,*) ii, kk, scont_new, scont_new/bnue2d(ii,kk), mint_vec(i), bnue2d(ii,kk), lambda_mat(i,i)
!      if(scont_new/bnue2d(ii,kk).gt.one) stop 'error in bla'
      scont_vec(i)=scont_new
      
   enddo

   do i=1, nx*nz-2*nz
      call conv_indx_1d_to_2d (i, nx, ii, kk)
      if(ii.eq.nx/2+1) write(*,*) ii, kk, scont_vec(i), mint_vec(i), phi_boundary(i)
   enddo

!   if(nit.ge.1) stop 'go on in testlambs_direct, step 40'
   
enddo
!
!transform solution to 2d array
do i=1, nx*nz-2*nz
   call conv_indx_1d_to_2d (i, nx-2, ii, kk)
   scont2d(ii,kk) = scont_vec(i)
enddo
!
do k=1, nz
   write(*,*) z(k), scont2d(1,k), scont2d(nx/2+1,k), scont2d(nx,k)
enddo
!
!dummy output
call output

stop 'testlambs_direct iteratively done, returning'
!
!
end subroutine testlambs_direct
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_grid2d
!
!test the upwind and downwind delta-tau steps
!such that alo-coefficients are not 'extreme'  
!
use prog_type
use fund_const
use dime2d, only: nx, nz, x, z, opac2d
use options, only: opt_method
use mod_integ1d, only: integ1d_tau_ud
use mod_interp1d, only: interpol_typ_quad2b, interpol_typ_quad2
use mod_math, only: calc_dtau_crit
!
implicit none
!
! ... local scalars
integer(i4b) :: iim2, iim1, i, iip1, kkm2, kkm1, k, kkp1
real(dp) :: delt_u, delt_d
real(dp) :: dels_xu, dels_zu, dels_u, dels_xd, dels_zd, dels_d
real(dp) :: nn_x, nn_z
real(dp) :: x_u, z_u, opac_u, opac_p, x_d, z_d, opac_d
real(dp) :: delt_crita, delt_critb, delt_critc
!
! ... local funcions
!
!linear interpolations give always stable alo-coefficients
if(opt_method.eq.4) return
!
do i=3, nx-2
   do k=3, nz-2
!
      opac_p=opac2d(1,k)
!
!-------------------check (n_x,n_z)=(0,1) direction---------------------
!      
      opac_u=opac2d(i,k-1)
      opac_d=opac2d(i,k+1)
      dels_u=(z(k)-z(k-1))
      dels_d=(z(k+1)-z(k))
      call integ1d_tau_ud(dels_u, dels_d, opac_u, opac_p, opac_d, delt_u, delt_d)
      call calc_dtau_crit(delt_u, delt_crita, delt_critb, delt_critc)
!      write(*,*) k, delt_u, delt_d, delt_crita, delt_critb, delt_critc
      if(delt_d.gt.delt_crita) stop 'error in check_grid1d: delt_d > delt_crita => increase resolution'
      if(delt_d.lt.delt_critb) stop 'error in check_grid1d: delt_d < delt_critb => decrease resolution'
      if(delt_d.lt.delt_critc) stop 'error in check_grid1d: delt_d < delt_critc => decrease resolution'
!
!-------------------check (n_x,n_z)=(1,0) direction---------------------
!      
      opac_u=opac2d(i-1,k)
      opac_d=opac2d(i+1,k)
      dels_u=(x(i)-x(i-1))
      dels_d=(x(i+1)-x(i))
      call integ1d_tau_ud(dels_u, dels_d, opac_u, opac_p, opac_d, delt_u, delt_d)
      call calc_dtau_crit(delt_u, delt_crita, delt_critb, delt_critc)
!      write(*,*) i, k, delt_u, delt_d, delt_crita, delt_critb, delt_critc
      if(delt_d.gt.delt_crita) stop 'error in check_grid1d: delt_d > delt_crita => increase resolution'
      if(delt_d.lt.delt_critb) stop 'error in check_grid1d: delt_d < delt_critb => decrease resolution'
      if(delt_d.lt.delt_critc) stop 'error in check_grid1d: delt_d < delt_critc => decrease resolution'
!
!-------------------check (n_x,n_z)=(1,1) direction---------------------
!
      nn_x=one/sqrt(two)
      nn_z=one/sqrt(two)

      iim1=i-1
      kkm1=k-1
      iip1=i+1
      kkp1=k+1
!
!calculate distance to intersection point with upwind x-grid-coordinate
      dels_xu=(x(i)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
      dels_zu=(z(k)-z(kkm1))/nn_z
      dels_u=min(dels_xu,dels_zu)
!
!calculate distance to intersection point with downwind x-grid-coordinate
      dels_xd=(x(iip1)-x(i))/nn_x
!calculate distance to instersection point with downwind z-grid-coordinate
      dels_zd=(z(kkp1)-z(k))/nn_z
      dels_d=min(dels_xd,dels_zd)
!
!upwind point
      if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
         x_u = x(i) - dels_u*nn_x
         z_u = z(k) - dels_u*nn_z
         kkm2=k-2
!standard version
         opac_u=interpol_typ_quad2b(opac2d(iim1,kkm2),opac2d(iim1,kkm1),opac2d(iim1,k), &
                                    z(kkm2), z(kkm1), z(k), z_u)
      elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
         x_u=x(i)-dels_u*nn_x
         z_u=z(k)-dels_u*nn_z
         iim2=i-2
         opac_u=interpol_typ_quad2b(opac2d(iim2,kkm1), opac2d(iim1,kkm1), opac2d(i,kkm1), &
                                    x(iim2), x(iim1), x(i), x_u)
      else
         write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
         stop 'error in check_grid2d: invalid dels_u'
      endif
!
!downwind point
!
      if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
         x_d=x(i)+dels_d*nn_x
         z_d=z(k)+dels_d*nn_z
!
         opac_d=interpol_typ_quad2b(opac2d(iip1,kkm1), opac2d(iip1,k), opac2d(iip1,kkp1), &
                                    z(kkm1), z(k), z(kkp1), z_d)
      elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
         x_d=x(i)+dels_d*nn_x
         z_d=z(k)+dels_d*nn_z
!
         opac_d=interpol_typ_quad2b(opac2d(iim1,kkp1), opac2d(i,kkp1), opac2d(iip1,kkp1), &
                                  x(iim1), x(i), x(iip1), x_d)
      else
         write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
         stop 'error in check_grid2d: invalid dels_d'
      endif
!     
      call integ1d_tau_ud(dels_u, dels_d, opac_u, opac_p, opac_d, delt_u, delt_d)
      call calc_dtau_crit(delt_u, delt_crita, delt_critb, delt_critc)
!      write(*,*) i, k, delt_u, delt_d, delt_crita, delt_critb, delt_critc
      if(delt_d.gt.delt_crita) stop 'error in check_grid1d: delt_d > delt_crita => increase resolution'
      if(delt_d.lt.delt_critb) stop 'error in check_grid1d: delt_d < delt_critb => decrease resolution'
      if(delt_d.lt.delt_critc) stop 'error in check_grid1d: delt_d < delt_critc => decrease resolution'

      
   enddo
enddo
!
!
!do i=1, nz
!   delt_u = 0.5d0*opac2d(1,i)*(x(nx)-x(1))
!   write(*,*) i, z(i), delt_u
!enddo
!stop
!
!
end subroutine check_grid2d

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_nodes_omega
!
!check if angular grid is reasonable:
!check for: opposite directions
!           normalization  
!
use prog_type
use fund_const
use angles, only: dim_omega, n_x, n_z, weight_omega
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j
real(dp) :: norm, norm_vec, opposite
!
! ... local logicals
logical :: lopposite
!
! ... local functions
!
write(*,*) '-------------------------checking the angular grid-----------------------------'
write(*,*)
!
!check normalization
!check if zero directions occurr
!check if vector is normalized
norm=0.d0
do i=1, dim_omega
!   write(*,*) n_x(i), n_z(i), weight_omega(i)
   norm=norm+weight_omega(i)
!
   norm_vec = n_x(i)**2 + n_z(i)**2
   if(abs(norm_vec-one).gt.1.d-15) stop 'error in check_nodes_omega: n_x^2+n_z^2 < 1 not allowed'
!
   if(n_x(i).eq.zero) stop 'error in check_nodes_omega: n_x = 0 not allowed'
   if(n_z(i).eq.zero) stop 'error in check_nodes_omega: n_z = 0 not allowed'   
!
enddo
if(abs(norm-one).gt.1.d-15) stop 'error in check_nodes_omega: normalization not one'
!
!
!
!check opposite directions
do i=1, dim_omega
   lopposite=.false.
   do j=1, dim_omega
      opposite = n_x(i)*n_x(j)+n_z(i)*n_z(j)
!      write(*,*) n_x(i), n_z(i), n_x(j), n_z(j), opposite, opposite+one
      if(abs(opposite+one).lt.1.d-15) then
         lopposite=.true.
         exit
      endif
   enddo
   if(.not.lopposite) stop 'error in check_nodes_omega: no opposite direction found'
enddo
!

!stop 'go on in check_nodes_omega'
!
end subroutine check_nodes_omega
