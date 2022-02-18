!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
program prog
!
use prog_type
use options, only: opt_method, opt_alo_cont, opt_ng_cont
use dime1d, only: nz, z, tau1d, bnue1d, scont1d, mint1d, scont1d_diff, &
                  scont1d_old, alo_diag, alo_subdiag, alo_supdiag, opac1d, t1d, imask1d
use iter, only: itmaxc, epsmaxc_arr, solm0, solm1, solm2, solm3, devmaxc
use params_input, only: eps_cont
use freq, only: xnue0
use mod_math, only: bnue
!
implicit none
!
! ... local scalars
integer(i4b) :: i, nconv, err
real(dp), parameter :: mu=1.d0/sqrt(3.d0)
real(dp) :: eps_max
integer(i4b) :: s1, s2, s3, s4, s5
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
!create z-grid
call gridz
!
!create model
call setup_mod1d
!
call setup_opacities1
!
!-----------------------------------------------------------------------
!
allocate(bnue1d(nz),stat=err)
do i=1, nz
   bnue1d(i) = bnue(xnue0,t1d(i))
enddo
!
!-----------------------------------------------------------------------
!
allocate(scont1d(nz),stat=err)
allocate(scont1d_old(nz),stat=err)
allocate(scont1d_diff(nz),stat=err)
allocate(mint1d(nz),stat=err)
scont1d_diff=0.d0
!
call calc_startval_cont
!
!-----------------------------------------------------------------------
!
!solution with diffusion equation
call diffusion(nz, eps_cont, bnue1d, tau1d, scont1d_diff)
!
!-----------------------------------------------------------------------
!
!solution with pp radiative transfer
!
!allocate all array
allocate(alo_diag(nz),stat=err)
allocate(alo_subdiag(nz),stat=err)
allocate(alo_supdiag(nz),stat=err)
allocate(solm0(nz),stat=err)
allocate(solm1(nz),stat=err)
allocate(solm2(nz),stat=err)
allocate(solm3(nz),stat=err)
!
s1=1
s2=2
s3=3
s4=4
s5=5
!
epsmaxc_arr=0.d0
!
!scont1d(nz-1)=scont1d_diff(nz-1)
!scont1d(nz)=scont1d_diff(nz)
!
do i=1, itmaxc
!
   scont1d_old=scont1d
!
!calculate mean intensity
   select case(opt_method)
      case(0) 
         call fvm1st(nz, mu, z, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d, imask1d)
      case(1)
         stop 'fvm2dn to be adapted'
         call fvm2nd(nz, mu, z, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d)
      case(2)
         stop 'diff1st to be adapted'
         call diff1st(nz, mu, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d)
      case(3) 
         stop 'diff2nd to be adapted'
         call diff2nd(nz, mu, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d)
      case(4)
         call sc1st(nz, mu, z, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d, imask1d)
      case(5)
         call sc2nd(nz, mu, z, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d, imask1d)
      case default
         stop 'error: opt_method not valid'
   end select
!
!calculate new source function (classical lambda iteration)
   select case(opt_alo_cont)
      case(0) 
         call calc_snew_classical(nz, eps_cont, mint1d, bnue1d, scont1d, imask1d)
      case(1) 
         call calc_snew_diag(nz, eps_cont, mint1d, bnue1d, alo_diag, scont1d, imask1d)
      case(2) 
         call calc_snew_tridiag(nz, eps_cont, mint1d, bnue1d, alo_diag, alo_subdiag, alo_supdiag, scont1d, imask1d)
      case default
         stop 'error: opt_alo not valid'
   end select

!ng extrapolation
   if(opt_ng_cont) then
     if(i.eq.s1) then
        solm3=scont1d
         s1=s1+8
      elseif(i.eq.s2) then
         solm2=scont1d
         s2=s2+8
      elseif(i.eq.s3) then
         solm1=scont1d
         s3=s3+8
      elseif(i.eq.s4) then
         solm0=scont1d
         s4=s4+8
      endif
!
      if(i.eq.s5) then
         call ng_accel(nz, solm0, solm1, solm2, solm3, scont1d)
         s5=s5+8
      endif
   endif
!
!calculate error
   call calc_err(nz, scont1d_old, scont1d, eps_max)
   write(*,*) 'iteration', i, 'max relative error', eps_max
   epsmaxc_arr(i)=eps_max
   nconv=i
   if(eps_max.lt.devmaxc) exit
!
enddo
!
!-----------------------------------------------------------------------
!
open(1, file='outputFILES/diff1d/output.dat', form='formatted')
   write(1, '(a5, 6a20)') 'i', 'z', 'opac', 'tau', 'scont', 'scont_diff', 'bnue'
   do i=1, nz
      write(1,'(i5, 6es20.8)') i, z(i), opac1d(i), tau1d(i), scont1d(i), scont1d_diff(i), bnue1d(i)
   enddo
close(1)
!
open(1, file='outputFILES/diff1d/convergence.dat', form='formatted')
   write(1,'(a5,a20)') 'i', 'rel error'
   do i=1, nconv
      write(1,'(i5, es20.8)') i, epsmaxc_arr(i)
   enddo
close(1)
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
!
implicit none
!
! ... local scalars
real(dp) :: mstar
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
namelist / dimensions_3dz / nz, zmin, zmax
namelist / input_options / opt_ng_cont, opt_ait_cont, opt_alo_cont
namelist / input_diff1d / opt_method
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
!--------read model parameters for modelling continuum and line---------
!
   rewind 1
   read(1, nml=input_cont)
!
!--------read model parameters for method to be used--------------------
!
   rewind 1
   read(1, nml=input_diff1d)
!
close(1)
!

!
write(*,*) 'input opt_method: ', opt_method
write(*,*) '   0 - 1st order fvm'
write(*,*) '   1 - 2nd order fvm'
write(*,*) '   2 - 1st order differences'
write(*,*) '   3 - 2nd order differences'
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
!calculation region is 1
imask1d=1
!inner boundary is 0
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
use dime1d, only: nz, z, rho1d, opac1d, tau1d, zmax, zmin
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
real(dp) :: rho, dtau, grad, a, b, c
!
! ... local logicals
!
! ... local functions
real(dp) :: opac_thomson
!
allocate(opac1d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in setup_opacities1'
allocate(tau1d(nz),stat=err)
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
!exponentially decreasing opacity
b=log(tau_max/tau_min)/(zmax-zmin)
a=b*tau_max*exp(b*zmin)
opac1d=a*exp(-b*z)
tau1d=a/b*exp(-b*z)
!
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
!
end subroutine setup_opacities1
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
use freq, only: xnue0
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
!   write(*,'(8es20.8,3l5)') z(i), tau, opac1, opac2, r1, r2, dtau, thdepth
   if(tau.gt.thdepth) exit
enddo
!
zshift=1.d0-z(2)
!
r_therm=z(nz+1-indx_therm1d)+zshift
t_therm=t1d(nz+1-indx_therm1d)
i_core=bnue(xnue0, t_therm)
!
!-----------------------------------------------------------------------
!
do i=3, nz
   if(z(i).le.r_therm) then
!when thermalized, source function = planck function
      scont1d(i)=bnue(xnue0, t1d(i))
   else
!when not thermalized, calculate source function from
!mean intensity (i_core * dilution factor)
      rad=z(i)+zshift
      dilfac = r_therm**2 / rad**2
      dilfac = 0.5d0*(1.d0 - sqrt(1.d0 - dilfac))
      dum_bnue=bnue(xnue0, t1d(i))
      scont1d(i)=(1.d0-eps_cont)*dilfac*i_core + eps_cont*dum_bnue
   endif

enddo
scont1d(2)=bnue(xnue0,t1d(2))
scont1d(1)=bnue(xnue0,t1d(1))
!
mint1d=0.d0
!
scont1d=0.d0
scont1d(2)=bnue(xnue0,t1d(2))
scont1d(1)=bnue(xnue0,t1d(1))
!
end subroutine calc_startval_cont
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine ng_accel(nz, solm0, solm1, solm2, solm3, snew)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nz
real(dp), dimension(nz), intent(in) :: solm0, solm1, solm2, solm3
real(dp), dimension(nz) :: snew
!
! ... local scalars
integer(i4b) :: i
real(dp) :: a1, a2, b1, b2, c1, c2, a, b
!
a1=0.d0
a2=0.d0
b1=0.d0
b2=0.d0
c1=0.d0
c2=0.d0

do i=1, nz
   a1 = a1+(solm0(i)-2.d0*solm1(i)+solm2(i))**2.d0
   b1 = b1+(solm0(i)-solm1(i)-solm2(i)+solm3(i))*(solm0(i)-2.d0*solm1(i)+solm2(i))
   c1 = c1+(solm0(i)-solm1(i))*(solm0(i)-2.d0*solm1(i)+solm2(i))
   a2 = a2+(solm0(i)-solm1(i)-solm2(i)+solm3(i))*(solm0(i)-2.d0*solm1(i)+solm2(i))
   b2 = b2+(solm0(i)-solm1(i)-solm2(i)+solm3(i))**2.d0
   c2 = c2+(solm0(i)-solm1(i))*(solm0(i)-solm1(i)-solm2(i)+solm3(i))
enddo
!
!
a=(b2*c1-b1*c2)/(a1*b2-a2*b1)
b=(a1*c2-a2*c1)/(a1*b2-a2*b1)
!
snew=(1.d0-a-b)*solm0 + a*solm1+b*solm2
!
end subroutine ng_accel
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_snew_classical(nz, eps, mint, bnue, scont, imask)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nz
integer(i4b), dimension(nz), intent(in) :: imask
real(dp), intent(in) :: eps
real(dp), dimension(nz), intent(in) :: mint, bnue
real(dp), dimension(nz) :: scont
!
! ... local scalars
integer(i4b) :: i
!
do i=1, nz
   select case(imask(i))
      case(1) 
         scont(i) = (1.d0-eps)*mint(i) + eps*bnue(i)
      case default
   end select
enddo
!
end subroutine calc_snew_classical
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_snew_diag(nz, eps, mint, bnue, diag, scont, imask)
!
use prog_type
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nz
integer(i4b), dimension(nz), intent(in) :: imask
real(dp), intent(in) :: eps
real(dp), dimension(nz), intent(in) :: mint, bnue, diag
real(dp), dimension(nz) :: scont
!
! ... local scalars
integer(i4b) :: i
!
do i=1, nz
   select case(imask(i))
      case(1) 
         scont(i) = (1.d0-eps)/(1.d0-(1.d0-eps)*diag(i))*mint(i) - &
           (1.d0-eps)*diag(i)/(1.d0-(1.d0-eps)*diag(i))*scont(i) + &
           eps/(1.d0-(1.d0-eps)*diag(i))*bnue(i)
      case default
   end select
enddo
!
end subroutine calc_snew_diag
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_snew_tridiag(nz, eps, mint, bnue, diag, sub_diag, sup_diag, scont, imask)
!
use prog_type
use mod_math, only: invtri_lev  
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: eps
integer(i4b), dimension(nz), intent(in) :: imask
real(dp), dimension(nz), intent(in) :: mint, bnue, diag, sub_diag, sup_diag
real(dp), dimension(nz) :: scont
!
! ... local scalars
integer(i4b) :: i
!
! ... local arrays
real(dp), dimension(nz) :: rhs, diag2, sub_diag2, sup_diag2
!
!prepare tridiagonal system
do i=3, nz-2
   rhs(i) = sub_diag(i)*scont(i-1) + diag(i)*scont(i) + sup_diag(i)*scont(i+1)
enddo
rhs = (1.d0-eps)*(mint-rhs) + eps*bnue
rhs(1) = bnue(1)
rhs(2) = bnue(2)
rhs(nz) = scont(nz)
rhs(nz-1) = scont(nz-1)
!
diag2 = 1.d0-(1.d0-eps)*diag
sub_diag2 = -(1.d0-eps)*sub_diag
sup_diag2 = -(1.d0-eps)*sup_diag


!write(*,*) sub_diag2
!write(*,*)
!write(*,*) diag2
!write(*,*)
!write(*,*) sup_diag2
!write(*,*)
!write(*,*) rhs
!stop
call invtri_lev(nz, sub_diag2, diag2, sup_diag2, scont, rhs)
!
!
end subroutine calc_snew_tridiag
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine diffusion(nz, eps, bnue1d, tau1d, scont1d)
!
! solve diffusion equation
!
use prog_type
use params_input, only:eps_cont
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: eps
real(dp), dimension(nz), intent(in) :: bnue1d, tau1d
real(dp), dimension(nz) :: scont1d
!
! ... local scalars
integer(i4b) :: i
real(dp) :: del_tau1, del_tau2, del_tau3
!
! ... local arrays
real(dp), dimension(nz) :: sub_diag, diag, sup_diag, rhs, mint1d
!
!
!
!matrix elements
do i=3,nz-2
   del_tau1=tau1d(i)-tau1d(i-1)
   del_tau2=tau1d(i+1)-tau1d(i)
   del_tau3=tau1d(i+1)-tau1d(i-1)
!
   sub_diag(i) = -2.d0/(del_tau1*del_tau3)
   diag(i) = 3.d0*eps + 2.d0/(del_tau1*del_tau3) + 2.d0/(del_tau2*del_tau3)
   sup_diag(i) = -2.d0/(del_tau2*del_tau3)
enddo
!
!right hand side
do i=3,nz-2
   rhs(i)   = 3.d0*eps*bnue1d(i)
enddo
!
!boundary condition left (j=b for optically thick)
sub_diag(1) = 0.d0
diag(1) = 1.d0
sup_diag(1) = 0.d0
rhs(1)   = bnue1d(1)
sub_diag(2) = 0.d0
diag(2) = 1.d0
sup_diag(2) = 0.d0
rhs(2)   = bnue1d(2)
!
!boundary condition right J=(dJ/dtau)/sqrt(3)
del_tau1=tau1d(nz-1)-tau1d(nz-2)
sub_diag(nz-1) = 1.d0/del_tau1/sqrt(3.d0)
diag(nz-1) = 1.d0-1.d0/del_tau1/sqrt(3.d0)
sup_diag(nz-1) = 0.d0
rhs(nz-1)   = 0.d0
del_tau1=tau1d(nz)-tau1d(nz-1)
sub_diag(nz) = 1.d0/del_tau1/sqrt(3.d0)
diag(nz) = 1.d0-1.d0/del_tau1/sqrt(3.d0)
sup_diag(nz) = 0.d0
rhs(nz)   = 0.d0
!!
!boundary condition right S=0
!sub_diag(nz-1) = 0.d0
!diag(nz-1) = 1.d0
!sup_diag(nz-1) = 0.d0
!rhs(nz-1) = -eps_cont*bnue1d(nz-1)/(1.d0-eps_cont)
!sub_diag(nz) = 0.d0
!diag(nz) = 1.d0
!sup_diag(nz) = 0.d0
!rhs(nz) = 0.d0


!
!solve for mean intensity
call tridag(sub_diag,diag,sup_diag,rhs,mint1d,nz)
!
!calculate source function
!
do i=1,nz
   scont1d(i) = (1.d0-eps)*mint1d(i) + eps*bnue1d(i)
enddo
!
end subroutine diffusion
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_err(nz, x_old, x_new, err)
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: nz
real(dp), dimension(nz), intent(in) :: x_old, x_new
real(dp) :: err
!
! ... local scalars
real(dp), dimension(nz) :: dev1d
!
dev1d=0.d0
where(x_new.ne.0.d0) dev1d = (x_new-x_old)/x_new
dev1d=abs(dev1d)
!
err=maxval(dev1d)
!
end subroutine calc_err
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fvm1st(nz, mu, z, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d, imask1d)
!
!1d first order scheme
!
use prog_type
!
implicit none
!
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: mu
integer(i4b), dimension(nz), intent(in) :: imask1d
real(dp), dimension(nz), intent(in) :: bnue1d, z, scont1d, tau1d
real(dp), dimension(nz) :: mint1d, alo_diag, alo_subdiag, alo_supdiag
!
! ... local scalars
integer(i4b) :: i
real(dp) :: del_tau, del_tau2, dum, bcoeff, ccoeff
!
! ... local arrays
real(dp), dimension(nz) :: int1d_plus, int1d_minus, alo_diag_plus, alo_diag_minus, &
                           alo_subdiag_plus, alo_subdiag_minus, alo_supdiag_plus, alo_supdiag_minus
!real(dp), dimension(nz) :: int1d_plus2, int1d_minus2, alo_diag_plus2, &
!                           alo_subdiag_plus2, alo_subdiag_minus2, alo_supdiag_plus2, alo_supdiag_minus2
!
!positive direction (intensity)
int1d_plus = 0.d0
alo_diag_plus = 0.d0
alo_subdiag_plus = 0.d0
alo_supdiag_plus = 0.d0
!
int1d_plus(1) = bnue1d(1)
int1d_plus(2) = bnue1d(2)
!
do i=3, nz-2
   del_tau = (tau1d(i-1)-tau1d(i+1))/2.d0/mu   !along direction 1.d0/sqrt(3) to be consistent with diffusion
   bcoeff=del_tau/(1.d0+del_tau)
   ccoeff=1.d0/(1.d0+del_tau)
!
   int1d_plus(i) = bcoeff*scont1d(i) + ccoeff*int1d_plus(i-1)
!   
   alo_subdiag_plus(i) = ccoeff*alo_diag_plus(i-1)
   alo_diag_plus(i) = bcoeff
   alo_supdiag_plus(i) = 0.d0
!
enddo

!write(*,*) 'test'
!write(*,*) alo_subdiag_plus
!write(*,*)
!write(*,*) alo_diag_plus
!write(*,*)
!write(*,*) alo_supdiag_plus
!write(*,*)
!stop
!
!negative direction (intensity)
int1d_minus=0.d0
alo_diag_minus = 0.d0
alo_subdiag_minus = 0.d0
alo_supdiag_minus = 0.d0
!
do i=nz-2, 3, -1
   del_tau = (tau1d(i-1)-tau1d(i+1))/2.d0/mu   !along direction 1.d0/sqrt(3) to be consistent with diffusion
   bcoeff=del_tau/(1.d0+del_tau)
   ccoeff=1.d0/(1.d0+del_tau)
!
   int1d_minus(i) = bcoeff*scont1d(i) + ccoeff*int1d_minus(i+1)

   alo_supdiag_minus(i) = ccoeff*alo_diag_minus(i+1)*imask1d(i+1)
   alo_diag_minus(i) = bcoeff
   alo_subdiag_minus(i) = 0.d0
enddo
!write(*,*) alo_subdiag_minus
!write(*,*)
!write(*,*) alo_diag_minus
!write(*,*)
!write(*,*) alo_supdiag_minus
!write(*,*)
!stop
!
!angular integration
mint1d=0.5d0*(int1d_minus+int1d_plus)
alo_diag=0.5d0*(alo_diag_minus+alo_diag_plus)
alo_subdiag=0.5d0*(alo_subdiag_minus+alo_subdiag_plus)
alo_supdiag=0.5d0*(alo_supdiag_minus+alo_supdiag_plus)
!
!write(*,*) alo_subdiag
!write(*,*)
!write(*,*) alo_diag
!write(*,*)
!write(*,*) alo_supdiag
!stop
!
end subroutine fvm1st

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fvm2nd(nz, mu, z, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d)
!
!2nd order fvm
!
use prog_type
!
implicit none
!
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: mu
real(dp), dimension(nz), intent(in) :: bnue1d, z, scont1d, tau1d
real(dp), dimension(nz) :: mint1d, alo_diag, alo_subdiag, alo_supdiag
!
! ... local scalars
integer(i4b) :: i
real(dp) :: del_tau, del_tau2
!
! ... local arrays
real(dp), dimension(nz) :: int1d_plus, int1d_minus, alo_diag_plus, alo_diag_minus, &
                           alo_subdiag_plus, alo_subdiag_minus, alo_supdiag_plus, alo_supdiag_minus
!
!positive direction (intensity)
int1d_plus(1) = bnue1d(1)
do i=2, nz
   del_tau = (tau1d(i-1)-tau1d(i))/2.d0/mu !along direction 1.d0/sqrt(3) to be consistent with diffusion
   int1d_plus(i) = del_tau/(1.d0+del_tau)*scont1d(i) + del_tau/(1.d0+del_tau)*scont1d(i-1) + (1.d0-del_tau)/(1.d0+del_tau) * int1d_plus(i-1)
enddo
!positive direction (alo)
alo_diag_plus = 0.d0
alo_subdiag_plus = 0.d0
alo_supdiag_plus = 0.d0
del_tau = (tau1d(1)-tau1d(2))/2.d0/mu
alo_diag_plus(2) = del_tau/(1.d0+del_tau)
alo_subdiag_plus(2) = del_tau/(1.d0+del_tau)
do i=3, nz
   del_tau = (tau1d(i-1)-tau1d(i))/2.d0/mu
   del_tau2 = (tau1d(i-2)-tau1d(i-1))/2.d0/mu
   alo_diag_plus(i) = del_tau/(1.d0+del_tau)
   alo_subdiag_plus(i) = (1.d0-del_tau)*del_tau2/(1.d0+del_tau)/(1.d0+del_tau2) + del_tau/(1.d0+del_tau)
   alo_supdiag_plus(i) = 0.d0
enddo
!
!negative direction (intensity)
int1d_minus(nz)=0.d0
do i=nz-1, 1, -1
   del_tau = (tau1d(i)-tau1d(i+1))/2.d0/mu !along direction 1.d0/sqrt(3) to be consistent with diffusion
   int1d_minus(i) = del_tau/(1.d0+del_tau)*scont1d(i) + del_tau/(1.d0+del_tau)*scont1d(i+1) + (1.d0-del_tau)/(1.d0+del_tau) * int1d_minus(i+1)
enddo
!negative direction (alo)
alo_diag_minus=0.d0
alo_subdiag_minus=0.d0
alo_supdiag_minus=0.d0
del_tau = (tau1d(nz-1)-tau1d(nz))/2.d0/mu
alo_diag_minus(nz-1) = del_tau/(1.d0+del_tau)
do i=nz-2, 1, -1
   del_tau = (tau1d(i)-tau1d(i+1))/2.d0/mu
   del_tau2 = (tau1d(i+1)-tau1d(i+2))/2.d0/mu
   alo_diag_minus(i) = del_tau/(1.d0+del_tau)
   alo_subdiag_minus(i) = 0.d0
   alo_supdiag_minus(i) = (1.d0-del_tau)*del_tau2/(1.d0+del_tau)/(1.d0+del_tau2) + del_tau/(1.d0+del_tau)
enddo
!
!angular integration
mint1d=0.5d0*(int1d_minus+int1d_plus)
alo_diag=0.5d0*(alo_diag_minus+alo_diag_plus)
alo_subdiag=0.5d0*(alo_subdiag_minus+alo_subdiag_plus)
alo_supdiag=0.5d0*(alo_supdiag_minus+alo_supdiag_plus)
!
end subroutine fvm2nd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine diff1st(nz, mu, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d)
!
!first order difference method
!
use prog_type
!
implicit none
!
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: mu
real(dp), dimension(nz), intent(in) :: bnue1d, scont1d, tau1d
real(dp), dimension(nz) :: mint1d, alo_diag, alo_subdiag, alo_supdiag
!
! ... local scalars
integer(i4b) :: i
real(dp) :: del_tau, del_tau2
!
! ... local arrays
real(dp), dimension(nz) :: int1d_plus, int1d_minus, alo_diag_plus, alo_diag_minus, &
                           alo_subdiag_plus, alo_subdiag_minus, alo_supdiag_plus, alo_supdiag_minus
!
!positive direction (intensity)
int1d_plus(1) = bnue1d(1)
do i=2, nz
   del_tau = (tau1d(i-1)-tau1d(i))/mu !along direction 1.d0/sqrt(3) to be consistent with diffusion
   int1d_plus(i) = del_tau/(1.d0+del_tau)*scont1d(i) + int1d_plus(i-1)/(1.d0+del_tau)
enddo
!positive direction (alo)
alo_diag_plus = 0.d0
del_tau = (tau1d(1)-tau1d(2))/mu
alo_diag_plus(2) = del_tau/(1.d0+del_tau)
alo_subdiag_plus = 0.d0
alo_supdiag_plus = 0.d0
do i=3, nz
   del_tau = (tau1d(i-1)-tau1d(i))/mu
   del_tau2 = (tau1d(i-2)-tau1d(i))/mu
   alo_diag_plus(i) = del_tau/(1.d0+del_tau)
   alo_subdiag_plus(i) = del_tau2/(1.d0+del_tau2)/(1.d0+del_tau)
   alo_supdiag_plus(i) = 0.d0
enddo
!
!negative direction (intensity)
int1d_minus(nz)=0.d0
do i=nz-1, 1, -1
   del_tau = (tau1d(i)-tau1d(i+1))/mu !along direction 1.d0/sqrt(3) to be consistent with diffusion
   int1d_minus(i) = del_tau/(1.d0+del_tau)*scont1d(i) + int1d_minus(i+1)/(1.d0+del_tau)
enddo
!negative direction (alo)
alo_diag_minus=0.d0
alo_subdiag_minus=0.d0
alo_supdiag_minus=0.d0
del_tau = (tau1d(nz-1)-tau1d(nz))/mu
alo_diag_minus(nz-1) = del_tau/(1.d0+del_tau)
do i=nz-2, 1, -1
   del_tau = (tau1d(i)-tau1d(i+1))/mu
   del_tau2 = (tau1d(i+1)-tau1d(i+2))/mu
   alo_diag_minus(i) = del_tau/(1.d0+del_tau)
   alo_subdiag_minus(i) = 0.d0
   alo_supdiag_minus(i) = del_tau2/(1.d0+del_tau2)/(1.d0+del_tau)
enddo
!
!angular integration
mint1d=0.5d0*(int1d_minus+int1d_plus)
alo_diag=0.5d0*(alo_diag_minus+alo_diag_plus)
alo_subdiag=0.5d0*(alo_subdiag_minus+alo_subdiag_plus)
alo_supdiag=0.5d0*(alo_supdiag_minus+alo_supdiag_plus)
!
end subroutine diff1st

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine diff2nd(nz, mu, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d)
!
!second order difference method
!
use prog_type
use mod_math, only: invtri_lev  
!
implicit none
!
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: mu
real(dp), dimension(nz), intent(in) :: bnue1d, scont1d, tau1d
real(dp), dimension(nz) :: mint1d, alo_diag, alo_subdiag, alo_supdiag
!
! ... local scalars
integer(i4b) :: i
real(dp) :: del_tau
!
! ... local arrays
real(dp), dimension(nz) :: int1d_plus, int1d_minus, alo_diag_plus, alo_diag_minus, &
                           alo_subdiag_plus, alo_subdiag_minus, alo_supdiag_plus, alo_supdiag_minus
real(dp), dimension(nz) :: rhs, a, b, c, x
!
!no alo
alo_diag_plus=0.d0
alo_subdiag_plus=0.d0
alo_supdiag_plus=0.d0
alo_diag_minus=0.d0
alo_subdiag_minus=0.d0
alo_supdiag_minus=0.d0
!
!positive direction (intensity)
rhs=0.d0
rhs(1)=bnue1d(1)
rhs(2:nz)=scont1d(2:nz)
a(1)=1.d0
b(1)=0.d0
c(1)=0.d0
do i=2, nz-1
   del_tau=(tau1d(i+1)-tau1d(i-1))/mu
   a(i) = 1.d0
   b(i) = 1.d0/del_tau
   c(i) = -1.d0/del_tau
enddo
del_tau=(tau1d(nz)-tau1d(nz-1))/mu
a(nz) = 1.d0-1.d0/del_tau
b(nz) = 1.d0/del_tau
c(nz) = 0.d0

call invtri_lev(nz, b, a, c, int1d_plus, rhs)
!
!negative direction
rhs=0.d0
rhs(1:nz-1)=scont1d(1:nz-1)
del_tau = (tau1d(1)-tau1d(2))/mu
a(1) = 1.d0+1.d0/del_tau
b(1) = 0.d0
c(1) = -1.d0/del_tau
do i=2, nz-1
   del_tau=(tau1d(i+1)-tau1d(i-1))/mu
   a(i) = 1.d0
   b(i) = -1.d0/del_tau
   c(i) = 1.d0/del_tau
enddo
a(nz)=1.d0
b(nz)=0.d0
c(nz)=0.d0
!
call invtri_lev(nz, b, a, c, int1d_minus, rhs)
!
!angular integration
mint1d=0.5d0*(int1d_minus+int1d_plus)
alo_diag=0.5d0*(alo_diag_minus+alo_diag_plus)
alo_subdiag=0.5d0*(alo_subdiag_minus+alo_subdiag_plus)
alo_supdiag=0.5d0*(alo_supdiag_minus+alo_supdiag_plus)
!
end subroutine diff2nd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine sc1st(nz, mu, z, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d, imask1d)
!
!1d sc scheme with linear source function
!
use prog_type
!
implicit none
!
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: mu
integer(i4b), dimension(nz), intent(in) :: imask1d
real(dp), dimension(nz), intent(in) :: bnue1d, z, scont1d, tau1d
real(dp), dimension(nz) :: mint1d, alo_diag, alo_subdiag, alo_supdiag
!
! ... local scalars
integer(i4b) :: i
real(dp) :: del_tau, acoeff, bcoeff, ccoeff
!
! ... local arrays
real(dp), dimension(nz) :: int1d_plus, int1d_minus, alo_diag_plus, alo_diag_minus, &
                           alo_subdiag_plus, alo_subdiag_minus, alo_supdiag_plus, alo_supdiag_minus
!real(dp), dimension(nz) :: int1d_plus2, int1d_minus2, alo_diag_plus2, &
!                           alo_subdiag_plus2, alo_subdiag_minus2, alo_supdiag_plus2, alo_supdiag_minus2
!
!positive direction (intensity)
int1d_plus = 0.d0
alo_diag_plus = 0.d0
alo_subdiag_plus = 0.d0
alo_supdiag_plus = 0.d0
!
int1d_plus(1) = bnue1d(1)
int1d_plus(2) = bnue1d(2)
do i=3, nz-2
!   del_tau = (tau1d(i-1)-tau1d(i+1))/2.d0/mu   !along direction 1.d0/sqrt(3) to be consistent with diffusion
   del_tau = (tau1d(i-1)-tau1d(i))/mu   !along direction 1.d0/sqrt(3) to be consistent with diffusion
!
   acoeff=1.d0/del_tau-exp(-del_tau)-exp(-del_tau)/del_tau
   bcoeff=1.d0-1.d0/del_tau+exp(-del_tau)/del_tau
   ccoeff=exp(-del_tau)
!
   int1d_plus(i) = acoeff*scont1d(i-1) + bcoeff*scont1d(i) + ccoeff*int1d_plus(i-1)
!
   alo_subdiag_plus(i) = acoeff + ccoeff*alo_diag_plus(i-1)
   alo_diag_plus(i) = bcoeff
   alo_supdiag_plus(i) = 0.d0
enddo
!write(*,*) 'test'
!write(*,*) alo_subdiag_plus
!write(*,*)
!write(*,*) alo_diag_plus
!write(*,*)
!write(*,*) alo_supdiag_plus
!write(*,*)
!stop
!
!negative direction (intensity)
int1d_minus=0.d0
alo_diag_minus = 0.d0
alo_subdiag_minus = 0.d0
alo_supdiag_minus = 0.d0
!
!
do i=nz-2, 3, -1
!   del_tau = (tau1d(i-1)-tau1d(i+1))/2.d0/mu   !along direction 1.d0/sqrt(3) to be consistent with diffusion
   del_tau = (tau1d(i)-tau1d(i+1))/mu   !along direction 1.d0/sqrt(3) to be consistent with diffusion
   acoeff=1.d0/del_tau-exp(-del_tau)-exp(-del_tau)/del_tau
   bcoeff=1.d0-1.d0/del_tau+exp(-del_tau)/del_tau
   ccoeff=exp(-del_tau)
   int1d_minus(i) = acoeff*scont1d(i+1) + bcoeff*scont1d(i) + ccoeff*int1d_minus(i+1)

   alo_supdiag_minus(i) = acoeff + ccoeff*alo_diag_minus(i+1)
   alo_diag_minus(i) = bcoeff
   alo_subdiag_minus(i) = 0.d0
enddo
!alo_diag_minus(1)=1.d0
!alo_diag_minus(2)=1.d0
!write(*,*) 'test'
!write(*,*) alo_subdiag_minus
!write(*,*)
!write(*,*) alo_diag_minus
!write(*,*)
!write(*,*) alo_supdiag_minus
!write(*,*)
!stop
!
!angular integration
mint1d=0.5d0*(int1d_minus+int1d_plus)
alo_diag=0.5d0*(alo_diag_minus+alo_diag_plus)
alo_subdiag=0.5d0*(alo_subdiag_minus+alo_subdiag_plus)
alo_supdiag=0.5d0*(alo_supdiag_minus+alo_supdiag_plus)
!

!do i=1, nz
!   write(*,*) alo_subdiag(i),alo_diag(i),alo_supdiag(i)
!enddo
!stop
!write(*,*) alo_subdiag
!write(*,*)
!write(*,*) alo_diag
!write(*,*)
!write(*,*) alo_supdiag
!stop
!
!
end subroutine sc1st
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine sc2nd(nz, mu, z, bnue1d, scont1d, tau1d, alo_diag, alo_subdiag, alo_supdiag, mint1d, imask1d)
!
!1d sc scheme with quadratic source function (tau needs to be equidistant)
!
use prog_type
!
implicit none
!
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: mu
integer(i4b), dimension(nz), intent(in) :: imask1d
real(dp), dimension(nz), intent(in) :: bnue1d, z, scont1d, tau1d
real(dp), dimension(nz) :: mint1d, alo_diag, alo_subdiag, alo_supdiag
!
! ... local scalars
integer(i4b) :: i
real(dp) :: del_tau, del_tau1, del_tau2, acoeff, bcoeff, ccoeff, dcoeff, e0, e1, e2
!
! ... local arrays
real(dp), dimension(nz) :: int1d_plus, int1d_minus, alo_diag_plus, alo_diag_minus, &
                           alo_subdiag_plus, alo_subdiag_minus, alo_supdiag_plus, alo_supdiag_minus
!real(dp), dimension(nz) :: int1d_plus2, int1d_minus2, alo_diag_plus2, &
!                           alo_subdiag_plus2, alo_subdiag_minus2, alo_supdiag_plus2, alo_supdiag_minus2
!
!
!
!positive direction (intensity)
int1d_plus = 0.d0
alo_diag_plus = 0.d0
alo_subdiag_plus = 0.d0
alo_supdiag_plus = 0.d0

int1d_plus(1) = bnue1d(1)
int1d_plus(2) = bnue1d(2)

do i=3, nz-2
!standard radiative transfer
   del_tau1 = (tau1d(i-1)-tau1d(i))/mu   !along direction 1.d0/sqrt(3) to be consistent with diffusion
   del_tau2 = (tau1d(i)-tau1d(i+1))/mu   !along direction 1.d0/sqrt(3) to be consistent with diffusion
   del_tau=del_tau1+del_tau2
   e0=1.d0-exp(-del_tau1)
   e1=del_tau1-e0
   e2=del_tau1**2-2.d0*e1
   acoeff=e0+(e2-(2.d0*del_tau1+del_tau2)*e1)/(del_tau1*del_tau)
   bcoeff=(del_tau*e1-e2)/(del_tau1*del_tau2)
   ccoeff=(e2-del_tau1*e1)/(del_tau2*del_tau)
   dcoeff=exp(-del_tau1)
   int1d_plus(i) = acoeff*scont1d(i-1) + bcoeff*scont1d(i) + ccoeff*scont1d(i+1) + dcoeff*int1d_plus(i-1)
!
   alo_subdiag_plus(i) = alo_diag_plus(i-1)*dcoeff + acoeff
   alo_diag_plus(i) = alo_supdiag_plus(i-1)*dcoeff + bcoeff
   alo_supdiag_plus(i) = ccoeff
enddo

!
!write(*,*) alo_subdiag_plus
!write(*,*)
!write(*,*) alo_diag_plus
!write(*,*)
!write(*,*) alo_supdiag_plus
!write(*,*)
!stop
!
!negative direction (intensity)
int1d_minus(nz)=0.d0
alo_diag_minus = 0.d0
alo_subdiag_minus = 0.d0
alo_supdiag_minus = 0.d0

do i=nz-2, 3, -1
   del_tau1 = (tau1d(i)-tau1d(i+1))/mu      !along direction 1.d0/sqrt(3) to be consistent with diffusion
   del_tau2 = (tau1d(i-1)-tau1d(i))/mu      !along direction 1.d0/sqrt(3) to be consistent with diffusion
   del_tau=del_tau1+del_tau2
   
   e0=1.d0-exp(-del_tau1)
   e1=del_tau1-e0
   e2=del_tau1**2-2.d0*e1
   acoeff=e0+(e2-(2.d0*del_tau1+del_tau2)*e1)/(del_tau1*del_tau)
   bcoeff=(del_tau*e1-e2)/(del_tau1*del_tau2)
   ccoeff=(e2-del_tau1*e1)/(del_tau2*del_tau)
   dcoeff=exp(-del_tau1)
!
   int1d_minus(i) = acoeff*scont1d(i+1) + bcoeff*scont1d(i) + ccoeff*scont1d(i-1) + dcoeff*int1d_minus(i+1)
!
   alo_supdiag_minus(i) = acoeff + dcoeff*alo_diag_minus(i+1)
   alo_subdiag_minus(i) = ccoeff
   alo_diag_minus(i) = bcoeff + dcoeff*alo_subdiag_minus(i+1)
!   write(*,*) del_tau1, del_tau2, del_tau, acoeff, bcoeff, ccoeff, dcoeff
enddo
!stop
!
!write(*,*) alo_subdiag_minus
!write(*,*)
!write(*,*) alo_diag_minus
!write(*,*)
!write(*,*) alo_supdiag_minus
!write(*,*)
!stop
!
!angular integration
mint1d=0.5d0*(int1d_minus+int1d_plus)
alo_diag=0.5d0*(alo_diag_minus+alo_diag_plus)
alo_subdiag=0.5d0*(alo_subdiag_minus+alo_subdiag_plus)
alo_supdiag=0.5d0*(alo_supdiag_minus+alo_supdiag_plus)
!
end subroutine sc2nd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine searchlight(nz, mu, tau1d, bnue1d)
!
!
use prog_type
!
implicit none
!
integer(i4b), intent(in) :: nz
real(dp), intent(in) :: mu
real(dp), dimension(nz), intent(in) :: tau1d, bnue1d
!
! ... local scalars
integer(i4b) :: i
real(dp) :: del_tau
!
! ... local arrays
real(dp), dimension(nz) :: int1d_1st, int1d_2nd, int1d_abs, scont1d, dtau1d
!
do i=1, nz
   scont1d(i)=tau1d(i)/maxval(tau1d)
!   scont1d(i)=0.d0
enddo
!
dtau1d(1) = tau1d(1)-tau1d(2)
do i=2, nz-1
   dtau1d(i) = (tau1d(i-1)-tau1d(i+1))/2.d0
enddo
dtau1d(nz) = tau1d(nz-1)-tau1d(nz)
!
!-------------------------positive direction----------------------------
!
int1d_1st(1) = bnue1d(1)
int1d_2nd(1) = bnue1d(1)
int1d_abs(1) = bnue1d(1)
do i=2, nz
!first order scheme
   del_tau = dtau1d(i)/mu !along direction 1.d0/sqrt(3) to be consistent with diffusion
   int1d_1st(i) = del_tau/(1.d0+del_tau)*scont1d(i) + int1d_1st(i-1)/(1.d0+del_tau)
!second order scheme
   del_tau = (tau1d(i-1)-tau1d(i))/2.d0/mu !along direction 1.d0/sqrt(3) to be consistent with diffusion
   int1d_2nd(i) = del_tau/(1.d0+del_tau)*scont1d(i) + del_tau/(1.d0+del_tau)*scont1d(i-1) + (1.d0-del_tau)/(1.d0+del_tau) * int1d_2nd(i-1)
!only absorption part
   int1d_abs(i) = exp(-2.d0*del_tau)*int1d_abs(i-1)
enddo
!
!-------------------------negative direction----------------------------
!
int1d_1st=0.d0
int1d_2nd=0.d0
int1d_1st(nz)=bnue1d(1)
int1d_2nd(nz)=bnue1d(1)
int1d_abs(nz)=bnue1d(1)
do i=nz-1, 1, -1
   del_tau = dtau1d(i)/mu !along direction 1.d0/sqrt(3) to be consistent with diffusion
   int1d_1st(i) = del_tau/(1.d0+del_tau)*scont1d(i) + int1d_1st(i+1)/(1.d0+del_tau)
   del_tau = (tau1d(i)-tau1d(i+1))/2.d0/mu
   int1d_2nd(i) = del_tau/(1.d0+del_tau)*scont1d(i) + del_tau/(1.d0+del_tau)*scont1d(i+1) + (1.d0-del_tau)/(1.d0+del_tau) * int1d_2nd(i+1)
   int1d_abs(i) = exp(-2.d0*del_tau)*int1d_abs(i+1)
enddo
!
!-----------------------------------------------------------------------
!
!
open(1, file='searchlight.dat', form='formatted')
   do i=1, nz
      write(1,'(4es20.8)') tau1d(i), int1d_1st(i), int1d_2nd(i), int1d_abs(i)
   enddo
close(1)
!
!
!
end subroutine searchlight
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine tridag(a,b,c,r,u,n)
!-------------------------------------------------------------------------
!                TRIDIAGONAL SOLVER OF "NUMERICAL RECIPES"
!-----------------------------------------------------------------------
  implicit none
  integer :: n
  double precision :: a(n),b(n),c(n),r(n),u(n)
  integer :: j
  double precision :: bet,gam(n)
  if(b(1).eq.0.d0) stop 'tridag: rewrite equations'
  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     if(bet.eq.0.) stop 'tridag failed'
     u(j)=(r(j)-a(j)*u(j-1))/bet
  enddo
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  enddo
  return
end subroutine tridag
