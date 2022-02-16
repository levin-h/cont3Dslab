!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_pp
!
!-------calculate model atmosphere from beta velocity law--------------
!
use prog_type
use fund_const
use dime_modext, only: nx, ny, nz
use model3d, only: x_modext3d, y_modext3d, z_modext3d, rho_modext3d, velx_modext3d, &
                   vely_modext3d, velz_modext3d, tgas_modext3d, trad_modext3d
use params_input, only: xmin, xmax, ymin, ymax, zmin, zmax, &
                        vmin, vmin_cgs, vmax, vmax_cgs, xmloss, xmloss_cgs, &
                        rstar, rstar_cgs, tmin, teff, beta, unit_length, unit_density, unit_velocity
use mod_directories, only: indat_file
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
real(dp) :: bconst, zshift
real(dp) :: velr, rad
!
! ... local characters
!
! ... local functions
real(dp) :: bvel
!
! ... namelists
namelist / input_usr / rstar, teff, tmin, xmloss, vmin, vmax, beta, xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz
!
write(*,*) '-----------------------creating 3d model atmosphere----------------------------'
write(*,*)
!
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
!   
!read 2d model parameters
   rewind 1
   read(1, nml=input_usr)
   rstar_cgs=rstar*rsu
   vmin_cgs=vmin*1.d5
   vmax_cgs=vmax*1.d5
   xmloss_cgs=xmloss*xmsu/yr
   tmin=teff*tmin
   
close(1)
!
!-----------------------allocate arrays---------------------------------
!
allocate(rho_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
allocate(velx_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
allocate(vely_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
allocate(velz_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
allocate(tgas_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
allocate(trad_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
allocate(x_modext3d(nx), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
allocate(y_modext3d(ny), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
allocate(z_modext3d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_pp'
!
!----------------------equidistant grid in x----------------------------
!
if(nx.le.3) stop 'error in calc_model3d_pp: nx needs to be larger than 3'   
do i=1, nx
   x_modext3d(i) = xmin + (i-1)*(xmax-xmin)/(nx-1)
enddo
!
!----------------------equidistant grid in y----------------------------
!
if(ny.le.3) stop 'error in calc_model3d_pp: ny needs to be larger than 3'   
do i=1, ny
   y_modext3d(i) = ymin + (i-1)*(ymax-ymin)/(ny-1)
enddo
!
!----------------------equidistant grid in z----------------------------
!
if(nz.le.5) stop 'error in calc_model3d_pp: nz needs to be larger than 5'   
do i=1, nz
   z_modext3d(i) = zmin + (i-1)*(zmax-zmin)/(nz-1)
enddo
!
!------------------calculate/define required constants------------------
!
!b-factor for beta-velocity-law
bconst=one-(vmin_cgs/vmax_cgs)**(one/beta)
!
if(z_modext3d(1).lt.0.) stop 'radial beta-velocity law not meaniningful for such z-range'
!
zshift=zero
if(z_modext3d(1).lt.1.) then
   write(*,*) 'shifting z-values to get a radial beta-velocity law starting at r=1.'
   write(*,*)
   zshift=one-z_modext3d(1)
endif
!
do i=1, nx
   do j=1, ny
      do k=1, nz
         rad = z_modext3d(k)+zshift
         velr = bvel(rad, vmax_cgs, bconst, beta)
         velz_modext3d(i,j,k) = velr
         velx_modext3d(i,j,k) = zero
         vely_modext3d(i,j,k) = zero
         rho_modext3d(i,j,k) = xmloss_cgs/four/pi/(rstar_cgs*rad)**2/velr
         trad_modext3d(i,j,k) = tmin
         tgas_modext3d(i,j,k) = tmin         
      enddo
   enddo
enddo
!
!transform everything to own units
x_modext3d = x_modext3d*rstar_cgs/unit_length/rsu
y_modext3d = y_modext3d*rstar_cgs/unit_length/rsu
z_modext3d = z_modext3d*rstar_cgs/unit_length/rsu
rho_modext3d = rho_modext3d/unit_density
velx_modext3d = velx_modext3d/unit_velocity
vely_modext3d = vely_modext3d/unit_velocity
velz_modext3d = velz_modext3d/unit_velocity
!
end subroutine calc_model3d_pp
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_nicowr2d
!
!-------calculate model atmosphere from beta velocity law--------------
!
use prog_type
use fund_const
use dime_modext, only: nx, ny, nz
use model3d, only: x_modext3d, y_modext3d, z_modext3d, rho_modext3d, velx_modext3d, &
                   vely_modext3d, velz_modext3d, tgas_modext3d, trad_modext3d
use params_input, only: xmin, xmax, ymin, ymax, zmin, zmax, &
                        vmin, vmin_cgs, vmax, vmax_cgs, xmloss, xmloss_cgs, &
                        rstar, rstar_cgs, tmin, teff, beta
use mod_directories, only: indat_file
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
real(dp) :: bconst, zshift
real(dp) :: velr, rad
!
! ... local characters
character(len=500) :: fname_model
!
! ... namelists
namelist / input_usr / fname_model
!
! ... for hdf5 file
integer(hid_t) :: file_id, group_id, dset_id
integer(hsize_t), dimension(1), parameter :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_x, dims_z
integer(hsize_t), dimension(2) :: dims_2d
integer(i4b) :: nx_modext2d, nz_modext2d
real(dp), dimension(:), allocatable :: x_modext2d, z_modext2d
real(dp), dimension(:,:), allocatable :: rho_modext2d, velx_modext2d, vely_modext2d, velz_modext2d, tgas_modext2d, trad_modext2d
!
! ... local functions
!
!
write(*,*) '---------------------read model from Nicos 2D WR models------------------------'
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
write(*,*) 'file name: ', trim(fname_model)
write(*,*)
!
call h5open_f(err)
call h5fopen_f(trim(fname_model), h5f_acc_rdonly_f, file_id, err)
!
!--------------------------read dimensions------------------------------
!
call h5gopen_f(file_id, 'dimensions', group_id, err)
!
   call h5aopen_f(group_id, 'nx', dset_id, err)
      call h5aread_f(dset_id, h5t_native_integer, nx_modext2d, dims_scalars, err)
   call h5aclose_f(dset_id, err)
   call h5aopen_f(group_id, 'nz', dset_id, err)
      call h5aread_f(dset_id, h5t_native_integer, nz_modext2d, dims_scalars, err)
   call h5aclose_f(dset_id, err)
!
call h5gclose_f(group_id, err)
!
!--------------------allocate coordinate and 2d arrays------------------
!
allocate(rho_modext2d(nx_modext2d,nz_modext2d), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(velx_modext2d(nx_modext2d,nz_modext2d), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(vely_modext2d(nx_modext2d,nz_modext2d), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(velz_modext2d(nx_modext2d,nz_modext2d), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(tgas_modext2d(nx_modext2d,nz_modext2d), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(trad_modext2d(nx_modext2d,nz_modext2d), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(x_modext2d(nx_modext2d), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(z_modext2d(nz_modext2d), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
!
!
!set dimension-arrays
dims_x = (/ nx_modext2d /)
dims_z = (/ nz_modext2d /)
dims_2d = (/ nx_modext2d, nz_modext2d /)
!
!-----------------------coordinates-------------------------------------
!
call h5gopen_f(file_id, 'coordinates', group_id, err)
!
   call h5dopen_f(group_id, 'x', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, x_modext2d, dims_x, err)
   call h5dclose_f(dset_id, err)
!
   call h5dopen_f(group_id, 'z', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, z_modext2d, dims_z, err)
   call h5dclose_f(dset_id, err)
!
call h5gclose_f(group_id, err)
!
!--------------------------model----------------------------------------
!
call h5gopen_f(file_id, 'model', group_id, err)
!
   call h5dopen_f(group_id, 'rho', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, rho_modext2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'tgas', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, tgas_modext2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'trad', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, trad_modext2d, dims_2d, err)
   call h5dclose_f(dset_id, err)      
   call h5dopen_f(group_id, 'velx', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velx_modext2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'vely', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, vely_modext2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
   call h5dopen_f(group_id, 'velz', dset_id, err)
      call h5dread_f(dset_id, h5t_native_double, velz_modext2d, dims_2d, err)
   call h5dclose_f(dset_id, err)
!
call h5gclose_f(group_id, err)
!
call h5fclose_f(file_id, err)
call h5close_f(err)
!
!------------allocate 3d arrays (including ghost zones)-----------------
!
nx = nx_modext2d
ny = nx_modext2d
nz = nz_modext2d
allocate(x_modext3d(nx), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(y_modext3d(ny), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(z_modext3d(nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(rho_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(velx_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(vely_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(velz_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(tgas_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
allocate(trad_modext3d(nx,ny,nz), stat=err)
   if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
!
!-----------------------create the grid---------------------------------
!
x_modext3d = x_modext2d
y_modext3d = x_modext2d
z_modext3d = z_modext2d
!
!------------------calculate/define required constants------------------
!
do i=1, nx
   do j=1, ny
      do k=1, nz
         rho_modext3d(i,j,k) = rho_modext2d(i,k)
         velx_modext3d(i,j,k) = velx_modext2d(i,k)
         vely_modext3d(i,j,k) = vely_modext2d(i,k)
         velz_modext3d(i,j,k) = velz_modext2d(i,k)
         tgas_modext3d(i,j,k) = tgas_modext2d(i,k)
         trad_modext3d(i,j,k) = trad_modext2d(i,k)         
      enddo
   enddo
enddo
!
!
end subroutine calc_model3d_nicowr2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_model3d_nicowr3d
!
!----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime_modext, only: nx, ny, nz
use model3d, only: x_modext3d, y_modext3d, z_modext3d, rho_modext3d, velx_modext3d, &
                   vely_modext3d, velz_modext3d, tgas_modext3d, trad_modext3d
use params_input, only: xmin, xmax, ymin, ymax, zmin, zmax, &
                        vmin_cgs, vmax, vmax_cgs, xmloss, xmloss_cgs, &
                        rstar, rstar_cgs, tmin, teff, beta, unit_length, unit_velocity, unit_density
use mod_directories, only: indat_file
use mod_amrvac_reader
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
integer(i4b) :: nw
integer(i4b) :: irho, ivelx, ively, ivelz, itgas, itrad
real(dp) :: bconst, zshift
real(dp) :: velr, rad, fdum, rho_avg, velr_avg, mdot_avg, tgas_avg, trad_avg, vmin, vinf, sr
integer(i4b) :: idum
!
! ... local characters
character(len=500)  :: fname_model
!
! ... namelists
logical :: opt_read_h5
integer(i4b) :: opt_bvel, max_refinement
integer(i4b), dimension(:), allocatable :: nd
real(dp), dimension(:), allocatable :: stretching
namelist / input_usr / fname_model, opt_bvel, max_refinement, nd, opt_read_h5
!
! ... for hdf5 file
integer(hid_t) :: file_id, group_id, dset_id
integer(hsize_t), dimension(1), parameter :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z
integer(hsize_t), dimension(3) :: dims_3d
!
! ... local functions
real(dp) :: bvel
!
! ... local logicals
logical :: check1
!
! ... local derived types
type(alldata) :: my_data
type(grid) :: my_grid
!
allocate(nd(3), stat=err)
!
write(*,*) '---------------------read model from Nicos 3D WR models------------------------'
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
open(1, file=trim(indat_file), status='old', form='formatted')
   rewind 1
   read(1, nml=input_usr)
close(1)
!
inquire(file=trim(fname_model), exist=check1)
if(.not.check1) then
   write(*,*) 'error in calc_model3d_nicowr3d: file "', trim(fname_model), '" does not exist'
   stop
endif
!
write(*,*) 'file name: ', trim(fname_model)
write(*,*)
!
!old version: read preprocessed h5 file
!
if(opt_read_h5) then
   call h5open_f(err)
   call h5fopen_f(trim(fname_model), h5f_acc_rdonly_f, file_id, err)
      call h5gopen_f(file_id, 'dimensions', group_id, err)

         call h5aopen_f(group_id, 'nx', dset_id, err)
            call h5aread_f(dset_id, h5t_native_integer, nx, dims_scalars, err)
         call h5aclose_f(dset_id, err)
         call h5aopen_f(group_id, 'ny', dset_id, err)
            call h5aread_f(dset_id, h5t_native_integer, ny, dims_scalars, err)
         call h5aclose_f(dset_id, err)      
         call h5aopen_f(group_id, 'nz', dset_id, err)
            call h5aread_f(dset_id, h5t_native_integer, nz, dims_scalars, err)
         call h5aclose_f(dset_id, err)
            
      call h5gclose_f(group_id, err)

      !allocate coordinate and 3d arrays
!
   allocate(x_modext3d(nx), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(y_modext3d(ny), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(z_modext3d(nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(rho_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(velx_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(vely_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(velz_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(tgas_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(trad_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
!
!
!set dimension-arrays
   dims_x = (/ nx /)
   dims_y = (/ ny /)   
   dims_z = (/ nz /)
   dims_3d = (/ nx, ny, nz /)
!
!coordinates
   call h5gopen_f(file_id, 'coordinates', group_id, err)
!
      call h5dopen_f(group_id, 'x', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, x_modext3d, dims_x, err)
      call h5dclose_f(dset_id, err)!!

      call h5dopen_f(group_id, 'y', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, y_modext3d, dims_y, err)
      call h5dclose_f(dset_id, err)      
!
      call h5dopen_f(group_id, 'z', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, z_modext3d, dims_z, err)
      call h5dclose_f(dset_id, err)
!
   call h5gclose_f(group_id, err)
!
!model
   call h5gopen_f(file_id, 'model', group_id, err)
!
      call h5dopen_f(group_id, 'rho', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, rho_modext3d, dims_3d, err)
      call h5dclose_f(dset_id, err)
      call h5dopen_f(group_id, 'tgas', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, tgas_modext3d, dims_3d, err)
      call h5dclose_f(dset_id, err)
      call h5dopen_f(group_id, 'trad', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, trad_modext3d, dims_3d, err)
      call h5dclose_f(dset_id, err)      
      call h5dopen_f(group_id, 'velx', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, velx_modext3d, dims_3d, err)
      call h5dclose_f(dset_id, err)
      call h5dopen_f(group_id, 'vely', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, vely_modext3d, dims_3d, err)
      call h5dclose_f(dset_id, err)
      call h5dopen_f(group_id, 'velz', dset_id, err)
         call h5dread_f(dset_id, h5t_native_double, velz_modext3d, dims_3d, err)
      call h5dclose_f(dset_id, err)
!
   call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
   call h5close_f(err)

else
!read directly from AMRVAC file------------

   allocate(stretching(3), stat=err)
   stretching=(/ one, one, one /)
!
!get the data and the grid for this snapshot
   my_data = get_data(fname_model, levmax_usr=max_refinement, stretching=stretching, nd=nd)
   my_grid = my_data%data%mesh
!
   nx = my_data%data%data_shape(3)
   ny = my_data%data%data_shape(2)
   nz = my_data%data%data_shape(1)
   nw = my_data%data%data_shape(4)
!
!--------------------allocate coordinate and 3d arrays------------------
!
   allocate(x_modext3d(nx), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(y_modext3d(ny), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(z_modext3d(nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(rho_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(velx_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(vely_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(velz_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(tgas_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
   allocate(trad_modext3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in calc_model3d_nicowr'
!
!find indices where data is stored
   do i=1, nw
!   write(*,*) fdata%w_names(i)
      if(trim(my_data%data%w_names(i)).eq.'rho') irho=i
      if(trim(my_data%data%w_names(i)).eq.'v1') ivelz=i
      if(trim(my_data%data%w_names(i)).eq.'v2') ively=i
      if(trim(my_data%data%w_names(i)).eq.'v3') ivelx=i
      if(trim(my_data%data%w_names(i)).eq.'Trad') itrad=i
      if(trim(my_data%data%w_names(i)).eq.'Tgas') itgas=i
   enddo
!
!
!store the data (all in units given by the input units
   x_modext3d = my_data%data%mesh%zgrid%coord
   y_modext3d = my_data%data%mesh%ygrid%coord
   z_modext3d = my_data%data%mesh%xgrid%coord
!
   do i=1, nx
      do j=1, ny
         do k=1, nz
            rho_modext3d(i,j,k) = my_data%data%data(k,j,i,irho)
            velx_modext3d(i,j,k) = my_data%data%data(k,j,i,ivelx)
            vely_modext3d(i,j,k) = my_data%data%data(k,j,i,ively)
            velz_modext3d(i,j,k) = my_data%data%data(k,j,i,ivelz)
            tgas_modext3d(i,j,k) = my_data%data%data(k,j,i,itgas)
            trad_modext3d(i,j,k) = my_data%data%data(k,j,i,itrad)                  
         enddo
      enddo
   enddo

endif
!
!-----overwrite everything with beta-velocity-law if option set---------
!
if(opt_bvel.eq.1) then
   write(*,*) 'calculating beta-velocity law'
   write(*,*)
!
!calculate lateral and azimuthal average velocities, average temperatures, average densities
   do k=1, nz
      trad_avg = 0.d0
      tgas_avg = 0.d0
      rho_avg = 0.d0
      velr_avg = 0.d0
!
      idum = 0
      do i=1, nx
         do j=1, ny
            trad_avg = trad_avg + trad_modext3d(i,j,k)
            tgas_avg = tgas_avg + tgas_modext3d(i,j,k)
            velr_avg = velr_avg + velz_modext3d(i,j,k)
            rho_avg = rho_avg + rho_modext3d(i,j,k)
            idum = idum+1
         enddo
      enddo
      
      rho_modext3d(:,:,k) = rho_avg/idum
      trad_modext3d(:,:,k) = tgas_avg/idum
      tgas_modext3d(:,:,k) = trad_avg/idum
      velz_modext3d(:,:,k) = velr_avg/idum
      velx_modext3d(:,:,k) = zero
      vely_modext3d(:,:,k) = zero
   enddo
!
!calculate radial average radial velociies and mass-loss rates (only outermost 40 grid points)
   idum=0
   mdot_avg = 0.d0
   velr_avg = 0.d0
   sr = unit_length*rsu
   do k=nz-40, nz
      velr_avg = velr_avg + velz_modext3d(1,1,k)
      mdot_avg = mdot_avg + four*pi*(z_modext3d(i)*sr)**2 * rho_modext3d(1,1,k) * velz_modext3d(1,1,k)
      idum=idum+1
   enddo
   velr_avg = velr_avg/idum
   mdot_avg = mdot_avg/idum

   beta=1.d0
   vmin=10.d5
   velr_avg = velr_avg*unit_velocity
   mdot_avg = mdot_avg*unit_velocity*unit_density

!calculate v_inf such that vel_r(beta-law) = vel_r(model) at r_max
   fdum = (vmin/velr_avg)**(one/beta)
   bconst = (one - fdum)/(one - fdum/z_modext3d(nz))
   vinf = vmin/(one-bconst)**beta
!
!or simply input vinf
!   vinf=4000.d5   
!   bconst = one - (vmin/vinf)**(one/beta)   

   write(*,*) 'using'
   write(*,*) 'beta', beta
   write(*,*) 'vmin [km/s]', vmin/1.d5
   write(*,*) 'v_r(r_max) [km/s]', velr_avg/1.d5   
   write(*,*) 'vinf [km/s]', vinf/1.d5
   write(*,*) 'mdot [Msun/yr]', mdot_avg*yr/xmsu
   
   do i=k, nz
      velz_modext3d(:,:,k) = bvel(z_modext3d(k), vinf, bconst, beta)/unit_velocity
      rho_modext3d(:,:,k) = mdot_avg/four/pi/velz_modext3d(1,1,k)/(z_modext3d(k)*sr)**2 / unit_density
   enddo

endif
!
!
!
end subroutine calc_model3d_nicowr3d
