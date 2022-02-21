
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_input
!
!-----------------------------------------------------------------------
!-------------------read in all input parameters------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const, only: rsu
use options_surfb, only: indat_file, input_file1, input_file2, output_file, input_mod, &
                         opt_opac, opt_scont, opt_sline, nxobs_surfb
use params_surfb, only: yhe, hei, vmicro, trad, eps_cont, kcont, eps_line, kline, iline, &
                        vth_fiducial, xnue0, na, gl, gu, flu, xic1, &
                        unit_length, unit_density, unit_velocity, unit_temperature
use mod_opal
use mod_math, only: bnue
!
implicit none
!
! ... local scalars
!
! ... local characters
!
! ... local functions
!
! ... namelist
namelist / input_options / input_file1, input_file2, output_file, input_mod, nxobs_surfb
namelist / input_model / yhe, hei, vmicro, trad
namelist / input_units / unit_length, unit_density, unit_velocity, unit_temperature, vth_fiducial
namelist / input_line / opt_sline, iline, eps_line, kline
namelist / input_cont / opt_opac, opt_scont, eps_cont, kcont
!
!-----------------------------------------------------------------------
!
write(*,*) '----------------------------read input-----------------------------------------'
write(*,*) 'input file name (*.nml) to define model-atmosphere'
read(*,*) indat_file
write(*,*) 'reading input from file: ', trim(indat_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
open(1, file=trim(indat_file), status='old', form='formatted')

!read options, i.e. which model shall be used
   rewind 1
   read(1, nml=input_options)

!read model parameters of star
   rewind 1
   read(1, nml=input_model)
!
!read line strength parameters etc
   rewind 1
   read(1, nml=input_line)
   call get_iline(iline, xnue0, na, gl, gu, flu)
!
!read continuum options
   rewind 1
   read(1, nml=input_cont)
!
!read units
   rewind 1
   read(1, nml=input_units)
!
close(1)
!
!
!velocities in cm/s
vth_fiducial = vth_fiducial * 1.d5
vmicro = vmicro * 1.d5

!unit length in cm
unit_length=unit_length*rsu
!
!core intensity
xic1 = bnue(xnue0,trad)
!
!opal tables
if(opt_opac.eq.1) call get_opal_table(yhe)      
!
end subroutine read_input
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine read_model3d
!
use prog_type
use fund_const
use options_surfb, only: input_file2, opt_opac, opt_sline, opt_scont
use dime_model3d, only: nx, ny, nz, x, y, z, tgas3d, trad3d, velx3d, vely3d, velz3d, &
                        rho3d, opac3d, opalbar3d, scont3d, sline3d
use params_surfb, only: trad, tmin, vth_min, vth_lim, xnue0, vth_fiducial, vmicro, na, vmax, yhe, hei, vmin, &
                        iline, kcont, eps_cont, eps_line, kline, xic1, &
                        unit_density, unit_length, unit_temperature, unit_velocity
use mod_surfb, only: del_vel, del_vel2
use mod_opal
use mod_math, only: bnue
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: i, j, k, err
real(dp) :: b2, b3, dilfac, rad
!
! ... local arrays
integer(i4b), dimension(3) :: indx_tmin
!
! ... local logicals
!
! ... for hdf5
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1), parameter :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z
integer(hsize_t), dimension(3) :: dims_3d
!
! ... local function
real(dp) :: get_opalbar, opac_thomson, opac_opal, sline_depcoeff, vthermal
!
write(*,*) '----------reading 3d model atmosphere from-------------'
write(*,*) 'file-name: ', input_file2
write(*,*)
!
!---------------------------read h5-file--------------------------------
!
!open hdf5-interface
call h5open_f(err)
   call h5fopen_f(trim(input_file2), h5f_acc_rdonly_f, file_id, err)
!
!--------------------------read dimensions------------------------------
!
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
!
!
!set dimension-arrays
      dims_x=(/ nx /)
      dims_y=(/ ny /)
      dims_z=(/ nz /)
      dims_3d=(/ nx, ny, nz /)
!
!------------------------read units and overwrite from------------------
!--------------------------subroutine read_input------------------------

      call h5gopen_f(file_id, 'units', group_id, err)
         call h5aopen_f(group_id, 'unit_length', dset_id, err)
            call h5aread_f(dset_id, h5t_native_double, unit_length, dims_scalars, err)
         call h5aclose_f(dset_id, err)
         call h5aopen_f(group_id, 'unit_velocity', dset_id, err)
            call h5aread_f(dset_id, h5t_native_double, unit_velocity, dims_scalars, err)
         call h5aclose_f(dset_id, err)
         call h5aopen_f(group_id, 'unit_density', dset_id, err)
            call h5aread_f(dset_id, h5t_native_double, unit_density, dims_scalars, err)
         call h5aclose_f(dset_id, err)
         call h5aopen_f(group_id, 'unit_temperature', dset_id, err)
            call h5aread_f(dset_id, h5t_native_double, unit_temperature, dims_scalars, err)
         call h5aclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!unit length in cm
      unit_length=unit_length*rsu      
!
!-----------------------coordinates-------------------------------------
!
!allocate arrays
      allocate(x(nx), stat=err)
      allocate(y(ny), stat=err)
      allocate(z(nz), stat=err)

      allocate(tgas3d(nx, ny, nz), stat=err)
      allocate(trad3d(nx, ny, nz), stat=err)            
      allocate(rho3d(nx, ny, nz), stat=err)
      allocate(velx3d(nx, ny, nz), stat=err)
      allocate(vely3d(nx, ny, nz), stat=err)
      allocate(velz3d(nx, ny, nz), stat=err)                  

      call h5gopen_f(file_id, 'coordinates', group_id, err)
         call h5dopen_f(group_id, 'x', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, x, dims_x, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'y', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, y, dims_y, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'z', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, z, dims_z, err)
         call h5dclose_f(dset_id, err)
      call h5gclose_f(group_id, err)
!
!--------------------------------3d model-------------------------------
!
      call h5gopen_f(file_id, 'model', group_id, err)
         call h5dopen_f(group_id, 'rho', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, rho3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velx', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'vely', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'velz', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'tgas', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, tgas3d, dims_3d, err)
         call h5dclose_f(dset_id, err)
         call h5dopen_f(group_id, 'trad', dset_id, err)
            call h5dread_f(dset_id, h5t_native_double, trad3d, dims_3d, err)
         call h5dclose_f(dset_id, err)            
      call h5gclose_f(group_id, err)
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!-----------------------------------------------------------------------
!
!transform density and temperature to cgs
rho3d = rho3d*unit_density
tgas3d = tgas3d*unit_temperature
trad3d = trad3d*unit_temperature
!trad3d = trad/2.d0
!
!measure velocities in vth_fiducial
velx3d = velx3d*unit_velocity/vth_fiducial
vely3d = vely3d*unit_velocity/vth_fiducial
velz3d = velz3d*unit_velocity/vth_fiducial
!
!minimum temperature, minimum thermal velocity, maximum allowed velocity steps in vth_fiducial
tmin=minval(tgas3d)
vth_min = max(vth_lim,vthermal(vmicro, tmin, na))   !thermal velocity not allowed to be smaller than 5 km/s for computational reasons
del_vel2=del_vel*vth_min/vth_fiducial
!
!maximum terminal velocity (to extend frequency range)
vmax = maxval(velz3d)*vth_fiducial
vmin = minval(velz3d)*vth_fiducial
!
!----------------model for source functions and opacities---------------
!
if(opt_scont.eq.4.or.opt_sline.eq.4) then
!read in source function from output-file of sc3d.eo
   stop 'error in read_model3d: opt_scont=4 and opt_sline=4 not implemented'
endif
!
!
allocate(scont3d(nx,ny,nz), stat=err)
allocate(sline3d(nx,ny,nz), stat=err)
allocate(opac3d(nx,ny,nz), stat=err)
allocate(opalbar3d(nx,ny,nz), stat=err)
!
do k=1, nz
   do j=1, ny
      do i=1, nx
!         
!         
!---------------------calculate continuum opacity-----------------------
!         
         select case(opt_opac)
            case(0)
               opac3d(i,j,k) = zero
            case(1)
!               opac3d(i,j,k) = opac_thomson(yhe,hei,rho3d(i,j,k),kcont) !in cgs
               opac3d(i,j,k) = opac_thomson(yhe,hei,rho3d(i,j,k),kcont)*unit_length !in 1/unit_length
            case(2)
               !in 1/unit_length               
               opac3d(i,j,k) = opac_opal(kcont, yhe, hei, log10(rho3d(i,j,k)),log10(tgas3d(i,j,k)), nrho_opal, ntemp_opal, rho_opal, temp_opal, kappa_opal)*unit_length
            case default
               stop 'error in read_model3d: opt_opac not specified'
         end select
!         
!---------------------calculate frequency integrated line opacity-------
!
         b2 = one
         b3 = one
!         opalbar3d(i,j,k) = get_opalbar(iline, kline, unit_length, yhe, hei, tgas3d(i,j,k), vth_fiducial, xnue0, b2, b3, rho3d(i,j,k))/unit_length !in cgs
         opalbar3d(i,j,k) = get_opalbar(iline, kline, unit_length, yhe, hei, tgas3d(i,j,k), vth_fiducial, xnue0, b2, b3, rho3d(i,j,k)) !in 1/unit_length
!         
!---------------------calculate continuum source function---------------
!         
         select case(opt_scont)
            case(0)
               scont3d(i,j,k) = zero
            case(1)
               scont3d(i,j,k) = bnue(xnue0,trad3d(i,j,k))
            case(2)
               rad = sqrt((x(i)**2+y(j)**2+z(k)**2))
               if(rad.lt.one) rad=one
               dilfac = 0.5d0*(one-sqrt(one-one/rad**2))
               scont3d(i,j,k) = xic1*dilfac
            case default
              stop 'error in read_model3d: opt_scont not specified'
         end select
!         
!---------------------calculate line source function--------------------
!         
         select case(opt_sline)
            case(0)
               sline3d(i,j,k) = zero
            case(1)
               sline3d(i,j,k) = bnue(xnue0,trad3d(i,j,k))
            case(2)
               rad = sqrt((x(i)**2+y(j)**2+z(k)**2))
               if(rad.lt.one) rad=one
               dilfac = half*(one-sqrt(one-one/rad**2))               
               sline3d(i,j,k) = xic1*dilfac
            case(3)
               b2 = one
               b3 = one
               sline3d(i,j,k) = sline_depcoeff(xnue0, tgas3d(i,j,k), b2, b3)
            case default
               stop 'error in read_model3d: opt_sline not specified'
         end select
!
!-----------------------------------------------------------------------
!         
      enddo
   enddo
enddo
!
!
!
end subroutine read_model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine print_model3d
!
use prog_type
use fund_const, only: rsu  
use dime_model3d, only: nx, ny, nz, x, y, z, &
                        opac3d, opalbar3d, velx3d, vely3d, velz3d, sline3d, scont3d, tgas3d, trad3d, rho3d
use params_surfb, only: xic1, vmicro, xnue0, vth_fiducial, na, tmin, vmin, vmax, unit_length, vth_min
use mod_surfb, only: del_vel, del_vel2
!
implicit none
!
! ... arguments
integer(i4b) :: i, j, k
!
write(*,*) '---------------------3d atomspheric structure along z--------------------------'
write(*,*)
!
write(*,'(a20, es20.8)') 'I_c', xic1
write(*,'(a20, es20.8)') 'xnue0', xnue0
write(*,'(a20, es20.8)') 'tmin [K]', tmin
write(*,'(a20, i20)') 'na', na
write(*,'(a20, es20.8)') 'vmin [km/s]', vmin/1.d5
write(*,'(a20, es20.8)') 'vmax [km/s]', vmax/1.d5
write(*,'(a20, es20.8)') 'v_micro [km/s]', vmicro/1.d5
write(*,'(a20, es20.8)') 'vth_fiducial [km/s]', vth_fiducial/1.d5
write(*,'(a20, es20.8)') 'vth_min [km/s]', vth_min/1.d5
write(*,'(a20, 2es20.8)') 'delta-v [vth],[vth_fiducial]', del_vel, del_vel2
write(*,'(a20, es20.8)') 'unit_length [r_sun]', unit_length/rsu
write(*,*)
!
i=nx/2+1
j=ny/2+1
!
write(*,*) 'at x, y', x(i), y(j)
write(*,'(a4, 11(a14))') '#', 'z [l_0]', 'opac [1/l_0]', 'opalbar [1/l_0]', 'velx [vth*]', &
                        'vely [vth*]', 'velz [vth*]', 't_gas [K]', 't_rad [K]', 'rho[g/cm^3]', 's_cont', 's_line'
do k=1, nz
   write(*,'(i4, 11(es14.6))')  k, z(k), opac3d(i,j,k), opalbar3d(i,j,k), velx3d(i,j,k), vely3d(i,j,k), &
                              velz3d(i,j,k), tgas3d(i,j,k), trad3d(i,j,k), rho3d(i,j,k), scont3d(i,j,k), sline3d(i,j,k)
enddo
write(*,*)
!
!stop 'go on in print_model'
!
end subroutine print_model3d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_surface(ixobs)
!
!-------------output of surface emergent intensity for given------------
!--------------xobs and direction specified by alpha, gamma-------------
!
use prog_type
use options_surfb, only: output_file, nxobs_surfb
use dime_surfb, only: nxp, nyp
use mod_surfb, only: xp, yp, xobs_surface, iem_surface, iemi_surface, iabs_surface, icont_surface
use params_surfb, only: xic1, vth_fiducial, xnue0
use hdf5
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: ixobs
!
! ... local characters
character(len=9) :: gname
!
! ... local scalars
integer(i4b) :: err, i, j
integer(hid_t) :: file_id, dspace_id, aspace_id, dset_id, group_id, group_id2, attr_id
!
! ... local arrays
integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_x , dims_y
integer(hsize_t), dimension(2) :: dims_2d


!i=nxp/2+1
!j=nyp/2+1
!!do i=1, nxp
!do j=1, nyp
!   write(*,*) iem_surface(i,j)/xic1
!enddo
!stop 'go on in iem_surface'

!
write(*,*) '------------------------------output to file-----------------------------------'
write(*,*) trim(output_file)
!
!-----------------------------------------------------------------------
!
dims_x = (/ nxp /)
dims_y = (/ nyp /)
dims_2d = (/ nxp, nyp /)
!
!-----------------------------------------------------------------------
!
!initialize fortran interface
call h5open_f(err)
!
call h5fopen_f(trim(output_file), h5f_acc_rdwr_f, file_id, err)
!-----------------------------------------------------------------------
!
if(ixobs.eq.1) then
   call h5gcreate_f(file_id, 'surfb', group_id, err)
!
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
         call h5acreate_f(group_id, 'nx', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, nxp, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5acreate_f(group_id, 'ny', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, nyp, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5acreate_f(group_id, 'nxobs', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, nxobs_surfb, dims_scalars, err)
         call h5aclose_f(attr_id, err)            
      call h5sclose_f(aspace_id, err)

      call h5screate_simple_f(1, dims_x, dspace_id, err)
         call h5dcreate_f(group_id, 'xp', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, xp, dims_x, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id, 'yp', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, yp, dims_y, err)
         call h5dclose_f(dset_id, err)
      call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)
endif
!
call h5gopen_f(file_id, 'surfb', group_id, err)
!
   write(gname,'(a4, i5.5)') 'xobs', ixobs
   call h5gcreate_f(group_id, gname, group_id2, err)
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
         call h5acreate_f(group_id2, 'xobs', h5t_native_double, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xobs_surface(ixobs), dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5sclose_f(aspace_id, err)
      call h5screate_simple_f(2, dims_2d, dspace_id, err)      
         call h5dcreate_f(group_id2, 'iem_surface', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iem_surface, dims_2d, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id2, 'iemi_surface', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iemi_surface, dims_2d, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id2, 'iabs_surface', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, iabs_surface, dims_2d, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id2, 'icont_surface', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, icont_surface, dims_2d, err)
         call h5dclose_f(dset_id, err)            
      call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id2, err)

call h5gclose_f(group_id, err)
!
!
call h5fclose_f(file_id, err)
!
!
!
end subroutine output_surface
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_fluxem(indx)
!
!-----------------output of emergent flux profile-----------------------
!
use prog_type
use fund_const, only: zero
use options_surfb, only: output_file
use fund_const, only: pi
use dime_surfb, only: nxobs_fs
use mod_surfb, only: flux_tot, flux_cont, flux_emi, flux_abs, xobs
use hdf5
!
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: indx
!
! ... local scalars
integer(i4b) :: i, err

! ... local characters
character(len=11) :: gname
!
!
! ... for hdf5
integer(hid_t) :: file_id, dspace_id, aspace_id, dset_id, group_id, group_id2, attr_id
integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_xobs
!
!reset small values
do i=1, nxobs_fs
   if(flux_tot(i).lt.1.d-15) flux_tot(i)=zero
   if(flux_cont(i).lt.1.d-15) flux_cont(i)=zero
   if(flux_emi(i).lt.1.d-80) flux_emi(i)=zero
   if(flux_abs(i).lt.1.d-80) flux_abs(i)=zero
enddo
!
write(*,*) '------------------------------output to file-----------------------------------'
write(*,*) trim(output_file)
!
dims_xobs = (/ nxobs_fs /)
!
!initialize fortran interface
call h5open_f (err)
!
call h5fopen_f(trim(output_file), h5f_acc_rdwr_f, file_id, err)
!
!------------------------------parameters-----------------------------------
if(indx.eq.1) then
   call h5gcreate_f(file_id, 'fluxem', group_id, err)
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
         call h5acreate_f(group_id, 'nxobs_fs', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, nxobs_fs, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5sclose_f(aspace_id, err)
   call h5gclose_f(group_id, err)
endif
!
call h5gopen_f(file_id, 'fluxem', group_id, err)
!
   write(gname,'(a6, i5.5)') 'fluxem', indx
   call h5gcreate_f(group_id, gname, group_id2, err)
      call h5screate_simple_f(1, dims_xobs, dspace_id, err)      
         call h5dcreate_f(group_id2, 'xobs', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, xobs, dims_xobs, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id2, 'ftot', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, flux_tot, dims_xobs, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id2, 'fcont', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, flux_cont, dims_xobs, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id2, 'femi', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, flux_emi, dims_xobs, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id2, 'fabs', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, flux_abs, dims_xobs, err)
         call h5dclose_f(dset_id, err)            
         call h5dcreate_f(group_id2, 'fnorm', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, flux_tot/flux_cont, dims_xobs, err)
         call h5dclose_f(dset_id, err)
      call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id2, err)
!
call h5gclose_f(group_id, err)
!
!
call h5fclose_f(file_id, err)



!
end subroutine output_fluxem
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_photprof
!
!-----------------output of emergent flux profile-----------------------
!
use prog_type
use fund_const, only: zero
use options_surfb, only: output_file
use fund_const, only: pi
use dime_surfb, only: nxobs_fs
use mod_surfb, only: xic_nue, xicc_nue, xobs
use params_surfb, only: vth_fiducial, xic1, xnue0
use hdf5
!
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, err

! ... local characters
character(len=10) :: gname
!
!
! ... for hdf5
integer(hid_t) :: file_id, dspace_id, aspace_id, dset_id, group_id, group_id2, attr_id
integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_xobs
!

!
write(*,*) '----------------------output photospheric profile to file----------------------'
write(*,*) trim(output_file)
write(*,*)
!
dims_xobs = (/ nxobs_fs /)
!
!initialize fortran interface
call h5open_f (err)
!
!create hdf5 file
call h5fcreate_f(trim(output_file), h5f_acc_trunc_f, file_id, err)
!
!------------------------------parameters-----------------------------------
!
   call h5gcreate_f(file_id, 'input_parameters', group_id, err)
!
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
         call h5acreate_f(group_id, 'xic1', h5t_native_double, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xic1, dims_scalars, err)
         call h5aclose_f(attr_id, err)
         call h5acreate_f(group_id, 'vth_fiducial', h5t_native_double, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, vth_fiducial, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5acreate_f(group_id, 'xnue0', h5t_native_double, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_double, xnue0, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5sclose_f(aspace_id, err)   
!
   call h5gclose_f(group_id, err)
!
!----------------------------photospheric profile---------------------------
!   
   call h5gcreate_f(file_id, 'photospheric_profile', group_id, err)
      call h5screate_simple_f(1, dims_scalars, aspace_id, err)
         call h5acreate_f(group_id, 'nxobs_fs', h5t_native_integer, aspace_id, attr_id, err)
            call h5awrite_f(attr_id, h5t_native_integer, nxobs_fs, dims_scalars, err)
         call h5aclose_f(attr_id, err)
      call h5sclose_f(aspace_id, err)
      call h5screate_simple_f(1, dims_xobs, dspace_id, err)      
         call h5dcreate_f(group_id, 'xobs', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, xobs, dims_xobs, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id, 'xic_nue', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, xic_nue, dims_xobs, err)
         call h5dclose_f(dset_id, err)
         call h5dcreate_f(group_id, 'xicc_nue', h5t_native_double, dspace_id, dset_id, err)
            call h5dwrite_f(dset_id, h5t_native_double, xicc_nue, dims_xobs, err)
         call h5dclose_f(dset_id, err)
      call h5sclose_f(dspace_id, err)
   call h5gclose_f(group_id, err)!
!
call h5fclose_f(file_id, err)
!
!
!
end subroutine output_photprof
