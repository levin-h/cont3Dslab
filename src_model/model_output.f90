!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_mod1d
!
!-----------------save 1d model atmosphere as h5-file--------------------
!
use prog_type
use mod_directories, only: model_dir, model1d_file
use dime_modext, only: nz
use model1d, only: z_modext1d, velz_modext1d, rho_modext1d, t_modext1d
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
! ... output to hdf5
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/ 1 /)
integer(hsize_t), dimension(1) :: dims_z
!
! ... local arrays
!
!-------------------------output to hdf5-file---------------------------
!
dims_z = (/ nz /)
!
write(*,*) 'save atmosphere to file ', trim(model_dir)//'/'//trim(model1d_file)
write(*,*)
 
call h5open_f (err)
   call h5fcreate_f(trim(model_dir)//'/'//trim(model1d_file), h5f_acc_trunc_f, file_id, err)
      CALL h5gcreate_f(file_id, 'dimensions', group_id, err)
         CALL h5screate_simple_f(1, dims_scalars, aspace_id, err)
            CALL h5acreate_f(group_id, 'nz', h5t_native_integer, aspace_id, &
                             attr_id, err)
               CALL h5awrite_f(attr_id, h5t_native_integer, nz, dims_scalars, err)
            CALL h5aclose_f(attr_id, err)
         CALL h5sclose_f(aspace_id, err)
      CALL h5gclose_f(group_id, err)
!
      CALL h5gcreate_f(file_id, 'coordinates', group_id, err)
         CALL h5screate_simple_f(1, dims_z, dspace_id, err)
            CALL h5dcreate_f(group_id, 'z', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, z_modext1d, dims_z, err)
            CALL h5dclose_f(dset_id, err)
         CALL h5sclose_f(dspace_id, err)
      CALL h5gclose_f(group_id, err)
!
      CALL h5gcreate_f(file_id, 'model', group_id, err)
         CALL h5screate_simple_f(1, dims_z, dspace_id, err)
            CALL h5dcreate_f(group_id, 'rho', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, rho_modext1d, dims_z, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'velz', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, velz_modext1d, dims_z, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'temperature', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, t_modext1d, dims_z, err)
            CALL h5dclose_f(dset_id, err)
         CALL h5sclose_f(dspace_id, err)
     CALL h5gclose_f(group_id, err)
!
   CALL h5fclose_f(file_id, err)
CALL h5close_f(err)
!
!
end subroutine output_mod1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_mod2d
!
!-----------------save 2d model atmosphere as h5-file--------------------
!
use prog_type
use mod_directories, only: model_dir, model2d_file
use dime_modext, only: nx, nz
use model2d, only: x_modext2d, z_modext2d, velx_modext2d, vely_modext2d, velz_modext2d, &
                   rho_modext2d, t_modext2d
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
! ... output to hdf5
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x, dims_z
integer(hsize_t), dimension(2) :: dims
!
! ... local arrays
!
!-------------------------output to hdf5-file---------------------------
!
dims_x = (/ nx /)
dims_z = (/ nz /)
dims = (/ nx, nz /)
!
write(*,*) 'save atmosphere to file ', trim(model_dir)//'/'//trim(model2d_file)
write(*,*)
 
call h5open_f (err)
   call h5fcreate_f(trim(model_dir)//'/'//trim(model2d_file), h5f_acc_trunc_f, file_id, err)
      CALL h5gcreate_f(file_id, 'dimensions', group_id, err)
         CALL h5screate_simple_f(1, dims_scalars, aspace_id, err)
            CALL h5acreate_f(group_id, 'nx', h5t_native_integer, aspace_id, &
                             attr_id, err)
               CALL h5awrite_f(attr_id, h5t_native_integer, nx, dims_scalars, err)
            CALL h5aclose_f(attr_id, err)
            CALL h5acreate_f(group_id, 'nz', h5t_native_integer, aspace_id, &
                             attr_id, err)
               CALL h5awrite_f(attr_id, h5t_native_integer, nz, dims_scalars, err)
            CALL h5aclose_f(attr_id, err)

         CALL h5sclose_f(aspace_id, err)
      CALL h5gclose_f(group_id, err)
!
      CALL h5gcreate_f(file_id, 'coordinates', group_id, err)
         CALL h5screate_simple_f(1, dims_x, dspace_id, err)
            CALL h5dcreate_f(group_id, 'x', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, x_modext2d, dims_x, err)
            CALL h5dclose_f(dset_id, err)
         CALL h5sclose_f(dspace_id, err)
         CALL h5screate_simple_f(1, dims_z, dspace_id, err)
            CALL h5dcreate_f(group_id, 'z', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, z_modext2d, dims_z, err)
            CALL h5dclose_f(dset_id, err)
         CALL h5sclose_f(dspace_id, err)
      CALL h5gclose_f(group_id, err)
!
      CALL h5gcreate_f(file_id, 'model', group_id, err)
         CALL h5screate_simple_f(2, dims, dspace_id, err)
            CALL h5dcreate_f(group_id, 'rho', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, rho_modext2d, dims, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'velx', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, velx_modext2d, dims, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'vely', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, vely_modext2d, dims, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'velz', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, velz_modext2d, dims, err)
            CALL h5dclose_f(dset_id, err)               
            CALL h5dcreate_f(group_id, 'temperature', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, t_modext2d, dims, err)
            CALL h5dclose_f(dset_id, err)
         CALL h5sclose_f(dspace_id, err)

     CALL h5gclose_f(group_id, err)
   CALL h5fclose_f(file_id, err)
CALL h5close_f(err)
!
!
end subroutine output_mod2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_mod3d
!
!-----------------save 3d model atmosphere as h5-file--------------------
!
use prog_type
use mod_directories, only: model_dir, model3d_file
use dime_modext, only: nx, ny, nz
use model3d, only: x_modext3d, y_modext3d, z_modext3d, &
                   velx_modext3d, vely_modext3d, velz_modext3d, &
                   rho_modext3d, tgas_modext3d, trad_modext3d
use params_input, only: unit_velocity, unit_length, unit_density, unit_temperature
use hdf5
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
! ... output to hdf5
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z
integer(hsize_t), dimension(3) :: dims
!
! ... local arrays
!
!-------------------------output to hdf5-file---------------------------
!
dims_x = (/ nx /)
dims_y = (/ ny /)
dims_z = (/ nz /)
dims = (/ nx, ny, nz /)
!
write(*,*) 'save atmosphere to file ', trim(model_dir)//'/'//trim(model3d_file)
write(*,*)
 
call h5open_f (err)
   call h5fcreate_f(trim(model_dir)//'/'//trim(model3d_file), h5f_acc_trunc_f, file_id, err)
      CALL h5gcreate_f(file_id, 'dimensions', group_id, err)
         CALL h5screate_simple_f(1, dims_scalars, aspace_id, err)
            CALL h5acreate_f(group_id, 'nx', h5t_native_integer, aspace_id, &
                             attr_id, err)
               CALL h5awrite_f(attr_id, h5t_native_integer, nx, dims_scalars, err)
            CALL h5aclose_f(attr_id, err)
            CALL h5acreate_f(group_id, 'ny', h5t_native_integer, aspace_id, &
                             attr_id, err)
               CALL h5awrite_f(attr_id, h5t_native_integer, ny, dims_scalars, err)
            CALL h5aclose_f(attr_id, err)
            CALL h5acreate_f(group_id, 'nz', h5t_native_integer, aspace_id, &
                             attr_id, err)
               CALL h5awrite_f(attr_id, h5t_native_integer, nz, dims_scalars, err)
            CALL h5aclose_f(attr_id, err)
         CALL h5sclose_f(aspace_id, err)
      CALL h5gclose_f(group_id, err)
!
      CALL h5gcreate_f(file_id, 'coordinates', group_id, err)
         CALL h5screate_simple_f(1, dims_x, dspace_id, err)
            CALL h5dcreate_f(group_id, 'x', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, x_modext3d, dims_x, err)
            CALL h5dclose_f(dset_id, err)
         CALL h5sclose_f(dspace_id, err)
         CALL h5screate_simple_f(1, dims_y, dspace_id, err)
            CALL h5dcreate_f(group_id, 'y', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, y_modext3d, dims_y, err)
            CALL h5dclose_f(dset_id, err)
         CALL h5sclose_f(dspace_id, err)
         CALL h5screate_simple_f(1, dims_z, dspace_id, err)
            CALL h5dcreate_f(group_id, 'z', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, z_modext3d, dims_z, err)
            CALL h5dclose_f(dset_id, err)
         CALL h5sclose_f(dspace_id, err)
      CALL h5gclose_f(group_id, err)
!
      CALL h5gcreate_f(file_id, 'model', group_id, err)
         CALL h5screate_simple_f(3, dims, dspace_id, err)
            CALL h5dcreate_f(group_id, 'rho', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, rho_modext3d, dims, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'velx', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, velx_modext3d, dims, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'vely', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, vely_modext3d, dims, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'velz', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, velz_modext3d, dims, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'tgas', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, tgas_modext3d, dims, err)
            CALL h5dclose_f(dset_id, err)
            CALL h5dcreate_f(group_id, 'trad', h5t_native_double, dspace_id, &
                             dset_id, err)
               CALL h5dwrite_f(dset_id, h5t_native_double, trad_modext3d, dims, err)
            CALL h5dclose_f(dset_id, err)               
         CALL h5sclose_f(dspace_id, err)
         CALL h5gclose_f(group_id, err)


      call h5gcreate_f(file_id, 'units', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'unit_length', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, unit_length, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'unit_velocity', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, unit_velocity, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
            call h5acreate_f(group_id, 'unit_density', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, unit_density, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'unit_temperature', h5t_native_double, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, unit_temperature, dims_scalars, err)
            call h5aclose_f(attr_id, err)               
         call h5sclose_f(aspace_id, err)
      call h5gclose_f(group_id, err)
         
   CALL h5fclose_f(file_id, err)
CALL h5close_f(err)
!
!
end subroutine output_mod3d
