subroutine output
!
use prog_type
use options, only: opt_alo_cont, opt_ng_cont, opt_ait_cont, opt_method
use angles, only: dim_mu, nodes_mu, n_x, n_z
use freq, only: nnue, nodes_nue, xic1_nue
use dime2d, only: x, z, nx, nz, t2d, velz2d, opac2d, scont2d, opac2d, imask2d, imaskb2d
use params_input, only: teff, trad, tmin, rstar, xmloss, vmin, vmax, beta, yhe, hei, kcont, eps_cont
use iter, only: itmaxc, nconvc, epsmaxc_arr, devmaxc
use hdf5
!
implicit none
!
! ... arguments
!
! ... local scalars
integer(i4b) :: i, j, k
integer(i4b) :: iopt_ng_cont, iopt_ait_cont
!
! ... local arrays
integer, dimension(:,:), allocatable :: mask2d
!
! ... local characters
character(len=100) :: output_file
!
! ... for output to hdf5
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x , dims_z, dims_mu, dims_nue, dims_itc
integer(hsize_t), dimension(2) :: dims_2d
!
!-----------------------------------------------------------------------
!
output_file='outputFILES/sc2d/model2d.h5'
write(*,*) '-------------------------output to directory-----------------------------------'
write(*,*) 'output to file ', trim(output_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
dims_x = (/ nx /)
dims_z = (/ nz /)
dims_mu = (/ dim_mu /)
dims_nue = (/ nnue /)
dims_itc = (/ itmaxc /)
dims_2d = (/ nx, nz /)


!write(*,*) itmaxc
!write(*,*) dims_itc
!write(*,*) epsmaxc_arr
!
!--------convert all logicals to integers (to be read in idl)-----------
!
!convert options
iopt_ng_cont=0
iopt_ait_cont=0
if(opt_ng_cont) iopt_ng_cont=1
if(opt_ait_cont) iopt_ait_cont=1
!
!-----------------------------------------------------------------------
! 
call h5open_f(err)
   call h5fcreate_f(trim(output_file), h5f_acc_trunc_f, file_id, err)
!
!------------------------------options----------------------------------
!
      call h5gcreate_f(file_id, 'options', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'opt_method', h5t_native_integer, aspace_id, attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, opt_method, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'opt_ng_cont', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, iopt_ng_cont, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'opt_ait_cont', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, iopt_ait_cont, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'opt_alo_cont', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, opt_alo_cont, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
      call h5gclose_f(group_id, err)
!
!----------------------------dimensions---------------------------------
!
      call h5gcreate_f(file_id, 'dimensions', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'ndxmax', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nx, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'ndzmax', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nz, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'dim_mu', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, dim_mu, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'nnue', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nnue, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
      call h5gclose_f(group_id, err)
!
!----------------------------input parameters---------------------------
!
      call h5gcreate_f(file_id, 'input_parameters', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'kcont', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, kcont, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'eps_cont', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, eps_cont, dims_scalars, err)
            call h5aclose_f(attr_id, err)

            call h5acreate_f(group_id, 'teff', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, teff, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'trad', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, trad, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'rstar', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, rstar, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'yhe', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, yhe, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'hei', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, hei, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'mdot', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, xmloss, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmin', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmin, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'vmax', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, vmax, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'beta', h5t_native_real, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_double, beta, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
     call h5gclose_f(group_id, err)
!
!------------------------boundary condition-----------------------------
!
!----------------------------coordinates--------------------------------
!
      call h5gcreate_f(file_id, 'coordinates', group_id, err)
         call h5screate_simple_f(1, dims_x, dspace_id, err)
            call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
         call h5screate_simple_f(1, dims_z, dspace_id, err)
            call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!---------------------angular grids-------------------------------------
!
      call h5gcreate_f(file_id, 'angles', group_id, err)
         call h5screate_simple_f(1, dims_mu, dspace_id, err)
            call h5dcreate_f(group_id, 'nodes_mu', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, nodes_mu, dims_mu, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'n_x', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, n_x, dims_mu, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'n_z', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, n_z, dims_mu, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!-----------------frequency grid----------------------------------------
!
      call h5gcreate_f(file_id, 'frequencies', group_id, err)
         call h5screate_simple_f(1, dims_nue, dspace_id, err)
            call h5dcreate_f(group_id, 'nodes_nue', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, nodes_nue, dims_nue, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------convergence behaviour----------------------------
!
     call h5gcreate_f(file_id, 'convergence_behaviour', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'itmaxc', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, itmaxc, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'devmaxc', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, devmaxc, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'nconvc', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nconvc, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
!
         call h5screate_simple_f(1, dims_itc, dspace_id, err)
            call h5dcreate_f(group_id, 'epsmaxc_arr', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, epsmaxc_arr, dims_itc, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------------2d solution--------------------------------
!
     call h5gcreate_f(file_id, 'solution2d', group_id, err)
         call h5screate_simple_f(2, dims_2d, dspace_id, err)
            call h5dcreate_f(group_id, 'scont2d', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, scont2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'mint2d', h5t_native_real, dspace_id, &
                             dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!--------------------------------3d model-------------------------------
!
     call h5gcreate_f(file_id, 'model2d', group_id, err)
         call h5screate_simple_f(2, dims_2d, dspace_id, err)
            call h5dcreate_f(group_id, 'mask2d', h5t_native_integer, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_integer, imask2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'maskb2d', h5t_native_integer, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_integer, imaskb2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
!
         call h5screate_simple_f(2, dims_2d, dspace_id, err)
            call h5dcreate_f(group_id, 't2d', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, t2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'opac2d', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, opac2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'velz2d', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, velz2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!------------------------------debugging--------------------------------
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
end subroutine output
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine output_benchmark01
!
use prog_type
use dime2d, only: x, z, nx, nz, int2d, opac2d, scont2d, imask2d, imaskb2d
use angles, only: dim_mu, nodes_mu, n_x, n_z
use mod_benchmark, only: int2d_theo, itmaxi, devmaxi, epsmaxi_arr, nconvi
use hdf5
!
implicit none
!
! ... arguments
!
! ... local scalars
!
! ... local arrays
!
! ... local characters
character(len=100) :: output_file
!
! ... for output to hdf5
integer(i4b) :: err
integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
integer(hsize_t), dimension(1) :: dims_x , dims_z, dims_mu, dims_nue, dims_itc
integer(hsize_t), dimension(2) :: dims_2d
!
!-----------------------------------------------------------------------
!
output_file='outputFILES/sc2d/searchlight2d.h5'
write(*,*) '-------------------------output to directory-----------------------------------'
write(*,*) 'output to file ', trim(output_file)
write(*,*)
!
!-----------------------------------------------------------------------
!
dims_x = (/ nx /)
dims_z = (/ nz /)
dims_2d = (/ nx, nz /)
dims_itc = (/ itmaxi /)
dims_mu = (/ dim_mu /)
!
!-----------------------------------------------------------------------
! 
call h5open_f(err)
   call h5fcreate_f(trim(output_file), h5f_acc_trunc_f, file_id, err)
!
!------------------------------options----------------------------------
!
!----------------------------dimensions---------------------------------
!
      call h5gcreate_f(file_id, 'dimensions', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'ndxmax', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nx, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'ndzmax', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nz, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
      call h5gclose_f(group_id, err)
!
!----------------------------input parameters---------------------------
!
!------------------------boundary condition-----------------------------
!
!----------------------------coordinates--------------------------------
!
      call h5gcreate_f(file_id, 'coordinates', group_id, err)
         call h5screate_simple_f(1, dims_x, dspace_id, err)
            call h5dcreate_f(group_id, 'x', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
         call h5screate_simple_f(1, dims_z, dspace_id, err)
            call h5dcreate_f(group_id, 'z', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!---------------------angular grids-------------------------------------
!
      call h5gcreate_f(file_id, 'angles', group_id, err)
         call h5screate_simple_f(1, dims_mu, dspace_id, err)
            call h5dcreate_f(group_id, 'nodes_mu', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, nodes_mu, dims_mu, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'n_x', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, n_x, dims_mu, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'n_z', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, n_z, dims_mu, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
         call h5gclose_f(group_id, err)
!
!-----------------frequency grid----------------------------------------
!
!----------------------convergence behaviour----------------------------
!
     call h5gcreate_f(file_id, 'convergence_behaviour', group_id, err)
         call h5screate_simple_f(1, dims_scalars, aspace_id, err)
            call h5acreate_f(group_id, 'itmaxi', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, itmaxi, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'devmaxi', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, devmaxi, dims_scalars, err)
            call h5aclose_f(attr_id, err)
            call h5acreate_f(group_id, 'nconvi', h5t_native_integer, aspace_id, &
                             attr_id, err)
               call h5awrite_f(attr_id, h5t_native_integer, nconvi, dims_scalars, err)
            call h5aclose_f(attr_id, err)
         call h5sclose_f(aspace_id, err)
!
         call h5screate_simple_f(1, dims_itc, dspace_id, err)
            call h5dcreate_f(group_id, 'epsmaxi_arr', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, epsmaxi_arr, dims_itc, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!----------------------------2d solution--------------------------------
!
     call h5gcreate_f(file_id, 'solution2d', group_id, err)
         call h5screate_simple_f(2, dims_2d, dspace_id, err)
            call h5dcreate_f(group_id, 'int2d', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, int2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'int2d_theo', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, int2d_theo, dims_2d, err)
            call h5dclose_f(dset_id, err)
     call h5gclose_f(group_id, err)
!
!--------------------------------3d model-------------------------------
!
     call h5gcreate_f(file_id, 'model2d', group_id, err)
         call h5screate_simple_f(2, dims_2d, dspace_id, err)
            call h5dcreate_f(group_id, 'mask2d', h5t_native_integer, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_integer, imask2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
            call h5dcreate_f(group_id, 'maskb2d', h5t_native_integer, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_integer, imaskb2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
!
         call h5screate_simple_f(2, dims_2d, dspace_id, err)
            call h5dcreate_f(group_id, 'opac2d', h5t_native_real, dspace_id, &
                             dset_id, err)
               call h5dwrite_f(dset_id, h5t_native_double, opac2d, dims_2d, err)
            call h5dclose_f(dset_id, err)
         call h5sclose_f(dspace_id, err)
     call h5gclose_f(group_id, err)
!
!------------------------------debugging--------------------------------
!
   call h5fclose_f(file_id, err)
call h5close_f(err)
!
!
end subroutine output_benchmark01
