module mod_io
  !
  use prog_type
  use fund_const
  !
  implicit none
  !
  !
contains
  !
  !-----------------------------------------------------------------------
  !-------------------read in all input parameters------------------------
  !-----------------------------------------------------------------------
  !
  subroutine read_input(input_terminal)
    !
    use options, only: opt_ng_cont, opt_ait_cont, opt_alo_cont, opt_method, opt_opac, opt_epsc
    use params_input, only: eps_cont, kcont, yhe, hei, verbose
    use mod_frequencies, only: xnue0, nnue, lfreqint, opt_grey
    use mod_grid3d, only: nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, opt_gridxyz
    use mod_bcondition3d, only: opt_bcondition, xic1_nue, xic2_nue
    use mod_angles, only: ntheta, opt_angint_method
    use mod_benchmark, only: n_z, n_y, thetab, phib, benchmark_mod
    use mod_directories, only: input_file
    use mod_model3d, only: model_dir
    use mod_opal, only: opal_dir
    !
    implicit none
    !
    ! ... arguments
    logical, intent(in) :: input_terminal   !decide whether indat file name is read from terminal or specified beforehand
    !
    ! ... local scalars
    integer(i4b) :: err
    real(dp) :: theta, phi
    integer(i4b) :: input_mod
    logical :: ver
    !
    ! ... local characters
    !
    ! ... local functions
    !
    ! ... namelist
    namelist / input_model / yhe, hei, input_mod
    namelist / input_cont / eps_cont, kcont
    namelist / dimensions_freq / nnue, xnue0
    namelist / dimensions_angles / ntheta
    namelist / dimensions_3dx / nx, xmin, xmax
    namelist / dimensions_3dy / ny, ymin, ymax
    namelist / dimensions_3dz / nz, zmin, zmax
    namelist / input_options / opt_ng_cont, opt_ait_cont, opt_alo_cont, opt_angint_method, model_dir, opal_dir, opt_grey, opt_opac, opt_epsc, opt_gridxyz, verbose
    namelist / input_sc3d / opt_method
    namelist / input_bcondition3d / opt_bcondition, xic1_nue, xic2_nue
    namelist / benchmark / benchmark_mod, theta, phi
    !
    !
    !-------------------set default values----------------------------------
    !
    opt_ng_cont=.true.
    opt_ait_cont=.false.
    opt_alo_cont=3
    opt_angint_method=-1
    !model_dir set in mod_model3d.f90
    !opal_dir set in mod_opal.f90
    opt_grey=-1
    opt_opac=-1
    opt_epsc=-1
    opt_gridxyz=1
    verbose=.true.
    !
    opt_method=-1
    !
    benchmark_mod=0
    theta=zero
    phi=zero
    !    
    yhe = 0.1
    hei = 2.
    input_mod = -1
    !
    eps_cont = zero
    kcont = one
    !
    nnue = -1
    xnue0 = zero
    !
    ntheta = -1
    !
    nx = -1
    ny = -1
    nz = -1
    xmin=zero
    xmax=zero
    ymin=zero
    ymax=zero
    zmin=zero
    zmax=zero
    !
    opt_bcondition = 0
    !
    !-----------------------------------------------------------------------
    !    
    !
    !
    write(*,*) '----------------------------read input-----------------------------------------'
    if(input_terminal) then
       write(*,*) 'input file name (*.nml) to define model-atmosphere'
       read(*,*) input_file
    endif
    write(*,*) 'reading input from file: ', trim(input_file)
    write(*,*)
    !
    !-----------------------------------------------------------------------
    !
    open(1, file=trim(input_file), status='old', form='formatted')
    !
    rewind 1
    read(1, nml=input_options)
    ver = verbose
    !
    !read 2d model parameters
    rewind 1
    read(1, nml=input_model)
    !
    !--------------------read dimensions for 2d grid------------------------
    !
    rewind 1
    read(1, nml=dimensions_3dx)
    !
    rewind 1
    read(1, nml=dimensions_3dy)
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
    read(1, nml=input_sc3d)
    !
    !----------------read boundary condition options etc. ------------------
    !
    allocate(xic1_nue(nnue), stat=err)
       if(err.ne.0) stop 'error in read_input: allocation'
    allocate(xic2_nue(nnue), stat=err)
       if(err.ne.0) stop 'error in read_input: allocation'       
    
    rewind 1
    read(1, nml=input_bcondition3d)    
    !
    !--------read benchmark models------------------------------------------
    !
    rewind 1
    read(1, nml=benchmark)
    thetab=theta*pi/180.d0
    phib=phi*pi/190.d0
    !
    close(1)
    !
    !
    select case(opt_method)
       case(4)
       case(5)
       case(6)
       case(7)
       case(14)
       case(15)
       case(16)
       case(17)
       case default
          opt_method=4
    end select
    !
    if(ver) write(*,*) 'input opt_method: ', opt_method
    if(ver) write(*,*) '   0 - '
    if(ver) write(*,*) '   1 - '
    if(ver) write(*,*) '   2 - '
    if(ver) write(*,*) '   3 - '
    if(ver) write(*,*) '   4 - 1st order sc'
    if(ver) write(*,*) '   5 - 2nd order sc'
    if(ver) write(*,*) '   6 - 1st order lc'
    if(ver) write(*,*) '   7 - 2nd order lc'
    if(ver) write(*,*) '   8 - '
    if(ver) write(*,*) '   9 - '
    if(ver) write(*,*) '  14 - LTE source function + 1st order sc radiation moments'
    if(ver) write(*,*) '  15 - LTE source function + 2nd order sc radiation moments'
    if(ver) write(*,*) '  16 - LTE source function + 1st order lc radiation moments'
    if(ver) write(*,*) '  17 - LTE source function + 2nd order lc radiation moments'
    if(ver) write(*,*)
    !
    if(ver) write(*,*) 'input opt_alo: ', opt_alo_cont
    if(ver) write(*,*) '   0 - classical lambda iteration (alo=0)'
    if(ver) write(*,*) '   1 - diagonal alo'
    if(ver) write(*,*) '   2 - tridiagonal alo'
    if(ver) write(*,*)
    !
    if(ver) write(*,*) 'input opt_ng: ', opt_ng_cont
    if(ver) write(*,*) '   0 - NG extrapolation off'
    if(ver) write(*,*) '   1 - NG extrapolation on'
    if(ver) write(*,*)
    !
    !
    !
    if(opt_grey.ne.0.and.&
         opt_grey.ne.1.and.&
         opt_grey.ne.2) opt_grey=1
    
    !
    if(ver) write(*,*) 'input opt_grey: ', opt_grey
    if(ver) write(*,*) '   0 - frequency dependent radiative transfer for multi frequencies'
    if(ver) write(*,*) '   1 - grey radiative transfer at one specific frequency bin'
    if(ver) write(*,*) '   2 - grey radiative transfer for frequency integrated quantities'
    if(ver) write(*,*)
    !
    select case(opt_grey)
       case(0,1)
          lfreqint = .false.
       case(2)
          lfreqint = .true.
       case default
          stop 'error in read_input: opt_ltec not properly specified'
    end select
    if(opt_grey.eq.0) stop 'error in read_input: opt_grey=0 to be implemented'
    !
    !
    !
    if(opt_opac.ne.0.and.&
         opt_opac.ne.1) opt_opac=0
    !
    if(ver) write(*,*) 'input opt_opac: ', opt_opac
    if(ver) write(*,*) '   0 - analytic opacity laws'
    if(ver) write(*,*) '   1 - OPAL opacities'
    if(ver) write(*,*)

    
    if(opt_epsc.ne.0.and.&
         opt_epsc.ne.1) opt_epsc=0
    !
    if(ver) write(*,*) 'input opt_epsc: ', opt_epsc
    if(ver) write(*,*) '   0 - constant eps_cont from input'
    if(ver) write(*,*) '   1 - calculating eps_cont from opacities'
    if(ver) write(*,*)
    !
    !
    if(opt_gridxyz.ne.0.and.&
         opt_gridxyz.ne.1) opt_gridxyz=1
    !
    if(ver) write(*,*) 'input opt_gridxyz: ', opt_gridxyz
    if(ver) write(*,*) '   0 - linear spacing for x, y, z'
    if(ver) write(*,*) '   1 - linear spacing for x, y, logarithmic spacing for z'
    if(ver) write(*,*)
    !
    !
    !
    select case(opt_bcondition)
       case(0)
       case default
          stop 'error in read_input: opt_bcondition3d not implemented yet'
    end select       
    if(ver) write(*,*) 'input opt_bcondition: ', opt_bcondition
    if(ver) write(*,*) '   0 - I_core(nue) = xic1(nue) - mu*xic2(nue)'
    if(ver) write(*,*) '   1 - set xic1= (1-eps)*B(Trad) + eps*B(Tgas), xic2=1/chi * dB/dz'
    if(ver) write(*,*) '   2 - read from user specified routine'
    !
  end subroutine read_input
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine print_model(nueindx, verbose)
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, &
         rho3d, opac3d, bnue3d, tgas3d, trad3d, velx3d, vely3d, velz3d, eps_cont3d, bnue3d
    use params_input, only: unit_length_cgs
    !
    implicit none
    !
    ! ... arguments
    integer(i4b), intent(in) :: nueindx
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... locals calars
    integer(i4b) :: i, j, k
    !
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '--------------------------model along z-axis-----------------------------------'
    if(ver) write(*,*)
    !
    i=nx/2+1
    j=ny/2+1
    if(ver) write(*,'(a20,i5)') 'frequency point', nueindx
    if(ver) write(*,'(a20, 2es20.8)') 'along z at (x,y)= ', x(i), y(j)
    if(ver) write(*,'(10a20)') 'z [r_star]', 'rho [g/cm^3]', 'opac [1/unit_length]', &
         'velx [km/s]', 'vely [km/s]', 'velz [km/s]', 'Tgas [K]', 'Trad [K]', 'B [cgs]', 'eps_cont'
    do k=1, nz
       if(ver) write(*,'(10es20.8)') z(k), rho3d(i,j,k), opac3d(i,j,k), velx3d(i,j,k)/1.d5, &
            vely3d(i,j,k)/1.d5, velz3d(i,j,k)/1.d5, &
            tgas3d(i,j,k), trad3d(i,j,k), bnue3d(i,j,k), eps_cont3d(i,j,k)
    enddo
    if(ver) write(*,*)
    !
    !
    !
  end subroutine print_model
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine print_solution(nueindx, verbose)
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, &
         rho3d, opac3d, bnue3d, eps_cont3d, mint3d, fcontx3d, fconty3d, fcontz3d, &
         kcontxx3d, kcontxy3d, kcontxz3d, kcontyy3d, kcontyz3d, kcontzz3d
    use params_input, only: unit_length_cgs
    !
    implicit none
    !
    ! ... arguments
    integer(i4b), intent(in) :: nueindx
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... locals calars
    integer(i4b) :: i, j, k
    !
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '------------------------solution along z-axis----------------------------------'
    if(ver) write(*,*)
    !
    i=nx/2+1
    j=ny/2+1
    if(ver) write(*,'(a20,i5)') 'frequency point', nueindx
    if(ver) write(*,'(a20, 2es20.8)') 'along z at (x,y)= ', x(i), y(j)
    if(ver) write(*,'(16a12)') 'z [r_star]', 'rho [cgs]', 'kappa [cgs]', &
         'B [cgs]', 'eps_cont', 'J [cgs]', 'H_x [cgs]', 'H_y [cgs]', 'H_z [cgs]', &
         'K_xx [cgs]', 'K_yy [cgs]', 'K_zz [cgs]', &
         'K_xy [cgs]', 'K_xz [cgs]', 'K_yz [cgs]', &
         'f = K_zz/J'
    do k=1, nz
       if(ver) write(*,'(16es12.4)') z(k), rho3d(i,j,k), opac3d(i,j,k)/rho3d(i,j,k)/unit_length_cgs, bnue3d(i,j,k), eps_cont3d(i,j,k), &
            mint3d(i,j,k), fcontx3d(i,j,k), fconty3d(i,j,k), fcontz3d(i,j,k), &
            kcontxx3d(i,j,k), kcontyy3d(i,j,k), kcontzz3d(i,j,k), &
            kcontxy3d(i,j,k), kcontxz3d(i,j,k), kcontyz3d(i,j,k), &
            kcontzz3d(i,j,k)/mint3d(i,j,k)
    enddo
    if(ver) write(*,*)
    !
  end subroutine print_solution
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine output(verbose)
    !
    use options, only: opt_alo_cont, opt_ng_cont, opt_ait_cont, opt_method, opt_opac
    use mod_angles, only: nomega, n_x, n_y, n_z, opt_angint_method
    use mod_frequencies, only: nnue, nodes_nue, opt_grey
    use mod_grid3d, only: x, y, z, nx, ny, nz, rho3d, tgas3d, trad3d, velx3d, vely3d, velz3d, opac3d, scont3d, imask3d, imaskb3d, &
         mint3d, eps_cont3d, fcontx3d, fconty3d, fcontz3d, &
         kcontxx3d, kcontxy3d, kcontxz3d, &
         kcontyy3d, kcontyz3d, &
         kcontzz3d
    use params_input, only: unit_density, unit_velocity, unit_temperature, unit_length, yhe, hei, kcont, eps_cont
    use mod_conttrans3d, only: itmaxc, nconvc, epsmaxc_arr, devmaxc
    use hdf5
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: iopt_ng_cont, iopt_ait_cont
    !
    ! ... local arrays
    integer, dimension(:,:,:), allocatable :: mask3d
    !
    ! ... local characters
    character(len=100) :: output_file
    !
    ! ... for output to hdf5
    integer(i4b) :: err
    integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
    integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
    integer(hsize_t), dimension(1) :: dims_x , dims_y, dims_z, dims_omega, dims_nue, dims_itc
    integer(hsize_t), dimension(3) :: dims_3d
    !
    !-----------------------------------------------------------------------
    !
    if(present(verbose)) ver = verbose
    !
    output_file='outputFILES/sc3d/output_sc3d.h5'
    !
    if(ver) write(*,*) '-------------------------output to directory-----------------------------------'
    if(ver) write(*,*) 'output to file ', trim(output_file)
    if(ver) write(*,*)
    !
    !-----------------------------------------------------------------------
    !
    dims_x = (/ nx /)
    dims_y = (/ ny /)
    dims_z = (/ nz /)
    dims_omega = (/ nomega /)
    dims_nue = (/ nnue /)
    dims_itc = (/ itmaxc /)
    dims_3d = (/ nx, ny, nz /)
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
                call h5acreate_f(group_id, 'opt_angint_method', h5t_native_integer, aspace_id, attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, opt_angint_method, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'opt_grey', h5t_native_integer, aspace_id, attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, opt_grey, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'opt_opac', h5t_native_integer, aspace_id, attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, opt_opac, dims_scalars, err)
                call h5aclose_f(attr_id, err)               
             call h5sclose_f(aspace_id, err)
          call h5gclose_f(group_id, err)
    !
    !----------------------------dimensions---------------------------------
    !
          call h5gcreate_f(file_id, 'dimensions', group_id, err)
             call h5screate_simple_f(1, dims_scalars, aspace_id, err)
                call h5acreate_f(group_id, 'nx', h5t_native_integer, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, nx, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'ny', h5t_native_integer, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, ny, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'nz', h5t_native_integer, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, nz, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'nomega', h5t_native_integer, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, nomega, dims_scalars, err)
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
                call h5acreate_f(group_id, 'unit_density', h5t_native_real, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_double, unit_density, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'unit_velocity', h5t_native_real, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_double, unit_velocity, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'unit_length', h5t_native_real, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_double, unit_length, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'unit_temperature', h5t_native_real, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_double, unit_temperature, dims_scalars, err)
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
             call h5screate_simple_f(1, dims_y, dspace_id, err)
                call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
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
             call h5screate_simple_f(1, dims_omega, dspace_id, err)
                call h5dcreate_f(group_id, 'n_x', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, n_x, dims_omega, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'n_y', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, n_y, dims_omega, err)
                call h5dclose_f(dset_id, err)               
                call h5dcreate_f(group_id, 'n_z', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, n_z, dims_omega, err)
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
    !----------------------------3d solution--------------------------------
    !
         call h5gcreate_f(file_id, 'solution3d', group_id, err)
             call h5screate_simple_f(3, dims_3d, dspace_id, err)
                call h5dcreate_f(group_id, 'scont3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, scont3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'mint3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, mint3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'fcontx3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, fcontx3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'fconty3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, fconty3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'fcontz3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, fcontz3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'kcontxx3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, kcontxx3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'kcontxy3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, kcontxy3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'kcontxz3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, kcontxz3d, dims_3d, err)
                call h5dclose_f(dset_id, err)                              
                call h5dcreate_f(group_id, 'kcontyy3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, kcontyy3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'kcontyz3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, kcontyz3d, dims_3d, err)
                call h5dclose_f(dset_id, err)      
                call h5dcreate_f(group_id, 'kcontzz3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, kcontzz3d, dims_3d, err)
                call h5dclose_f(dset_id, err)               
             call h5sclose_f(dspace_id, err)
         call h5gclose_f(group_id, err)
    !
    !--------------------------------3d model-------------------------------
    !
         call h5gcreate_f(file_id, 'model3d', group_id, err)
             call h5screate_simple_f(3, dims_3d, dspace_id, err)
                call h5dcreate_f(group_id, 'mask3d', h5t_native_integer, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_integer, imask3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'maskb3d', h5t_native_integer, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_integer, imaskb3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
             call h5sclose_f(dspace_id, err)
    !
             call h5screate_simple_f(3, dims_3d, dspace_id, err)
                call h5dcreate_f(group_id, 'eps_cont3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, eps_cont3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'rho3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, rho3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'tgas3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, tgas3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'trad3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, trad3d, dims_3d, err)
                call h5dclose_f(dset_id, err)               
                call h5dcreate_f(group_id, 'opac3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'velx3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, velx3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'vely3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, vely3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'velz3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, velz3d, dims_3d, err)
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
  subroutine output_benchmark01(verbose)
    !
    use mod_grid3d, only: x, y, z, nx, ny, nz, int3d, opac3d, scont3d, imask3d, imaskb3d
    use mod_angles, only: nomega, n_x, n_y, n_z
    use mod_benchmark, only: int3d_theo
    use hdf5
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local ligcals
    logical :: ver=.false.
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
    integer(hsize_t), dimension(1) :: dims_x , dims_y, dims_z, dims_omega, dims_nue
    integer(hsize_t), dimension(3) :: dims_3d
    !
    !-----------------------------------------------------------------------
    !
    if(present(verbose)) ver=verbose
    !
    output_file='outputFILES/sc3d/searchlight3d.h5'
    if(ver) write(*,*) '-------------------------output to directory-----------------------------------'
    if(ver) write(*,*) 'output to file ', trim(output_file)
    if(ver) write(*,*)
    !
    !-----------------------------------------------------------------------
    !
    dims_x = (/ nx /)
    dims_y = (/ ny /)
    dims_z = (/ nz /)
    dims_3d = (/ nx, ny, nz /)
    dims_omega = (/ nomega /)
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
                call h5acreate_f(group_id, 'ndymax', h5t_native_integer, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, ny, dims_scalars, err)
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
             call h5screate_simple_f(1, dims_y, dspace_id, err)
                call h5dcreate_f(group_id, 'y', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
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
             call h5screate_simple_f(1, dims_omega, dspace_id, err)
                call h5dcreate_f(group_id, 'n_x', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, n_x, dims_omega, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'n_y', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, n_y, dims_omega, err)
                call h5dclose_f(dset_id, err)               
                call h5dcreate_f(group_id, 'n_z', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, n_z, dims_omega, err)
                call h5dclose_f(dset_id, err)
             call h5sclose_f(dspace_id, err)
             call h5gclose_f(group_id, err)
    !
    !-----------------frequency grid----------------------------------------
    !
    !----------------------convergence behaviour----------------------------
    !
    !----------------------------3d solution--------------------------------
    !
         call h5gcreate_f(file_id, 'solution3d', group_id, err)
             call h5screate_simple_f(3, dims_3d, dspace_id, err)
                call h5dcreate_f(group_id, 'int3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, int3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'int3d_theo', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, int3d_theo, dims_3d, err)
                call h5dclose_f(dset_id, err)
         call h5gclose_f(group_id, err)
    !
    !--------------------------------3d model-------------------------------
    !
         call h5gcreate_f(file_id, 'model3d', group_id, err)
             call h5screate_simple_f(3, dims_3d, dspace_id, err)
                call h5dcreate_f(group_id, 'mask3d', h5t_native_integer, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_integer, imask3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'maskb3d', h5t_native_integer, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_integer, imaskb3d, dims_3d, err)
                call h5dclose_f(dset_id, err)
             call h5sclose_f(dspace_id, err)
    !
             call h5screate_simple_f(3, dims_3d, dspace_id, err)
                call h5dcreate_f(group_id, 'opac3d', h5t_native_real, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, opac3d, dims_3d, err)
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


end module mod_io
    
    
    
    
    
