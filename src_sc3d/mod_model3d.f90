module mod_model3d
  
  use prog_type
  use fund_const

  implicit none

  !modelxd_file:   file where model-atmosphere is stored
  !model_dir:       directory of input models  
  character(len=300) :: model_dir='inputFILES'
  character(len=10), parameter :: model1d_file='model1d.h5'
  character(len=10), parameter :: model2d_file='model2d.h5'
  character(len=10), parameter :: model3d_file='model3d.h5'
  !
  !dimensions of the external model
  integer(i4b) :: nx_modext, ny_modext, nz_modext
  !
  !data of external model
  real(dp), dimension(:), allocatable :: x_modext, y_modext, z_modext
  real(dp), dimension(:,:), allocatable :: velx_modext2d, vely_modext2d, velz_modext2d, &
       rho_modext2d, tgas_modext2d, trad_modext2d, vth_modext2d  
  real(dp), dimension(:,:,:), allocatable :: velx_modext3d, vely_modext3d, velz_modext3d, &
       rho_modext3d, tgas_modext3d, trad_modext3d, vth_modext3d
  
contains  
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine setup_mod3d(imodel, verbose)
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, rho3d, tgas3d, trad3d, &
    velx3d, vely3d, velz3d, xmin, xmax, ymin, ymax, zmin, zmax
    use mod_interp1d, only: find_index
    use mod_interp3d, only: coeff3d_8p_lin
    !
    ! ... arguments
    integer(i4b), intent(in) :: imodel
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: err
    integer(i4b) :: i, j, k, iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
    !
    ! ... local logicals
    !
    ! ... local functions
    !
    !-----------------------------------------------------------------------
    !
    if(present(verbose)) ver=verbose
    !
    !read 3d model
    select case(imodel)
       case(1)
          call read_mod3d_hdf5(verbose=ver)    !read from hdf5 file
       case(2)
          call read_mod3d_amrvac2d(verbose=ver)  !take directly from within amrvac simulations
       case default
          stop 'error in setup_mod3d: read_mod3d method not specified'
    end select
    !
    if(xmin.lt.minval(x_modext)) stop 'error in setup_mod3d: input xmin smaller than min(x_modext)'
    if(xmax.gt.maxval(x_modext)) stop 'error in setup_mod3d: input xmax greater than max(x_modext)'
    if(ymin.lt.minval(y_modext)) stop 'error in setup_mod3d: input ymin smaller than min(y_modext)'
    if(ymax.gt.maxval(y_modext)) stop 'error in setup_mod3d: input ymax greater than max(y_modext)'
    if(zmin.lt.minval(z_modext)) stop 'error in setup_mod3d: input zmin smaller than min(z_modext)'
    if(zmax.gt.maxval(z_modext)) stop 'error in setup_mod3d: input zmax greater than max(z_modext)'
    !
    !--------calculate calculation volume and perform interpolations--------
    !
    if(allocated(rho3d)) deallocate(rho3d)
    if(allocated(velx3d)) deallocate(velx3d)
    if(allocated(vely3d)) deallocate(vely3d)
    if(allocated(velz3d)) deallocate(velz3d)    
    if(allocated(tgas3d)) deallocate(tgas3d)
    if(allocated(trad3d)) deallocate(trad3d)    
    !    
    allocate(rho3d(nx,ny,nz),stat=err)
       if(err.ne.0) stop 'error: allocation in setup_mod3d'
    allocate(velx3d(nx,ny,nz),stat=err)
       if(err.ne.0) stop 'error: allocation in setup_mod3d'
    allocate(vely3d(nx,ny,nz),stat=err)
       if(err.ne.0) stop 'error: allocation in setup_mod3d'
    allocate(velz3d(nx,ny,nz),stat=err)
       if(err.ne.0) stop 'error: allocation in setup_mod3d'
    allocate(tgas3d(nx,ny,nz),stat=err)
       if(err.ne.0) stop 'error: allocation in setup_mod3d'
    allocate(trad3d(nx,ny,nz),stat=err)
       if(err.ne.0) stop 'error: allocation in setup_mod3d'
    !
    !
    rho3d = zero
    velx3d = zero
    vely3d = zero
    velz3d = zero
    tgas3d = zero
    trad3d = zero
    !
    !interpolate everything outside the ghost zones
    do i=1, nx-1
      call find_index(x(i), x_modext, nx_modext, iim2, iim1, ii, iip1)
      do j=1, ny-1
        call find_index(y(j), y_modext, ny_modext, jjm2, jjm1, jj, jjp1)
        do k=1, nz-2
          call find_index(z(k), z_modext, nz_modext, kkm2, kkm1, kk, kkp1)

          call coeff3d_8p_lin(x_modext(iim1), x_modext(ii), &
               y_modext(jjm1), y_modext(jj), &
               z_modext(kkm1), z_modext(kk), &
               x(i), y(j), z(k), &
               acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff)

          rho3d(i,j,k) = acoeff*rho_modext3d(iim1,jjm1,kkm1) + bcoeff*rho_modext3d(ii,jjm1,kkm1) + &
               ccoeff*rho_modext3d(iim1,jj,kkm1)   + dcoeff*rho_modext3d(ii,jj,kkm1) + &
               ecoeff*rho_modext3d(iim1,jjm1,kk)   + fcoeff*rho_modext3d(ii,jjm1,kk) + &
               gcoeff*rho_modext3d(iim1,jj,kk)     + hcoeff*rho_modext3d(ii,jj,kk)
          
          tgas3d(i,j,k) = acoeff*tgas_modext3d(iim1,jjm1,kkm1) + bcoeff*tgas_modext3d(ii,jjm1,kkm1) + &
               ccoeff*tgas_modext3d(iim1,jj,kkm1)   + dcoeff*tgas_modext3d(ii,jj,kkm1) + &
               ecoeff*tgas_modext3d(iim1,jjm1,kk)   + fcoeff*tgas_modext3d(ii,jjm1,kk) + &
               gcoeff*tgas_modext3d(iim1,jj,kk)     + hcoeff*tgas_modext3d(ii,jj,kk)

          trad3d(i,j,k) = acoeff*trad_modext3d(iim1,jjm1,kkm1) + bcoeff*trad_modext3d(ii,jjm1,kkm1) + &
               ccoeff*trad_modext3d(iim1,jj,kkm1)   + dcoeff*trad_modext3d(ii,jj,kkm1) + &
               ecoeff*trad_modext3d(iim1,jjm1,kk)   + fcoeff*trad_modext3d(ii,jjm1,kk) + &
               gcoeff*trad_modext3d(iim1,jj,kk)     + hcoeff*trad_modext3d(ii,jj,kk)

          velx3d(i,j,k) = acoeff*velx_modext3d(iim1,jjm1,kkm1) + bcoeff*velx_modext3d(ii,jjm1,kkm1) + &
               ccoeff*velx_modext3d(iim1,jj,kkm1)   + dcoeff*velx_modext3d(ii,jj,kkm1) + &
               ecoeff*velx_modext3d(iim1,jjm1,kk)   + fcoeff*velx_modext3d(ii,jjm1,kk) + &
               gcoeff*velx_modext3d(iim1,jj,kk)     + hcoeff*velx_modext3d(ii,jj,kk)

          vely3d(i,j,k) = acoeff*vely_modext3d(iim1,jjm1,kkm1) + bcoeff*vely_modext3d(ii,jjm1,kkm1) + &
               ccoeff*vely_modext3d(iim1,jj,kkm1)   + dcoeff*vely_modext3d(ii,jj,kkm1) + &
               ecoeff*vely_modext3d(iim1,jjm1,kk)   + fcoeff*vely_modext3d(ii,jjm1,kk) + &
               gcoeff*vely_modext3d(iim1,jj,kk)     + hcoeff*vely_modext3d(ii,jj,kk)
          
          velz3d(i,j,k) = acoeff*velz_modext3d(iim1,jjm1,kkm1) + bcoeff*velz_modext3d(ii,jjm1,kkm1) + &
               ccoeff*velz_modext3d(iim1,jj,kkm1)   + dcoeff*velz_modext3d(ii,jj,kkm1) + &
               ecoeff*velz_modext3d(iim1,jjm1,kk)   + fcoeff*velz_modext3d(ii,jjm1,kk) + &
               gcoeff*velz_modext3d(iim1,jj,kk)     + hcoeff*velz_modext3d(ii,jj,kk)
        enddo
      enddo
    enddo
    !
    !
    !------------------------set the ghost zones----------------------------
    !
    !ghost zones
    do j=1, ny-1
      do k=3, nz-2
        !ghost zones right side
        rho3d(nx,j,k) = rho3d(1,j,k)
        tgas3d(nx,j,k) = tgas3d(1,j,k)
        trad3d(nx,j,k) = trad3d(1,j,k)
        velx3d(nx,j,k) = velx3d(1,j,k)
        vely3d(nx,j,k) = vely3d(1,j,k)
        velz3d(nx,j,k) = velz3d(1,j,k)
      enddo
    enddo
    !
    do i=1, nx-1
      do k=3, nz-2
        !ghost zones back side
        rho3d(i,ny,k) = rho3d(i,1,k)
        tgas3d(i,ny,k) = tgas3d(i,1,k)
        trad3d(i,ny,k) = trad3d(i,1,k)
        velx3d(i,ny,k) = velx3d(i,1,k)
        vely3d(i,ny,k) = vely3d(i,1,k)
        velz3d(i,ny,k) = velz3d(i,1,k)
      enddo
    enddo
    !
    do k=3, nz-2
      !ghost zones right back edge
      rho3d(nx,ny,k) = rho3d(1,1,k)
      tgas3d(nx,ny,k) = tgas3d(1,1,k)
      trad3d(nx,ny,k) = trad3d(1,1,k)
      velx3d(nx,ny,k) = velx3d(1,1,k)
      vely3d(nx,ny,k) = vely3d(1,1,k)
      velz3d(nx,ny,k) = velz3d(1,1,k)
    enddo
    !
    !
    !
    !ghost zones at bottom and top
    do i=1, nx
      do j=1, ny
        velx3d(i,j,2) = velx3d(i,j,3)
        vely3d(i,j,2) = vely3d(i,j,3)
        velz3d(i,j,2) = velz3d(i,j,3)
        rho3d(i,j,2) = rho3d(i,j,3)
        tgas3d(i,j,2) = tgas3d(i,j,3)
        trad3d(i,j,2) = trad3d(i,j,3)
        velx3d(i,j,1) = velx3d(i,j,2)
        vely3d(i,j,1) = vely3d(i,j,2)
        velz3d(i,j,1) = velz3d(i,j,2)
        rho3d(i,j,1) = rho3d(i,j,2)
        tgas3d(i,j,1) = tgas3d(i,j,2)
        trad3d(i,j,1) = trad3d(i,j,2)
        !
        velx3d(i,j,nz-1) = velx3d(i,j,nz-2)
        vely3d(i,j,nz-1) = vely3d(i,j,nz-2)
        velz3d(i,j,nz-1) = velz3d(i,j,nz-2)
        rho3d(i,j,nz-1) = rho3d(i,j,nz-2)
        tgas3d(i,j,nz-1) = tgas3d(i,j,nz-2)
        trad3d(i,j,nz-1) = trad3d(i,j,nz-2)
        velx3d(i,j,nz) = velx3d(i,j,nz-1)
        vely3d(i,j,nz) = vely3d(i,j,nz-1)
        velz3d(i,j,nz) = velz3d(i,j,nz-1)
        rho3d(i,j,nz) = rho3d(i,j,nz-1)
        tgas3d(i,j,nz) = tgas3d(i,j,nz-1)
        trad3d(i,j,nz) = trad3d(i,j,nz-1)
      enddo
    enddo
    !
    deallocate(rho_modext3d)
    deallocate(tgas_modext3d)
    deallocate(trad_modext3d)
    deallocate(velx_modext3d)
    deallocate(vely_modext3d)
    deallocate(velz_modext3d)                
    !
  end subroutine setup_mod3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine setup_opacities1(nueindx, verbose)
    !
    use params_input, only: kcont, yhe, hei, unit_length_cgs
    use mod_grid3d, only: nx, ny, nz, z, rho3d, opac3d, zmax, zmin
    !
    ! ... arguments
    integer(i4b), intent(in) :: nueindx    
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    real(dp) :: rho, dtau, a, b
    !
    ! ... local logicals
    !
    ! ... local functions
    real(dp) :: opac_thomson
    !
    if(present(verbose)) ver=verbose
    !
    if(ver) write(*,*) '-------------------setting up opacities: analytical law------------------------'
    if(ver) write(*,*)
    !
    do i=1, nx
      do j=1, ny
        do k=1, nz
          rho=rho3d(i,j,k)
          !in units 1/unit_length
          opac3d(i,j,k) = opac_thomson(yhe, hei, rho, kcont)*unit_length_cgs
        enddo
      enddo
    enddo
    !
    !
    !use constant opacity (k=1 corresponds to tau=1) for test cases
    !dtau=one
    !opac3d=kcont*dtau/(zmax-zmin)
    !!
    !!exponentially decreasing opacity
    !b=log(tau_max/tau_min)/(zmax-zmin)
    !a=b*tau_max*exp(b*zmin)
    !do i=1, nx
    !   do j=1, ny
    !      do k=1, nz
    !         opac3d(i,j,k)=a*exp(-b*z(k))
    !      enddo
    !   enddo
    !enddo
    !
  end subroutine setup_opacities1
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine setup_opacities2(verbose)
    !
    use params_input, only: kcont
    use params_input, only: yhe, hei, unit_length_cgs
    use mod_grid3d, only: nx, ny, nz, z, rho3d, opac3d, tgas3d, trad3d
    use mod_opal
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    real(dp) :: rho, temp
    !
    ! ... local logicals
    !
    ! ... local functions
    real(dp) :: opac_opal
    !
    if(present(verbose)) ver=verbose
    !
    if(ver) write(*,*) '--------------------setting up opacities: OPAL tables--------------------------'
    if(ver) write(*,*)
    !
    call get_opal_table(yhe, verbose=ver)
    !
    do i=1, nx
      do j=1, ny
        !do i=nx/2+1, nx/2+1
        !   do j=ny/2+1, ny/2+1
        do k=1, nz
          rho=rho3d(i,j,k)
          !         temp=trad3d(i,j,k)
          temp=tgas3d(i,j,k)
          !in units 1/unit_length_cgs
          opac3d(i,j,k) = opac_opal(kcont, yhe, hei, log10(rho),log10(temp), nrho_opal, ntemp_opal, rho_opal, temp_opal, kappa_opal)*unit_length_cgs

        enddo
      enddo
    enddo
    !
    !write(*,*) yhe, hei
    !stop 'go on in setup_opacities'
    !
    !
  end subroutine setup_opacities2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine setup_epsc1(verbose)
    !
    use params_input, only: kcont, eps_cont
    use mod_grid3d, only: eps_cont3d
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    real(dp) :: rho, temp
    !
    ! ... local logicals
    !
    ! ... local functions
    real(dp) :: opac_opal
    !
    if(present(verbose)) ver=verbose
    !
    if(ver) write(*,*) '--------------------setting up 3d thermalization parameter: constant-----------'
    if(ver) write(*,*)
    !
    !
    eps_cont3d = eps_cont
    !
    !
  end subroutine setup_epsc1
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine setup_epsc2(verbose)
    !
    use params_input, only: kcont, eps_cont, unit_length_cgs, yhe, hei
    use mod_grid3d, only: nx, ny, nz, eps_cont3d, rho3d, opac3d
    !
    ! ... arguments    
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    real(dp) :: rho, opac_sc, opac_tot
    !
    ! ... local logicals
    !
    ! ... local functions
    real(dp) :: opac_thomson
    !
    if(present(verbose)) ver=verbose
    !
    if(ver) write(*,*) '--------------------setting up 3d thermalization parameter: constant-----------'
    if(ver) write(*,*)
    !
    if(ver) write(*,*) yhe, hei, kcont
    !
    do k=1, nz
      do j=1, ny
        do i=1, nx
          rho=rho3d(i,j,k)
          !in units 1/unit_length
          opac_sc = opac_thomson(yhe, hei, rho, kcont)*unit_length_cgs
          opac_tot = opac3d(i,j,k)
          eps_cont3d(i,j,k) = (opac_tot-opac_sc)/opac_tot
          if(ver) write(*,*) opac_tot/rho/unit_length_cgs, opac_sc/rho/unit_length_cgs, eps_cont3d(i,j,k)
          if(eps_cont3d(i,j,k).lt.zero) stop 'error in setup_epsc2: eps_cont3d(i,j,k) < 0'
        enddo
      enddo
    enddo

    stop
    !
    !
    !
  end subroutine setup_epsc2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine setup_bnue(nueindx, verbose)
    !
    use mod_grid3d, only: nx, ny, nz, tgas3d, trad3d, bnue3d, eps_cont3d
    use mod_frequencies, only: nodes_nue, lfreqint, opt_grey
    use mod_math, only: bnue2
    !
    ! ... arguments
    integer(i4b), intent(in) :: nueindx    
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    !
    ! ... local logicals
    !
    ! ... local functions
    !
    if(present(verbose)) ver=verbose
    !
    if(ver) write(*,*) '-----------------------setting up planck function------------------------------'
    if(ver) write(*,*)
    !
    if(.not.allocated(bnue3d)) then
      allocate(bnue3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error: allocation in setup_bnue'
    endif
    !
    do i=1, nx
      do j=1, ny
        do k=1, nz
          !         bnue3d(i,j,k)=bnue2(nodes_nue(nueindx),trad3d(i,j,k),lfreqint)
          !         bnue3d(i,j,k)=bnue2(nodes_nue(nueindx),tgas3d(i,j,k),lfreqint)
          bnue3d(i,j,k) = &
          (one-eps_cont3d(i,j,k))*bnue2(nodes_nue(nueindx),trad3d(i,j,k),lfreqint) &
          + eps_cont3d(i,j,k)*bnue2(nodes_nue(nueindx),tgas3d(i,j,k),lfreqint)
        enddo
      enddo
    enddo

    !
  end subroutine setup_bnue
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine read_mod3d_hdf5(verbose)
    !
    use params_input, only: unit_velocity, unit_temperature, unit_density, unit_length, unit_length_cgs
    use hdf5
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: i, err
    ! ... for hdf5-file
    integer(hid_t) :: file_id, group_id, attr_id, dset_id
    integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z, dims_scalars
    integer(hsize_t), dimension(3) :: dims
    !
    !----------------read in radius and theta grid from 3d file-------------
    !
    if(present(verbose)) ver=verbose
    !
    if(ver) write(*,*) '-----------------read model atmosphere from 3d input file----------------------'
    if(ver) write(*,*) 'file name: ', trim(model_dir)//'/'//model3d_file
    if(ver) write(*,*)
    !
    if(allocated(x_modext)) deallocate(x_modext)
    if(allocated(y_modext)) deallocate(y_modext)
    if(allocated(z_modext)) deallocate(z_modext)
    if(allocated(rho_modext3d)) deallocate(rho_modext3d)
    if(allocated(tgas_modext3d)) deallocate(tgas_modext3d)
    if(allocated(trad_modext3d)) deallocate(trad_modext3d)
    if(allocated(velx_modext3d)) deallocate(velx_modext3d)
    if(allocated(vely_modext3d)) deallocate(vely_modext3d)
    if(allocated(velz_modext3d)) deallocate(velz_modext3d)
    !
    dims_scalars=(/ 1 /)
    !
    call h5open_f(err)
    call h5fopen_f(trim(model_dir)//'/'//model3d_file, h5f_acc_rdonly_f, file_id, err)
    !
    !read units
    call h5gopen_f(file_id, 'units', group_id, err)
    call h5aopen_f(group_id, 'unit_length', attr_id, err)
    call h5aread_f(attr_id, h5t_native_double, unit_length, dims_scalars, err)
    call h5aclose_f(attr_id, err)
    call h5aopen_f(group_id, 'unit_density', attr_id, err)
    call h5aread_f(attr_id, h5t_native_double, unit_density, dims_scalars, err)
    call h5aclose_f(attr_id, err)
    call h5aopen_f(group_id, 'unit_temperature', attr_id, err)
    call h5aread_f(attr_id, h5t_native_double, unit_temperature, dims_scalars, err)
    call h5aclose_f(attr_id, err)
    call h5aopen_f(group_id, 'unit_velocity', attr_id, err)
    call h5aread_f(attr_id, h5t_native_double, unit_velocity, dims_scalars, err)
    call h5aclose_f(attr_id, err)
    call h5gclose_f(group_id, err)
    unit_length_cgs = unit_length*rsu
    !
    !read dimensions
    call h5gopen_f(file_id, 'dimensions', group_id, err)
    call h5aopen_f(group_id, 'nx', attr_id, err)
    call h5aread_f(attr_id, h5t_native_integer, nx_modext, dims_scalars, err)
    call h5aclose_f(attr_id, err)
    call h5aopen_f(group_id, 'ny', attr_id, err)
    call h5aread_f(attr_id, h5t_native_integer, ny_modext, dims_scalars, err)
    call h5aclose_f(attr_id, err)
    call h5aopen_f(group_id, 'nz', attr_id, err)
    call h5aread_f(attr_id, h5t_native_integer, nz_modext, dims_scalars, err)
    call h5aclose_f(attr_id, err)

    call h5gclose_f(group_id, err)
    !
    !read coordinates
    dims_x=(/nx_modext/)
    dims_y=(/ny_modext/)
    dims_z=(/nz_modext/)
    dims=(/ nx_modext, ny_modext, nz_modext /)
    !
    allocate(x_modext(nx_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    allocate(y_modext(ny_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    allocate(z_modext(nz_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    allocate(rho_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    allocate(tgas_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    allocate(trad_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    allocate(velx_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    allocate(vely_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    allocate(velz_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
    if(err.ne.0) stop 'error in read_mod3d: allocation'
    !
    call h5gopen_f(file_id, 'coordinates', group_id, err)
    call h5dopen_f(group_id, 'x', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, x_modext, dims_x, err)
    call h5dclose_f(dset_id, err)
    call h5dopen_f(group_id, 'y', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, y_modext, dims_y, err)
    call h5dclose_f(dset_id, err)
    call h5dopen_f(group_id, 'z', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, z_modext, dims_z, err)
    call h5dclose_f(dset_id, err)
    call h5gclose_f(group_id, err)
    !
    call h5gopen_f(file_id, 'model', group_id, err)
    call h5dopen_f(group_id, 'rho', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, rho_modext3d, dims, err)
    call h5dclose_f(dset_id, err)
    call h5dopen_f(group_id, 'velx', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, velx_modext3d, dims, err)
    call h5dclose_f(dset_id, err)
    call h5dopen_f(group_id, 'vely', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, vely_modext3d, dims, err)
    call h5dclose_f(dset_id, err)
    call h5dopen_f(group_id, 'velz', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, velz_modext3d, dims, err)
    call h5dclose_f(dset_id, err)
    call h5dopen_f(group_id, 'tgas', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, tgas_modext3d, dims, err)
    call h5dclose_f(dset_id, err)
    call h5dopen_f(group_id, 'trad', dset_id, err)
    call h5dread_f(dset_id, h5t_native_double, trad_modext3d, dims, err)
    call h5dclose_f(dset_id, err)
    call h5gclose_f(group_id, err)
    !
    call h5fclose_f(file_id, err)
    call h5close_f(err)
    !
    !transform everything to cgs
    rho_modext3d = rho_modext3d*unit_density
    velx_modext3d = velx_modext3d*unit_velocity
    vely_modext3d = vely_modext3d*unit_velocity
    velz_modext3d = velz_modext3d*unit_velocity
    tgas_modext3d = tgas_modext3d*unit_temperature
    trad_modext3d = trad_modext3d*unit_temperature
    !
    !
  end subroutine read_mod3d_hdf5
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine read_mod3d_amrvac2d(verbose)
    !
    !note: global 2d arrays have to be set within amrvac
    !
    use params_input, only: unit_length, unit_density, unit_temperature, unit_velocity
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver=.false.
    !
    ! ... local scalars
    integer(i4b) :: j, err
    !
    if(present(verbose)) ver=verbose    
    !
    !----------------read in radius and theta grid from 3d file-------------
    !
    if(ver) write(*,*) '-----------------read model atmosphere from 2d amrvac model--------------------'
    if(ver) write(*,*)
    !
    !
    if(allocated(rho_modext3d)) deallocate(rho_modext3d)
    if(allocated(tgas_modext3d)) deallocate(tgas_modext3d)
    if(allocated(trad_modext3d)) deallocate(trad_modext3d)
    if(allocated(velx_modext3d)) deallocate(velx_modext3d)
    if(allocated(vely_modext3d)) deallocate(vely_modext3d)
    if(allocated(velz_modext3d)) deallocate(velz_modext3d)
    !
    !
    allocate(rho_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in read_mod3d_amrvac: allocation'
    allocate(tgas_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in read_mod3d_amrvac: allocation'
    allocate(trad_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in read_mod3d_amrvac: allocation'
    allocate(velx_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in read_mod3d_amrvac: allocation'
    allocate(vely_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in read_mod3d_amrvac: allocation'
    allocate(velz_modext3d(nx_modext, ny_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in read_mod3d_amrvac: allocation'    
    !
    !
    !transform everything to cgs
    do j=1, ny_modext
       rho_modext3d(:,j,:) = rho_modext2d
       velx_modext3d(:,j,:) = velx_modext2d
       vely_modext3d(:,j,:) = vely_modext2d
       velz_modext3d(:,j,:) = velz_modext2d
       tgas_modext3d(:,j,:) = tgas_modext2d
       trad_modext3d(:,j,:) = trad_modext2d
    enddo

    rho_modext3d = rho_modext3d*unit_density
    velx_modext3d = velx_modext3d*unit_velocity
    vely_modext3d = vely_modext3d*unit_velocity
    velz_modext3d = velz_modext3d*unit_velocity
    tgas_modext3d = tgas_modext3d*unit_temperature
    trad_modext3d = trad_modext3d*unit_temperature
    
  end subroutine read_mod3d_amrvac2d
  

end module mod_model3d
