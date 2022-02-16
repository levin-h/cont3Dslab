module mod_sc3damrvac

  use prog_type
  use fund_const
  
  implicit none

contains

  subroutine sc3d_init(fname, unit_l, unit_d, unit_t, unit_v)
    !
    !on input (from amrvac main file, unit_length, unit_density, unit_temperature, unit_velocity
    !
    use mod_directories, only: input_file
    use mod_angles, only: calcnodes_omega
    use mod_frequencies, only: calcnodes_nue
    use mod_grid3d, only: gridxyz
    use mod_io, only: read_input
    use params_input, only: verbose
    !
    ! ... arguments
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: unit_l, unit_d, unit_t, unit_v
    !
    !
    !store parameter file in global variable
    input_file=trim(fname)

    !print the logo
    call sc3d_logo

    !read the input from namelist
    call read_input(.false.)

    !create the coordinate grid
    call gridxyz(verbose=verbose)

    !create the angular integration grid
    call calcnodes_omega(verbose=verbose)

    !create the frequency grid
    call calcnodes_nue(verbose=verbose)
    !
    !set the units
    call set_units(unit_l, unit_d, unit_t, unit_v)

  end subroutine sc3d_init
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine sc3d_logo

    write(*,*)
    write(*,*) '-----------------------------------------------------------------------------' 
    write(*,*) '-----------------------------------------------------------------------------'
    write(*,*) '|        _____  _____     ____  _____     _____ _               ____        |'
    write(*,*) '|       / ____|/ ____|   |___ \|  __ \   / ____| |        /\   |  _ \       |'
    write(*,*) '|      | (___ | |   ______ __) | |  | | | (___ | |       /  \  | |_) |      |'
    write(*,*) '|       \___ \| |  |______|__ <| |  | |  \___ \| |      / /\ \ |  _ <       |'
    write(*,*) '|       ____) | |____     ___) | |__| |  ____) | |____ / ____ \| |_) |      |'
    write(*,*) '|      |_____/ \_____|   |____/|_____/  |_____/|______/_/    \_\____/       |'
    write(*,*) '-----------------------------------------------------------------------------'
    write(*,*) '-----------------------------------------------------------------------------'    
    write(*,*)
    

  end subroutine sc3d_logo
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine set_units(unit_l, unit_d, unit_t, unit_v)
    !
    use params_input, only: unit_length, unit_length_cgs, unit_density, unit_temperature, unit_velocity
    !
    real(dp), intent(in) :: unit_l, unit_d, unit_t, unit_v
    !
    write(*,*) '-------------------------setting the units to sc3d scheme----------------------'
    write(*,*)
    unit_length=unit_l/rsu
    unit_length_cgs = unit_l
    unit_density=unit_d
    unit_temperature=unit_t
    unit_velocity=unit_v
    !
  end subroutine set_units
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine set_model2d(nx, ny, x, y, rho2d, velx2d, vely2d, tgas2d, trad2d)
    !
    use mod_model3d, only: nx_modext, ny_modext, nz_modext, x_modext, y_modext, z_modext, &
         rho_modext2d, velx_modext2d, vely_modext2d, velz_modext2d, &
         tgas_modext2d, trad_modext2d
    use mod_bcondition3d, only: xic1_nue, xic2_nue
    !
    ! ... arguments
    integer(i4b), intent(in) :: nx, ny
    real(dp), dimension(nx), intent(in) :: x
    real(dp), dimension(ny), intent(in) :: y
    real(dp), dimension(nx,ny), intent(in) :: rho2d, velx2d, vely2d, tgas2d, trad2d
    !
    ! ... local scalars
    integer(i4b) :: i, j, err
    !    
    if(allocated(x_modext)) deallocate(x_modext)
    if(allocated(y_modext)) deallocate(y_modext)
    if(allocated(z_modext)) deallocate(z_modext)    
    if(allocated(rho_modext2d)) deallocate(rho_modext2d)
    if(allocated(tgas_modext2d)) deallocate(tgas_modext2d)
    if(allocated(trad_modext2d)) deallocate(trad_modext2d)
    if(allocated(velx_modext2d)) deallocate(velx_modext2d)
    if(allocated(vely_modext2d)) deallocate(vely_modext2d)
    if(allocated(velz_modext2d)) deallocate(velz_modext2d)

    !note that in own code, z corresponds to x and x corresponds to y in amrvac
    nx_modext = ny
    nz_modext = nx
    ny_modext = ny
    allocate(x_modext(nx_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    allocate(y_modext(ny_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    allocate(z_modext(nz_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    
    allocate(rho_modext2d(nx_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    allocate(tgas_modext2d(nx_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    allocate(trad_modext2d(nx_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    allocate(velx_modext2d(nx_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    allocate(vely_modext2d(nx_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    allocate(velz_modext2d(nx_modext, nz_modext), stat=err)
       if(err.ne.0) stop 'error in set_model2d: allocation'
    !
    x_modext = y
    y_modext = y
    z_modext = x

    do i=1, nx_modext
       do j=1, nz_modext
          rho_modext2d(i,j) = rho2d(j,i)
          tgas_modext2d(i,j) = tgas2d(j,i)
          trad_modext2d(i,j) = trad2d(j,i)
          velx_modext2d(i,j) = vely2d(j,i)   !again, x in sc3d corresponds to y in amrac
          vely_modext2d(i,j) = zero          !for the moment, no y velocity component
          velz_modext2d(i,j) = velx2d(j,i)
       enddo
    enddo

  end subroutine set_model2d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine get_model2d(nx_amrvac, ny_amrvac, x_amrvac, y_amrvac, rho2d_amrvac, kappa2d_amrvac, erad2d_amrvac, fcontx2d_amrvac, fconty2d_amrvac, kcontxx2d_amrvac, kcontyy2d_amrvac, kcontxy2d_amrvac)
    !
    !rho2d only to test interpolations
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, rho3d, opac3d, mint3d, fcontx3d, fconty3d, fcontz3d, &
                          kcontxx3d, kcontxy3d, kcontxz3d, &
                          kcontyy3d, kcontyz3d, kcontzz3d
    use params_input, only: unit_length_cgs
    use mod_interp1d, only: find_index
    use mod_interp2d, only: coeff2d_4p_lin
    !
    ! ... arguments
    integer(i4b), intent(in) :: nx_amrvac, ny_amrvac
    real(dp), dimension(nx_amrvac), intent(in) :: x_amrvac
    real(dp), dimension(ny_amrvac), intent(in) :: y_amrvac
    real(dp), dimension(nx_amrvac,ny_amrvac), intent(inout) :: kappa2d_amrvac, fcontx2d_amrvac, fconty2d_amrvac, erad2d_amrvac, rho2d_amrvac, kcontxx2d_amrvac, kcontyy2d_amrvac, kcontxy2d_amrvac
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: iim2, iim1, ii, iip1, &
         kkm2, kkm1, kk, kkp1, jj

    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
    !
    !note that in own code, z corresponds to x and x corresponds to y in amrvac
    !
    jj=ny/2
    do i=1, nx_amrvac
       call find_index(x_amrvac(i), z, nz, kkm2, kkm1, kk, kkp1)
       do j=1, ny_amrvac
          call find_index(y_amrvac(j), x, nx, iim2, iim1, ii, iip1)

          call coeff2d_4p_lin(x(iim1), x(ii), z(kkm1), z(kk), y_amrvac(j), x_amrvac(i), &
               acoeff, bcoeff, ccoeff, dcoeff)

          rho2d_amrvac(i,j) = acoeff*rho3d(iim1,jj,kkm1) + bcoeff*rho3d(ii,jj,kkm1) + ccoeff*rho3d(iim1,jj,kk) + dcoeff*rho3d(ii,jj,kk)
          kappa2d_amrvac(i,j) = acoeff*opac3d(iim1,jj,kkm1) + bcoeff*opac3d(ii,jj,kkm1) + ccoeff*opac3d(iim1,jj,kk) + dcoeff*opac3d(ii,jj,kk)          
          erad2d_amrvac(i,j) = acoeff*mint3d(iim1,jj,kkm1) + bcoeff*mint3d(ii,jj,kkm1) + ccoeff*mint3d(iim1,jj,kk) + dcoeff*mint3d(ii,jj,kk)          
          fcontx2d_amrvac(i,j) = acoeff*fcontz3d(iim1,jj,kkm1) + bcoeff*fcontz3d(ii,jj,kkm1) + ccoeff*fcontz3d(iim1,jj,kk) + dcoeff*fcontz3d(ii,jj,kk)
          fconty2d_amrvac(i,j) = acoeff*fcontx3d(iim1,jj,kkm1) + bcoeff*fcontx3d(ii,jj,kkm1) + ccoeff*fcontx3d(iim1,jj,kk) + dcoeff*fcontx3d(ii,jj,kk)
          kcontxx2d_amrvac(i,j) = acoeff*kcontzz3d(iim1,jj,kkm1) + bcoeff*kcontzz3d(ii,jj,kkm1) + ccoeff*kcontzz3d(iim1,jj,kk) + dcoeff*kcontzz3d(ii,jj,kk)
          kcontyy2d_amrvac(i,j) = acoeff*kcontxx3d(iim1,jj,kkm1) + bcoeff*kcontxx3d(ii,jj,kkm1) + ccoeff*kcontxx3d(iim1,jj,kk) + dcoeff*kcontxx3d(ii,jj,kk)
          kcontxy2d_amrvac(i,j) = acoeff*kcontxz3d(iim1,jj,kkm1) + bcoeff*kcontxz3d(ii,jj,kkm1) + ccoeff*kcontxz3d(iim1,jj,kk) + dcoeff*kcontxz3d(ii,jj,kk)          
       enddo
    enddo

    !transform opacity to kappa in cgs units (opacity was given in 1/unit_length_cgs)
    kappa2d_amrvac = kappa2d_amrvac/unit_length_cgs/rho2d_amrvac
    !
    !transform intensity moments to energy, flux, pressure tensor
    erad2d_amrvac = erad2d_amrvac * four*pi/cgs_clight
    kcontxx2d_amrvac = kcontxx2d_amrvac * four*pi/cgs_clight
    kcontyy2d_amrvac = kcontyy2d_amrvac * four*pi/cgs_clight
    kcontxy2d_amrvac = kcontxy2d_amrvac * four*pi/cgs_clight    
    fcontx2d_amrvac = fcontx2d_amrvac * four*pi
    fconty2d_amrvac = fconty2d_amrvac * four*pi
    !
    
!    stop 'go on in get_model2d'

  end subroutine get_model2d  
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine run_sc3d(f_bound)
    !
    use mod_model3d, only: setup_mod3d
    use mod_grid3d, only: allocate_global3d, check_grid3d, tgas3d
    use options, only: opt_opac, opt_epsc, opt_method
    use mod_model3d, only: setup_opacities1, setup_opacities2, setup_bnue, setup_epsc1, setup_epsc2
    use mod_io, only: print_model, print_solution, output
    use mod_conttrans3d, only: conttrans_sc3d, conttrans_lte3d
    use mod_bcondition3d, only: setup_bcondition3d
    use params_input, only: verbose
    !
    integer(i4b) :: nueindx
    !
    real(dp), intent(in) :: f_bound
    !
    !
    !setup the model from amrvac input (option 2)
    call setup_mod3d(2, verbose=verbose)
    !
    !allocate all global arrays
    call allocate_global3d(verbose=verbose)
    !
    !set the frequency index
    nueindx=1
    !
    !set up opacities
    select case(opt_opac)
       case(0)
          call setup_opacities1(nueindx, verbose=verbose)
       case(1)
          call setup_opacities2(verbose=verbose)
       case default
         stop 'error in run_sc3d: opt_opac not set'
   end select
   !
   !set up 3d thermalization parameter
   select case(opt_epsc)
      case(0)
         call setup_epsc1(verbose=verbose)
      case(1)
         call setup_epsc2(verbose=verbose)
      case default
         stop 'error in main: opt_opac not set'
   end select
   !
   !set up planck function
   call setup_bnue(nueindx, verbose=verbose)
   !
   call setup_bcondition3d_amrvac(nueindx, f_bound, verbose=verbose)   
   !
   !print out the model
   call print_model(nueindx, verbose=verbose)
   !
   !run short characteristics method
   select case(opt_method)
      case(4,5,6,7)
      !check the grid and perform radiative transfer
         call check_grid3d(nueindx, verbose=verbose)
         call conttrans_sc3d(nueindx, verbose=verbose)
      case(14,15,16,17)
         !simply set the source function to lte value
         call conttrans_lte3d(nueindx, verbose=verbose)
      case default
         stop 'error in run_sc3d: opt_method not properly defined'
   end select
   !
   !
   call print_solution(nueindx, verbose=verbose)

!   call output(verbose=verbose)
   
!   stop 'go on in run_sc3d'
   

 end subroutine run_sc3d


 subroutine setup_bcondition3d_amrvac(nueindx, f_bound, verbose)
   !
   use mod_grid3d, only: nx, ny, bnue3d, x, y, z
   use mod_frequencies, only: nnue
   use mod_bcondition3d
   !
   ! ... arguments
   integer(i4b), intent(in) :: nueindx
   real(dp), intent(in) :: f_bound
   logical, intent(in), optional :: verbose
   !
   ! ... local logicals
   logical :: ver = .false.
   !
   ! ... locals calars
   integer(i4b) :: i, j, err
   !
   !
   if(present(verbose)) ver = verbose
   !
   
   if(ver) write(*,*) '-----------------------setting up boundary conditions--------------------------'
   if(ver) write(*,*)
   !
   if(allocated(xic1_3d)) deallocate(xic1_3d)
   if(allocated(xic2_3d)) deallocate(xic2_3d)
   allocate(xic1_3d(nx,ny,3), stat=err)  !only 3 inner points
      if(err.ne.0) stop 'error in setup_bcondition3d: allocation'
   allocate(xic2_3d(nx,ny,3), stat=err)  !only 3 inner points
      if(err.ne.0) stop 'error in setup_bcondition3d: allocation'

   !
   if(nnue.gt.1) stop 'error in setup_bcondition3d_amrvac: nnue>1 not implemented yet with f_bound)'   
   !set xic1 = bnue (already set beforehand at this frequency index)
   do i=1, nx
      do j=1, ny
         xic1_3d(i,j,1) = bnue3d(i,j,1)
         xic1_3d(i,j,2) = bnue3d(i,j,2)

         xic2_3d(i,j,1) =  - three * f_bound/four/pi
         xic2_3d(i,j,2) =  - three * f_bound/four/pi                  
      enddo
   enddo
   !
   i=nx/2+1
   j=ny/2+1
   if(ver) write(*,'(a20,i5)') 'frequency point', nueindx
   if(ver) write(*,'(a20, 2es20.8)') 'at (x,y)= ', x(i), y(j)
   if(ver) write(*,'(10a20)') 'z [r_star]', 'xic1 [cgs]', 'xic2 [cgs]'
   if(ver) write(*,'(10es20.8)') z(1), xic1_3d(i,j,1), xic2_3d(i,j,1)
   if(ver) write(*,'(10es20.8)') z(2), xic1_3d(i,j,2), xic2_3d(i,j,2)
   if(ver) write(*,*)

 end subroutine setup_bcondition3d_amrvac
 !
end module mod_sc3damrvac  
