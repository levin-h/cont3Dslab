!
!-----------------------------------------------------------------------
!----------------------------for benchmarks-----------------------------
!-----------------------------------------------------------------------
!
module mod_benchmark3d
  !
  use prog_type
  use fund_const
  
  implicit none
  !
  integer(i4b) :: benchmark_mod
  !
  !for searchlight beam test
  real(dp) :: n_y, n_z, nn_x, nn_y, nn_z, thetab, phib
  real(dp), dimension(:,:,:), allocatable :: int3d_theo
  !
  !for plane-parallel diffusion
  real(dp), parameter :: tau_min=1.d-3, tau_max=1.d6

contains
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine make_benchmark(verbose)
    !
    use mod_benchmark, only: benchmark_mod
    use options, only: opt_method
    use mod_io, only: output_benchmark01
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i
    real(dp) :: ts, te
    !
    !calculate delx, dely, delz grids if not done yet
    !
    if(present(verbose)) ver=verbose
    !
    select case(benchmark_mod)
       case(1)
          if(ver) write(*,*) '-----------performing benchmark model 1: searchligh beam test 3d---------------'
          if(ver) write(*,*)
          !calculating solution
          call cpu_time(ts)
          call benchmark01_solution(verbose=ver)
          call cpu_time(te)
          if(ver) write(*,*) 'computation time benchmark01', te-ts
          if(ver) write(*,*)
          !output to file
          call output_benchmark01(verbose=ver)
          !
       case(2)
          if(ver) write(*,*) '-----performing benchmark model 2: LC test for periodic boundary condition-----'
          if(ver) write(*,*)
          !calculating solution
          call cpu_time(ts)
          call benchmark02_solution(verbose=ver)
          call cpu_time(te)
          if(ver) write(*,*) 'computation time benchmark02', te-ts
          if(ver) write(*,*)
          !output to file
          !call output_benchmark02(verbose=ver)
          !
       case default
          if(ver) write(*,*) '----------------------no benchmark is being performed--------------------------'
          if(ver) write(*,*)
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
  subroutine benchmark01_solution(verbose)
    !
    use mod_grid3d, only: z, int3d, scont3d, opac3d, alocont_o_nn3d, alocont_nn3d_tmp, &
         normalization3d_tmp, mint3d_tmp, fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
         nx, ny, nz, tgas3d, trad3d, bnue3d, zmin, zmax
    use mod_angles, only: nomega, n_z, n_y, n_x
    use options, only: opt_method
    use mod_frequencies, only: nodes_nue, xnue0
    use mod_benchmark, only: int3d_theo, nn_x, nn_y, nn_z, thetab, phib
    use params_input, only: kcont
    use mod_model3d, only: setup_bnue
    use mod_conttrans3d, only: formallc_cont3d, formalsc_cont3d_lin, formallc_cont3d_lin, formalsc_cont3d
    !
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    
    ! ... local scalars
    integer(i4b) :: i, j, k, err, ix, iz, oindx, nueindx
    integer(i4b) :: startz, endz, gamma
    integer(i4b) :: s1, s2, s3, s4, s4b
    real(dp) :: mu, eps_max, dtau, ts, te, norm
    !
    ! ... local arrays
    real(dp), dimension(:,:), allocatable :: eps2d
    real(dp), dimension(:,:), allocatable :: intbound2d_ng, intbound2d_single
    logical(dp), dimension(:,:), allocatable :: imaskbound2d
    !
    ! ... local functions
    real(dp) :: bnue
    !
    ! ... local characters
    !
    if(present(verbose)) ver=verbose
    !
    call cpu_time(ts)
    !
    !use constant opacity (k=1 corresponds to tau=1) for test cases
    dtau=one
    opac3d=kcont*dtau/(zmax-zmin)
    !
    nn_x = sin(thetab)*cos(phib)
    nn_y = sin(thetab)*sin(phib)
    nn_z = cos(thetab)
    norm = sqrt(nn_x**2+nn_y**2+nn_z**2)
    nn_x = nn_x/norm
    nn_y = nn_y/norm
    nn_z = nn_z/norm
    n_x=nn_x
    n_z=nn_z
    n_y=nn_y
    oindx=1
    !
    nodes_nue=xnue0
    nueindx=1

    if(.not.allocated(int3d)) allocate(int3d(nx,ny,nz))
    if(.not.allocated(alocont_o_nn3d)) allocate(alocont_o_nn3d(nx,ny,nz,27))
    if(.not.allocated(alocont_nn3d_tmp)) allocate(alocont_nn3d_tmp(nx,ny,nz,27))
    if(.not.allocated(normalization3d_tmp)) allocate(normalization3d_tmp(nx,ny,nz))
    if(.not.allocated(mint3d_tmp)) allocate(mint3d_tmp(nx,ny,nz))
    if(.not.allocated(fcontx3d_tmp)) allocate(fcontx3d_tmp(nx,ny,nz))
    if(.not.allocated(fconty3d_tmp)) allocate(fconty3d_tmp(nx,ny,nz))
    if(.not.allocated(fcontz3d_tmp)) allocate(fcontz3d_tmp(nx,ny,nz))
    !
    !scont3d=zero
    !calculating planck function
    call setup_bnue(nueindx, verbose=ver)
    !
    select case(opt_method)
       case(4)
          call formalsc_cont3d_lin(oindx,nueindx)
       case(5)
          call formalsc_cont3d(oindx,nueindx)
       case(6)
          call formallc_cont3d_lin(oindx,nueindx)
       case(7)
          call formallc_cont3d(oindx,nueindx)
       case default
          stop 'error in benchmark01_solution: opt_method not specified'
    end select
    !
    deallocate(alocont_nn3d_tmp)
    deallocate(normalization3d_tmp)
    deallocate(mint3d_tmp)
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
    
    allocate(int3d_theo(nx,ny,nz),stat=err)
    int3d_theo=zero
    int3d_theo(:,:,startz-gamma)=int3d(:,:,startz-gamma)
    int3d_theo(:,:,startz-2*gamma)=int3d(:,:,startz-2*gamma)
    !
    i=nx/2+1
    j=ny/2+1
    do k=startz, endz, gamma
       dtau = half*(opac3d(i,j,k)+opac3d(i,j,k-gamma))*(z(k)-z(k-gamma))/n_z(oindx)
       int3d_theo(:,:,k)=int3d_theo(:,:,k-gamma)*exp(-dtau)
    enddo
    !
  end subroutine benchmark01_solution
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine benchmark02_solution(verbose)
    !
    use mod_grid3d, only: x, y, z, int3d, scont3d, opac3d, alocont_o_nn3d, alocont_nn3d_tmp, &
         normalization3d_tmp, mint3d_tmp, fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
         nx, ny, nz, tgas3d, trad3d, bnue3d, zmin, zmax
    use mod_angles, only: nomega, n_z, n_y, n_x
    use options, only: opt_method
    use mod_frequencies, only: nodes_nue, xnue0
    use mod_benchmark, only: int3d_theo, nn_x, nn_y, nn_z, thetab, phib
    use params_input, only: kcont
    use mod_conttrans3d, only: formallc_cont3d, formalsc_cont3d_lin, formallc_cont3d_lin, formalsc_cont3d
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b), parameter :: nslab = 16
    integer(i4b) :: startx, endx, startxb, endxb, alpha, &
         starty, endy, startyb, endyb, beta, &
         startz, endz, gamma
    integer(i4b) :: k, ii, jj, kk, iip1, jjp1, kkp1
    integer(i4b) :: xindx, yindx, zindx, oindx, nueindx
    real(dp) :: x_p, y_p, z_p, x_d, y_d, z_d
    real(dp) :: int_p, opac_p, scont_p, opac_d, scont_d
    real(dp) :: dels_xyd, dels_xzd, dels_yzd, dels_d
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff
    real(dp) :: dtau, opacmax, norm
    !
    ! ... local arrays
    !
    ! ... local functions
    real(dp) :: bnue
    !
    ! ... local characters
    !
    write(*,*) 'remember to uncomment some output in flc_cont3d'
    stop

    
    if(.not.allocated(int3d)) allocate(int3d(nx,ny,nz))
    int3d=zero
    !
    nn_x = sin(thetab)*cos(phib)
    nn_y = sin(thetab)*sin(phib)
    nn_z = cos(thetab)
    norm = sqrt(nn_x**2+nn_y**2+nn_z**2)
    nn_x = nn_x/norm
    nn_y = nn_y/norm
    nn_z = nn_z/norm
    n_x=nn_x
    n_z=nn_z
    n_y=nn_y
    !
    oindx=1
    !
    nodes_nue=xnue0
    nueindx=1
    !
    !for a point
    xindx=nx!1!nx/2+1!nx-2
    yindx=ny/2+1!ny-2
    zindx=nz/2+1
    
    xindx=3
    yindx=1
    zindx=13
    !
    !zero source function
    scont3d=zero
    !
    !-----------------------------opacity law-------------------------------
    !
    !use constant opacity (k=1 corresponds to tau=1) for test cases
    !dtau=one
    !opac3d=kcont*dtau/(zmax-zmin)
    !
    !linear opacity law along z, such that delta_tau along ray is dtau within the z-level
    dtau=10.d0!one
    opac3d=zero
    if(nn_z.gt.zero) then
       opacmax = nn_z*two*dtau*(z(nz)-z(1))/(two*z(nz)*(z(zindx)-z(zindx-1)) - z(zindx)**2 + z(zindx-1)**2)
    elseif(nn_z.lt.zero) then
       !   stop 'to test opacmax for negative n_z'
       opacmax = nn_z*two*dtau*(z(nz)-z(1))/(two*z(nz)*(z(zindx)-z(zindx+1)) - z(zindx)**2 + z(zindx+1)**2)
    endif
    do k=1, nz
       opac3d(:,:,k)=opacmax-opacmax*(z(k)-z(1))/(z(nz)-z(1))
    enddo
    !
    !------------------------------------------------------------------------
    !
    int3d=one
    !
    !
    if(nn_x.gt.zero) then
       startx = 3
       endx = nx-2
       startxb = nx
       endxb = 1
       alpha=  1
    elseif(nn_x.lt.zero) then
       startx = nx-2
       endx = 3
       startxb = 1
       endxb = nx
       alpha=-1
    endif
    !
    if(nn_y.gt.zero) then
       starty = 3
       endy = ny-2
       startyb = ny
       endyb = 1
       beta =  1
    elseif(nn_y.lt.zero) then
       starty = ny-2
       endy = 3
       startyb = 1
       endyb = ny
       beta =-1
    endif
    !
    if(nn_z.gt.zero) then
       startz = 3
       endz = nz-2
       gamma=  1
    elseif(nn_z.lt.zero) then
       startz = nz-2
       endz = 3
       gamma=-1
    endif
    !
    !
    !
    call flc_cont3d(oindx, nueindx, xindx, yindx, zindx, startxb, startyb, endxb, endyb, &
         alpha, beta, gamma, nn_x, nn_y, nn_z)
    !
    if(ver) write(*,*) 'solution from benchmark/theoretical solution', int3d(xindx,yindx,zindx)/exp(-dtau), (int3d(xindx,yindx,zindx)-exp(-dtau))/exp(-dtau)

  end subroutine benchmark02_solution

  
end module mod_benchmark3d
