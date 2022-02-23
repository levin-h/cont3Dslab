module mod_angles
  !
  use prog_type
  use fund_const
  !
  implicit none
  !
  integer(i4b) :: opt_angint_method
  !opt_angint_method=0 if angular integration is used with trapezoidal rule (nodes equidistant in theta and phi)
  !opt_angint_method=1 if angular integration is used with trapezoidal rule (nodes from Lobell & Blomme 2008)
  !opt_angint_method=2 if angular integration is used with simpsons rule (nodes equidistant in theta and phi)
  !                       (note: mu-grid and phi-grid will be made equidistant for three subsequent points)
  !opt_angint_method=3 if angular integration is used with simpson rule corrected for the error
  !                       from a grid with half resolution (also known as boole's rule)
  !opt_angint_method=4 if angular integration is used with cubic splines (catmull-rom-spline, nodes equidistant in theta and phi)
  !opt_angint_method=5 if angular integration is used with gauss-legendre-integration (for each octant)
  !opt_angint_method=6 if angular integration is used with gauss-chebyshev-integration (for each octant)
  !opt_angint_method=7 if angular integration is used with triangulation (linear integrals)
  !opt_angint_method=8 if angular integration is used with triangulation ('pseudo'-gauss integrals per triangle)
  !opt_angint_method=9 if angular integration is used with lebedev interpolation (optimized nodes on the sphere)  
  !
  !data
  integer(i4b) :: ntheta
  integer(i4b) :: dim_mu, dim_phi, nomega
  !
  real(dp), dimension(:), allocatable :: n_x, n_y, n_z, weight_omega
  !
  integer, dimension(:,:), allocatable :: q_alo
  !
  !
  !dim_mu: total number of mu-integration-nodes
  !nodes_mu: mu-integration-nodes
  !weight_mu: mu-integration-weight
  !
  !n_theta: number of integration nodes for theta=[0,pi/2]
  !
  !q_alo: indices for alo-calculations for each direction
  !       to store nearest neighbours correctly


contains
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega(verbose)
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, err    
    !
    if(present(verbose)) ver=verbose    
    !
    if(ver) write(*,*) '---------------------------creating angular grid-------------------------------'
    if(ver) write(*,*)
    !
    select case (opt_angint_method)
      !
    case(0)
      !trapezoidal rule
      call calcnodes_omega0(verbose=ver)
    case(1)
      !trapezoidal rule with direcitonal distribution from lobell&blomme
      stop 'opt_angint_method=1 needs to be tested for slab'
      call calcnodes_omega1(verbose=ver)
    case(2)
      !simpsons rule
      stop 'opt_angint_method=2 needs to be tested for slab'
      call calcnodes_omega2(verbose=ver)
    case(3)
      !boole's rule
      stop 'opt_angint_method=3 needs to be tested for slab'
      call calcnodes_omega3(verbose=ver)
    case(4)
      !hermite spline
      stop 'opt_angint_method=4 needs to be tested for slab'
      call calcnodes_omega4(verbose=ver)
    case(5)
      !gauss-legendre integration
      stop 'opt_angint_method=5 needs to be tested for slab'
      call calcnodes_omega5(verbose=ver)
    case(6)
      !gauss-chebyshev itnegration
      stop 'opt_angint_method=6 needs to be tested for slab'
      call calcnodes_omega6(verbose=ver)
    case(7)
      !triangulation with triangles
      stop 'opt_angint_method=7 needs to be tested for slab'
      call calcnodes_omega7(verbose=ver)
    case(8)
      !triangulation with triangles (pseudo gauss)
      stop 'opt_angint_method=8 needs to be tested for slab'
      call calcnodes_omega8(verbose=ver)
    case(9)
      !lebedev quadrature
      stop 'opt_angint_method=9 needs to be tested for slab'
      call calcnodes_omega9(verbose=ver)
    case default
      stop 'error in calcnodes_omega: opt_angint_method not specified'
      !
    end select
    !
    !
    !set indices for nearest neighbour alo-calculations (direction dependent)
    allocate(q_alo(nomega,27), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega'
    q_alo=0
    do i=1, nomega
      if(n_x(i).gt.zero.and.n_y(i).gt.zero.and.n_z(i).gt.zero) then
        q_alo(i,:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, &
        10, 11, 12, 13, 14, 15, 16, 17, 18, &
        19, 20, 21, 22, 23, 24, 25, 26, 27 /)
      elseif(n_x(i).gt.zero.and.n_y(i).gt.zero.and.n_z(i).lt.zero) then
        q_alo(i,:) = (/ 19, 20, 21, 22, 23, 24, 25, 26, 27, &
        10, 11, 12, 13, 14, 15, 16, 17, 18, &
        1, 2, 3, 4, 5, 6, 7, 8, 9 /)
      elseif(n_x(i).gt.zero.and.n_y(i).lt.zero.and.n_z(i).gt.zero) then
        q_alo(i,:) = (/ 7, 8, 9, 4, 5, 6, 1, 2, 3, &
        16, 17, 18, 13, 14, 15, 10, 11, 12, &
        25, 26, 27, 22, 23, 24, 19, 20, 21 /)
      elseif(n_x(i).lt.zero.and.n_y(i).gt.zero.and.n_z(i).gt.zero) then
        q_alo(i,:) = (/ 3, 2, 1, 6, 5, 4, 9, 8, 7, &
        12, 11, 10, 15, 14, 13, 18, 17, 16, &
        21, 20, 19, 24, 23, 22, 27, 26, 25 /)
      elseif(n_x(i).gt.zero.and.n_y(i).lt.zero.and.n_z(i).lt.zero) then
        q_alo(i,:) = (/ 25, 26, 27, 22, 23, 24, 19, 20, 21, &
        16, 17, 18, 13, 14, 15, 10, 11, 12, &
        7, 8, 9, 4, 5, 6, 1, 2, 3 /)
      elseif(n_x(i).lt.zero.and.n_y(i).gt.zero.and.n_z(i).lt.zero) then
        q_alo(i,:) = (/ 21, 20, 19, 24, 23, 22, 27, 26, 25, &
        12, 11, 10, 15, 14, 13, 18, 17, 16, &
        3, 2, 1, 6, 5, 4, 9, 8, 7 /)
      elseif(n_x(i).lt.zero.and.n_y(i).lt.zero.and.n_z(i).gt.zero) then
        q_alo(i,:) = (/ 9, 8, 7, 6, 5, 4, 3, 2, 1, &
        18, 17, 16, 15, 14, 13, 12, 11, 10, &
        27, 26, 25, 24, 23, 22, 21, 20, 19 /)
      elseif(n_x(i).lt.zero.and.n_y(i).lt.zero.and.n_z(i).lt.zero) then
        q_alo(i,:) = (/ 27, 26, 25, 24, 23, 22, 21, 20, 19, &
        18, 17, 16, 15, 14, 13, 12, 11, 10, &
        9, 8, 7, 6, 5, 4, 3, 2, 1 /)
      endif
    enddo
    !
    if(ver) write(*,*) 'n_omega', nomega
    if(ver) write(*,*)
    !
  end subroutine calcnodes_omega
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega0(verbose)
    !
    !-----------------------------------------------------------------------
    !   equidistant theta-grid
    !   equidistant phi-grid, with delta_phi=delta-theta
    !   including exact opposite directions
    !   trapezoidal rule
    !-----------------------------------------------------------------------
    !
    use mod_integ1d, only: precalc_oweight_trapez
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: dim_mu, dim_phi
    real(dp) :: theta_min, theta_max
    real(dp) :: phi_min, phi_max
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_mu, weight_mu, nodes_phi, weight_phi, fdum_arr
    !
    if(present(verbose)) ver=verbose
    !
    dim_mu=2*ntheta
    dim_phi=2*ntheta
    nomega=2*dim_mu*dim_phi
    !
    !----------------------------mu grid------------------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    allocate(weight_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
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
    do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.small_number) nodes_mu(i)=zero
    enddo
    !
    call precalc_oweight_trapez(nodes_mu, dim_mu, -one, one, weight_mu)
    !
    !normalization
    weight_mu=weight_mu/two
    !
    !----------------------------phi grid-----------------------------------
    !
    allocate(nodes_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    allocate(weight_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    !
    if(dim_phi.eq.2) then
      !to be consistent with diffusion equation (two-stream approximation)
      phi_min=pi/four
      phi_max=three*pi/four
    else
      phi_min=1.d-6
      phi_max=pi-1.d-6
    endif
    !
    do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
    enddo
    !
    call precalc_oweight_trapez(nodes_phi, dim_phi, zero, pi, weight_phi)
    !
    !write(*,*) nodes_phi*180./pi
    !stop
    !normalization
    weight_phi=weight_phi/pi
    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        n_x(k)=sqrt(one-nodes_mu(i)**2)*cos(nodes_phi(j))
        n_y(k)=sqrt(one-nodes_mu(i)**2)*sin(nodes_phi(j))
        n_z(k)=nodes_mu(i)
        if(abs(n_x(k)).lt.1.d-14) n_x(k)=zero
        if(abs(n_y(k)).lt.1.d-14) n_y(k)=zero
        if(abs(n_z(k)).lt.1.d-14) n_z(k)=zero
        !
        !opposite direction
        n_x(k+1)=-n_x(k)
        n_y(k+1)=-n_y(k)
        n_z(k+1)=-n_z(k)
        if(abs(n_x(k+1)).lt.1.d-14) n_x(k)=zero
        if(abs(n_y(k+1)).lt.1.d-14) n_y(k)=zero
        if(abs(n_z(k+1)).lt.1.d-14) n_z(k)=zero
        !
        weight_omega(k) = weight_mu(i)*weight_phi(j)/two
        weight_omega(k+1) = weight_mu(i)*weight_phi(j)/two
        k=k+2
      enddo
    enddo
    !
    !allocate(fdum_arr(nomega))
    !fdum_arr=n_y
    !n_y=n_x
    !n_x=fdum_arr
    !write(*,*) 'NOTE*************: UPDATE CALCNODES_OMEGA0 FOR FDUM_ARR'
    !
    write(*,*) n_x
    write(*,*) n_y
    write(*,*) n_z
    stop 'go on in calcnodes_omega0'
    !
  end subroutine calcnodes_omega0
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega1(verbose)
    !
    !-----------------------------------------------------------------------
    !   equidistant theta-grid
    !   equidistant phi-grid, calculated for each theta-level
    !      such that domega=const (see Lobell&Blomme 2008)
    !   trapezoidal rule
    !-----------------------------------------------------------------------
    !
    use mod_integ1d, only: precalc_oweight_trapez
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err, n_phi
    integer(i4b) :: dim_mu, dim_phi
    real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
    !                                                             and slightly asymmetric, to avoid mu=0.

    real(dp), parameter :: phi_min=1.d-6, phi_max=two*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
    !                                                             !and slightly asymmetric therefore as well

    real(dp) :: del_theta, del_phi, del_omega
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_mu, weight_mu, nodes_phi, weight_phi, nodes_theta
    integer(i4b), dimension(:), allocatable :: dim_phi_arr
    real(dp), dimension(:,:), allocatable :: nodes_phi_pair, weight_phi_pair
    logical, dimension(:,:), allocatable :: phi_mask
    !
    ! ... local functions
    !
    if(present(verbose)) ver=verbose
    !
    dim_mu=2*ntheta-1
    dim_phi=2*dim_mu-1
    !nomega=dim_mu*dim_phi
    !
    !----------------------------mu grid------------------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    allocate(weight_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    !
    do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
        nodes_mu(i)=zero
      endif
    enddo
    !
    call precalc_oweight_trapez(nodes_mu, dim_mu, -one, one, weight_mu)
    !
    !normalization
    weight_mu=weight_mu/two
    !
    !----------------calculate theta-array from (given) nodes_mu-------------
    !
    allocate(nodes_theta(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error calcnodes_omega1'
    !
    !also reverse array
    do i=1, dim_mu
      nodes_theta(i) = acos(nodes_mu(dim_mu+1-i))
    enddo
    !
    !-------------------calculate equidistant phi array---------------------
    !----------------with different n_phi on each theta-level---------------
    !
    allocate(nodes_phi_pair(dim_mu,dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    allocate(weight_phi_pair(dim_mu,dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    allocate(phi_mask(dim_mu,dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    allocate(dim_phi_arr(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    !
    phi_mask=.false.
    nodes_phi_pair=zero
    weight_phi_pair=zero
    !
    !del_theta, del_phi and del_omega in equatorial plane
    del_theta=nodes_theta(ntheta)-nodes_theta(ntheta-1)
    del_phi=two*pi/(dim_phi-1)
    del_omega=del_theta*del_phi
    !
    nomega=0
    !
    do i=1, dim_mu
      !
      !-----------------calculate phi-grid for given theta--------------------
      !-------new version: del_theta is allowed to be latitude dependent------
      !
      !calculate del_theta (of theta)
      if(i.eq.1) then
        del_theta=nodes_theta(i+1)-nodes_theta(i)
      else if(i.eq.dim_mu) then
        del_theta=nodes_theta(i)-nodes_theta(i-1)
      else
        !use central differences in order to provide symmetry
        del_theta=(nodes_theta(i+1)-nodes_theta(i-1))/two
      endif
      !
      !n_phi minimum to be 2 and has to be odd
      n_phi=n_phi_lobl(del_omega, del_theta, nodes_theta(i))
      n_phi=maxval((/n_phi,2/))
      n_phi=n_phi/2
      n_phi=2*n_phi+1
      !
      if(n_phi.gt.dim_phi) stop 'error calcnodes_omega1: n_phi gt dim_phi'
      !
      dim_phi_arr(i)=n_phi
      !
      !calculate phi grid
      do j=1, n_phi
        nodes_phi_pair(i,j)=phi_min + float(j-1)*(phi_max-phi_min)/(n_phi-1)
        phi_mask(i,j) = .true.
        nomega=nomega+1
      enddo
      !
      !calculate weights
      allocate(nodes_phi(n_phi), stat=err)
      if(err.ne.0) stop 'allocation error in calcnodes_omega1'
      allocate(weight_phi(n_phi), stat=err)
      if(err.ne.0) stop 'allocation error in calcnodes_omega1'
      nodes_phi=nodes_phi_pair(i,1:n_phi)
      call precalc_oweight_trapez(nodes_phi, n_phi, zero, two*pi, weight_phi)
      weight_phi_pair(i,1:n_phi)=weight_phi
      deallocate(nodes_phi)
      deallocate(weight_phi)
      !
    enddo
    !
    weight_phi_pair=weight_phi_pair/two/pi
    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega1'
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        if(phi_mask(i,j)) then
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi_pair(i,j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi_pair(i,j))
          n_z(k)=nodes_mu(i)

          !         write(*,*) nodes_phi_pair(i,j), n_x(k)
          !
          if(abs(n_x(k)).lt.1.d-14) n_x(k)=zero
          if(abs(n_y(k)).lt.1.d-14) n_y(k)=zero
          if(abs(n_z(k)).lt.1.d-14) n_z(k)=zero
          !
          weight_omega(k) = weight_mu(i)*weight_phi_pair(i,j)
          k=k+1
        endif
      enddo
    enddo
    !
    !
    !
  end subroutine calcnodes_omega1
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega2(verbose)
    !
    !-----------------------------------------------------------------------
    !   equidistant theta-grid
    !   equidistant phi-grid, with delta_phi=delta-theta
    !   simpsons rule
    !-----------------------------------------------------------------------
    !
    use mod_integ1d, only: precalc_oweight_simps
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: dim_mu, dim_phi
    real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
    !                                                             and slightly asymmetric, to avoid mu=0.
    real(dp), parameter :: phi_min=1.d-6, phi_max=two*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
    !                                                             !and slightly asymmetric therefore as well
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_phi, weight_phi, nodes_mu, weight_mu

    if(present(verbose)) ver=verbose
    
    !
    dim_mu=2*ntheta-1
    dim_phi=2*dim_mu-1
    nomega=dim_mu*dim_phi
    !
    !----------------------------mu grid------------------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega2'
    allocate(weight_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega2'
    !
    do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
        nodes_mu(i)=zero
      endif
    enddo
    !
    call precalc_oweight_simps(nodes_mu, dim_mu, -one, one, weight_mu)
    !
    !normalization
    weight_mu=weight_mu/two
    !
    !----------------------------phi grid-----------------------------------
    !
    allocate(nodes_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega2'
    allocate(weight_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega2'
    !
    do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
    enddo
    !
    call precalc_oweight_simps(nodes_phi, dim_phi, zero, two*pi, weight_phi)
    !
    !normalization
    weight_phi=weight_phi/two/pi
    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega2'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega2'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega2'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega2'
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
        n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
        n_z(k)=nodes_mu(i)
        !
        if(abs(n_x(k)).lt.1.d-14) n_x(k)=zero
        if(abs(n_y(k)).lt.1.d-14) n_y(k)=zero
        if(abs(n_z(k)).lt.1.d-14) n_z(k)=zero
        !
        weight_omega(k) = weight_mu(i)*weight_phi(j)
        k=k+1
      enddo
    enddo
    !
    !
    !
  end subroutine calcnodes_omega2
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega3(verbose)
    !
    !-----------------------------------------------------------------------
    !   equidistant theta-grid
    !   equidistant phi-grid, with delta_phi=delta-theta
    !   boole's rule
    !-----------------------------------------------------------------------
    !
    use mod_integ1d, only: precalc_oweight_boole
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: dim_mu, dim_phi
    real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
    !                                                             and slightly asymmetric, to avoid mu=0.
    real(dp), parameter :: phi_min=1.d-6, phi_max=two*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
    !                                                             !and slightly asymmetric therefore as well
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_phi, weight_phi, nodes_mu, weight_mu

    if(present(verbose)) ver=verbose
    
    !
    dim_mu=2*ntheta-1
    dim_phi=2*dim_mu-1
    nomega=dim_mu*dim_phi
    !
    !----------------------------mu grid------------------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega3'
    allocate(weight_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega3'
    !
    do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
        nodes_mu(i)=zero
      endif
    enddo
    !
    call precalc_oweight_boole(nodes_mu, dim_mu, -one, one, weight_mu)
    !
    !normalization
    weight_mu=weight_mu/two
    !
    !----------------------------phi grid-----------------------------------
    !
    allocate(nodes_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega3'
    allocate(weight_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega3'
    !
    do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
    enddo
    !
    call precalc_oweight_boole(nodes_phi, dim_phi, zero, two*pi, weight_phi)
    !
    !normalization
    weight_phi=weight_phi/two/pi
    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega3'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega3'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega3'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega3'
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
        n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
        n_z(k)=nodes_mu(i)
        !
        if(abs(n_x(k)).lt.1.d-14) n_x(k)=zero
        if(abs(n_y(k)).lt.1.d-14) n_y(k)=zero
        if(abs(n_z(k)).lt.1.d-14) n_z(k)=zero
        !
        weight_omega(k) = weight_mu(i)*weight_phi(j)
        k=k+1
      enddo
    enddo
    !
    !
    !
  end subroutine calcnodes_omega3
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega4(verbose)
    !
    !-----------------------------------------------------------------------
    !   equidistant theta-grid
    !   equidistant phi-grid, with delta_phi=delta-theta
    !   hermite spline integration
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: dim_mu, dim_phi
    real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
    !                                                             and slightly asymmetric, to avoid mu=0.
    real(dp), parameter :: phi_min=1.d-6, phi_max=two*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
    !                                                             !and slightly asymmetric therefore as well
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_phi, weight_phi, nodes_mu, weight_mu
    !
    if(present(verbose)) ver=verbose        
    !
    dim_mu=2*ntheta-1
    dim_phi=2*dim_mu-1
    nomega=dim_mu*dim_phi
    !
    !----------------------------mu grid------------------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega4'
    allocate(weight_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega4'
    !
    do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
      nodes_mu(i) = cos(nodes_mu(i))
      if(abs(nodes_mu(i)).lt.1.d-15) then
        nodes_mu(i)=zero
      endif
    enddo
    !
    stop 'precalc_weight_spline needs to be rewritten for open boundaries'
    call precalc_weight_spline(nodes_mu, dim_mu, weight_mu, .false.)
    !
    !normalization
    weight_mu=weight_mu/two
    !
    !----------------------------phi grid-----------------------------------
    !
    allocate(nodes_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega4'
    allocate(weight_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega4'
    !
    do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
    enddo
    !
    stop 'precalc_weight_spline needs to be rewritten for open boundaries'
    call precalc_weight_spline(nodes_phi, dim_phi, weight_phi, .true.)
    !
    !normalization
    weight_phi=weight_phi/two/pi
    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega4'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega4'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega4'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega4'
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
        n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
        n_z(k)=nodes_mu(i)
        !
        if(abs(n_x(k)).lt.1.d-14) n_x(k)=zero
        if(abs(n_y(k)).lt.1.d-14) n_y(k)=zero
        if(abs(n_z(k)).lt.1.d-14) n_z(k)=zero
        !
        weight_omega(k) = weight_mu(i)*weight_phi(j)
        k=k+1
      enddo
    enddo
    !
    !
    !
  end subroutine calcnodes_omega4
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega5(verbose)
    !
    !-----------------------------------------------------------------------
    !   gauss legendre integration with ntheta points in [zero,pi/2.d0]
    !      (i.e., same integration in each octant)
    !-----------------------------------------------------------------------
    !
    use mod_integ1d, only: precalc_weight_legendre
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: dim_mu, dim_phi
    real(dp) :: theta_min, theta_max, phi_min, phi_max, sum
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_dum, weight_dum, nodes_phi, weight_phi, &
    nodes_rev, weight_rev, nodes_mu, weight_mu
    !

    if(present(verbose)) ver=verbose
    
    dim_mu=2*ntheta
    dim_phi=4*ntheta
    nomega=dim_mu*dim_phi
    !
    !---------------------standard grid in [-1,d0,1.d0]---------------------
    !
    allocate(nodes_dum(ntheta), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    allocate(weight_dum(ntheta), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    !
    call precalc_weight_legendre(ntheta, -one, one, nodes_dum, weight_dum)
    !
    !-----------------------------mu grid-----------------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    allocate(weight_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    !
    theta_min=zero
    theta_max=pi/two
    nodes_mu(1:ntheta)=(theta_max-theta_min)/two * nodes_dum + (theta_min+theta_max)/two
    weight_mu(1:ntheta)=weight_dum*(theta_max-theta_min)/two*sin(nodes_mu(1:ntheta))
    !
    theta_min=pi/two
    theta_max=pi
    nodes_mu(ntheta+1:dim_mu) = (theta_max-theta_min)/two * nodes_dum + (theta_min+theta_max)/two
    weight_mu(ntheta+1:dim_mu)=weight_dum*sin(nodes_mu(ntheta+1:dim_mu))*(theta_max-theta_min)/two
    !
    nodes_mu=cos(nodes_mu)
    !
    weight_mu = weight_mu/two
    !
    allocate(nodes_rev(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error calcnodes_omega5'
    allocate(weight_rev(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error calcnodes_omega5'
    !
    do i=1, dim_mu
      nodes_rev(i) = nodes_mu(dim_mu+1-i)
      weight_rev(i) = weight_mu(dim_mu+1-i)
    enddo
    !
    nodes_mu=nodes_rev
    weight_mu=weight_rev
    !
    !-------------------------------phi grid--------------------------------
    !
    allocate(nodes_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    allocate(weight_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    !
    phi_min=zero
    phi_max=pi/two
    nodes_phi(1:ntheta) = (phi_max-phi_min)/two * nodes_dum + (phi_min+phi_max)/two
    weight_phi(1:ntheta) = weight_dum*(phi_max-phi_min)/two
    !
    phi_min=pi/two
    phi_max=pi
    nodes_phi(ntheta+1:2*ntheta) = (phi_max-phi_min)/two * nodes_dum + (phi_min+phi_max)/two
    weight_phi(ntheta+1:2*ntheta) = weight_dum*(phi_max-phi_min)/two
    !
    phi_min=pi
    phi_max=three*pi/two
    nodes_phi(2*ntheta+1:3*ntheta) = (phi_max-phi_min)/two * nodes_dum + (phi_min+phi_max)/two
    weight_phi(2*ntheta+1:3*ntheta) = weight_dum*(phi_max-phi_min)/two
    !
    phi_min=three*pi/two
    phi_max=two*pi
    nodes_phi(3*ntheta+1:4*ntheta) = (phi_max-phi_min)/two * nodes_dum + (phi_min+phi_max)/two
    weight_phi(3*ntheta+1:4*ntheta) = weight_dum*(phi_max-phi_min)/two
    !
    weight_phi=weight_phi/two/pi
    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega5'
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
        n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
        n_z(k)=nodes_mu(i)
        !
        if(abs(n_x(k)).lt.1.d-14) n_x(k)=zero
        if(abs(n_y(k)).lt.1.d-14) n_y(k)=zero
        if(abs(n_z(k)).lt.1.d-14) n_z(k)=zero
        !
        weight_omega(k) = weight_mu(i)*weight_phi(j)
        k=k+1
      enddo
    enddo
    !
    !
    !
  end subroutine calcnodes_omega5
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega6(verbose)
    !
    !-----------------------------------------------------------------------
    !   gauss chebyshev integration with ntheta points in [0.d0,pi/2.d0]
    !      (i.e., same integration in each octant)
    !
    ! NOTE: EVENTUALLY TO BE DEBUGGED BECAUSE WEIGHTS DO NOT SUM TO 1!!!!
    !-----------------------------------------------------------------------
    !
    use mod_integ1d, only: precalc_weight_chebyshev
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: dim_mu, dim_phi
    real(dp) :: theta_min, theta_max, phi_min, phi_max, sum
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_dum, weight_dum, nodes_phi, weight_phi, &
    nodes_rev, weight_rev, nodes_mu, weight_mu
    !
    if(present(verbose)) ver=verbose
    
    dim_mu=2*ntheta
    dim_phi=4*ntheta
    nomega=dim_mu*dim_phi
    !
    !---------------------standard grid in [-1,d0,1.d0]---------------------
    !
    allocate(nodes_dum(ntheta), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    allocate(weight_dum(ntheta), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    !
    call precalc_weight_chebyshev(ntheta, -one, one, nodes_dum, weight_dum)
    !
    !-----------------------------mu grid-----------------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    allocate(weight_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    !
    theta_min=zero
    theta_max=pi/two
    nodes_mu(1:ntheta)=(theta_max-theta_min)/two * nodes_dum + (theta_min+theta_max)/two
    weight_mu(1:ntheta)=weight_dum*(theta_max-theta_min)/two*sin(nodes_mu(1:ntheta))
    !
    theta_min=pi/two
    theta_max=pi
    nodes_mu(ntheta+1:dim_mu) = (theta_max-theta_min)/two * nodes_dum + (theta_min+theta_max)/two
    weight_mu(ntheta+1:dim_mu)=weight_dum*sin(nodes_mu(ntheta+1:dim_mu))*(theta_max-theta_min)/two
    !
    nodes_mu=cos(nodes_mu)
    !
    weight_mu = weight_mu/two
    !
    allocate(nodes_rev(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error calcnodes_omega6'
    allocate(weight_rev(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error calcnodes_omega6'
    !
    do i=1, dim_mu
      nodes_rev(i) = nodes_mu(dim_mu+1-i)
      weight_rev(i) = weight_mu(dim_mu+1-i)
    enddo
    !
    nodes_mu=nodes_rev
    weight_mu=weight_rev
    !
    !-------------------------------phi grid--------------------------------
    !
    allocate(nodes_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    allocate(weight_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    !
    phi_min=zero
    phi_max=pi/two
    nodes_phi(1:ntheta) = (phi_max-phi_min)/two * nodes_dum + (phi_min+phi_max)/two
    weight_phi(1:ntheta) = weight_dum*(phi_max-phi_min)/two
    !
    phi_min=pi/two
    phi_max=pi
    nodes_phi(ntheta+1:2*ntheta) = (phi_max-phi_min)/two * nodes_dum + (phi_min+phi_max)/two
    weight_phi(ntheta+1:2*ntheta) = weight_dum*(phi_max-phi_min)/two
    !
    phi_min=pi
    phi_max=three*pi/two
    nodes_phi(2*ntheta+1:3*ntheta) = (phi_max-phi_min)/two * nodes_dum + (phi_min+phi_max)/two
    weight_phi(2*ntheta+1:3*ntheta) = weight_dum*(phi_max-phi_min)/two
    !
    phi_min=three*pi/two
    phi_max=two*pi
    nodes_phi(3*ntheta+1:4*ntheta) = (phi_max-phi_min)/two * nodes_dum + (phi_min+phi_max)/two
    weight_phi(3*ntheta+1:4*ntheta) = weight_dum*(phi_max-phi_min)/two
    !
    weight_phi=weight_phi/two/pi
    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega6'
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
        n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
        n_z(k)=nodes_mu(i)
        !
        if(abs(n_x(k)).lt.1.d-14) n_x(k)=zero
        if(abs(n_y(k)).lt.1.d-14) n_y(k)=zero
        if(abs(n_z(k)).lt.1.d-14) n_z(k)=zero
        !
        weight_omega(k) = weight_mu(i)*weight_phi(j)
        k=k+1
      enddo
    enddo
    !
    !
  end subroutine calcnodes_omega6
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega7(verbose)
    !
    !
    !-----------------------------------------------------------------------
    !   split 2d integration domain in triangles
    !   integration with barycentric interpolation in each triangle
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: dim_mu, dim_phi
    real(dp) :: del_im, del_ip, phi1, phi2, mu1, mu2, w1, w2, w3, sum
    real(dp), parameter :: phi_min=1.d-6, phi_max=two*pi-1.d-5   !open interval, n_x=0, n_y=0, n_z=0 not allowed
    !                                                             !and slightly asymmetric therefore as well
    real(dp), parameter :: theta_min=1.d-6, theta_max=pi-1.d-5   !open interval, since n_x=0 or n_z=0 not allowed,
    !                                                             and slightly asymmetric, to avoid mu=0.
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_phi, nodes_mu
    !
    if(present(verbose)) ver=verbose        
    !
    dim_mu=2*ntheta-1
    dim_phi=2*dim_mu-1
    !
    if(dim_mu.lt.3) stop 'error in calcnodes_omega7: dim_mu has to be ge 3'
    !if(mod(dim_mu,2).eq.0) stop 'error in calcnodes_omega7: dim_mu has to be odd'
    if(dim_phi.lt.3) stop 'error in calcnodes_omega7:: dim_phi has to be ge 3'
    !if(mod(dim_phi,2).eq.0) stop 'error in calcnodes_omega7: dim_phi has to be odd'
    !
    !----------------------------mu grid------------------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega7'
    !
    do i=1, dim_mu
      nodes_mu(i) = theta_max - (i-1)*(theta_max-theta_min)/(dim_mu-1)
    enddo
    nodes_mu=cos(nodes_mu)
    !
    !make mu-grid equidistant for three subsequent nodes
    do i=2, dim_mu, 2
      nodes_mu(i) = (nodes_mu(i+1)+nodes_mu(i-1))/two
    enddo
    !
    !check if three subsequent x-points are really equidistant
    do i=2, dim_mu-1, 2
      del_im=nodes_mu(i)-nodes_mu(i-1)
      del_ip=nodes_mu(i+1)-nodes_mu(i)
      if(abs(del_im-del_ip).gt.1.d-15) then
        stop 'error in calcnodes_omega7: mu-grid not equidistant'
      endif
    enddo
    !
    !----------------------------phi grid-----------------------------------
    !
    allocate(nodes_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega7'
    !
    do i=1, dim_phi
      nodes_phi(i) = phi_min + (i-1)*(phi_max-phi_min)/(dim_phi-1)
    enddo
    !
    nomega=ceiling(dim_mu/2.)*ceiling(dim_phi/2.) + floor(dim_mu/2.)*floor(dim_phi/2.)
    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega7'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega7'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega7'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega7'
    !
    weight_omega=zero
    !
    phi1=zero
    phi2=two*pi
    mu1=-one
    mu2=one
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        if(i.eq.1.and.j.eq.1) then   !weight at (1,1)
          w1=((nodes_mu(i+1)-nodes_mu(i))*(nodes_phi(j+2)-nodes_phi(j)) + &
          (nodes_mu(i+2)-nodes_mu(i))*(nodes_phi(j+1)-nodes_phi(j)))/six
          w2=((nodes_mu(i+2)-nodes_mu(i))*(nodes_phi(j)-phi1) + &
          (nodes_mu(i)-mu1)*(nodes_phi(j+2)-nodes_phi(j)))/two
          w3=(nodes_mu(i)-mu1)*(nodes_phi(j)-phi1)
          weight_omega(k)=w1+w2+w3
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(i.eq.dim_mu.and.j.eq.1) then   !weight at (dim_mu,1)
          w1=((nodes_mu(i)-nodes_mu(i-2))*(nodes_phi(j+1)-nodes_phi(j)) + &
          (nodes_mu(i)-nodes_mu(i-1))*(nodes_phi(j+2)-nodes_phi(j)))/six
          w2=((nodes_mu(i)-nodes_mu(i-2))*(nodes_phi(j)-phi1) + &
          (mu2-nodes_mu(i))*(nodes_phi(j+2)-nodes_phi(j)))/two
          w3=(mu2-nodes_mu(i))*(nodes_phi(j)-phi1)
          weight_omega(k)=w1+w2+w3
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(i.eq.1.and.j.eq.dim_phi) then   !weight at (1,dim_phi)
          w1=((nodes_mu(i+1)-nodes_mu(i))*(nodes_phi(j)-nodes_phi(j-2)) + &
          (nodes_mu(i+2)-nodes_mu(i))*(nodes_phi(j)-nodes_phi(j-1)))/six
          w2=((nodes_mu(i+2)-nodes_mu(i))*(phi2-nodes_phi(j)) + &
          (nodes_mu(i)-mu1)*(nodes_phi(j)-nodes_phi(j-2)))/two
          w3=(nodes_mu(i)-mu1)*(phi2-nodes_phi(j))
          weight_omega(k)=w1+w2+w3
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(i.eq.dim_mu.and.j.eq.dim_phi) then   !weight at (dim_mu,dim_phi)
          w1=((nodes_mu(i)-nodes_mu(i-2))*(nodes_phi(j)-nodes_phi(j-1)) + &
          (nodes_mu(i)-nodes_mu(i-1))*(nodes_phi(j)-nodes_phi(j-2)))/six
          w2=((nodes_mu(i)-nodes_mu(i-2))*(phi2-nodes_phi(j)) + &
          (mu2-nodes_mu(i))*(nodes_phi(j)-nodes_phi(j-2)))/two
          w3=(mu2-nodes_mu(i))*(phi2-nodes_phi(j))
          weight_omega(k)=w1+w2+w3
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(j.eq.1.and.mod(i,2).ne.0) then   !weights along edge at j=1 (only odd indices contribute)
          w1=((nodes_mu(i+2)-nodes_mu(i-2))*(nodes_phi(j+1)-nodes_phi(j)) + &
          (nodes_mu(i+1)-nodes_mu(i-1))*(nodes_phi(j+2)-nodes_phi(j)))/six
          w2=(nodes_mu(i+2)-nodes_mu(i-2))*(nodes_phi(j)-phi1)/two
          weight_omega(k)=w1+w2
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(j.eq.dim_phi.and.mod(i,2).ne.0) then   !weights along edge at j=dim_phi (only odd indices contribute)
          w1=((nodes_mu(i+2)-nodes_mu(i-2))*(nodes_phi(j)-nodes_phi(j-1)) + &
          (nodes_mu(i+1)-nodes_mu(i-1))*(nodes_phi(j)-nodes_phi(j-2)))/six
          w2=(nodes_mu(i+2)-nodes_mu(i-2))*(phi2-nodes_phi(j))/two
          weight_omega(k)=w1+w2
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(i.eq.1.and.mod(j,2).ne.0) then   !weights along edge at i=1 (only odd indices contribute)
          w1=((nodes_mu(i+1)-nodes_mu(i))*(nodes_phi(j+2)-nodes_phi(j-2)) + &
          (nodes_mu(i+2)-nodes_mu(i))*(nodes_phi(j+1)-nodes_phi(j-1)))/six
          w2=(nodes_mu(i)-mu1)*(nodes_phi(j+2)-nodes_phi(j-2))/two
          weight_omega(k)=w1+w2
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(i.eq.dim_mu.and.mod(j,2).ne.0) then   !weights along edge at i=dim_mu (only odd indices contribute)
          w1=((nodes_mu(i)-nodes_mu(i-1))*(nodes_phi(j+2)-nodes_phi(j-2)) + &
          (nodes_mu(i)-nodes_mu(i-2))*(nodes_phi(j+1)-nodes_phi(j-1)))/six
          w2=(mu2-nodes_mu(i))*(nodes_phi(j+2)-nodes_phi(j-2))/two
          weight_omega(k)=w1+w2
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(mod(i,2).eq.0.and.mod(j,2).eq.0) then   !weights at each rectangle center (only even indices)
          weight_omega(k)=(nodes_mu(i+1)-nodes_mu(i-1))*(nodes_phi(j+1)-nodes_phi(j-1))/three
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        elseif(mod(i,2).ne.0.and.mod(j,2).ne.0) then   !weights at each rectangle corner (only odd indices)
          w1=(nodes_mu(i+2)-nodes_mu(i-2))*(nodes_phi(j+1)-nodes_phi(j-1))
          w2=(nodes_mu(i+1)-nodes_mu(i-1))*(nodes_phi(j+2)-nodes_phi(j-2))
          weight_omega(k)=(w1+w2)/six
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          k=k+1
        endif
      enddo
    enddo
    !
    !normalize weights
    weight_omega=weight_omega/four/pi
    !
    !
  end subroutine calcnodes_omega7
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega8(verbose)
    !
    !-----------------------------------------------------------------------
    !   split 2d integration domain in triangles
    !   integration with 7 point 'pseudo gauss' in each triangle
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, kx, ky, err
    integer(i4b) :: dim_mu, dim_phi
    integer(i4b) :: dim_mu2, dim_phi2
    real(dp) :: jaca, jacb, jacc, jacd, w1, w2
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_mu, nodes_mu2, nodes_phi2, nodes_phi
    real(dp), dimension(:,:), allocatable :: weight_omega2
    logical, dimension(:,:), allocatable :: phi_mask
    !
    if(present(verbose)) ver=verbose
    !
    dim_mu=7*ntheta
    dim_phi=2*dim_mu
    !
    if(dim_mu.lt.7) stop 'error in calcnodes_omega8: dim_mu has to be ge 7'
    if(dim_phi.lt.7) stop 'error in calcnodes_omega8: dim_phi has to be ge 7'
    !
    dim_mu2=dim_mu/7 + 1
    dim_phi2=dim_phi/7 + 1
    !
    !actual grid, for which gauss points are distributed
    allocate(nodes_mu2(dim_mu2), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega8'
    allocate(nodes_phi2(dim_phi2), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega8'
    !
    !-------------------------equidistant theta grid-------------------------
    !
    do i=1, dim_mu2
      nodes_mu2(i) = pi - (i-1)*pi/(dim_mu2-1)
    enddo
    nodes_mu2=cos(nodes_mu2)
    !
    !-------------------------equidistant phi grid--------------------------
    !
    do i=1, dim_phi2
      nodes_phi2(i) = (i-1)*two*pi/(dim_phi2-1)
    enddo
    !
    !-----------------------distribute gauss points-------------------------
    !
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega8'
    allocate(nodes_phi(dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega8'
    !
    nodes_mu=zero
    k=1
    do i=2, dim_mu2
      nodes_mu(k)   = (9*nodes_mu2(i-1)+nodes_mu2(i))/ten
      nodes_mu(k+1) = (5*nodes_mu2(i-1)+nodes_mu2(i))/six
      nodes_mu(k+2) = (7*nodes_mu2(i-1)+3*nodes_mu2(i))/ten
      nodes_mu(k+3) = (nodes_mu2(i-1)+nodes_mu2(i))/two
      nodes_mu(k+4) = (7*nodes_mu2(i)+3*nodes_mu2(i-1))/ten
      nodes_mu(k+5) = (5*nodes_mu2(i)+nodes_mu2(i-1))/six
      nodes_mu(k+6) = (9*nodes_mu2(i)+nodes_mu2(i-1))/ten
      k=k+7
    enddo
    !
    k=1
    do i=2, dim_phi2
      nodes_phi(k)   = (9*nodes_phi2(i-1)+nodes_phi2(i))/ten
      nodes_phi(k+1) = (5*nodes_phi2(i-1)+nodes_phi2(i))/six
      nodes_phi(k+2) = (7*nodes_phi2(i-1)+3*nodes_phi2(i))/ten
      nodes_phi(k+3) = (nodes_phi2(i-1)+nodes_phi2(i))/two
      nodes_phi(k+4) = (7*nodes_phi2(i)+3*nodes_phi2(i-1))/ten
      nodes_phi(k+5) = (5*nodes_phi2(i)+nodes_phi2(i-1))/six
      nodes_phi(k+6) = (9*nodes_phi2(i)+nodes_phi2(i-1))/ten
      k=k+7
    enddo
    !
    !--------------------calculate weights and mask-------------------------
    !
    allocate(phi_mask(dim_mu,dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega8'
    allocate(weight_omega2(dim_mu,dim_phi), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega8'
    !
    phi_mask=.false.
    weight_omega2=zero
    !
    !standard weights (on unit triangle)
    w1=twentyfive/fortyeight/two
    w2=-twentyseven/fortyeight/two
    !
    nomega=0
    !
    ky=1
    do j=2, dim_phi2
      kx=1
      do i=2, dim_mu2
        !jacobians for each triangle
        jaca=(nodes_mu2(i)-nodes_mu2(i-1))*(nodes_phi2(j)-nodes_phi2(j-1))/two
        jacb=(nodes_mu2(i)-nodes_mu2(i-1))*(nodes_phi2(j)-nodes_phi2(j-1))/two
        jacc=(nodes_mu2(i)-nodes_mu2(i-1))*(nodes_phi2(j)-nodes_phi2(j-1))/two
        jacd=(nodes_mu2(i)-nodes_mu2(i-1))*(nodes_phi2(j)-nodes_phi2(j-1))/two
        !triangle a
        weight_omega2(kx+2,ky)=w1*jaca
        weight_omega2(kx+4,ky)=w1*jaca
        weight_omega2(kx+3,ky+1)=w2*jaca
        weight_omega2(kx+3,ky+2)=w1*jaca
        phi_mask(kx+2,ky)=.true.
        phi_mask(kx+4,ky)=.true.
        phi_mask(kx+3,ky+1)=.true.
        phi_mask(kx+3,ky+2)=.true.
        !triangle b
        weight_omega2(kx+6,ky+2)=w1*jacb
        weight_omega2(kx+4,ky+3)=w1*jacb
        weight_omega2(kx+5,ky+3)=w2*jacb
        weight_omega2(kx+6,ky+4)=w1*jacb
        phi_mask(kx+6,ky+2)=.true.
        phi_mask(kx+4,ky+3)=.true.
        phi_mask(kx+5,ky+3)=.true.
        phi_mask(kx+6,ky+4)=.true.
        !triangle c
        weight_omega2(kx+3,ky+4)=w1*jacc
        weight_omega2(kx+3,ky+5)=w2*jacc
        weight_omega2(kx+2,ky+6)=w1*jacc
        weight_omega2(kx+4,ky+6)=w1*jacc
        phi_mask(kx+3,ky+4)=.true.
        phi_mask(kx+3,ky+5)=.true.
        phi_mask(kx+2,ky+6)=.true.
        phi_mask(kx+4,ky+6)=.true.
        !triangle c
        weight_omega2(kx,ky+2)=w1*jacd
        weight_omega2(kx+1,ky+3)=w2*jacd
        weight_omega2(kx+2,ky+3)=w1*jacd
        weight_omega2(kx,ky+4)=w1*jacd
        phi_mask(kx,ky+2)=.true.
        phi_mask(kx+1,ky+3)=.true.
        phi_mask(kx+2,ky+3)=.true.
        phi_mask(kx,ky+4)=.true.
        nomega=nomega+16
        kx=kx+7
      enddo
      ky=ky+7
    enddo

    !
    !------------------------all directions and weights---------------------
    !
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega0'
    !
    k=1
    do i=1, dim_mu
      do j=1, dim_phi
        if(phi_mask(i,j)) then
          n_x(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*cos(nodes_phi(j))
          n_y(k)=sqrt(one-nodes_mu(i)*nodes_mu(i))*sin(nodes_phi(j))
          n_z(k)=nodes_mu(i)
          !
          if(abs(n_x(k)).lt.1.d-14) n_x(k)=zero
          if(abs(n_y(k)).lt.1.d-14) n_y(k)=zero
          if(abs(n_z(k)).lt.1.d-14) n_z(k)=zero
          !
          weight_omega(k) = weight_omega2(i,j)
          k=k+1
        endif
      enddo
    enddo
    !
    !normalize weights
    weight_omega=weight_omega/four/pi
    !
    !write(*,'(10es20.8)') weight_omega
    !stop
    !
  end subroutine calcnodes_omega8
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calcnodes_omega9(verbose)
    !
    !-----------------------------------------------------------------------
    !   lebedev quadrature with routines from
    !   https://people.sc.fsu.edu/~jburkardt/f_src/sphere_lebedev_rule/sphere_lebedev_rule.html
    !-----------------------------------------------------------------------
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, err, available, n_omega
    integer(i4b) :: dim_mu, dim_phi
    real(dp) :: cb, sb, ca, sa, cc, sc, n_x2, n_y2, n_z2
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: nodes_mu, weight_mu
    !
    ! ... local functions
    integer(i4b) :: order_table, available_table
    !
    if(present(verbose)) ver=verbose
    !
    !only rough estimate for number of integration points
    dim_mu=2*ntheta-1
    dim_phi=2*dim_mu-1
    n_omega=dim_mu*dim_phi
    !
    !allocate nodes_mu and weight mu (although not required, but for output-routines...)
    allocate(nodes_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega9'
    allocate(weight_mu(dim_mu), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega9'
    !
    !
    !check if rule is available and calculate the order of integration scheme
    do i=1, 200
      available = available_table (i)
      nomega=order_table(i)
      if(available.eq.1.and.nomega.ge.n_omega) exit
    enddo
    !
    !calculate lebedev nodes and weights
    allocate(n_x(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega9'
    allocate(n_y(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega9'
    allocate(n_z(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega9'
    allocate(weight_omega(nomega), stat=err)
    if(err.ne.0) stop 'allocation error in calcnodes_omega9'
    !
    !calculate nodes and weights
    call ld_by_order (nomega, n_x, n_y, n_z, weight_omega)
    !
    !!exclude zero directions (by rotation of all directions around x, y, and z-axis)
    ca=0.999999999d0
    cb=ca
    cc=ca
    sa=sqrt(one-ca**2)
    sb=sa
    sc=sa
    !
    do i=1, nomega
      n_x2=cc*cb*n_x(i)+(cc*sb*sa-sc*ca)*n_y(i)+(cc*sb*ca+sc*sa)*n_z(i)
      n_y2=sc*cb*n_x(i)+(sc*sb*sa+cc*ca)*n_y(i)+(sc*sb*ca-cc*sa)*n_z(i)
      n_z2=-sb*n_x(i)+cb*sa*n_y(i)+cb*ca*n_z(i)
      n_x(i)=n_x2
      n_y(i)=n_y2
      n_z(i)=n_z2
    enddo
    !
    !
    !
  end subroutine calcnodes_omega9
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  function n_phi_lobl(del_omega, del_theta, theta)
    !
    !--------calculates number of phi integration nodes for a given---------
    !--(constant) del_omega at mu-level theta with corresponding del_theta--
    !
    !
    ! ... arguments
    integer(i4b) :: n_phi_lobl
    real(dp), intent(in) :: del_omega, del_theta, theta
    !
    ! ... local scalars
    real(dp) :: n_phi_real
    !
    n_phi_real=two*pi*sin(theta)*del_theta / del_omega + one
    n_phi_lobl=nint(n_phi_real)
    !
  end function n_phi_lobl
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine check_nodes_omega(verbose)
    !
    !check if angular grid is reasonable:
    !check for: opposite directions (not necessarily required if many directions are used)
    !           zero directions
    !           normalization
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
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
    if(present(verbose)) ver=verbose
    !
    if(ver) write(*,*) '-------------------------checking the angular grid-----------------------------'
    if(ver) write(*,*)
    !
    !check normalization
    norm=sum(weight_omega)
    if(abs(norm-one).gt.1.d-14) stop 'error in check_nodes_omega: normalization not one'
    !
    !check if only possitive weights occurr
    !check if zero directions occurr
    !check if vector is normalized
    !norm=0.d0
    do i=1, nomega
      !   write(*,*) n_x(i), n_z(i), weight_omega(i)
      !
      norm_vec = n_x(i)**2 + n_y(i)**2 + n_z(i)**2
      if(abs(norm_vec-one).gt.1.d-15) stop 'error in check_nodes_omega: n_x^2+n_y^2+n_z^2 < 1 not allowed'
      !
      if(n_x(i).eq.zero) stop 'error in check_nodes_omega: n_x = 0 not allowed'
      if(n_y(i).eq.zero) stop 'error in check_nodes_omega: n_y = 0 not allowed'
      if(n_z(i).eq.zero) stop 'error in check_nodes_omega: n_z = 0 not allowed'
      !
      if(weight_omega(i).lt.zero) stop 'error in check_nodes_omega: negative weights are not allowed'
      !
    enddo
    !
    !
    !
    !check opposite directions (not so important for many directions)
    do i=1, nomega
      lopposite=.false.
      do j=1, nomega
        opposite = n_x(i)*n_x(j)+n_y(i)*n_y(j)+n_z(i)*n_z(j)
        if(abs(opposite+one).lt.1.d-15) then
          lopposite=.true.
          exit
        endif
      enddo
      if(.not.lopposite) stop 'error in check_nodes_omega: no opposite direction found'
    enddo
    !
    !do i=1, nomega
    !   write(*,*) n_x(i), n_y(i), n_z(i), weight_omega(i)
    !enddo
    !
  end subroutine check_nodes_omega



end module mod_angles
