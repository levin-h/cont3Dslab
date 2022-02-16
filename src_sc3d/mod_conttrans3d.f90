module mod_conttrans3d
  !
  !modules for continuum transport
  !
  use prog_type
  use fund_const
  !
  implicit none
  !
  !temporary output
  character(len=21), parameter :: tempfiles_dir='outputFILES_TEMP/sc3d'
  !
  !-----------------------------------------------------------------------
  !-----------------------iteration-variables-----------------------------
  !-----------------------------------------------------------------------
  !
  integer(i4b), parameter :: itmaxc=1350
  integer(i4b) :: it_start_cont, nconvc
  real(dp), parameter :: devmaxc=1.d-5
  real(dp), dimension(itmaxc) :: epsmaxc_arr
  !devmax(l): maximum required percentage-deviation (continuum, line)
  !
  !-----------------------------------------------------------------------
  !----------------------------ng-extrapolation---------------------------
  !-----------------------------------------------------------------------
  !
  integer(i4b), parameter :: ng_const=6
  !ng_const: ng (or aitken)-extrapolation is performed at each ng_const iterate
  !
contains
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calc_startval_cont(nueindx)
    !
    !-----------------------------------------------------------------------
    !---------calculates start value of continuum source function-----------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, scont3d, opac3d, bnue3d
    use params_input, only: eps_cont
    use mod_frequencies, only: nodes_nue, weight_nue, nnue, lfreqint
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nueindx
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    real(dp) :: opac1, opac2, z1, z2, dtau, tau, thdepth
    real(dp) :: r_therm, t_therm, i_core, dum_bnue, dilfac, zshift
    !
    ! ... local arrays
    !
    ! ... local functions
    real(dp) :: bnue2
    !
    !-----------------------calculate thermalization depth------------------
    !
    if(eps_cont.eq.zero) then
      !set thermalization depth to arbitrary value
      thdepth = 10.d0
    else
      thdepth = one/sqrt(eps_cont)
    endif
    !
    !--------------calculate tau in z-direction (for each frequency)--------
    !
    scont3d=zero
    !
    do i=1, nx
      do j=1, ny
        !
        !set planck function at inner boundary
        scont3d(i,j,1)=bnue3d(i,j,1)
        scont3d(i,j,2)=bnue3d(i,j,2)
        !
        !set zero at outer boundary
        scont3d(i,j,nz)=zero
        scont3d(i,j,nz-1)=zero
        !
        !calculate thermalization in between
        !start at outer boundary with tau=0
        tau=zero
        do k=nz-2, 3, -1
          opac1=opac3d(i,j,k+1)
          opac2=opac3d(i,j,k)
          z1=z(k+1)
          z2=z(k)
          dtau=(opac2+opac1)*(z1-z2)/two
          tau=tau+dtau
          if(tau.gt.thdepth) then
            scont3d(i,j,k)=bnue3d(i,j,k)
          else
            scont3d(i,j,k)=zero
          endif
        enddo
      enddo
    enddo
    !
    !-----------------------------------------------------------------------
    !
    !
    !or just zero initial condition
    scont3d=zero
    do i=1, nx
      do j=1, ny
        scont3d(i,j,1)=bnue3d(i,j,1)
        scont3d(i,j,2)=bnue3d(i,j,2)
      enddo
    enddo
    !
    !
  end subroutine calc_startval_cont
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine conttrans_sc3d(nueindx, verbose)
    !
    !------------------continuum transport for a frequency bin--------------
    !
    !---------calculating intensities for different angles mu---------------
    !-------------------in plane parallel symmetry--------------------------
    !--------------performing angle integration on z-axis-------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, scont3d, mint3d, normalization3d, alocont_nn3d, imask3d, imaskb3d, bnue3d
    use params_input, only: eps_cont
    use mod_frequencies, only: nodes_nue
    use options, only: opt_ng_cont, opt_ait_cont, opt_method
    use mod_interp2d, only: wp_integ1d, wpa_integ1d, wpb_integ1d, wp_interp1d, wpa_interp1d, wpb_interp1d, wp_interp2d, &
         wpa_interp2d, wpb_interp2d, interpolation_threshold
    use mod_timing, only: ts_formal, te_formal, ttot_formal, ts_snew, te_snew, ttot_snew
    use omp_lib
    use mod_math, only: calc_dev3d
    use mod_scont_new3d, only: scont_new3d
    use mod_ng_extrapol, only: store_ng3d, ng_expol3d, ait_expol3d
    !
    ! ... arguments
    integer(i4b), intent(in) :: nueindx
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, l, err, ix, iy, iz, indx_threshold
    integer(i4b) :: s1, s2, s3, s4, s4b
    integer(i4b) :: iim2, iim1, ii, iip1
    real(dp) :: dummy1, dummy2, rad
    real(dp) :: eps_max
    !
    ! ... local arrays
    real(dp), dimension(:,:,:), allocatable :: eps3d
    real(dp), dimension(:,:), allocatable :: scont3d_ng
    !
    !... local characters
    !
    ! ... local functions
    !
    ! ... local logicals
    !
    !-----------------------------------------------------------------------
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '-------------------------continuum transport (3d sc)---------------------------'
    if(ver) write(*,*)
    !
    if(allocated(eps3d)) deallocate(eps3d)
    if(allocated(scont3d_ng)) deallocate(scont3d_ng)
    !
    allocate(eps3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error in conttrans_sc3d'
    allocate(scont3d_ng(4,nx*ny*nz), stat=err)
       if(err.ne.0) stop 'allocation error in conttrans_sc3d'
    !
    !-----------------------------------------------------------------------
    !
    !initialisation of iteration step at which the old solution has to be
    !stored (for ng-extrapolation/aitken extrapolation)
    !
    s1=1
    s2=2
    s3=3
    s4=4
    s4b=4
    !
    !setting index for threshold
    indx_threshold=5
    !
    !-----------------------------------------------------------------------
    !
    !calculating start values of continuum source function
    call calc_startval_cont(nueindx)
    !
    mint3d=zero
    normalization3d=zero
    alocont_nn3d=zero
    !
    eps3d=zero
    scont3d_ng=zero
    !
    ttot_formal=zero
    ttot_snew=zero
    !
    !************************start iteration scheme*************************
    !
    do i=1, itmaxc
      !
      !****************************output*************************************
      !
      ix=nx/2+1
      iy=ny/2+1
      if(ver) write(*,fmt='(a5, 6(a20))') '#', 'z', 'j', 's_c', 'deviation(old-new)', 'normalization', 'alo_diag'
      do k=1, nz
        if(ver) write(*,fmt='(i5, 10(e20.8))') k, z(k), mint3d(ix,iy,k), scont3d(ix,iy,k), eps3d(ix,iy,k), normalization3d(ix,iy,k), alocont_nn3d(ix,iy,k,14)
      end do
      if(ver) write(*,*)
      !
      !*************************output end************************************
      !
      eps3d=mint3d
      !
      call output_itstep_cont(i, verbose=ver)
      !
      !----------------calculating mean intensities on central ray------------
      !-------------------------------step i----------------------------------
      !
      if(ver) write(*,*) '--------------calculating angle integrated mean intensities in 1d--------------'
      if(ver) write(*,*) 'step', i
      if(ver) write(*,*)
      !
      ts_formal=omp_get_wtime()
      !   if(i.le.s4b) then
      !use linear interpolations for the first s4b iteration steps
      !in order to get a first guess of a 'smooth' solution.
      !otherwise: convergence behaviour might be oscillatory
      !      call mint_sc3d(4)
      !   else
      call mint_sc3d(opt_method, nueindx, verbose=ver)
      !   endif
      te_formal=omp_get_wtime()
      ttot_formal=ttot_formal+te_formal-ts_formal
      !
      !-------------calculating percentage-error of mean intensities----------
      !
      if(ver) write(*,*) '-----------------------------calculating deviation-----------------------------'
      if(ver) write(*,*)
      call calc_dev3d(eps3d, mint3d, imask3d, nx, ny, nz, eps_max, ix, iy, iz)
      epsmaxc_arr(i)=eps_max
      if(ver) write(*,'(a30, 3i4, 3f14.4, es18.8)') 'max (dev) at grid-point:', ix, iy, iz, x(ix), y(iy), z(iz) , eps_max
      if(ver) write(*,*)
      !
      nconvc=i
      !
      if(abs(eps_max).lt.devmaxc) then
        if(ver) write(*,*) "convergence after iteration no. ", i
        if(ver) write(*,*) "max (dev): ", eps_max
        if(ver) write(*,*)
        !for thin continua, ensure that higher order interpolation scheme
        !   has been used, and not only linear approach
        if(i.gt.s4b) exit
      else if(i.eq.itmaxc) then
        if(ver) write(*,*) "no convergence after iteration no. ", i
        if(ver) write(*,*)
      end if
      !
      if(i.gt.indx_threshold) then
         if(abs(epsmaxc_arr(i)).ge.abs(epsmaxc_arr(i-1)).and. &
              abs(epsmaxc_arr(i-1)).le.abs(epsmaxc_arr(i-2)).and. &
              abs(epsmaxc_arr(i-2)).ge.abs(epsmaxc_arr(i-3)).and. &
              abs(epsmaxc_arr(i-3)).le.abs(epsmaxc_arr(i-4)).and. &
              abs(epsmaxc_arr(i-4)).ge.abs(epsmaxc_arr(i-5))) then
            !
            if(ver) write(*,*) 'error in iteration scheme: oscillations!!!'
            if(ver) write(*,*) 'max deviation at iteration i-5', epsmaxc_arr(i-5)
            if(ver) write(*,*) 'max deviation at iteration i-4', epsmaxc_arr(i-4)
            if(ver) write(*,*) 'max deviation at iteration i-3', epsmaxc_arr(i-3)
            if(ver) write(*,*) 'max deviation at iteration i-2', epsmaxc_arr(i-2)
            if(ver) write(*,*) 'max deviation at iteration i-1', epsmaxc_arr(i-1)
            if(ver) write(*,*) 'max deviation at iteration i  ', epsmaxc_arr(i)
            if(ver) write(*,*) 'possible solutions: '
            if(ver) write(*,*) '   1. use linear interpolations for upwind/downwind source function'
            if(ver) write(*,*) '      and for upwind intensities to have same solution procedure in'
            if(ver) write(*,*) '      each iteration step (independent of the source function and intensities'
            if(ver) write(*,*) '      themselves.'
            if(ver) write(*,*) '   2. avoid monotonicity constraint in integration of source contribution'
            if(ver) write(*,*) '      (try linear approximation of source function along ray)'
            if(ver) write(*,*) '   3. increase the spatial grid resolution, in order to avoid'
            if(ver) write(*,*) '      monotonicity constraints in quadratic interpolation procedures'
            if(ver) write(*,*) '   4. increasing interpolation threshold in quadratic interpolation procedure'
            interpolation_threshold=min(interpolation_threshold+0.1d0,one)
            if(ver) write(*,*) 'setting interpolation threshold to', interpolation_threshold
            wpa_interp2d=wpa_interp2d+one
            wpb_interp2d=wpb_interp2d+one
            wp_interp2d=wpa_interp2d/wpb_interp2d
            wpa_interp1d=wpa_interp1d+one
            wpb_interp1d=wpb_interp1d+one
            wp_interp1d=wpa_interp1d/wpb_interp1d
            wpa_integ1d=wpa_integ1d+one
            wpb_integ1d=wpb_integ1d+one
            wp_integ1d=wpa_integ1d/wpb_integ1d
            if(ver) write(*,*) 'setting derivative weights for 2d bezier interpolation from, to', (wpa_interp2d-one)/(wpb_interp2d-one), wp_interp2d
            if(ver) write(*,*) 'setting derivative weights for 1d bezier interpolation from, to', (wpa_interp1d-one)/(wpb_interp1d-one), wp_interp1d
            if(ver) write(*,*) 'setting derivative weights for 1d bezier integration   from, to', (wpa_integ1d-one)/(wpb_integ1d-one), wp_integ1d
            !
            !start ng-extrapolation from beginning
            s1=i+1
            s2=i+2
            s3=i+3
            s4=i+4
            indx_threshold=i+5
            !
         endif
      endif
      !
      !-------calculating alo-corrected source-functions on central ray-------
      !
      ts_snew = omp_get_wtime()
      call scont_new3d(verbose=ver)
      !
      !------------extrapolation of old subsequent source functions-----------
      !
      if(opt_ng_cont.or.opt_ait_cont) then
         !
         !----------storing old source-functions for ng-extrapolation------------
         !
         if(i.eq.s1) then
            if(ver) write(*,*) '----------------------storing source fct. at step n-3--------------------------'
            if(ver) write(*,*)
            call store_ng3d(1,scont3d_ng,nx,ny,nz,scont3d, verbose=ver)
            s1=s1+ng_const
         elseif(i.eq.s2) then
            if(ver) write(*,*) '----------------------storing source fct. at step n-2--------------------------'
            if(ver) write(*,*)
            call store_ng3d(2,scont3d_ng,nx,ny,nz,scont3d, verbose=ver)
            s2=s2+ng_const
         elseif(i.eq.s3) then
            if(ver) write(*,*) '----------------------storing source fct. at step n-1--------------------------'
            if(ver) write(*,*)
            call store_ng3d(3,scont3d_ng,nx,ny,nz,scont3d, verbose=ver)
            s3=s3+ng_const
         elseif(i.eq.s4) then
            if(ver) write(*,*) '----------------------storing source fct. at step n----------------------------'
            if(ver) write(*,*)
            call store_ng3d(4,scont3d_ng,nx,ny,nz,scont3d, verbose=ver)
            s4=s4+ng_const
            !
            if(opt_ng_cont) call ng_expol3d(scont3d_ng,nx,ny,nz,scont3d, verbose=ver)
            if(opt_ait_cont) call ait_expol3d(scont3d_ng,nx,ny,nz,scont3d, verbose=ver)
         endif
         !
      endif

      te_snew = omp_get_wtime()
      ttot_snew = ttot_snew+te_snew-ts_snew
      !
      if(ver) write(*,*) 'computation time for formal solution', te_formal-ts_formal
      if(ver) write(*,*) 'computation time for new source fct ', te_snew-ts_snew
      if(ver) write(*,*)
      !   if(i.ge.1) stop 'go on in conttrans'
    enddo
    !
    !make average for computation time
    ttot_snew=ttot_snew/float(nconvc-1)
    ttot_formal=ttot_formal/float(nconvc)
    !
    if(ver) write(*,*)
    if(ver) write(*,*) 'finally used derivative weights for 2d bezier interpolation', wp_interp2d
    if(ver) write(*,*) 'finally used derivative weights for 1d bezier interpolation', wp_interp1d
    if(ver) write(*,*) 'finally used derivative weights for 1d bezier integration  ', wp_integ1d
    if(ver) write(*,*)
    if(ver) write(*,*) 'average computation time for formal solution', ttot_formal
    if(ver) write(*,*) 'average computation time for new source fct ', ttot_snew
    if(ver) write(*,*)
    !
    !
  end subroutine conttrans_sc3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine conttrans_lte3d(nueindx, verbose)
    !
    !----------------set source function to LTE source function-------------
    !-------------------------for one frequency bin-------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, scont3d, imask3d, imaskb3d, bnue3d
    use params_input, only: eps_cont
    use mod_frequencies, only: nodes_nue
    use options, only: opt_method
    use omp_lib
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: nueindx
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    !
    ! ... local arrays
    !
    !... local characters
    !
    ! ... local functions
    !
    ! ... local logicals
    !
    !-----------------------------------------------------------------------
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '-------------------------continuum transport (3d lte)--------------------------'
    if(ver) write(*,*)
    !
    !set up the source function
    scont3d=zero
    !
    do i=1, nx
       do j=1, ny
          do k=1, nz
             scont3d(i,j,k)=bnue3d(i,j,k)
          enddo
       enddo
    enddo
    !
    !calculate radiation moments
    call mint_sc3d(opt_method,nueindx, verbose=ver)
    !
  end subroutine conttrans_lte3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine output_itstep_cont(itnr, verbose)
    !
    !------------------output after each iteration step---------------------
    !input: itnr: current iteration step
    !
    use mod_grid3d, only: scont3d, mint3d
    !
    ! ... arguments
    integer(i4b), intent(in) :: itnr
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver
    !
    ! ... local scalars
    !
    ! ... local arrays
    !
    !-----------------------------------------------------------------------
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '-----------------------storing current iterate in------------------------------'
    if(ver) write(*,*) trim(tempfiles_dir)//'/iterparam_cont.dat'
    if(ver) write(*,*) trim(tempfiles_dir)//'/iterepsmax_cont.dat'
    if(ver) write(*,*) trim(tempfiles_dir)//'/solution_scont3d.dat'
    if(ver) write(*,*) trim(tempfiles_dir)//'/solution_mint3d.dat'
    if(ver) write(*,*)
    !
    !----------------------output iteration parameter-----------------------
    !
    open(1, file=trim(tempfiles_dir)//'/iterparam_cont.dat', form='formatted')
       write(1,'(3(a20))') 'itmaxc', 'current it', 'devmax'
       write(1,'(2i20, e20.8)') itmaxc, itnr, devmaxc
    close(1)
    !
    open(1, file=trim(tempfiles_dir)//'/iterepsmax_cont.dat', form='unformatted')
       write(1) epsmaxc_arr
    close(1)
    !
    !-----------------output 3-d solution grids-----------------------------
    !
    open(1, file=trim(tempfiles_dir)//'/solution_scont3d.dat', form='unformatted')
       write(1) scont3d
    close(1)
    !
    open(1, file=trim(tempfiles_dir)//'/solution_mint3d.dat', form='unformatted')
       write(1) mint3d
    close(1)
    !
    !
    !
  end subroutine output_itstep_cont
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine mint_sc3d(opt_method, nueindx, verbose)
    !
    !-----------------------------------------------------------------------
    !-----------calculates mean intensities at all grid point---------------
    !-----------------------------------------------------------------------
    !
    !   formal solution from 3d short characteristics
    !
    !                  mu-integration between [-1,1]
    !                 phi-integration between [0,2*pi]
    !
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: int3d, alocont_o_nn3d, alocont_nn3d, mint3d, fcontx3d, fconty3d, fcontz3d, normalization3d, &
         alocont_nn3d_tmp, mint3d_tmp, fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp, &
         kcontxx3d, kcontxy3d, kcontxz3d, &
         kcontyy3d, kcontyz3d, &
         kcontzz3d, &
         normalization3d_tmp, bnue3d, tgas3d, trad3d, &
         nx, ny, nz, imask3d, scont3d
    use mod_angles, only: nomega
    use mod_frequencies, only: nodes_nue, nnue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: opt_method, nueindx
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: oindx
    integer(i4b) :: err
    !
    ! ... local arrays
    !
    ! ... local functions
    real(dp) :: bnue
    !
    ! ... local characters
    !
    if(present(verbose)) ver = verbose
    !
    alocont_nn3d=zero
    normalization3d=zero
    mint3d=zero
    fcontx3d=zero
    fconty3d=zero
    fcontz3d=zero
    kcontxx3d=zero
    kcontxy3d=zero
    kcontxz3d=zero
    kcontyy3d=zero
    kcontyz3d=zero
    kcontzz3d=zero
    !
    !--------------deallocation of global (threadprivate) arrays------------
    !
    if(allocated(int3d)) deallocate(int3d)
    if(allocated(alocont_o_nn3d)) deallocate(alocont_o_nn3d)
    !
    !
    !-----------------------begin of parallel region------------------------
    !
    !$omp parallel &
    !$omp private(err, oindx)
    !
    !-------------allocation of global (threadprivate) arrays---------------
    !
    allocate(alocont_o_nn3d(nx,ny,nz,27), stat=err)
       if(err.ne.0) stop 'allocation error in mint_sc3d: alocont_o_nn3d'
    allocate(alocont_nn3d_tmp(nx,ny,nz,27), stat=err)
       if(err.ne.0) stop 'allocation error in mint_sc3d: alocont_nn3d_tmp'
    allocate(normalization3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: normalization3d_tmp'
    allocate(mint3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: mint3d_tmp'
    allocate(fcontx3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: fcontx3d_tmp'
    allocate(fconty3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: fconty3d_tmp'
    allocate(fcontz3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: fcontz3d_tmp'
    allocate(kcontxx3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: kcontxx3d_tmp'
    allocate(kcontxy3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: kcontxy3d_tmp'
    allocate(kcontxz3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: kcontxz3d_tmp'
    allocate(kcontyy3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: kcontyy3d_tmp'
    allocate(kcontyz3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: kcontyz3d_tmp'
    allocate(kcontzz3d_tmp(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d: kcontzz3d_tmp'
    allocate(int3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error mint_sc3d:: int3d'
    !
    alocont_nn3d_tmp=zero
    normalization3d_tmp=zero
    mint3d_tmp=zero
    fcontx3d_tmp=zero
    fconty3d_tmp=zero
    fcontz3d_tmp=zero
    kcontxx3d_tmp=zero
    kcontxy3d_tmp=zero
    kcontxz3d_tmp=zero
    kcontyy3d_tmp=zero
    kcontyz3d_tmp=zero
    kcontzz3d_tmp=zero
    !
    !$omp do schedule(dynamic)
    do oindx=nomega, 1, -1
       !      if(ver) write(*,'(a30,2i10,a5,i5)') 'calculating omega', oindx, nomega, 'bez'
       select case(opt_method)
          case(4,14)
             call formalsc_cont3d_lin(oindx,nueindx)
          case(5,15)
             call formalsc_cont3d(oindx,nueindx)
          case(6,16)
             call formallc_cont3d_lin(oindx,nueindx)
          case(7,17)
             call formallc_cont3d(oindx,nueindx)
          case default
             stop 'error in mint_sc3d: opt_method not defined'
       end select
    enddo
    !$omp enddo
    !
    !------------------------add up temporary arrays------------------------
    !
    !$omp critical
    alocont_nn3d = alocont_nn3d+alocont_nn3d_tmp
    normalization3d = normalization3d + normalization3d_tmp
    mint3d = mint3d + mint3d_tmp
    fcontx3d = fcontx3d + fcontx3d_tmp
    fconty3d = fconty3d + fconty3d_tmp
    fcontz3d = fcontz3d + fcontz3d_tmp

    kcontxx3d = kcontxx3d + kcontxx3d_tmp
    kcontxy3d = kcontxy3d + kcontxy3d_tmp
    kcontxz3d = kcontxz3d + kcontxz3d_tmp
    kcontyy3d = kcontyy3d + kcontyy3d_tmp
    kcontyz3d = kcontyz3d + kcontyz3d_tmp
    kcontzz3d = kcontzz3d + kcontzz3d_tmp
    !$omp end critical
    !
    !------------deallocation of global (threadprivate) arrays--------------
    !
    deallocate(alocont_nn3d_tmp)
    deallocate(normalization3d_tmp)
    deallocate(mint3d_tmp)
    deallocate(fcontx3d_tmp)
    deallocate(fconty3d_tmp)
    deallocate(fcontz3d_tmp)
    deallocate(kcontxx3d_tmp)
    deallocate(kcontxy3d_tmp)
    deallocate(kcontxz3d_tmp)
    deallocate(kcontyy3d_tmp)
    deallocate(kcontyz3d_tmp)
    deallocate(kcontzz3d_tmp)
    deallocate(alocont_o_nn3d)
    deallocate(int3d)
    !
    !$omp end parallel
    !
    !renormalize (note that mean intensities also calculated at ghost points on left and right boundary
    do i=1, nx
       do j=1, ny
          do k=1, nz
             select case(imask3d(i,j,k))
                case(1,2)
                   if(abs(normalization3d(i,j,k)-one).gt.small_number) then
                      write(*,'(a50, 2i5, es20.8)') 'normalization error in mint_sc3d:', i, k, normalization3d(i,j,k)
                      stop
                   endif
                   mint3d(i,j,k) = mint3d(i,j,k)/normalization3d(i,j,k)
                   alocont_nn3d(i,j,k,:) = alocont_nn3d(i,j,k,:)/normalization3d(i,j,k)
                   normalization3d(i,j,k) = normalization3d(i,j,k)/normalization3d(i,j,k)
                   fcontx3d(i,j,k) = fcontx3d(i,j,k)/normalization3d(i,j,k)
                   fconty3d(i,j,k) = fconty3d(i,j,k)/normalization3d(i,j,k)
                   fcontz3d(i,j,k) = fcontz3d(i,j,k)/normalization3d(i,j,k)
                   kcontxx3d(i,j,k) = kcontxx3d(i,j,k)/normalization3d(i,j,k)
                   kcontxy3d(i,j,k) = kcontxy3d(i,j,k)/normalization3d(i,j,k)
                   kcontxz3d(i,j,k) = kcontxz3d(i,j,k)/normalization3d(i,j,k)
                   kcontyy3d(i,j,k) = kcontyy3d(i,j,k)/normalization3d(i,j,k)
                   kcontyz3d(i,j,k) = kcontyz3d(i,j,k)/normalization3d(i,j,k)
                   kcontzz3d(i,j,k) = kcontzz3d(i,j,k)/normalization3d(i,j,k)
                case default
             end select
          enddo
       enddo
    enddo
    !
    !
    !
  end subroutine mint_sc3d


  !
  !***********************************************************************
  !***********************************************************************
  !
  !        SOME ROUTINES USED FOR 3D CONTINUUM AND LINE TRANSPORT
  !
  !***********************************************************************
  !***********************************************************************
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine set_boundary3d(nue, nn_x, nn_y, nn_z)
    !
    use prog_type
    use fund_const
    use mod_grid3d, only: nx, ny, nz, bnue3d, int3d
    use mod_bcondition3d, only: xic1_3d, xic2_3d
    !
    implicit none
    !
    ! ... arguments
    real(dp), intent(in) :: nue, nn_x, nn_y, nn_z
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    !
    ! ... local functions
    !
    int3d=zero
    !
    do i=1, nx
      do j=1, ny
!        int3d(i,j,1) = bnue3d(i,j,1)
         !        int3d(i,j,2) = bnue3d(i,j,2)
        int3d(i,j,1) = xic1_3d(i,j,1) - nn_z * xic2_3d(i,j,1)
        int3d(i,j,2) = xic1_3d(i,j,2) - nn_z * xic2_3d(i,j,2)
      enddo
   enddo
    !
    !
  end subroutine set_boundary3d
  !
  !***********************************************************************
  !***********************************************************************
  !
  !                 CONTINUUM TRASNFER ROUTINES
  !
  !***********************************************************************
  !***********************************************************************
  !
  subroutine formalsc_cont3d(oindx,nueindx)
    !
    !-----------------------------------------------------------------------
    !------short characteristics for continuum radiative transfer in 3d-----
    !-----------calculating intensties for given mu,phi specified-----------
    !---------------------------by input oindx------------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, int3d, opac3d, scont3d, imaskb3d, &
         alocont_o_nn3d, alocont_nn3d_tmp, mint3d_tmp, normalization3d_tmp, &
         fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, bnue3d, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp
    use mod_angles, only: n_x, n_y, n_z, weight_omega, q_alo
    use mod_frequencies, only: nodes_nue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: oindx, nueindx
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: alpha, beta, gamma
    integer(i4b) :: startx, starty, startz, endx, endy, endz, &
    startxb, startyb, endxb, endyb
    integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
    integer(i4b) :: aloindx_im1, aloindx_ip1, aloindx_jm1, aloindx_jp1, aloindx_km1, aloindx_kp1
    real(dp) :: nn_x, nn_y, nn_z, xnue, wall
    real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
    real(dp) :: dels_xyd, dels_xzd, dels_yzd, dels_d
    integer :: q1, q2, q3, q4, q5, q6, q7, q8, q9, &
    q10, q11, q12, q13, q14, q15, q16, q17, q18, &
    q19, q20, q21, q22, q23, q24, q25, q26, q27
    !
    !for debugging
    real(dp) :: int_u2, scont_u2, opac_u2, int_u3, scont_u3, opac_u3, scont_d2, scont_d3, opac_d2, opac_d3, interpol2d_9p_quad, interpol2d_9p_bez
    !
    ! ... local arrays
    !
    ! ... local functions
    real(dp) :: dist_ellipsoid, calc_icore_gdark
    !
    ! ... local characters
    !character(len=50) :: enter
    !
    ! ... local logicals
    !
    !frequency
    xnue=nodes_nue(nueindx)
    !
    !directions
    nn_x=n_x(oindx)
    nn_y=n_y(oindx)
    nn_z=n_z(oindx)
    !
    !angulare integration weight
    wall=weight_omega(oindx)
    !
    !indices for nearest neighbour alo
    q1=q_alo(oindx,1)
    q2=q_alo(oindx,2)
    q3=q_alo(oindx,3)
    q4=q_alo(oindx,4)
    q5=q_alo(oindx,5)
    q6=q_alo(oindx,6)
    q7=q_alo(oindx,7)
    q8=q_alo(oindx,8)
    q9=q_alo(oindx,9)
    q10=q_alo(oindx,10)
    q11=q_alo(oindx,11)
    q12=q_alo(oindx,12)
    q13=q_alo(oindx,13)
    q14=q_alo(oindx,14)
    q15=q_alo(oindx,15)
    q16=q_alo(oindx,16)
    q17=q_alo(oindx,17)
    q18=q_alo(oindx,18)
    q19=q_alo(oindx,19)
    q20=q_alo(oindx,20)
    q21=q_alo(oindx,21)
    q22=q_alo(oindx,22)
    q23=q_alo(oindx,23)
    q24=q_alo(oindx,24)
    q25=q_alo(oindx,25)
    q26=q_alo(oindx,26)
    q27=q_alo(oindx,27)
    !
    !-----------------------------------------------------------------------
    !
    !set directional index-parameter (phantom points are excluded from calculation)
    !
    !index parameter:
    !         if n_x >= 0                 if n_x < 0
    !                startx = 2                  startx = ndxmax-1
    !                endx = ndxmax-1             endx = 2
    !                alpha=  1                   alpha=-1
    !
    !         if n_y >= 0                 if n_y < 0
    !                starty = 2                  starty = ndymax-1
    !                endy = ndymax-1             endy = 2
    !                beta =  1                   beta =-1
    !
    !         if n_z >= 0                 if n_z < 0
    !                startz = 2                  startz = ndzmax-1
    !                endz = ndzmax-1             endz = 2
    !                gamma = 1                   gamma = -1
    !
    if(nn_x.gt.zero) then
      startx = 1
      endx = nx
      startxb = nx
      endxb = 1
      alpha=  1
      aloindx_im1 = nx-1
      aloindx_ip1 = 2
    elseif(nn_x.lt.zero) then
      startx = nx
      endx = 1
      startxb = 1
      endxb = nx
      alpha=-1
      aloindx_im1 = 2
      aloindx_ip1 = nx-1
    else
      stop 'error in formalsc_cont3d: n_x = 0 not allowed'
    endif
    !
    if(nn_y.gt.zero) then
      starty = 1
      endy = ny
      startyb = ny
      endyb = 1
      beta =  1
      aloindx_jm1 = ny-1
      aloindx_jp1 = 2
    elseif(nn_y.lt.zero) then
      starty = ny
      endy = 1
      startyb = 1
      endyb = ny
      beta =-1
      aloindx_jm1 = 2
      aloindx_jp1 = ny-1
    else
      stop 'error in formalsc_cont3d: n_y = 0 not allowed'
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
    else
      stop 'error in formalsc_cont3d: n_z = 0 not allowed'
    endif
    !
    !--------------------reset the intensities and alo----------------------
    !
    call set_boundary3d(xnue, nn_x, nn_y, nn_z)
    !
    !-----------------------------------------------------------------------
    !
    alocont_o_nn3d=zero
    !
    !-----------------------------------------------------------------------
    !
    do k=startz, endz, gamma
      do j=starty, endy, beta
        do i=startx, endx, alpha
          !
          select case(imaskb3d(i,j,k))
            !
            !*************left and right boundary with long characteristics*********
            !
          case(12,13)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
            !*************front and back boundary with long characteristics*********
            !
          case(14,15)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
            !********(left,front), (right,front), (left,back), (right,back)*********
            !*******************boundary edges with long characteristics************
            !
          case(18,21,24,27)
            !
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
            !***********************adjacent to left boundary***********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(4)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on left towards right
            alocont_o_nn3d(nx,j-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,j-1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,j-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,j+1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx,j+1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx,j+1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,j  ,k+1, 6) = alocont_o_nn3d(i-1,j  ,k+1, 6)
            alocont_o_nn3d(nx,j  ,k  ,15) = alocont_o_nn3d(i-1,j  ,k  ,15)
            alocont_o_nn3d(nx,j  ,k-1,24) = alocont_o_nn3d(i-1,j  ,k-1,24)
            !
            !***********************adjacent to right boundary**********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(5)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on right towards left
            alocont_o_nn3d(1,j-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,j-1,k,  16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,j-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,j+1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,j+1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,j+1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,j  ,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,j  ,k,  13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,j  ,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            !
            !***********************adjacent to front boundary**********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(6)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards front
            alocont_o_nn3d(i-1,ny,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(i-1,ny,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(i-1,ny,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(i  ,ny,k+1, 8) = alocont_o_nn3d(i  ,j-1,k+1, 8)
            alocont_o_nn3d(i  ,ny,k  ,17) = alocont_o_nn3d(i  ,j-1,k  ,17)
            alocont_o_nn3d(i  ,ny,k-1,26) = alocont_o_nn3d(i  ,j-1,k-1,26)
            alocont_o_nn3d(i+1,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(i+1,ny,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(i+1,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            !
            !***********************adjacent to back boundary***********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(7)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards front
            alocont_o_nn3d(i+1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(i+1,1,k,  10) = alocont_o_nn3d(i+1,j+1,k,  10)
            alocont_o_nn3d(i+1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(i  ,1,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(i  ,1,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(i  ,1,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(i-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(i-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(i-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            !
            !***********************adjacent to left front boundary*****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(8)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(nx,1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,2,k+1, 6) = alocont_o_nn3d(i-1,j  ,k+1, 6)
            alocont_o_nn3d(nx,2,k  ,15) = alocont_o_nn3d(i-1,j  ,k  ,15)
            alocont_o_nn3d(nx,2,k-1,24) = alocont_o_nn3d(i-1,j  ,k-1,24)
            alocont_o_nn3d(nx,j+1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx,j+1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx,j+1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(1,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(1,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(1,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(2,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(2,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(2,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(3,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(3,ny,k,  16) = alocont_o_nn3d(i+1,j-1,k,  16)
            alocont_o_nn3d(3,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx ,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx ,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx ,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            !
            !***********************adjacent to right front boundary****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(9)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)

            alocont_o_nn3d(1,2,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,2,k  ,13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,2,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            alocont_o_nn3d(1,j+1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,j+1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,j+1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,1,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,ny,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx,ny ,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(nx,ny ,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(nx,ny ,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx-2,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx-2,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx-2,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx-1,ny ,k+1, 8) = alocont_o_nn3d(i  ,j-1,k+1, 8)
            alocont_o_nn3d(nx-1,ny ,k  ,17) = alocont_o_nn3d(i  ,j-1,k  ,17)
            alocont_o_nn3d(nx-1,ny ,k-1,26) = alocont_o_nn3d(i  ,j-1,k-1,26)

            !
            !***********************adjacent to left back boundary******************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(10)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1  ,1  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(1  ,1  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(1  ,1  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(2  ,1  ,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(2  ,1  ,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(2  ,1  ,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(nx  ,1  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx  ,1  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx  ,1  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,j-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,j-1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,j-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx  ,ny  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx  ,ny  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx  ,ny  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx ,ny-1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx ,ny-1,k,  15) = alocont_o_nn3d(i-1,j,k,  15)
            alocont_o_nn3d(nx ,ny-1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(i+1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(i+1,1,k,  10) = alocont_o_nn3d(i+1,j+1,k,  10)
            alocont_o_nn3d(i+1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            !
            !***********************adjacent to right back boundary*****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(11)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(nx-1,1,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(nx-1,1,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(nx-1,1,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(nx,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(nx,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(nx,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,j-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,j-1,k,  16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,j-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,j  ,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,j  ,k  ,13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,j  ,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            alocont_o_nn3d(1,ny,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,ny,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,ny,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(i-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(i-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(i-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            !
            !**********************left boundary with long characteristics**********
            !
          case(16)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
            alocont_o_nn3d(2,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(2,ny,k,  16) = alocont_o_nn3d(i+1,j-1,k,  16)
            alocont_o_nn3d(2,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(nx,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(nx,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(1,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(1,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(1,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            !
            !*********************front boundary with long characteristics**********
            !
          case(17)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            alocont_o_nn3d(nx,1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,1,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,2,k+1,3) = alocont_o_nn3d(i-1,j+1,k+1,3)
            alocont_o_nn3d(nx,2,k,12) = alocont_o_nn3d(i-1,j+1,k,12)
            alocont_o_nn3d(nx,2,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,ny,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,ny,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,ny,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            !
            !*********************right boundary with long characteristics**********
            !
          case(19)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            alocont_o_nn3d(nx-1,ny,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx-1,ny,k,  18) = alocont_o_nn3d(i-1,j-1,k,  18)
            alocont_o_nn3d(nx-1,ny,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(nx,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(nx,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(1,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(1,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(1,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            !
            !*********************front boundary with long characteristics**********
            !
          case(20)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,1,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,1,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            alocont_o_nn3d(1,2,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,2,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,2,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,ny,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,ny,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,ny,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            !
            !**********************left boundary with long characteristics**********
            !
          case(22)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(1,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(1,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            alocont_o_nn3d(2,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(2,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(2,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(nx,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(nx,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(nx,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            !
            !*********************back boundary with long characteristics***********
            !
          case(23)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            alocont_o_nn3d(nx,1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,1,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,ny,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,ny,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,ny,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,ny-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,ny-1,k,  18) = alocont_o_nn3d(i-1,j-1,k,  18)
            alocont_o_nn3d(nx,ny-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            !
            !**********************right boundary with long characteristics*********
            !
          case(25)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(1,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(1,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            alocont_o_nn3d(nx-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(nx,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(nx,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            !
            !**********************back boundary with long characteristics**********
            !
          case(26)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,1,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,1,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            alocont_o_nn3d(1,ny-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,ny-1,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,ny-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,ny,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,ny,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,ny,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            !
            !************************standard radiative transfer********************
            !
          case(3)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            !
            call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
          end select
          !
          !
        enddo
      enddo
    enddo
    !
    !---------------------perform angular integration-----------------------
    !NOTE: due to periodic boundary conditions, and copying of ALO coefficients
    !      angular integration cannot (easily) be performed within loop from above
    !
    do k=1, nz
      do j=1, ny
        do i=1, nx
          mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
          fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
          fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
          fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
          kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
          kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
          kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
          kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
          kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
          kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
          normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall

          alocont_nn3d_tmp(i,j,k,1) = alocont_nn3d_tmp(i,j,k,1) + wall*alocont_o_nn3d(i,j,k,1)
          alocont_nn3d_tmp(i,j,k,2) = alocont_nn3d_tmp(i,j,k,2) + wall*alocont_o_nn3d(i,j,k,2)
          alocont_nn3d_tmp(i,j,k,3) = alocont_nn3d_tmp(i,j,k,3) + wall*alocont_o_nn3d(i,j,k,3)
          alocont_nn3d_tmp(i,j,k,4) = alocont_nn3d_tmp(i,j,k,4) + wall*alocont_o_nn3d(i,j,k,4)
          alocont_nn3d_tmp(i,j,k,5) = alocont_nn3d_tmp(i,j,k,5) + wall*alocont_o_nn3d(i,j,k,5)
          alocont_nn3d_tmp(i,j,k,6) = alocont_nn3d_tmp(i,j,k,6) + wall*alocont_o_nn3d(i,j,k,6)
          alocont_nn3d_tmp(i,j,k,7) = alocont_nn3d_tmp(i,j,k,7) + wall*alocont_o_nn3d(i,j,k,7)
          alocont_nn3d_tmp(i,j,k,8) = alocont_nn3d_tmp(i,j,k,8) + wall*alocont_o_nn3d(i,j,k,8)
          alocont_nn3d_tmp(i,j,k,9) = alocont_nn3d_tmp(i,j,k,9) + wall*alocont_o_nn3d(i,j,k,9)
          alocont_nn3d_tmp(i,j,k,10) = alocont_nn3d_tmp(i,j,k,10) + wall*alocont_o_nn3d(i,j,k,10)
          alocont_nn3d_tmp(i,j,k,11) = alocont_nn3d_tmp(i,j,k,11) + wall*alocont_o_nn3d(i,j,k,11)
          alocont_nn3d_tmp(i,j,k,12) = alocont_nn3d_tmp(i,j,k,12) + wall*alocont_o_nn3d(i,j,k,12)
          alocont_nn3d_tmp(i,j,k,13) = alocont_nn3d_tmp(i,j,k,13) + wall*alocont_o_nn3d(i,j,k,13)
          alocont_nn3d_tmp(i,j,k,14) = alocont_nn3d_tmp(i,j,k,14) + wall*alocont_o_nn3d(i,j,k,14)
          alocont_nn3d_tmp(i,j,k,15) = alocont_nn3d_tmp(i,j,k,15) + wall*alocont_o_nn3d(i,j,k,15)
          alocont_nn3d_tmp(i,j,k,16) = alocont_nn3d_tmp(i,j,k,16) + wall*alocont_o_nn3d(i,j,k,16)
          alocont_nn3d_tmp(i,j,k,17) = alocont_nn3d_tmp(i,j,k,17) + wall*alocont_o_nn3d(i,j,k,17)
          alocont_nn3d_tmp(i,j,k,18) = alocont_nn3d_tmp(i,j,k,18) + wall*alocont_o_nn3d(i,j,k,18)
          alocont_nn3d_tmp(i,j,k,19) = alocont_nn3d_tmp(i,j,k,19) + wall*alocont_o_nn3d(i,j,k,19)
          alocont_nn3d_tmp(i,j,k,20) = alocont_nn3d_tmp(i,j,k,20) + wall*alocont_o_nn3d(i,j,k,20)
          alocont_nn3d_tmp(i,j,k,21) = alocont_nn3d_tmp(i,j,k,21) + wall*alocont_o_nn3d(i,j,k,21)
          alocont_nn3d_tmp(i,j,k,22) = alocont_nn3d_tmp(i,j,k,22) + wall*alocont_o_nn3d(i,j,k,22)
          alocont_nn3d_tmp(i,j,k,23) = alocont_nn3d_tmp(i,j,k,23) + wall*alocont_o_nn3d(i,j,k,23)
          alocont_nn3d_tmp(i,j,k,24) = alocont_nn3d_tmp(i,j,k,24) + wall*alocont_o_nn3d(i,j,k,24)
          alocont_nn3d_tmp(i,j,k,25) = alocont_nn3d_tmp(i,j,k,25) + wall*alocont_o_nn3d(i,j,k,25)
          alocont_nn3d_tmp(i,j,k,26) = alocont_nn3d_tmp(i,j,k,26) + wall*alocont_o_nn3d(i,j,k,26)
          alocont_nn3d_tmp(i,j,k,27) = alocont_nn3d_tmp(i,j,k,27) + wall*alocont_o_nn3d(i,j,k,27)
        enddo
      enddo
    enddo
    !
  end subroutine formalsc_cont3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine formalsc_cont3d_lin(oindx,nueindx)
    !
    !-----------------------------------------------------------------------
    !------short characteristics for continuum radiative transfer in 3d-----
    !-----------calculating intensties for given mu,phi specified-----------
    !---------------------------by input oindx------------------------------
    !-----------------------------------------------------------------------
    !
    use prog_type
    use fund_const
    use mod_grid3d, only: nx, ny, nz, x, y, z, int3d, opac3d, scont3d, imaskb3d, &
         alocont_o_nn3d, alocont_nn3d_tmp, mint3d_tmp, normalization3d_tmp, &
         fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, bnue3d, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp
    use mod_angles, only: n_x, n_y, n_z, weight_omega, q_alo
    use mod_frequencies, only: nodes_nue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: oindx, nueindx
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: alpha, beta, gamma
    integer(i4b) :: startx, starty, startz, endx, endy, endz, &
    startxb, startyb, endxb, endyb
    integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
    integer(i4b) :: aloindx_im1, aloindx_ip1, aloindx_jm1, aloindx_jp1, aloindx_km1, aloindx_kp1
    real(dp) :: nn_x, nn_y, nn_z, xnue, wall
    real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
    real(dp) :: dels_xyd, dels_xzd, dels_yzd, dels_d
    integer :: q14, q15, q17, q18, q23, q24, q26, q27
    !
    !for debugging
    real(dp) :: int_u2, scont_u2, opac_u2, int_u3, scont_u3, opac_u3, scont_d2, scont_d3, opac_d2, opac_d3, interpol2d_9p_quad, interpol2d_9p_bez
    !
    ! ... local arrays
    !
    ! ... local functions
    real(dp) :: dist_ellipsoid, calc_icore_gdark
    !
    ! ... local characters
    !character(len=50) :: enter
    !
    ! ... local logicals
    !
    !frequency
    xnue=nodes_nue(nueindx)
    !
    !directions
    nn_x=n_x(oindx)
    nn_y=n_y(oindx)
    nn_z=n_z(oindx)
    !
    !angulare integration weight
    wall=weight_omega(oindx)
    !
    !indices for nearest neighbour alo
    q14=q_alo(oindx,14)
    q15=q_alo(oindx,15)
    q17=q_alo(oindx,17)
    q18=q_alo(oindx,18)
    q23=q_alo(oindx,23)
    q24=q_alo(oindx,24)
    q26=q_alo(oindx,26)
    q27=q_alo(oindx,27)
    !
    !-----------------------------------------------------------------------
    !
    !set directional index-parameter (phantom points are excluded from calculation)
    !
    !index parameter:
    !         if n_x >= 0                 if n_x < 0
    !                startx = 2                  startx = ndxmax-1
    !                endx = ndxmax-1             endx = 2
    !                alpha=  1                   alpha=-1
    !
    !         if n_y >= 0                 if n_y < 0
    !                starty = 2                  starty = ndymax-1
    !                endy = ndymax-1             endy = 2
    !                beta =  1                   beta =-1
    !
    !         if n_z >= 0                 if n_z < 0
    !                startz = 2                  startz = ndzmax-1
    !                endz = ndzmax-1             endz = 2
    !                gamma = 1                   gamma = -1
    !
    if(nn_x.gt.zero) then
      startx = 1
      endx = nx
      startxb = nx
      endxb = 1
      alpha=  1
      aloindx_im1 = nx-1
      aloindx_ip1 = 2
    elseif(nn_x.lt.zero) then
      startx = nx
      endx = 1
      startxb = 1
      endxb = nx
      alpha=-1
      aloindx_im1 = 2
      aloindx_ip1 = nx-1
    else
      stop 'error in formalsc_cont3d: n_x = 0 not allowed'
    endif
    !
    if(nn_y.gt.zero) then
      starty = 1
      endy = ny
      startyb = ny
      endyb = 1
      beta =  1
      aloindx_jm1 = ny-1
      aloindx_jp1 = 2
    elseif(nn_y.lt.zero) then
      starty = ny
      endy = 1
      startyb = 1
      endyb = ny
      beta =-1
      aloindx_jm1 = 2
      aloindx_jp1 = ny-1
    else
      stop 'error in formalsc_cont3d: n_y = 0 not allowed'
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
    else
      stop 'error in formalsc_cont3d: n_z = 0 not allowed'
    endif
    !
    !--------------------reset the intensities and alo----------------------
    !
    call set_boundary3d(xnue, nn_x, nn_y, nn_z)
    !
    !-----------------------------------------------------------------------
    !
    alocont_o_nn3d=zero
    !
    !-----------------------------------------------------------------------
    !
    do k=startz, endz, gamma
      do j=starty, endy, beta
        do i=startx, endx, alpha
          !
          select case(imaskb3d(i,j,k))
            !
            !*************left and right boundary with long characteristics*********
            !
          case(12,13)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
            !*************front and back boundary with long characteristics*********
            !
          case(14,15)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
            !********(left,front), (right,front), (left,back), (right,back)*********
            !*******************boundary edges with long characteristics************
            !
          case(18,21,24,27)
            !
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
            !***********************adjacent to left boundary***********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(4)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on left towards right
            alocont_o_nn3d(nx,j-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,j-1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,j-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,j+1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx,j+1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx,j+1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,j  ,k+1, 6) = alocont_o_nn3d(i-1,j  ,k+1, 6)
            alocont_o_nn3d(nx,j  ,k  ,15) = alocont_o_nn3d(i-1,j  ,k  ,15)
            alocont_o_nn3d(nx,j  ,k-1,24) = alocont_o_nn3d(i-1,j  ,k-1,24)
            !
            !***********************adjacent to right boundary**********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(5)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on right towards left
            alocont_o_nn3d(1,j-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,j-1,k,  16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,j-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,j+1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,j+1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,j+1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,j  ,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,j  ,k,  13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,j  ,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            !
            !***********************adjacent to front boundary**********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(6)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards front
            alocont_o_nn3d(i-1,ny,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(i-1,ny,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(i-1,ny,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(i  ,ny,k+1, 8) = alocont_o_nn3d(i  ,j-1,k+1, 8)
            alocont_o_nn3d(i  ,ny,k  ,17) = alocont_o_nn3d(i  ,j-1,k  ,17)
            alocont_o_nn3d(i  ,ny,k-1,26) = alocont_o_nn3d(i  ,j-1,k-1,26)
            alocont_o_nn3d(i+1,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(i+1,ny,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(i+1,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            !
            !***********************adjacent to back boundary***********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(7)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards front
            alocont_o_nn3d(i+1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(i+1,1,k,  10) = alocont_o_nn3d(i+1,j+1,k,  10)
            alocont_o_nn3d(i+1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(i  ,1,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(i  ,1,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(i  ,1,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(i-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(i-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(i-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            !
            !***********************adjacent to left front boundary*****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(8)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(nx,1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,2,k+1, 6) = alocont_o_nn3d(i-1,j  ,k+1, 6)
            alocont_o_nn3d(nx,2,k  ,15) = alocont_o_nn3d(i-1,j  ,k  ,15)
            alocont_o_nn3d(nx,2,k-1,24) = alocont_o_nn3d(i-1,j  ,k-1,24)
            alocont_o_nn3d(nx,j+1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx,j+1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx,j+1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(1,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(1,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(1,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(2,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(2,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(2,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(3,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(3,ny,k,  16) = alocont_o_nn3d(i+1,j-1,k,  16)
            alocont_o_nn3d(3,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx ,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx ,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx ,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            !
            !***********************adjacent to right front boundary****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(9)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)

            alocont_o_nn3d(1,2,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,2,k  ,13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,2,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            alocont_o_nn3d(1,j+1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,j+1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,j+1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,1,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,ny,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx,ny ,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(nx,ny ,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(nx,ny ,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx-2,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx-2,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx-2,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx-1,ny ,k+1, 8) = alocont_o_nn3d(i  ,j-1,k+1, 8)
            alocont_o_nn3d(nx-1,ny ,k  ,17) = alocont_o_nn3d(i  ,j-1,k  ,17)
            alocont_o_nn3d(nx-1,ny ,k-1,26) = alocont_o_nn3d(i  ,j-1,k-1,26)

            !
            !***********************adjacent to left back boundary******************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(10)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1  ,1  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(1  ,1  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(1  ,1  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(2  ,1  ,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(2  ,1  ,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(2  ,1  ,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(nx  ,1  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx  ,1  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx  ,1  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,j-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,j-1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,j-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx  ,ny  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx  ,ny  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx  ,ny  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx ,ny-1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx ,ny-1,k,  15) = alocont_o_nn3d(i-1,j,k,  15)
            alocont_o_nn3d(nx ,ny-1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(i+1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(i+1,1,k,  10) = alocont_o_nn3d(i+1,j+1,k,  10)
            alocont_o_nn3d(i+1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            !
            !***********************adjacent to right back boundary*****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(11)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(nx-1,1,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(nx-1,1,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(nx-1,1,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(nx,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(nx,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(nx,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,j-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,j-1,k,  16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,j-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,j  ,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,j  ,k  ,13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,j  ,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            alocont_o_nn3d(1,ny,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,ny,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,ny,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(i-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(i-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(i-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            !
            !**********************left boundary with long characteristics**********
            !
          case(16)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
            alocont_o_nn3d(2,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(2,ny,k,  16) = alocont_o_nn3d(i+1,j-1,k,  16)
            alocont_o_nn3d(2,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(nx,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(nx,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(1,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(1,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(1,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            !
            !*********************front boundary with long characteristics**********
            !
          case(17)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            alocont_o_nn3d(nx,1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,1,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,2,k+1,3) = alocont_o_nn3d(i-1,j+1,k+1,3)
            alocont_o_nn3d(nx,2,k,12) = alocont_o_nn3d(i-1,j+1,k,12)
            alocont_o_nn3d(nx,2,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,ny,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,ny,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,ny,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            !
            !*********************right boundary with long characteristics**********
            !
          case(19)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            alocont_o_nn3d(nx-1,ny,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx-1,ny,k,  18) = alocont_o_nn3d(i-1,j-1,k,  18)
            alocont_o_nn3d(nx-1,ny,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(nx,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(nx,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(1,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(1,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(1,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            !
            !*********************front boundary with long characteristics**********
            !
          case(20)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,1,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,1,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            alocont_o_nn3d(1,2,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,2,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,2,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,ny,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,ny,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,ny,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            !
            !**********************left boundary with long characteristics**********
            !
          case(22)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(1,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(1,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            alocont_o_nn3d(2,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(2,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(2,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(nx,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(nx,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(nx,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            !
            !*********************back boundary with long characteristics***********
            !
          case(23)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            alocont_o_nn3d(nx,1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,1,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,ny,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,ny,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,ny,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,ny-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,ny-1,k,  18) = alocont_o_nn3d(i-1,j-1,k,  18)
            alocont_o_nn3d(nx,ny-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            !
            !**********************right boundary with long characteristics*********
            !
          case(25)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(1,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(1,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            alocont_o_nn3d(nx-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(nx,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(nx,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            !
            !**********************back boundary with long characteristics**********
            !
          case(26)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,1,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,1,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            alocont_o_nn3d(1,ny-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,ny-1,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,ny-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,ny,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,ny,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,ny,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            !
            !************************standard radiative transfer********************
            !
          case(3)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            !
            call fsc_cont3d_lin(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
          end select
          !
          !
        enddo
      enddo
    enddo
    !
    !---------------------perform angular integration-----------------------
    !NOTE: due to periodic boundary conditions, and copying of ALO coefficients
    !      angular integration cannot (easily) be performed within loop from above
    !
    do k=1, nz
      do j=1, ny
        do i=1, nx
          mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
          fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
          fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
          fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
          kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
          kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
          kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
          kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
          kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
          kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
          normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall

          alocont_nn3d_tmp(i,j,k,1) = alocont_nn3d_tmp(i,j,k,1) + wall*alocont_o_nn3d(i,j,k,1)
          alocont_nn3d_tmp(i,j,k,2) = alocont_nn3d_tmp(i,j,k,2) + wall*alocont_o_nn3d(i,j,k,2)
          alocont_nn3d_tmp(i,j,k,3) = alocont_nn3d_tmp(i,j,k,3) + wall*alocont_o_nn3d(i,j,k,3)
          alocont_nn3d_tmp(i,j,k,4) = alocont_nn3d_tmp(i,j,k,4) + wall*alocont_o_nn3d(i,j,k,4)
          alocont_nn3d_tmp(i,j,k,5) = alocont_nn3d_tmp(i,j,k,5) + wall*alocont_o_nn3d(i,j,k,5)
          alocont_nn3d_tmp(i,j,k,6) = alocont_nn3d_tmp(i,j,k,6) + wall*alocont_o_nn3d(i,j,k,6)
          alocont_nn3d_tmp(i,j,k,7) = alocont_nn3d_tmp(i,j,k,7) + wall*alocont_o_nn3d(i,j,k,7)
          alocont_nn3d_tmp(i,j,k,8) = alocont_nn3d_tmp(i,j,k,8) + wall*alocont_o_nn3d(i,j,k,8)
          alocont_nn3d_tmp(i,j,k,9) = alocont_nn3d_tmp(i,j,k,9) + wall*alocont_o_nn3d(i,j,k,9)
          alocont_nn3d_tmp(i,j,k,10) = alocont_nn3d_tmp(i,j,k,10) + wall*alocont_o_nn3d(i,j,k,10)
          alocont_nn3d_tmp(i,j,k,11) = alocont_nn3d_tmp(i,j,k,11) + wall*alocont_o_nn3d(i,j,k,11)
          alocont_nn3d_tmp(i,j,k,12) = alocont_nn3d_tmp(i,j,k,12) + wall*alocont_o_nn3d(i,j,k,12)
          alocont_nn3d_tmp(i,j,k,13) = alocont_nn3d_tmp(i,j,k,13) + wall*alocont_o_nn3d(i,j,k,13)
          alocont_nn3d_tmp(i,j,k,14) = alocont_nn3d_tmp(i,j,k,14) + wall*alocont_o_nn3d(i,j,k,14)
          alocont_nn3d_tmp(i,j,k,15) = alocont_nn3d_tmp(i,j,k,15) + wall*alocont_o_nn3d(i,j,k,15)
          alocont_nn3d_tmp(i,j,k,16) = alocont_nn3d_tmp(i,j,k,16) + wall*alocont_o_nn3d(i,j,k,16)
          alocont_nn3d_tmp(i,j,k,17) = alocont_nn3d_tmp(i,j,k,17) + wall*alocont_o_nn3d(i,j,k,17)
          alocont_nn3d_tmp(i,j,k,18) = alocont_nn3d_tmp(i,j,k,18) + wall*alocont_o_nn3d(i,j,k,18)
          alocont_nn3d_tmp(i,j,k,19) = alocont_nn3d_tmp(i,j,k,19) + wall*alocont_o_nn3d(i,j,k,19)
          alocont_nn3d_tmp(i,j,k,20) = alocont_nn3d_tmp(i,j,k,20) + wall*alocont_o_nn3d(i,j,k,20)
          alocont_nn3d_tmp(i,j,k,21) = alocont_nn3d_tmp(i,j,k,21) + wall*alocont_o_nn3d(i,j,k,21)
          alocont_nn3d_tmp(i,j,k,22) = alocont_nn3d_tmp(i,j,k,22) + wall*alocont_o_nn3d(i,j,k,22)
          alocont_nn3d_tmp(i,j,k,23) = alocont_nn3d_tmp(i,j,k,23) + wall*alocont_o_nn3d(i,j,k,23)
          alocont_nn3d_tmp(i,j,k,24) = alocont_nn3d_tmp(i,j,k,24) + wall*alocont_o_nn3d(i,j,k,24)
          alocont_nn3d_tmp(i,j,k,25) = alocont_nn3d_tmp(i,j,k,25) + wall*alocont_o_nn3d(i,j,k,25)
          alocont_nn3d_tmp(i,j,k,26) = alocont_nn3d_tmp(i,j,k,26) + wall*alocont_o_nn3d(i,j,k,26)
          alocont_nn3d_tmp(i,j,k,27) = alocont_nn3d_tmp(i,j,k,27) + wall*alocont_o_nn3d(i,j,k,27)
        enddo
      enddo
    enddo
    !
  end subroutine formalsc_cont3d_lin
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine formallc_cont3d(oindx,nueindx)
    !
    !-----------------------------------------------------------------------
    !------long characteristics for continuum radiative transfer in 3d------
    !-----------calculating intensties for given mu,phi specified-----------
    !---------------------------by input oindx------------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, int3d, opac3d, scont3d, imaskb3d, &
         alocont_o_nn3d, alocont_nn3d_tmp, mint3d_tmp, normalization3d_tmp, &
         fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, bnue3d, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp
    use mod_angles, only: n_x, n_y, n_z, weight_omega, q_alo
    use mod_frequencies, only: nodes_nue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: oindx, nueindx
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: alpha, beta, gamma
    integer(i4b) :: startx, starty, startz, endx, endy, endz, &
    startxb, startyb, endxb, endyb
    integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
    integer(i4b) :: aloindx_im1, aloindx_ip1, aloindx_jm1, aloindx_jp1, aloindx_km1, aloindx_kp1
    real(dp) :: nn_x, nn_y, nn_z, xnue, wall
    real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
    real(dp) :: dels_xyd, dels_xzd, dels_yzd, dels_d
    integer :: q1, q2, q3, q4, q5, q6, q7, q8, q9, &
    q10, q11, q12, q13, q14, q15, q16, q17, q18, &
    q19, q20, q21, q22, q23, q24, q25, q26, q27
    !
    !for debugging
    real(dp) :: int_u2, scont_u2, opac_u2, int_u3, scont_u3, opac_u3, scont_d2, scont_d3, opac_d2, opac_d3, interpol2d_9p_quad, interpol2d_9p_bez
    !
    ! ... local arrays
    !
    ! ... local functions
    real(dp) :: dist_ellipsoid, calc_icore_gdark
    !
    ! ... local characters
    !character(len=50) :: enter
    !
    ! ... local logicals
    !
    !frequency
    xnue=nodes_nue(nueindx)
    !
    !directions
    nn_x=n_x(oindx)
    nn_y=n_y(oindx)
    nn_z=n_z(oindx)
    !
    !angulare integration weight
    wall=weight_omega(oindx)
    !
    !indices for nearest neighbour alo
    q1=q_alo(oindx,1)
    q2=q_alo(oindx,2)
    q3=q_alo(oindx,3)
    q4=q_alo(oindx,4)
    q5=q_alo(oindx,5)
    q6=q_alo(oindx,6)
    q7=q_alo(oindx,7)
    q8=q_alo(oindx,8)
    q9=q_alo(oindx,9)
    q10=q_alo(oindx,10)
    q11=q_alo(oindx,11)
    q12=q_alo(oindx,12)
    q13=q_alo(oindx,13)
    q14=q_alo(oindx,14)
    q15=q_alo(oindx,15)
    q16=q_alo(oindx,16)
    q17=q_alo(oindx,17)
    q18=q_alo(oindx,18)
    q19=q_alo(oindx,19)
    q20=q_alo(oindx,20)
    q21=q_alo(oindx,21)
    q22=q_alo(oindx,22)
    q23=q_alo(oindx,23)
    q24=q_alo(oindx,24)
    q25=q_alo(oindx,25)
    q26=q_alo(oindx,26)
    q27=q_alo(oindx,27)
    !
    !-----------------------------------------------------------------------
    !
    !set directional index-parameter (phantom points are excluded from calculation)
    !
    !index parameter:
    !         if n_x >= 0                 if n_x < 0
    !                startx = 2                  startx = ndxmax-1
    !                endx = ndxmax-1             endx = 2
    !                alpha=  1                   alpha=-1
    !
    !         if n_y >= 0                 if n_y < 0
    !                starty = 2                  starty = ndymax-1
    !                endy = ndymax-1             endy = 2
    !                beta =  1                   beta =-1
    !
    !         if n_z >= 0                 if n_z < 0
    !                startz = 2                  startz = ndzmax-1
    !                endz = ndzmax-1             endz = 2
    !                gamma = 1                   gamma = -1
    !
    if(nn_x.gt.zero) then
      startx = 1
      endx = nx
      startxb = nx
      endxb = 1
      alpha=  1
      aloindx_im1 = nx-1
      aloindx_ip1 = 2
    elseif(nn_x.lt.zero) then
      startx = nx
      endx = 1
      startxb = 1
      endxb = nx
      alpha=-1
      aloindx_im1 = 2
      aloindx_ip1 = nx-1
    else
      stop 'error in formallc_cont3d: n_x = 0 not allowed'
    endif
    !
    if(nn_y.gt.zero) then
      starty = 1
      endy = ny
      startyb = ny
      endyb = 1
      beta =  1
      aloindx_jm1 = ny-1
      aloindx_jp1 = 2
    elseif(nn_y.lt.zero) then
      starty = ny
      endy = 1
      startyb = 1
      endyb = ny
      beta =-1
      aloindx_jm1 = 2
      aloindx_jp1 = ny-1
    else
      stop 'error in formallc_cont3d: n_y = 0 not allowed'
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
    else
      stop 'error in formallc_cont3d: n_z = 0 not allowed'
    endif
    !
    !--------------------reset the intensities and alo----------------------
    !
    call set_boundary3d(xnue, nn_x, nn_y, nn_z)
    !
    !-----------------------------------------------------------------------
    !
    alocont_o_nn3d=zero
    !
    !-----------------------------------------------------------------------
    !
    do k=startz, endz, gamma
      do j=starty, endy, beta
        do i=startx, endx, alpha
          !
          select case(imaskb3d(i,j,k))
            !
            !*************left and right boundary with long characteristics*********
            !
          case(12,13)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
            !*************front and back boundary with long characteristics*********
            !
          case(14,15)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
            !********(left,front), (right,front), (left,back), (right,back)*********
            !*******************boundary edges with long characteristics************
            !
          case(18,21,24,27)
            !
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
            !***********************adjacent to left boundary***********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(4)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on left towards right
            alocont_o_nn3d(nx,j-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,j-1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,j-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,j+1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx,j+1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx,j+1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,j  ,k+1, 6) = alocont_o_nn3d(i-1,j  ,k+1, 6)
            alocont_o_nn3d(nx,j  ,k  ,15) = alocont_o_nn3d(i-1,j  ,k  ,15)
            alocont_o_nn3d(nx,j  ,k-1,24) = alocont_o_nn3d(i-1,j  ,k-1,24)
            !
            !               alocont_nn3d_tmp(nx,j-1,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !               alocont_nn3d_tmp(nx,j-1,k  ,18) = alocont_nn3d_tmp(i-1,j-1,k  ,18)
            !               alocont_nn3d_tmp(nx,j-1,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !               alocont_nn3d_tmp(nx,j+1,k+1, 3) = alocont_nn3d_tmp(i-1,j+1,k+1, 3)
            !               alocont_nn3d_tmp(nx,j+1,k  ,12) = alocont_nn3d_tmp(i-1,j+1,k  ,12)
            !               alocont_nn3d_tmp(nx,j+1,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !               alocont_nn3d_tmp(nx,j  ,k+1, 6) = alocont_nn3d_tmp(i-1,j  ,k+1, 6)
            !               alocont_nn3d_tmp(nx,j  ,k  ,15) = alocont_nn3d_tmp(i-1,j  ,k  ,15)
            !               alocont_nn3d_tmp(nx,j  ,k-1,24) = alocont_nn3d_tmp(i-1,j  ,k-1,24)
            !
            !***********************adjacent to right boundary**********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(5)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on right towards left
            alocont_o_nn3d(1,j-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,j-1,k,  16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,j-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,j+1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,j+1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,j+1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,j  ,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,j  ,k,  13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,j  ,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)

            !               alocont_nn3d_tmp(1,j-1,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !               alocont_nn3d_tmp(1,j-1,k,  16) = alocont_nn3d_tmp(i+1,j-1,k  ,16)
            !               alocont_nn3d_tmp(1,j-1,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !               alocont_nn3d_tmp(1,j+1,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !               alocont_nn3d_tmp(1,j+1,k  ,10) = alocont_nn3d_tmp(i+1,j+1,k  ,10)
            !               alocont_nn3d_tmp(1,j+1,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !               alocont_nn3d_tmp(1,j  ,k+1, 4) = alocont_nn3d_tmp(i+1,j  ,k+1, 4)
            !               alocont_nn3d_tmp(1,j  ,k,  13) = alocont_nn3d_tmp(i+1,j  ,k  ,13)
            !               alocont_nn3d_tmp(1,j  ,k-1,22) = alocont_nn3d_tmp(i+1,j  ,k-1,22)
            !
            !***********************adjacent to front boundary**********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(6)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards front
            alocont_o_nn3d(i-1,ny,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(i-1,ny,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(i-1,ny,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(i  ,ny,k+1, 8) = alocont_o_nn3d(i  ,j-1,k+1, 8)
            alocont_o_nn3d(i  ,ny,k  ,17) = alocont_o_nn3d(i  ,j-1,k  ,17)
            alocont_o_nn3d(i  ,ny,k-1,26) = alocont_o_nn3d(i  ,j-1,k-1,26)
            alocont_o_nn3d(i+1,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(i+1,ny,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(i+1,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)

            !               alocont_nn3d_tmp(i-1,ny,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !               alocont_nn3d_tmp(i-1,ny,k  ,18) = alocont_nn3d_tmp(i-1,j-1,k  ,18)
            !               alocont_nn3d_tmp(i-1,ny,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !               alocont_nn3d_tmp(i  ,ny,k+1, 8) = alocont_nn3d_tmp(i  ,j-1,k+1, 8)
            !               alocont_nn3d_tmp(i  ,ny,k  ,17) = alocont_nn3d_tmp(i  ,j-1,k  ,17)
            !               alocont_nn3d_tmp(i  ,ny,k-1,26) = alocont_nn3d_tmp(i  ,j-1,k-1,26)
            !               alocont_nn3d_tmp(i+1,ny,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !               alocont_nn3d_tmp(i+1,ny,k  ,16) = alocont_nn3d_tmp(i+1,j-1,k  ,16)
            !               alocont_nn3d_tmp(i+1,ny,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !
            !***********************adjacent to back boundary***********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(7)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards front
            alocont_o_nn3d(i+1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(i+1,1,k,  10) = alocont_o_nn3d(i+1,j+1,k,  10)
            alocont_o_nn3d(i+1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(i  ,1,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(i  ,1,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(i  ,1,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(i-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(i-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(i-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)

            !              alocont_nn3d_tmp(i+1,1,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !              alocont_nn3d_tmp(i+1,1,k,  10) = alocont_nn3d_tmp(i+1,j+1,k,  10)
            !              alocont_nn3d_tmp(i+1,1,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !              alocont_nn3d_tmp(i  ,1,k+1, 2) = alocont_nn3d_tmp(i  ,j+1,k+1, 2)
            !              alocont_nn3d_tmp(i  ,1,k  ,11) = alocont_nn3d_tmp(i  ,j+1,k  ,11)
            !              alocont_nn3d_tmp(i  ,1,k-1,20) = alocont_nn3d_tmp(i  ,j+1,k-1,20)
            !              alocont_nn3d_tmp(i-1,1,k+1, 3) = alocont_nn3d_tmp(i-1,j+1,k+1, 3)
            !              alocont_nn3d_tmp(i-1,1,k  ,12) = alocont_nn3d_tmp(i-1,j+1,k  ,12)
            !              alocont_nn3d_tmp(i-1,1,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !
            !***********************adjacent to left front boundary*****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(8)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(nx,1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,2,k+1, 6) = alocont_o_nn3d(i-1,j  ,k+1, 6)
            alocont_o_nn3d(nx,2,k  ,15) = alocont_o_nn3d(i-1,j  ,k  ,15)
            alocont_o_nn3d(nx,2,k-1,24) = alocont_o_nn3d(i-1,j  ,k-1,24)
            alocont_o_nn3d(nx,j+1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx,j+1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx,j+1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(1,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(1,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(1,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(2,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(2,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(2,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(3,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(3,ny,k,  16) = alocont_o_nn3d(i+1,j-1,k,  16)
            alocont_o_nn3d(3,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx ,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx ,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx ,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)

            !               alocont_nn3d_tmp(nx,1,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !               alocont_nn3d_tmp(nx,1,k  ,18) = alocont_nn3d_tmp(i-1,j-1,k  ,18)
            !               alocont_nn3d_tmp(nx,1,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !               alocont_nn3d_tmp(nx,2,k+1, 6) = alocont_nn3d_tmp(i-1,j  ,k+1, 6)
            !               alocont_nn3d_tmp(nx,2,k  ,15) = alocont_nn3d_tmp(i-1,j  ,k  ,15)
            !               alocont_nn3d_tmp(nx,2,k-1,24) = alocont_nn3d_tmp(i-1,j  ,k-1,24)
            !               alocont_nn3d_tmp(nx,j+1,k+1, 3) = alocont_nn3d_tmp(i-1,j+1,k+1, 3)
            !               alocont_nn3d_tmp(nx,j+1,k  ,12) = alocont_nn3d_tmp(i-1,j+1,k  ,12)
            !               alocont_nn3d_tmp(nx,j+1,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !               alocont_nn3d_tmp(1,ny ,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !               alocont_nn3d_tmp(1,ny ,k  ,18) = alocont_nn3d_tmp(i-1,j-1,k  ,18)
            !               alocont_nn3d_tmp(1,ny ,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !               alocont_nn3d_tmp(2,ny,k+1, 8) = alocont_nn3d_tmp(i,j-1,k+1, 8)
            !               alocont_nn3d_tmp(2,ny,k,  17) = alocont_nn3d_tmp(i,j-1,k,  17)
            !               alocont_nn3d_tmp(2,ny,k-1,26) = alocont_nn3d_tmp(i,j-1,k-1,26)
            !               alocont_nn3d_tmp(3,ny,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !               alocont_nn3d_tmp(3,ny,k,  16) = alocont_nn3d_tmp(i+1,j-1,k,  16)
            !               alocont_nn3d_tmp(3,ny,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !               alocont_nn3d_tmp(nx ,ny ,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !               alocont_nn3d_tmp(nx ,ny ,k  ,18) = alocont_nn3d_tmp(i-1,j-1,k  ,18)
            !               alocont_nn3d_tmp(nx ,ny ,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !
            !***********************adjacent to right front boundary****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(9)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)

            alocont_o_nn3d(1,2,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,2,k  ,13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,2,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            alocont_o_nn3d(1,j+1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,j+1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,j+1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,1,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,ny,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx,ny ,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(nx,ny ,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(nx,ny ,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx-2,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx-2,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx-2,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx-1,ny ,k+1, 8) = alocont_o_nn3d(i  ,j-1,k+1, 8)
            alocont_o_nn3d(nx-1,ny ,k  ,17) = alocont_o_nn3d(i  ,j-1,k  ,17)
            alocont_o_nn3d(nx-1,ny ,k-1,26) = alocont_o_nn3d(i  ,j-1,k-1,26)

            !              alocont_nn3d_tmp(1,2,k+1, 4) = alocont_nn3d_tmp(i+1,j  ,k+1, 4)
            !              alocont_nn3d_tmp(1,2,k  ,13) = alocont_nn3d_tmp(i+1,j  ,k  ,13)
            !              alocont_nn3d_tmp(1,2,k-1,22) = alocont_nn3d_tmp(i+1,j  ,k-1,22)
            !              alocont_nn3d_tmp(1,j+1,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !              alocont_nn3d_tmp(1,j+1,k  ,10) = alocont_nn3d_tmp(i+1,j+1,k  ,10)
            !              alocont_nn3d_tmp(1,j+1,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !              alocont_nn3d_tmp(1,1,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !              alocont_nn3d_tmp(1,1,k  ,16) = alocont_nn3d_tmp(i+1,j-1,k  ,16)
            !              alocont_nn3d_tmp(1,1,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !              alocont_nn3d_tmp(1,ny,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !              alocont_nn3d_tmp(1,ny,k  ,16) = alocont_nn3d_tmp(i+1,j-1,k  ,16)
            !              alocont_nn3d_tmp(1,ny,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !              alocont_nn3d_tmp(nx,ny ,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !              alocont_nn3d_tmp(nx,ny ,k  ,16) = alocont_nn3d_tmp(i+1,j-1,k  ,16)
            !              alocont_nn3d_tmp(nx,ny ,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !              alocont_nn3d_tmp(nx-2,ny ,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !              alocont_nn3d_tmp(nx-2,ny ,k  ,18) = alocont_nn3d_tmp(i-1,j-1,k  ,18)
            !              alocont_nn3d_tmp(nx-2,ny ,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !              alocont_nn3d_tmp(nx-1,ny ,k+1, 8) = alocont_nn3d_tmp(i  ,j-1,k+1, 8)
            !              alocont_nn3d_tmp(nx-1,ny ,k  ,17) = alocont_nn3d_tmp(i  ,j-1,k  ,17)
            !              alocont_nn3d_tmp(nx-1,ny ,k-1,26) = alocont_nn3d_tmp(i  ,j-1,k-1,26)

            !
            !***********************adjacent to left back boundary******************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(10)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1  ,1  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(1  ,1  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(1  ,1  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(2  ,1  ,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(2  ,1  ,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(2  ,1  ,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(nx  ,1  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx  ,1  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx  ,1  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,j-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,j-1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,j-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx  ,ny  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx  ,ny  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx  ,ny  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx ,ny-1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx ,ny-1,k,  15) = alocont_o_nn3d(i-1,j,k,  15)
            alocont_o_nn3d(nx ,ny-1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(i+1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(i+1,1,k,  10) = alocont_o_nn3d(i+1,j+1,k,  10)
            alocont_o_nn3d(i+1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)

            !               alocont_nn3d_tmp(1  ,1  ,k+1, 3) = alocont_nn3d_tmp(i-1,j+1,k+1, 3)
            !               alocont_nn3d_tmp(1  ,1  ,k  ,12) = alocont_nn3d_tmp(i-1,j+1,k  ,12)
            !               alocont_nn3d_tmp(1  ,1  ,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !               alocont_nn3d_tmp(2  ,1  ,k+1, 2) = alocont_nn3d_tmp(i  ,j+1,k+1, 2)
            !               alocont_nn3d_tmp(2  ,1  ,k  ,11) = alocont_nn3d_tmp(i  ,j+1,k  ,11)
            !               alocont_nn3d_tmp(2  ,1  ,k-1,20) = alocont_nn3d_tmp(i  ,j+1,k-1,20)
            !               alocont_nn3d_tmp(nx  ,1  ,k+1, 3) = alocont_nn3d_tmp(i-1,j+1,k+1, 3)
            !               alocont_nn3d_tmp(nx  ,1  ,k  ,12) = alocont_nn3d_tmp(i-1,j+1,k  ,12)
            !               alocont_nn3d_tmp(nx  ,1  ,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !               alocont_nn3d_tmp(nx,j-1,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !               alocont_nn3d_tmp(nx,j-1,k  ,18) = alocont_nn3d_tmp(i-1,j-1,k  ,18)
            !               alocont_nn3d_tmp(nx,j-1,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !               alocont_nn3d_tmp(nx  ,ny  ,k+1, 3) = alocont_nn3d_tmp(i-1,j+1,k+1, 3)
            !               alocont_nn3d_tmp(nx  ,ny  ,k  ,12) = alocont_nn3d_tmp(i-1,j+1,k  ,12)
            !               alocont_nn3d_tmp(nx  ,ny  ,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !               alocont_nn3d_tmp(nx ,ny-1,k+1, 6) = alocont_nn3d_tmp(i-1,j,k+1, 6)
            !               alocont_nn3d_tmp(nx ,ny-1,k,  15) = alocont_nn3d_tmp(i-1,j,k,  15)
            !               alocont_nn3d_tmp(nx ,ny-1,k-1,24) = alocont_nn3d_tmp(i-1,j,k-1,24)
            !               alocont_nn3d_tmp(i+1,1,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !               alocont_nn3d_tmp(i+1,1,k,  10) = alocont_nn3d_tmp(i+1,j+1,k,  10)
            !               alocont_nn3d_tmp(i+1,1,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !
            !***********************adjacent to right back boundary*****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(11)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(nx-1,1,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(nx-1,1,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(nx-1,1,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(nx,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(nx,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(nx,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,j-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,j-1,k,  16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,j-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,j  ,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,j  ,k  ,13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,j  ,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            alocont_o_nn3d(1,ny,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,ny,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,ny,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(i-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(i-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(i-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)

            !               alocont_nn3d_tmp(1,1,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !               alocont_nn3d_tmp(1,1,k  ,10) = alocont_nn3d_tmp(i+1,j+1,k  ,10)
            !               alocont_nn3d_tmp(1,1,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !               alocont_nn3d_tmp(nx-1,1,k+1, 2) = alocont_nn3d_tmp(i  ,j+1,k+1, 2)
            !               alocont_nn3d_tmp(nx-1,1,k  ,11) = alocont_nn3d_tmp(i  ,j+1,k  ,11)
            !               alocont_nn3d_tmp(nx-1,1,k-1,20) = alocont_nn3d_tmp(i  ,j+1,k-1,20)
            !               alocont_nn3d_tmp(nx,1,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !               alocont_nn3d_tmp(nx,1,k  ,10) = alocont_nn3d_tmp(i+1,j+1,k  ,10)
            !               alocont_nn3d_tmp(nx,1,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !               alocont_nn3d_tmp(1,j-1,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !               alocont_nn3d_tmp(1,j-1,k,  16) = alocont_nn3d_tmp(i+1,j-1,k  ,16)
            !               alocont_nn3d_tmp(1,j-1,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !               alocont_nn3d_tmp(1,j  ,k+1, 4) = alocont_nn3d_tmp(i+1,j  ,k+1, 4)
            !               alocont_nn3d_tmp(1,j  ,k  ,13) = alocont_nn3d_tmp(i+1,j  ,k  ,13)
            !               alocont_nn3d_tmp(1,j  ,k-1,22) = alocont_nn3d_tmp(i+1,j  ,k-1,22)
            !               alocont_nn3d_tmp(1,ny,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !               alocont_nn3d_tmp(1,ny,k  ,10) = alocont_nn3d_tmp(i+1,j+1,k  ,10)
            !               alocont_nn3d_tmp(1,ny,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !               alocont_nn3d_tmp(i-1,1,k+1, 3) = alocont_nn3d_tmp(i-1,j+1,k+1, 3)
            !               alocont_nn3d_tmp(i-1,1,k  ,12) = alocont_nn3d_tmp(i-1,j+1,k  ,12)
            !               alocont_nn3d_tmp(i-1,1,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !
            !**********************left boundary with long characteristics**********
            !
          case(16)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !
            alocont_o_nn3d(2,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(2,ny,k,  16) = alocont_o_nn3d(i+1,j-1,k,  16)
            alocont_o_nn3d(2,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(nx,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(nx,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(1,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(1,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(1,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)

            !               alocont_nn3d_tmp(2,ny,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !               alocont_nn3d_tmp(2,ny,k,  16) = alocont_nn3d_tmp(i+1,j-1,k,  16)
            !               alocont_nn3d_tmp(2,ny,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !               alocont_nn3d_tmp(nx,ny,k+1, 8) = alocont_nn3d_tmp(i,j-1,k+1, 8)
            !               alocont_nn3d_tmp(nx,ny,k,  17) = alocont_nn3d_tmp(i,j-1,k,  17)
            !               alocont_nn3d_tmp(nx,ny,k-1,26) = alocont_nn3d_tmp(i,j-1,k-1,26)
            !               alocont_nn3d_tmp(1,ny,k+1, 8) = alocont_nn3d_tmp(i,j-1,k+1, 8)
            !               alocont_nn3d_tmp(1,ny,k,  17) = alocont_nn3d_tmp(i,j-1,k,  17)
            !               alocont_nn3d_tmp(1,ny,k-1,26) = alocont_nn3d_tmp(i,j-1,k-1,26)
            !
            !*********************front boundary with long characteristics**********
            !
          case(17)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            alocont_o_nn3d(nx,1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,1,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,2,k+1,3) = alocont_o_nn3d(i-1,j+1,k+1,3)
            alocont_o_nn3d(nx,2,k,12) = alocont_o_nn3d(i-1,j+1,k,12)
            alocont_o_nn3d(nx,2,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,ny,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,ny,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,ny,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)

            !               alocont_nn3d_tmp(nx,1,k+1, 6) = alocont_nn3d_tmp(i-1,j,k+1, 6)
            !               alocont_nn3d_tmp(nx,1,k,15) = alocont_nn3d_tmp(i-1,j,k,15)
            !               alocont_nn3d_tmp(nx,1,k-1,24) = alocont_nn3d_tmp(i-1,j,k-1,24)
            !               alocont_nn3d_tmp(nx,2,k+1,3) = alocont_nn3d_tmp(i-1,j+1,k+1,3)
            !               alocont_nn3d_tmp(nx,2,k,12) = alocont_nn3d_tmp(i-1,j+1,k,12)
            !               alocont_nn3d_tmp(nx,2,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !               alocont_nn3d_tmp(nx,ny,k+1, 6) = alocont_nn3d_tmp(i-1,j,k+1, 6)
            !               alocont_nn3d_tmp(nx,ny,k,15) = alocont_nn3d_tmp(i-1,j,k,15)
            !               alocont_nn3d_tmp(nx,ny,k-1,24) = alocont_nn3d_tmp(i-1,j,k-1,24)
            !
            !*********************right boundary with long characteristics**********
            !
          case(19)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            alocont_o_nn3d(nx-1,ny,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx-1,ny,k,  18) = alocont_o_nn3d(i-1,j-1,k,  18)
            alocont_o_nn3d(nx-1,ny,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(nx,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(nx,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(1,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(1,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(1,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)

            !               alocont_nn3d_tmp(nx-1,ny,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !               alocont_nn3d_tmp(nx-1,ny,k,  18) = alocont_nn3d_tmp(i-1,j-1,k,  18)
            !               alocont_nn3d_tmp(nx-1,ny,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !               alocont_nn3d_tmp(nx,ny,k+1, 8) = alocont_nn3d_tmp(i,j-1,k+1, 8)
            !               alocont_nn3d_tmp(nx,ny,k,  17) = alocont_nn3d_tmp(i,j-1,k,  17)
            !               alocont_nn3d_tmp(nx,ny,k-1,26) = alocont_nn3d_tmp(i,j-1,k-1,26)
            !               alocont_nn3d_tmp(1,ny,k+1, 8) = alocont_nn3d_tmp(i,j-1,k+1, 8)
            !               alocont_nn3d_tmp(1,ny,k,  17) = alocont_nn3d_tmp(i,j-1,k,  17)
            !               alocont_nn3d_tmp(1,ny,k-1,26) = alocont_nn3d_tmp(i,j-1,k-1,26)
            !
            !*********************front boundary with long characteristics**********
            !
          case(20)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,1,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,1,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            alocont_o_nn3d(1,2,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,2,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,2,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,ny,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,ny,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,ny,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)

            !               alocont_nn3d_tmp(1,1,k+1, 4) = alocont_nn3d_tmp(i+1,j,k+1, 4)
            !               alocont_nn3d_tmp(1,1,k  ,13) = alocont_nn3d_tmp(i+1,j,k  ,13)
            !               alocont_nn3d_tmp(1,1,k-1,22) = alocont_nn3d_tmp(i+1,j,k-1,22)
            !               alocont_nn3d_tmp(1,2,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !               alocont_nn3d_tmp(1,2,k  ,10) = alocont_nn3d_tmp(i+1,j+1,k  ,10)
            !               alocont_nn3d_tmp(1,2,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !               alocont_nn3d_tmp(1,ny,k+1, 4) = alocont_nn3d_tmp(i+1,j,k+1, 4)
            !               alocont_nn3d_tmp(1,ny,k  ,13) = alocont_nn3d_tmp(i+1,j,k  ,13)
            !               alocont_nn3d_tmp(1,ny,k-1,22) = alocont_nn3d_tmp(i+1,j,k-1,22)
            !
            !**********************left boundary with long characteristics**********
            !
          case(22)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(1,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(1,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            alocont_o_nn3d(2,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(2,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(2,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(nx,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(nx,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(nx,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)

            !               alocont_nn3d_tmp(1,1,k+1, 2) = alocont_nn3d_tmp(i,j+1,k+1, 2)
            !               alocont_nn3d_tmp(1,1,k  ,11) = alocont_nn3d_tmp(i,j+1,k  ,11)
            !               alocont_nn3d_tmp(1,1,k-1,20) = alocont_nn3d_tmp(i,j+1,k-1,20)
            !               alocont_nn3d_tmp(2,1,k+1, 1) = alocont_nn3d_tmp(i+1,j+1,k+1, 1)
            !               alocont_nn3d_tmp(2,1,k  ,10) = alocont_nn3d_tmp(i+1,j+1,k  ,10)
            !               alocont_nn3d_tmp(2,1,k-1,19) = alocont_nn3d_tmp(i+1,j+1,k-1,19)
            !               alocont_nn3d_tmp(nx,1,k+1, 2) = alocont_nn3d_tmp(i,j+1,k+1, 2)
            !               alocont_nn3d_tmp(nx,1,k  ,11) = alocont_nn3d_tmp(i,j+1,k  ,11)
            !               alocont_nn3d_tmp(nx,1,k-1,20) = alocont_nn3d_tmp(i,j+1,k-1,20)
            !
            !*********************back boundary with long characteristics***********
            !
          case(23)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            alocont_o_nn3d(nx,1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,1,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,ny,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,ny,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,ny,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,ny-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,ny-1,k,  18) = alocont_o_nn3d(i-1,j-1,k,  18)
            alocont_o_nn3d(nx,ny-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)

            !               alocont_nn3d_tmp(nx,1,k+1, 6) = alocont_nn3d_tmp(i-1,j,k+1, 6)
            !               alocont_nn3d_tmp(nx,1,k,15) = alocont_nn3d_tmp(i-1,j,k,15)
            !               alocont_nn3d_tmp(nx,1,k-1,24) = alocont_nn3d_tmp(i-1,j,k-1,24)
            !               alocont_nn3d_tmp(nx,ny,k+1, 6) = alocont_nn3d_tmp(i-1,j,k+1, 6)
            !               alocont_nn3d_tmp(nx,ny,k,15) = alocont_nn3d_tmp(i-1,j,k,15)
            !               alocont_nn3d_tmp(nx,ny,k-1,24) = alocont_nn3d_tmp(i-1,j,k-1,24)
            !               alocont_nn3d_tmp(nx,ny-1,k+1, 9) = alocont_nn3d_tmp(i-1,j-1,k+1, 9)
            !               alocont_nn3d_tmp(nx,ny-1,k,  18) = alocont_nn3d_tmp(i-1,j-1,k,  18)
            !               alocont_nn3d_tmp(nx,ny-1,k-1,27) = alocont_nn3d_tmp(i-1,j-1,k-1,27)
            !
            !**********************right boundary with long characteristics*********
            !
          case(25)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(1,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(1,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            alocont_o_nn3d(nx-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(nx,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(nx,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)

            !               alocont_nn3d_tmp(1,1,k+1, 2) = alocont_nn3d_tmp(i,j+1,k+1, 2)
            !               alocont_nn3d_tmp(1,1,k  ,11) = alocont_nn3d_tmp(i,j+1,k  ,11)
            !               alocont_nn3d_tmp(1,1,k-1,20) = alocont_nn3d_tmp(i,j+1,k-1,20)
            !               alocont_nn3d_tmp(nx-1,1,k+1, 3) = alocont_nn3d_tmp(i-1,j+1,k+1, 3)
            !               alocont_nn3d_tmp(nx-1,1,k  ,12) = alocont_nn3d_tmp(i-1,j+1,k  ,12)
            !               alocont_nn3d_tmp(nx-1,1,k-1,21) = alocont_nn3d_tmp(i-1,j+1,k-1,21)
            !               alocont_nn3d_tmp(nx,1,k+1, 2) = alocont_nn3d_tmp(i,j+1,k+1, 2)
            !               alocont_nn3d_tmp(nx,1,k  ,11) = alocont_nn3d_tmp(i,j+1,k  ,11)
            !               alocont_nn3d_tmp(nx,1,k-1,20) = alocont_nn3d_tmp(i,j+1,k-1,20)
            !
            !**********************back boundary with long characteristics**********
            !
          case(26)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,1,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,1,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            alocont_o_nn3d(1,ny-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,ny-1,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,ny-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,ny,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,ny,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,ny,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)

            !               alocont_nn3d_tmp(1,1,k+1, 4) = alocont_nn3d_tmp(i+1,j,k+1, 4)
            !               alocont_nn3d_tmp(1,1,k  ,13) = alocont_nn3d_tmp(i+1,j,k  ,13)
            !               alocont_nn3d_tmp(1,1,k-1,22) = alocont_nn3d_tmp(i+1,j,k-1,22)
            !               alocont_nn3d_tmp(1,ny-1,k+1, 7) = alocont_nn3d_tmp(i+1,j-1,k+1, 7)
            !               alocont_nn3d_tmp(1,ny-1,k  ,16) = alocont_nn3d_tmp(i+1,j-1,k  ,16)
            !               alocont_nn3d_tmp(1,ny-1,k-1,25) = alocont_nn3d_tmp(i+1,j-1,k-1,25)
            !               alocont_nn3d_tmp(1,ny,k+1, 4) = alocont_nn3d_tmp(i+1,j,k+1, 4)
            !               alocont_nn3d_tmp(1,ny,k  ,13) = alocont_nn3d_tmp(i+1,j,k  ,13)
            !               alocont_nn3d_tmp(1,ny,k-1,22) = alocont_nn3d_tmp(i+1,j,k-1,22)
            !
            !************************standard radiative transfer********************
            !
          case(3)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            q25, q26, q27)

            !               call fsc_cont3d(iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            !                               nn_x, nn_y, nn_z, wall, &
            !                               q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
            !                               q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
            !                               q25, q26, q27)
            !
          end select
          !
          !
        enddo
      enddo
    enddo
    !
    !
    !---------------------perform angular integration-----------------------
    !NOTE: due to periodic boundary conditions, and copying of ALO coefficients
    !      angular integration cannot (easily) be performed within loop from above
    !
    do k=1, nz
      do j=1, ny
        do i=1, nx
          mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
          fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
          fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
          fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
          kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
          kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
          kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
          kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
          kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
          kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
          normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall

          alocont_nn3d_tmp(i,j,k,1) = alocont_nn3d_tmp(i,j,k,1) + wall*alocont_o_nn3d(i,j,k,1)
          alocont_nn3d_tmp(i,j,k,2) = alocont_nn3d_tmp(i,j,k,2) + wall*alocont_o_nn3d(i,j,k,2)
          alocont_nn3d_tmp(i,j,k,3) = alocont_nn3d_tmp(i,j,k,3) + wall*alocont_o_nn3d(i,j,k,3)
          alocont_nn3d_tmp(i,j,k,4) = alocont_nn3d_tmp(i,j,k,4) + wall*alocont_o_nn3d(i,j,k,4)
          alocont_nn3d_tmp(i,j,k,5) = alocont_nn3d_tmp(i,j,k,5) + wall*alocont_o_nn3d(i,j,k,5)
          alocont_nn3d_tmp(i,j,k,6) = alocont_nn3d_tmp(i,j,k,6) + wall*alocont_o_nn3d(i,j,k,6)
          alocont_nn3d_tmp(i,j,k,7) = alocont_nn3d_tmp(i,j,k,7) + wall*alocont_o_nn3d(i,j,k,7)
          alocont_nn3d_tmp(i,j,k,8) = alocont_nn3d_tmp(i,j,k,8) + wall*alocont_o_nn3d(i,j,k,8)
          alocont_nn3d_tmp(i,j,k,9) = alocont_nn3d_tmp(i,j,k,9) + wall*alocont_o_nn3d(i,j,k,9)
          alocont_nn3d_tmp(i,j,k,10) = alocont_nn3d_tmp(i,j,k,10) + wall*alocont_o_nn3d(i,j,k,10)
          alocont_nn3d_tmp(i,j,k,11) = alocont_nn3d_tmp(i,j,k,11) + wall*alocont_o_nn3d(i,j,k,11)
          alocont_nn3d_tmp(i,j,k,12) = alocont_nn3d_tmp(i,j,k,12) + wall*alocont_o_nn3d(i,j,k,12)
          alocont_nn3d_tmp(i,j,k,13) = alocont_nn3d_tmp(i,j,k,13) + wall*alocont_o_nn3d(i,j,k,13)
          alocont_nn3d_tmp(i,j,k,14) = alocont_nn3d_tmp(i,j,k,14) + wall*alocont_o_nn3d(i,j,k,14)
          alocont_nn3d_tmp(i,j,k,15) = alocont_nn3d_tmp(i,j,k,15) + wall*alocont_o_nn3d(i,j,k,15)
          alocont_nn3d_tmp(i,j,k,16) = alocont_nn3d_tmp(i,j,k,16) + wall*alocont_o_nn3d(i,j,k,16)
          alocont_nn3d_tmp(i,j,k,17) = alocont_nn3d_tmp(i,j,k,17) + wall*alocont_o_nn3d(i,j,k,17)
          alocont_nn3d_tmp(i,j,k,18) = alocont_nn3d_tmp(i,j,k,18) + wall*alocont_o_nn3d(i,j,k,18)
          alocont_nn3d_tmp(i,j,k,19) = alocont_nn3d_tmp(i,j,k,19) + wall*alocont_o_nn3d(i,j,k,19)
          alocont_nn3d_tmp(i,j,k,20) = alocont_nn3d_tmp(i,j,k,20) + wall*alocont_o_nn3d(i,j,k,20)
          alocont_nn3d_tmp(i,j,k,21) = alocont_nn3d_tmp(i,j,k,21) + wall*alocont_o_nn3d(i,j,k,21)
          alocont_nn3d_tmp(i,j,k,22) = alocont_nn3d_tmp(i,j,k,22) + wall*alocont_o_nn3d(i,j,k,22)
          alocont_nn3d_tmp(i,j,k,23) = alocont_nn3d_tmp(i,j,k,23) + wall*alocont_o_nn3d(i,j,k,23)
          alocont_nn3d_tmp(i,j,k,24) = alocont_nn3d_tmp(i,j,k,24) + wall*alocont_o_nn3d(i,j,k,24)
          alocont_nn3d_tmp(i,j,k,25) = alocont_nn3d_tmp(i,j,k,25) + wall*alocont_o_nn3d(i,j,k,25)
          alocont_nn3d_tmp(i,j,k,26) = alocont_nn3d_tmp(i,j,k,26) + wall*alocont_o_nn3d(i,j,k,26)
          alocont_nn3d_tmp(i,j,k,27) = alocont_nn3d_tmp(i,j,k,27) + wall*alocont_o_nn3d(i,j,k,27)
        enddo
      enddo
    enddo
    !
  end subroutine formallc_cont3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine formallc_cont3d_lin(oindx,nueindx)
    !
    !-----------------------------------------------------------------------
    !------long characteristics for continuum radiative transfer in 3d------
    !-----------calculating intensties for given mu,phi specified-----------
    !---------------------------by input oindx------------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, int3d, opac3d, scont3d, imaskb3d, &
         alocont_o_nn3d, alocont_nn3d_tmp, mint3d_tmp, normalization3d_tmp, &
         fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, bnue3d, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp
    use mod_angles, only: n_x, n_y, n_z, weight_omega, q_alo
    use mod_frequencies, only: nodes_nue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: oindx, nueindx
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: alpha, beta, gamma
    integer(i4b) :: startx, starty, startz, endx, endy, endz, &
    startxb, startyb, endxb, endyb
    integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
    integer(i4b) :: aloindx_im1, aloindx_ip1, aloindx_jm1, aloindx_jp1, aloindx_km1, aloindx_kp1
    real(dp) :: nn_x, nn_y, nn_z, xnue, wall
    real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, dels_r
    real(dp) :: dels_xyd, dels_xzd, dels_yzd, dels_d
    integer :: q14, q15, q17, q18, q23, q24, q26, q27
    !
    !for debugging
    real(dp) :: int_u2, scont_u2, opac_u2, int_u3, scont_u3, opac_u3, scont_d2, scont_d3, opac_d2, opac_d3, interpol2d_9p_quad, interpol2d_9p_bez
    !
    ! ... local arrays
    !
    ! ... local functions
    real(dp) :: dist_ellipsoid, calc_icore_gdark
    !
    ! ... local characters
    !character(len=50) :: enter
    !
    ! ... local logicals
    !
    !frequency
    xnue=nodes_nue(nueindx)
    !
    !directions
    nn_x=n_x(oindx)
    nn_y=n_y(oindx)
    nn_z=n_z(oindx)
    !
    !angulare integration weight
    wall=weight_omega(oindx)
    !
    !indices for nearest neighbour alo
    q14=q_alo(oindx,14)
    q15=q_alo(oindx,15)
    q17=q_alo(oindx,17)
    q18=q_alo(oindx,18)
    q23=q_alo(oindx,23)
    q24=q_alo(oindx,24)
    q26=q_alo(oindx,26)
    q27=q_alo(oindx,27)
    !
    !-----------------------------------------------------------------------
    !
    !set directional index-parameter (phantom points are excluded from calculation)
    !
    !index parameter:
    !         if n_x >= 0                 if n_x < 0
    !                startx = 2                  startx = ndxmax-1
    !                endx = ndxmax-1             endx = 2
    !                alpha=  1                   alpha=-1
    !
    !         if n_y >= 0                 if n_y < 0
    !                starty = 2                  starty = ndymax-1
    !                endy = ndymax-1             endy = 2
    !                beta =  1                   beta =-1
    !
    !         if n_z >= 0                 if n_z < 0
    !                startz = 2                  startz = ndzmax-1
    !                endz = ndzmax-1             endz = 2
    !                gamma = 1                   gamma = -1
    !
    if(nn_x.gt.zero) then
      startx = 1
      endx = nx
      startxb = nx
      endxb = 1
      alpha=  1
      aloindx_im1 = nx-1
      aloindx_ip1 = 2
    elseif(nn_x.lt.zero) then
      startx = nx
      endx = 1
      startxb = 1
      endxb = nx
      alpha=-1
      aloindx_im1 = 2
      aloindx_ip1 = nx-1
    else
      stop 'error in formalc_cont3d_lin: n_x = 0 not allowed'
    endif
    !
    if(nn_y.gt.zero) then
      starty = 1
      endy = ny
      startyb = ny
      endyb = 1
      beta =  1
      aloindx_jm1 = ny-1
      aloindx_jp1 = 2
    elseif(nn_y.lt.zero) then
      starty = ny
      endy = 1
      startyb = 1
      endyb = ny
      beta =-1
      aloindx_jm1 = 2
      aloindx_jp1 = ny-1
    else
      stop 'error in formallc_cont3d_lin: n_y = 0 not allowed'
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
    else
      stop 'error in formallc_cont3d_lin: n_z = 0 not allowed'
    endif
    !
    !--------------------reset the intensities and alo----------------------
    !
    call set_boundary3d(xnue, nn_x, nn_y, nn_z)
    !
    !-----------------------------------------------------------------------
    !
    alocont_o_nn3d=zero
    !
    !-----------------------------------------------------------------------
    !
    do k=startz, endz, gamma
      do j=starty, endy, beta
        do i=startx, endx, alpha
          !
          select case(imaskb3d(i,j,k))
            !
            !*************left and right boundary with long characteristics*********
            !
          case(12,13)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
            !*************front and back boundary with long characteristics*********
            !
          case(14,15)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
            !********(left,front), (right,front), (left,back), (right,back)*********
            !*******************boundary edges with long characteristics************
            !
          case(18,21,24,27)
            !
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
            !***********************adjacent to left boundary***********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(4)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on left towards right
            alocont_o_nn3d(nx,j-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,j-1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,j-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,j+1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx,j+1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx,j+1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,j  ,k+1, 6) = alocont_o_nn3d(i-1,j  ,k+1, 6)
            alocont_o_nn3d(nx,j  ,k  ,15) = alocont_o_nn3d(i-1,j  ,k  ,15)
            alocont_o_nn3d(nx,j  ,k-1,24) = alocont_o_nn3d(i-1,j  ,k-1,24)
            !
            !***********************adjacent to right boundary**********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(5)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on right towards left
            alocont_o_nn3d(1,j-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,j-1,k,  16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,j-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,j+1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,j+1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,j+1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,j  ,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,j  ,k,  13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,j  ,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            !
            !***********************adjacent to front boundary**********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(6)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards front
            alocont_o_nn3d(i-1,ny,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(i-1,ny,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(i-1,ny,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(i  ,ny,k+1, 8) = alocont_o_nn3d(i  ,j-1,k+1, 8)
            alocont_o_nn3d(i  ,ny,k  ,17) = alocont_o_nn3d(i  ,j-1,k  ,17)
            alocont_o_nn3d(i  ,ny,k-1,26) = alocont_o_nn3d(i  ,j-1,k-1,26)
            alocont_o_nn3d(i+1,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(i+1,ny,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(i+1,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            !
            !***********************adjacent to back boundary***********************
            !********(standard radiative transfer + copy ALO to opposite side)******
            !
          case(7)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards front
            alocont_o_nn3d(i+1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(i+1,1,k,  10) = alocont_o_nn3d(i+1,j+1,k,  10)
            alocont_o_nn3d(i+1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(i  ,1,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(i  ,1,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(i  ,1,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(i-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(i-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(i-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            !
            !***********************adjacent to left front boundary*****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(8)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(nx,1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,2,k+1, 6) = alocont_o_nn3d(i-1,j  ,k+1, 6)
            alocont_o_nn3d(nx,2,k  ,15) = alocont_o_nn3d(i-1,j  ,k  ,15)
            alocont_o_nn3d(nx,2,k-1,24) = alocont_o_nn3d(i-1,j  ,k-1,24)
            alocont_o_nn3d(nx,j+1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx,j+1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx,j+1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(1,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(1,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(1,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(2,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(2,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(2,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(3,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(3,ny,k,  16) = alocont_o_nn3d(i+1,j-1,k,  16)
            alocont_o_nn3d(3,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx ,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx ,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx ,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            !
            !***********************adjacent to right front boundary****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(9)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)

            alocont_o_nn3d(1,2,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,2,k  ,13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,2,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            alocont_o_nn3d(1,j+1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,j+1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,j+1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,1,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,ny,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx,ny ,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(nx,ny ,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(nx,ny ,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx-2,ny ,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx-2,ny ,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx-2,ny ,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx-1,ny ,k+1, 8) = alocont_o_nn3d(i  ,j-1,k+1, 8)
            alocont_o_nn3d(nx-1,ny ,k  ,17) = alocont_o_nn3d(i  ,j-1,k  ,17)
            alocont_o_nn3d(nx-1,ny ,k-1,26) = alocont_o_nn3d(i  ,j-1,k-1,26)

            !
            !***********************adjacent to left back boundary******************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(10)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1  ,1  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(1  ,1  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(1  ,1  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(2  ,1  ,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(2  ,1  ,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(2  ,1  ,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(nx  ,1  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx  ,1  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx  ,1  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,j-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,j-1,k  ,18) = alocont_o_nn3d(i-1,j-1,k  ,18)
            alocont_o_nn3d(nx,j-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx  ,ny  ,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx  ,ny  ,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx  ,ny  ,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx ,ny-1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx ,ny-1,k,  15) = alocont_o_nn3d(i-1,j,k,  15)
            alocont_o_nn3d(nx ,ny-1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(i+1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(i+1,1,k,  10) = alocont_o_nn3d(i+1,j+1,k,  10)
            alocont_o_nn3d(i+1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            !
            !***********************adjacent to right back boundary*****************
            !*******(standard radiative transfer + copy ALO to opposite sides*******
            !
          case(11)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(nx-1,1,k+1, 2) = alocont_o_nn3d(i  ,j+1,k+1, 2)
            alocont_o_nn3d(nx-1,1,k  ,11) = alocont_o_nn3d(i  ,j+1,k  ,11)
            alocont_o_nn3d(nx-1,1,k-1,20) = alocont_o_nn3d(i  ,j+1,k-1,20)
            alocont_o_nn3d(nx,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(nx,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(nx,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,j-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,j-1,k,  16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,j-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,j  ,k+1, 4) = alocont_o_nn3d(i+1,j  ,k+1, 4)
            alocont_o_nn3d(1,j  ,k  ,13) = alocont_o_nn3d(i+1,j  ,k  ,13)
            alocont_o_nn3d(1,j  ,k-1,22) = alocont_o_nn3d(i+1,j  ,k-1,22)
            alocont_o_nn3d(1,ny,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,ny,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,ny,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(i-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(i-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(i-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            !
            !**********************left boundary with long characteristics**********
            !
          case(16)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            alocont_o_nn3d(2,ny,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(2,ny,k,  16) = alocont_o_nn3d(i+1,j-1,k,  16)
            alocont_o_nn3d(2,ny,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(nx,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(nx,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(nx,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(1,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(1,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(1,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            !
            !*********************front boundary with long characteristics**********
            !
          case(17)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            alocont_o_nn3d(nx,1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,1,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,2,k+1,3) = alocont_o_nn3d(i-1,j+1,k+1,3)
            alocont_o_nn3d(nx,2,k,12) = alocont_o_nn3d(i-1,j+1,k,12)
            alocont_o_nn3d(nx,2,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,ny,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,ny,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,ny,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            !
            !*********************right boundary with long characteristics**********
            !
          case(19)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            alocont_o_nn3d(nx-1,ny,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx-1,ny,k,  18) = alocont_o_nn3d(i-1,j-1,k,  18)
            alocont_o_nn3d(nx-1,ny,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            alocont_o_nn3d(nx,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(nx,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(nx,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            alocont_o_nn3d(1,ny,k+1, 8) = alocont_o_nn3d(i,j-1,k+1, 8)
            alocont_o_nn3d(1,ny,k,  17) = alocont_o_nn3d(i,j-1,k,  17)
            alocont_o_nn3d(1,ny,k-1,26) = alocont_o_nn3d(i,j-1,k-1,26)
            !
            !*********************front boundary with long characteristics**********
            !
          case(20)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,1,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,1,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            alocont_o_nn3d(1,2,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(1,2,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(1,2,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(1,ny,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,ny,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,ny,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            !
            !**********************left boundary with long characteristics**********
            !
          case(22)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(1,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(1,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            alocont_o_nn3d(2,1,k+1, 1) = alocont_o_nn3d(i+1,j+1,k+1, 1)
            alocont_o_nn3d(2,1,k  ,10) = alocont_o_nn3d(i+1,j+1,k  ,10)
            alocont_o_nn3d(2,1,k-1,19) = alocont_o_nn3d(i+1,j+1,k-1,19)
            alocont_o_nn3d(nx,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(nx,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(nx,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            !
            !*********************back boundary with long characteristics***********
            !
          case(23)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            alocont_o_nn3d(nx,1,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,1,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,1,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,ny,k+1, 6) = alocont_o_nn3d(i-1,j,k+1, 6)
            alocont_o_nn3d(nx,ny,k,15) = alocont_o_nn3d(i-1,j,k,15)
            alocont_o_nn3d(nx,ny,k-1,24) = alocont_o_nn3d(i-1,j,k-1,24)
            alocont_o_nn3d(nx,ny-1,k+1, 9) = alocont_o_nn3d(i-1,j-1,k+1, 9)
            alocont_o_nn3d(nx,ny-1,k,  18) = alocont_o_nn3d(i-1,j-1,k,  18)
            alocont_o_nn3d(nx,ny-1,k-1,27) = alocont_o_nn3d(i-1,j-1,k-1,27)
            !
            !**********************right boundary with long characteristics*********
            !
          case(25)
            !appropriate neighbours
            iim1=aloindx_im1
            jjm1=j-beta
            kkm1=k-gamma
            iip1=aloindx_ip1
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(1,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(1,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            alocont_o_nn3d(nx-1,1,k+1, 3) = alocont_o_nn3d(i-1,j+1,k+1, 3)
            alocont_o_nn3d(nx-1,1,k  ,12) = alocont_o_nn3d(i-1,j+1,k  ,12)
            alocont_o_nn3d(nx-1,1,k-1,21) = alocont_o_nn3d(i-1,j+1,k-1,21)
            alocont_o_nn3d(nx,1,k+1, 2) = alocont_o_nn3d(i,j+1,k+1, 2)
            alocont_o_nn3d(nx,1,k  ,11) = alocont_o_nn3d(i,j+1,k  ,11)
            alocont_o_nn3d(nx,1,k-1,20) = alocont_o_nn3d(i,j+1,k-1,20)
            !
            !**********************back boundary with long characteristics**********
            !
          case(26)
            !
            !appropriate neighbours
            iim1=i-alpha
            jjm1=aloindx_jm1
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=aloindx_jp1
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !copy alo on back side towards opposite sides
            alocont_o_nn3d(1,1,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,1,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,1,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            alocont_o_nn3d(1,ny-1,k+1, 7) = alocont_o_nn3d(i+1,j-1,k+1, 7)
            alocont_o_nn3d(1,ny-1,k  ,16) = alocont_o_nn3d(i+1,j-1,k  ,16)
            alocont_o_nn3d(1,ny-1,k-1,25) = alocont_o_nn3d(i+1,j-1,k-1,25)
            alocont_o_nn3d(1,ny,k+1, 4) = alocont_o_nn3d(i+1,j,k+1, 4)
            alocont_o_nn3d(1,ny,k  ,13) = alocont_o_nn3d(i+1,j,k  ,13)
            alocont_o_nn3d(1,ny,k-1,22) = alocont_o_nn3d(i+1,j,k-1,22)
            !
            !************************standard radiative transfer********************
            !
          case(3)
            iim1=i-alpha
            jjm1=j-beta
            kkm1=k-gamma
            iip1=i+alpha
            jjp1=j+beta
            kkp1=k+gamma
            call flc_cont3d_lin(oindx, nueindx, iim1, i, iip1, jjm1, j, jjp1, kkm1, k, kkp1, &
            startxb, startyb, endxb, endyb, &
            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
            q14, q15, q17, q18, q23, q24, q26, q27)
            !
          end select
          !
          !
        enddo
      enddo
    enddo
    !
    !
    !---------------------perform angular integration-----------------------
    !NOTE: due to periodic boundary conditions, and copying of ALO coefficients
    !      angular integration cannot (easily) be performed within loop from above
    !
    do k=1, nz
      do j=1, ny
        do i=1, nx
          mint3d_tmp(i,j,k) = mint3d_tmp(i,j,k) + int3d(i,j,k)*wall
          fcontx3d_tmp(i,j,k) = fcontx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*wall
          fconty3d_tmp(i,j,k) = fconty3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*wall
          fcontz3d_tmp(i,j,k) = fcontz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*wall
          kcontxx3d_tmp(i,j,k) = kcontxx3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_x*wall
          kcontxy3d_tmp(i,j,k) = kcontxy3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_y*wall
          kcontxz3d_tmp(i,j,k) = kcontxz3d_tmp(i,j,k) + int3d(i,j,k)*nn_x*nn_z*wall
          kcontyy3d_tmp(i,j,k) = kcontyy3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_y*wall
          kcontyz3d_tmp(i,j,k) = kcontyz3d_tmp(i,j,k) + int3d(i,j,k)*nn_y*nn_z*wall
          kcontzz3d_tmp(i,j,k) = kcontzz3d_tmp(i,j,k) + int3d(i,j,k)*nn_z*nn_z*wall
          normalization3d_tmp(i,j,k) = normalization3d_tmp(i,j,k) + wall

          alocont_nn3d_tmp(i,j,k,1) = alocont_nn3d_tmp(i,j,k,1) + wall*alocont_o_nn3d(i,j,k,1)
          alocont_nn3d_tmp(i,j,k,2) = alocont_nn3d_tmp(i,j,k,2) + wall*alocont_o_nn3d(i,j,k,2)
          alocont_nn3d_tmp(i,j,k,3) = alocont_nn3d_tmp(i,j,k,3) + wall*alocont_o_nn3d(i,j,k,3)
          alocont_nn3d_tmp(i,j,k,4) = alocont_nn3d_tmp(i,j,k,4) + wall*alocont_o_nn3d(i,j,k,4)
          alocont_nn3d_tmp(i,j,k,5) = alocont_nn3d_tmp(i,j,k,5) + wall*alocont_o_nn3d(i,j,k,5)
          alocont_nn3d_tmp(i,j,k,6) = alocont_nn3d_tmp(i,j,k,6) + wall*alocont_o_nn3d(i,j,k,6)
          alocont_nn3d_tmp(i,j,k,7) = alocont_nn3d_tmp(i,j,k,7) + wall*alocont_o_nn3d(i,j,k,7)
          alocont_nn3d_tmp(i,j,k,8) = alocont_nn3d_tmp(i,j,k,8) + wall*alocont_o_nn3d(i,j,k,8)
          alocont_nn3d_tmp(i,j,k,9) = alocont_nn3d_tmp(i,j,k,9) + wall*alocont_o_nn3d(i,j,k,9)
          alocont_nn3d_tmp(i,j,k,10) = alocont_nn3d_tmp(i,j,k,10) + wall*alocont_o_nn3d(i,j,k,10)
          alocont_nn3d_tmp(i,j,k,11) = alocont_nn3d_tmp(i,j,k,11) + wall*alocont_o_nn3d(i,j,k,11)
          alocont_nn3d_tmp(i,j,k,12) = alocont_nn3d_tmp(i,j,k,12) + wall*alocont_o_nn3d(i,j,k,12)
          alocont_nn3d_tmp(i,j,k,13) = alocont_nn3d_tmp(i,j,k,13) + wall*alocont_o_nn3d(i,j,k,13)
          alocont_nn3d_tmp(i,j,k,14) = alocont_nn3d_tmp(i,j,k,14) + wall*alocont_o_nn3d(i,j,k,14)
          alocont_nn3d_tmp(i,j,k,15) = alocont_nn3d_tmp(i,j,k,15) + wall*alocont_o_nn3d(i,j,k,15)
          alocont_nn3d_tmp(i,j,k,16) = alocont_nn3d_tmp(i,j,k,16) + wall*alocont_o_nn3d(i,j,k,16)
          alocont_nn3d_tmp(i,j,k,17) = alocont_nn3d_tmp(i,j,k,17) + wall*alocont_o_nn3d(i,j,k,17)
          alocont_nn3d_tmp(i,j,k,18) = alocont_nn3d_tmp(i,j,k,18) + wall*alocont_o_nn3d(i,j,k,18)
          alocont_nn3d_tmp(i,j,k,19) = alocont_nn3d_tmp(i,j,k,19) + wall*alocont_o_nn3d(i,j,k,19)
          alocont_nn3d_tmp(i,j,k,20) = alocont_nn3d_tmp(i,j,k,20) + wall*alocont_o_nn3d(i,j,k,20)
          alocont_nn3d_tmp(i,j,k,21) = alocont_nn3d_tmp(i,j,k,21) + wall*alocont_o_nn3d(i,j,k,21)
          alocont_nn3d_tmp(i,j,k,22) = alocont_nn3d_tmp(i,j,k,22) + wall*alocont_o_nn3d(i,j,k,22)
          alocont_nn3d_tmp(i,j,k,23) = alocont_nn3d_tmp(i,j,k,23) + wall*alocont_o_nn3d(i,j,k,23)
          alocont_nn3d_tmp(i,j,k,24) = alocont_nn3d_tmp(i,j,k,24) + wall*alocont_o_nn3d(i,j,k,24)
          alocont_nn3d_tmp(i,j,k,25) = alocont_nn3d_tmp(i,j,k,25) + wall*alocont_o_nn3d(i,j,k,25)
          alocont_nn3d_tmp(i,j,k,26) = alocont_nn3d_tmp(i,j,k,26) + wall*alocont_o_nn3d(i,j,k,26)
          alocont_nn3d_tmp(i,j,k,27) = alocont_nn3d_tmp(i,j,k,27) + wall*alocont_o_nn3d(i,j,k,27)
        enddo
      enddo
    enddo
    !
  end subroutine formallc_cont3d_lin
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine coeff2d_contu_lin(opac_im1jm1, opac_ijm1, opac_im1j, opac_ij, &
                               scont_im1jm1, scont_ijm1, scont_im1j, scont_ij, &
                               int_im1jm1, int_ijm1, int_im1j, int_ij, &
                               x_im1, x_i, y_jm1, y_j, x_p, y_p, &
                               a_scont, b_scont, c_scont, d_scont, &
                               a_inten, b_inten, c_inten, d_inten, &
                               opac_p, scont_p, int_p)
    !
    !         interpolates opacity, continuum source function and intensity
    !               values given on a 2d grid onto point x_p, y_p
    !
    !on input (f_* stands for opac_*, scont_* and int_*, respectivly):
    !
    ! y_j      f_im1j--------------f_ij
    !  |          |                  |
    !  |          |                  |
    !  |          |         x        |
    !  |          |     (x_p,y_p)    |
    !  |          |                  |
    !y_jm1    f_im1jm1------------f_ijm1
    !  |
    !  --------x_im2-----------x_im1---------------x_i
    !
    !        x_p, y_p: coordinates of point onto which shall be interpolated
    !
    !on output:
    !   1. interpolation coefficients for source function and intensity
    !      (required for ALO calculations):
    !         a_scont, b_scont, c_scont, d_scont
    !         a_inten, b_inten, c_inten, d_inten
    !
    !      such that:
    !         f_p = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
    !
    !   2. interpolated values at point p: opac_p, scont_p, int_p
    !
    !
    ! ... arguments
    real(dp), intent(in) :: opac_im1jm1, opac_ijm1, opac_im1j,   opac_ij, &
    scont_im1jm1, scont_ijm1, scont_im1j, scont_ij, &
    int_im1jm1, int_ijm1, int_im1j, int_ij, &
    x_im1, x_i, y_jm1, y_j, x_p, y_p
    real(dp), intent(out) :: a_scont, b_scont, c_scont, d_scont, &
    a_inten, b_inten, c_inten, d_inten, &
    opac_p, int_p, scont_p
    !
    ! ... local scalars
    real(dp) :: dxi, tx, dyj, ty, rdxdy
    !
    !
    !define deltax, deltay
    dxi = x_i-x_im1
    dyj = y_j-y_jm1
    !
    !define deltax, deltay-ratios
    tx = (x_p-x_im1)/dxi
    ty = (y_p-y_jm1)/dyj
    !
    !-------------------------bilinear interpolation------------------------
    !
    rdxdy=tx*ty
    !
    a_scont=1.d0-tx-ty+rdxdy
    b_scont=tx-rdxdy
    c_scont=ty-rdxdy
    d_scont=rdxdy
    !
    a_inten=a_scont
    b_inten=b_scont
    c_inten=c_scont
    d_inten=d_scont
    !
    opac_p = a_scont*opac_im1jm1 + b_scont*opac_ijm1 + c_scont*opac_im1j + d_scont*opac_ij
    scont_p = a_scont*scont_im1jm1 + b_scont*scont_ijm1 + c_scont*scont_im1j + d_scont*scont_ij
    int_p = a_scont*int_im1jm1 + b_scont*int_ijm1 + c_scont*int_im1j + d_scont*int_ij
    return
    !
  end subroutine coeff2d_contu_lin
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine coeff2d_contd_lin(opac_im1jm1, opac_ijm1, opac_im1j,   opac_ij, &
                               scont_im1jm1, scont_ijm1, scont_im1j,   scont_ij, &
                               x_im1, x_i, y_jm1, y_j, x_p, y_p, &
                               a_scont, b_scont, c_scont, d_scont, &
                               opac_p, scont_p)
    !
    !            interpolates opacity and continuum source function
    !               values given on a 2d grid onto point x_p, y_p
    !
    !on input (f_* stands for opac_* and scont_* respectivly):
    !
    ! y_j      f_im1j--------------f_ij
    !  |       |                  |
    !  |       |                  |
    !  |       |         x        |
    !  |       |     (x_p,y_p)    |
    !  |       |                  |
    !y_jm1  f_im1jm1------------f_ijm1
    !  |
    !  ------x_im1---------------x_i
    !
    !        x_p, y_p: coordinates of point onto which shall be interpolated
    !
    !on output:
    !   1. interpolation coefficients for source function
    !      (required for ALO calculations):
    !         a_scont, b_scont, c_scont, d_scont
    !
    !      such that:
    !         f_p = a*f_im1jm1 + b*f_ijm1 + c*f_im1j + d*f_ij
    !
    !   2. interpolated values at point p: opac_p, scont_p
    !
    !
    ! ... arguments
    real(dp), intent(in) :: opac_im1jm1, opac_ijm1, opac_im1j, opac_ij, &
    scont_im1jm1, scont_ijm1, scont_im1j, scont_ij, &
    x_im1, x_i, y_jm1, y_j, x_p, y_p
    real(dp), intent(out) :: a_scont, b_scont, c_scont, d_scont, opac_p, scont_p
    !
    ! ... local scalars
    real(dp) :: dxi, tx, dyj, ty, rdxdy
    !
    !define deltax, deltay
    dxi = x_i-x_im1
    dyj = y_j-y_jm1
    !
    !define deltax, deltay-ratios
    tx = (x_p-x_im1)/dxi
    ty = (y_p-y_jm1)/dyj
    !
    !-------------------------bilinear interpolation------------------------
    !
    rdxdy=tx*ty
    !
    a_scont=1.d0-tx-ty+rdxdy
    b_scont=tx-rdxdy
    c_scont=ty-rdxdy
    d_scont=rdxdy
    !
    opac_p = a_scont*opac_im1jm1 + b_scont*opac_ijm1 + c_scont*opac_im1j + d_scont*opac_ij
    scont_p = a_scont*scont_im1jm1 + b_scont*scont_ijm1 + c_scont*scont_im1j + d_scont*scont_ij
    !
    !
  end subroutine coeff2d_contd_lin
  !
  !***********************************************************************
  !***********************************************************************
  !
  !                 SHORT CHARACTERISTICS FOR A CELL
  !
  !***********************************************************************
  !***********************************************************************
  !
  subroutine fsc_cont3d(iim1, ii, iip1, jjm1, jj, jjp1, kkm1, kk, kkp1, &
                        nn_x, nn_y, nn_z, wall, &
                        q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
                        q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
                        q25, q26, q27)
    !
    !-----------------------------------------------------------------------
    !------short characteristics for continuum radiative transfer in 3d-----
    !------------------at any given point x_p, y_p, z_p---------------------
    !-----------calculating intensty for given mu,phi specified-------------
    !----------------------by input nn_x, nn_y, nn_z------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, opac3d, scont3d, int3d, alocont_o_nn3d, &
         fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
         mint3d_tmp, alocont_nn3d_tmp, normalization3d_tmp, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp
    use mod_angles, only: n_x, n_y, n_z
    use mod_frequencies, only: nodes_nue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: iim1, ii, iip1, jjm1, jj, jjp1, kkm1, kk, kkp1, &
    q1, q2, q3, q4, q5, q6, q7, q8, q9, &
    q10, q11, q12, q13, q14, q15, q16, q17, q18, &
    q19, q20, q21, q22, q23, q24, q25, q26, q27
    real(dp), intent(in) :: nn_x, nn_y, nn_z, wall
    !
    ! ... local scalars
    real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, &
    dels_xyd, dels_xzd, dels_yzd, dels_d
    real(dp) :: opac_p, scont_p
    real(dp) :: x_u, y_u, z_u, int_u, opac_u, scont_u
    real(dp) :: x_d, y_d, z_d, opac_d, scont_d
    real(dp) :: abs_sc, int_sc, contr_sc
    real(dp) :: alo_u, alo_p, alo_d
    real(dp) :: c14_scontu, c15_scontu, c17_scontu, c18_scontu, c23_scontu, c24_scontu, c26_scontu, &
    c14_opacu,  c15_opacu,  c17_opacu,  c18_opacu,  c23_opacu,  c24_opacu,  c26_opacu, &
    c14_intu,   c15_intu,   c17_intu,   c18_intu,   c23_intu,   c24_intu,   c26_intu, &
    c15_scontd, c17_scontd, c18_scontd, c23_scontd, c24_scontd, c26_scontd, c27_scontd, &
    c15_opacd,  c17_opacd,  c18_opacd,  c23_opacd,  c24_opacd,  c26_opacd,  c27_opacd
    !
    !calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
    dels_xyu=(z(kk)-z(kkm1))/nn_z
    dels_xzu=(y(jj)-y(jjm1))/nn_y
    dels_yzu=(x(ii)-x(iim1))/nn_x
    dels_u=min(dels_xyu,dels_xzu,dels_yzu)
    !
    !calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
    dels_xyd=(z(kkp1)-z(kk))/nn_z
    dels_xzd=(y(jjp1)-y(jj))/nn_y
    dels_yzd=(x(iip1)-x(ii))/nn_x
    dels_d=min(dels_xyd,dels_xzd,dels_yzd)
    !
    !----------------------------local point--------------------------------
    !
    scont_p=scont3d(ii,jj,kk)
    opac_p=opac3d(ii,jj,kk)
    !
    !----------------------------upwind point-------------------------------
    !
    if(dels_xyu.eq.dels_u) then
      !intersection with x-y plane on level k-gamma
      x_u = x(ii) - dels_u*nn_x
      y_u = y(jj) - dels_u*nn_y
      z_u = z(kkm1)
      !
      call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(ii,jjm1,kkm1), &
      opac3d(iim1,jj,kkm1), opac3d(ii,jj,kkm1), &
      scont3d(iim1,jjm1,kkm1), scont3d(ii,jjm1,kkm1), &
      scont3d(iim1,jj,kkm1), scont3d(ii,jj,kkm1), &
      int3d(iim1,jjm1,kkm1), int3d(ii,jjm1,kkm1), &
      int3d(iim1,jj,kkm1), int3d(ii,jj,kkm1), &
      x(iim1), x(ii), y(jjm1), y(jj), x_u, y_u, &
      c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
      c14_intu, c15_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
      !
      !set interpolation coefficients that are not used to zero
      c23_scontu = zero
      c24_scontu = zero
      c26_scontu = zero
      c23_intu = zero
      c24_intu = zero
      c26_intu = zero
      !
    elseif(dels_xzu.eq.dels_u) then
      !intersection with x-z plane on level j-beta
      x_u = x(ii) - dels_u*nn_x
      y_u = y(jjm1)
      z_u = z(kk) - dels_u*nn_z
      !
      call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(ii,jjm1,kkm1), &
      opac3d(iim1,jjm1,kk), opac3d(ii,jjm1,kk), &
      scont3d(iim1,jjm1,kkm1), scont3d(ii,jjm1,kkm1), &
      scont3d(iim1,jjm1,kk), scont3d(ii,jjm1,kk), &
      int3d(iim1,jjm1,kkm1), int3d(ii,jjm1,kkm1), &
      int3d(iim1,jjm1,kk), int3d(ii,jjm1,kk), &
      x(iim1), x(ii), z(kkm1), z(kk), x_u, z_u, &
      c14_scontu, c15_scontu, c23_scontu, c24_scontu, &
      c14_intu, c15_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
      !
      !set interpolation coefficients that are not used to zero
      c17_scontu = zero
      c18_scontu = zero
      c26_scontu = zero
      c17_intu = zero
      c18_intu = zero
      c26_intu = zero
      !
    elseif(dels_yzu.eq.dels_u) then
      !intersection with y-z plane on level i-alpha
      x_u = x(iim1)
      y_u = y(jj) - dels_u*nn_y
      z_u = z(kk) - dels_u*nn_z
      !
      call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,jj,kkm1), &
      opac3d(iim1,jjm1,kk), opac3d(iim1,jj,kk), &
      scont3d(iim1,jjm1,kkm1), scont3d(iim1,jj,kkm1), &
      scont3d(iim1,jjm1,kk), scont3d(iim1,jj,kk), &
      int3d(iim1,jjm1,kkm1), int3d(iim1,jj,kkm1), &
      int3d(iim1,jjm1,kk), int3d(iim1,jj,kk), &
      y(jjm1), y(jj), z(kkm1), z(kk), y_u, z_u, &
      c14_scontu, c17_scontu, c23_scontu, c26_scontu, &
      c14_intu, c17_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontu = zero
      c18_scontu = zero
      c24_scontu = zero
      c15_intu = zero
      c18_intu = zero
      c24_intu = zero
      !
    else
      write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
      stop 'error in fsc_cont3d: invalid dels_u'
    endif
    !
    !---------------------------downwind point------------------------------
    !
    if(dels_xyd.eq.dels_d) then
      !intersection with x-y plane on level k+gamma
      x_d = x(ii) + dels_d*nn_x
      y_d = y(jj) + dels_d*nn_y
      z_d = z(kkp1)
      !
      call coeff2d_contd_lin(opac3d(ii,jj,kkp1), opac3d(iip1,jj,kkp1), &
      opac3d(ii,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(ii,jj,kkp1), scont3d(iip1,jj,kkp1), &
      scont3d(ii,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
      x(ii), x(iip1), y(jj), y(jjp1), x_d, y_d, &
      c23_scontd, c24_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontd = zero
      c17_scontd = zero
      c18_scontd = zero
      !
    elseif(dels_xzd.eq.dels_d) then
      !intersection with x-z plane on level j+beta
      x_d = x(ii) + dels_d*nn_x
      y_d = y(jjp1)
      z_d = z(kk) + dels_d*nn_z
      !
      call coeff2d_contd_lin(opac3d(ii,jjp1,kk), opac3d(iip1,jjp1,kk), &
      opac3d(ii,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(ii,jjp1,kk), scont3d(iip1,jjp1,kk), &
      scont3d(ii,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
      x(ii), x(iip1), z(kk), z(kkp1), x_d, z_d, &
      c17_scontd, c18_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontd = zero
      c23_scontd = zero
      c24_scontd = zero
      !
    elseif(dels_yzd.eq.dels_d) then
      !intersection with y-z plane on level i+alpha
      x_d = x(iip1)
      y_d = y(jj) + dels_d*nn_y
      z_d = z(kk) + dels_d*nn_z
      !
      call coeff2d_contd_lin(opac3d(iip1,jj,kk), opac3d(iip1,jjp1,kk), &
      opac3d(iip1,jj,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(iip1,jj,kk), scont3d(iip1,jjp1,kk), &
      scont3d(iip1,jj,kkp1), scont3d(iip1,jjp1,kkp1), &
      y(jj), y(jjp1), z(kk), z(kkp1), y_d, z_d, &
      c15_scontd, c18_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c17_scontd = zero
      c23_scontd = zero
      c26_scontd = zero
      !
    else
      write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
      stop 'error in fsc_cont3d: invalid dels_d'
    endif
    !
    !--------------------------------radiative transfer---------------------
    !
    call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
    dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
    !call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
    !                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
    !alo_d=zero
    !
    int3d(ii,jj,kk) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p + alo_d*scont_d
    !
    alocont_o_nn3d(iip1,jjp1,kkp1,q1) = alo_d*c27_scontd
    alocont_o_nn3d(ii,jjp1,kkp1,q2) = (alo_d*c26_scontd + abs_sc*c26_intu*alocont_o_nn3d(ii,jjp1,kkp1,q1))
    alocont_o_nn3d(iim1,jjp1,kkp1,q3) = (abs_sc*c26_intu*alocont_o_nn3d(iim1,jjp1,kkp1,q2))
    alocont_o_nn3d(iip1,jj,kkp1,q4) = (alo_d*c24_scontd + abs_sc*c24_intu*alocont_o_nn3d(iip1,jj,kkp1,q1))
    alocont_o_nn3d(ii,jj,kkp1,q5) = (alo_d*c23_scontd + abs_sc*(c23_intu*alocont_o_nn3d(ii,jj,kkp1,q1) + &
    c24_intu*alocont_o_nn3d(ii,jj,kkp1,q2) + &
    c26_intu*alocont_o_nn3d(ii,jj,kkp1,q4)))
    alocont_o_nn3d(iim1,jj,kkp1,q6) = (abs_sc*(c23_intu*alocont_o_nn3d(iim1,jj,kkp1,q2) + &
    c24_intu*alocont_o_nn3d(iim1,jj,kkp1,q3) + &
    c26_intu*alocont_o_nn3d(iim1,jj,kkp1,q5)))
    alocont_o_nn3d(iip1,jjm1,kkp1,q7) = (abs_sc*c24_intu*alocont_o_nn3d(iip1,jjm1,kkp1,q4))
    alocont_o_nn3d(ii,jjm1,kkp1,q8) = (abs_sc*(c23_intu*alocont_o_nn3d(ii,jjm1,kkp1,q4) + &
    c24_intu*alocont_o_nn3d(ii,jjm1,kkp1,q5) + &
    c26_intu*alocont_o_nn3d(ii,jjm1,kkp1,q7)))
    alocont_o_nn3d(iim1,jjm1,kkp1,q9) = (abs_sc*(c23_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q5) + &
    c24_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q6) + &
    c26_intu*alocont_o_nn3d(iim1,jjm1,kkp1,q8)))
    alocont_o_nn3d(iip1,jjp1,kk,q10) = (alo_d*c18_scontd + abs_sc*c18_intu*alocont_o_nn3d(iip1,jjp1,kk,q1))
    alocont_o_nn3d(ii,jjp1,kk,q11) = (alo_d*c17_scontd + abs_sc*(c17_intu*alocont_o_nn3d(ii,jjp1,kk,q1) + &
    c18_intu*alocont_o_nn3d(ii,jjp1,kk,q2) + &
    c26_intu*alocont_o_nn3d(ii,jjp1,kk,q10)))
    alocont_o_nn3d(iim1,jjp1,kk,q12) = (abs_sc*(c17_intu*alocont_o_nn3d(iim1,jjp1,kk,q2) + &
    c18_intu*alocont_o_nn3d(iim1,jjp1,kk,q3) + &
    c26_intu*alocont_o_nn3d(iim1,jjp1,kk,q11)))
    alocont_o_nn3d(iip1,jj,kk,q13) = (alo_d*c15_scontd + abs_sc*(c15_intu*alocont_o_nn3d(iip1,jj,kk,q1) + &
    c18_intu*alocont_o_nn3d(iip1,jj,kk,q4) + &
    c24_intu*alocont_o_nn3d(iip1,jj,kk,q10)))
    alocont_o_nn3d(ii,jj,kk,q14) = (alo_p + abs_sc*(c14_intu*alocont_o_nn3d(ii,jj,kk,q1) + &
    c15_intu*alocont_o_nn3d(ii,jj,kk,q2) + &
    c17_intu*alocont_o_nn3d(ii,jj,kk,q4) + &
    c26_intu*alocont_o_nn3d(ii,jj,kk,q13) + &
    c18_intu*alocont_o_nn3d(ii,jj,kk,q5) + &
    c23_intu*alocont_o_nn3d(ii,jj,kk,q10) + &
    c24_intu*alocont_o_nn3d(ii,jj,kk,q11)))
    alocont_o_nn3d(iim1,jj,kk,q15) = (alo_u*c26_scontu + abs_sc*(c14_intu*alocont_o_nn3d(iim1,jj,kk,q2) + &
    c15_intu*alocont_o_nn3d(iim1,jj,kk,q3) + &
    c17_intu*alocont_o_nn3d(iim1,jj,kk,q5) + &
    c18_intu*alocont_o_nn3d(iim1,jj,kk,q6) + &
    c23_intu*alocont_o_nn3d(iim1,jj,kk,q11) + &
    c24_intu*alocont_o_nn3d(iim1,jj,kk,q12) + &
    c26_intu*alocont_o_nn3d(iim1,jj,kk,q14)))
    alocont_o_nn3d(iip1,jjm1,kk,q16) = (abs_sc*(c15_intu*alocont_o_nn3d(iip1,jjm1,kk,q4) + &
    c18_intu*alocont_o_nn3d(iip1,jjm1,kk,q7) + &
    c24_intu*alocont_o_nn3d(iip1,jjm1,kk,q13)))
    alocont_o_nn3d(ii,jjm1,kk,q17) = (alo_u*c24_scontu + abs_sc*(c14_intu*alocont_o_nn3d(ii,jjm1,kk,q4) + &
    c15_intu*alocont_o_nn3d(ii,jjm1,kk,q5) + &
    c17_intu*alocont_o_nn3d(ii,jjm1,kk,q7) + &
    c18_intu*alocont_o_nn3d(ii,jjm1,kk,q8) + &
    c23_intu*alocont_o_nn3d(ii,jjm1,kk,q13) + &
    c24_intu*alocont_o_nn3d(ii,jjm1,kk,q14) + &
    c26_intu*alocont_o_nn3d(ii,jjm1,kk,q16)))
    alocont_o_nn3d(iim1,jjm1,kk,q18) = (alo_u*c23_scontu + abs_sc*(c14_intu*alocont_o_nn3d(iim1,jjm1,kk,q5) + &
    c15_intu*alocont_o_nn3d(iim1,jjm1,kk,q6) + &
    c17_intu*alocont_o_nn3d(iim1,jjm1,kk,q8) + &
    c18_intu*alocont_o_nn3d(iim1,jjm1,kk,q9) + &
    c23_intu*alocont_o_nn3d(iim1,jjm1,kk,q14) + &
    c24_intu*alocont_o_nn3d(iim1,jjm1,kk,q15) + &
    c26_intu*alocont_o_nn3d(iim1,jjm1,kk,q17)))
    alocont_o_nn3d(iip1,jjp1,kkm1,q19) = (abs_sc*c18_intu*alocont_o_nn3d(iip1,jjp1,kkm1,q10))
    alocont_o_nn3d(ii,jjp1,kkm1,q20) = (abs_sc*(c17_intu*alocont_o_nn3d(ii,jjp1,kkm1,q10) + &
    c18_intu*alocont_o_nn3d(ii,jjp1,kkm1,q11) + &
    c26_intu*alocont_o_nn3d(ii,jjp1,kkm1,q19)))
    alocont_o_nn3d(iim1,jjp1,kkm1,q21) = (abs_sc*(c17_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q11) + &
    c18_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q12) + &
    c26_intu*alocont_o_nn3d(iim1,jjp1,kkm1,q20)))
    alocont_o_nn3d(iip1,jj,kkm1,q22) = (abs_sc*(c15_intu*alocont_o_nn3d(iip1,jj,kkm1,q10) + &
    c18_intu*alocont_o_nn3d(iip1,jj,kkm1,q13) + &
    c24_intu*alocont_o_nn3d(iip1,jj,kkm1,q19)))
    alocont_o_nn3d(ii,jj,kkm1,q23) = (alo_u*c18_scontu + abs_sc*(c14_intu*alocont_o_nn3d(ii,jj,kkm1,q10) + &
    c15_intu*alocont_o_nn3d(ii,jj,kkm1,q11) + &
    c17_intu*alocont_o_nn3d(ii,jj,kkm1,q13) + &
    c18_intu*alocont_o_nn3d(ii,jj,kkm1,q14) + &
    c23_intu*alocont_o_nn3d(ii,jj,kkm1,q19) + &
    c24_intu*alocont_o_nn3d(ii,jj,kkm1,q20) + &
    c26_intu*alocont_o_nn3d(ii,jj,kkm1,q22)))
    alocont_o_nn3d(iim1,jj,kkm1,q24) = (alo_u*c17_scontu + abs_sc*(c14_intu*alocont_o_nn3d(iim1,jj,kkm1,q11) + &
    c15_intu*alocont_o_nn3d(iim1,jj,kkm1,q12) + &
    c17_intu*alocont_o_nn3d(iim1,jj,kkm1,q14) + &
    c18_intu*alocont_o_nn3d(iim1,jj,kkm1,q15) + &
    c23_intu*alocont_o_nn3d(iim1,jj,kkm1,q20) + &
    c24_intu*alocont_o_nn3d(iim1,jj,kkm1,q21) + &
    c26_intu*alocont_o_nn3d(iim1,jj,kkm1,q23)))
    alocont_o_nn3d(iip1,jjm1,kkm1,q25) = (abs_sc*(c15_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q13) + &
    c18_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q16) + &
    c24_intu*alocont_o_nn3d(iip1,jjm1,kkm1,q22)))
    alocont_o_nn3d(ii,jjm1,kkm1,q26) = (alo_u*c15_scontu + abs_sc*(c14_intu*alocont_o_nn3d(ii,jjm1,kkm1,q13) + &
    c15_intu*alocont_o_nn3d(ii,jjm1,kkm1,q14) + &
    c17_intu*alocont_o_nn3d(ii,jjm1,kkm1,q16) + &
    c18_intu*alocont_o_nn3d(ii,jjm1,kkm1,q17) + &
    c23_intu*alocont_o_nn3d(ii,jjm1,kkm1,q22) + &
    c24_intu*alocont_o_nn3d(ii,jjm1,kkm1,q23) + &
    c26_intu*alocont_o_nn3d(ii,jjm1,kkm1,q25)))
    alocont_o_nn3d(iim1,jjm1,kkm1,q27) = (alo_u*c14_scontu + abs_sc*(c14_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q14) + &
    c15_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q15) + &
    c17_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q17) + &
    c18_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q18) + &
    c23_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q23) + &
    c24_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q24) + &
    c26_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q26)))
    !
    !if(ii.eq.10.and.jj.eq.10.and.kk.eq.10) then
    !   write(*,*) nn_x, nn_y, nn_z, alocont_o_nn3d(i,j,k,14), int3d(ii,jj,kk)
    !endif
    !
    !mint3d_tmp(ii,jj,kk) = mint3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*wall
    !fcontx3d_tmp(ii,jj,kk) = fcontx3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_x*wall
    !fconty3d_tmp(ii,jj,kk) = fconty3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_y*wall
    !fcontz3d_tmp(ii,jj,kk) = fcontz3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_z*wall
    !kcontxx3d_tmp(ii,jj,kk) = kcontxx3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_x*nn_x*wall
    !kcontxy3d_tmp(ii,jj,kk) = kcontxy3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_x*nn_y*wall
    !kcontxz3d_tmp(ii,jj,kk) = kcontxz3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_x*nn_z*wall
    !kcontyy3d_tmp(ii,jj,kk) = kcontyy3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_y*nn_y*wall
    !kcontyz3d_tmp(ii,jj,kk) = kcontyz3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_y*nn_z*wall
    !kcontzz3d_tmp(ii,jj,kk) = kcontzz3d_tmp(ii,jj,kk) + int3d(ii,jj,kk)*nn_z*nn_z*wall
    !normalization3d_tmp(ii,jj,kk) = normalization3d_tmp(ii,jj,kk) + wall
    !
    !alocont_nn3d_tmp(iip1,jjp1,kkp1,q1) = alocont_nn3d_tmp(iip1,jjp1,kkp1,q1) + wall*alocont_o_nn3d(iip1,jjp1,kkp1,q1)
    !alocont_nn3d_tmp(ii,jjp1,kkp1,q2) = alocont_nn3d_tmp(ii,jjp1,kkp1,q2) + wall*alocont_o_nn3d(ii,jjp1,kkp1,q2)
    !alocont_nn3d_tmp(iim1,jjp1,kkp1,q3) = alocont_nn3d_tmp(iim1,jjp1,kkp1,q3) + wall*alocont_o_nn3d(iim1,jjp1,kkp1,q3)
    !alocont_nn3d_tmp(iip1,jj,kkp1,q4) = alocont_nn3d_tmp(iip1,jj,kkp1,q4) + wall*alocont_o_nn3d(iip1,jj,kkp1,q4)
    !alocont_nn3d_tmp(ii,jj,kkp1,q5) = alocont_nn3d_tmp(ii,jj,kkp1,q5) + wall*alocont_o_nn3d(ii,jj,kkp1,q5)
    !alocont_nn3d_tmp(iim1,jj,kkp1,q6) = alocont_nn3d_tmp(iim1,jj,kkp1,q6) + wall*alocont_o_nn3d(iim1,jj,kkp1,q6)
    !alocont_nn3d_tmp(iip1,jjm1,kkp1,q7) = alocont_nn3d_tmp(iip1,jjm1,kkp1,q7) + wall*alocont_o_nn3d(iip1,jjm1,kkp1,q7)
    !alocont_nn3d_tmp(ii,jjm1,kkp1,q8) = alocont_nn3d_tmp(ii,jjm1,kkp1,q8) + wall*alocont_o_nn3d(ii,jjm1,kkp1,q8)
    !alocont_nn3d_tmp(iim1,jjm1,kkp1,q9) = alocont_nn3d_tmp(iim1,jjm1,kkp1,q9) + wall*alocont_o_nn3d(iim1,jjm1,kkp1,q9)
    !alocont_nn3d_tmp(iip1,jjp1,kk,q10) = alocont_nn3d_tmp(iip1,jjp1,kk,q10) + wall*alocont_o_nn3d(iip1,jjp1,kk,q10)
    !alocont_nn3d_tmp(ii,jjp1,kk,q11) = alocont_nn3d_tmp(ii,jjp1,kk,q11) + wall*alocont_o_nn3d(ii,jjp1,kk,q11)
    !alocont_nn3d_tmp(iim1,jjp1,kk,q12) = alocont_nn3d_tmp(iim1,jjp1,kk,q12) + wall*alocont_o_nn3d(iim1,jjp1,kk,q12)
    !alocont_nn3d_tmp(iip1,jj,kk,q13) = alocont_nn3d_tmp(iip1,jj,kk,q13) + wall*alocont_o_nn3d(iip1,jj,kk,q13)
    !alocont_nn3d_tmp(ii,jj,kk,q14) = alocont_nn3d_tmp(ii,jj,kk,q14) + wall*alocont_o_nn3d(ii,jj,kk,q14)
    !alocont_nn3d_tmp(iim1,jj,kk,q15) = alocont_nn3d_tmp(iim1,jj,kk,q15) + wall*alocont_o_nn3d(iim1,jj,kk,q15)
    !alocont_nn3d_tmp(iip1,jjm1,kk,q16) = alocont_nn3d_tmp(iip1,jjm1,kk,q16) + wall*alocont_o_nn3d(iip1,jjm1,kk,q16)
    !alocont_nn3d_tmp(ii,jjm1,kk,q17) = alocont_nn3d_tmp(ii,jjm1,kk,q17) + wall*alocont_o_nn3d(ii,jjm1,kk,q17)
    !alocont_nn3d_tmp(iim1,jjm1,kk,q18) = alocont_nn3d_tmp(iim1,jjm1,kk,q18) + wall*alocont_o_nn3d(iim1,jjm1,kk,q18)
    !alocont_nn3d_tmp(iip1,jjp1,kkm1,q19) = alocont_nn3d_tmp(iip1,jjp1,kkm1,q19) + wall*alocont_o_nn3d(iip1,jjp1,kkm1,q19)
    !alocont_nn3d_tmp(ii,jjp1,kkm1,q20) = alocont_nn3d_tmp(ii,jjp1,kkm1,q20) + wall*alocont_o_nn3d(ii,jjp1,kkm1,q20)
    !alocont_nn3d_tmp(iim1,jjp1,kkm1,q21) = alocont_nn3d_tmp(iim1,jjp1,kkm1,q21) + wall*alocont_o_nn3d(iim1,jjp1,kkm1,q21)
    !alocont_nn3d_tmp(iip1,jj,kkm1,q22) = alocont_nn3d_tmp(iip1,jj,kkm1,q22) + wall*alocont_o_nn3d(iip1,jj,kkm1,q22)
    !alocont_nn3d_tmp(ii,jj,kkm1,q23) = alocont_nn3d_tmp(ii,jj,kkm1,q23) + wall*alocont_o_nn3d(ii,jj,kkm1,q23)
    !alocont_nn3d_tmp(iim1,jj,kkm1,q24) = alocont_nn3d_tmp(iim1,jj,kkm1,q24) + wall*alocont_o_nn3d(iim1,jj,kkm1,q24)
    !alocont_nn3d_tmp(iip1,jjm1,kkm1,q25) = alocont_nn3d_tmp(iip1,jjm1,kkm1,q25) + wall*alocont_o_nn3d(iip1,jjm1,kkm1,q25)
    !alocont_nn3d_tmp(ii,jjm1,kkm1,q26) = alocont_nn3d_tmp(ii,jjm1,kkm1,q26) + wall*alocont_o_nn3d(ii,jjm1,kkm1,q26)
    !alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) = alocont_nn3d_tmp(iim1,jjm1,kkm1,q27) + wall*alocont_o_nn3d(iim1,jjm1,kkm1,q27)
    !
    !
    !
    !
  end subroutine fsc_cont3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine fsc_cont3d_lin(iim1, ii, iip1, jjm1, jj, jjp1, kkm1, kk, kkp1, &
                            nn_x, nn_y, nn_z, wall, &
                            q14, q15, q17, q18, q23, q24, q26, q27)
    !
    !-----------------------------------------------------------------------
    !------short characteristics for continuum radiative transfer in 3d-----
    !------------------at any given point x_p, y_p, z_p---------------------
    !-----------calculating intensty for given mu,phi specified-------------
    !----------------------by input nn_x, nn_y, nn_z------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, opac3d, scont3d, int3d, alocont_o_nn3d, &
         fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
         mint3d_tmp, alocont_nn3d_tmp, normalization3d_tmp, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp
    use mod_angles, only: n_x, n_y, n_z
    use mod_frequencies, only: nodes_nue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: iim1, ii, iip1, jjm1, jj, jjp1, kkm1, kk, kkp1, &
    q14, q15, q17, q18, q23, q24, q26, q27
    real(dp), intent(in) :: nn_x, nn_y, nn_z, wall
    !
    ! ... local scalars
    real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_u, &
    dels_xyd, dels_xzd, dels_yzd, dels_d
    real(dp) :: opac_p, scont_p
    real(dp) :: x_u, y_u, z_u, int_u, opac_u, scont_u
    real(dp) :: x_d, y_d, z_d, opac_d, scont_d
    real(dp) :: abs_sc, int_sc, contr_sc
    real(dp) :: alo_u, alo_p
    real(dp) :: c14_scontu, c15_scontu, c17_scontu, c18_scontu, c23_scontu, c24_scontu, c26_scontu, &
    c14_opacu,  c15_opacu,  c17_opacu,  c18_opacu,  c23_opacu,  c24_opacu,  c26_opacu, &
    c14_intu,   c15_intu,   c17_intu,   c18_intu,   c23_intu,   c24_intu,   c26_intu, &
    c15_scontd, c17_scontd, c18_scontd, c23_scontd, c24_scontd, c26_scontd, c27_scontd, &
    c15_opacd,  c17_opacd,  c18_opacd,  c23_opacd,  c24_opacd,  c26_opacd,  c27_opacd
    !
    !calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
    dels_xyu=(z(kk)-z(kkm1))/nn_z
    dels_xzu=(y(jj)-y(jjm1))/nn_y
    dels_yzu=(x(ii)-x(iim1))/nn_x
    dels_u=min(dels_xyu,dels_xzu,dels_yzu)
    !
    !calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
    dels_xyd=(z(kkp1)-z(kk))/nn_z
    dels_xzd=(y(jjp1)-y(jj))/nn_y
    dels_yzd=(x(iip1)-x(ii))/nn_x
    dels_d=min(dels_xyd,dels_xzd,dels_yzd)
    !
    !----------------------------local point--------------------------------
    !
    scont_p=scont3d(ii,jj,kk)
    opac_p=opac3d(ii,jj,kk)
    !
    !----------------------------upwind point-------------------------------
    !
    if(dels_xyu.eq.dels_u) then
      !intersection with x-y plane on level k-gamma
      x_u = x(ii) - dels_u*nn_x
      y_u = y(jj) - dels_u*nn_y
      z_u = z(kkm1)
      !
      call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(ii,jjm1,kkm1), &
      opac3d(iim1,jj,kkm1), opac3d(ii,jj,kkm1), &
      scont3d(iim1,jjm1,kkm1), scont3d(ii,jjm1,kkm1), &
      scont3d(iim1,jj,kkm1), scont3d(ii,jj,kkm1), &
      int3d(iim1,jjm1,kkm1), int3d(ii,jjm1,kkm1), &
      int3d(iim1,jj,kkm1), int3d(ii,jj,kkm1), &
      x(iim1), x(ii), y(jjm1), y(jj), x_u, y_u, &
      c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
      c14_intu, c15_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
      !
      !set interpolation coefficients that are not used to zero
      c23_scontu = zero
      c24_scontu = zero
      c26_scontu = zero
      c23_intu = zero
      c24_intu = zero
      c26_intu = zero
      !
    elseif(dels_xzu.eq.dels_u) then
      !intersection with x-z plane on level j-beta
      x_u = x(ii) - dels_u*nn_x
      y_u = y(jjm1)
      z_u = z(kk) - dels_u*nn_z
      !
      call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(ii,jjm1,kkm1), &
      opac3d(iim1,jjm1,kk), opac3d(ii,jjm1,kk), &
      scont3d(iim1,jjm1,kkm1), scont3d(ii,jjm1,kkm1), &
      scont3d(iim1,jjm1,kk), scont3d(ii,jjm1,kk), &
      int3d(iim1,jjm1,kkm1), int3d(ii,jjm1,kkm1), &
      int3d(iim1,jjm1,kk), int3d(ii,jjm1,kk), &
      x(iim1), x(ii), z(kkm1), z(kk), x_u, z_u, &
      c14_scontu, c15_scontu, c23_scontu, c24_scontu, &
      c14_intu, c15_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
      !
      !set interpolation coefficients that are not used to zero
      c17_scontu = zero
      c18_scontu = zero
      c26_scontu = zero
      c17_intu = zero
      c18_intu = zero
      c26_intu = zero
      !
    elseif(dels_yzu.eq.dels_u) then
      !intersection with y-z plane on level i-alpha
      x_u = x(iim1)
      y_u = y(jj) - dels_u*nn_y
      z_u = z(kk) - dels_u*nn_z
      !
      call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,jj,kkm1), &
      opac3d(iim1,jjm1,kk), opac3d(iim1,jj,kk), &
      scont3d(iim1,jjm1,kkm1), scont3d(iim1,jj,kkm1), &
      scont3d(iim1,jjm1,kk), scont3d(iim1,jj,kk), &
      int3d(iim1,jjm1,kkm1), int3d(iim1,jj,kkm1), &
      int3d(iim1,jjm1,kk), int3d(iim1,jj,kk), &
      y(jjm1), y(jj), z(kkm1), z(kk), y_u, z_u, &
      c14_scontu, c17_scontu, c23_scontu, c26_scontu, &
      c14_intu, c17_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontu = zero
      c18_scontu = zero
      c24_scontu = zero
      c15_intu = zero
      c18_intu = zero
      c24_intu = zero
      !
    else
      write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
      stop 'error in fsc_cont3d: invalid dels_u'
    endif
    !
    !---------------------------downwind point------------------------------
    !
    if(dels_xyd.eq.dels_d) then
      !intersection with x-y plane on level k+gamma
      x_d = x(ii) + dels_d*nn_x
      y_d = y(jj) + dels_d*nn_y
      z_d = z(kkp1)
      !
      call coeff2d_contd_lin(opac3d(ii,jj,kkp1), opac3d(iip1,jj,kkp1), &
      opac3d(ii,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(ii,jj,kkp1), scont3d(iip1,jj,kkp1), &
      scont3d(ii,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
      x(ii), x(iip1), y(jj), y(jjp1), x_d, y_d, &
      c23_scontd, c24_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontd = zero
      c17_scontd = zero
      c18_scontd = zero
      !
    elseif(dels_xzd.eq.dels_d) then
      !intersection with x-z plane on level j+beta
      x_d = x(ii) + dels_d*nn_x
      y_d = y(jjp1)
      z_d = z(kk) + dels_d*nn_z
      !
      call coeff2d_contd_lin(opac3d(ii,jjp1,kk), opac3d(iip1,jjp1,kk), &
      opac3d(ii,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(ii,jjp1,kk), scont3d(iip1,jjp1,kk), &
      scont3d(ii,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
      x(ii), x(iip1), z(kk), z(kkp1), x_d, z_d, &
      c17_scontd, c18_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontd = zero
      c23_scontd = zero
      c24_scontd = zero
      !
    elseif(dels_yzd.eq.dels_d) then
      !intersection with y-z plane on level i+alpha
      x_d = x(iip1)
      y_d = y(jj) + dels_d*nn_y
      z_d = z(kk) + dels_d*nn_z
      !
      call coeff2d_contd_lin(opac3d(iip1,jj,kk), opac3d(iip1,jjp1,kk), &
      opac3d(iip1,jj,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(iip1,jj,kk), scont3d(iip1,jjp1,kk), &
      scont3d(iip1,jj,kkp1), scont3d(iip1,jjp1,kkp1), &
      y(jj), y(jjp1), z(kk), z(kkp1), y_d, z_d, &
      c15_scontd, c18_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c17_scontd = zero
      c23_scontd = zero
      c26_scontd = zero
      !
    else
      write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
      stop 'error in fsc_cont3d: invalid dels_d'
    endif
    !
    !--------------------------------radiative transfer---------------------
    !
    call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
    dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
    !
    int3d(ii,jj,kk) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p
    !
    alocont_o_nn3d(ii,jj,kk,q14) = alo_p
    alocont_o_nn3d(iim1,jj,kk,q15) = (alo_u*c26_scontu + abs_sc*c26_intu*alocont_o_nn3d(iim1,jj,kk,q14))
    alocont_o_nn3d(ii,jjm1,kk,q17) = (alo_u*c24_scontu + abs_sc*c24_intu*alocont_o_nn3d(ii,jjm1,kk,q14))
    alocont_o_nn3d(iim1,jjm1,kk,q18) = (alo_u*c23_scontu + abs_sc*(c23_intu*alocont_o_nn3d(iim1,jjm1,kk,q14) + &
    c24_intu*alocont_o_nn3d(iim1,jjm1,kk,q15) + &
    c26_intu*alocont_o_nn3d(iim1,jjm1,kk,q17)))
    alocont_o_nn3d(ii,jj,kkm1,q23) = (alo_u*c18_scontu + abs_sc*c18_intu*alocont_o_nn3d(ii,jj,kkm1,q14))
    alocont_o_nn3d(iim1,jj,kkm1,q24) = (alo_u*c17_scontu + abs_sc*(c17_intu*alocont_o_nn3d(iim1,jj,kkm1,q14) + &
    c18_intu*alocont_o_nn3d(iim1,jj,kkm1,q15) + &
    c26_intu*alocont_o_nn3d(iim1,jj,kkm1,q23)))
    alocont_o_nn3d(ii,jjm1,kkm1,q26) = (alo_u*c15_scontu + abs_sc*(c15_intu*alocont_o_nn3d(ii,jjm1,kkm1,q14) + &
    c18_intu*alocont_o_nn3d(ii,jjm1,kkm1,q17) + &
    c24_intu*alocont_o_nn3d(ii,jjm1,kkm1,q23)))
    alocont_o_nn3d(iim1,jjm1,kkm1,q27) = (alo_u*c14_scontu + abs_sc*(c14_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q14) + &
    c15_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q15) + &
    c17_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q17) + &
    c18_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q18) + &
    c23_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q23) + &
    c24_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q24) + &
    c26_intu*alocont_o_nn3d(iim1,jjm1,kkm1,q26)))
    !
    !
    !
    !
  end subroutine fsc_cont3d_lin
  !
  !***********************************************************************
  !***********************************************************************
  !
  !                 LONG CHARACTERISTICS FOR A CELL
  !
  !***********************************************************************
  !***********************************************************************
  !
  subroutine flc_cont3d(oindx, nueindx, xindxm1, xindx, xindxp1, yindxm1, yindx, yindxp1, zindxm1, zindx, zindxp1, &
                        startxb, startyb, endxb, endyb, &
                        alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
                        q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, &
                        q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, &
                        q25, q26, q27)
    !
    !-----------------------------------------------------------------------
    !------long characteristics for continuum radiative transfer in 3d------
    !------------------at any given point x_p, y_p, z_p---------------------
    !-----------calculating intensty for given mu,phi specified-------------
    !---------------------------by input oindx------------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, opac3d, scont3d, int3d, alocont_o_nn3d, &
         fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
         mint3d_tmp, alocont_nn3d_tmp, normalization3d_tmp, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp
    use mod_angles, only: n_x, n_y, n_z
    use mod_frequencies, only: nodes_nue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: oindx, nueindx, &
    xindxm1, xindx, xindxp1, &
    yindxm1, yindx, yindxp1, &
    zindxm1, zindx, zindxp1, &
    startxb, startyb, endxb, endyb, &
    alpha, beta, gamma, &
    q1, q2, q3, q4, q5, q6, q7, q8, q9, &
    q10, q11, q12, q13, q14, q15, q16, q17, q18, &
    q19, q20, q21, q22, q23, q24, q25, q26, q27
    real(dp), intent(in) :: nn_x, nn_y, nn_z, wall
    !
    ! ... local scalars
    integer(i4b) :: nmax, nslab
    integer(i4b) :: i, j, ii, jj, kk, iim1, jjm1, kkm1, iip1, jjp1, kkp1, iid, jjd
    integer(i1b) :: ixz, iyz, ixy, ixz1, ixz2, ixz3, ixz4, iyz1, iyz2, iyz3, iyz4
    real(dp) :: opac_u, opac_p, opac_d, &
    scont_u, scont_p, scont_d, &
    x_u, x_p, x_d, &
    y_u, y_p, y_d, &
    z_u, z_p, z_d, &
    alo_u, alo_p, alo_d, &
    dels_u, dels_xyu, dels_xzu, dels_yzu, &
    dels_d, dels_xyd, dels_xzd, dels_yzd
    real(dp) :: int_u, int_sc, abs_sc, contr_sc, abs_lc, contr_lc
    real(dp) :: c14_intu, c15_intu, c17_intu, c18_intu, c23_intu, c24_intu, c26_intu, &
    c14_scontu, c15_scontu, c17_scontu, c18_scontu, c23_scontu, c24_scontu, c26_scontu
    real(dp) :: c15_scontd, c17_scontd, c18_scontd, c23_scontd, c24_scontd, c26_scontd, c27_scontd
    real(dp) :: c14_scontu1, c15_scontu1, c17_scontu1, c18_scontu1, c23_scontu1, c24_scontu1, c26_scontu1, &
    c15_scontu2, c17_scontu2, c18_scontu2, c24_scontu2, c26_scontu2, &
    c18_scontu3
    real(dp) :: a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, &
    c1, c2, c3, c4, c5, d1, d2, d3, d4, d5, &
    alo_u1, alo_u2, alo_u3, alo_p1, alo_d1, abs_scc, &
    iu10, iu11, iu12, iu13, iu14, iu15, iu16, &
    iu17, iu18, iu19, iu20, iu21, iu22, iu23, &
    iu24, iu25, iu26, iu27
    !
    ! ... local logicals
    logical :: lboundz
    !
    ! ... local characters
    !
    !define total abosorption and source contribution
    abs_lc=one
    contr_lc=zero
    !
    !define logical to check if previous z-layer is hit
    lboundz=.false.
    !
    !
    !set local indices
    ii=xindx
    jj=yindx
    kk=zindx
    iip1=ii+alpha
    jjp1=jj+beta
    kkp1=kk+gamma
    !
    !-------------------------local point----------------------------------
    !
    x_p = x(ii)
    y_p = y(jj)
    z_p = z(kk)
    scont_p = scont3d(ii,jj,kk)
    opac_p = opac3d(ii,jj,kk)
    !
    !---------------------------downwind point------------------------------
    !
    !set downwind point to opposite side if outside of boundary
    if(iip1.gt.nx) then
      iid = 1
      iip1 = iid+alpha
    elseif(iip1.lt.1) then
      iid = nx
      iip1 = iid+alpha
    else
      iid = ii
    endif
    !
    if(jjp1.gt.ny) then
      jjd = 1
      jjp1 = jjd+beta
    elseif(jjp1.lt.1) then
      jjd = ny
      jjp1 = jjd+beta
    else
      jjd = jj
    endif
    !
    !calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
    dels_xyd=(z(kkp1)-z(kk))/nn_z
    dels_xzd=(y(jjp1)-y(jjd))/nn_y
    dels_yzd=(x(iip1)-x(iid))/nn_x
    dels_d=min(dels_xyd,dels_xzd,dels_yzd)
    !
    !
    !
    if(dels_xyd.eq.dels_d) then
      !intersection with x-y plane on level k+gamma
      x_d = x(iid) + dels_d*nn_x
      y_d = y(jjd) + dels_d*nn_y
      z_d = z(kkp1)
      !
      call coeff2d_contd_lin(opac3d(iid,jjd,kkp1), opac3d(iip1,jjd,kkp1), &
      opac3d(iid,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(iid,jjd,kkp1), scont3d(iip1,jjd,kkp1), &
      scont3d(iid,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
      x(iid), x(iip1), y(jjd), y(jjp1), x_d, y_d, &
      c23_scontd, c24_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontd = zero
      c17_scontd = zero
      c18_scontd = zero
      !
    elseif(dels_xzd.eq.dels_d) then
      !intersection with x-z plane on level j+beta
      x_d = x(iid) + dels_d*nn_x
      y_d = y(jjp1)
      z_d = z(kk) + dels_d*nn_z
      !
      call coeff2d_contd_lin(opac3d(iid,jjp1,kk), opac3d(iip1,jjp1,kk), &
      opac3d(iid,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(iid,jjp1,kk), scont3d(iip1,jjp1,kk), &
      scont3d(iid,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
      x(iid), x(iip1), z(kk), z(kkp1), x_d, z_d, &
      c17_scontd, c18_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontd = zero
      c23_scontd = zero
      c24_scontd = zero
      !
    elseif(dels_yzd.eq.dels_d) then
      !intersection with y-z plane on level i+alpha
      x_d = x(iip1)
      y_d = y(jjd) + dels_d*nn_y
      z_d = z(kk) + dels_d*nn_z
      !
      call coeff2d_contd_lin(opac3d(iip1,jjd,kk), opac3d(iip1,jjp1,kk), &
      opac3d(iip1,jjd,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(iip1,jjd,kk), scont3d(iip1,jjp1,kk), &
      scont3d(iip1,jjd,kkp1), scont3d(iip1,jjp1,kkp1), &
      y(jjd), y(jjp1), z(kk), z(kkp1), y_d, z_d, &
      c15_scontd, c18_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c17_scontd = zero
      c23_scontd = zero
      c26_scontd = zero
    endif
    !
    !-----------------------------------------------------------------------
    !
    !maximum number of loops through the slab
    nslab=16
    !
    !maximum number of grid points until a boundary is hit
    nmax=nslab*2*max(nx,ny)
    !
    !for benchmark02
    !open(1, file='TRASH/benchmark02_intersections.dat', form='formatted')
    !
    !--------------default values to calculate alo-coefficients-------------
    !
    a1=one
    a2=one
    a3=one
    a4=one
    a5=one
    b1=zero
    b2=zero
    b3=zero
    b4=zero
    b5=zero
    c1=zero
    c2=zero
    c3=zero
    c4=zero
    c5=zero
    d1=zero
    d2=zero
    d3=zero
    d4=zero
    d5=zero
    !
    alo_u1=zero
    alo_u2=zero
    alo_u3=zero
    alo_p1=zero
    alo_d1=zero
    !
    ixz1=0
    ixz2=0
    ixz3=0
    ixz4=0
    iyz1=0
    iyz2=0
    iyz3=0
    iyz4=0
    !
    c14_scontu1 = zero
    c15_scontu1 = zero
    c17_scontu1 = zero
    c18_scontu1 = zero
    c23_scontu1 = zero
    c24_scontu1 = zero
    c26_scontu1 = zero
    c15_scontu2 = zero
    c17_scontu2 = zero
    c18_scontu2 = zero
    c26_scontu2 = zero
    c18_scontu3 = zero
    !
    !-----------------------------------------------------------------------
    !
    do i=1, nmax

      if(ii.eq.endxb) then
        !reset grid position if at boundary
        ii=startxb
        x_p=x(ii)
      endif
      !reset grid position if at boundary
      if(jj.eq.endyb) then
        jj=startyb
        y_p=y(jj)
      endif
      !
      iim1=ii-alpha
      jjm1=jj-beta
      kkm1=kk-gamma
      !
      !calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
      dels_xyu=(z_p-z(kkm1))/nn_z
      dels_xzu=(y_p-y(jjm1))/nn_y
      dels_yzu=(x_p-x(iim1))/nn_x
      dels_u=min(dels_xyu,dels_xzu,dels_yzu)

      !for benchmark02
      !   write(1,*) x_p, y_p, z_p
      !
      !----------------------------local point--------------------------------
      !
      !has been set at beginning of complete procedure, and will be updated
      !using the upwind point
      !
      !----------------------------upwind point-------------------------------
      !
      if(dels_xyu.eq.dels_u) then
        !intersection with x-y plane on level k-gamma
        x_u = x_p - dels_u*nn_x
        y_u = y_p - dels_u*nn_y
        z_u = z(kkm1)
        !
        call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(ii,jjm1,kkm1), &
        opac3d(iim1,jj,kkm1), opac3d(ii,jj,kkm1), &
        scont3d(iim1,jjm1,kkm1), scont3d(ii,jjm1,kkm1), &
        scont3d(iim1,jj,kkm1), scont3d(ii,jj,kkm1), &
        int3d(iim1,jjm1,kkm1), int3d(ii,jjm1,kkm1), &
        int3d(iim1,jj,kkm1), int3d(ii,jj,kkm1), &
        x(iim1), x(ii), y(jjm1), y(jj), x_u, y_u, &
        c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
        c14_intu, c15_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
        kk=kkm1
        lboundz=.true.
        !
        !set interpolation coefficients that are not used to zero
        c23_scontu = zero
        c24_scontu = zero
        c26_scontu = zero
        c23_intu = zero
        c24_intu = zero
        c26_intu = zero
        !
        !set integer values for intersection plane
        ixy=1
        ixz=0
        iyz=0
        !
      elseif(dels_xzu.eq.dels_u) then
        !intersection with x-z plane on level j-beta
        x_u = x_p - dels_u*nn_x
        y_u = y(jjm1)
        z_u = z_p - dels_u*nn_z
        !
        call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(ii,jjm1,kkm1), &
        opac3d(iim1,jjm1,kk), opac3d(ii,jjm1,kk), &
        scont3d(iim1,jjm1,kkm1), scont3d(ii,jjm1,kkm1), &
        scont3d(iim1,jjm1,kk), scont3d(ii,jjm1,kk), &
        int3d(iim1,jjm1,kkm1), int3d(ii,jjm1,kkm1), &
        int3d(iim1,jjm1,kk), int3d(ii,jjm1,kk), &
        x(iim1), x(ii), z(kkm1), z(kk), x_u, z_u, &
        c14_scontu, c15_scontu, c23_scontu, c24_scontu, &
        c14_intu, c15_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
        jj=jjm1
        !
        !set interpolation coefficients that are not used to zero
        c17_scontu = zero
        c18_scontu = zero
        c26_scontu = zero
        c17_intu = zero
        c18_intu = zero
        c26_intu = zero

        !set integer values for intersection plane (iyz=1 describing previous y-level)
        ixy=0
        ixz=0
        iyz=1
        !
      elseif(dels_yzu.eq.dels_u) then
        !intersection with y-z plane on level i-alpha
        x_u = x(iim1)
        y_u = y_p - dels_u*nn_y
        z_u = z_p - dels_u*nn_z
        !
        call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,jj,kkm1), &
        opac3d(iim1,jjm1,kk), opac3d(iim1,jj,kk), &
        scont3d(iim1,jjm1,kkm1), scont3d(iim1,jj,kkm1), &
        scont3d(iim1,jjm1,kk), scont3d(iim1,jj,kk), &
        int3d(iim1,jjm1,kkm1), int3d(iim1,jj,kkm1), &
        int3d(iim1,jjm1,kk), int3d(iim1,jj,kk), &
        y(jjm1), y(jj), z(kkm1), z(kk), y_u, z_u, &
        c14_scontu, c17_scontu, c23_scontu, c26_scontu, &
        c14_intu, c17_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
        ii=iim1
        !
        !set interpolation coefficients that are not used to zero
        c15_scontu = zero
        c18_scontu = zero
        c24_scontu = zero
        c15_intu = zero
        c18_intu = zero
        c24_intu = zero
        !
        !set integer values for intersection plane (ixz=1 describing previous x-level)
        ixy=0
        ixz=1
        iyz=0
        !
      else
        write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
        stop 'error in flc_cont3d: invalid dels_u'
      endif
      !
      !--------------------------------radiative transfer---------------------
      !
      !calculate absorption and source contribution
      call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
      dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
      !   call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
      !                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
      !   alo_d=zero
      !add up source contribution and absorption (first source-contribution, then absorption!!!!)
      contr_lc = contr_lc + abs_lc*contr_sc
      abs_lc = abs_lc*abs_sc
      !
      !set alo-coefficients
      if(i.eq.1) then
        a1 = abs_sc
        b1 = alo_u
        c1 = alo_p
        d1 = alo_d
        ixz1 = ixz
        iyz1 = iyz
        c14_scontu1 = c14_scontu
        c15_scontu1 = c15_scontu
        c17_scontu1 = c17_scontu
        c18_scontu1 = c18_scontu
        c23_scontu1 = c23_scontu
        c24_scontu1 = c24_scontu
        c26_scontu1 = c26_scontu
        iu27 = c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q14) + &
        c15_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q15) + &
        c17_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q17) + &
        c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q18)
        iu26 = c14_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q13) + &
        c15_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q14) + &
        c17_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q16) + &
        c18_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q17)
        iu25 = c15_intu * alocont_o_nn3d(xindxp1,yindxm1,zindxm1,q13) + &
        c18_intu * alocont_o_nn3d(xindxp1,yindxm1,zindxm1,q16)
        iu24 = c14_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q11) + &
        c15_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q12) + &
        c17_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q14) + &
        c18_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q15)
        iu23 = c14_intu * alocont_o_nn3d(xindx,yindx,zindxm1,q10) + &
        c15_intu * alocont_o_nn3d(xindx,yindx,zindxm1,q11) + &
        c17_intu * alocont_o_nn3d(xindx,yindx,zindxm1,q13) + &
        c18_intu * alocont_o_nn3d(xindx,yindx,zindxm1,q14)
        iu22 = c15_intu * alocont_o_nn3d(xindxp1,yindx,zindxm1,q10) + &
        c18_intu * alocont_o_nn3d(xindxp1,yindx,zindxm1,q13)
        iu21 = c17_intu * alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q11) + &
        c18_intu * alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q12)
        iu20 = c17_intu * alocont_o_nn3d(xindx,yindxp1,zindxm1,q10) + &
        c18_intu * alocont_o_nn3d(xindx,yindxp1,zindxm1,q11)
        iu19 = c18_intu * alocont_o_nn3d(xindxp1,yindxp1,zindxm1,q10)
        iu18 = c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q5) + &
        c15_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q6) + &
        c17_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q8) + &
        c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q9)
        iu17 = c14_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q4) + &
        c15_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q5) + &
        c17_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q7) + &
        c18_intu * alocont_o_nn3d(xindxp1,yindxm1,zindx,q8)
        iu16 = c15_intu * alocont_o_nn3d(xindxp1,yindxm1,zindx,q4) + &
        c18_intu * alocont_o_nn3d(xindxp1,yindxm1,zindx,q7)
        iu15 = c14_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q2) + &
        c15_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q3) + &
        c17_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q5) + &
        c18_intu * alocont_o_nn3d(xindx,yindx,zindx,q6)
        iu14 = c14_intu * alocont_o_nn3d(xindx,yindx,zindx,q1) + &
        c15_intu * alocont_o_nn3d(xindx,yindx,zindx,q2) + &
        c17_intu * alocont_o_nn3d(xindx,yindx,zindx,q4) + &
        c18_intu * alocont_o_nn3d(xindx,yindx,zindx,q5)
        iu13 = c15_intu * alocont_o_nn3d(xindxp1,yindx,zindx,q1) + &
        c18_intu * alocont_o_nn3d(xindxp1,yindx,zindx,q4)
        iu12 = c17_intu * alocont_o_nn3d(xindxm1,yindxp1,zindx,q2) + &
        c18_intu * alocont_o_nn3d(xindxm1,yindxp1,zindx,q3)
        iu11 = c17_intu * alocont_o_nn3d(xindx,yindxp1,zindx,q1) + &
        c18_intu * alocont_o_nn3d(xindx,yindxp1,zindx,q2)
        iu10 = c18_intu * alocont_o_nn3d(xindxp1,yindxp1,zindx,q1)
      elseif(i.eq.2) then
        a2 = abs_sc
        b2 = alo_u
        c2 = alo_p
        d2 = alo_d
        ixz2 = ixz
        iyz2 = iyz
        c15_scontu2 = c15_scontu
        c17_scontu2 = c17_scontu
        c18_scontu2 = c18_scontu
        c24_scontu2 = c24_scontu
        c26_scontu2 = c26_scontu
        iu27 = ixz1*c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q13) + &
        iyz1*c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q11) + &
        ixz1*c17_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q16) + &
        iyz1*c15_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q12) + &
        ixz1*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q17) + &
        iyz1*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q15) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q14)
        iu26 = iyz1*c14_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q10) + &
        iyz1*c15_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q11) + &
        ixz1*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q16) + &
        iyz1*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q14) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindx,yindxm1,zindxm1,q13)
        iu25 = iyz1*c15_intu * alocont_o_nn3d(xindxp1,yindxm1,zindxm1,q10) + &
        iyz1*c18_intu * alocont_o_nn3d(xindxp1,yindxm1,zindxm1,q13)
        iu24 = ixz1*c14_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q10) + &
        ixz1*c17_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q13) + &
        ixz1*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q14) + &
        iyz1*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q12) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindxm1,yindx,zindxm1,q11)
        iu23 = ixz1*c18_intu * alocont_o_nn3d(xindx,yindx,zindxm1,q13) + &
        iyz1*c18_intu * alocont_o_nn3d(xindx,yindx,zindxm1,q11) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindx,yindx,zindxm1,q10)
        iu22 = iyz1*c18_intu * alocont_o_nn3d(xindxp1,yindx,zindxm1,q10)
        iu21 = ixz1*c18_intu * alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q11) + &
        ixz1*c17_intu * alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q10)
        !      if(xindx.eq.nx-1.and.yindx.eq.11.and.zindx.eq.18) then
        !         write(*,*) 't1', iu21, q11, q10
        !         write(*,*) 't1', int3d(nx,11,17), int3d(1,11,17), int3d(2,11,17)
        !         write(*,*) 't1', alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q11), alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q10)
        !      endif
        iu20 = ixz1*c18_intu * alocont_o_nn3d(xindx,yindxp1,zindxm1,q10)
        iu19 = zero
        iu18 = ixz1*c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q4) + &
        iyz1*c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q2) + &
        ixz1*c17_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q7) + &
        iyz1*c15_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q3) + &
        ixz1*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q8) + &
        iyz1*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q6) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindxm1,yindxm1,zindx,q5)
        iu17 = iyz1*c14_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q1) + &
        iyz1*c15_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q2) + &
        ixz1*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q7) + &
        iyz1*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q5) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindx,yindxm1,zindx,q4)
        iu16 = iyz1*c15_intu * alocont_o_nn3d(xindxp1,yindxm1,zindx,q1) + &
        iyz1*c18_intu * alocont_o_nn3d(xindxp1,yindxm1,zindx,q4)
        iu15 = ixz1*c14_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q1) + &
        ixz1*c17_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q4) + &
        ixz1*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q5) + &
        iyz1*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q3) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindxm1,yindx,zindx,q2)
        iu14 = ixz1*c18_intu * alocont_o_nn3d(xindx,yindx,zindx,q4) + &
        iyz1*c18_intu * alocont_o_nn3d(xindx,yindx,zindx,q2) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindx,yindx,zindx,q1)
        iu13 = iyz1*c18_intu * alocont_o_nn3d(xindxp1,yindx,zindx,q1)
        !      if(xindx.eq.nx.and.yindx.eq.ny.and.zindx.eq.17) then
        !         write(*,*) 't1', iu13, iyz1, c18_intu, alocont_o_nn3d(xindxp1,yindx,zindx,q1), int3d(nx-1,ny-1,16)
        !      endif
        iu12 = ixz1*c17_intu * alocont_o_nn3d(xindxm1,yindxp1,zindx,q1) + &
        ixz1*c18_intu * alocont_o_nn3d(xindxm1,yindxp1,zindx,q2)
        iu11 = ixz1*c18_intu * alocont_o_nn3d(xindx,yindxp1,zindx,q1)
        iu10 = zero
        !      if(xindx.eq.1.and.yindx.eq.1.and.zindx.eq.17) then
        !         write(*,*) 't1', iu13, iyz1, c18_intu, alocont_o_nn3d(xindxp1,yindx,zindx,q1), xindxp1, yindx, zindx, q1
        !      endif
      elseif(i.eq.3) then
        a3 = abs_sc
        b3 = alo_u
        c3 = alo_p
        d3 = alo_d
        ixz3 = ixz
        iyz3 = iyz
        c18_scontu3 = c18_scontu
        iu27 = (ixz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q16) + &
        (iyz1*iyz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q12) + &
        (ixz1*iyz2+iyz1*ixz2)*c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q10) + &
        (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q14) + &
        ((ixz1*iyz2+iyz1*ixz2)*c17_intu + (ixz1*ixz2)*c15_intu)* alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q13) + &
        ((ixz1*iyz2+iyz1*ixz2)*c15_intu + (iyz1*iyz2)*c17_intu)* alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q11)
        iu26 = (iyz1*iyz2)*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q11) + &
        (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q13) + &
        ((ixz1*iyz2+iyz1*ixz2)*c15_intu + (iyz1*iyz2)*c17_intu) * alocont_o_nn3d(xindx,yindxm1,zindxm1,q10)
        iu25 = (iyz1*iyz2)*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q10)
        iu24 = (ixz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q13) + &
        (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q11) + &
        ((ixz1*iyz2+iyz1*ixz2)*c17_intu + (ixz1*ixz2)*c15_intu)* alocont_o_nn3d(xindxm1,yindx,zindxm1,q10)
        iu23 = (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindx,yindx,zindxm1,q10)
        iu22 = zero
        iu21 = (ixz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q10)
        iu20 = zero
        iu18 = (ixz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q7) + &
        (iyz1*iyz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q3) + &
        (ixz1*iyz2+iyz1*ixz2)*c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q1) + &
        (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q5) + &
        ((ixz1*iyz2+iyz1*ixz2)*c17_intu + (ixz1*ixz2)*c15_intu)* alocont_o_nn3d(xindxm1,yindxm1,zindx,q4) + &
        ((ixz1*iyz2+iyz1*ixz2)*c15_intu + (iyz1*iyz2)*c17_intu)* alocont_o_nn3d(xindxm1,yindxm1,zindx,q2)
        iu17 = (iyz1*iyz2)*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q2) + &
        (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q4) + &
        ((ixz1*iyz2+iyz1*ixz2)*c15_intu + (iyz1*iyz2)*c17_intu)* alocont_o_nn3d(xindx,yindxm1,zindx,q1)
        iu16 = (iyz1*iyz2)*c18_intu * alocont_o_nn3d(xindxp1,yindxm1,zindx,q1)
        iu15 = (ixz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q4) + &
        (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q2) + &
        ((ixz1*iyz2+iyz1*ixz2)*c17_intu + (ixz1*ixz2)*c15_intu)* alocont_o_nn3d(xindxm1,yindx,zindx,q1)
        iu14 = (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindx,yindx,zindx,q1)
        iu13 = zero
        iu12 = (ixz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxp1,zindx,q1)
        iu11 = zero
      elseif(i.eq.4) then
        a4 = abs_sc
        b4 = alo_u
        c4 = alo_p
        d4 = alo_d
        ixz4 = ixz
        iyz4 = iyz
        iu27 = (ixz1*ixz2*iyz3 + ixz1*iyz2*ixz3 + iyz1*ixz2*ixz3)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q13) + &
        (ixz1*iyz2*iyz3 + iyz1*iyz2*ixz3 + iyz1*ixz2*iyz3)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q11) + &
        ((ixz1*ixz2*iyz3 + ixz1*iyz2*ixz3 + iyz1*ixz2*ixz3)*c15_intu + &
        (ixz1*iyz2*iyz3 + iyz1*iyz2*ixz3 + iyz1*ixz2*iyz3)*c17_intu) * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q10)
        iu26 = ((ixz1*iyz2*iyz3 + iyz1*iyz2*ixz3 + iyz1*ixz2*iyz3)*c18_intu) * alocont_o_nn3d(xindx,yindxm1,zindxm1,q10)
        iu25 = zero
        iu24 = (ixz1*ixz2*iyz3 + ixz1*iyz2*ixz3 + iyz1*ixz2*ixz3)*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q10)
        iu23 = zero
        iu21 = zero
        iu18 = (ixz1*ixz2*iyz3 + ixz1*iyz2*ixz3 + iyz1*ixz2*ixz3)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q4) + &
        (ixz1*iyz2*iyz3 + iyz1*iyz2*ixz3 + iyz1*ixz2*iyz3)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindx,q2) + &
        ((ixz1*ixz2*iyz3 + ixz1*iyz2*ixz3 + iyz1*ixz2*ixz3)*c15_intu + &
        (ixz1*iyz2*iyz3 + iyz1*iyz2*ixz3 + iyz1*ixz2*iyz3)*c17_intu) * alocont_o_nn3d(xindxm1,yindxm1,zindx,q1)
        iu17 = (ixz1*iyz2*iyz3 + iyz1*iyz2*ixz3 + iyz1*ixz2*iyz3)*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindx,q1)
        iu16 = zero
        iu15 = (ixz1*ixz2*iyz3 + ixz1*iyz2*ixz3 + iyz1*ixz2*ixz3)*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindx,q1)
        iu14 = zero
        iu12 = zero
      elseif(i.eq.5) then
        a5 = abs_sc
        b5 = alo_u
        c5 = alo_p
        d5 = alo_d
        iu27 = ((ixz1*ixz2*iyz3*iyz4) + (iyz1*iyz2*ixz1*ixz2) + &
        (ixz1*iyz2*iyz3*ixz4) + (iyz1*ixz2*ixz3*iyz4) + &
        (ixz1*iyz2*ixz3*iyz4) + (iyz1*ixz2*iyz3*ixz4)) * c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q1)
        iu26 = zero
        iu24 = zero
        iu18 = zero
        iu17 = zero
        iu15 = zero
      elseif(i.eq.6) then
        iu27 = zero
      endif
      !
      !update the downwind and current point for next step
      opac_d=opac_p
      dels_d=dels_u
      scont_d=scont_p
      x_p=x_u
      y_p=y_u
      z_p=z_u
      opac_p=opac_u
      scont_p=scont_u
      !   write(*,*) 't1', i, ixz1, ixz2, ixz3, iyz1, iyz2, iyz3
      !exit z-loop if z-layer is hit
      if(lboundz) exit
      !
    enddo
    !
    if(.not.lboundz) stop 'error in flc_cont3d: no previous z-layer has been found, increase nslab'
    !
    !for benchmark02
    !close(1)
    !
    !
    !
    !set the intensity + alo
    int3d(xindx,yindx,zindx) = int_u*abs_lc + contr_lc
    !
    alo_u1 = b1+a1*c2+a1*a2*d3
    alo_p1 = c1+a1*d2
    alo_d1 = d1
    alo_u2 = a1*b2+a1*a2*c3+a1*a2*a3*d4
    alo_u3 = a1*a2*b3 + a1*a2*a3*c4 + a1*a2*a3*a4*d5
    abs_scc = a1*a2*a3*a4*a5
    !
    alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q27) = alo_u1*c14_scontu1 + alo_u2*(c15_scontu2*ixz1+c17_scontu2*iyz1) + &
    alo_u3*c18_scontu3*(ixz1*iyz2+iyz1*ixz2) + abs_scc*iu27
    alocont_o_nn3d(xindx,yindxm1,zindxm1,q26) = alo_u1*c15_scontu1 + alo_u2*c18_scontu2*iyz1 + abs_scc*iu26
    alocont_o_nn3d(xindxp1,yindxm1,zindxm1,q25) = abs_scc*iu25
    alocont_o_nn3d(xindxm1,yindx,zindxm1,q24) = alo_u1*c17_scontu1 + alo_u2*c18_scontu2*ixz1 + abs_scc*iu24
    alocont_o_nn3d(xindx,yindx,zindxm1,q23) = alo_u1*c18_scontu1 + abs_scc*iu23
    alocont_o_nn3d(xindxp1,yindx,zindxm1,q22) = abs_scc*iu22
    alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q21) = abs_scc*iu21
    alocont_o_nn3d(xindx,yindxp1,zindxm1,q20) = abs_scc*iu20
    alocont_o_nn3d(xindxp1,yindxp1,zindxm1,q19) = abs_scc*iu19
    alocont_o_nn3d(xindxm1,yindxm1,zindx,q18) = alo_u1*c23_scontu1 + alo_u2*(iyz1*c26_scontu2+ixz1*c24_scontu2) + abs_scc*iu18
    alocont_o_nn3d(xindx,yindxm1,zindx,q17) = alo_u1*c24_scontu1 + abs_scc*iu17
    alocont_o_nn3d(xindxp1,yindxm1,zindx,q16) = abs_scc*iu16
    alocont_o_nn3d(xindxm1,yindx,zindx,q15) = alo_u1*c26_scontu1 + abs_scc*iu15
    alocont_o_nn3d(xindx,yindx,zindx,q14) = alo_p1 + abs_scc*iu14
    alocont_o_nn3d(xindxp1,yindx,zindx,q13) = alo_d1*c15_scontd + abs_scc*iu13
    alocont_o_nn3d(xindxm1,yindxp1,zindx,q12) = abs_scc*iu12
    alocont_o_nn3d(xindx,yindxp1,zindx,q11) = alo_d1*c17_scontd + abs_scc*iu11
    alocont_o_nn3d(xindxp1,yindxp1,zindx,q10) = alo_d1*c18_scontd + abs_scc*iu10
    alocont_o_nn3d(xindxm1,yindxm1,zindxp1,q9) = zero
    alocont_o_nn3d(xindx,yindxm1,zindxp1,q8) = zero
    alocont_o_nn3d(xindxp1,yindxm1,zindxp1,q7) = zero
    alocont_o_nn3d(xindxm1,yindx,zindxp1,q6) = zero
    alocont_o_nn3d(xindx,yindx,zindxp1,q5) = alo_d1*c23_scontd
    alocont_o_nn3d(xindxp1,yindx,zindxp1,q4) = alo_d1*c24_scontd
    alocont_o_nn3d(xindxm1,yindxp1,zindxp1,q3) = zero
    alocont_o_nn3d(xindx,yindxp1,zindxp1,q2) = alo_d1*c26_scontd
    alocont_o_nn3d(xindxp1,yindxp1,zindxp1,q1) = alo_d1*c27_scontd
    !
    !if(xindx.eq.2.and.yindx.eq.2.and.zindx.eq.16) then
    !   write(*,*) 't2', q18, scont3d(xindxm1,yindxm1,zindx), alocont_o_nn3d(xindxm1,yindxm1,zindx,q18), int3d(xindx,yindx,zindx), iu18, int_u, &
    !        alo_u1*c23_scontu1 + alo_u2*(iyz1*c26_scontu2+ixz1*c24_scontu2), contr_lc, iu18, int_u, alo_u1*c23_scontu1, alo_u2*(iyz1*c26_scontu2+ixz1*c24_scontu2), alo_u1, contr_sc
    !!   , int3d(nx-1,yindx-1,zindx+1), alocont_o_nn3d(nx,yindx,zindx,19), alocont_o_nn3d(xindx,yindx,zindx,q1)!
    !endif

    !
    !mint3d_tmp(xindx,yindx,zindx) = mint3d_tmp(xindx,yindx,zindx) + int3d(xindx,yindx,zindx)*wall
    !fcontx3d_tmp(xindx,yindx,zindx) = fcontx3d_tmp(xindx,yindx,zindx) + int3d(xindx,yindx,zindx)*nn_x*wall
    !fconty3d_tmp(xindx,yindx,zindx) = fconty3d_tmp(xindx,yindx,zindx) + int3d(xindx,yindx,zindx)*nn_y*wall
    !fcontz3d_tmp(xindx,yindx,zindx) = fcontz3d_tmp(xindx,yindx,zindx) + int3d(xindx,yindx,zindx)*nn_z*wall
    !normalization3d_tmp(xindx,yindx,zindx) = normalization3d_tmp(xindx,yindx,zindx) + wall
    !
    !alocont_nn3d_tmp(xindxp1,yindxp1,zindxp1,q1) = alocont_nn3d_tmp(xindxp1,yindxp1,zindxp1,q1) + wall*alocont_o_nn3d(xindxp1,yindxp1,zindxp1,q1)
    !alocont_nn3d_tmp(xindx,yindxp1,zindxp1,q2) = alocont_nn3d_tmp(xindx,yindxp1,zindxp1,q2) + wall*alocont_o_nn3d(xindx,yindxp1,zindxp1,q2)
    !alocont_nn3d_tmp(xindxm1,yindxp1,zindxp1,q3) = alocont_nn3d_tmp(xindxm1,yindxp1,zindxp1,q3) + wall*alocont_o_nn3d(xindxm1,yindxp1,zindxp1,q3)
    !alocont_nn3d_tmp(xindxp1,yindx,zindxp1,q4) = alocont_nn3d_tmp(xindxp1,yindx,zindxp1,q4) + wall*alocont_o_nn3d(xindxp1,yindx,zindxp1,q4)
    !alocont_nn3d_tmp(xindx,yindx,zindxp1,q5) = alocont_nn3d_tmp(xindx,yindx,zindxp1,q5) + wall*alocont_o_nn3d(xindx,yindx,zindxp1,q5)
    !alocont_nn3d_tmp(xindxm1,yindx,zindxp1,q6) = alocont_nn3d_tmp(xindxm1,yindx,zindxp1,q6) + wall*alocont_o_nn3d(xindxm1,yindx,zindxp1,q6)
    !alocont_nn3d_tmp(xindxp1,yindxm1,zindxp1,q7) = alocont_nn3d_tmp(xindxp1,yindxm1,zindxp1,q7) + wall*alocont_o_nn3d(xindxp1,yindxm1,zindxp1,q7)
    !alocont_nn3d_tmp(xindx,yindxm1,zindxp1,q8) = alocont_nn3d_tmp(xindx,yindxm1,zindxp1,q8) + wall*alocont_o_nn3d(xindx,yindxm1,zindxp1,q8)
    !alocont_nn3d_tmp(xindxm1,yindxm1,zindxp1,q9) = alocont_nn3d_tmp(xindxm1,yindxm1,zindxp1,q9) + wall*alocont_o_nn3d(xindxm1,yindxm1,zindxp1,q9)
    !alocont_nn3d_tmp(xindxp1,yindxp1,zindx,q10) = alocont_nn3d_tmp(xindxp1,yindxp1,zindx,q10) + wall*alocont_o_nn3d(xindxp1,yindxp1,zindx,q10)
    !alocont_nn3d_tmp(xindx,yindxp1,zindx,q11) = alocont_nn3d_tmp(xindx,yindxp1,zindx,q11) + wall*alocont_o_nn3d(xindx,yindxp1,zindx,q11)
    !alocont_nn3d_tmp(xindxm1,yindxp1,zindx,q12) = alocont_nn3d_tmp(xindxm1,yindxp1,zindx,q12) + wall*alocont_o_nn3d(xindxm1,yindxp1,zindx,q12)
    !alocont_nn3d_tmp(xindxp1,yindx,zindx,q13) = alocont_nn3d_tmp(xindxp1,yindx,zindx,q13) + wall*alocont_o_nn3d(xindxp1,yindx,zindx,q13)
    !alocont_nn3d_tmp(xindx,yindx,zindx,q14) = alocont_nn3d_tmp(xindx,yindx,zindx,q14) + wall*alocont_o_nn3d(xindx,yindx,zindx,q14)
    !alocont_nn3d_tmp(xindxm1,yindx,zindx,q15) = alocont_nn3d_tmp(xindxm1,yindx,zindx,q15) + wall*alocont_o_nn3d(xindxm1,yindx,zindx,q15)
    !alocont_nn3d_tmp(xindxp1,yindxm1,zindx,q16) = alocont_nn3d_tmp(xindxp1,yindxm1,zindx,q16) + wall*alocont_o_nn3d(xindxp1,yindxm1,zindx,q16)
    !alocont_nn3d_tmp(xindx,yindxm1,zindx,q17) = alocont_nn3d_tmp(xindx,yindxm1,zindx,q17) + wall*alocont_o_nn3d(xindx,yindxm1,zindx,q17)
    !alocont_nn3d_tmp(xindxm1,yindxm1,zindx,q18) = alocont_nn3d_tmp(xindxm1,yindxm1,zindx,q18) + wall*alocont_o_nn3d(xindxm1,yindxm1,zindx,q18)
    !alocont_nn3d_tmp(xindxp1,yindxp1,zindxm1,q19) = alocont_nn3d_tmp(xindxp1,yindxp1,zindxm1,q19) + wall*alocont_o_nn3d(xindxp1,yindxp1,zindxm1,q19)
    !alocont_nn3d_tmp(xindx,yindxp1,zindxm1,q20) = alocont_nn3d_tmp(xindx,yindxp1,zindxm1,q20) + wall*alocont_o_nn3d(xindx,yindxp1,zindxm1,q20)
    !alocont_nn3d_tmp(xindxm1,yindxp1,zindxm1,q21) = alocont_nn3d_tmp(xindxm1,yindxp1,zindxm1,q21) + wall*alocont_o_nn3d(xindxm1,yindxp1,zindxm1,q21)
    !alocont_nn3d_tmp(xindxp1,yindx,zindxm1,q22) = alocont_nn3d_tmp(xindxp1,yindx,zindxm1,q22) + wall*alocont_o_nn3d(xindxp1,yindx,zindxm1,q22)
    !alocont_nn3d_tmp(xindx,yindx,zindxm1,q23) = alocont_nn3d_tmp(xindx,yindx,zindxm1,q23) + wall*alocont_o_nn3d(xindx,yindx,zindxm1,q23)
    !alocont_nn3d_tmp(xindxm1,yindx,zindxm1,q24) = alocont_nn3d_tmp(xindxm1,yindx,zindxm1,q24) + wall*alocont_o_nn3d(xindxm1,yindx,zindxm1,q24)
    !alocont_nn3d_tmp(xindxp1,yindxm1,zindxm1,q25) = alocont_nn3d_tmp(xindxp1,yindxm1,zindxm1,q25) + wall*alocont_o_nn3d(xindxp1,yindxm1,zindxm1,q25)
    !alocont_nn3d_tmp(xindx,yindxm1,zindxm1,q26) = alocont_nn3d_tmp(xindx,yindxm1,zindxm1,q26) + wall*alocont_o_nn3d(xindx,yindxm1,zindxm1,q26)
    !alocont_nn3d_tmp(xindxm1,yindxm1,zindxm1,q27) = alocont_nn3d_tmp(xindxm1,yindxm1,zindxm1,q27) + wall*alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q27)




    !if(xindx.eq.9.and.yindx.eq.10.and.zindx.eq.4) then
    !   write(*,*) xindx, yindx, zindx, alocont_o_nn3d(xindxp1,yindx,zindxp1,q4), int_u*abs_lc + contr_lc, alocont_o_nn3d(10,10,5,q4)
    !stop
    !endif

    !write(*,*) int3d(xindx,yindx,zindx), abs_lc, int_u
    !if(int3d(xindx,yindx,zindx).lt.zero) stop

    !
    !for benchmark02
    !open(1, file='TRASH/benchmark02_x.dat', form='formatted')
    !   do i=1, nx
    !      write(1,'(es20.8)') x(i)
    !   enddo
    !close(1)!
    !
    !open(1, file='TRASH/benchmark02_y.dat', form='formatted')
    !   do i=1, ny
    !      write(1,'(es20.8)') y(i)
    !   enddo
    !close(1)
    !


    !
  end subroutine flc_cont3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine flc_cont3d_lin(oindx, nueindx, xindxm1, xindx, xindxp1, yindxm1, yindx, yindxp1, zindxm1, zindx, zindxp1, &
                            startxb, startyb, endxb, endyb, &
                            alpha, beta, gamma, nn_x, nn_y, nn_z, wall, &
                            q14, q15, q17, q18, q23, q24, q26, q27)
    !
    !-----------------------------------------------------------------------
    !------long characteristics for continuum radiative transfer in 3d------
    !------------------at any given point x_p, y_p, z_p---------------------
    !-----------calculating intensty for given mu,phi specified-------------
    !---------------------------by input oindx------------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, x, y, z, opac3d, scont3d, int3d, alocont_o_nn3d, &
         fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
         mint3d_tmp, alocont_nn3d_tmp, normalization3d_tmp, &
         kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
         kcontyy3d_tmp, kcontyz3d_tmp, &
         kcontzz3d_tmp
    use mod_angles, only: n_x, n_y, n_z
    use mod_frequencies, only: nodes_nue
    !
    !
    ! ... arguments
    integer(i4b), intent(in) :: oindx, nueindx, &
    xindxm1, xindx, xindxp1, &
    yindxm1, yindx, yindxp1, &
    zindxm1, zindx, zindxp1, &
    startxb, startyb, endxb, endyb, &
    alpha, beta, gamma, &
    q14, q15, q17, q18, q23, q24, q26, q27
    real(dp), intent(in) :: nn_x, nn_y, nn_z, wall
    !
    ! ... local scalars
    integer(i4b) :: nmax, nslab
    integer(i4b) :: i, j, ii, jj, kk, iim1, jjm1, kkm1, iip1, jjp1, kkp1, iid, jjd
    integer(i1b) :: ixz, iyz, ixy, ixz1, ixz2, ixz3, ixz4, iyz1, iyz2, iyz3, iyz4
    real(dp) :: opac_u, opac_p, opac_d, &
    scont_u, scont_p, scont_d, &
    x_u, x_p, x_d, &
    y_u, y_p, y_d, &
    z_u, z_p, z_d, &
    alo_u, alo_p, alo_d, &
    dels_u, dels_xyu, dels_xzu, dels_yzu, &
    dels_d, dels_xyd, dels_xzd, dels_yzd
    real(dp) :: int_u, int_sc, abs_sc, contr_sc, abs_lc, contr_lc
    real(dp) :: c14_intu, c15_intu, c17_intu, c18_intu, c23_intu, c24_intu, c26_intu, &
    c14_scontu, c15_scontu, c17_scontu, c18_scontu, c23_scontu, c24_scontu, c26_scontu
    real(dp) :: c15_scontd, c17_scontd, c18_scontd, c23_scontd, c24_scontd, c26_scontd, c27_scontd
    real(dp) :: c14_scontu1, c15_scontu1, c17_scontu1, c18_scontu1, c23_scontu1, c24_scontu1, c26_scontu1, &
    c15_scontu2, c17_scontu2, c18_scontu2, c24_scontu2, c26_scontu2, &
    c18_scontu3
    real(dp) :: a1, a2, a3, a4, b1, b2, b3, b4, &
    c1, c2, c3, c4, &
    alo_u1, alo_u2, alo_u3, alo_p1, abs_scc, &
    iu14, iu15, iu17, iu18, iu23, iu24, iu26, iu27
    !
    ! ... local logicals
    logical :: lboundz
    !
    ! ... local characters
    !
    !define total abosorption and source contribution
    abs_lc=one
    contr_lc=zero
    !
    !define logical to check if previous z-layer is hit
    lboundz=.false.
    !
    !
    !set local indices
    ii=xindx
    jj=yindx
    kk=zindx
    iip1=ii+alpha
    jjp1=jj+beta
    kkp1=kk+gamma
    !
    !-------------------------local point----------------------------------
    !
    x_p = x(ii)
    y_p = y(jj)
    z_p = z(kk)
    scont_p = scont3d(ii,jj,kk)
    opac_p = opac3d(ii,jj,kk)
    !
    !---------------------------downwind point------------------------------
    !
    !set downwind point to opposite side if outside of boundary
    if(iip1.gt.nx) then
      iid = 1
      iip1 = iid+alpha
    elseif(iip1.lt.1) then
      iid = nx
      iip1 = iid+alpha
    else
      iid = ii
    endif
    !
    if(jjp1.gt.ny) then
      jjd = 1
      jjp1 = jjd+beta
    elseif(jjp1.lt.1) then
      jjd = ny
      jjp1 = jjd+beta
    else
      jjd = jj
    endif
    !
    !calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
    dels_xyd=(z(kkp1)-z(kk))/nn_z
    dels_xzd=(y(jjp1)-y(jjd))/nn_y
    dels_yzd=(x(iip1)-x(iid))/nn_x
    dels_d=min(dels_xyd,dels_xzd,dels_yzd)
    !
    !
    !
    if(dels_xyd.eq.dels_d) then
      !intersection with x-y plane on level k+gamma
      x_d = x(iid) + dels_d*nn_x
      y_d = y(jjd) + dels_d*nn_y
      z_d = z(kkp1)
      !
      call coeff2d_contd_lin(opac3d(iid,jjd,kkp1), opac3d(iip1,jjd,kkp1), &
      opac3d(iid,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(iid,jjd,kkp1), scont3d(iip1,jjd,kkp1), &
      scont3d(iid,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
      x(iid), x(iip1), y(jjd), y(jjp1), x_d, y_d, &
      c23_scontd, c24_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontd = zero
      c17_scontd = zero
      c18_scontd = zero
      !
    elseif(dels_xzd.eq.dels_d) then
      !intersection with x-z plane on level j+beta
      x_d = x(iid) + dels_d*nn_x
      y_d = y(jjp1)
      z_d = z(kk) + dels_d*nn_z
      !
      call coeff2d_contd_lin(opac3d(iid,jjp1,kk), opac3d(iip1,jjp1,kk), &
      opac3d(iid,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(iid,jjp1,kk), scont3d(iip1,jjp1,kk), &
      scont3d(iid,jjp1,kkp1), scont3d(iip1,jjp1,kkp1), &
      x(iid), x(iip1), z(kk), z(kkp1), x_d, z_d, &
      c17_scontd, c18_scontd, c26_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c15_scontd = zero
      c23_scontd = zero
      c24_scontd = zero
      !
    elseif(dels_yzd.eq.dels_d) then
      !intersection with y-z plane on level i+alpha
      x_d = x(iip1)
      y_d = y(jjd) + dels_d*nn_y
      z_d = z(kk) + dels_d*nn_z
      !
      call coeff2d_contd_lin(opac3d(iip1,jjd,kk), opac3d(iip1,jjp1,kk), &
      opac3d(iip1,jjd,kkp1), opac3d(iip1,jjp1,kkp1), &
      scont3d(iip1,jjd,kk), scont3d(iip1,jjp1,kk), &
      scont3d(iip1,jjd,kkp1), scont3d(iip1,jjp1,kkp1), &
      y(jjd), y(jjp1), z(kk), z(kkp1), y_d, z_d, &
      c15_scontd, c18_scontd, c24_scontd, c27_scontd, opac_d, scont_d)
      !
      !set interpolation coefficients that are not used to zero
      c17_scontd = zero
      c23_scontd = zero
      c26_scontd = zero
    endif
    !
    !-----------------------------------------------------------------------
    !
    !maximum number of loops through the slab
    nslab=16
    !
    !maximum number of grid points until a boundary is hit
    nmax=nslab*2*max(nx,ny)
    !
    !for benchmark02
    !open(1, file='TRASH/benchmark02_intersections.dat', form='formatted')
    !
    !--------------default values to calculate alo-coefficients-------------
    !
    a1=one
    a2=one
    a3=one
    a4=one
    b1=zero
    b2=zero
    b3=zero
    b4=zero
    c1=zero
    c2=zero
    c3=zero
    c4=zero
    !
    alo_u1=zero
    alo_u2=zero
    alo_u3=zero
    alo_p1=zero
    !
    ixz1=0
    ixz2=0
    ixz3=0
    ixz4=0
    iyz1=0
    iyz2=0
    iyz3=0
    iyz4=0
    !
    c14_scontu1 = zero
    c15_scontu1 = zero
    c17_scontu1 = zero
    c18_scontu1 = zero
    c23_scontu1 = zero
    c24_scontu1 = zero
    c26_scontu1 = zero
    c15_scontu2 = zero
    c17_scontu2 = zero
    c18_scontu2 = zero
    c26_scontu2 = zero
    c18_scontu3 = zero
    !
    !-----------------------------------------------------------------------
    !
    do i=1, nmax

      if(ii.eq.endxb) then
        !reset grid position if at boundary
        ii=startxb
        x_p=x(ii)
      endif
      !reset grid position if at boundary
      if(jj.eq.endyb) then
        jj=startyb
        y_p=y(jj)
      endif
      !
      iim1=ii-alpha
      jjm1=jj-beta
      kkm1=kk-gamma
      !
      !calculate distance of upwind point to (previous) xy-plane, xz-plane and yz-plane
      dels_xyu=(z_p-z(kkm1))/nn_z
      dels_xzu=(y_p-y(jjm1))/nn_y
      dels_yzu=(x_p-x(iim1))/nn_x
      dels_u=min(dels_xyu,dels_xzu,dels_yzu)

      !for benchmark02
      !   write(1,*) x_p, y_p, z_p
      !
      !----------------------------local point--------------------------------
      !
      !has been set at beginning of complete procedure, and will be updated
      !using the upwind point
      !
      !----------------------------upwind point-------------------------------
      !
      if(dels_xyu.eq.dels_u) then
        !intersection with x-y plane on level k-gamma
        x_u = x_p - dels_u*nn_x
        y_u = y_p - dels_u*nn_y
        z_u = z(kkm1)
        !
        call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(ii,jjm1,kkm1), &
        opac3d(iim1,jj,kkm1), opac3d(ii,jj,kkm1), &
        scont3d(iim1,jjm1,kkm1), scont3d(ii,jjm1,kkm1), &
        scont3d(iim1,jj,kkm1), scont3d(ii,jj,kkm1), &
        int3d(iim1,jjm1,kkm1), int3d(ii,jjm1,kkm1), &
        int3d(iim1,jj,kkm1), int3d(ii,jj,kkm1), &
        x(iim1), x(ii), y(jjm1), y(jj), x_u, y_u, &
        c14_scontu, c15_scontu, c17_scontu, c18_scontu, &
        c14_intu, c15_intu, c17_intu, c18_intu, opac_u, scont_u, int_u)
        kk=kkm1
        lboundz=.true.
        !
        !set interpolation coefficients that are not used to zero
        c23_scontu = zero
        c24_scontu = zero
        c26_scontu = zero
        c23_intu = zero
        c24_intu = zero
        c26_intu = zero
        !
        !set integer values for intersection plane
        ixy=1
        ixz=0
        iyz=0
        !
      elseif(dels_xzu.eq.dels_u) then
        !intersection with x-z plane on level j-beta
        x_u = x_p - dels_u*nn_x
        y_u = y(jjm1)
        z_u = z_p - dels_u*nn_z
        !
        call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(ii,jjm1,kkm1), &
        opac3d(iim1,jjm1,kk), opac3d(ii,jjm1,kk), &
        scont3d(iim1,jjm1,kkm1), scont3d(ii,jjm1,kkm1), &
        scont3d(iim1,jjm1,kk), scont3d(ii,jjm1,kk), &
        int3d(iim1,jjm1,kkm1), int3d(ii,jjm1,kkm1), &
        int3d(iim1,jjm1,kk), int3d(ii,jjm1,kk), &
        x(iim1), x(ii), z(kkm1), z(kk), x_u, z_u, &
        c14_scontu, c15_scontu, c23_scontu, c24_scontu, &
        c14_intu, c15_intu, c23_intu, c24_intu, opac_u, scont_u, int_u)
        jj=jjm1
        !
        !set interpolation coefficients that are not used to zero
        c17_scontu = zero
        c18_scontu = zero
        c26_scontu = zero
        c17_intu = zero
        c18_intu = zero
        c26_intu = zero

        !set integer values for intersection plane (iyz=1 describing previous y-level)
        ixy=0
        ixz=0
        iyz=1
        !
      elseif(dels_yzu.eq.dels_u) then
        !intersection with y-z plane on level i-alpha
        x_u = x(iim1)
        y_u = y_p - dels_u*nn_y
        z_u = z_p - dels_u*nn_z
        !
        call coeff2d_contu_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,jj,kkm1), &
        opac3d(iim1,jjm1,kk), opac3d(iim1,jj,kk), &
        scont3d(iim1,jjm1,kkm1), scont3d(iim1,jj,kkm1), &
        scont3d(iim1,jjm1,kk), scont3d(iim1,jj,kk), &
        int3d(iim1,jjm1,kkm1), int3d(iim1,jj,kkm1), &
        int3d(iim1,jjm1,kk), int3d(iim1,jj,kk), &
        y(jjm1), y(jj), z(kkm1), z(kk), y_u, z_u, &
        c14_scontu, c17_scontu, c23_scontu, c26_scontu, &
        c14_intu, c17_intu, c23_intu, c26_intu, opac_u, scont_u, int_u)
        ii=iim1
        !
        !set interpolation coefficients that are not used to zero
        c15_scontu = zero
        c18_scontu = zero
        c24_scontu = zero
        c15_intu = zero
        c18_intu = zero
        c24_intu = zero
        !
        !set integer values for intersection plane (ixz=1 describing previous x-level)
        ixy=0
        ixz=1
        iyz=0
        !
      else
        write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
        stop 'error in flc_cont3d: invalid dels_u'
      endif
      !
      !--------------------------------radiative transfer---------------------
      !
      !calculate absorption and source contribution
      call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
      dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
      !
      !add up source contribution and absorption (first source-contribution, then absorption!!!!)
      contr_lc = contr_lc + abs_lc*contr_sc
      abs_lc = abs_lc*abs_sc
      !
      !set alo-coefficients
      if(i.eq.1) then
        a1 = abs_sc
        b1 = alo_u
        c1 = alo_p
        ixz1 = ixz
        iyz1 = iyz
        c14_scontu1 = c14_scontu
        c15_scontu1 = c15_scontu
        c17_scontu1 = c17_scontu
        c18_scontu1 = c18_scontu
        c23_scontu1 = c23_scontu
        c24_scontu1 = c24_scontu
        c26_scontu1 = c26_scontu
        iu27 = c14_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q14) + &
        c15_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q15) + &
        c17_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q17) + &
        c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q18)
        iu26 = c15_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q14) + &
        c18_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q17)
        iu24 = c17_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q14) + &
        c18_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q15)
        iu23 = c18_intu * alocont_o_nn3d(xindx,yindx,zindxm1,q14)
        iu18 = zero
        iu17 = zero
        iu15 = zero
        iu14 = zero
      elseif(i.eq.2) then
        a2 = abs_sc
        b2 = alo_u
        c2 = alo_p
        ixz2 = ixz
        iyz2 = iyz
        c15_scontu2 = c15_scontu
        c17_scontu2 = c17_scontu
        c18_scontu2 = c18_scontu
        c24_scontu2 = c24_scontu
        c26_scontu2 = c26_scontu
        iu27 = ixz1*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q17) + &
        iyz1*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q15) + &
        (ixz1*c15_intu+iyz1*c17_intu) * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q14)
        iu26 = iyz1*c18_intu * alocont_o_nn3d(xindx,yindxm1,zindxm1,q14)
        iu24 = ixz1*c18_intu * alocont_o_nn3d(xindxm1,yindx,zindxm1,q14)
        iu23 = zero
      elseif(i.eq.3) then
        a3 = abs_sc
        b3 = alo_u
        c3 = alo_p
        ixz3 = ixz
        iyz3 = iyz
        c18_scontu3 = c18_scontu
        iu27 = (ixz1*iyz2+iyz1*ixz2)*c18_intu * alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q14)
        iu26 = zero
        iu24 = zero
        iu23 = zero
      elseif(i.eq.4) then
        a4 = abs_sc
        b4 = alo_u
        c4 = alo_p
        ixz4 = ixz
        iyz4 = iyz
        iu27 = zero
      endif
      !
      !update the downwind and current point for next step
      opac_d=opac_p
      dels_d=dels_u
      scont_d=scont_p
      x_p=x_u
      y_p=y_u
      z_p=z_u
      opac_p=opac_u
      scont_p=scont_u
      !   write(*,*) 't1', i, ixz1, ixz2, ixz3, iyz1, iyz2, iyz3
      !exit z-loop if z-layer is hit
      if(lboundz) exit
      !
    enddo
    !
    if(.not.lboundz) stop 'error in flc_cont3d: no previous z-layer has been found, increase nslab'
    !
    !for benchmark02
    !close(1)
    !
    !
    !
    !set the intensity + alo
    int3d(xindx,yindx,zindx) = int_u*abs_lc + contr_lc
    !
    alo_u1 = b1+a1*c2
    alo_p1 = c1
    alo_u2 = a1*b2+a1*a2*c3
    alo_u3 = a1*a2*b3 + a1*a2*a3*c4
    abs_scc = a1*a2*a3*a4
    !
    alocont_o_nn3d(xindxm1,yindxm1,zindxm1,q27) = alo_u1*c14_scontu1 + alo_u2*(c15_scontu2*ixz1+c17_scontu2*iyz1) + &
    alo_u3*c18_scontu3*(ixz1*iyz2+iyz1*ixz2) + abs_scc*iu27
    alocont_o_nn3d(xindx,yindxm1,zindxm1,q26) = alo_u1*c15_scontu1 + alo_u2*c18_scontu2*iyz1 + abs_scc*iu26
    alocont_o_nn3d(xindxm1,yindx,zindxm1,q24) = alo_u1*c17_scontu1 + alo_u2*c18_scontu2*ixz1 + abs_scc*iu24
    alocont_o_nn3d(xindx,yindx,zindxm1,q23) = alo_u1*c18_scontu1 + abs_scc*iu23
    !
    alocont_o_nn3d(xindxm1,yindxm1,zindx,q18) = alo_u1*c23_scontu1 + alo_u2*(iyz1*c26_scontu2+ixz1*c24_scontu2) + abs_scc*iu18
    alocont_o_nn3d(xindx,yindxm1,zindx,q17) = alo_u1*c24_scontu1 + abs_scc*iu17
    alocont_o_nn3d(xindxm1,yindx,zindx,q15) = alo_u1*c26_scontu1 + abs_scc*iu15
    alocont_o_nn3d(xindx,yindx,zindx,q14) = alo_p1 + abs_scc*iu14
    !

    !if(xindx.eq.2.and.yindx.eq.2.and.zindx.eq.16) then
    !   write(*,*) 't2', q18, scont3d(xindxm1,yindxm1,zindx), alocont_o_nn3d(xindxm1,yindxm1,zindx,q18), int3d(xindx,yindx,zindx), iu18, int_u, &
    !        alo_u1*c23_scontu1 + alo_u2*(iyz1*c26_scontu2+ixz1*c24_scontu2), contr_lc, iu18, int_u, alo_u1*c23_scontu1, alo_u2*(iyz1*c26_scontu2+ixz1*c24_scontu2), alo_u1, contr_sc
    !!   , int3d(nx-1,yindx-1,zindx+1), alocont_o_nn3d(nx,yindx,zindx,19), alocont_o_nn3d(xindx,yindx,zindx,q1)!
    !endif

    !
    !
  end subroutine flc_cont3d_lin


  !
end module mod_conttrans3d
