module mod_scont_new3d

  use prog_type
  use fund_const

  implicit none

contains

  subroutine scont_new3d(verbose)
    !
    !-----------------------------------------------------------------------
    !---------calculates new iterate of continuum source function-----------
    !   different options: classical lambda-iteration
    !                      diagonal of lambda-matrix
    !                      nearest neighbours
    !-----------------------------------------------------------------------
    !
    use options, only: opt_alo_cont
    !
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: err
    !
    if(present(verbose)) ver=verbose
    !
    if(ver) write(*,*) '--------------calculating new iterate for source function (alo)----------------'
    if(ver) write(*,*)
    !
    !-----------------------------------------------------------------------
    !
    select case(opt_alo_cont)
    case(0)
      call scont_new3d_classic(verbose=ver)
    case(1)
      call scont_new3d_diag(verbose=ver)
    case(2)
      call scont_new3d_dn(verbose=ver)
    case(3)
      call scont_new3d_nn(verbose=ver)
    case default
      stop 'set option opt_alo_cont'
    end select
    !
    !calculate new source function on ghost points
    !
  end subroutine scont_new3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine scont_new3d_classic(verbose)
    !
    !-----------------------------------------------------------------------
    !--------calculates new iterate of continuum source function------------
    !-------------------classical lambda-iteration--------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, imaskb3d, scont3d, bnue3d, mint3d
    use params_input, only: eps_cont
    !
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    !
    ! ... local functions
    !
    !-----------------------------------------------------------------------
    !
    if(present(verbose)) ver=verbose
    !
    do k=1, nz
      do j=1, ny
        do i=1, nx
          select case(imaskb3d(i,j,k))
          case(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
            scont3d(i,j,k) = (1.d0-eps_cont) * mint3d(i,j,k) + eps_cont*bnue3d(i,j,k)
          case default
          endselect
        enddo
      enddo
    enddo
    !
    !
  end subroutine scont_new3d_classic
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine scont_new3d_diag(verbose)
    !
    !-----------------------------------------------------------------------
    !--------calculates new iterate of continuum source function------------
    !---------approximate lambda iteration using only diagonal--------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, scont3d, alocont_nn3d, mint3d, imaskb3d, bnue3d
    use params_input, only: eps_cont
    !
    !
    !! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: indx_1d
    real(dp) :: dummy1, dummy2, scont_new
    real(dp) :: wor
    !
    ! ... local functions
    !
    if(present(verbose)) ver=verbose
    !
    !----------------calculating snew directly for 2-d arrays---------------
    !
    dummy2=1.d0-eps_cont
    !wor=0.99999d0
    !alocont_nn3d(:,:,14)=alocont_nn3d(:,:,14)*wor
    !
    do k=1, nz
      do j=1, ny
        do i=1, nx
          select case(imaskb3d(i,j,k))
          case(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
            dummy1=1.d0-(1.d0-eps_cont)*alocont_nn3d(i,j,k,14)
            scont_new=(dummy2/dummy1) * mint3d(i,j,k) - &
            (dummy2/dummy1) * alocont_nn3d(i,j,k,14) * scont3d(i,j,k) + &
            (eps_cont/dummy1)*bnue3d(i,j,k)
            if(scont_new.ge.0.d0) scont3d(i,j,k) = scont_new
          case default
          endselect
        enddo
      enddo
    enddo
    !

    !stop 'go on in scont_new3d_diag'
    !
  end subroutine scont_new3d_diag
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine scont_new3d_dn(verbose)
    !
    !
    !-----------------------------------------------------------------------
    !------------calculates new iterate of source function------------------
    !---approximate lambda iteration using nearest neighbour contribution---
    !----------alo has to be stored in sparse-matrix formats----------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, scont3d, mint3d, imaskb3d, bnue3d, &
    alocont_nn3d, alocont_rowindx, alocont_colindx, &
    alocont_data, alocont_data_diag
    use params_input, only: eps_cont
    use mod_math, only: conv_indx_1d_to_3d, conv_indx_3d_to_1d
    use mod_sparse, only: matmul_coo, jsor_coo
    !
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: indx_1d, indx_x, indx_y, indx_z, err, indx_rspec
    real(dp) :: dummy1, dummy2, rspec
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: dummy_vec, sum_vec, scont_vec, mint_vec, bnue_vec
    logical, dimension(:), allocatable :: mask_vec
    !
    ! ... local functions
    !
    if(present(verbose)) ver=verbose
    !
    !-----------------------------------------------------------------------
    !---------------------allocation of local arrays------------------------
    !-----------------------------------------------------------------------
    !
    allocate(bnue_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_dn'
    allocate(scont_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_dn'
    allocate(mint_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_dn'
    allocate(dummy_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_dn'
    allocate(mask_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_dn'
    allocate(sum_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_dn'
    !
    !----------------transform 3d-arrays to 1-d array-----------------------
    !
    mint_vec=zero
    scont_vec=zero
    bnue_vec=zero
    dummy_vec=zero
    mask_vec=.false.
    !
    do i=1, nx
      do j=1, ny
        do k=1, nz
          select case(imaskb3d(i,j,k))
          case(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d)
            mask_vec(indx_1d)=.true.
            mint_vec(indx_1d)=mint3d(i,j,k)
            scont_vec(indx_1d)=scont3d(i,j,k)
            bnue_vec(indx_1d)=bnue3d(i,j,k)
          case default
          endselect
        enddo
      enddo
    enddo
    !
    call calc_alocont_dn3d_coo
    !
    !------------set up linear system that needs to be solved---------------
    !-------------------in order to obtain s_new----------------------------
    !
    call matmul_coo(alocont_data, alocont_colindx, alocont_rowindx, scont_vec, dummy_vec, nx*ny*nz, 7*nx*ny*nz)
    !
    scont_vec=(1.d0-eps_cont)*dummy_vec
    !
    bnue_vec=eps_cont*bnue_vec
    !
    mint_vec=(1.d0-eps_cont)*mint_vec
    !
    dummy_vec=mint_vec-scont_vec+bnue_vec
    !
    alocont_data=(1.d0-eps_cont)*alocont_data
    alocont_data_diag=1.d0-(1.d0-eps_cont)*alocont_data_diag
    !
    do i=1, 7*nx*ny*nz
      if(alocont_colindx(i).eq.alocont_rowindx(i)) then
        alocont_data(i)=1.d0-alocont_data(i)
      else
        alocont_data(i)=-1.d0*alocont_data(i)
      endif
    enddo
    !
    !--------------check if matrix is diagonal dominant---------------------
    !-----------------and estimate spectral radius--------------------------
    !
    sum_vec=0.d0
    !
    do i=1, 7*nx*ny*nz
      if(alocont_rowindx(i).ne.alocont_colindx(i)) then
        sum_vec(alocont_rowindx(i)) = sum_vec(alocont_rowindx(i)) + abs(alocont_data(i))
      endif
    enddo
    !
    rspec=maxval(sum_vec/alocont_data_diag)
    indx_rspec=maxloc(sum_vec/alocont_data_diag, 1)
    call conv_indx_1d_to_3d(indx_rspec, nx, ny, indx_x, indx_y, indx_z)
    if(ver) write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
    if(ver) write(*,'(a30,i10,3i5)') 'at 1d, 3d indices', indx_rspec, indx_x, indx_y, indx_z
    if(ver) write(*,*)
    !
    !--------------------linear system now reads:---------------------------
    !              alocont_matrix * s_new = scont_vec
    !
    if(ver) write(*,*) '-----------------------inverting nearest neighbour alo-------------------------'
    !
    dummy_vec=-dummy_vec
    !
    !call output_alocont_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, nx*ny*nz, 7*nx*ny*nz)
    !
    call jsor_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, 1.d0, nx*ny*nz, 7*nx*ny*nz, .false., scont_vec, verbose=ver)
    !
    !---------back-transformation of source function on 3-d grid------------
    !
    do i=1, nx*ny*nz
      if(mask_vec(i)) then
        call conv_indx_1d_to_3d(i, nx, ny, indx_x, indx_y, indx_z)
        !      write(*,*) indx_x, indx_z, scont_vec(i)
        if(scont_vec(i).ge.zero) then
          scont3d(indx_x,indx_y,indx_z)=scont_vec(i)
        endif
      endif
    enddo
    !
    !update periodic boundary condition
    do k=1, nz
      !right boundary = left boundary
      do j=2,ny-1
        scont3d(nx,j,k)=scont3d(1,j,k)
      enddo
      !back boundary = front boundary
      do i=2, nx-1
        scont3d(i,ny,k)=scont3d(i,1,k)
      enddo
      !right front boundary = left front boundary
      scont3d(nx,1,k)=scont3d(1,1,k)
      !left back boundary = left front boundary
      scont3d(1,ny,k)=scont3d(1,1,k)
      !right back boundary = left front boundary
      scont3d(nx,ny,k)=scont3d(1,1,k)
    enddo
    !
    if(ver) write(*,*) '-------------------------------------------------------------------------------'
    if(ver) write(*,*)

    !write(*,*) scont3d(1,1)
    !stop 'go on in scont_new'
    !
    !-----------------------------------------------------------------------
    !
  end subroutine scont_new3d_dn
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine scont_new3d_nn(verbose)
    !
    !
    !-----------------------------------------------------------------------
    !------------calculates new iterate of source function------------------
    !---approximate lambda iteration using nearest neighbour contribution---
    !----------alo has to be stored in sparse-matrix formats----------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, scont3d, mint3d, imaskb3d, bnue3d, x, y, z, &
    alocont_nn3d, alocont_rowindx, alocont_colindx, alocont_data, alocont_data_diag
    use params_input, only: eps_cont
    use mod_math, only: conv_indx_1d_to_3d, conv_indx_3d_to_1d
    use mod_sparse, only: matmul_coo, jsor_coo    
    !
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k
    integer(i4b) :: indx_1d, indx_x, indx_y, indx_z, err, indx_rspec
    real(dp) :: dummy1, dummy2, rspec
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: mint_vec, scont_vec, bnue_vec, dummy_vec, sum_vec
    logical, dimension(:), allocatable :: mask_vec
    !
    ! ... local functions
    !
    if(present(verbose)) ver=verbose
    !
    !-----------------------------------------------------------------------
    !---------------------allocation of local arrays------------------------
    !-----------------------------------------------------------------------
    !
    allocate(bnue_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_nn'
    allocate(scont_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_nn'
    allocate(mint_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_nn'
    allocate(dummy_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_nn'
    allocate(mask_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_nn'
    allocate(sum_vec(nx*ny*nz), stat=err)
    if(err.ne.0) stop 'allocation error scont_new3d_nn'
    !
    !----------------transform 3d-arrays to 1-d array-----------------------
    !
    mint_vec=zero
    scont_vec=zero
    bnue_vec=zero
    dummy_vec=zero
    mask_vec=.false.
    !
    do i=1, nx
      do j=1, ny
        do k=1, nz
          select case(imaskb3d(i,j,k))
          case(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d)
            mask_vec(indx_1d)=.true.
            mint_vec(indx_1d)=mint3d(i,j,k)
            scont_vec(indx_1d)=scont3d(i,j,k)
            bnue_vec(indx_1d)=bnue3d(i,j,k)
          case default
          endselect
        enddo
      enddo
    enddo
    !
    call calc_alocont_nn3d_coo
    !
    !------------set up linear system that needs to be solved---------------
    !-------------------in order to obtain s_new----------------------------
    !
    call matmul_coo(alocont_data, alocont_colindx, alocont_rowindx, scont_vec, dummy_vec, nx*ny*nz, 27*nx*ny*nz)
    !
    scont_vec=(one-eps_cont)*dummy_vec
    !
    bnue_vec=eps_cont*bnue_vec
    !
    mint_vec=(one-eps_cont)*mint_vec
    !
    dummy_vec=mint_vec-scont_vec+bnue_vec
    !
    alocont_data=(one-eps_cont)*alocont_data
    alocont_data_diag=one-(one-eps_cont)*alocont_data_diag
    !
    do i=1, 27*nx*ny*nz
      if(alocont_colindx(i).eq.alocont_rowindx(i)) then
        alocont_data(i)=one-alocont_data(i)
      else
        alocont_data(i)=-one*alocont_data(i)
      endif
    enddo
    !
    !--------------check if matrix is diagonal dominant---------------------
    !-----------------and estimate spectral radius--------------------------
    !
    sum_vec=0.d0
    !
    do i=1, 27*nx*ny*nz
      if(alocont_rowindx(i).ne.alocont_colindx(i)) then
        sum_vec(alocont_rowindx(i)) = sum_vec(alocont_rowindx(i)) + abs(alocont_data(i))
      endif
    enddo
    !
    rspec=maxval(sum_vec/alocont_data_diag)
    indx_rspec=maxloc(sum_vec/alocont_data_diag, 1)
    call conv_indx_1d_to_3d(indx_rspec, nx, ny, indx_x, indx_y, indx_z)
    if(ver) write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
    if(ver) write(*,'(a30,i10,3i5)') 'at 1d, 3d indices', indx_rspec, indx_x, indx_y, indx_z
    if(ver) write(*,*)
    !
    !--------------------linear system now reads:---------------------------
    !              alocont_matrix * s_new = scont_vec
    !
    if(ver) write(*,*) '-----------------------inverting nearest neighbour alo-------------------------'
    !
    dummy_vec=-dummy_vec
    !
    !call output_alocont_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, nx*ny*nz, 27*nx*ny*nz)
    !
    call jsor_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, one, nx*ny*nz, 27*nx*ny*nz, .false., scont_vec, verbose=ver)
    !
    !
    !---------back-transformation of source function on 3-d grid------------
    !
    do i=1, nx*ny*nz
      if(mask_vec(i)) then
        call conv_indx_1d_to_3d(i, nx, ny, indx_x, indx_y, indx_z)
        if(scont_vec(i).ge.zero) then
          scont3d(indx_x,indx_y,indx_z)=scont_vec(i)
        endif
        !    write(*,*) scont3d(indx_x,indx_z), imaskb3d(indx_x,indx_z)
      endif
    enddo
    !
    !update periodic boundary condition
    do k=1, nz
      !right boundary = left boundary
      do j=2,ny-1
        scont3d(nx,j,k)=scont3d(1,j,k)
      enddo
      !back boundary = front boundary
      do i=2, nx-1
        scont3d(i,ny,k)=scont3d(i,1,k)
      enddo
      !right front boundary = left front boundary
      scont3d(nx,1,k)=scont3d(1,1,k)
      !left back boundary = left front boundary
      scont3d(1,ny,k)=scont3d(1,1,k)
      !right back boundary = left front boundary
      scont3d(nx,ny,k)=scont3d(1,1,k)
    enddo
    !
    if(ver) write(*,*) '-------------------------------------------------------------------------------'
    if(ver) write(*,*)
    !
    !-----------------------------------------------------------------------
    !
  end subroutine scont_new3d_nn
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calc_alocont_dn3d_coo(verbose)
    !
    !-----------------------------------------------------------------------
    !------------calculates approximate lambda operator---------------------
    !-------using only direct neighbours (6 neighbours + local point)-------
    !----------------------coo storage format-------------------------------
    !-----------------------------------------------------------------------
    !
    !         IMPORTANT NOTE: FOR PERIODIC BOUNDARY CONDITIONS
    !              ALO IS ONLY CALCULATED ON LEFT AND FRONT BOUNDARY;
    !              RIGHT AND BACK BOUNDARY SOURCE FUNCTION IS UPDATED AT THE END
    !              OTHERWISE: INVERSION BECOMES ILL-CONDITIONED (SPECTRAL RADIUS > 1)
    !
    !
    use mod_grid3d, only: nx, ny, nz, alocont_nn3d, imaskb3d, imask3d, &
    alocont_data, alocont_data_diag, alocont_colindx, alocont_rowindx
    use mod_math, only: conv_indx_3d_to_1d
    !
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: iindx, indx_1d_row, indx_1d_col, indx_x, indx_z
    integer :: ndiags_ijk, ndiags_max, ndiags_min
    real(dp) :: rspec, rspec_max
    !
    ! ... local arrays
    !
    ! ... local functions
    !
    ! ... local characters
    character(len=50) :: enter
    !
    if(present(verbose)) ver=verbose
    !
    !-------------------allocation of alo-matrix----------------------------
    !
    if(.not.allocated(alocont_data)) then
      allocate(alocont_data(7*nx*ny*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn3d_coo: alocont_data'
    endif

    if(.not.allocated(alocont_data_diag)) then
      allocate(alocont_data_diag(nx*ny*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn3d_coo: alocont_data_diag'
    endif
    !
    if(.not.allocated(alocont_colindx)) then
      allocate(alocont_colindx(7*nx*ny*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn3d_coo: alocont_col_indx'
    endif
    !
    if(.not.allocated(alocont_rowindx)) then
      allocate(alocont_rowindx(7*nx*ny*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn3d_coo: alocont_row_indx'
    endif
    !
    !-----------------calculate alo in coo storage format-------------------
    !
    !define maximum allowed spectral radius
    rspec_max=one
    !
    alocont_data=zero
    alocont_data_diag=zero
    alocont_rowindx=1
    alocont_colindx=1
    !
    iindx=1
    ndiags_max = 0
    ndiags_min = 7
    !
    !stop 'go on in alocont_dn3d'

    do k=1, nz
      do j=1, ny
        do i=1, nx
          !
          select case(imaskb3d(i,j,k))
            !
          case(3,4,6,8)
            !
            !-------------------------standard points-------------------------------
            !-----------------points adjacent to left boundary----------------------
            !----------------points adjacent to front boundary----------------------
            !--------------points adjacent to left front boundary ------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
          case(5,9)
            !
            !-----------------points adjacent to right boundary---------------------
            !-----------------points adjacent to right front boundary---------------
            !---------------------(influence left boundary)-------------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
          case(7,10)
            !
            !-----------------points adjacent to back boundary----------------------
            !-----------------points adjacent to left back boundary-----------------
            !---------------------(influence front boundary)------------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
          case(11)
            !
            !-----------------points adjacent to right back boundary----------------
            !--------------------(influence left front boundary)--------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
          case(12,16)
            !
            !---------------------left boundary (xz-plane)--------------------------
            !--------------left boundary adjacent to front boundary-----------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !!include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction (periodic boundary condition)
            !               stop 'this element causes osciallations'
            call conv_indx_3d_to_1d(nx-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !               write(*,*) i, j, k, alocont_nn3d(i,j,k,13), alocont_nn3d(nx,j,k,13)
            !
            !!include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !!
          case(22)
            !
            !----------------left boundary adjecent to back boundary----------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !!
            !!include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
          case(14,17)
            !
            !---------------------front boundary (yz-plane)-------------------------
            !-----------------front boundary adjacent to left boundary--------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
          case(20)
            !
            !-----------------front boundary adjacent to right boundary-------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
          case(18)
            !
            !--------------------left front boundary (edge)-------------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(nx-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !!include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
            !---------------------right and back boundaries-------------------------
            !-----------------(only diagonal part, will be updated------------------
            !------------according to periodic boundary condition anyways)----------
            !
          case(13,15,19,21,23,24,25,26,27)
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
          case(1,41,51,61,71,81,91,101,111,121,141,161,171,181,201)!,211,221,191,131,251,271,261,151,231,241)
            !
            !--------------------bottom boundary (non-periodic)--------------------
            !-----------------------(only in z-direction)---------------------------
            !
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
          case(2,42,52,62,72,82,92,102,112,122,142,162,172,182,202)!,212,222,192,132,252,272,262,152,232,242)
            !
            !---------------------top boundary (non-periodic)-----------------------
            !-----------------------(only in z-direction)---------------------------
            !
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
          case default
          end select
          !
          iindx=iindx+7
          !
        enddo
      enddo
    enddo
    !
    !
    !write(*,*)
    !write(*,*) 'maximum number of used neighbours for ALO calculation', ndiags_max-1
    !write(*,*) 'minimum number of used neighbours for ALO calculation', ndiags_min-1
    !write(*,*)

    !open(1, file='TRASH/alo_sc3d.dat')
    !do i=1, 5*nx*nz
    !   write(1,'(2i10,es20.8)'), alocont_rowindx(i), alocont_colindx(i), alocont_data(i)
    !   write(*,'(2i10,es20.8)'), alocont_rowindx(i), alocont_colindx(i), alocont_data(i)
    !enddo
    !
    !stop 'go on in calc_alocont_dn3d_coo'

  end subroutine calc_alocont_dn3d_coo

  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calc_alocont_nn3d_coo(verbose)
    !
    !-----------------------------------------------------------------------
    !------------calculates approximate lambda operator---------------------
    !-------using only nearest neighbours (26 neighbours + local point)-----
    !----------------------coo storage format-------------------------------
    !-----------------------------------------------------------------------
    !
    use mod_grid3d, only: nx, ny, nz, alocont_nn3d, &
    alocont_data, alocont_data_diag, alocont_colindx, alocont_rowindx, &
    imaskb3d, imask3d, x, z
    use mod_math, only: conv_indx_3d_to_1d
    !
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local logicals
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, err
    integer(i4b) :: iindx, indx_1d_col, indx_1d_row, indx_x, indx_y, indx_z
    integer :: ndiags_ijk, ndiags_max, ndiags_min
    real(dp) :: sum, eps, rspec, rspec_max
    !
    ! ... local arrays
    !
    ! ... local functions
    !
    ! ... local characters
    character(len=50) :: enter
    !
    ! ... for debug
    real(dp), dimension(:), allocatable :: alocont_data2, alocont_data_diag2
    integer(i4b), dimension(:), allocatable :: alocont_rowindx2, alocont_colindx2
    !
    if(present(verbose)) ver=verbose
    !
    !-------------------allocation of alo-matrix----------------------------
    !
    if(.not.allocated(alocont_data)) then
      allocate(alocont_data(27*nx*ny*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn3d_coo: alocont_data'
    endif

    if(.not.allocated(alocont_data_diag)) then
      allocate(alocont_data_diag(nx*ny*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn3d_coo: alocont_data_diag'
    endif
    !
    if(.not.allocated(alocont_colindx)) then
      allocate(alocont_colindx(27*nx*ny*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn2_coo: alocont_col_indx'
    endif
    !
    if(.not.allocated(alocont_rowindx)) then
      allocate(alocont_rowindx(27*nx*ny*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn3d_coo: alocont_row_indx'
    endif
    !
    !-----------------calculate alo in coo storage format-------------------
    !
    !define maximum allowed spectral radius
    rspec_max=one
    !
    alocont_data=zero
    alocont_data_diag=zero
    alocont_rowindx=1
    alocont_colindx=1
    !
    iindx=1
    ndiags_max = 0
    ndiags_min = 27


    do k=1, nz
      do j=1, ny
        do i=1, nx
          !
          select case(imaskb3d(i,j,k))
            !
          case(3,4,6,8)
            !
            !-------------------------standard points-------------------------------
            !-----------------points adjacent to left boundary----------------------
            !----------------points adjacent to front boundary----------------------
            !--------------points adjacent to left front boundary ------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
            !include neighbour in positive x- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn3d(i,j,k,24)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
            !
            !include neighbour in positive x- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn3d(i,j,k,6)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
            !
            !include neighbour in negative x- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+9) = alocont_nn3d(i,j,k,22)
            alocont_colindx(iindx+9)=indx_1d_col
            alocont_rowindx(iindx+9)=indx_1d_row
            !
            !include neighbour in negative x- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+10) = alocont_nn3d(i,j,k,4)
            alocont_colindx(iindx+10)=indx_1d_col
            alocont_rowindx(iindx+10)=indx_1d_row
            !
            !include neighbour in positive y- and positive z-direction
            call conv_indx_3d_to_1d(i,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+11) = alocont_nn3d(i,j,k,26)
            alocont_colindx(iindx+11)=indx_1d_col
            alocont_rowindx(iindx+11)=indx_1d_row
            !
            !include neighbour in positive y- and negative z-direction
            call conv_indx_3d_to_1d(i,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+12) = alocont_nn3d(i,j,k,8)
            alocont_colindx(iindx+12)=indx_1d_col
            alocont_rowindx(iindx+12)=indx_1d_row
            !
            !include neighbour in negative y- and positive z-direction
            call conv_indx_3d_to_1d(i,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+13) = alocont_nn3d(i,j,k,20)
            alocont_colindx(iindx+13)=indx_1d_col
            alocont_rowindx(iindx+13)=indx_1d_row
            !
            !include neighbour in negative y- and negative z-direction
            call conv_indx_3d_to_1d(i,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+14) = alocont_nn3d(i,j,k,2)
            alocont_colindx(iindx+14)=indx_1d_col
            alocont_rowindx(iindx+14)=indx_1d_row
            !
            !include neighbour in positive x- and positive y-direction
            call conv_indx_3d_to_1d(i+1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+15) = alocont_nn3d(i,j,k,18)
            alocont_colindx(iindx+15)=indx_1d_col
            alocont_rowindx(iindx+15)=indx_1d_row
            !
            !include neighbour in positive x- and negative y-direction
            call conv_indx_3d_to_1d(i+1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+16) = alocont_nn3d(i,j,k,12)
            alocont_colindx(iindx+16)=indx_1d_col
            alocont_rowindx(iindx+16)=indx_1d_row
            !
            !include neighbour in negative x- and positive y-direction
            call conv_indx_3d_to_1d(i-1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+17) = alocont_nn3d(i,j,k,16)
            alocont_colindx(iindx+17)=indx_1d_col
            alocont_rowindx(iindx+17)=indx_1d_row
            !
            !include neighbour in negative x- and negative y-direction
            call conv_indx_3d_to_1d(i-1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+18) = alocont_nn3d(i,j,k,10)
            alocont_colindx(iindx+18)=indx_1d_col
            alocont_rowindx(iindx+18)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+19) = alocont_nn3d(i,j,k,27)
            alocont_colindx(iindx+19)=indx_1d_col
            alocont_rowindx(iindx+19)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+20) = alocont_nn3d(i,j,k,25)
            alocont_colindx(iindx+20)=indx_1d_col
            alocont_rowindx(iindx+20)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+21) = alocont_nn3d(i,j,k,21)
            alocont_colindx(iindx+21)=indx_1d_col
            alocont_rowindx(iindx+21)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+22) = alocont_nn3d(i,j,k,19)
            alocont_colindx(iindx+22)=indx_1d_col
            alocont_rowindx(iindx+22)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+23) = alocont_nn3d(i,j,k,9)
            alocont_colindx(iindx+23)=indx_1d_col
            alocont_rowindx(iindx+23)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+24) = alocont_nn3d(i,j,k,7)
            alocont_colindx(iindx+24)=indx_1d_col
            alocont_rowindx(iindx+24)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+25) = alocont_nn3d(i,j,k,3)
            alocont_colindx(iindx+25)=indx_1d_col
            alocont_rowindx(iindx+25)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+26) = alocont_nn3d(i,j,k,1)
            alocont_colindx(iindx+26)=indx_1d_col
            alocont_rowindx(iindx+26)=indx_1d_row
            !
          case(5,9)
            !
            !-----------------points adjacent to right boundary---------------------
            !-----------------points adjacent to right front boundary---------------
            !---------------------(influence left boundary)-------------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
            !include neighbour in positive x- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn3d(i,j,k,24)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
            !
            !include neighbour in positive x- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn3d(i,j,k,6)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
            !
            !include neighbour in negative x- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+9) = alocont_nn3d(i,j,k,22)
            alocont_colindx(iindx+9)=indx_1d_col
            alocont_rowindx(iindx+9)=indx_1d_row
            !
            !include neighbour in negative x- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+10) = alocont_nn3d(i,j,k,4)
            alocont_colindx(iindx+10)=indx_1d_col
            alocont_rowindx(iindx+10)=indx_1d_row
            !
            !include neighbour in positive y- and positive z-direction
            call conv_indx_3d_to_1d(i,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+11) = alocont_nn3d(i,j,k,26)
            alocont_colindx(iindx+11)=indx_1d_col
            alocont_rowindx(iindx+11)=indx_1d_row
            !
            !include neighbour in positive y- and negative z-direction
            call conv_indx_3d_to_1d(i,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+12) = alocont_nn3d(i,j,k,8)
            alocont_colindx(iindx+12)=indx_1d_col
            alocont_rowindx(iindx+12)=indx_1d_row
            !
            !include neighbour in negative y- and positive z-direction
            call conv_indx_3d_to_1d(i,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+13) = alocont_nn3d(i,j,k,20)
            alocont_colindx(iindx+13)=indx_1d_col
            alocont_rowindx(iindx+13)=indx_1d_row
            !
            !include neighbour in negative y- and negative z-direction
            call conv_indx_3d_to_1d(i,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+14) = alocont_nn3d(i,j,k,2)
            alocont_colindx(iindx+14)=indx_1d_col
            alocont_rowindx(iindx+14)=indx_1d_row
            !
            !include neighbour in positive x- and positive y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+15) = alocont_nn3d(i,j,k,18)
            alocont_colindx(iindx+15)=indx_1d_col
            alocont_rowindx(iindx+15)=indx_1d_row
            !
            !include neighbour in positive x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+16) = alocont_nn3d(i,j,k,12)
            alocont_colindx(iindx+16)=indx_1d_col
            alocont_rowindx(iindx+16)=indx_1d_row
            !
            !include neighbour in negative x- and positive y-direction
            call conv_indx_3d_to_1d(i-1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+17) = alocont_nn3d(i,j,k,16)
            alocont_colindx(iindx+17)=indx_1d_col
            alocont_rowindx(iindx+17)=indx_1d_row
            !
            !include neighbour in negative x- and negative y-direction
            call conv_indx_3d_to_1d(i-1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+18) = alocont_nn3d(i,j,k,10)
            alocont_colindx(iindx+18)=indx_1d_col
            alocont_rowindx(iindx+18)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+19) = alocont_nn3d(i,j,k,27)
            alocont_colindx(iindx+19)=indx_1d_col
            alocont_rowindx(iindx+19)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+20) = alocont_nn3d(i,j,k,25)
            alocont_colindx(iindx+20)=indx_1d_col
            alocont_rowindx(iindx+20)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+21) = alocont_nn3d(i,j,k,21)
            alocont_colindx(iindx+21)=indx_1d_col
            alocont_rowindx(iindx+21)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+22) = alocont_nn3d(i,j,k,19)
            alocont_colindx(iindx+22)=indx_1d_col
            alocont_rowindx(iindx+22)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+23) = alocont_nn3d(i,j,k,9)
            alocont_colindx(iindx+23)=indx_1d_col
            alocont_rowindx(iindx+23)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+24) = alocont_nn3d(i,j,k,7)
            alocont_colindx(iindx+24)=indx_1d_col
            alocont_rowindx(iindx+24)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+25) = alocont_nn3d(i,j,k,3)
            alocont_colindx(iindx+25)=indx_1d_col
            alocont_rowindx(iindx+25)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+26) = alocont_nn3d(i,j,k,1)
            alocont_colindx(iindx+26)=indx_1d_col
            alocont_rowindx(iindx+26)=indx_1d_row
            !
          case(7,10)
            !
            !-----------------points adjacent to back boundary----------------------
            !-----------------points adjacent to left back boundary-----------------
            !---------------------(influence front boundary)------------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
            !include neighbour in positive x- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn3d(i,j,k,24)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
            !
            !include neighbour in positive x- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn3d(i,j,k,6)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
            !
            !include neighbour in negative x- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+9) = alocont_nn3d(i,j,k,22)
            alocont_colindx(iindx+9)=indx_1d_col
            alocont_rowindx(iindx+9)=indx_1d_row
            !
            !include neighbour in negative x- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+10) = alocont_nn3d(i,j,k,4)
            alocont_colindx(iindx+10)=indx_1d_col
            alocont_rowindx(iindx+10)=indx_1d_row
            !
            !include neighbour in positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+11) = alocont_nn3d(i,j,k,26)
            alocont_colindx(iindx+11)=indx_1d_col
            alocont_rowindx(iindx+11)=indx_1d_row
            !
            !include neighbour in positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+12) = alocont_nn3d(i,j,k,8)
            alocont_colindx(iindx+12)=indx_1d_col
            alocont_rowindx(iindx+12)=indx_1d_row
            !
            !include neighbour in negative y- and positive z-direction
            call conv_indx_3d_to_1d(i,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+13) = alocont_nn3d(i,j,k,20)
            alocont_colindx(iindx+13)=indx_1d_col
            alocont_rowindx(iindx+13)=indx_1d_row
            !
            !include neighbour in negative y- and negative z-direction
            call conv_indx_3d_to_1d(i,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+14) = alocont_nn3d(i,j,k,2)
            alocont_colindx(iindx+14)=indx_1d_col
            alocont_rowindx(iindx+14)=indx_1d_row
            !
            !include neighbour in positive x- and positive y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+15) = alocont_nn3d(i,j,k,18)
            alocont_colindx(iindx+15)=indx_1d_col
            alocont_rowindx(iindx+15)=indx_1d_row
            !
            !include neighbour in positive x- and negative y-direction
            call conv_indx_3d_to_1d(i+1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+16) = alocont_nn3d(i,j,k,12)
            alocont_colindx(iindx+16)=indx_1d_col
            alocont_rowindx(iindx+16)=indx_1d_row
            !
            !include neighbour in negative x- and positive y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+17) = alocont_nn3d(i,j,k,16)
            alocont_colindx(iindx+17)=indx_1d_col
            alocont_rowindx(iindx+17)=indx_1d_row
            !
            !include neighbour in negative x- and negative y-direction
            call conv_indx_3d_to_1d(i-1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+18) = alocont_nn3d(i,j,k,10)
            alocont_colindx(iindx+18)=indx_1d_col
            alocont_rowindx(iindx+18)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+19) = alocont_nn3d(i,j,k,27)
            alocont_colindx(iindx+19)=indx_1d_col
            alocont_rowindx(iindx+19)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+20) = alocont_nn3d(i,j,k,25)
            alocont_colindx(iindx+20)=indx_1d_col
            alocont_rowindx(iindx+20)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+21) = alocont_nn3d(i,j,k,21)
            alocont_colindx(iindx+21)=indx_1d_col
            alocont_rowindx(iindx+21)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+22) = alocont_nn3d(i,j,k,19)
            alocont_colindx(iindx+22)=indx_1d_col
            alocont_rowindx(iindx+22)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+23) = alocont_nn3d(i,j,k,9)
            alocont_colindx(iindx+23)=indx_1d_col
            alocont_rowindx(iindx+23)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+24) = alocont_nn3d(i,j,k,7)
            alocont_colindx(iindx+24)=indx_1d_col
            alocont_rowindx(iindx+24)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+25) = alocont_nn3d(i,j,k,3)
            alocont_colindx(iindx+25)=indx_1d_col
            alocont_rowindx(iindx+25)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+26) = alocont_nn3d(i,j,k,1)
            alocont_colindx(iindx+26)=indx_1d_col
            alocont_rowindx(iindx+26)=indx_1d_row
            !
          case(11)
            !
            !-----------------points adjacent to right back boundary----------------
            !--------------------(influence left front boundary)--------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
            !include neighbour in positive x- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn3d(i,j,k,24)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
            !
            !include neighbour in positive x- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn3d(i,j,k,6)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
            !
            !include neighbour in negative x- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+9) = alocont_nn3d(i,j,k,22)
            alocont_colindx(iindx+9)=indx_1d_col
            alocont_rowindx(iindx+9)=indx_1d_row
            !
            !include neighbour in negative x- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+10) = alocont_nn3d(i,j,k,4)
            alocont_colindx(iindx+10)=indx_1d_col
            alocont_rowindx(iindx+10)=indx_1d_row
            !
            !include neighbour in positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+11) = alocont_nn3d(i,j,k,26)
            alocont_colindx(iindx+11)=indx_1d_col
            alocont_rowindx(iindx+11)=indx_1d_row
            !
            !include neighbour in positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+12) = alocont_nn3d(i,j,k,8)
            alocont_colindx(iindx+12)=indx_1d_col
            alocont_rowindx(iindx+12)=indx_1d_row
            !
            !include neighbour in negative y- and positive z-direction
            call conv_indx_3d_to_1d(i,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+13) = alocont_nn3d(i,j,k,20)
            alocont_colindx(iindx+13)=indx_1d_col
            alocont_rowindx(iindx+13)=indx_1d_row
            !
            !include neighbour in negative y- and negative z-direction
            call conv_indx_3d_to_1d(i,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+14) = alocont_nn3d(i,j,k,2)
            alocont_colindx(iindx+14)=indx_1d_col
            alocont_rowindx(iindx+14)=indx_1d_row
            !
            !include neighbour in positive x- and positive y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+15) = alocont_nn3d(i,j,k,18)
            alocont_colindx(iindx+15)=indx_1d_col
            alocont_rowindx(iindx+15)=indx_1d_row
            !
            !include neighbour in positive x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+16) = alocont_nn3d(i,j,k,12)
            alocont_colindx(iindx+16)=indx_1d_col
            alocont_rowindx(iindx+16)=indx_1d_row
            !
            !include neighbour in negative x- and positive y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+17) = alocont_nn3d(i,j,k,16)
            alocont_colindx(iindx+17)=indx_1d_col
            alocont_rowindx(iindx+17)=indx_1d_row
            !
            !include neighbour in negative x- and negative y-direction
            call conv_indx_3d_to_1d(i-1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+18) = alocont_nn3d(i,j,k,10)
            alocont_colindx(iindx+18)=indx_1d_col
            alocont_rowindx(iindx+18)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+19) = alocont_nn3d(i,j,k,27)
            alocont_colindx(iindx+19)=indx_1d_col
            alocont_rowindx(iindx+19)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+20) = alocont_nn3d(i,j,k,25)
            alocont_colindx(iindx+20)=indx_1d_col
            alocont_rowindx(iindx+20)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+21) = alocont_nn3d(i,j,k,21)
            alocont_colindx(iindx+21)=indx_1d_col
            alocont_rowindx(iindx+21)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+22) = alocont_nn3d(i,j,k,19)
            alocont_colindx(iindx+22)=indx_1d_col
            alocont_rowindx(iindx+22)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+23) = alocont_nn3d(i,j,k,9)
            alocont_colindx(iindx+23)=indx_1d_col
            alocont_rowindx(iindx+23)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+24) = alocont_nn3d(i,j,k,7)
            alocont_colindx(iindx+24)=indx_1d_col
            alocont_rowindx(iindx+24)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+25) = alocont_nn3d(i,j,k,3)
            alocont_colindx(iindx+25)=indx_1d_col
            alocont_rowindx(iindx+25)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and negative z-directiony
            call conv_indx_3d_to_1d(i-1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+26) = alocont_nn3d(i,j,k,1)
            alocont_colindx(iindx+26)=indx_1d_col
            alocont_rowindx(iindx+26)=indx_1d_row
            !
          case(12,16)
            !
            !---------------------left boundary (xz-plane)--------------------------
            !--------------left boundary adjacent to front boundary-----------------
            !-----------------(influence right adjacent points)
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !!include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !!include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
            !include neighbour in positive x- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn3d(i,j,k,24)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
            !
            !include neighbour in positive x- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn3d(i,j,k,6)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
            !
            !include neighbour in negative x- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+9) = alocont_nn3d(i,j,k,22)
            alocont_colindx(iindx+9)=indx_1d_col
            alocont_rowindx(iindx+9)=indx_1d_row
            !
            !include neighbour in negative x- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+10) = alocont_nn3d(i,j,k,4)
            alocont_colindx(iindx+10)=indx_1d_col
            alocont_rowindx(iindx+10)=indx_1d_row
            !
            !include neighbour in positive y- and positive z-direction
            call conv_indx_3d_to_1d(i,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+11) = alocont_nn3d(i,j,k,26)
            alocont_colindx(iindx+11)=indx_1d_col
            alocont_rowindx(iindx+11)=indx_1d_row
            !
            !include neighbour in positive y- and negative z-direction
            call conv_indx_3d_to_1d(i,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+12) = alocont_nn3d(i,j,k,8)
            alocont_colindx(iindx+12)=indx_1d_col
            alocont_rowindx(iindx+12)=indx_1d_row
            !
            !include neighbour in negative y- and positive z-direction
            call conv_indx_3d_to_1d(i,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+13) = alocont_nn3d(i,j,k,20)
            alocont_colindx(iindx+13)=indx_1d_col
            alocont_rowindx(iindx+13)=indx_1d_row
            !
            !include neighbour in negative y- and negative z-direction
            call conv_indx_3d_to_1d(i,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+14) = alocont_nn3d(i,j,k,2)
            alocont_colindx(iindx+14)=indx_1d_col
            alocont_rowindx(iindx+14)=indx_1d_row
            !
            !include neighbour in positive x- and positive y-direction
            call conv_indx_3d_to_1d(i+1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+15) = alocont_nn3d(i,j,k,18)
            alocont_colindx(iindx+15)=indx_1d_col
            alocont_rowindx(iindx+15)=indx_1d_row
            !
            !include neighbour in positive x- and negative y-direction
            call conv_indx_3d_to_1d(i+1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+16) = alocont_nn3d(i,j,k,12)
            alocont_colindx(iindx+16)=indx_1d_col
            alocont_rowindx(iindx+16)=indx_1d_row
            !
            !include neighbour in negative x- and positive y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+17) = alocont_nn3d(i,j,k,16)
            alocont_colindx(iindx+17)=indx_1d_col
            alocont_rowindx(iindx+17)=indx_1d_row
            !
            !include neighbour in negative x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+18) = alocont_nn3d(i,j,k,10)
            alocont_colindx(iindx+18)=indx_1d_col
            alocont_rowindx(iindx+18)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+19) = alocont_nn3d(i,j,k,27)
            alocont_colindx(iindx+19)=indx_1d_col
            alocont_rowindx(iindx+19)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+20) = alocont_nn3d(i,j,k,25)
            alocont_colindx(iindx+20)=indx_1d_col
            alocont_rowindx(iindx+20)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+21) = alocont_nn3d(i,j,k,21)
            alocont_colindx(iindx+21)=indx_1d_col
            alocont_rowindx(iindx+21)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+22) = alocont_nn3d(i,j,k,19)
            alocont_colindx(iindx+22)=indx_1d_col
            alocont_rowindx(iindx+22)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+23) = alocont_nn3d(i,j,k,9)
            alocont_colindx(iindx+23)=indx_1d_col
            alocont_rowindx(iindx+23)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+24) = alocont_nn3d(i,j,k,7)
            alocont_colindx(iindx+24)=indx_1d_col
            alocont_rowindx(iindx+24)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+25) = alocont_nn3d(i,j,k,3)
            alocont_colindx(iindx+25)=indx_1d_col
            alocont_rowindx(iindx+25)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+26) = alocont_nn3d(i,j,k,1)
            alocont_colindx(iindx+26)=indx_1d_col
            alocont_rowindx(iindx+26)=indx_1d_row


            !
          case(22)
            !
            !----------------left boundary adjecent to back boundary----------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction
            call conv_indx_3d_to_1d(i,j-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !!
            !!include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
          case(14,17)
            !
            !---------------------front boundary (yz-plane)-------------------------
            !-----------------front boundary adjacent to left boundary--------------
            !-------------------(influence back adjacent points)--------------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
            !include neighbour in positive x- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn3d(i,j,k,24)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
            !
            !include neighbour in positive x- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn3d(i,j,k,6)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
            !
            !include neighbour in negative x- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+9) = alocont_nn3d(i,j,k,22)
            alocont_colindx(iindx+9)=indx_1d_col
            alocont_rowindx(iindx+9)=indx_1d_row
            !
            !include neighbour in negative x- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+10) = alocont_nn3d(i,j,k,4)
            alocont_colindx(iindx+10)=indx_1d_col
            alocont_rowindx(iindx+10)=indx_1d_row
            !
            !include neighbour in positive y- and positive z-direction
            call conv_indx_3d_to_1d(i,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+11) = alocont_nn3d(i,j,k,26)
            alocont_colindx(iindx+11)=indx_1d_col
            alocont_rowindx(iindx+11)=indx_1d_row
            !
            !include neighbour in positive y- and negative z-direction
            call conv_indx_3d_to_1d(i,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+12) = alocont_nn3d(i,j,k,8)
            alocont_colindx(iindx+12)=indx_1d_col
            alocont_rowindx(iindx+12)=indx_1d_row
            !
            !include neighbour in negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+13) = alocont_nn3d(i,j,k,20)
            alocont_colindx(iindx+13)=indx_1d_col
            alocont_rowindx(iindx+13)=indx_1d_row
            !
            !include neighbour in negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+14) = alocont_nn3d(i,j,k,2)
            alocont_colindx(iindx+14)=indx_1d_col
            alocont_rowindx(iindx+14)=indx_1d_row
            !
            !include neighbour in positive x- and positive y-direction
            call conv_indx_3d_to_1d(i+1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+15) = alocont_nn3d(i,j,k,18)
            alocont_colindx(iindx+15)=indx_1d_col
            alocont_rowindx(iindx+15)=indx_1d_row
            !
            !include neighbour in positive x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+16) = alocont_nn3d(i,j,k,12)
            alocont_colindx(iindx+16)=indx_1d_col
            alocont_rowindx(iindx+16)=indx_1d_row
            !
            !include neighbour in negative x- and positive y-direction
            call conv_indx_3d_to_1d(i-1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+17) = alocont_nn3d(i,j,k,16)
            alocont_colindx(iindx+17)=indx_1d_col
            alocont_rowindx(iindx+17)=indx_1d_row
            !
            !include neighbour in negative x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+18) = alocont_nn3d(i,j,k,10)
            alocont_colindx(iindx+18)=indx_1d_col
            alocont_rowindx(iindx+18)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+19) = alocont_nn3d(i,j,k,27)
            alocont_colindx(iindx+19)=indx_1d_col
            alocont_rowindx(iindx+19)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+20) = alocont_nn3d(i,j,k,25)
            alocont_colindx(iindx+20)=indx_1d_col
            alocont_rowindx(iindx+20)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+21) = alocont_nn3d(i,j,k,21)
            alocont_colindx(iindx+21)=indx_1d_col
            alocont_rowindx(iindx+21)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+22) = alocont_nn3d(i,j,k,19)
            alocont_colindx(iindx+22)=indx_1d_col
            alocont_rowindx(iindx+22)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+23) = alocont_nn3d(i,j,k,9)
            alocont_colindx(iindx+23)=indx_1d_col
            alocont_rowindx(iindx+23)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+24) = alocont_nn3d(i,j,k,7)
            alocont_colindx(iindx+24)=indx_1d_col
            alocont_rowindx(iindx+24)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+25) = alocont_nn3d(i,j,k,3)
            alocont_colindx(iindx+25)=indx_1d_col
            alocont_rowindx(iindx+25)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+26) = alocont_nn3d(i,j,k,1)
            alocont_colindx(iindx+26)=indx_1d_col
            alocont_rowindx(iindx+26)=indx_1d_row
            !
            !
          case(20)
            !
            !-----------------front boundary adjacent to right boundary-------------
            !---------------(influences left and back adjacent boundary)------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction
            call conv_indx_3d_to_1d(i-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !include direct neighbour in negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row



            !
            !include neighbour in positive x- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn3d(i,j,k,24)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
            !
            !include neighbour in positive x- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn3d(i,j,k,6)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
            !
            !include neighbour in negative x- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+9) = alocont_nn3d(i,j,k,22)
            alocont_colindx(iindx+9)=indx_1d_col
            alocont_rowindx(iindx+9)=indx_1d_row
            !
            !include neighbour in negative x- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+10) = alocont_nn3d(i,j,k,4)
            alocont_colindx(iindx+10)=indx_1d_col
            alocont_rowindx(iindx+10)=indx_1d_row
            !
            !include neighbour in positive y- and positive z-direction
            call conv_indx_3d_to_1d(i,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+11) = alocont_nn3d(i,j,k,26)
            alocont_colindx(iindx+11)=indx_1d_col
            alocont_rowindx(iindx+11)=indx_1d_row
            !
            !include neighbour in positive y- and negative z-direction
            call conv_indx_3d_to_1d(i,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+12) = alocont_nn3d(i,j,k,8)
            alocont_colindx(iindx+12)=indx_1d_col
            alocont_rowindx(iindx+12)=indx_1d_row
            !
            !include neighbour in negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+13) = alocont_nn3d(i,j,k,20)
            alocont_colindx(iindx+13)=indx_1d_col
            alocont_rowindx(iindx+13)=indx_1d_row
            !
            !include neighbour in negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+14) = alocont_nn3d(i,j,k,2)
            alocont_colindx(iindx+14)=indx_1d_col
            alocont_rowindx(iindx+14)=indx_1d_row
            !
            !include neighbour in positive x- and positive y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+15) = alocont_nn3d(i,j,k,18)
            alocont_colindx(iindx+15)=indx_1d_col
            alocont_rowindx(iindx+15)=indx_1d_row
            !
            !include neighbour in positive x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+16) = alocont_nn3d(i,j,k,12)
            alocont_colindx(iindx+16)=indx_1d_col
            alocont_rowindx(iindx+16)=indx_1d_row
            !
            !include neighbour in negative x- and positive y-direction
            call conv_indx_3d_to_1d(i-1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+17) = alocont_nn3d(i,j,k,16)
            alocont_colindx(iindx+17)=indx_1d_col
            alocont_rowindx(iindx+17)=indx_1d_row
            !
            !include neighbour in negative x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+18) = alocont_nn3d(i,j,k,10)
            alocont_colindx(iindx+18)=indx_1d_col
            alocont_rowindx(iindx+18)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+19) = alocont_nn3d(i,j,k,27)
            alocont_colindx(iindx+19)=indx_1d_col
            alocont_rowindx(iindx+19)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and positive z-direction
            call conv_indx_3d_to_1d(i-1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+20) = alocont_nn3d(i,j,k,25)
            alocont_colindx(iindx+20)=indx_1d_col
            alocont_rowindx(iindx+20)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+21) = alocont_nn3d(i,j,k,21)
            alocont_colindx(iindx+21)=indx_1d_col
            alocont_rowindx(iindx+21)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+22) = alocont_nn3d(i,j,k,19)
            alocont_colindx(iindx+22)=indx_1d_col
            alocont_rowindx(iindx+22)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+23) = alocont_nn3d(i,j,k,9)
            alocont_colindx(iindx+23)=indx_1d_col
            alocont_rowindx(iindx+23)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and negative z-direction
            call conv_indx_3d_to_1d(i-1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+24) = alocont_nn3d(i,j,k,7)
            alocont_colindx(iindx+24)=indx_1d_col
            alocont_rowindx(iindx+24)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(1,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+25) = alocont_nn3d(i,j,k,3)
            alocont_colindx(iindx+25)=indx_1d_col
            alocont_rowindx(iindx+25)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i-1,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+26) = alocont_nn3d(i,j,k,1)
            alocont_colindx(iindx+26)=indx_1d_col
            alocont_rowindx(iindx+26)=indx_1d_row
            !
          case(18)
            !
            !--------------------left front boundary (edge)-------------------------
            !------------(influences right and top adjacent points)-----------------
            !
            !include diagonal part
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
            !include direct neighbour in positive x- direction
            call conv_indx_3d_to_1d(i+1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn3d(i,j,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
            !
            !include direct neighbour in negative x-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j,k,nx,ny,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn3d(i,j,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
            !
            !include direct neighbour in positive y- direction
            call conv_indx_3d_to_1d(i,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,17)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
            !!include direct neighbour in negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,11)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
            !
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
            !
            !include neighbour in positive x- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn3d(i,j,k,24)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
            !
            !include neighbour in positive x- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn3d(i,j,k,6)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
            !
            !include neighbour in negative x- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+9) = alocont_nn3d(i,j,k,22)
            alocont_colindx(iindx+9)=indx_1d_col
            alocont_rowindx(iindx+9)=indx_1d_row
            !
            !include neighbour in negative x- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+10) = alocont_nn3d(i,j,k,4)
            alocont_colindx(iindx+10)=indx_1d_col
            alocont_rowindx(iindx+10)=indx_1d_row
            !
            !include neighbour in positive y- and positive z-direction
            call conv_indx_3d_to_1d(i,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+11) = alocont_nn3d(i,j,k,26)
            alocont_colindx(iindx+11)=indx_1d_col
            alocont_rowindx(iindx+11)=indx_1d_row
            !
            !include neighbour in positive y- and negative z-direction
            call conv_indx_3d_to_1d(i,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+12) = alocont_nn3d(i,j,k,8)
            alocont_colindx(iindx+12)=indx_1d_col
            alocont_rowindx(iindx+12)=indx_1d_row
            !
            !include neighbour in negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+13) = alocont_nn3d(i,j,k,20)
            alocont_colindx(iindx+13)=indx_1d_col
            alocont_rowindx(iindx+13)=indx_1d_row
            !
            !include neighbour in negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+14) = alocont_nn3d(i,j,k,2)
            alocont_colindx(iindx+14)=indx_1d_col
            alocont_rowindx(iindx+14)=indx_1d_row
            !
            !include neighbour in positive x- and positive y-direction
            call conv_indx_3d_to_1d(i+1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+15) = alocont_nn3d(i,j,k,18)
            alocont_colindx(iindx+15)=indx_1d_col
            alocont_rowindx(iindx+15)=indx_1d_row
            !
            !include neighbour in positive x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+16) = alocont_nn3d(i,j,k,12)
            alocont_colindx(iindx+16)=indx_1d_col
            alocont_rowindx(iindx+16)=indx_1d_row
            !
            !include neighbour in negative x- and positive y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j+1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+17) = alocont_nn3d(i,j,k,16)
            alocont_colindx(iindx+17)=indx_1d_col
            alocont_rowindx(iindx+17)=indx_1d_row
            !
            !include neighbour in negative x- and negative y-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,ny-1,k,nx,ny,indx_1d_row)
            alocont_data(iindx+18) = alocont_nn3d(i,j,k,10)
            alocont_colindx(iindx+18)=indx_1d_col
            alocont_rowindx(iindx+18)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and positive z-direction
            call conv_indx_3d_to_1d(i+1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+19) = alocont_nn3d(i,j,k,27)
            alocont_colindx(iindx+19)=indx_1d_col
            alocont_rowindx(iindx+19)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j+1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+20) = alocont_nn3d(i,j,k,25)
            alocont_colindx(iindx+20)=indx_1d_col
            alocont_rowindx(iindx+20)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+21) = alocont_nn3d(i,j,k,21)
            alocont_colindx(iindx+21)=indx_1d_col
            alocont_rowindx(iindx+21)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and positive z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,ny-1,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+22) = alocont_nn3d(i,j,k,19)
            alocont_colindx(iindx+22)=indx_1d_col
            alocont_rowindx(iindx+22)=indx_1d_row
            !
            !include neighbour in positive x-, positive y- and negative z-direction
            call conv_indx_3d_to_1d(i+1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+23) = alocont_nn3d(i,j,k,9)
            alocont_colindx(iindx+23)=indx_1d_col
            alocont_rowindx(iindx+23)=indx_1d_row
            !
            !include neighbour in negative x-, positive y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,j+1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+24) = alocont_nn3d(i,j,k,7)
            alocont_colindx(iindx+24)=indx_1d_col
            alocont_rowindx(iindx+24)=indx_1d_row
            !
            !include neighbour in positive x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(i+1,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+25) = alocont_nn3d(i,j,k,3)
            alocont_colindx(iindx+25)=indx_1d_col
            alocont_rowindx(iindx+25)=indx_1d_row
            !
            !include neighbour in negative x-, negative y- and negative z-direction (periodic boundary condition)
            call conv_indx_3d_to_1d(nx-1,ny-1,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+26) = alocont_nn3d(i,j,k,1)
            alocont_colindx(iindx+26)=indx_1d_col
            alocont_rowindx(iindx+26)=indx_1d_row
            !
            !---------------------right and back boundaries-------------------------
            !-----------------(only diagonal part, will be updated------------------
            !------------according to periodic boundary condition anyways)----------
            !
          case(13,15,19,21,23,24,25,26,27)
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            alocont_data(iindx) = alocont_nn3d(i,j,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn3d(i,j,k,14)
            !
          case(1,41,51,61,71,81,91,101,111,121,141,161,171,181,201)!,211,221,191,131,251,271,261,151,231,241)
            !
            !--------------------bottom boundary (non-periodic)--------------------
            !-----------------------(only in z-direction)---------------------------
            !
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            !include direct neighbour in positive z- direction
            call conv_indx_3d_to_1d(i,j,k+1,nx,ny,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn3d(i,j,k,23)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
            !
          case(2,42,52,62,72,82,92,102,112,122,142,162,172,182,202)!,212,222,192,132,252,272,262,152,232,242)
            !
            !---------------------top boundary (non-periodic)-----------------------
            !-----------------------(only in z-direction)---------------------------
            !
            call conv_indx_3d_to_1d(i,j,k,nx,ny,indx_1d_col)
            !include direct neighbour in negative z-direction
            call conv_indx_3d_to_1d(i,j,k-1,nx,ny,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn3d(i,j,k,5)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            !
          case default
          end select
          !
          iindx=iindx+27
          !
        enddo
      enddo
    enddo
    !
    !write(*,*)
    !write(*,*) 'maximum number of used neighbours for ALO calculation', ndiags_max-1
    !write(*,*) 'minimum number of used neighbours for ALO calculation', ndiags_min-1
    !write(*,*)
    !return

    !i=2
    !do k=1, nz
    !   write(*,*) alocont_nn3d(2,k,22), alocont_nn3d(3,k,22), alocont_nn3d(4,k,22), alocont_nn3d(5,k,22), alocont_nn3d(nx-2,k,22), alocont_nn3d(nx-1,k,22), alocont_nn3d(nx,k,22)
    !enddo
    !
    !stop 'go on in calc_alocont'


  end subroutine calc_alocont_nn3d_coo


end module mod_scont_new3d
