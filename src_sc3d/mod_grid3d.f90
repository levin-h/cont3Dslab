module mod_grid3d
  !
  use prog_type
  use fund_const
  use omp_lib
  !
  implicit none
  !
  !----------------------options to set up the grid-----------------------
  !
  integer(i4b) :: opt_gridxyz
  !opt_gridxyz = 0   linear spacing in x,y,z
  !opt_gridxyz = 1   logarithmic spacing in z, linear in x,y
  !
  !
  !-----------------------------data--------------------------------------
  !
  integer(i4b) :: nx, ny, nz
  !
  real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
  !
  real(dp), dimension(:), allocatable :: x, y, z
  integer, dimension(:,:,:), allocatable :: imask3d, imaskb3d
  !
  real(dp), dimension(:,:,:), allocatable :: rho3d, tgas3d, trad3d, velx3d, vely3d, velz3d, &
       eps_cont3d, opac3d, bnue3d, scont3d, normalization3d
  real(dp), dimension(:,:,:), allocatable :: int3d, mint3d, mintbar3d, fcontx3d, fconty3d, fcontz3d, &
                                             kcontxx3d, kcontxy3d, kcontxz3d, &
                                             kcontyy3d, kcontyz3d, kcontzz3d
  !
  !for all alo terms
  real(dp), dimension(:,:,:,:), allocatable :: alocont_nn3d, alocont_o_nn3d
  real(dp), dimension(:), allocatable :: alocont_data, alocont_data_diag, aloline_data, aloline_data_diag
  integer(i4b), dimension(:), allocatable :: alocont_colindx, alocont_rowindx, aloline_colindx, aloline_rowindx
  !
  !for paralleliztion
  real(dp), dimension(:,:,:,:), allocatable :: alocont_nn3d_tmp
  real(dp), dimension(:,:,:), allocatable :: mint3d_tmp, normalization3d_tmp, opac3d_tmp, &
       fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp, &
       kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp, &
       kcontyy3d_tmp, kcontyz3d_tmp, &
       kcontzz3d_tmp
  !$omp threadprivate(int3d, alocont_o_nn3d, alocont_nn3d_tmp, normalization3d_tmp, mint3d_tmp)
  !$omp threadprivate(fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp)
  !$omp threadprivate(kcontxx3d_tmp, kcontxy3d_tmp, kcontxz3d_tmp)
  !$omp threadprivate(kcontyy3d_tmp, kcontyz3d_tmp)
  !$omp threadprivate(kcontzz3d_tmp)
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
contains
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine allocate_global3d(verbose)
    !
    ! ... arguments
    logical, intent(in), optional :: verbose
    !
    ! ... local characters
    logical :: ver = .false.
    !
    integer(i4b) :: err
    !
    !
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '-------------------------allocating all global 3d arrays-----------------------'
    if(ver) write(*,*)
    !
    if(allocated(mint3d)) deallocate(mint3d)
    if(allocated(scont3d)) deallocate(scont3d)
    if(allocated(normalization3d)) deallocate(normalization3d)
    if(allocated(alocont_nn3d)) deallocate(alocont_nn3d)
    if(allocated(bnue3d)) deallocate(bnue3d)
    if(allocated(opac3d)) deallocate(opac3d)
    if(allocated(eps_cont3d)) deallocate(eps_cont3d)
    if(allocated(fcontx3d)) deallocate(fcontx3d)
    if(allocated(fconty3d)) deallocate(fconty3d)
    if(allocated(fcontz3d)) deallocate(fcontz3d)
    if(allocated(kcontxx3d)) deallocate(kcontxx3d)
    if(allocated(kcontxy3d)) deallocate(kcontxy3d)
    if(allocated(kcontxz3d)) deallocate(kcontxz3d)
    if(allocated(kcontyy3d)) deallocate(kcontyy3d)
    if(allocated(kcontyz3d)) deallocate(kcontyz3d)
    if(allocated(kcontzz3d)) deallocate(kcontzz3d)

    
    allocate(mint3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error in allocate_global3d'
    allocate(scont3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error in allocate_global3d'
    allocate(normalization3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error in allocate_global3d'
    allocate(alocont_nn3d(nx,ny,nz,27), stat=err)
       if(err.ne.0) stop 'allocation error in allocate_global3d'
    allocate(bnue3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error in allocate_global3d'
    allocate(opac3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error in allocate_global3d'
    allocate(eps_cont3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'allocation error in allocate_global3d'
    allocate(fcontx3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
    allocate(fconty3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
    allocate(fcontz3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
    allocate(kcontxx3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
    allocate(kcontxy3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
    allocate(kcontxz3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
    allocate(kcontyy3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
    allocate(kcontyz3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
    allocate(kcontzz3d(nx,ny,nz), stat=err)
       if(err.ne.0) stop 'error in allocate_global3d: allocation'
!
end subroutine allocate_global3d
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine gridxyz(verbose)
    !
    !input grid has information (along x and y direction)
    !
    !   y
    !   |  
    !   |
    !   ------x  
    !
    !   41----42----43----44
    !   |     |     |     |
    !   31----32----33----34
    !   |     |     |     |
    !   21----22----23----24
    !   |     |     |     |
    !   11----12----13----14  
    !
    !for efficient radiative transfer scheme with periodic boundary conditions
    !   RT grid requires information to be stored on
    !
    !   11----12----13----14----11
    !   |     |     |     |     |
    !   41----42----43----44----41
    !   |     |     |     |     |
    !   31----32----33----34----31
    !   |     |     |     |     |
    !   21----22----23----24----21
    !   |     |     |     |     |
    !   11----12----13----14----11
    !
    !thus, need to assume spatial distance between x_11-x_14 calculated as 1/2*(x_14-x_13)+1/2*(x_12-x_11)
    !                     spatial distance between y_11-x_41 calculated as 1/2*(y_41-y_31)+1/2*(y_21-y_11)    
    !
    !since the x and y grids are then not exactly symmetric, nx and ny shall be even numbers (in order to include grid point
    !at (xmax-xmin)/2  and (ymax-ymin)/2
    !
    logical, intent(in), optional :: verbose
    !
    ! ... local characters
    logical :: ver = .false.
    !
    integer(i4b) :: i, j, k, err
    !
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '-----------------------creating x, y and z grid---------------------------------'
    if(ver) write(*,*)
    !
    select case(opt_gridxyz)
       case(0)
          call calc_gridxyz00(verbose=ver)
       case(1)
         call calc_gridxyz01(verbose=ver)
      case default
         stop 'error in gridxyz: opt_gridxyz not properly specified'
   end select
   !
   !------------------------masks for computational domain-----------------
   !
   allocate(imask3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error in gridxyz: allocation'
   allocate(imaskb3d(nx,ny,nz), stat=err)
      if(err.ne.0) stop 'error in gridxyz: allocation'
   !
   !standard points within computational domain are labelled 1
   imask3d=1
   do i=1, nx
      do j=1, ny
         imask3d(i,j,1)=0
         imask3d(i,j,2)=0
         imask3d(i,j,nz-1)=0
         imask3d(i,j,nz)=0
      enddo
   enddo
   !
   !standard points within computational domain are labelled 3
   imaskb3d=3
   !
   !
   !
   do i=1, nx
      do j=1, ny
         imaskb3d(i,j,1)=0     !lower most boundary is labelled 0
         imaskb3d(i,j,nz)=0    !upper most boundary is labelled 0
         imaskb3d(i,j,2)=1     !next lower most boundary is labelled 1
         imaskb3d(i,j,nz-1)=2  !next upper most boundary is labelled 2
      enddo
   enddo
   !
   !
   !
   do k=2, nz-1
      !
      !left boundary (xz-plane) is labelled 12 (on lower and upper z-boundary: 121,122)
      !right boundary (xz-plane) is labelled 13 (on lower and upper z-boundary: 131,132)
      !left adjacent to boundary is labelled 4 (on lower and upper z-boundary: 41,42)
      !right adjacent to boundary is labelled 5 (on lower and upper z-boundary: 51,52)
      do j=3, ny-2
         if(k.eq.2) then
            imaskb3d(1,j,k)=121
            imaskb3d(nx,j,k)=131
            imaskb3d(2,j,k)=41
            imaskb3d(nx-1,j,k)=51
         elseif(k.eq.nz-1) then
            imaskb3d(1,j,k)=122
            imaskb3d(1,j,k)=132
            imaskb3d(2,j,k)=42
            imaskb3d(nx-1,j,k)=52
         else
            imaskb3d(1,j,k)=12
            imaskb3d(nx,j,k)=13
            imaskb3d(2,j,k)=4
            imaskb3d(nx-1,j,k)=5
         endif
      enddo
      !front boundary (yz-plane) is labelled 14 (on lower and upper z-boundary: 141,142)
      !back boundary (yz-plane) is labelled 15 (on lower and upper z-boundary: 151,152)
      !front adjacent to boundary is labelled 6 (on lower and upper z-boundary: 61,62)
      !back adjacent to boundary is labelled 7 (on lower and upper z-boundary: 71,72)
      do i=3, nx-2
         if(k.eq.2) then
            imaskb3d(i,1,k)=141
            imaskb3d(i,ny,k)=151
            imaskb3d(i,2,k)=61
            imaskb3d(i,ny-1,k)=71
         elseif(k.eq.nz-1) then
            imaskb3d(i,1,k)=142
            imaskb3d(i,ny,k)=152
            imaskb3d(i,2,k)=62
            imaskb3d(i,ny-1,k)=72
         else
            imaskb3d(i,1,k)=14
            imaskb3d(i,ny,k)=15
            imaskb3d(i,2,k)=6
            imaskb3d(i,ny-1,k)=7
         endif
      enddo
      !
      !left adjacent front adjacent to boundary is labelled 8  (on lower and upper z-boundary: 81, 82)
      !right adjacent front adjacent to boundary is labelled 9  (on lower and upper z-boundary: 91, 92)
      !left adjacent back adjacent to boundary is labelled 10  (on lower and upper z-boundary: 101, 102)
      !right adjacent back adjacent to boundary is labelled 11  (on lower and upper z-boundary: 111, 112)
      !left and front adjacent to boundary is labelled 16  (on lower and upper z-boundary: 161, 162)
      !left adjacent and front is labelled 17  (on lower and upper z-boundary: 171, 172)
      !left and front is labelled 18  (on lower and upper z-boundary: 181, 182)
      !right and front adjacent to boundary is labelled 19  (on lower and upper z-boundary: 191, 192)
      !left adjacent and front is labelled 20  (on lower and upper z-boundary: 201, 202)
      !right and front is labelled 21  (on lower and upper z-boundary: 211, 212)
      !left and back adjacent to boundary is labelled 22  (on lower and upper z-boundary: 221, 222)
      !left adjacent and back is labelled 23  (on lower and upper z-boundary: 231, 232)
      !left and back is labelled 24  (on lower and upper z-boundary: 241, 242)
      !right and back adjacent to boundary is labelled 25  (on lower and upper z-boundary: 251, 252)
      !right adjacent and back is labelled 26  (on lower and upper z-boundary: 261, 262)
      !right and back is labelled 27  (on lower and upper z-boundary: 271, 272)   
      if(k.eq.2) then
         imaskb3d(2,2,k)=81
         imaskb3d(nx-1,2,k)=91
         imaskb3d(2,ny-1,k)=101
         imaskb3d(nx-1,ny-1,k)=111
         imaskb3d(1,2,k)=161
         imaskb3d(2,1,k)=171
         imaskb3d(1,1,k)=181
         imaskb3d(nx,2,k)=191
         imaskb3d(nx-1,1,k)=201
         imaskb3d(nx,1,k)=211
         imaskb3d(1,ny-1,k)=221
         imaskb3d(2,ny,k)=231
         imaskb3d(1,ny,k)=241
         imaskb3d(nx,ny-1,k)=251
         imaskb3d(nx-1,ny,k)=261
         imaskb3d(nx,ny,k)=271
      elseif(k.eq.nz-1) then
         imaskb3d(2,2,k)=82
         imaskb3d(nx-1,2,k)=92
         imaskb3d(2,ny-1,k)=102
         imaskb3d(nx-1,ny-1,k)=112
         imaskb3d(1,2,k)=162
         imaskb3d(2,1,k)=172
         imaskb3d(1,1,k)=182
         imaskb3d(nx,2,k)=192
         imaskb3d(nx-1,1,k)=202
         imaskb3d(nx,1,k)=212
         imaskb3d(1,ny-1,k)=222
         imaskb3d(2,ny,k)=232
         imaskb3d(1,ny,k)=242
         imaskb3d(nx,ny-1,k)=252
         imaskb3d(nx-1,ny,k)=262
         imaskb3d(nx,ny,k)=272
      else
         imaskb3d(2,2,k)=8
         imaskb3d(nx-1,2,k)=9
         imaskb3d(2,ny-1,k)=10
         imaskb3d(nx-1,ny-1,k)=11
         imaskb3d(1,2,k)=16
         imaskb3d(2,1,k)=17
         imaskb3d(1,1,k)=18
         imaskb3d(nx,2,k)=19
         imaskb3d(nx-1,1,k)=20
         imaskb3d(nx,1,k)=21
         imaskb3d(1,ny-1,k)=22
         imaskb3d(2,ny,k)=23
         imaskb3d(1,ny,k)=24
         imaskb3d(nx,ny-1,k)=25
         imaskb3d(nx-1,ny,k)=26
         imaskb3d(nx,ny,k)=27
      endif
      !
   enddo!

   !
   !
   !
 end subroutine gridxyz
 !
 !-----------------------------------------------------------------------
 !-----------------------------------------------------------------------
 !-----------------------------------------------------------------------
 !
 !
 subroutine calc_gridxyz00(verbose)
    !
    !input grid has information (along x and y direction)
    !
    !   y
    !   |  
    !   |
    !   ------x  
    !
    !   41----42----43----44
    !   |     |     |     |
    !   31----32----33----34
    !   |     |     |     |
    !   21----22----23----24
    !   |     |     |     |
    !   11----12----13----14  
    !
    !for efficient radiative transfer scheme with periodic boundary conditions
    !   RT grid requires information to be stored on
    !
    !   11----12----13----14----11
    !   |     |     |     |     |
    !   41----42----43----44----41
    !   |     |     |     |     |
    !   31----32----33----34----31
    !   |     |     |     |     |
    !   21----22----23----24----21
    !   |     |     |     |     |
    !   11----12----13----14----11
    !
    !thus, need to assume spatial distance between x_11-x_14 calculated as 1/2*(x_14-x_13)+1/2*(x_12-x_11)
    !                     spatial distance between y_11-x_41 calculated as 1/2*(y_41-y_31)+1/2*(y_21-y_11)
    !
    !since the x and y grids are then not exactly symmetric, nx and ny shall be even numbers (in order to include grid point!at (xmax-xmin)/2  and (ymax-ymin)/2
    !
    !
    logical, intent(in), optional :: verbose
    !
    ! ... local characters
    logical :: ver = .false.
    !
    integer(i4b) :: i, j, k, err
    !
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) 'linear spacing in x, y, z'
    if(ver) write(*,*)
    !
    !ensure even number of nx and ny
    if(modulo(nx,2).ne.0) nx=nx+1
    if(modulo(ny,2).ne.0) ny=ny+1
    !
    !
    !
    allocate(x(nx), stat=err)
       if(err.ne.0) stop 'error in calc_gridxyz00: allocation'
    allocate(y(ny), stat=err)
       if(err.ne.0) stop 'error in calc_gridxyz00: allocation'
    allocate(z(nz), stat=err)
       if(err.ne.0) stop 'error in calc_gridxyz00: allocation'
    !
    !----------------------equidistant grid in x----------------------------
    !
    if(nx.le.3) stop 'error in calc_gridxyz00: nx needs to be larger than 3'
    do i=1, nx-1
       x(i) = xmin + (i-1)*(xmax-xmin)/(nx-2)
    enddo
    !ghost zones
    x(nx) = x(nx-1) + half*(x(nx-1)-x(nx-2)) + half*(x(2)-x(1))
    !
    !----------------------equidistant grid in y----------------------------
    !
    if(ny.le.3) stop 'error in calc_gridxyz00: ny needs to be larger than 3'   
    !
    do i=1, ny-1
       y(i) = ymin + (i-1)*(ymax-ymin)/(ny-2)
    enddo
    !ghost zones
    y(ny) = y(ny-1) + half*(y(ny-1)-y(ny-2)) + half*(y(2)-y(1))
    !
    !
    !do i=2, ny-1
    !   y(i) = ymin + (i-2)*(ymax-ymin)/(ny-3)
    !enddo
    !!ghost zones
    !y(1)=two*y(2)-y(3)
    !y(ny)=two*y(ny-1)-y(ny-2)
    !
    !----------------------equidistant grid in z----------------------------
    !
    if(nz.le.5) stop 'error in calc_gridxyz00: nz needs to be larger than 5'   
    do i=3, nz-2
       z(i) = zmin + (i-3)*(zmax-zmin)/(nz-5)
    enddo
    !ghost zones
    z(2)=two*z(3)-z(4)
    z(1)=two*z(2)-z(3)
    z(nz-1)=two*z(nz-2)-z(nz-3)
    z(nz)=two*z(nz-1)-z(nz-2)
    !
  end subroutine calc_gridxyz00
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine calc_gridxyz01(verbose)
    !
    !input grid has information (along x and y direction)
    !
    !   y
    !   |
    !   |
    !   ------x
    !  
    !   41----42----43----44
    !   |     |     |     |
    !   31----32----33----34
    !   |     |     |     |
    !   21----22----23----24
    !   |     |     |     |
    !   11----12----13----14
    !
    !for efficient radiative transfer scheme with periodic boundary conditions
    !   RT grid requires information to be stored on  
    !
    !   11----12----13----14----11
    !   |     |     |     |     |
    !   41----42----43----44----41
    !   |     |     |     |     |
    !   31----32----33----34----31
    !   |     |     |     |     |
    !   21----22----23----24----21
    !   |     |     |     |     |
    !   11----12----13----14----11
    !
    !thus, need to assume spatial distance between x_11-x_14 calculated as 1/2*(x_14-x_13)+1/2*(x_12-x_11)
    !                     spatial distance between y_11-x_41 calculated as 1/2*(y_41-y_31)+1/2*(y_21-y_11)
    !
    !since the x and y grids are then not exactly symmetric, nx and ny shall be even numbers (in order to include grid point
    !at (xmax-xmin)/2  and (ymax-ymin)/2
    !
    logical, intent(in), optional :: verbose
    !
    ! ... local characters
    logical :: ver = .false.
    !
    integer(i4b) :: i, j, k, err
    real(dp) :: zmin_log, zmax_log, zshift, del
    !
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) 'linear spacing in x, y'
    if(ver) write(*,*) 'logarithmic spacing in z in [1,zmax-zmin]'
    if(ver) write(*,*)
    !
    !ensure even number of nx and ny
    if(modulo(nx,2).ne.0) nx=nx+1
    if(modulo(ny,2).ne.0) ny=ny+1
    !
    !
    !
    allocate(x(nx), stat=err)
       if(err.ne.0) stop 'error in calc_gridxyz00: allocation'
    allocate(y(ny), stat=err)
       if(err.ne.0) stop 'error in calc_gridxyz00: allocation'
    allocate(z(nz), stat=err)
       if(err.ne.0) stop 'error in calc_gridxyz00: allocation'
    !
    !----------------------equidistant grid in x----------------------------
    !
    if(nx.le.3) stop 'error in calc_gridxyz00: nx needs to be larger than 3'
    do i=1, nx-1
       x(i) = xmin + (i-1)*(xmax-xmin)/(nx-2)
    enddo
    !ghost zones
    x(nx) = x(nx-1) + half*(x(nx-1)-x(nx-2)) + half*(x(2)-x(1))
    !
    !----------------------equidistant grid in y----------------------------
    !
    if(ny.le.3) stop 'error in calc_gridxyz00: ny needs to be larger than 3'   
    !
    do i=1, ny-1
       y(i) = ymin + (i-1)*(ymax-ymin)/(ny-2)
    enddo
    !ghost zones
    y(ny) = y(ny-1) + half*(y(ny-1)-y(ny-2)) + half*(y(2)-y(1))
    !
    !
    !do i=2, ny-1
    !   y(i) = ymin + (i-2)*(ymax-ymin)/(ny-3)
    !enddo
    !!ghost zones
    !y(1)=two*y(2)-y(3)
    !y(ny)=two*y(ny-1)-y(ny-2)
    !
    !----------------------logarithmic grid in z----------------------------
    !
    if(nz.le.5) stop 'error in calc_gridxyz00: nz needs to be larger than 5'
    !
    zshift = one-zmin
    zmin_log = zmin + zshift
    zmax_log = zmax + zshift
    !
    del=log10(zmax_log/zmin_log)/float(nz-5)
    z(3) = zmin_log
    do i=4, nz-2
       z(i) = z(i-1)*10.d0**del
    enddo
    !shift everything back to original zmin, zmax
    z = z - zshift
    z(3) = zmin
    z(nz-2) = zmax
    !ghost zones
    z(2)=two*z(3)-z(4)
    z(1)=two*z(2)-z(3)
    z(nz-1)=two*z(nz-2)-z(nz-3)
    z(nz)=two*z(nz-1)-z(nz-2)
    !
    !
    !
  end subroutine calc_gridxyz01
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  subroutine check_grid3d(nueindx, verbose)
    !
    !test the upwind and downwind delta-tau steps
    !such that alo-coefficients are not 'extreme'
    !needs to be done for each nue  
    !
    use options, only: opt_method
    use mod_integ1d, only: integ1d_tau_ud
    use mod_interp2d, only: interpol2d_4p_lin
    use mod_interp1d, only: interpol_yp, find_index
    use mod_grid, only: grid_log
    use mod_math, only: calc_dtau_crit
    !
    ! ... arguments
    integer(i4b), intent(in) :: nueindx
    logical, intent(in), optional :: verbose
    !
    ! ... local characters
    logical :: ver = .false.
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, oindx
    integer(i4b) :: iim1, iip1, jjm1, jjp1, kkm1, kkp1, llm2, llm1, ll, llp1
    real(dp) :: dels_u, dels_d, delt_u, delt_d
    real(dp) :: dels_xyu, dels_xzu, dels_yzu, dels_xyd, dels_xzd, dels_yzd
    real(dp) :: nn_x, nn_y, nn_z
    real(dp) :: x_u, y_u, z_u, x_d, y_d, z_d
    real(dp) :: opac_u, opac_p, opac_d
    real(dp) :: delt_crita, delt_critb, delt_critc
    real(dp) :: eps, dtumax, dtumin, dtdmax

    !
    ! ... local arrays
    integer(i4b), parameter :: nomega=7, ntest=101
    real(dp), dimension(nomega) :: n_x, n_y, n_z, norm
    real(dp), dimension(ntest) :: delt_crita_arr, delt_critb_arr, delt_critc_arr, delt_u_arr
    !
    ! ... local funcions
    !
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '-----------------------testing the delta-tau steps-----------------------------'
    if(ver) write(*,*)
    !
    !linear interpolations give always stable alo-coefficients
    if(opt_method.eq.4.or.&
         opt_method.eq.6) return
    !
    !----------------------directions to be tested--------------------------
    !
    eps = 1.d-5
    n_x = (/ one, eps, eps, one, eps, one, one /)
    n_y = (/ eps, one, eps, one, one, eps, one /)
    n_z = (/ eps, eps, one, eps, one, one, one /)
    norm = sqrt(n_x**2+n_y**2+n_z**2)
    n_x = n_x/norm
    n_y = n_y/norm
    n_z = n_z/norm
    !
    !-----------------grid for critical delta-tau steps---------------------
    !
    dtumin=1.d-4
    dtumax=1.d7
    !
    call grid_log(dtumin, dtumax, ntest, delt_u_arr)
    !
    do i=1, ntest   
       call calc_dtau_crit(delt_u_arr(i), delt_crita, delt_critb, delt_critc)
       delt_crita_arr(i)=delt_crita
       delt_critb_arr(i)=delt_critb
       delt_critc_arr(i)=delt_critc
    enddo
    !
    open(1,file='outputFILES_TEST/delt_crit.dat', form='formatted')
       if(ver) write(1,'(4a20)') 'dtau_u', 'dtaud_crita', 'dtaud_critb', 'dtaud_critc'
       do i=1, ntest
          if(ver) write(1,'(4es20.8)') delt_u_arr(i), delt_crita_arr(i), delt_critb_arr(i), delt_critc_arr(i)
       enddo
    close(1)
    !
    !-----------------------------------------------------------------------
    !
    do oindx=1, nomega
       !
       nn_x=n_x(oindx)
       nn_y=n_y(oindx)
       nn_z=n_z(oindx)
       !
       do i=3, nx-2
          do j=3, ny-2
             do k=3, nz-2
                !
                !-------------------------check direction-------------------------------
                !
                iim1=i-1
                jjm1=j-1
                kkm1=k-1
                iip1=i+1
                jjp1=j+1
                kkp1=k+1
                !
                !calculate distance of upwind point to previous xy-plane, xz-plane and yz-plane
                dels_xyu=(z(k)-z(kkm1))/nn_z
                dels_xzu=(y(j)-y(jjm1))/nn_y
                dels_yzu=(x(i)-x(iim1))/nn_x
                dels_u=min(dels_xyu,dels_xzu,dels_yzu)
                !
                !calculate distance of downwind point to (next) xy-plane, xz-plane and yz-plane
                dels_xyd=(z(kkp1)-z(k))/nn_z
                dels_xzd=(y(jjp1)-y(j))/nn_y
                dels_yzd=(x(iip1)-x(i))/nn_x
                dels_d=min(dels_xyd,dels_xzd,dels_yzd)
                !
                !----------------------------local point--------------------------------
                !
                opac_p=opac3d(i,j,k)
                !
                !----------------------------upwind point-------------------------------
                
                if(dels_xyu.eq.dels_u) then
                   !intersection with x-y plane on level k-gamma
                   x_u = x(i) - dels_u*nn_x
                   y_u = y(j) - dels_u*nn_y
                   z_u = z(kkm1)
                   !
                   opac_u = interpol2d_4p_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                        opac3d(iim1,j,kkm1), opac3d(i,j,kkm1), &
                        x(iim1), x(i), y(jjm1), y(j), x_u, y_u)
                   !
                elseif(dels_xzu.eq.dels_u) then
                   !intersection with x-z plane on level j-beta
                   x_u = x(i) - dels_u*nn_x
                   y_u = y(jjm1)
                   z_u = z(k) - dels_u*nn_z
                   !
                   opac_u = interpol2d_4p_lin(opac3d(iim1,jjm1,kkm1), opac3d(i,jjm1,kkm1), &
                        opac3d(iim1,jjm1,k), opac3d(i,jjm1,k), &
                        x(iim1), x(i), z(kkm1), z(k), x_u, z_u)
                   !
                elseif(dels_yzu.eq.dels_u) then
                   !intersection with y-z plane on level i-alpha
                   x_u = x(iim1)
                   y_u = y(j) - dels_u*nn_y
                   z_u = z(k) - dels_u*nn_z
                   !
                   opac_u = interpol2d_4p_lin(opac3d(iim1,jjm1,kkm1), opac3d(iim1,j,kkm1), &
                        opac3d(iim1,jjm1,k), opac3d(iim1,j,k), &
                        y(jjm1), y(j), z(kkm1), z(k), y_u, z_u)
                   !
                else
                   if(ver) write(*,'(4es20.8,2l4)') dels_u, dels_xzu, dels_xyu, dels_yzu, dels_u.eq.dels_xzu, dels_xyu.eq.dels_u
                   stop 'error in check_grid3d: invalid dels_u'
                endif
                !
                !---------------------------downwind point------------------------------
                !
                if(dels_xyd.eq.dels_d) then
                   !intersection with x-y plane on level k+gamma
                   x_d = x(i) + dels_d*nn_x
                   y_d = y(j) + dels_d*nn_y
                   z_d = z(kkp1)
                   !
                   opac_d = interpol2d_4p_lin(opac3d(i,j,kkp1), opac3d(iip1,j,kkp1), &
                        opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                        x(i), x(iip1), y(j), y(jjp1), x_d, y_d)
                   !
                elseif(dels_xzd.eq.dels_d) then
                   !intersection with x-z plane on level j+beta
                   x_d = x(i) + dels_d*nn_x
                   y_d = y(jjp1)
                   z_d = z(k) + dels_d*nn_z
                   !
                   opac_d = interpol2d_4p_lin(opac3d(i,jjp1,k), opac3d(iip1,jjp1,k), &
                        opac3d(i,jjp1,kkp1), opac3d(iip1,jjp1,kkp1), &
                        x(i), x(iip1), z(k), z(kkp1), x_d, z_d)
                   !
                elseif(dels_yzd.eq.dels_d) then
                   !intersection with y-z plane on level i+alpha
                   x_d = x(iip1)
                   y_d = y(j) + dels_d*nn_y
                   z_d = z(k) + dels_d*nn_z
                   !
                   opac_d = interpol2d_4p_lin(opac3d(iip1,j,k), opac3d(iip1,jjp1,k), &
                        opac3d(iip1,j,kkp1), opac3d(iip1,jjp1,kkp1), &
                        y(j), y(jjp1), z(k), z(kkp1), y_d, z_d)
                   !
                else
                   if(ver) write(*,'(4es20.8,l4)') dels_d, dels_xzd, dels_xyd, dels_yzd, dels_d.eq.dels_xzd
                   stop 'error in check_grid3d: invalid dels_d'
                endif
                !
                !calculate upwind and downdwind delta-tau steps and compare with critical value
                call integ1d_tau_ud(dels_u, dels_d, opac_u, opac_p, opac_d, delt_u, delt_d)
                call find_index(delt_u, delt_u_arr, ntest, llm2, llm1, ll, llp1)
                delt_crita = interpol_yp(delt_u_arr(llm1), delt_u_arr(ll), delt_crita_arr(llm1), delt_crita_arr(ll), delt_u)
                delt_critb = interpol_yp(delt_u_arr(llm1), delt_u_arr(ll), delt_critb_arr(llm1), delt_critb_arr(ll), delt_u)
                delt_critc = interpol_yp(delt_u_arr(llm1), delt_u_arr(ll), delt_critc_arr(llm1), delt_critc_arr(ll), delt_u)
                !            call calc_dtau_crit(delt_u, delt_crita, delt_critb, delt_critc)
                !            if(ver) write(*,*) i, k, delt_u, delt_d, delt_crita, delt_critb, delt_critcy
                if(delt_d.gt.delt_crita) stop 'error in check_grid3d: delt_d > delt_crita => increase resolution'
                if(delt_d.lt.delt_critb) stop 'error in check_grid3d: delt_d < delt_critb => decrease resolution'
                if(delt_d.lt.delt_critc) stop 'error in check_grid3d: delt_d < delt_critc => decrease resolution'
                !            
             enddo
          enddo
       enddo
    enddo
    !
    !
    !
  end subroutine check_grid3d
  !

end module mod_grid3d
