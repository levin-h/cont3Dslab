module mod_frequencies
  !
  use prog_type
  use fund_const
  !
  implicit none
  !
  integer(i4b) :: opt_grey
  !opt_grey = 0   if a number of frequency points shall be used (to be implemented)
  !opt_grey = 1   if frequency dependent transfer at one specific frequency bin
  !opt_grey = 2   if grey approximation with frequency integrated values
  !
  !data
  integer(i4b) :: nnue
  real(dp) :: xnue0
  !
  ! ... logicals
  logical :: lfreqint
  !set to true if frequency integrated radiative transfer used
  !
  ! ... arrays
  real(dp), dimension(:), allocatable :: nodes_nue, weight_nue, xic1_nue
  !
  !-----------------------------------------------------------------------
  !
contains
  !
  subroutine calcnodes_nue(verbose)
    !
    logical, intent(in), optional :: verbose
    !
    ! ... local characters
    logical :: ver = .false.
    
    ! ... local scalars
    integer(i4b) :: i, err
    real(dp) :: rho, dtau
    real(dp) :: theta_min, theta_max
    !
    ! ... local logicals
    !
    ! ... local functions
    real(dp) :: opac_thomson
    !
    !
    if(present(verbose)) ver = verbose
    !
    if(ver) write(*,*) '---------------------------creating frequency grid-----------------------------'
    if(ver) write(*,*)
    !
    allocate(nodes_nue(nnue), stat=err)
       if(err.ne.0) stop 'error: allocation in calcnodes_nue'
    allocate(weight_nue(nnue), stat=err)
       if(err.ne.0) stop 'error: allocation in calcnodes_nue'
    !
    if(nnue.eq.1) then
       nodes_nue=xnue0
       weight_nue=one
    else
       stop 'todo: calcnodes_nue'
    endif
    !
  end subroutine calcnodes_nue


end module mod_frequencies

