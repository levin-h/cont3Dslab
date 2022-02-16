module mod_bcondition3d
  use prog_type
  use fund_const

  implicit none

  integer(i4b) :: opt_bcondition
  real(dp), dimension(:), allocatable :: xic1_nue, xic2_nue
  real(dp), dimension(:,:,:), allocatable :: xic1_3d, xic2_3d

  contains


    subroutine setup_bcondition3d(nueindx, verbose)
      !
      use mod_grid3d, only: nx, ny, bnue3d, x, y, z
      ! ... arguments
      integer(i4b), intent(in) :: nueindx
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



      select case(opt_bcondition)
         case(0)
           !set xic1 and xic2 to input parameters
            xic1_3d = xic1_nue(nueindx)
            xic2_3d = xic2_nue(nueindx)
         case(1)
            !set xic1 = bnue (already set beforehand at this frequency index)
            do i=1, nx
               do j=1, ny
                  xic1_3d(i,j,1) = bnue3d(i,j,1)
                  xic1_3d(i,j,2) = bnue3d(i,j,2)
               enddo
            enddo
            stop 'todo: set xic2_3d in setup_bcondition3d'
         case default
            stop 'error in setup_bcondition3d: opt_bcondition not properly set'
      end select


      !
      i=nx/2+1
      j=ny/2+1
      if(ver) write(*,'(a20,i5)') 'frequency point', nueindx
      if(ver) write(*,'(a20, 2es20.8)') 'at (x,y)= ', x(i), y(j)
      if(ver) write(*,'(10a20)') 'z [r_star]', 'xic1 [cgs]', 'xic2 [cgs]'
      if(ver) write(*,'(10es20.8)') z(1), xic1_3d(i,j,1), xic2_3d(i,j,1)
      if(ver) write(*,'(10es20.8)') z(2), xic1_3d(i,j,2), xic2_3d(i,j,2)
      if(ver) write(*,*)
!     stop

    end subroutine setup_bcondition3d
    






end module mod_bcondition3d
