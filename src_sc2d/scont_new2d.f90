subroutine scont_new2d
!
!-----------------------------------------------------------------------
!---------calculates new iterate of continuum source function-----------
!   different options: classical lambda-iteration
!                      diagonal of lambda-matrix
!                      nearest neighbours
!-----------------------------------------------------------------------
!
use prog_type
use options, only: opt_alo_cont
!
implicit none
!
! ... local scalars
integer(i4b) :: err
!
write(*,*) '--------------calculating new iterate for source function (alo)----------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
select case(opt_alo_cont)
   case(0)
      call scont_new2d_classic
   case(1)
      call scont_new2d_diag
   case(2)
!      stop 'include boundary alo in new2d_dn'
      call scont_new2d_dn
   case(3)
!     stop 'include boundary alo in new2d_nn'
      call scont_new2d_nn
   case default
      stop 'set option opt_alo_cont'
end select
!
!calculate new source function on ghost points
!
end subroutine scont_new2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new2d_classic
!
!-----------------------------------------------------------------------
!--------calculates new iterate of continuum source function------------
!-------------------classical lambda-iteration--------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime2d, only: nx, nz, imaskb2d, scont2d, bnue2d, mint2d
use params_input, only: eps_cont
!
implicit none
!
! ... local scalars
integer(i4b) :: i, k
!
! ... local functions
!
!-----------------------------------------------------------------------
!
do i=1, nx
   do k=1, nz
      select case(imaskb2d(i,k))
         case(3,4,5,6,7,8,9)
            scont2d(i,k) = (1.d0-eps_cont) * mint2d(i,k) + eps_cont*bnue2d(i,k)
         case default
      endselect
   enddo
enddo
!
!
end subroutine scont_new2d_classic
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new2d_diag
!
!-----------------------------------------------------------------------
!--------calculates new iterate of continuum source function------------
!---------approximate lambda iteration using only diagonal--------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime2d, only: nx, nz, scont2d, alocont_nn2d, mint2d, imaskb2d, bnue2d
use params_input, only: eps_cont
!
implicit none
!
! ... local scalars
integer(i4b) :: i, k
integer(i4b) :: indx_1d
real(dp) :: dummy1, dummy2, scont_new
real(dp) :: wor
!
! ... local functions
!
!----------------calculating snew directly for 2-d arrays---------------
!
dummy2=1.d0-eps_cont
!wor=0.99999d0
!alocont_nn2d(:,:,14)=alocont_nn2d(:,:,14)*wor
!
do k=1, nz
   do i=1, nx
      select case(imaskb2d(i,k))
        case(4,5,6,7,9)
            dummy1=1.d0-(1.d0-eps_cont)*alocont_nn2d(i,k,14)
            scont_new=(dummy2/dummy1) * mint2d(i,k) - &
                      (dummy2/dummy1) * alocont_nn2d(i,k,14) * scont2d(i,k) + &
                      (eps_cont/dummy1)*bnue2d(i,k)
            if(scont_new.ge.0.d0) scont2d(i,k) = scont_new
         case(3)
            !use alo on right side
            dummy1=1.d0-(1.d0-eps_cont)*alocont_nn2d(nx-1,k,14)
            scont_new=(dummy2/dummy1) * mint2d(i,k) - &
                      (dummy2/dummy1) * alocont_nn2d(nx-1,k,14) * scont2d(i,k) + &
                      (eps_cont/dummy1)*bnue2d(i,k)
            if(scont_new.ge.0.d0) scont2d(i,k) = scont_new            
         case(8)
            !use alo on left side
            dummy1=1.d0-(1.d0-eps_cont)*alocont_nn2d(2,k,14)
            scont_new=(dummy2/dummy1) * mint2d(i,k) - &
                      (dummy2/dummy1) * alocont_nn2d(2,k,14) * scont2d(i,k) + &
                      (eps_cont/dummy1)*bnue2d(i,k)
!            write(*,*) mint2d(i,k), alocont_nn2d(i,k,14), scont2d(i,k), i, k, scont_new, dummy1/dummy2, bnue2d(2,k)
            if(scont_new.ge.0.d0) scont2d(i,k) = scont_new            
         case default
      endselect
   enddo
enddo
!

!stop 'go on in scont_new2d_diag'
!
end subroutine scont_new2d_diag
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new2d_dn
!
!
!-----------------------------------------------------------------------
!------------calculates new iterate of source function------------------
!---approximate lambda iteration using nearest neighbour contribution---
!----------alo has to be stored in sparse-matrix formats----------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime2d, only: nx, nz, scont2d, mint2d, imaskb2d, bnue2d, &
                  alocont_nn2d, alocont_rowindx, alocont_colindx, &
                  alocont_data, alocont_data_diag
use params_input, only: eps_cont
use mod_math, only: conv_indx_2d_to_1d, conv_indx_1d_to_2d
use mod_sparse, only: matmul_coo, jsor_coo
!
implicit none
!
! ... local scalars
integer(i4b) :: i, k
integer(i4b) :: indx_1d, indx_x, indx_z, err, indx_rspec
real(dp) :: dummy1, dummy2, rspec
!
! ... local arrays
real(dp), dimension(:), allocatable :: dummy_vec, sum_vec, scont_vec, mint_vec, bnue_vec
logical, dimension(:), allocatable :: mask_vec
!
! ... local functions
!
!-----------------------------------------------------------------------
!---------------------allocation of local arrays------------------------
!-----------------------------------------------------------------------
!
allocate(bnue_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_dn'
allocate(scont_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_dn'
allocate(mint_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_dn'
allocate(dummy_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_dn'
allocate(mask_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_dn'
allocate(sum_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_dn'
!
!----------------transform 2d-arrays to 1-d array-----------------------
!
mint_vec=0.d0
scont_vec=0.d0
bnue_vec=0.d0
dummy_vec=0.d0
mask_vec=.false.
!
do i=1, nx
   do k=1, nz
      select case(imaskb2d(i,k))
         case(3,4,5,6,9)
         call conv_indx_2d_to_1d(i,k,nx,indx_1d)
         mask_vec(indx_1d)=.true.
         mint_vec(indx_1d)=mint2d(i,k)
         scont_vec(indx_1d)=scont2d(i,k)
         bnue_vec(indx_1d)=bnue2d(i,k)
      case default
      endselect
   enddo
enddo
!
call calc_alocont_dn2d_coo
!
!------------set up linear system that needs to be solved---------------
!-------------------in order to obtain s_new----------------------------
!
call matmul_coo(alocont_data, alocont_colindx, alocont_rowindx, scont_vec, dummy_vec, nx*nz, 5*nx*nz)
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
do i=1, 5*nx*nz
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
do i=1, 5*nx*nz
   if(alocont_rowindx(i).ne.alocont_colindx(i)) then
      sum_vec(alocont_rowindx(i)) = sum_vec(alocont_rowindx(i)) + abs(alocont_data(i))
   endif
enddo
!
rspec=maxval(sum_vec/alocont_data_diag)
indx_rspec=maxloc(sum_vec/alocont_data_diag, 1)
call conv_indx_1d_to_2d(indx_rspec, nx, indx_x, indx_z)
write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
write(*,'(a30,i10,2i5)') 'at 1d, 2d indices', indx_rspec, indx_x, indx_z
write(*,*)
!
!--------------------linear system now reads:---------------------------
!              alocont_matrix * s_new = scont_vec
!
write(*,*) '-----------------------inverting nearest neighbour alo-------------------------'
!
dummy_vec=-dummy_vec
!
!call output_alocont_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, nz, 3*nz)
!
!do i=1, 3*nz
!   write(*,*) alocont_data(i), alocont_colindx(i), alocont_rowindx(i)
!enddo
!stop 'go on in scont_new1d_dn'
!
call jsor_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, 1.d0, nx*nz, 5*nx*nz, .false., scont_vec)
!stop 'go on in scont_new1d_dn'

!
!
!---------back-transformation of source function on 3-d grid------------
!
do i=1, nx*nz
   if(mask_vec(i)) then
      call conv_indx_1d_to_2d(i, nx, indx_x, indx_z)
!      write(*,*) indx_x, indx_z, scont_vec(i)
      if(scont_vec(i).ge.0.d0) then
         scont2d(indx_x,indx_z)=scont_vec(i)
      endif
   endif
enddo
!
!update periodic boundary condition
do i=1, nz
   scont2d(nx-1,i)=scont2d(1,i)
   scont2d(nx,i)=scont2d(2,i)
enddo
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)

!write(*,*) scont2d(1,1)
!stop 'go on in scont_new'
!
!-----------------------------------------------------------------------
!
end subroutine scont_new2d_dn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new2d_nn
!
!
!-----------------------------------------------------------------------
!------------calculates new iterate of source function------------------
!---approximate lambda iteration using nearest neighbour contribution---
!----------alo has to be stored in sparse-matrix formats----------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime2d, only: nx, nz, scont2d, mint2d, imaskb2d, t2d, bnue2d, x, z, &
     alocont_nn2d, alocont_rowindx, alocont_colindx, alocont_data, alocont_data_diag
use params_input, only: eps_cont
use mod_math, only: conv_indx_2d_to_1d, conv_indx_1d_to_2d
use mod_sparse, only: matmul_coo, jsor_coo
!
implicit none
!
! ... local scalars
integer(i4b) :: i, k
integer(i4b) :: indx_1d, indx_x, indx_z, err, indx_rspec
real(dp) :: dummy1, dummy2, rspec
!
! ... local arrays
real(dp), dimension(:), allocatable :: mint_vec, scont_vec, bnue_vec, dummy_vec, sum_vec
logical, dimension(:), allocatable :: mask_vec
!
! ... local functions
!
!-----------------------------------------------------------------------
!---------------------allocation of local arrays------------------------
!-----------------------------------------------------------------------
!
allocate(bnue_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_nn: bnue_vec'
!
allocate(scont_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_nn: scont_vec'
!
allocate(mint_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_nn: mint_vec'
!
allocate(dummy_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_nn: dummy_vec'
!
allocate(mask_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_nn: mask_vec'
!
allocate(sum_vec(nx*nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new2d_nn: sum_vec'
!
!----------------transform 2d-arrays to 1-d array-----------------------
!
mint_vec=0.d0
scont_vec=0.d0
bnue_vec=0.d0
dummy_vec=0.d0
mask_vec=.false.
!
do i=1, nx
   do k=1, nz
      select case(imaskb2d(i,k))
         case(3,4,5,6,9)
            call conv_indx_2d_to_1d(i,k,nx,indx_1d)
            mint_vec(indx_1d)=mint2d(i,k)
            scont_vec(indx_1d)=scont2d(i,k)
            mask_vec(indx_1d)=.true.
            bnue_vec(indx_1d)=bnue2d(i,k)
         case default
      endselect
   enddo
enddo
!
call calc_alocont_nn2d_coo
!
!------------set up linear system that needs to be solved---------------
!-------------------in order to obtain s_new----------------------------
!
call matmul_coo(alocont_data, alocont_colindx, alocont_rowindx, scont_vec, dummy_vec, nx*nz, 9*nx*nz)
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
do i=1, 9*nx*nz
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
do i=1, 9*nx*nz
   if(alocont_rowindx(i).ne.alocont_colindx(i)) then
      sum_vec(alocont_rowindx(i)) = sum_vec(alocont_rowindx(i)) + abs(alocont_data(i))
   endif
enddo
!
rspec=maxval(sum_vec/alocont_data_diag)
indx_rspec=maxloc(sum_vec/alocont_data_diag, 1)
call conv_indx_1d_to_2d(indx_rspec, nx, indx_x, indx_z)
write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
write(*,'(a30,i10,2i5)') 'at 1d, 2d indices', indx_rspec, indx_x, indx_z
write(*,*)
!
!--------------------linear system now reads:---------------------------
!              alocont_matrix * s_new = scont_vec
!
write(*,*) '-----------------------inverting nearest neighbour alo-------------------------'
!
dummy_vec=-dummy_vec
!
!call output_alocont_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, nx*nz, 9*nx*nz)
!
call jsor_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, 1.d0, nx*nz, 9*nx*nz, .false., scont_vec)
!
!
!---------back-transformation of source function on 3-d grid------------
!
do i=1, nx*nz
   if(mask_vec(i)) then
      call conv_indx_1d_to_2d(i, nx, indx_x, indx_z)
      if(scont_vec(i).ge.0.d0) then
         scont2d(indx_x,indx_z)=scont_vec(i)
      endif
!    write(*,*) scont2d(indx_x,indx_z), imaskb2d(indx_x,indx_z)
   endif
enddo
!
!update periodic boundary condition
do i=1, nz
   scont2d(nx-1,i)=scont2d(1,i)
   scont2d(nx,i)=scont2d(2,i)
enddo
!
!stop 'go on in scont_new2d_nn'
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
end subroutine scont_new2d_nn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_alocont_dn2d_coo
!
!-----------------------------------------------------------------------
!------------calculates approximate lambda operator---------------------
!-------using only direct neighbours (8 neighbours + local point)-------
!----------------------coo storage format-------------------------------
!-----------------------------------------------------------------------
!
!         IMPORTANT NOTE: FOR PERIODIC BOUNDARY CONDITIONS
!              ALO IS ONLY CALCULATED ON LEFT BOUNDARY,
!              AND RIGHT BOUNDARY SOURCE FUNCTION IS UPDATED AT THE END
!              OTHERWISE: INVERSION BECOMES ILL-CONDITIONED (SPECTRAL RADIUS > 1)  
!
!  
use prog_type
use fund_const
use dime2d, only: nx, nz, alocont_nn2d, imaskb2d, imask2d, &
     alocont_data, alocont_data_diag, alocont_colindx, alocont_rowindx
use mod_math, only: conv_indx_2d_to_1d, conv_indx_1d_to_2d
!
implicit none
!
! ... local scalars
integer(i4b) :: i, k, err
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
!-------------------allocation of alo-matrix----------------------------
!
if(.not.allocated(alocont_data)) then
   allocate(alocont_data(5*nx*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn2d_coo: alocont_data'
endif

if(.not.allocated(alocont_data_diag)) then
   allocate(alocont_data_diag(nx*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn2d_coo: alocont_data_diag'
endif
!
if(.not.allocated(alocont_colindx)) then
   allocate(alocont_colindx(5*nx*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn2d_coo: alocont_col_indx'
endif
!
if(.not.allocated(alocont_rowindx)) then
   allocate(alocont_rowindx(5*nx*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_dn2d_coo: alocont_row_indx'
endif
!
!-----------------calculate alo in coo storage format-------------------
!
!define maximum allowed spectral radius
rspec_max=1.d0
!
alocont_data=0.d0
alocont_data_diag=0.d0
alocont_rowindx=1
alocont_colindx=1
!
iindx=1
ndiags_max = 0
ndiags_min = 5
!
do k=1, nz
   do i=1, nx

      select case(imaskb2d(i,k))
!
         case(9)
!         
!---------------------direct neighbour alo for--------------------------
!--------------------------standard points------------------------------
!
!include diagonal part
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction
             call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_row)
             alocont_data(iindx+1) = alocont_nn2d(i,k,15)
             alocont_colindx(iindx+1)=indx_1d_col
             alocont_rowindx(iindx+1)=indx_1d_row

!include direct neighbour in negative x-direction
             call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_row)
             alocont_data(iindx+2) = alocont_nn2d(i,k,13)
             alocont_colindx(iindx+2)=indx_1d_col
             alocont_rowindx(iindx+2)=indx_1d_row
!
!include direct neighbour in positive z- direction
             call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
             alocont_data(iindx+3) = alocont_nn2d(i,k,23)
             alocont_colindx(iindx+3)=indx_1d_col
             alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
             call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
             alocont_data(iindx+4) = alocont_nn2d(i,k,5)
             alocont_colindx(iindx+4)=indx_1d_col
             alocont_rowindx(iindx+4)=indx_1d_row
!
         case(3)
!         
!---------------------direct neighbour alo for--------------------------
!------------------------left left boundary-----------------------------
!
!include diagonal part
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(nx-1,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(nx-1,k,14)            
!
!include direct neighbour in positive x- direction
            call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
!
!include direct neighbour in negative x-direction (note: periodic boundary condition!!!)
            call conv_indx_2d_to_1d(nx-2,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(nx-1,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
!
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
!
         case(4)
!         
!---------------------direct neighbour alo for--------------------------
!---------------------------left boundary-------------------------------
!
!include diagonal part
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction
            call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
!
!include direct neighbour in negative x-direction (note: periodic boundary condition!!!)
            call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(i,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
!
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)!*imask2d(i,k+1)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)!*imask2d(i,k-1)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
            
         case(5)
!            
!---------------------direct neighbour alo for--------------------------
!--------------left point directly adjacent to left boundary------------
!(note: is treated analogously to 9, but for debugging and consistency reasons
!   we treat this point here explicitly

            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction
            call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
!
!include direct neighbour in negative x- direction
            call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(i,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
!            
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)!*imask2d(i,k+1)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)!*imask2d(i,k-1)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
!
         case(6)
!            
!---------------------direct neighbour alo for--------------------------
!-------------right point directly adjacent to right boundary-----------
!
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction (periodic boundary condition)
            call conv_indx_2d_to_1d(1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
!
!include direct neighbour in negative x- direction (periodic boundary condition)
            call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(i,k,13)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
!            
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)!*imask2d(i,k+1)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)!*imask2d(i,k-1)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
!
         case(1)
!            
!---------------------direct neighbour alo for--------------------------
!---------------------lower boundary condition--------------------------
!-----------------------(only in z-direction)---------------------------
!
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
         case(2)
!            
!---------------------direct neighbour alo for--------------------------
!---------------------upper boundary condition--------------------------
!-----------------------(only in z-direction)---------------------------
!
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
!
         case default
      end select
!
      iindx=iindx+5
!

!**********************old version (maybe for debugging)****************
!
!      rspec=0.d0
!      select case(imaskb2d(i,k))
!
!         case(3,4,5,6,7,8,9)
!---------------------direct neighbour alo for--------------------------
!-------------standard points directly adjacent boundary points---------
!
!include diagonal part
!            rspec=rspec+abs(alocont_nn2d(i,k,14))
!            call conv_indx_2d_to_1d(i,k, nx, indx_1d_row)
!            alocont_data(iindx) = alocont_nn2d(i,k,14)
!            alocont_colindx(iindx)=indx_1d_row
!            alocont_rowindx(iindx)=indx_1d_row
!            alocont_data_diag(indx_1d_row) = alocont_nn2d(i,k,14)
!            ndiags_ijk=1
!
!include direct neighbour in positive x- direction
! (if included, then need also to include at boundary points
!            rspec=rspec+abs(alocont_nn2d(i+1,k,13))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_col)
!               alocont_data(iindx+1) = alocont_nn2d(i+1,k,13)
!               alocont_colindx(iindx+1)=indx_1d_col
!               alocont_rowindx(iindx+1)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!
!include direct neighbour in negative x-direction
!            rspec=rspec+abs(alocont_nn2d(i-1,k,15))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_col)
!               alocont_data(iindx+2) = alocont_nn2d(i-1,k,15)
!               alocont_colindx(iindx+2)=indx_1d_col
!               alocont_rowindx(iindx+2)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!
!include direct neighbour in positive z-direction
!            rspec=rspec+abs(alocont_nn2d(i,k+1,5))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_col)
!               alocont_data(iindx+3) = alocont_nn2d(i,k+1,5)
!               alocont_colindx(iindx+3)=indx_1d_col
!               alocont_rowindx(iindx+3)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!            write(*,*) indx_1d_col, indx_1d_row, alocont_data(iindx+3)
!
!include direct neighbour in negative z-direction
!            rspec=rspec+abs(alocont_nn2d(i,k-1,23))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_col)
!               alocont_data(iindx+4) = alocont_nn2d(i,k-1,23)
!               alocont_colindx(iindx+4)=indx_1d_col
!               alocont_rowindx(iindx+4)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!
!         case(10)
!
!------------alo procedure on left boundary (point bl)------------------
!
!!include diagonal part
!            rspec=rspec+abs(alocont_nn2d(i,k,14))
!            call conv_indx_2d_to_1d(i,k, nx, indx_1d_row)
!            alocont_data(iindx) = alocont_nn2d(i,k,14)
!            alocont_colindx(iindx)=indx_1d_row
!            alocont_rowindx(iindx)=indx_1d_row
!            alocont_data_diag(indx_1d_row) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction
!            rspec=rspec+abs(alocont_nn2d(i+1,k,13))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_col)
!               alocont_data(iindx+1) = alocont_nn2d(i+1,k,13)
!               alocont_colindx(iindx+1)=indx_1d_col
!               alocont_rowindx(iindx+1)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!            write(*,*) indx_1d_col, indx_1d_row, alocont_data(iindx+1)
!
!include direct neighbour in negative x-direction (from right boundary)
!            rspec=rspec+abs(alocont_nn2d(nx-2,k,15))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(nx-2,k,nx,indx_1d_col)
!               alocont_data(iindx+2) = alocont_nn2d(nx-2,k,15)
!               alocont_colindx(iindx+2)=indx_1d_col
!               alocont_rowindx(iindx+2)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!              ndiags_min = min(ndiags_ijk,ndiags_min)
!              iindx=iindx+5
!              cycle
!           endif
!            write(*,*) indx_1d_col, indx_1d_row, alocont_data(iindx+2), alocont_nn2d(nx-2,k,15)
!
!include direct neighbour in positive z-direction
!            rspec=rspec+abs(alocont_nn2d(i,k+1,5))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_col)
!               alocont_data(iindx+3) = alocont_nn2d(i,k+1,5)
!               alocont_colindx(iindx+3)=indx_1d_col
!               alocont_rowindx(iindx+3)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!            write(*,*) indx_1d_col, indx_1d_row, alocont_data(iindx+2)
!
!include direct neighbour in negative z-direction
!            rspec=rspec+abs(alocont_nn2d(i,k-1,23))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_col)
!               alocont_data(iindx+4) = alocont_nn2d(i,k-1,23)
!               alocont_colindx(iindx+4)=indx_1d_col
!               alocont_rowindx(iindx+4)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!
!
!         case(11)
!
!----------------alo procedure on right boundary-------------------------
!
!include diagonal part
!            rspec=rspec+abs(alocont_nn2d(i,k,14))
!            call conv_indx_2d_to_1d(i,k, nx, indx_1d_row)
!            alocont_data(iindx) = alocont_nn2d(i,k,14)
!            alocont_colindx(iindx)=indx_1d_row
!            alocont_rowindx(iindx)=indx_1d_row
!            alocont_data_diag(indx_1d_row) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction (from left boundary)
!            rspec=rspec+abs(alocont_nn2d(3,k,13))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(3,k,nx,indx_1d_col)
!               alocont_data(iindx+1) = alocont_nn2d(3,k,13)
!               alocont_colindx(iindx+1)=indx_1d_col
!               alocont_rowindx(iindx+1)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!
!include direct neighbour in negative x-direction
!            rspec=rspec+abs(alocont_nn2d(i-1,k,15))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_col)
!               alocont_data(iindx+2) = alocont_nn2d(i-1,k,15)
!               alocont_colindx(iindx+2)=indx_1d_col
!               alocont_rowindx(iindx+2)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!            write(*,*) indx_1d_col, indx_1d_row, alocont_data(iindx+2)
!
!include direct neighbour in positive z-direction
!            rspec=rspec+abs(alocont_nn2d(i,k+1,5))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_col)
!               alocont_data(iindx+3) = alocont_nn2d(i,k+1,5)
!               alocont_colindx(iindx+3)=indx_1d_col
!               alocont_rowindx(iindx+3)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!
!include direct neighbour in negative z-direction
!            rspec=rspec+abs(alocont_nn2d(i,k-1,23))
!            if(rspec.lt.rspec_max) then
!               call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_col)
!               alocont_data(iindx+4) = alocont_nn2d(i,k-1,23)
!               alocont_colindx(iindx+4)=indx_1d_col
!               alocont_rowindx(iindx+4)=indx_1d_row
!               ndiags_ijk=ndiags_ijk+1
!               ndiags_max = max(ndiags_ijk,ndiags_max)
!            else
!               ndiags_min = min(ndiags_ijk,ndiags_min)
!               iindx=iindx+5
!               cycle
!            endif
!         case default
!      end select
!!
!
!      write(*,*) rspec
!      iindx=iindx+5
!
!******************************************************************      
!      
   enddo
enddo
!!
!!
!write(*,*)
!write(*,*) 'maximum number of used neighbours for ALO calculation', ndiags_max-1
!write(*,*) 'minimum number of used neighbours for ALO calculation', ndiags_min-1
!write(*,*)

!open(1, file='TRASH/alo_sc2d.dat')
!do i=1, 5*nx*nz
!   write(1,'(2i10,es20.8)'), alocont_rowindx(i), alocont_colindx(i), alocont_data(i)
!   write(*,'(2i10,es20.8)'), alocont_rowindx(i), alocont_colindx(i), alocont_data(i)
!enddo
!
!stop 'go on in calc_alocont_dn2d_coo'

end subroutine calc_alocont_dn2d_coo

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_alocont_nn2d_coo
!
!-----------------------------------------------------------------------
!------------calculates approximate lambda operator---------------------
!-------using only nearest neighbours (26 neighbours + local point)-----
!----------------------coo storage format-------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime2d, only: nx, nz, alocont_nn2d, &
                  alocont_data, alocont_data_diag, alocont_colindx, alocont_rowindx, &
                  imaskb2d, imask2d, x, z
use mod_math, only: conv_indx_2d_to_1d, conv_indx_1d_to_2d
!
implicit none
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
!-------------------allocation of alo-matrix----------------------------
!
if(.not.allocated(alocont_data)) then
   allocate(alocont_data(9*nx*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn2d_coo: alocont_data'
endif

if(.not.allocated(alocont_data_diag)) then
   allocate(alocont_data_diag(nx*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn2d_coo: alocont_data_diag'
endif
!
if(.not.allocated(alocont_colindx)) then
   allocate(alocont_colindx(9*nx*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn2_coo: alocont_col_indx'
endif
!
if(.not.allocated(alocont_rowindx)) then
   allocate(alocont_rowindx(9*nx*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont_nn2d_coo: alocont_row_indx'
endif
!
!-----------------calculate alo in coo storage format-------------------
!
!define maximum allowed spectral radius
rspec_max=1.d0
!
alocont_data=0.d0
alocont_data_diag=0.d0
alocont_rowindx=1
alocont_colindx=1
!
iindx=1
ndiags_max = 0
ndiags_min = 9
!
do i=1, nx
   do k=1, nz


      select case(imaskb2d(i,k))
!
         case(9)
!         
!---------------------nearest neighbour alo for-------------------------
!--------------------------standard points------------------------------
!
!include diagonal part
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction
            call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)!*imask2d(i+1,k)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row

!include direct neighbour in negative x-direction
            call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(i,k,13)!*imask2d(i-1,k)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
!
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)!*imask2d(i,k+1)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)!*imask2d(i,k-1)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row

!include direct neighbour in positive x-direction and positive z-direction
            call conv_indx_2d_to_1d(i+1,k+1,nx,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn2d(i,k,24)!*imask2d(i+1,k+1)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
!
!include direct neighbour in positive x-direction and negative z-direction
            call conv_indx_2d_to_1d(i+1,k-1,nx,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn2d(i,k,6)!*imask2d(i+1,k-1)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
!
!include direct neighbour in negative x-direction and positive z-direction
            call conv_indx_2d_to_1d(i-1,k+1,nx,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn2d(i,k,22)!*imask2d(i-1,k+1)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
!
!include direct neighbour in negative x-direction and negative z-direction
            call conv_indx_2d_to_1d(i-1,k-1,nx,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn2d(i,k,4)!*imask2d(i-1,k-1)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
!
         case(3)
!         
!--------------------nearest neighbour alo for--------------------------
!------------------------left left boundary-----------------------------
!
!include diagonal part
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction
            call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)!*imask2d(i+1,k)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
!
!include direct neighbour in negative x-direction (note: periodic boundary condition!!!)
            call conv_indx_2d_to_1d(nx-2,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(nx-1,k,13)!*imask2d(nx-2,k)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
!
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)!*imask2d(i,k+1)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)!*imask2d(i,k-1)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
!
!include direct neighbour in positive x-direction and positive z-direction
            call conv_indx_2d_to_1d(i+1,k+1,nx,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn2d(i,k,24)!*imask2d(i+1,k+1)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
!
!include direct neighbour in positive x-direction and negative z-direction
            call conv_indx_2d_to_1d(i+1,k-1,nx,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn2d(nx-1,k,6)!*imask2d(i+1,k-1)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
!
!include direct neighbour in negative x-direction and positive z-direction
!(in case of quadratic interpolations at boundary points, 22 can case problems)            
!           call conv_indx_2d_to_1d(nx-2,k+1,nx,indx_1d_row)
!           alocont_data(iindx+7) = alocont_nn2d(nx-1,k,22)!*imask2d(nx-2,k+1)
!           alocont_colindx(iindx+7)=indx_1d_col
!           alocont_rowindx(iindx+7)=indx_1d_row
!
!include direct neighbour in negative x-direction and negative z-direction
!            call conv_indx_2d_to_1d(nx-2,k-1,nx,indx_1d_row)
!            alocont_data(iindx+8) = alocont_nn2d(nx-1,k,4)!*imask2d(nx-2,k-1)
!            alocont_colindx(iindx+8)=indx_1d_col
!            alocont_rowindx(iindx+8)=indx_1d_row
!
         case(4)
!         
!--------------------nearest neighbour alo for--------------------------
!---------------------------left boundary-------------------------------
!
!include diagonal part
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction
            call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)!*imask2d(i+1,k)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
!
!include direct neighbour in negative x-direction (note: periodic boundary condition!!!)
            call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(i,k,13)!*imask2d(i-1,k)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row

!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)!*imask2d(i,k+1)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)!*imask2d(i,k-1)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row

!include direct neighbour in positive x-direction and positive z-direction
            call conv_indx_2d_to_1d(i+1,k+1,nx,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn2d(i,k,24)!*imask2d(i+1,k+1)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
!
!include direct neighbour in positive x-direction and negative z-direction
            call conv_indx_2d_to_1d(i+1,k-1,nx,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn2d(i,k,6)!*imask2d(i+1,k-1)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
!
!include direct neighbour in negative x-direction and positive z-direction
            call conv_indx_2d_to_1d(i-1,k+1,nx,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn2d(i,k,22)!*imask2d(i-1,k+1)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
!
!include direct neighbour in negative x-direction and negative z-direction
            call conv_indx_2d_to_1d(i-1,k-1,nx,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn2d(i,k,4)!*imask2d(i-1,k-1)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
!            
         case(5)
!            
!--------------------nearest neighbour alo for--------------------------
!--------------left point directly adjacent to left boundary------------
!(note: is treated analogously to 9, but for debugging and consistency reasons
!   we treat this point here explicitly

            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction
            call conv_indx_2d_to_1d(i+1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)!*imask2d(i+1,k)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
!
!include direct neighbour in negative x- direction
            call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(i,k,13)!*imask2d(i-1,k)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
!            
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)!*imask2d(i,k+1)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)!*imask2d(i,k-1)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row

!include direct neighbour in positive x-direction and positive z-direction
            call conv_indx_2d_to_1d(i+1,k+1,nx,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn2d(i,k,24)!*imask2d(i+1,k+1)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
!
!include direct neighbour in positive x-direction and negative z-direction
            call conv_indx_2d_to_1d(i+1,k-1,nx,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn2d(i,k,6)!*imask2d(i+1,k-1)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
!
!include direct neighbour in negative x-direction and positive z-direction
            call conv_indx_2d_to_1d(i-1,k+1,nx,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn2d(i,k,22)!*imask2d(i-1,k+1)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
!
!include direct neighbour in negative x-direction and negative z-direction
            call conv_indx_2d_to_1d(i-1,k-1,nx,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn2d(i,k,4)!*imask2d(i-1,k-1)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
!
         case(6)
!            
!--------------------nearest neighbour alo for--------------------------
!-------------right point directly adjacent to right boundary-----------
!
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
            alocont_data(iindx) = alocont_nn2d(i,k,14)
            alocont_colindx(iindx)=indx_1d_col
            alocont_rowindx(iindx)=indx_1d_col
            alocont_data_diag(indx_1d_col) = alocont_nn2d(i,k,14)
!
!include direct neighbour in positive x- direction (periodic boundary condition)
            call conv_indx_2d_to_1d(1,k,nx,indx_1d_row)
            alocont_data(iindx+1) = alocont_nn2d(i,k,15)!*imask2d(1,k)
            alocont_colindx(iindx+1)=indx_1d_col
            alocont_rowindx(iindx+1)=indx_1d_row
!
!include direct neighbour in negative x- direction (periodic boundary condition)
            call conv_indx_2d_to_1d(i-1,k,nx,indx_1d_row)
            alocont_data(iindx+2) = alocont_nn2d(i,k,13)!*imask2d(i-1,k)
            alocont_colindx(iindx+2)=indx_1d_col
            alocont_rowindx(iindx+2)=indx_1d_row
!            
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)!*imask2d(i,k+1)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)!*imask2d(i,k-1)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row

!include direct neighbour in positive x-direction and positive z-direction
            call conv_indx_2d_to_1d(1,k+1,nx,indx_1d_row)
            alocont_data(iindx+5) = alocont_nn2d(1,k,24)!*imask2d(1,k+1)
            alocont_colindx(iindx+5)=indx_1d_col
            alocont_rowindx(iindx+5)=indx_1d_row
!
!include direct neighbour in positive x-direction and negative z-direction
            call conv_indx_2d_to_1d(1,k-1,nx,indx_1d_row)
            alocont_data(iindx+6) = alocont_nn2d(1,k,6)!*imask2d(1,k-1)
            alocont_colindx(iindx+6)=indx_1d_col
            alocont_rowindx(iindx+6)=indx_1d_row
!
!include direct neighbour in negative x-direction and positive z-direction
            call conv_indx_2d_to_1d(i-1,k+1,nx,indx_1d_row)
            alocont_data(iindx+7) = alocont_nn2d(i,k,22)!*imask2d(i-1,k+1)
            alocont_colindx(iindx+7)=indx_1d_col
            alocont_rowindx(iindx+7)=indx_1d_row
!
!include direct neighbour in negative x-direction and negative z-direction
            call conv_indx_2d_to_1d(i-1,k-1,nx,indx_1d_row)
            alocont_data(iindx+8) = alocont_nn2d(i,k,4)!*imask2d(i-1,k-1)
            alocont_colindx(iindx+8)=indx_1d_col
            alocont_rowindx(iindx+8)=indx_1d_row
!
         case(1)
!            
!--------------------nearest neighbour alo for--------------------------
!---------------------lower boundary condition--------------------------
!-----------------------(only in z-direction)---------------------------
!
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
!include direct neighbour in positive z- direction
            call conv_indx_2d_to_1d(i,k+1,nx,indx_1d_row)
            alocont_data(iindx+3) = alocont_nn2d(i,k,23)
            alocont_colindx(iindx+3)=indx_1d_col
            alocont_rowindx(iindx+3)=indx_1d_row
!include direct neighbour in positive x-direction and positive z-direction
!            call conv_indx_2d_to_1d(i+1,k+1,nx,indx_1d_row)
!            alocont_data(iindx+3) = alocont_nn2d(i,k,24)
!            alocont_colindx(iindx+3)=indx_1d_col
!            alocont_rowindx(iindx+3)=indx_1d_row
!
         case(2)
!            
!--------------------nearest neighbour alo for--------------------------
!---------------------upper boundary condition--------------------------
!-----------------------(only in z-direction)---------------------------
!
            call conv_indx_2d_to_1d(i,k, nx, indx_1d_col)
!include direct neighbour in negative z-direction
            call conv_indx_2d_to_1d(i,k-1,nx,indx_1d_row)
            alocont_data(iindx+4) = alocont_nn2d(i,k,5)
            alocont_colindx(iindx+4)=indx_1d_col
            alocont_rowindx(iindx+4)=indx_1d_row
!include direct neighbour in positive x-direction and negative z-direction
!            call conv_indx_2d_to_1d(i+1,k-1,nx,indx_1d_row)
!            alocont_data(iindx+3) = alocont_nn2d(i,k,6)
!            alocont_colindx(iindx+3)=indx_1d_col
!            alocont_rowindx(iindx+3)=indx_1d_row
!
         case default
      end select
!
!-----------------------------------------------------------------------
!
      iindx=iindx+9
!
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
!   write(*,*) alocont_nn2d(2,k,22), alocont_nn2d(3,k,22), alocont_nn2d(4,k,22), alocont_nn2d(5,k,22), alocont_nn2d(nx-2,k,22), alocont_nn2d(nx-1,k,22), alocont_nn2d(nx,k,22)
!enddo
!
!stop 'go on in calc_alocont'


end subroutine calc_alocont_nn2d_coo
