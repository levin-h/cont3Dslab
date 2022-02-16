subroutine scont_new1d
!
!-----------------------------------------------------------------------
!---------calculates new iterate of continuum source function-----------
!   different options: classical lambda-iteration
!                      diagonal of lambda-matrix
!                      nearest neighbours (= direct neighbour in 1D)
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
      call scont_new1d_classic
   case(1)
      call scont_new1d_diag
   case(2)
      call scont_new1d_dn
   case(3)
      call scont_new1d_dn
   case default
      stop 'set option opt_alo_cont'
end select
!
!
end subroutine scont_new1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new1d_classic
!
!-----------------------------------------------------------------------
!--------calculates new iterate of continuum source function------------
!-------------------classical lambda-iteration--------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime1d, only: nz, imask1d, scont1d, bnue1d, mint1d
use params_input, only: eps_cont
!
implicit none
!
! ... local scalars
integer(i4b) :: i
!
! ... local functions
!
!-----------------------------------------------------------------------
!
do i=1, nz
   select case(imask1d(i))
      case(1)
         scont1d(i) = (1.d0-eps_cont) * mint1d(i) + eps_cont*bnue1d(i)
      case default
   endselect
enddo
!
end subroutine scont_new1d_classic
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new1d_diag
!
!-----------------------------------------------------------------------
!--------calculates new iterate of continuum source function------------
!---------approximate lambda iteration using only diagonal--------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime1d, only: nz, scont1d, alocont_nn1d, mint1d, imask1d, bnue1d
use params_input, only: eps_cont
!
implicit none
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: indx_1d
real(dp) :: dummy1, dummy2, scont_new
!
! ... local functions
!
!----------------calculating snew directly for 1-d arrays---------------
!
dummy2=1.d0-eps_cont
!
do i=1,nz
   select case(imask1d(i))
      case(1)
         dummy1=1.d0-(1.d0-eps_cont)*alocont_nn1d(i,2)
         scont_new =  (dummy2/dummy1) * mint1d(i) - &
                      (dummy2/dummy1) * alocont_nn1d(i,2) * scont1d(i) + &
                      (eps_cont/dummy1)*bnue1d(i)
!         write(*,*) i, alocont_nn1d(i,2), scont_new
         if(scont_new.ge.0.d0) then
            scont1d(i) = scont_new
         endif
      case default
   endselect
enddo
!
!stop 'go on in scont_new1d_diag'
!
!
end subroutine scont_new1d_diag
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine scont_new1d_dn
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
use dime1d, only: nz, scont1d, mint1d, imask1d, bnue1d, &
                  alocont_nn1d, alocont_rowindx, alocont_colindx, alocont_data, alocont_data_diag
use params_input, only: eps_cont
use mod_sparse, only: matmul_coo, jsor_coo
!
implicit none
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: indx_1d, indx_x, indx_y, indx_z, err, indx_rspec
real(dp) :: dummy1, dummy2, rspec
!
! ... local arrays
real(dp), dimension(:), allocatable :: dummy_vec, sum_vec, scont_vec, mint_vec, bnue_vec
logical, dimension(:), allocatable :: mask_vec
!
! ... local functions
real(dp) :: bnue
!
!-----------------------------------------------------------------------
!---------------------allocation of local arrays------------------------
!-----------------------------------------------------------------------
!
allocate(scont_vec(nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new1d_dn'
allocate(bnue_vec(nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new1d_dn'
allocate(mint_vec(nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new1d_dn'
allocate(dummy_vec(nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new1d_dn: dummy_vec'
allocate(mask_vec(nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new1d_dn: mask_vec'
allocate(sum_vec(nz), stat=err)
   if(err.ne.0) stop 'allocation error scont_new1d_dn: sum_vec'
!
!----------------transform 3d-arrays to 1-d array-----------------------
!   
scont_vec=scont1d
mint_vec=mint1d
bnue_vec=bnue1d
dummy_vec=0.d0
mask_vec=.false.
!
do i=1, nz
   select case(imask1d(i))
      case(1)
         mask_vec(i)=.true.
      case default
   endselect
enddo
!
call calc_alocont1d_dn_coo

!do i=1, 3*nz
!   write(*,*) alocont_data(i), alocont_colindx(i), alocont_rowindx(i)
!enddo
!stop 'go on in scont_new1d_dn'

!
!------------set up linear system that needs to be solved---------------
!-------------------in order to obtain s_new----------------------------
!
call matmul_coo(alocont_data, alocont_colindx, alocont_rowindx, scont_vec, dummy_vec, nz, 3*nz)
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
do i=1, 3*nz
   if(alocont_colindx(i).eq.alocont_rowindx(i)) then
      alocont_data(i)=1.d0-alocont_data(i)
   else
      alocont_data(i)=-1.d0*alocont_data(i)
   endif
enddo
!
!write(*,*)
!do i=1, 3*nz
!   write(*,*) alocont_data(i), alocont_colindx(i), alocont_rowindx(i)
!enddo
!stop 'go on in scont_new1d_dn'
!
!--------------check if matrix is diagonal dominant---------------------
!-----------------and estimate spectral radius--------------------------
!
sum_vec=0.d0
!
do i=1, 3*nz
   if(alocont_rowindx(i).ne.alocont_colindx(i)) then
      sum_vec(alocont_rowindx(i)) = sum_vec(alocont_rowindx(i)) + abs(alocont_data(i))
   endif
enddo
!
rspec=maxval(sum_vec/alocont_data_diag)
indx_rspec=maxloc(sum_vec/alocont_data_diag, 1)
write(*,'(a30,es20.8)') 'estimate of spectral radius', rspec
write(*,'(a30,i10,3i5)') 'at 1d index', indx_rspec
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

call jsor_coo(alocont_data, alocont_colindx, alocont_rowindx, alocont_data_diag, dummy_vec, 1.d0, nz, 3*nz, .false., scont_vec)
!stop 'go on in scont_new1d_dn'

!
!
!---------back-transformation of source function on 3-d grid------------
!
do i=1, nz
   if(mask_vec(i)) then
      if(scont_vec(i).ge.0.d0) then
         scont1d(i)=scont_vec(i)
      endif
   endif
enddo
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)
!
!-----------------------------------------------------------------------
!
end subroutine scont_new1d_dn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_alocont1d_dn_coo
!
!-----------------------------------------------------------------------
!------------calculates approximate lambda operator---------------------
!-------using only direct neighbours (2 neighbours + local point)-------
!----------------------coo storage format-------------------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime1d, only: nz, alocont_nn1d, &
                  alocont_data, alocont_data_diag, alocont_colindx, alocont_rowindx
!
implicit none
!
! ... local scalars
integer(i4b) :: i, err
integer(i4b) :: iindx, indx_1d, indx_1d_col, indx_x, indx_y, indx_z
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
   allocate(alocont_data(3*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont1d_dn_coo: alocont_data'
endif

if(.not.allocated(alocont_data_diag)) then
   allocate(alocont_data_diag(nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont1d_dn_coo: alocont_data_diag'
endif
!
if(.not.allocated(alocont_colindx)) then
   allocate(alocont_colindx(3*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont1d_dn_coo: alocont_col_indx'
endif
!
if(.not.allocated(alocont_rowindx)) then
   allocate(alocont_rowindx(3*nz), stat=err)
      if(err.ne.0) stop 'allocation error calc_alocont1d_dn_coo: alocont_row_indx'
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
!
do i=1, nz
!
!   write(*,*) alocont_nn1d(i,1), alocont_nn1d(i,2), alocont_nn1d(i,3)
   
   if(i.eq.1) then
      alocont_data_diag(i) = 1.d0      
      alocont_data(iindx) = 1.d0
      alocont_rowindx(iindx) = i
      alocont_colindx(iindx) = i
      iindx=iindx+1
      alocont_data(iindx) = 0.d0
      alocont_rowindx(iindx) = i
      alocont_colindx(iindx) = i+1
      iindx=iindx+1
      alocont_data(iindx) = 0.d0
      alocont_rowindx(iindx) = i
      alocont_colindx(iindx) = i+2
      iindx=iindx+1
   elseif(i.eq.nz) then 
      alocont_data_diag(i) = 1.d0      
      alocont_data(iindx) = 1.d0
      alocont_rowindx(iindx) = i
      alocont_colindx(iindx) = i
      iindx=iindx+1
      alocont_data(iindx) = 0.d0
      alocont_rowindx(iindx) = i
      alocont_colindx(iindx) = i-1
      iindx=iindx+1
      alocont_data(iindx) = 0.d0
      alocont_rowindx(iindx) = i
      alocont_colindx(iindx) = i-2
      iindx=iindx+1
   else
      
   rspec=0.d0
!
!include diagonal part
   rspec=rspec+abs(alocont_nn1d(i,2))
   indx_1d_col=i
   alocont_data(iindx) = alocont_nn1d(i,2)
   alocont_colindx(iindx)=indx_1d_col
   alocont_rowindx(iindx)=indx_1d_col
   alocont_data_diag(indx_1d_col) = alocont_nn1d(i,2)
!
!include direct neighbour in positive z-direction
   rspec=rspec+abs(alocont_nn1d(i+1,1))
   indx_1d=i+1
   alocont_data(iindx+1) = alocont_nn1d(i+1,1)
   alocont_colindx(iindx+1)=indx_1d
   alocont_rowindx(iindx+1)=indx_1d_col
!
!include direct neighbour in negative z-direction
   rspec=rspec+abs(alocont_nn1d(i-1,3))
   indx_1d=i-1
   alocont_data(iindx+2) = alocont_nn1d(i-1,3)
   alocont_colindx(iindx+2)=indx_1d
   alocont_rowindx(iindx+2)=indx_1d_col
!
!-----------------------------------------------------------------------
!
   iindx=iindx+3
!
   endif
!
enddo
!


!stop 'go on in calc_alocont1d_dn_coo'
!
end subroutine calc_alocont1d_dn_coo
