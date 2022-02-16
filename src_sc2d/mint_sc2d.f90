subroutine mint_sc2d(opt_method)
!
!-----------------------------------------------------------------------
!-----------calculates mean intensities at all grid point---------------
!-----------------------------------------------------------------------
!
!   formal solution from 2d short characteristics
!
!                  mu-integration between [-1,1]
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime2d, only: int2d, alocont_o_nn2d, alocont_nn2d, mint2d, normalization2d, &
                  alocont_nn2d_tmp, mint2d_tmp, normalization2d_tmp, bnue2d, t2d, &
                  nx, nz, imask2d, scont2d
use angles, only: dim_omega
use freq, only: nodes_nue
!
implicit none
!
! ... local scalars
integer(i4b) :: i, k
integer(i4b) :: nuindx
integer(i4b) :: oindx
integer(i4b) :: err
integer(i4b), intent(in) :: opt_method
!
! ... local arrays
!
! ... local functions
real(dp) :: bnue
!
! ... local characters
!
alocont_nn2d=0.d0
normalization2d=0.d0
mint2d=0.d0
!
!--------------deallocation of global (threadprivate) arrays------------
!
if(allocated(int2d)) deallocate(int2d)
if(allocated(alocont_o_nn2d)) deallocate(alocont_o_nn2d)
!
!
!-----------------------begin of parallel region------------------------
!
!$omp parallel &
!$omp private(err, oindx)
!
!-------------allocation of global (threadprivate) arrays---------------
!
   allocate(alocont_o_nn2d(nx,nz,27), stat=err)
      if(err.ne.0) stop 'allocation error in mint_sc2d: alocont_o_nn2d'
!
   allocate(alocont_nn2d_tmp(nx,nz,27), stat=err)
      if(err.ne.0) stop 'allocation error in mint_sc2d: alocont_nn2d_tmp'
!
   allocate(normalization2d_tmp(nx,nz), stat=err)
      if(err.ne.0) stop 'allocation error mint_sc2d: normalization2d_tmp'
!
   allocate(mint2d_tmp(nx,nz), stat=err)
      if(err.ne.0) stop 'allocation error mint_sc2d: mint2d_tmp'
!
   allocate(int2d(nx,nz), stat=err)
      if(err.ne.0) stop 'allocation error mint_sc2d:: int2d'
!
   alocont_nn2d_tmp=0.d0
   normalization2d_tmp=0.d0
   mint2d_tmp=0.d0
!
   !$omp do schedule(dynamic)
      do oindx=dim_omega, 1, -1
!         write(*,'(a30,2i10,a5,i5)') 'calculating omega', oindx, dim_omega, 'bez'
         if(opt_method.eq.4) then
            call fsc_cont2d_lin(oindx,1)
         elseif(opt_method.eq.5) then
            call fsc_cont2d(oindx,1)
         endif
      enddo
   !$omp enddo
!
!------------------------add up temporary arrays------------------------
!
   !$omp critical
      alocont_nn2d = alocont_nn2d+alocont_nn2d_tmp
      normalization2d = normalization2d + normalization2d_tmp
      mint2d = mint2d + mint2d_tmp
   !$omp end critical
!
!------------deallocation of global (threadprivate) arrays--------------
!
   deallocate(alocont_nn2d_tmp)
   deallocate(normalization2d_tmp)
   deallocate(mint2d_tmp)
   deallocate(alocont_o_nn2d)
   deallocate(int2d)
!
!$omp end parallel
!
!renormalize (note that mean intensities also calculated at ghost points on left and right boundary   
do i=1, nx
   do k=1, nz
!do i=3, nx-2
!   do k=3, nz-2
      select case(imask2d(i,k))
         case(1,2)
!            if(abs(normalization2d(i,k)-0.d0).lt.1.d-6) then
            if(abs(normalization2d(i,k)-one).gt.1.d-15) then
            write(*,'(a50, 2i5, es20.8)') 'normalization error in mint_sc2d:', i, k, normalization2d(i,k)
               stop
            endif
            mint2d(i,k) = mint2d(i,k)/normalization2d(i,k)
            alocont_nn2d(i,k,:) = alocont_nn2d(i,k,:)/normalization2d(i,k)
            normalization2d(i,k) = normalization2d(i,k)/normalization2d(i,k)
         case default
      end select
   enddo
enddo
!
!do k=1, nz
!   write(*,*) k, scont2d(nx,k), normalization2d(nx,k), mint2d(nx,k)
!enddo
!write(*,*) mint2d(nx-1,4), alocont_nn2d(nx-2,4,15)
!write(*,*) mint2d(nx-1,5), alocont_nn2d(nx-2,5,15)
!write(*,*) mint2d(nx-1,6), alocont_nn2d(nx-2,6,15)
!write(*,*) mint2d(nx-1,7), alocont_nn2d(nx-2,7,15)
!stop
!
end subroutine mint_sc2d
