subroutine mint_sc1d
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
use prog_type
use fund_const
use dime1d, only: int1d, alocont_o_nn1d, alocont_nn1d, mint1d, normalization1d, &
                  alocont_nn1d_tmp, mint1d_tmp, normalization1d_tmp, bnue1d, t1d, &
                  nz, imask1d
use angles, only: dim_mu
use freq, only: nodes_nue
use options, only: opt_method
!
implicit none
!
! ... local scalars
integer(i4b) :: i
integer(i4b) :: nuindx
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
alocont_nn1d=0.d0
normalization1d=0.d0
mint1d=0.d0
!
!--------------deallocation of global (threadprivate) arrays------------
!
if(allocated(int1d)) deallocate(int1d)
if(allocated(alocont_o_nn1d)) deallocate(alocont_o_nn1d)
!
!
!-----------------------begin of parallel region------------------------
!
!$omp parallel &
!$omp private(err, oindx)
!
!-------------allocation of global (threadprivate) arrays---------------
!
   allocate(alocont_o_nn1d(nz,3), stat=err)
      if(err.ne.0) stop 'allocation error in mint_sc1d: alocont_o_nn1d'
!
   allocate(alocont_nn1d_tmp(nz,3), stat=err)
      if(err.ne.0) stop 'allocation error in mint_sc1d: alocont_nn1d_tmp'
!
   allocate(normalization1d_tmp(nz), stat=err)
      if(err.ne.0) stop 'allocation error mint_sc1d: normalization1d_tmp'
!
   allocate(mint1d_tmp(nz), stat=err)
      if(err.ne.0) stop 'allocation error mint_sc1d: mint1d_tmp'
!
   allocate(int1d(nz), stat=err)
      if(err.ne.0) stop 'allocation error ffvm_line1d: int1d'
!
   alocont_nn1d_tmp=0.d0
   normalization1d_tmp=0.d0
   mint1d_tmp=0.d0
!
   !$omp do schedule(dynamic)
      do oindx=1, dim_mu
!         write(*,'(a30,2i10,a5,i5)') 'calculating omega', oindx, dim_omega, 'bez'
         if(opt_method.eq.4) then
            call fsc_cont1d_lin(oindx,1)
         elseif(opt_method.eq.5) then
            call fsc_cont1d(oindx,1)
         endif
      enddo
   !$omp enddo
!   stop 'go on in mint_sc1d'
!
!------------------------add up temporary arrays------------------------
!
   !$omp critical
      alocont_nn1d = alocont_nn1d+alocont_nn1d_tmp
      normalization1d = normalization1d + normalization1d_tmp
      mint1d = mint1d + mint1d_tmp
   !$omp end critical
!
!------------deallocation of global (threadprivate) arrays--------------
!
   deallocate(alocont_nn1d_tmp)
   deallocate(normalization1d_tmp)
   deallocate(mint1d_tmp)
   deallocate(alocont_o_nn1d)
   deallocate(int1d)
!
!$omp end parallel
!
!renormalize
do i=1, nz
   select case(imask1d(i))
      case(1)
         if(abs(normalization1d(i)-0.d0).lt.1.d-6) then
            write(*,'(a50, i5, es20.8)') 'normalization error in mint_sc1d:', i, normalization1d(i)
            stop
         endif
         mint1d(i) = mint1d(i)/normalization1d(i)
         alocont_nn1d(i,:) = alocont_nn1d(i,:)/normalization1d(i)
         normalization1d(i) = normalization1d(i)/normalization1d(i)
      case default
   end select
enddo
!
!
end subroutine mint_sc1d
