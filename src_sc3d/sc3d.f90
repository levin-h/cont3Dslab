!
!continuum transport for slab with periodic boundary conditions
!                                  at left,right,back and front boundary
! vNEW: including modularized interpolation routines
!       including modularized routines for most other stuff as well
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
program prog
!
use prog_type
use fund_const
use omp_lib
use mod_timing
use options, only: opt_method, opt_opac, opt_epsc
use mod_grid3d, only: check_grid3d, gridxyz, allocate_global3d
use mod_angles, only: calcnodes_omega, check_nodes_omega
use mod_frequencies, only: calcnodes_nue
use mod_model3d, only: setup_mod3d, setup_opacities1, setup_opacities2, setup_bnue, setup_epsc1, setup_epsc2
use mod_io, only: read_input, print_model, print_solution, output
use mod_conttrans3d, only: conttrans_sc3d, conttrans_lte3d
use params_input, only: verbose
use mod_benchmark3d, only: make_benchmark
use mod_bcondition3d, only: setup_bcondition3d
!
implicit none
!
! ... local scalars
integer(i4b) :: nueindx
!
! ... local arrays
!
! ... local characters
!
! ... local logicals
!
! ... local functions
!
!-----------------------------------------------------------------------
!
ts_tot = omp_get_wtime()
!
!read input
call read_input(.true.)
!
!create x,y,z-grid
call gridxyz(verbose=verbose)
!
!create model
call setup_mod3d(1, verbose=verbose)
!
!
!create angular grid
call calcnodes_omega(verbose=verbose)
call check_nodes_omega(verbose=verbose)
!
!create frequency grid
call calcnodes_nue(verbose=verbose)
!
call allocate_global3d(verbose=verbose)
!
!
call make_benchmark
!
!check if alo coefficients are correct
!call check_alo3d_o
!call check_alo3d
!
!call testlambs_direct
!
!---------perform radiative transfer for a given frequency point--------
!
nueindx=1
!
!set up opacities
select case(opt_opac)
   case(0)
      call setup_opacities1(nueindx, verbose=verbose)
   case(1)
      call setup_opacities2(verbose=verbose)
   case default
      stop 'error in main: opt_opac not set'
end select
!
!   
!
!set up 3d thermalization parameter
select case(opt_epsc)
   case(0)
      call setup_epsc1(verbose=verbose)
   case(1)
      call setup_epsc2(verbose=verbose)
   case default
      stop 'error in main: opt_opac not set'
end select
!
!set up planck function
call setup_bnue(nueindx, verbose=verbose)
!
call setup_bcondition3d(nueindx, verbose=verbose)
!
call print_model(nueindx, verbose=verbose)
!
select case(opt_method)
   case(4,5,6,7)
!check the grid and perform radiative transfer
      call check_grid3d(nueindx, verbose=verbose)   
      call conttrans_sc3d(nueindx, verbose=verbose)
   case(14,15,16,17)
!simply set the source function to lte value
      call conttrans_lte3d(nueindx, verbose=verbose)
   case default
      stop 'error in main: opt_method not properly defined'
end select
!
!
call print_solution(nueindx, verbose=verbose)
!
!
call output(verbose=verbose)
!
te_tot=omp_get_wtime()
!
write(*,*) '-------------------------------------------------------------------------------'
write(*,*)
!
write(*,*) 'total computation time:', te_tot-ts_tot
write(*,*)
!
write(*,*) '---------------------------main program done-----------------------------------'
write(*,*)
!
end program prog
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine make_benchmark
!!
!use prog_type
!use mod_benchmark, only: benchmark_mod
!use options, only: opt_method
!!
!implicit none
!!
!! ... local scalars
!integer(i4b) :: i
!!
!!calculate delx, dely, delz grids if not done yet
!!
!select case(benchmark_mod)
!   case(1)
!      write(*,*) '-----------performing benchmark model 1: searchligh beam test 2d---------------'
!      write(*,*)
!!calculating solution
!      call benchmark01_solution
!!output to file
!      call output_benchmark01
!
!!
!   case default
!      write(*,*) '----------------------no benchmark is being performed--------------------------'
!      write(*,*)
!      return
!end select
!!
!!stop program, because opacities and source functions have been overwritten
!if(benchmark_mod.ne.0) then
!   write(*,*)
!   stop '----------------------benchmark and main program done--------------------------'
!endif
!!
!!
!end subroutine make_benchmark
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine benchmark01_solution
!!
!use prog_type
!use fund_const
!use dime2d, only: z, int2d, scont2d, opac2d, alocont_o_nn2d, alocont_nn2d_tmp, &
!                  normalization2d_tmp, mint2d_tmp, nx, nz, t2d, bnue2d, intbound2d, zmin, zmax
!use mod_angles, only: nomega, nodes_mu, n_z, n_x
!use mod_frequencies, only: nodes_nue
!use ng_extra, only: ng_const
!use options, only: opt_ng_cont, opt_ait_cont, opt_method
!use mod_benchmark, only: int2d_theo, itmaxi, nconvi, devmaxi, epsmaxi_arr, nn_z, nn_x
!use params_input, only: kcont
!!
!implicit none
!!
!! ... local scalars
!integer(i4b) :: i, j, k, err, ix, iz, oindx
!integer(i4b) :: startz, endz, gamma
!integer(i4b) :: s1, s2, s3, s4, s4b
!real(dp) :: mu, eps_max, dtau, ts, te
!!
!! ... local arrays
!real(dp), dimension(:,:), allocatable :: eps2d
!real(dp), dimension(:,:), allocatable :: intbound2d_ng, intbound2d_single
!logical(dp), dimension(:,:), allocatable :: imaskbound2d
!!
!! ... local functions
!real(dp) :: bnue
!!
!! ... local characters
!!
!call cpu_time(ts)
!!
!!use constant opacity (k=1 corresponds to tau=1) for test cases
!dtau=one
!opac2d=kcont*dtau/(zmax-zmin)
!!
!nodes_mu=nn_z
!n_z=nn_z
!n_x=-sqrt(one-n_z**2)
!oindx=1
!
!if(.not.allocated(int2d)) allocate(int2d(nx,nz))
!if(.not.allocated(alocont_o_nn2d)) allocate(alocont_o_nn2d(nx,nz,27))
!if(.not.allocated(alocont_nn2d_tmp)) allocate(alocont_nn2d_tmp(nx,nz,27))
!if(.not.allocated(normalization2d_tmp)) allocate(normalization2d_tmp(nx,nz))
!if(.not.allocated(mint2d_tmp)) allocate(mint2d_tmp(nx,nz))
!allocate(eps2d(2,nz), stat=err)
!allocate(imaskbound2d(2,nz), stat=err)
!allocate(intbound2d_single(2,nz), stat=err)
!allocate(intbound2d_ng(4,2*nz), stat=err)
!eps2d=zero
!intbound2d_ng=zero
!imaskbound2d=1
!!
!scont2d=zero
!!calculating planck function
!do i=1, nx
!   do j=1, nz
!      bnue2d(i,j) = bnue(nodes_nue(1), t2d(i,j))
!   enddo
!enddo
!!
!!-----------------------------------------------------------------------
!!
!!initialisation of iteration step at which the old solution has to be 
!!stored (for ng-extrapolation/aitken extrapolation)
!!
!s1=1
!s2=2
!s3=3
!s4=4
!s4b=4
!
!epsmaxi_arr=zero
!!
!!-----------------------------------------------------------------------
!!
!!************************start iteration scheme*************************
!!
!do i=1, itmaxi
!!
!   if(.not.allocated(alocont_nn2d_tmp)) allocate(alocont_nn2d_tmp(nx,nz,27))
!   if(.not.allocated(normalization2d_tmp)) allocate(normalization2d_tmp(nx,nz))
!   if(.not.allocated(mint2d_tmp)) allocate(mint2d_tmp(nx,nz))
!!
!   eps2d=intbound2d(:,:,oindx)
!!
!   write(*,*) '--------------calculating intensity for a given ray----------------------------'
!   write(*,*) 'step', i
!!   write(*,*)
!!   write(*,fmt='(a5, 4(a20))') '#', 'z', 'intboundary(1)', 'intboundary(2)', 'int_theo (chi=const)'
!!   do j=1, nz
!!      write(*,fmt='(i5, 10(e20.8))') j, z(j), intbound2d(1,j,oindx), intbound2d(2,j,oindx), int2d(nx/2+1,2)*exp(-opac2d(nx/2+1,2)*(z(j)-z(2))/n_z(oindx))
!!   end do
!!   write(*,*)
!
!   if(opt_method.eq.4) then
!      call fsc_cont2d_lin(oindx,1)
!   elseif(opt_method.eq.5) then
!      call fsc_cont2d(oindx,1)
!   endif
!   intbound2d_single=intbound2d(:,:,oindx)
!!
!!-------------calculating percentage-error of mean intensities----------
!!
!   write(*,*) '-----------------------------calculating deviation-----------------------------'
!   write(*,*)
!   call calc_dev2d(eps2d, intbound2d_single, imaskbound2d, 2, nz, eps_max, ix, iz)
!   epsmaxi_arr(i)=eps_max
!   write(*,'(a30, 2i4, f8.4, es18.8)') 'max (dev) at grid-point:', ix, iz, z(iz) , eps_max
!   write(*,*)
!!
!   nconvi=i
!!
!   if(abs(eps_max).lt.devmaxi) then
!      write(*,*) "convergence after iteration no. ", i
!      write(*,*) "max (dev): ", eps_max
!      write(*,*)
!      exit
!   else if(i.eq.itmaxi) then
!      write(*,*) "no convergence after iteration no. ", i
!      write(*,*)
!   end if
!!
!!------------extrapolation of old subsequent source functions-----------
!!
!   if(opt_ng_cont.or.opt_ait_cont) then
!!
!!----------storing old source-functions for ng-extrapolation------------
!!
!      if(i.eq.s1) then
!         write(*,*) '----------------------storing source fct. at step n-3--------------------------'
!         write(*,*)
!         call store_ng2d(1,intbound2d_ng,2,nz,intbound2d_single)
!         s1=s1+ng_const
!      elseif(i.eq.s2) then
!         write(*,*) '----------------------storing source fct. at step n-2--------------------------'
!         write(*,*)
!         call store_ng2d(2,intbound2d_ng,2,nz,intbound2d_single)
!         s2=s2+ng_const
!      elseif(i.eq.s3) then
!         write(*,*) '----------------------storing source fct. at step n-1--------------------------'
!         write(*,*)
!         call store_ng2d(3,intbound2d_ng,2,nz,intbound2d_single)
!         s3=s3+ng_const
!      elseif(i.eq.s4) then
!         write(*,*) '----------------------storing source fct. at step n----------------------------'
!         write(*,*)
!         call store_ng2d(4,intbound2d_ng,2,nz,intbound2d_single)
!         s4=s4+ng_const
!
!         
!         if(opt_ng_cont) call ng_expol2d(intbound2d_ng,2,nz,intbound2d_single)
!         if(opt_ait_cont) call ait_expol2d(intbound2d_ng,2,nz,intbound2d_single)
!         intbound2d(:,:,oindx)=intbound2d_single
!      endif
!   endif
!   
!enddo
!!
!deallocate(alocont_nn2d_tmp)
!deallocate(normalization2d_tmp)
!deallocate(mint2d_tmp)
!!
!call cpu_time(te)
!write(*,*)
!write(*,*) 'total computation time', te-ts
!write(*,*)
!
!-----------------set theoretical intensities (for plane-parallel atmosphere)--------------
!
!if(n_z(oindx).gt.zero) then
!   startz = 3
!   endz = nz-2
!   gamma =  1
!elseif(n_z(oindx).lt.zero) then
!   startz = nz-2
!   endz = 3
!endif
!
!allocate(int2d_theo(nx,nz),stat=err)
!int2d_theo=zero
!int2d_theo(:,startz-gamma)=int2d(:,startz-gamma)
!int2d_theo(:,startz-2*gamma)=int2d(:,startz-2*gamma)
!!
!i=nx/2+1
!do k=startz, endz, gamma
!   dtau = half*(opac2d(i,k)+opac2d(i,k-gamma))*(z(k)-z(k-gamma))/n_z(oindx)
!   int2d_theo(:,k)=int2d_theo(:,k-gamma)*exp(-dtau)
!enddo
!!
!end subroutine benchmark01_solution
!!
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!!
!subroutine testlambs_direct
!!
!use prog_type
!use fund_const
!use dime2d, only: int2d, alocont_o_nn2d, alocont_nn2d, mint2d, normalization2d, &
!                  alocont_nn2d_tmp, mint2d_tmp, normalization2d_tmp, bnue2d, t2d, &
!                  nx, nz, imask2d, scont2d, imaskb2d, z
!use mod_angles, only: nomega, n_z
!use mod_frequencies, only: nodes_nue
!use mod_debug, only: iindx, kindx
!use options, only: opt_method
!use params_input, only: eps_cont
!use iter, only: itmaxc, devmaxc, epsmaxc_arr, nconvc
!
!!
!implicit none
!!
!!local scalars
!integer(i4b) :: nit, j, i, ii, k, kk, oindx, indx_1d_col, indx_1d_row, err
!real(dp) :: eps_max, fdum1, fdum2, fdum3, scont_new, wor, sum, sum2
!!
!!local arrays
!real(dp), dimension(:,:), allocatable :: lambda_mat
!real(dp), dimension(:), allocatable :: phi_boundary
!real(dp), dimension(:,:), allocatable :: a_mat, unity_mat, eps_mat, zeta_mat
!real(dp), dimension(:), allocatable :: bnue_vec
!real(dp), dimension(:), allocatable :: scont_vec, b_vec, mint_vec, eps_vec
!!
!!local functions
!real(dp) :: bnue
!
!!
!allocate(lambda_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
!allocate(a_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
!allocate(zeta_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
!allocate(unity_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
!allocate(eps_mat(nx*nz-2*nz,nx*nz-2*nz), stat=err)
!allocate(bnue_vec(nx*nz-2*nz), stat=err)
!allocate(scont_vec(nx*nz-2*nz), stat=err)
!allocate(b_vec(nx*nz-2*nz), stat=err)
!allocate(phi_boundary(nx*nz-2*nz), stat=err)
!allocate(mint_vec(nx*nz-2*nz), stat=err)
!allocate(eps_vec(nx*nz-2*nz), stat=err)
!!
!!
!!prepare bnue array
!do i=1, nx
!   do k=1, nz
!!      bnue2d(i,k) = one!bnue(nodes_nue(1), t2d(i,k))
!      bnue2d(i,k) = bnue(nodes_nue(1), t2d(i,k))
!   enddo
!enddo
!!
!alocont_nn2d=zero
!normalization2d=zero
!mint2d=zero
!!
!if(allocated(int2d)) deallocate(int2d)
!if(allocated(alocont_o_nn2d)) deallocate(alocont_o_nn2d)
!if(allocated(alocont_nn2d_tmp)) deallocate(alocont_nn2d_tmp)
!if(allocated(normalization2d_tmp)) deallocate(normalization2d_tmp)
!if(allocated(mint2d_tmp)) deallocate(mint2d_tmp)
!!
!allocate(int2d(nx,nz))
!allocate(alocont_o_nn2d(nx,nz,27))
!allocate(alocont_nn2d_tmp(nx,nz,27))
!allocate(normalization2d_tmp(nx,nz))
!allocate(mint2d_tmp(nx,nz))
!
!alocont_nn2d_tmp=zero
!normalization2d_tmp=zero
!mint2d_tmp=zero
!
!goto 40
!!
!!--------------------calculate lambda matrix-----------------------
!!
!10 continue
!!note: set boundary to zero
!do i=1, nx-2
!   do k=1, nz
!
!     alocont_nn2d_tmp=zero
!     normalization2d_tmp=zero
!     mint2d_tmp=zero
!     scont2d=zero
!!     scont2d(i,k)=one
!     call conv_indx_2d_to_1d(i,k, nx-2, indx_1d_col)      
!!
!     select case(imaskb2d(i,k))
!         case(5,6,9)
!            scont2d(i,k)=one
!         case(3,7)
!            scont2d(1,k)=one
!            scont2d(nx-1,k)=one
!!            scont2d(i,k)=one            
!         case(4,8)
!            scont2d(2,k)=one
!            scont2d(nx,k)=one
!!            scont2d(i,k)=one            
!         case(1,2)
!            scont2d(i,k)=one
!            if(i.eq.1) scont2d(nx-1,k)=one
!            if(i.eq.2) scont2d(nx,k)=one
!      end select
!
!      do oindx=1, nomega
!         if(opt_method.eq.4) then
!            call fsc_cont2d_lin(oindx,1)
!         elseif(opt_method.eq.5) then
!            call fsc_cont2d(oindx,1)
!         endif
!      enddo
!
!!      write(*,*) i, k, mint2d_tmp(i,k), alocont_nn2d_tmp(i,k,14), normalization2d_tmp(i,k)
!!
!!store corresponding matrix elements      
!      do ii=1, nx-2
!         do kk=1, nz
!            call conv_indx_2d_to_1d(ii,kk, nx-2, indx_1d_row)
!            lambda_mat(indx_1d_row,indx_1d_col)=mint2d_tmp(ii,kk)
!         enddo
!      enddo
!   enddo
!enddo
!
!open(1,file='TRASH/lambda_mat', form='unformatted')
!   write(1) lambda_mat
!close(1)
!!
!stop 'lambda matrix calculated, returning'
!!
!!--------------------calculate boundary contribution---------------
!!
!20 continue
!!note: set source function to zero, and only boundary intensities
!scont2d=zero
!
!do oindx=1, nomega
!   if(opt_method.eq.4) then
!      call fsc_cont2d_lin(oindx,1)
!   elseif(opt_method.eq.5) then
!      call fsc_cont2d(oindx,1)
!   endif
!enddo
!      
!do ii=1, nx-2
!   do kk=1, nz
!      call conv_indx_2d_to_1d(ii,kk, nx-2, indx_1d_row)
!      phi_boundary(indx_1d_row)=mint2d_tmp(ii,kk)
!      bnue_vec(indx_1d_row)=bnue2d(ii,kk)
!      if(kk.eq.nz) bnue_vec(indx_1d_row)=zero
!      if(kk.eq.nz-1) bnue_vec(indx_1d_row)=zero
!      if(kk.eq.1) phi_boundary(indx_1d_row)=bnue2d(ii,kk)
!      if(kk.eq.2) phi_boundary(indx_1d_row)=bnue2d(ii,kk)
!!      write(*,*) ii, kk, bnue_vec(indx_1d_row), bnue2d(ii,kk)
!   enddo
!enddo
!!
!open(1,file='TRASH/phi_boundary', form='unformatted')
!   write(1) phi_boundary
!close(1)
!!
!open(1,file='TRASH/bnue_vec', form='unformatted')
!   write(1) bnue_vec
!close(1)
!
!stop 'boundary contribution calculated, returning'
!!
!!------------------set up matrix system----------------------------
!!
!30 continue
!!
!!
!!
!open(1,file='TRASH/lambda_mat', form='unformatted')
!   read(1) lambda_mat
!close(1)
!open(1,file='TRASH/phi_boundary', form='unformatted')
!   read(1) phi_boundary
!close(1)
!open(1,file='TRASH/bnue_vec', form='unformatted')
!   read(1) bnue_vec
!close(1)
!!
!!
!!
!eps_mat=zero
!unity_mat=zero
!do i=1, nx*nz-2*nz
!   eps_mat(i,i)=eps_cont
!   unity_mat(i,i)=one
!enddo
!zeta_mat=zero
!zeta_mat=unity_mat-eps_mat
!!
!a_mat=matmul(zeta_mat,lambda_mat)
!a_mat=unity_mat-a_mat
!!
!bnue_vec=matmul(eps_mat,bnue_vec)
!phi_boundary=matmul(zeta_mat,phi_boundary)
!b_vec=phi_boundary+bnue_vec
!!
!note: bnue_vec here serves only as dummy argument
!call dgesv(nx*nz-2*nz, 1, a_mat, nx*nz-2*nz, bnue_vec, b_vec, nx*nz-2*nz, err)
!scont_vec=b_vec
!!if(err.ne.0) stop 'error in testlambs_direct: no inversion possible'
!!
!!transform solution to 2d array
!do i=1, nx*nz-2*nz
!   call conv_indx_1d_to_2d (i, nx-2, ii, kk)
!   scont2d(ii,kk) = scont_vec(i)
!enddo
!!
!do k=1, nz
!   write(*,*) z(k), scont2d(1,k), scont2d(nx/2+1,k), scont2d(nx,k), scont2d(nx/2+1,k)/bnue2d(1,1)
!enddo
!!
!!dummy output
!epsmaxc_arr=one
!nconvc=10
!
!write(*,*) bnue2d(1,:)
!call output
!
!stop 'testlambs_direct done, returning'
!!
!!------------------solve iteratively-------------------------------
!!
!40 continue
!!
!open(1,file='TRASH/lambda_mat', form='unformatted')
!   read(1) lambda_mat
!close(1)
!open(1,file='TRASH/phi_boundary', form='unformatted')
!   read(1) phi_boundary
!close(1)
!open(1,file='TRASH/bnue_vec', form='unformatted')
!   read(1) bnue_vec
!close(1)
!
!check the lambda-matrix
!do i=1, nx*nz-2*nz
!   sum=zero
!   do j=1, nx*nz-2*nz
!      sum=sum+lambda_mat(i,j)
!   enddo
!   call conv_indx_1d_to_2d (i, nx-2, ii, kk)
!   write(*,*) i, ii, kk, sum, lambda_mat(i,i)
!enddo
!!stop
!!
!!start value
!scont_vec=zero
!do i=1, nx-2
!   do k=1, nz
!      call conv_indx_2d_to_1d(i,k,nx-2,indx_1d_row)
!      scont_vec(indx_1d_row)=bnue2d(i,k)
!   enddo
!enddo
!!
!
!!stop
!ii=1
!kk=17
!call conv_indx_2d_to_1d(ii,kk, nx-2, indx_1d_row)
!sum=zero
!sum2=zero
!do j=1, nx*nz-2*nz
!   call conv_indx_1d_to_2d (j, nx-2, ii, kk)
!   write(*,*) ii, kk, lambda_mat(indx_1d_row,j), scont_vec(j), lambda_mat(indx_1d_row,j)*scont_vec(j), phi_boundary(j)
!   sum=sum+lambda_mat(indx_1d_row,j)
!   sum2=sum2+lambda_mat(indx_1d_row,j)*scont_vec(j)
!enddo
!write(*,*) 'sum', sum, sum2, bnue2d(ii,kk), phi_boundary(indx_1d_row), sum*bnue2d(ii,kk)
!!stop
!
!!
!
!mint_vec=zero
!eps_vec=zero
!!
!!for underestimation of diagonal
!wor=one!.99d0
!!
!do nit=1, itmaxc
!   mint_vec=matmul(lambda_mat,scont_vec)+phi_boundary
!
!   call calc_dev(eps_vec, mint_vec, nx*nz-2*nz, eps_max)
!   eps_vec=mint_vec
!   write(*,*) 'maximum deviation', nit, eps_max
!   epsmaxc_arr(nit)=eps_max
!   nconvc=nit
!   if(abs(eps_max).lt.devmaxc) exit
!   
!   do i=1, nx*nz-2*nz
!      fdum1=(one-eps_cont)
!      fdum2=one-fdum1*lambda_mat(i,i)
!      scont_new=fdum1/fdum2*mint_vec(i)-fdum1/fdum2*lambda_mat(i,i)*scont_vec(i) + eps_cont/fdum2*bnue_vec(i)
!
!      fdum3=one/(one/fdum1-lambda_mat(i,i))
!      scont_new=fdum3*mint_vec(i)-fdum3*lambda_mat(i,i)*scont_vec(i) + eps_cont/fdum2*bnue_vec(i)
!      call conv_indx_1d_to_2d (i, nx-2, ii, kk)
!
!!      write(*,*) ii, kk, fdum1/fdum2, fdum3, lambda_mat(i,i), mint_vec(i), phi_boundary(i), scont_new/bnue2d(ii,kk), one-lambda_mat(i,i)
!!      if(scont_new/bnue2d(ii,kk).gt.one) write(*,*) ii, kk, scont_new, scont_new/bnue2d(ii,kk), mint_vec(i), bnue2d(ii,kk), lambda_mat(i,i)
!!      if(scont_new/bnue2d(ii,kk).gt.one) stop 'error in bla'
!      scont_vec(i)=scont_new
!      
!   enddo
!
!   do i=1, nx*nz-2*nz
!      call conv_indx_1d_to_2d (i, nx, ii, kk)
!      if(ii.eq.nx/2+1) write(*,*) ii, kk, scont_vec(i), mint_vec(i), phi_boundary(i)
!   enddo
!
!!   if(nit.ge.1) stop 'go on in testlambs_direct, step 40'
!   
!enddo
!!
!!transform solution to 2d array
!do i=1, nx*nz-2*nz
!   call conv_indx_1d_to_2d (i, nx-2, ii, kk)
!   scont2d(ii,kk) = scont_vec(i)
!enddo
!!
!do k=1, nz
!   write(*,*) z(k), scont2d(1,k), scont2d(nx/2+1,k), scont2d(nx,k)
!enddo
!!
!!dummy output
!call output
!
!stop 'testlambs_direct iteratively done, returning'
!!
!!
!end subroutine testlambs_direct
