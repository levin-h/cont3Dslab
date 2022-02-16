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
subroutine set_boundary1d(nue,n_z)
!
use prog_type
use dime1d, only: z, int1d, t1d, scont1d, bnue1d
use mod_math, only: bnue
!
implicit none
!
! ... arguments
real(dp), intent(in) :: nue, n_z
!
! ... local scalars
integer(i4b) :: i
!
! ... local functions
!
int1d=0.d0
!
int1d(1)=bnue(nue,t1d(1))
int1d(2)=bnue(nue,t1d(2))
!
!
!scont1d=1.d0
!int1d(1)=1.d0
!bnue1d=1.d0
!
end subroutine set_boundary1d
!
!***********************************************************************
!***********************************************************************
!
!                 CONTINUUM TRASNFER ROUTINES
!
!***********************************************************************
!***********************************************************************
!
subroutine fsc_cont1d(oindx,nueindx)
!
!-----------------------------------------------------------------------
!------short characteristics for continuum radiative transfer in 1d-----
!-----------calculating intensties for given mu,nue specified-----------
!--------------------by input oindx and nueindx-------------------------
!-----------------------------------------------------------------------
!
use prog_type

use fund_const
use dime1d, only: nz, z
use dime1d, only: int1d, opac1d, scont1d, alocont_o_nn1d, alocont_nn1d_tmp, &
     mint1d_tmp, normalization1d_tmp, imask1d
use angles, only: nodes_mu, weight_mu, q_alo
use freq, only: nodes_nue, xic1_nue
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: oindx, nueindx
!
! ... local scalars
integer(i4b) :: k, kkm1, kkp1
integer(i4b) :: gamma
integer(i4b) :: startz, endz
real(dp) :: nn_z, xnue, wall
real(dp) :: opac_p, scont_p
real(dp) :: dels_u, int_u, opac_u, scont_u
real(dp) :: dels_d, opac_d, scont_d
real(dp) :: abs_sc, int_sc, contr_sc
real(dp) :: alo_u, alo_p, alo_d
integer :: q5, q14, q23
!
!for debugging
!
! ... local arrays
!
! ... local functions
!
! ... local characters
!
! ... local logicals
!
!frequency
xnue=nodes_nue(nueindx)
!
!directions
nn_z=nodes_mu(oindx)
!
!angulare integration weight
wall=weight_mu(oindx)
!
!indices for nearest neighbour alo
q5=q_alo(oindx,1)
q14=q_alo(oindx,2)
q23=q_alo(oindx,3)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
if(nn_z.gt.0.d0) then
   startz = 3
   endz = nz-2
   gamma=  1
elseif(nn_z.lt.0.d0) then
   startz = nz-2
   endz = 3
   gamma=-1
endif
!
!--------------------reset the intensities and alo----------------------
!
call set_boundary1d(xnue, nn_z)
!
alocont_o_nn1d=0.d0
!
!-----------------------------------------------------------------------
!
do k=startz, endz, gamma
!
   kkm1=k-gamma
   kkp1=k+gamma
   
   int_u = int1d(kkm1)
   scont_u=scont1d(kkm1)
   opac_u=opac1d(kkm1)
   
   scont_p=scont1d(k)
   opac_p=opac1d(k)

   scont_d=scont1d(kkp1)
   opac_d=opac1d(kkp1)

   dels_u=(z(k)-z(kkm1))/nn_z
   dels_d=(z(kkp1)-z(k))/nn_z


   call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                   dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
   int1d(k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p + alo_d*scont_d
!   call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!   int1d(k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p
!   alo_d=0.d0

   alocont_o_nn1d(kkp1,q5) = alo_d*imask1d(kkp1)
   alocont_o_nn1d(k,q14) = (alo_p + abs_sc*alocont_o_nn1d(k,q5))*imask1d(k)
   alocont_o_nn1d(kkm1,q23) = (alo_u + abs_sc*alocont_o_nn1d(kkm1,q14))*imask1d(kkm1)

!   write(*,*) k, alocont_o_nn1d(k,q14), alocont_o_nn1d(k,q5), abs_sc, alo_u, alo_p, alo_d!, alo_u+alo_p+alo_d
!   if(kkm1.eq. nz-2) then
!      write(*,*) q5, q14, q23, alocont_o_nn1d(kkp1, q14), alocont_o_nn1d(kkp1,q5), alocont_o_nn1d(kkp1,q23)
!      stop
   !   endif
!   write(*,*) k, alo_d!, alocont_o
   
   
!perform angular integration
!*********debug start
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
!               if(dels_r.ne.1.d10) then
!                  int3d(i,j,k)=xic1
!               else
!                  int3d(i,j,k)=0.d0
!               endif
!*********debug end
!
   mint1d_tmp(k) = mint1d_tmp(k) + int1d(k)*wall
   normalization1d_tmp(k) = normalization1d_tmp(k) + wall
   alocont_nn1d_tmp(kkp1,q5) = alocont_nn1d_tmp(kkp1,q5) + wall*alocont_o_nn1d(kkp1,q5)
   alocont_nn1d_tmp(k,q14) = alocont_nn1d_tmp(k,q14) + wall*alocont_o_nn1d(k,q14)
   alocont_nn1d_tmp(kkm1,q23) = alocont_nn1d_tmp(kkm1,q23) + wall*alocont_o_nn1d(kkm1,q23)
!
!   write(*,'(i5, 10es20.8)') k, dels_u, dels_d, alo_u, alo_p, alo_d, alo_u+alo_p+alo_d, abs_sc, alocont_nn1d_tmp(k,q14), mint1d_tmp(k)
!   
enddo
!
!write(*,*) nn_z, q5, q14, q23
!do k=1, nz
!   write(*,*) alocont_o_nn1d(k,q5), alocont_o_nn1d(k,q14), alocont_o_nn1d(k,q23)
!enddo
!stop
!write(*,*) nn_z
!write(*,'(10es20.8)') alocont_nn1d_tmp(1,2), alocont_nn1d_tmp(1,3), 0., 0.
!write(*,'(10es20.8)') alocont_nn1d_tmp(2,1), alocont_nn1d_tmp(2,2), alocont_nn1d_tmp(2,3), 0.
!write(*,'(10es20.8)') 0.d0, alocont_nn1d_tmp(3,1), alocont_nn1d_tmp(3,2), alocont_nn1d_tmp(3,3)
!write(*,'(10es20.8)') 0.d0, 0.d0, alocont_nn1d_tmp(4,1), alocont_nn1d_tmp(4,2)
!write(*,*)
!stop
!
end subroutine fsc_cont1d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_cont1d_lin(oindx,nueindx)
!
!-----------------------------------------------------------------------
!------short characteristics for continuum radiative transfer in 1d-----
!-----------calculating intensties for given mu,nue specified-----------
!--------------------by input oindx and nueindx-------------------------
!-----------------------------------------------------------------------
!
use prog_type

use fund_const
use dime1d, only: nz, z
use dime1d, only: int1d, opac1d, scont1d, alocont_o_nn1d, alocont_nn1d_tmp, &
     mint1d_tmp, normalization1d_tmp, imask1d
use angles, only: nodes_mu, weight_mu, q_alo
use freq, only: nodes_nue, xic1_nue
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: oindx, nueindx
!
! ... local scalars
integer(i4b) :: k, kkm1, kkp1
integer(i4b) :: gamma
integer(i4b) :: startz, endz
real(dp) :: nn_z, xnue, wall
real(dp) :: opac_p, scont_p
real(dp) :: dels_u, int_u, opac_u, scont_u
real(dp) :: dels_d, opac_d, scont_d
real(dp) :: abs_sc, int_sc, contr_sc
real(dp) :: alo_u, alo_p, alo_d
integer :: q5, q14, q23
!
!for debugging
!
! ... local arrays
!
! ... local functions
!
! ... local characters
!
! ... local logicals
!
!frequency
xnue=nodes_nue(nueindx)
!
!directions
nn_z=nodes_mu(oindx)
!
!angulare integration weight
wall=weight_mu(oindx)
!
!indices for nearest neighbour alo
q5=q_alo(oindx,1)
q14=q_alo(oindx,2)
q23=q_alo(oindx,3)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
if(nn_z.gt.0.d0) then
   startz = 3
   endz = nz-2
   gamma=  1
elseif(nn_z.lt.0.d0) then
   startz = nz-2
   endz = 3
   gamma=-1
endif
!
!--------------------reset the intensities and alo----------------------
!
call set_boundary1d(xnue, nn_z)
!
alocont_o_nn1d=0.d0
!
!-----------------------------------------------------------------------
!
do k=startz, endz, gamma
!
   kkm1=k-gamma
   kkp1=k+gamma
   
   int_u = int1d(kkm1)
   scont_u=scont1d(kkm1)
   opac_u=opac1d(kkm1)
   
   scont_p=scont1d(k)
   opac_p=opac1d(k)

   scont_d=scont1d(kkp1)
   opac_d=opac1d(kkp1)

   dels_u=(z(k)-z(kkm1))/nn_z
   dels_d=(z(kkp1)-z(k))/nn_z


!   call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                   dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!   int1d(k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p + alo_d*scont_d
   call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
   int1d(k) = abs_sc*int_u + alo_u*scont_u + alo_p*scont_p
   alo_d=0.d0

   alocont_o_nn1d(kkp1,q5) = alo_d*imask1d(kkp1)
   alocont_o_nn1d(k,q14) = (alo_p + abs_sc*alocont_o_nn1d(k,q5))*imask1d(k)
   alocont_o_nn1d(kkm1,q23) = (alo_u + abs_sc*alocont_o_nn1d(kkm1,q14))*imask1d(kkm1)

   
!perform angular integration
!*********debug start
!               dels_r=dist_sphere(nn_x,nn_y,nn_z,x(i),y(j),z(k))
!               if(dels_r.ne.1.d10) then
!                  int3d(i,j,k)=xic1
!               else
!                  int3d(i,j,k)=0.d0
!               endif
!*********debug end
!
   mint1d_tmp(k) = mint1d_tmp(k) + int1d(k)*wall
   normalization1d_tmp(k) = normalization1d_tmp(k) + wall
   alocont_nn1d_tmp(kkp1,q5) = alocont_nn1d_tmp(kkp1,q5) + wall*alocont_o_nn1d(kkp1,q5)
   alocont_nn1d_tmp(k,q14) = alocont_nn1d_tmp(k,q14) + wall*alocont_o_nn1d(k,q14)
   alocont_nn1d_tmp(kkm1,q23) = alocont_nn1d_tmp(kkm1,q23) + wall*alocont_o_nn1d(kkm1,q23)
!
!   write(*,'(i5, 13es20.8)') k, alocont_nn1d_tmp(kkp1,q5), alocont_nn1d_tmp(k,q14), alocont_nn1d_tmp(kkm1,q23)
!   
enddo
!

!stop 'go on in fsc_cont1d_lin'

!write(*,*) nn_z
!write(*,'(10es20.8)') alocont_nn1d_tmp(1,2), alocont_nn1d_tmp(1,3), 0., 0.
!write(*,'(10es20.8)') alocont_nn1d_tmp(2,1), alocont_nn1d_tmp(2,2), alocont_nn1d_tmp(2,3), 0.
!write(*,'(10es20.8)') 0.d0, alocont_nn1d_tmp(3,1), alocont_nn1d_tmp(3,2), alocont_nn1d_tmp(3,3)
!write(*,'(10es20.8)') 0.d0, 0.d0, alocont_nn1d_tmp(4,1), alocont_nn1d_tmp(4,2)
!write(*,*)
!stop
!
end subroutine fsc_cont1d_lin
