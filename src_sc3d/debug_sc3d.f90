!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_alo3d_o
!
use prog_type
use fund_const
use mod_grid3d, only: int3d, alocont_o_nn3d, alocont_nn3d, mint3d, normalization3d, &
                  alocont_nn3d_tmp, mint3d_tmp, normalization3d_tmp, bnue3d, tgas3d, trad3d, &
                  nx, ny, nz, imask3d, scont3d, imaskb3d, &
                  fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp
use mod_angles, only: nomega, n_x, n_y, n_z
use mod_frequencies, only: nodes_nue
use options, only: opt_method
use mod_model3d, only: setup_opacities1
use mod_conttrans3d, only: formallc_cont3d, formalsc_cont3d_lin, formallc_cont3d_lin, formalsc_cont3d
!
implicit none
!
! ... local scalars
integer(i4b) :: istart, jstart, kstart, iend, jend, kend, i, j, k, oindx, nueindx
!
! ... local arrays
!
alocont_nn3d=zero
normalization3d=zero
mint3d=zero
bnue3d=zero
!
if(allocated(int3d)) deallocate(int3d)
if(allocated(alocont_o_nn3d)) deallocate(alocont_o_nn3d)
if(allocated(alocont_nn3d_tmp)) deallocate(alocont_nn3d_tmp)
if(allocated(normalization3d_tmp)) deallocate(normalization3d_tmp)
if(allocated(mint3d_tmp)) deallocate(mint3d_tmp)
if(allocated(fcontx3d_tmp)) deallocate(fcontx3d_tmp)
if(allocated(fconty3d_tmp)) deallocate(fconty3d_tmp)
if(allocated(fcontz3d_tmp)) deallocate(fcontz3d_tmp)
!
allocate(int3d(nx,ny,nz))
allocate(alocont_o_nn3d(nx,ny,nz,27))
allocate(alocont_nn3d_tmp(nx,ny,nz,27))
allocate(normalization3d_tmp(nx,ny,nz))
allocate(mint3d_tmp(nx,ny,nz))
allocate(fcontx3d_tmp(nx,ny,nz))
allocate(fconty3d_tmp(nx,ny,nz))
allocate(fcontz3d_tmp(nx,ny,nz))
!
nueindx=1
!
!set up opacities
call setup_opacities1(nueindx)
!
istart=1
iend=nx
jstart=1
jend=ny
kstart=2
kend=nz-1


istart=1
iend=1
jstart=1
jend=1
kstart=16
kend=16

!istart=1
!iend=nx
!jstart=1
!jend=ny
!kstart=5
!kend=5

write(*,*) 'possible directions'
do oindx=1, nomega
   write(*,*) oindx, n_x(oindx), n_y(oindx), n_z(oindx)
enddo
write(*,*)

do i=istart, iend
   do j=jstart, jend
      do k=kstart, kend
!
         scont3d=zero
         alocont_nn3d_tmp=zero
         mint3d_tmp=zero
         normalization3d_tmp=zero
         select case(imaskb3d(i,j,k))
            case(1,2)
               scont3d(i,j,k)=one
            case(3,4,5,6,7,8,9,10,11)
               scont3d(i,j,k)=one
            case(12,13,16,19,22,25)
               scont3d(1,j,k)=one
               scont3d(nx,j,k)=one
            case(14,15,17,20,23,26)
               scont3d(i,1,k)=one
               scont3d(i,ny,k)=one
            case(18,21,24,27)
               scont3d(1,1,k)=one
               scont3d(1,ny,k)=one
               scont3d(nx,1,k)=one
               scont3d(nx,ny,k)=one
            case default
         end select
!
!         do oindx=nomega, 1, -1
         do oindx=5,5
            write(*,*) 'oindx, n_x, n_y, n_z', oindx, n_x(oindx), n_y(oindx), n_z(oindx)
            select case(opt_method)
               case(4)
                  call formalsc_cont3d_lin(oindx,nueindx)
               case(5)
                  call formalsc_cont3d(oindx,nueindx)
               case(6)
                  call formallc_cont3d_lin(oindx,nueindx)
               case(7)
                  call formallc_cont3d(oindx,nueindx)
            end select
!            
            select case(imaskb3d(i,j,k))
               case(3,4,5,6,7,8,9,10,11) !standard alo
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
                  
               case(12,16,22) !left boundary
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(nx-1,j-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(nx-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j+1,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(nx-1,j+1,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j-1,k), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,10)-int3d(nx-1,j-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(nx-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j+1,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(nx-1,j+1,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(nx-1,j-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(nx-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j+1,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(nx-1,j+1,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(nx-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(nx-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(nx-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(nx-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(nx-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(nx-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(nx-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(nx-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(nx-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
                  
               case(13,19,25) !right boundary
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(2,j-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(2,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j+1,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(2,j+1,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(2,j-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(2,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j+1,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(2,j+1,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(2,j-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(2,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j+1,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(2,j+1,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(2,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(2,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(2,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(2,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(2,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(2,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(2,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(2,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(2,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
!
               case(14,17,20) !front boundary
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,ny-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(i-1,ny-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,ny-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,ny-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(i+1,ny-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,ny-1,k), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,10)-int3d(i-1,ny-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,ny-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,ny-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(i+1,ny-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,ny-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(i-1,ny-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,ny-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,ny-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(i+1,ny-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(i-1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(i+1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(i-1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(i+1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(i-1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(i+1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
!
               case(15,23,26) !back boundary
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,2,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(i-1,2,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,2,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,2,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(i+1,2,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,2,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(i-1,2,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,2,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,2,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(i+1,2,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,2,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(i-1,2,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,2,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,2,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(i+1,2,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(i-1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,2,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(i+1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(i-1,2,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,2,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(i+1,2,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(i-1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,2,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(i+1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
!
               case(18) !left front
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,ny-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(nx-1,ny-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,ny-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,ny-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(i+1,ny-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(nx-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j+1,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(nx-1,j+1,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,ny-1,k), alocont_o_nn3d(i,i,k,10), alocont_o_nn3d(i,j,k,10)-int3d(nx-1,ny-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,ny-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,ny-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(i+1,ny-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(nx-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j+1,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(nx-1,j+1,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,ny-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(nx-1,ny-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,ny-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,ny-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(i+1,ny-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(nx-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j+1,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(nx-1,j+1,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j+1,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(nx-1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(i+1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(nx-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(nx-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(i+1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(nx-1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(i+1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(nx-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(nx-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(i+1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(nx-1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(i+1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(nx-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(nx-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(i+1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
!
               case(21) !right front
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,ny-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(i-1,ny-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,ny-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,ny-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(2,ny-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(2,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j+1,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(2,j+1,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,ny-1,k), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,10)-int3d(i-1,ny-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,ny-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,ny-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(2,ny-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(2,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j+1,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(2,j+1,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,ny-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(i-1,ny-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,ny-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,ny-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,ny-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(2,ny-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(2,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j+1,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j+1,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j+1,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(2,j+1,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(i-1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(2,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(2,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(i-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(2,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(i-1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(2,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(2,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(i-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(2,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(i-1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(2,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(2,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(i-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(2,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
               case(24) !left back boundary
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(nx-1,j-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(nx-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,2,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(nx-1,2,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,2,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,2,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(i+1,2,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j-1,k), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,10)-int3d(nx-1,j-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(nx-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,2,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(nx-1,2,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,2,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,2,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(i+1,2,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(nx-1,j-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(nx-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(nx-1,2,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(nx-1,2,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,2,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(i+1,2,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(i+1,2,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(nx-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(i+1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(nx-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(nx-1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,2,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(i+1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(nx-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(i+1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(nx-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(nx-1,2,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,2,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(i+1,2,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(nx-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(i+1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(nx-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(nx-1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,2,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(i+1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
               case(27) !right back boundary
                  write(*,*) 1, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k-1), alocont_o_nn3d(i,j,k,1), alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)
                  write(*,*) 2, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k-1), alocont_o_nn3d(i,j,k,2), alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)
                  write(*,*) 3, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j-1,k-1), alocont_o_nn3d(i,j,k,3), alocont_o_nn3d(i,j,k,3)-int3d(2,j-1,k-1)
                  write(*,*) 4, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k-1), alocont_o_nn3d(i,j,k,4), alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)
                  write(*,*) 5, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k-1), alocont_o_nn3d(i,j,k,5), alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)
                  write(*,*) 6, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k-1), alocont_o_nn3d(i,j,k,6), alocont_o_nn3d(i,j,k,6)-int3d(2,j,k-1)
                  write(*,*) 7, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,2,k-1), alocont_o_nn3d(i,j,k,7), alocont_o_nn3d(i,j,k,7)-int3d(i-1,2,k-1)
                  write(*,*) 8, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k-1), alocont_o_nn3d(i,j,k,8), alocont_o_nn3d(i,j,k,8)-int3d(i,2,k-1)
                  write(*,*) 9, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,2,k-1), alocont_o_nn3d(i,j,k,9), alocont_o_nn3d(i,j,k,9)-int3d(2,2,k-1)
                  write(*,*) 10, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k), alocont_o_nn3d(i,j,k,10), alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)
                  write(*,*) 11, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k), alocont_o_nn3d(i,j,k,11), alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)
                  write(*,*) 12, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j-1,k), alocont_o_nn3d(i,j,k,12), alocont_o_nn3d(i,j,k,12)-int3d(2,j-1,k)
                  write(*,*) 13, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k), alocont_o_nn3d(i,j,k,13), alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)
                  write(*,*) 14, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k), alocont_o_nn3d(i,j,k,14), alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)
                  write(*,*) 15, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k), alocont_o_nn3d(i,j,k,15), alocont_o_nn3d(i,j,k,15)-int3d(2,j,k)
                  write(*,*) 16, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,2,k), alocont_o_nn3d(i,j,k,16), alocont_o_nn3d(i,j,k,16)-int3d(i-1,2,k)
                  write(*,*) 17, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k), alocont_o_nn3d(i,j,k,17), alocont_o_nn3d(i,j,k,17)-int3d(i,2,k)
                  write(*,*) 18, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,2,k), alocont_o_nn3d(i,j,k,18), alocont_o_nn3d(i,j,k,18)-int3d(2,2,k)
                  write(*,*) 19, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j-1,k+1), alocont_o_nn3d(i,j,k,19), alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)
                  write(*,*) 20, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j-1,k+1), alocont_o_nn3d(i,j,k,20), alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)
                  write(*,*) 21, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j-1,k+1), alocont_o_nn3d(i,j,k,21), alocont_o_nn3d(i,j,k,21)-int3d(2,j-1,k+1)
                  write(*,*) 22, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,j,k+1), alocont_o_nn3d(i,j,k,22), alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)
                  write(*,*) 23, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,j,k+1), alocont_o_nn3d(i,j,k,23), alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)
                  write(*,*) 24, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,j,k+1), alocont_o_nn3d(i,j,k,24), alocont_o_nn3d(i,j,k,24)-int3d(2,j,k+1)
                  write(*,*) 25, i, j, k, oindx, imaskb3d(i,j,k), int3d(i-1,2,k+1), alocont_o_nn3d(i,j,k,25), alocont_o_nn3d(i,j,k,25)-int3d(i-1,2,k+1)
                  write(*,*) 26, i, j, k, oindx, imaskb3d(i,j,k), int3d(i,2,k+1), alocont_o_nn3d(i,j,k,26), alocont_o_nn3d(i,j,k,26)-int3d(i,2,k+1)
                  write(*,*) 27, i, j, k, oindx, imaskb3d(i,j,k), int3d(2,2,k+1), alocont_o_nn3d(i,j,k,27), alocont_o_nn3d(i,j,k,27)-int3d(2,2,k+1)
                  write(*,*)
                  
                  if(abs(alocont_o_nn3d(i,j,k,1)-int3d(i-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
                  if(abs(alocont_o_nn3d(i,j,k,2)-int3d(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
                  if(abs(alocont_o_nn3d(i,j,k,3)-int3d(2,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
                  if(abs(alocont_o_nn3d(i,j,k,4)-int3d(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
                  if(abs(alocont_o_nn3d(i,j,k,5)-int3d(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
                  if(abs(alocont_o_nn3d(i,j,k,6)-int3d(2,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
                  if(abs(alocont_o_nn3d(i,j,k,7)-int3d(i-1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
                  if(abs(alocont_o_nn3d(i,j,k,8)-int3d(i,2,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
                  if(abs(alocont_o_nn3d(i,j,k,9)-int3d(2,2,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
                  if(abs(alocont_o_nn3d(i,j,k,10)-int3d(i-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
                  if(abs(alocont_o_nn3d(i,j,k,11)-int3d(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
                  if(abs(alocont_o_nn3d(i,j,k,12)-int3d(2,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
                  if(abs(alocont_o_nn3d(i,j,k,13)-int3d(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
                  if(abs(alocont_o_nn3d(i,j,k,14)-int3d(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
                  if(abs(alocont_o_nn3d(i,j,k,15)-int3d(2,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
                  if(abs(alocont_o_nn3d(i,j,k,16)-int3d(i-1,2,k)).gt.small_number) stop 'error in check_alo2d: 16'
                  if(abs(alocont_o_nn3d(i,j,k,17)-int3d(i,2,k)).gt.small_number) stop 'error in check_alo2d: 17'
                  if(abs(alocont_o_nn3d(i,j,k,18)-int3d(2,2,k)).gt.small_number) stop 'error in check_alo2d: 18'
                  if(abs(alocont_o_nn3d(i,j,k,19)-int3d(i-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
                  if(abs(alocont_o_nn3d(i,j,k,20)-int3d(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
                  if(abs(alocont_o_nn3d(i,j,k,21)-int3d(2,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
                  if(abs(alocont_o_nn3d(i,j,k,22)-int3d(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
                  if(abs(alocont_o_nn3d(i,j,k,23)-int3d(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
                  if(abs(alocont_o_nn3d(i,j,k,24)-int3d(2,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
                  if(abs(alocont_o_nn3d(i,j,k,25)-int3d(i-1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
                  if(abs(alocont_o_nn3d(i,j,k,26)-int3d(i,2,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
                  if(abs(alocont_o_nn3d(i,j,k,27)-int3d(2,2,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
               case default
            end select
         enddo
      enddo
   enddo
enddo
!
!
stop 'go on in check_alo3d_o'
!
!
end subroutine check_alo3d_o
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine check_alo3d
!
use prog_type
use fund_const
use mod_grid3d, only: int3d, alocont_o_nn3d, alocont_nn3d, mint3d, normalization3d, &
                  alocont_nn3d_tmp, mint3d_tmp, normalization3d_tmp, bnue3d, tgas3d, trad3d, &
                  nx, ny, nz, imask3d, scont3d, imaskb3d, &
                  fcontx3d_tmp, fconty3d_tmp, fcontz3d_tmp
use mod_angles, only: nomega, n_x, n_y, n_z
use mod_frequencies, only: nodes_nue
use options, only: opt_method
use mod_model3d, only: setup_opacities1
use mod_conttrans3d, only: formallc_cont3d, formalsc_cont3d_lin, formallc_cont3d_lin, formalsc_cont3d
!
implicit none
!
integer(i4b) :: istart, jstart, kstart, iend, jend, kend, i, j, k, oindx, nueindx
!
alocont_nn3d=zero
normalization3d=zero
mint3d=zero
bnue3d=zero
!
if(allocated(int3d)) deallocate(int3d)
if(allocated(alocont_o_nn3d)) deallocate(alocont_o_nn3d)
if(allocated(alocont_nn3d_tmp)) deallocate(alocont_nn3d_tmp)
if(allocated(normalization3d_tmp)) deallocate(normalization3d_tmp)
if(allocated(mint3d_tmp)) deallocate(mint3d_tmp)
if(allocated(fcontx3d_tmp)) deallocate(fcontx3d_tmp)
if(allocated(fconty3d_tmp)) deallocate(fconty3d_tmp)
if(allocated(fcontz3d_tmp)) deallocate(fcontz3d_tmp)
!
allocate(int3d(nx,ny,nz))
allocate(alocont_o_nn3d(nx,ny,nz,27))
allocate(alocont_nn3d_tmp(nx,ny,nz,27))
allocate(normalization3d_tmp(nx,ny,nz))
allocate(mint3d_tmp(nx,ny,nz))
allocate(fcontx3d_tmp(nx,ny,nz))
allocate(fconty3d_tmp(nx,ny,nz))
allocate(fcontz3d_tmp(nx,ny,nz))
!
nueindx=1
!
!set up opacities
call setup_opacities1(nueindx)
!
istart=1
iend=nx
jstart=1
jend=ny
kstart=2
kend=nz-1


!istart=10
!iend=10
!jstart=ny
!jend=ny
!kstart=5
!kend=5

istart=1
iend=nx
jstart=1
jend=1
kstart=15
kend=15


write(*,*) 'possible directions'
do oindx=1, nomega
   write(*,*) oindx, n_x(oindx), n_y(oindx), n_z(oindx)
enddo
write(*,*)

do i=istart, iend
   do j=jstart, jend
      do k=kstart, kend
!
         scont3d=zero
         alocont_nn3d_tmp=zero
         mint3d_tmp=zero
         normalization3d_tmp=zero
         select case(imaskb3d(i,j,k))
            case(1,2)
               scont3d(i,j,k)=one
            case(3,4,5,6,7,8,9,10,11)
               scont3d(i,j,k)=one
            case(12,13,16,19,22,25)
               scont3d(1,j,k)=one
               scont3d(nx,j,k)=one
            case(14,15,17,20,23,26)
               scont3d(i,1,k)=one
               scont3d(i,ny,k)=one
            case(18,21,24,27)
               scont3d(1,1,k)=one
               scont3d(1,ny,k)=one
               scont3d(nx,1,k)=one
               scont3d(nx,ny,k)=one
            case default
         end select
!
         do oindx=nomega, 1, -1
            select case(opt_method)
               case(4)
                  call formalsc_cont3d_lin(oindx,nueindx)
               case(5)
                  call formalsc_cont3d(oindx,nueindx)
               case(6)
                  call formallc_cont3d_lin(oindx,nueindx)
               case(7)
                  call formallc_cont3d(oindx,nueindx)
            end select
         enddo
!
         select case(imaskb3d(i,j,k))
            case(3,4,5,6,7,8,9,10,11) !standard alo
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,j-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,j-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,j+1,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,j+1,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k), alocont_nn3d_tmp(i,j,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,j-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,j-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,j+1,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,j+1,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,j-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,j-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,j+1,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,j+1,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
               
            case(12,16,22) !left boundary
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(nx-1,j-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,j-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(nx-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j+1,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(nx-1,j+1,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,j+1,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j-1,k), alocont_nn3d_tmp(i,j,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(nx-1,j-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,j-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(nx-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j+1,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(nx-1,j+1,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,j+1,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(nx-1,j-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,j-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(nx-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j+1,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(nx-1,j+1,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,j+1,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(nx-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(nx-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(nx-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(nx-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(nx-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(nx-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(nx-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(nx-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(nx-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
               
            case(13,19,25) !right boundary
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,j-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(2,j-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(2,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,j+1,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j+1,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(2,j+1,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k), alocont_nn3d_tmp(i,j,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,j-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(2,j-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(2,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,j+1,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j+1,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(2,j+1,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,j-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(2,j-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(2,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,j+1,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j+1,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(2,j+1,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(2,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(2,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(2,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(2,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(2,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(2,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(2,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(2,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(2,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'

            case(14,17,20) !front boundary
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,ny-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,ny-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,ny-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,ny-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,ny-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,j+1,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,j+1,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,ny-1,k), alocont_nn3d_tmp(i,j,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,ny-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,ny-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,ny-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,ny-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,j+1,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,j+1,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,ny-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,ny-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,ny-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,ny-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,ny-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,j+1,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,j+1,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'

            case(15,23,26) !back boundary
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,j-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,j-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,2,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,2,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,2,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,2,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,2,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k), alocont_nn3d_tmp(i,j,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,j-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,j-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,2,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,2,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,2,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,2,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,2,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,j-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,j-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,2,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,2,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,2,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,2,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,2,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,2,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,2,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,2,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,2,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,2,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 27'

            case(18) !left front
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,ny-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(nx-1,ny-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,ny-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,ny-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,ny-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(nx-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j+1,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(nx-1,j+1,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,j+1,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,ny-1,k), alocont_nn3d_tmp(i,i,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(nx-1,ny-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,ny-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,ny-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,ny-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(nx-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j+1,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(nx-1,j+1,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,j+1,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,ny-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(nx-1,ny-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,ny-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,ny-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,ny-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(nx-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j+1,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(nx-1,j+1,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j+1,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,j+1,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(nx-1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(nx-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(nx-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(nx-1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(nx-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(nx-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(nx-1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(nx-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(nx-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'

            case(21) !right front
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,ny-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,ny-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,ny-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,ny-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(2,ny-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(2,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,j+1,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j+1,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(2,j+1,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,ny-1,k), alocont_nn3d_tmp(i,j,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,ny-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,ny-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,ny-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(2,ny-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(2,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,j+1,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j+1,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(2,j+1,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,ny-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,ny-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,ny-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,ny-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,ny-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(2,ny-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(2,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j+1,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,j+1,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j+1,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j+1,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(2,j+1,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(2,ny-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(2,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(2,j+1,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(2,ny-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(2,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,j+1,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,j+1,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(2,j+1,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(2,ny-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(2,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(2,j+1,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
            case(24) !left back boundary
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(nx-1,j-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,j-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(nx-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,2,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(nx-1,2,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,2,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,2,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,2,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j-1,k), alocont_nn3d_tmp(i,j,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(nx-1,j-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,j-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(nx-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,2,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(nx-1,2,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,2,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,2,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,2,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(nx-1,j-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,j-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(nx-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(nx-1,2,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(nx-1,2,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,2,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i+1,2,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,2,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(nx-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(i+1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(nx-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(i+1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(nx-1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,2,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(i+1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(nx-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(i+1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(nx-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(i+1,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(nx-1,2,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,2,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(i+1,2,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(nx-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(i+1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(nx-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(i+1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(nx-1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,2,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(i+1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
            case(27) !right back boundary
               write(*,*) 1, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k-1), alocont_nn3d_tmp(i,j,k,1), alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,j-1,k-1)
               write(*,*) 2, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k-1), alocont_nn3d_tmp(i,j,k,2), alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)
               write(*,*) 3, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j-1,k-1), alocont_nn3d_tmp(i,j,k,3), alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(2,j-1,k-1)
               write(*,*) 4, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k-1), alocont_nn3d_tmp(i,j,k,4), alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)
               write(*,*) 5, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k-1), alocont_nn3d_tmp(i,j,k,5), alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)
               write(*,*) 6, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k-1), alocont_nn3d_tmp(i,j,k,6), alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(2,j,k-1)
               write(*,*) 7, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,2,k-1), alocont_nn3d_tmp(i,j,k,7), alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,2,k-1)
               write(*,*) 8, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k-1), alocont_nn3d_tmp(i,j,k,8), alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,2,k-1)
               write(*,*) 9, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,2,k-1), alocont_nn3d_tmp(i,j,k,9), alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(2,2,k-1)
               write(*,*) 10, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k), alocont_nn3d_tmp(i,j,k,10), alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,j-1,k)
               write(*,*) 11, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k), alocont_nn3d_tmp(i,j,k,11), alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)
               write(*,*) 12, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j-1,k), alocont_nn3d_tmp(i,j,k,12), alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(2,j-1,k)
               write(*,*) 13, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k), alocont_nn3d_tmp(i,j,k,13), alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)
               write(*,*) 14, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k), alocont_nn3d_tmp(i,j,k,14), alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)
               write(*,*) 15, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k), alocont_nn3d_tmp(i,j,k,15), alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(2,j,k)
               write(*,*) 16, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,2,k), alocont_nn3d_tmp(i,j,k,16), alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,2,k)
               write(*,*) 17, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k), alocont_nn3d_tmp(i,j,k,17), alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,2,k)
               write(*,*) 18, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,2,k), alocont_nn3d_tmp(i,j,k,18), alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(2,2,k)
               write(*,*) 19, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j-1,k+1), alocont_nn3d_tmp(i,j,k,19), alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,j-1,k+1)
               write(*,*) 20, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j-1,k+1), alocont_nn3d_tmp(i,j,k,20), alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)
               write(*,*) 21, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j-1,k+1), alocont_nn3d_tmp(i,j,k,21), alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(2,j-1,k+1)
               write(*,*) 22, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,j,k+1), alocont_nn3d_tmp(i,j,k,22), alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)
               write(*,*) 23, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,j,k+1), alocont_nn3d_tmp(i,j,k,23), alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)
               write(*,*) 24, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,j,k+1), alocont_nn3d_tmp(i,j,k,24), alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(2,j,k+1)
               write(*,*) 25, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i-1,2,k+1), alocont_nn3d_tmp(i,j,k,25), alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,2,k+1)
               write(*,*) 26, i, j, k, imaskb3d(i,j,k), mint3d_tmp(i,2,k+1), alocont_nn3d_tmp(i,j,k,26), alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,2,k+1)
               write(*,*) 27, i, j, k, imaskb3d(i,j,k), mint3d_tmp(2,2,k+1), alocont_nn3d_tmp(i,j,k,27), alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(2,2,k+1)
               write(*,*)
               
               if(abs(alocont_nn3d_tmp(i,j,k,1)-mint3d_tmp(i-1,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 1'
               if(abs(alocont_nn3d_tmp(i,j,k,2)-mint3d_tmp(i,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 2'
               if(abs(alocont_nn3d_tmp(i,j,k,3)-mint3d_tmp(2,j-1,k-1)).gt.small_number) stop 'error in check_alo2d: 3'
               if(abs(alocont_nn3d_tmp(i,j,k,4)-mint3d_tmp(i-1,j,k-1)).gt.small_number) stop 'error in check_alo2d: 4'
               if(abs(alocont_nn3d_tmp(i,j,k,5)-mint3d_tmp(i,j,k-1)).gt.small_number) stop 'error in check_alo2d: 5'
               if(abs(alocont_nn3d_tmp(i,j,k,6)-mint3d_tmp(2,j,k-1)).gt.small_number) stop 'error in check_alo2d: 6'
               if(abs(alocont_nn3d_tmp(i,j,k,7)-mint3d_tmp(i-1,2,k-1)).gt.small_number) stop 'error in check_alo2d: 7'
               if(abs(alocont_nn3d_tmp(i,j,k,8)-mint3d_tmp(i,2,k-1)).gt.small_number) stop 'error in check_alo2d: 8'
               if(abs(alocont_nn3d_tmp(i,j,k,9)-mint3d_tmp(2,2,k-1)).gt.small_number) stop 'error in check_alo2d: 9'
               if(abs(alocont_nn3d_tmp(i,j,k,10)-mint3d_tmp(i-1,j-1,k)).gt.small_number) stop 'error in check_alo2d: 10'
               if(abs(alocont_nn3d_tmp(i,j,k,11)-mint3d_tmp(i,j-1,k)).gt.small_number) stop 'error in check_alo2d: 11'
               if(abs(alocont_nn3d_tmp(i,j,k,12)-mint3d_tmp(2,j-1,k)).gt.small_number) stop 'error in check_alo2d: 12'
               if(abs(alocont_nn3d_tmp(i,j,k,13)-mint3d_tmp(i-1,j,k)).gt.small_number) stop 'error in check_alo2d: 13'
               if(abs(alocont_nn3d_tmp(i,j,k,14)-mint3d_tmp(i,j,k)).gt.small_number) stop 'error in check_alo2d: 14'
               if(abs(alocont_nn3d_tmp(i,j,k,15)-mint3d_tmp(2,j,k)).gt.small_number) stop 'error in check_alo2d: 15'
               if(abs(alocont_nn3d_tmp(i,j,k,16)-mint3d_tmp(i-1,2,k)).gt.small_number) stop 'error in check_alo2d: 16'
               if(abs(alocont_nn3d_tmp(i,j,k,17)-mint3d_tmp(i,2,k)).gt.small_number) stop 'error in check_alo2d: 17'
               if(abs(alocont_nn3d_tmp(i,j,k,18)-mint3d_tmp(2,2,k)).gt.small_number) stop 'error in check_alo2d: 18'
               if(abs(alocont_nn3d_tmp(i,j,k,19)-mint3d_tmp(i-1,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 19'
               if(abs(alocont_nn3d_tmp(i,j,k,20)-mint3d_tmp(i,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 20'
               if(abs(alocont_nn3d_tmp(i,j,k,21)-mint3d_tmp(2,j-1,k+1)).gt.small_number) stop 'error in check_alo2d: 21'
               if(abs(alocont_nn3d_tmp(i,j,k,22)-mint3d_tmp(i-1,j,k+1)).gt.small_number) stop 'error in check_alo2d: 22'
               if(abs(alocont_nn3d_tmp(i,j,k,23)-mint3d_tmp(i,j,k+1)).gt.small_number) stop 'error in check_alo2d: 23'
               if(abs(alocont_nn3d_tmp(i,j,k,24)-mint3d_tmp(2,j,k+1)).gt.small_number) stop 'error in check_alo2d: 24'
               if(abs(alocont_nn3d_tmp(i,j,k,25)-mint3d_tmp(i-1,2,k+1)).gt.small_number) stop 'error in check_alo2d: 25'
               if(abs(alocont_nn3d_tmp(i,j,k,26)-mint3d_tmp(i,2,k+1)).gt.small_number) stop 'error in check_alo2d: 26'
               if(abs(alocont_nn3d_tmp(i,j,k,27)-mint3d_tmp(2,2,k+1)).gt.small_number) stop 'error in check_alo2d: 27'
            case default
         end select
         
      enddo
   enddo
enddo
!
!
stop 'go on in check_alo3d'
!
!
end subroutine check_alo3d

