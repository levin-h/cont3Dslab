!
!***********************************************************************
!***********************************************************************
!
!        SOME ROUTINES USED FOR 2D CONTINUUM AND LINE TRANSPORT
!
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine set_boundary2d(nue,nn_x,nn_z)
!
use prog_type
use fund_const 
use dime2d, only: nx, nz, x, z, imask2d, int2d, bnue2d, alocont_o_nn2d, scont2d
!
implicit none
!
! ... arguments
real(dp), intent(in) :: nue, nn_x, nn_z
!
! ... local scalars
integer(i4b) :: i, k
!
! ... local characters
!
!for debugging
real(dp), dimension(51) :: int1d_plus, int1d_minus, scont1d
!
int2d=zero
alocont_o_nn2d=zero
!
do i=1, nx
   int2d(i,1)=bnue2d(i,1)
   int2d(i,2)=bnue2d(i,2)
enddo
!
!
!DEBUG (for eps=1.d-3)
!write(*,*) 'debug'
!do i=1, nx
!   int2d(i,1)=one
!   int2d(i,2)=one
!enddo
!int2d=zero
!
!int1d_plus=(/6.2055404011444005d-3, 6.2055404011444005d-3, 6.2055404011423397d-3,  &
!6.2055404011443112d-3, 6.2055404011439738d-3, 6.2055404011436138d-3,  &
!6.2055404011441412d-3, 6.2055404011433744d-3, 6.2055404011414289d-3,  &
!6.2055404010916675d-3, 6.2055403982367602d-3, 6.2055402973796303d-3,  &
!6.2055379805864615d-3, 6.2055030209482694d-3, 6.2051512766667979d-3,  &
!6.2027420358412343d-3, 6.1911967625580665d-3, 6.1511620371428053d-3,  &
!6.0467711056277230d-3, 5.8336610217539489d-3, 5.4794162395446584d-3,  &
!4.9824705233615565d-3, 4.3755491377665424d-3, 3.7131467029983601d-3,  &
!3.0530784163283319d-3, 2.4420639139706448d-3, 1.9093371640856876d-3,  &
!1.4670875616556950d-3, 1.1145802459983506d-3, 8.4328041093139998d-4,  &
!6.4123487023452612d-4, 4.9561460096960808d-4, 3.9394096568582933d-4,  &
!3.2489786386051116d-4, 2.7900364791260844d-4, 2.4894681645690680d-4,  &
!2.2945419944642159d-4, 2.1689241726356240d-4, 2.0882977170174933d-4,  &
!2.0366815444863767d-4, 2.0036915557748112d-4, 1.9826283095735568d-4,  &
!1.9691889116591699d-4, 1.9606175370491939d-4, 1.9551523627495010d-4,  &
!1.9516683243421515d-4, 1.9494474995377376d-4, 1.9480319816446459d-4,  &
!1.9471297934948203d-4, 0.d0, 0.d0/)
!
!
!int1d_minus=(/6.2055404011444005d-3, 6.2055404011444005d-3, 6.2055404011423397d-3, &
!6.2055404011443121d-3, 6.2055404011439746d-3, 6.2055404011436130d-3, &
!6.2055404011441403d-3, 6.2055404011433736d-3, 6.2055404011414281d-3, &
!6.2055404010915591d-3, 6.2055403982307277d-3, 6.2055402971601227d-3, &
!6.2055379753161227d-3, 6.2055029360644142d-3, 6.2051503389847955d-3, &
!6.2027347281956103d-3, 6.1911552217874012d-3, 6.1509833749334034d-3, &
!6.0461675873517877d-3, 5.8320010616479592d-3, 5.4755733488954664d-3, &
!4.9747586197114212d-3, 4.3617881602751743d-3, 3.6908395550058516d-3, &
!3.0196367861604206d-3, 2.3950171832003843d-3, 1.8464830691501278d-3, &
!1.3865732447399753d-3, 1.0149346915656779d-3, 7.2359387527619909d-4, &
!5.0197744911162801d-4, 3.3944729098410301d-4, 2.2491323536243717d-4, &
!1.4691218920830283d-4, 9.5064930196600081d-5, 6.1139302225244192d-5, &
!3.9157919019417402d-5, 2.5002622031345590d-5, 1.5921941496846711d-5, &
!1.0110656892045707d-5, 6.3973049777276812d-6, 4.0267916948905593d-6, &
!2.5144370143423697d-6, 1.5499497496229929d-6, 9.3501051446542064d-7, &
!5.4299819915566460d-7, 2.9312263476773896d-7, 1.3385768819840271d-7, &
!3.2349988520750620d-8, 0.d0, 0.d0 /)
!
!scont1d=(/6.2055404011444005d-3, 6.2055404011444005d-3, 6.2055404011423397d-3,  &
!6.2055404011443112d-3, 6.2055404011439746d-3, 6.2055404011436138d-3,  &
!6.2055404011441403d-3, 6.2055404011433744d-3, 6.2055404011414298d-3,  &
!6.2055404010916658d-3, 6.2055403982366544d-3, 6.2055402973737496d-3,  &
!6.2055379803744844d-3, 6.2055030159289789d-3, 6.2051511974191157d-3,  &
!6.2027411840375477d-3, 6.1911903565817061d-3, 6.1511271737332129d-3,  &
!6.0466284175444193d-3, 5.8332037510604364d-3, 5.4782228398271094d-3,  &
!4.9798414973664263d-3, 4.3705055207737542d-3, 3.7044966762758780d-3,  &
!3.0395267840473048d-3, 2.4223275484433647d-3, 1.8822377469099550d-3,  &
!1.4316091131297013d-3, 1.0698982495058332d-3, 7.8885922654363557d-4,  &
!5.7724001423403119d-4, 4.2331876385169457d-4, 3.1532287385390638d-4,  &
!2.4187416271903095d-4, 1.9305215602209472d-4, 1.6109281519680369d-4,  &
!1.4037648714531158d-4, 1.2703126868179041d-4, 1.1846815693331636d-4,  &
!1.1298718163711620d-4, 1.0948450666256632d-4, 1.0724832300570374d-4,  &
!1.0582160217175601d-4, 1.0491169961597317d-4, 1.0433155141643713d-4,  &
!1.0396171318605180d-4, 1.0372596992078429d-4, 1.0357571182309240d-4,  &
!1.0347994400883367d-4, 0.d0, 0.d0/)
!
!if(n_z.gt.0.) then
!   int2d(1,:)=int1d_plus
!   int2d(2,:)=int1d_plus
!   int2d(nx-1,:)=int1d_plus
!   int2d(nx,:)=int1d_plus
!else
!   int2d(1,:)=int1d_minus
!   int2d(2,:)=int1d_minus
!   int2d(nx-1,:)=int1d_minus
!   int2d(nx,:)=int1d_minus
!endif
!scont2d(1,:)=scont1d
!scont2d(2,:)=scont1d
!scont2d(nx-1,:)=scont1d
!scont2d(nx,:)=scont1d
!
!
end subroutine set_boundary2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine set_boundaryp2d_linb(kk,oindx,bindxa,bindxb,startx,alpha)
!
!          sets the intensities on int2d to values stored in
!           intbound2d on a z-level k (boundary intensities)
!
use prog_type
use dime2d, only: int2d, intbound2d
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: kk,oindx,bindxa,bindxb,startx,alpha
!
int2d(startx-2*alpha,kk)=intbound2d(bindxa,kk,oindx)
int2d(startx-alpha,kk)=intbound2d(bindxb,kk,oindx)
!
end subroutine set_boundaryp2d_linb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine calc_boundaryp2d_linb(startx,endx,alpha,gamma,kk,oindx,bindxa,bindxb,nn_x,nn_z,wall,q14)
!
!           calculates intensity at left and right boundary
!         for periodic boundary conditions using the SC method  
!
use prog_type
use dime2d, only: int2d, opac2d, scont2d, x, z, mint2d_tmp, normalization2d_tmp, alocont_o_nn2d, alocont_nn2d_tmp, intbound2d
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: startx, endx, alpha, gamma, kk, oindx, bindxa, bindxb, q14
real(dp), intent(in) :: nn_x, nn_z, wall
!
! ... local scalars
integer(i4b) :: ii, iim1, iip1, ii_op, kkm1, kkp1
real(dp) :: x_u, x_p, x_d, z_u, z_p, z_d
real(dp) :: scont_u, scont_p, scont_d
real(dp) :: opac_u, opac_p, opac_d
real(dp) :: alo_u, alo_p, alo_pp, alo_d
real(dp) :: int_u
real(dp) :: dels_xu, dels_zu, dels_u, dels_xd, dels_zd, dels_d
real(dp) :: abs_sc, contr_sc, int_sc
real(dp) :: abs_sct, contr_sct
real(dp) :: c03_intu, c04_intu, c05_intu, &
            c03_scontu, c04_scontu, c05_scontu, &
            c02_scontd, c04_scontd, c05_scontd
!
! ... local logicals
!
! ... local characters
!
!
!----------------------directly adjacent ghost point---------------------
!
!set indices
ii=endx+alpha
iim1=ii-alpha
iip1=ii+alpha
!
kkm1=kk-gamma
kkp1=kk+gamma
!
!index on opposite boundary
ii_op=startx-2*alpha
!
!dummy values
scont_d=zero
opac_d=zero
dels_d=zero
!
!calculate distance to intersection point with upwind x-grid-coordinate
dels_xu=(x(ii)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
dels_zu=(z(kk)-z(kkm1))/nn_z
dels_u=min(dels_xu,dels_zu)
!
!calculate distance to intersection point with downwind x-grid-coordinate
dels_xd=(x(iip1)-x(ii))/nn_x
!calculate distance to instersection point with downwind z-grid-coordinate
dels_zd=(z(kkp1)-z(kk))/nn_z
dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
opac_p=opac2d(ii,kk)
scont_p=scont2d(ii,kk)
!
!----------------------------upwind point-------------------------------
!
if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
   x_u = x(ii) - dels_u*nn_x
   z_u = z(kk) - dels_u*nn_z
 
   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(iim1,kk), &
                          scont2d(iim1,kkm1), scont2d(iim1,kk), &
                          int2d(iim1,kkm1), int2d(iim1,kk), &
                          z(kkm1), z(kk), z_u, &
                          c03_scontu, c05_scontu, &
                          c03_intu, c05_intu, &
                          opac_u, scont_u, int_u)
!
elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
   x_u=x(ii)-dels_u*nn_x
   z_u=z(kk)-dels_u*nn_z
!
   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(ii,kkm1), &
                          scont2d(iim1,kkm1), scont2d(ii,kkm1), &
                          int2d(iim1,kkm1), int2d(ii,kkm1), &
                          x(iim1), x(ii), x_u, &
                          c03_scontu, c04_scontu, &
                          c03_intu, c04_intu, &
                          opac_u, scont_u, int_u)              
else
   write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
   stop 'error in set_boundaryp2d_linb: invalid dels_u'
endif
!
!--------------------------downwind point-------------------------------
!
if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd_lin(opac2d(iip1,kk), opac2d(iip1,kkp1), &
                          scont2d(iip1,kk), scont2d(iip1,kkp1), &
                          z(kk), z(kkp1), z_d, &
                          c02_scontd, c05_scontd, &
                          opac_d, scont_d)
   c04_scontd=zero
elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd_lin(opac2d(ii,kkp1), opac2d(iip1,kkp1), &
                          scont2d(ii,kkp1), scont2d(iip1,kkp1), &
                          x(ii), x(iip1), x_d, &
                          c04_scontd, c05_scontd, &
                          opac_d, scont_d)
   c02_scontd=zero
else
   write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
   stop 'error in set_boundaryp2d_linb: invalid dels_d'
endif
!
!-----------------------radiative transfer------------------------------
!
call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
intbound2d(bindxa,kk,oindx)=int_sc
int2d(ii,kk)=int_sc
alocont_o_nn2d(ii,kk,q14)=alo_p
!
!set to opposite side
int2d(ii_op,kk)=int_sc
alocont_o_nn2d(ii_op,kk,14)=alo_p
!
!perform angular integration
mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int_sc*wall
normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alo_p
!
!perform angular integration on opposite side
mint2d_tmp(ii_op,kk) = mint2d_tmp(ii_op,kk) + int_sc*wall
normalization2d_tmp(ii_op,kk) = normalization2d_tmp(ii_op,kk) + wall
alocont_nn2d_tmp(ii_op,kk,q14) = alocont_nn2d_tmp(ii_op,kk,q14) + wall*alo_p
!
!--------------------------next ghost point------------------------------
!
ii=ii_op+alpha
iim1=ii-alpha
iip1=ii+alpha
ii_op=endx+2*alpha
!
!calculate distance to intersection point with upwind x-grid-coordinate
dels_xu=(x(ii)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
dels_zu=(z(kk)-z(kkm1))/nn_z
dels_u=min(dels_xu,dels_zu)
!
!calculate distance to intersection point with downwind x-grid-coordinate
dels_xd=(x(iip1)-x(ii))/nn_x
!calculate distance to instersection point with downwind z-grid-coordinate
dels_zd=(z(kkp1)-z(kk))/nn_z
dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
opac_p=opac2d(ii,kk)
scont_p=scont2d(ii,kk)
!
!----------------------------upwind point-------------------------------
!
if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
   x_u = x(ii) - dels_u*nn_x
   z_u = z(kk) - dels_u*nn_z
 
   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(iim1,kk), &
                          scont2d(iim1,kkm1), scont2d(iim1,kk), &
                          int2d(iim1,kkm1), int2d(iim1,kk), &
                          z(kkm1), z(kk), z_u, &
                          c03_scontu, c05_scontu, &
                          c03_intu, c05_intu, &
                          opac_u, scont_u, int_u)
!
elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
   x_u=x(ii)-dels_u*nn_x
   z_u=z(kk)-dels_u*nn_z
!
   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(ii,kkm1), &
                          scont2d(iim1,kkm1), scont2d(ii,kkm1), &
                          int2d(iim1,kkm1), int2d(ii,kkm1), &
                          x(iim1), x(ii), x_u, &
                          c03_scontu, c04_scontu, &
                          c03_intu, c04_intu, &
                          opac_u, scont_u, int_u)              
else
   write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
   stop 'error in set_boundaryp2d_linb: invalid dels_u'
endif
!
!--------------------------downwind point-------------------------------
!
if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd_lin(opac2d(iip1,kk), opac2d(iip1,kkp1), &
                          scont2d(iip1,kk), scont2d(iip1,kkp1), &
                          z(kk), z(kkp1), z_d, &
                          c02_scontd, c05_scontd, &
                          opac_d, scont_d)
   c04_scontd=zero
elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd_lin(opac2d(ii,kkp1), opac2d(iip1,kkp1), &
                          scont2d(ii,kkp1), scont2d(iip1,kkp1), &
                          x(ii), x(iip1), x_d, &
                          c04_scontd, c05_scontd, &
                          opac_d, scont_d)
   c02_scontd=zero
else
   write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
   stop 'error in set_boundaryp2d_linb: invalid dels_u'
endif
!
!-----------------------radiative transfer------------------------------
!
call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
intbound2d(bindxb,kk,oindx)=int_sc
int2d(ii,kk)=int_sc
alocont_o_nn2d(ii,kk,q14)=alo_p
!
!set to opposite side
int2d(ii_op,kk)=int_sc
alocont_o_nn2d(ii_op,kk,14)=alo_p
!
!perform angular integration
mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int_sc*wall
alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alo_p
normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
!
!perform angular integration on opposite side
mint2d_tmp(ii_op,kk) = mint2d_tmp(ii_op,kk) + int_sc*wall
normalization2d_tmp(ii_op,kk) = normalization2d_tmp(ii_op,kk) + wall
alocont_nn2d_tmp(ii_op,kk,q14) = alocont_nn2d_tmp(ii_op,kk,q14) + wall*alo_p
!
!
!
end subroutine calc_boundaryp2d_linb
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine set_boundaryp2d_lin(nn_x,nn_z,startxb,endxb,alphab,gammab,kk,bindxa,bindxb,oindx,wall,q14,q15,q23,q24)
!
!           calculates intensity at left and right boundary
!         for periodic boundary conditions using the LC method  
!
use prog_type
use dime2d, only: int2d, opac2d, scont2d, x, z, mint2d_tmp, normalization2d_tmp, alocont_o_nn2d, alocont_nn2d_tmp, intbound2d, imask2d
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: startxb,endxb,alphab,gammab,kk,bindxa,bindxb,oindx
integer, intent(in) :: q14, q15, q23, q24
real(dp), intent(in) :: nn_x, nn_z, wall
!
! ... local scalars
integer(i4b) :: i, j, ii, iim1, kkm1, iip1, kkp1, nslab
real(dp) :: x_u, x_p, x_d, z_u, z_p, z_d
real(dp) :: scont_u, scont_p, scont_d
real(dp) :: opac_u, opac_p, opac_d
real(dp) :: alo_u, alo_p, alo_d
real(dp) :: alo_uu, alo_pp, alo_pp2, abs_scc, cc03_scontu, cc04_scontu, cc05_scontu, cc03_intu, cc04_intu, cc05_intu
real(dp) :: int_u
real(dp) :: dels_xu, dels_zu, dels_u, dels_xd, dels_zd, dels_d
real(dp) :: abs_sc, contr_sc, int_sc
real(dp) :: abs_sct, contr_sct
real(dp) :: c03_intu, c04_intu, c05_intu, &
            c03_scontu, c04_scontu, c05_scontu, &
            c02_scontd, c04_scontd, c05_scontd
real(dp) :: a2, b2, c2, iu2, su2, sp2, a1, b1, c1, iu1, su1, sp1
!
! ... local logicals
logical :: lbound
!
! ... local characters
!
!-------------------------if zero boundary condition--------------------
!
!ii=startxb
!!
!int2d(ii,kk)=zero
!int2d(endxb+alphab,kk)=zero
!alocont_o_nn2d(ii,kk,q14) = zero
!alocont_o_nn2d(endxb+alphab,kk,q14) = zero
!mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int2d(ii,kk)*wall
!mint2d_tmp(endxb+alphab,kk) = mint2d_tmp(endxb+alphab,kk) + int2d(endxb+alphab,kk)*wall
!normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
!normalization2d_tmp(endxb+alphab,kk) = normalization2d_tmp(endxb+alphab,kk) + wall
!alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alocont_o_nn2d(ii,kk,q14)
!alocont_nn2d_tmp(endxb+alphab,kk,q14) = alocont_nn2d_tmp(endxb+alphab,kk,q14) + wall*alocont_o_nn2d(endxb+alphab,kk,q14)
!!
!ii=startxb-alphab
!int2d(ii,kk) = int_sc
!int2d(endxb,kk) = int_sc
!alocont_o_nn2d(ii,kk,q14) = zero
!alocont_o_nn2d(endxb,kk,q14) = zero
!mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int2d(ii,kk)*wall
!mint2d_tmp(endxb,kk) = mint2d_tmp(endxb,kk) + int2d(endxb,kk)*wall
!normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
!normalization2d_tmp(endxb,kk) = normalization2d_tmp(endxb,kk) + wall
!alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alocont_o_nn2d(ii,kk,q14)
!alocont_nn2d_tmp(endxb,kk,q14) = alocont_nn2d_tmp(endxb,kk,q14) + wall*alocont_o_nn2d(endxb,kk,q14)
!
!return
!
!-------------first step: calculate intensities at ghost point----------
!-----------------directly adjacent to computational region-------------
!---------------(point G2r for n_x>0 and point G1l for n_x < 0)---------
!--------------------------using LC method------------------------------
!
!G2l----------G1l--------x--------------x----------x-------G2r--------G1r
!
!
!set index where boundary intensity shall be calculated
ii=startxb
iip1=ii-alphab
!
!define total abosorption and source contribution
abs_sct=one
contr_sct=zero
!
!define logical to check if previous z-layer is hit
lbound=.false.
!
!calculate number of required calculations to sweep trough complete slab
nslab = floor((z(kk)-z(kk+gammab))/nn_z/abs(x(ii)-x(endxb+alphab)))+1


!stop
!
!define current z-coordinate and previous z-index
z_p = z(kk)
kkm1=kk+gammab
kkp1=kk-gammab
!
!--------------------current point--------------------------------------
!
opac_p=opac2d(ii,kk)
scont_p=scont2d(ii,kk)
!
!-----------------downwind point (required for dtau steps)--------------
!
dels_xd=(x(iip1)-x(ii))/nn_x
dels_zd=(z(kkp1)-z(kk))/nn_z
dels_d=min(dels_xd,dels_zd)
!
if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd_lin(opac2d(iip1,kk), opac2d(iip1,kkp1), &
                          scont2d(iip1,kk), scont2d(iip1,kkp1), &
                          z(kk), z(kkp1), z_d, &
                          c02_scontd, c05_scontd, &
                          opac_d, scont_d)
   c04_scontd=zero
elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd_lin(opac2d(ii,kkp1), opac2d(iip1,kkp1), &
                          scont2d(ii,kkp1), scont2d(iip1,kkp1), &
                          x(ii), x(iip1), x_d, &
                          c04_scontd, c05_scontd, &
                          opac_d, scont_d)
   c02_scontd=zero
else
   write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
   stop 'error in set_boundaryp2d_lin: invalid dels_d'
endif
!
!-----------------------------------------------------------------------
!
!sweep trough the number of slabs
do j=1, nslab
!
!reset the x_p coordinate to beginnning of the slab
   x_p=x(ii)
!   write(*,*) j
!
!sweep through one slab
   do i=startxb, endxb, alphab
      
      iim1=i+alphab
!
!calculate distance to intersection point with upwind x-grid-coordinate
      dels_xu=(x_p-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
      dels_zu=(z_p-z(kkm1))/nn_z
      dels_u=min(dels_xu,dels_zu)
!      write(*,*) nslab, i, dels_xu, dels_zu, dels_u
!      write(*,*) dels_xu, dels_zu, nn_z, opac_p
!
!----------------------------local point--------------------------------
!
!has been set at beginning of complete procedure, and will be updated
!using the upwind point
!
!----------------------------upwind point-------------------------------
!
      if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
         x_u = x_p - dels_u*nn_x
         z_u = z_p - dels_u*nn_z
!               
         call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(iim1,kk), &
                                scont2d(iim1,kkm1), scont2d(iim1,kk), &
                                int2d(iim1,kkm1), int2d(iim1,kk), &
                                z(kkm1), z(kk), z_u, &
                                c03_scontu, c05_scontu, &
                                c03_intu, c05_intu, &
                                opac_u, scont_u, int_u)
         c04_scontu=zero
         c04_intu=zero
!calculate absorption and source contribution
         call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                           dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!add up source contribution and absorption (first source-contribution, then absorption!!!!)
         contr_sct = contr_sct + abs_sct*contr_sc
         abs_sct = abs_sct*abs_sc
!set local alo at first intersection point
         if(j.eq.1) then
            if(i.eq.startxb) then
               alo_pp=alo_p
               alo_uu=alo_u
               abs_scc=abs_sc
               cc03_scontu=c03_scontu
               cc04_scontu=c04_scontu
               cc05_scontu=c05_scontu
               cc03_intu=c03_intu
               cc04_intu=c04_intu
               cc05_intu=zero       !note: within LC method, no interpolaton of upwind intensities if previous x-level is hit
               a1=abs_scc
               b1=alo_uu
               su1=cc03_scontu
            endif
            if(i.eq.startxb+alphab) then
               alo_pp2=alo_p
               iu2=zero
               a2=abs_sc
               b2=alo_u
               c2=alo_p
               su2=zero
               sp2=su1
               iu1=iu2*a2+b2*su2+c2*sp2
            endif            
         endif
!update the downwind and current point for next step
         x_d=x_p
         z_d=z_p
         opac_d=opac_d
         scont_d=scont_d
         x_p=x_u
         z_p=z_u
         opac_p=opac_u
         scont_p=scont_u
         
      elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
         x_u=x_p - dels_u*nn_x
         z_u=z_p - dels_u*nn_z
!
         call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(i,kkm1), &
                                scont2d(iim1,kkm1), scont2d(i,kkm1), &
                                int2d(iim1,kkm1), int2d(i,kkm1), &
                                x(iim1), x(i), x_u, &
                                c03_scontu, c04_scontu, &
                                c03_intu, c04_intu, &
                                opac_u, scont_u, int_u)
         c05_scontu=zero
         c05_intu=zero
!calculate absorption and source contribution
         call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                           dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!add upt source contribution and absorption (first source-contribution, then absorption!!!!)
         contr_sct = contr_sct + abs_sct*contr_sc
         abs_sct = abs_sct*abs_sc
!set local alo at first intersection point
         if(j.eq.1) then
            if(i.eq.startxb) then
               alo_pp=alo_p
               alo_uu=alo_u
               abs_scc=abs_sc
               cc03_scontu=c03_scontu
               cc04_scontu=c04_scontu
               cc05_scontu=c05_scontu
               cc03_intu=c03_intu
               cc04_intu=c04_intu
               cc05_intu=c05_intu          !if previous z-level is hit, cc05_intu=0 anyways!!
               iu1=cc03_intu*alocont_o_nn2d(ii+alphab,kk+gammab,q14) + &
                   cc04_intu*alocont_o_nn2d(ii+alphab,kk+gammab,q15)
               a1=abs_scc
               b1=alo_uu
               su1=cc03_scontu
            endif
            if(i.eq.startxb+alphab) then
               iu2=c04_intu*alocont_o_nn2d(ii+alphab,kk+gammab,q14)
               alo_pp2=alo_p
               a2=abs_sc
               b2=alo_u
               c2=alo_p
               su2=c04_scontu
               sp2=su1
               iu1=a2*iu2+b2*su2+c2*sp2
            endif
         endif
!set logical and exit loop (note: outer loop is automatically at nslab)
         lbound=.true.
         exit
!
      else
         write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
         stop 'error in set_boundaryp2d_lin: invalid dels_u'
      endif
   enddo
!
enddo
!
!set the intensity + alo
int2d(ii,kk)=int_u*abs_sct + contr_sct
int2d(endxb+alphab,kk)=int2d(ii,kk)
intbound2d(bindxa,kk,oindx)=int2d(ii,kk)

alocont_o_nn2d(ii,kk,q14) = alo_pp
alocont_o_nn2d(ii+alphab,kk,q15) = (alo_uu*cc05_scontu + alo_pp2*cc05_scontu*abs_scc)!*imask2d(ii+alphab,kk)            
alocont_o_nn2d(ii,kk+gammab,q23) = (alo_uu*cc04_scontu + abs_scc*(cc04_intu*alocont_o_nn2d(ii,kk+gammab,q14)))!*imask2d(ii,kk+gammab)
alocont_o_nn2d(ii+alphab,kk+gammab,q24) = alo_uu*cc03_scontu + abs_scc*iu1

!write(*,*) ii, kk, int2d(ii,kk), alocont_o_nn2d(ii+alphab,kk,q15), alo_uu, cc05_scontu, alo_uu*cc05_scontu
!write(*,*) alo_pp2, cc05_scontu, abs_scc!, alo_pp2*cc05_scontu*abs_scc

!
!copy everything on left side (not possible for q15 and q24, because index of previous grid point at left side not defined;
!This point will be explicitly set when setting up the coo-alo
alocont_o_nn2d(endxb+alphab,kk,q14) = alocont_o_nn2d(ii,kk,q14)
alocont_o_nn2d(endxb+alphab,kk+gammab,q23) = alocont_o_nn2d(ii,kk+gammab,q23)
!
!perform angular integration
mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int2d(ii,kk)*wall
mint2d_tmp(endxb+alphab,kk) = mint2d_tmp(endxb+alphab,kk) + int2d(endxb+alphab,kk)*wall
normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
normalization2d_tmp(endxb+alphab,kk) = normalization2d_tmp(endxb+alphab,kk) + wall

alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alocont_o_nn2d(ii,kk,q14)
alocont_nn2d_tmp(ii,kk+gammab,q23) = alocont_nn2d_tmp(ii,kk+gammab,q23) + wall*alocont_o_nn2d(ii,kk+gammab,q23)
alocont_nn2d_tmp(ii+alphab,kk,q15) = alocont_nn2d_tmp(ii+alphab,kk,q15) + wall*alocont_o_nn2d(ii+alphab,kk,q15)
alocont_nn2d_tmp(ii+alphab,kk+gammab,q24) = alocont_nn2d_tmp(ii+alphab,kk+gammab,q24) + wall*alocont_o_nn2d(ii+alphab,kk+gammab,q24)

!again, left side for q15 and q24 not possible
alocont_nn2d_tmp(endxb+alphab,kk,q14) = alocont_nn2d_tmp(endxb+alphab,kk,q14) + wall*alocont_o_nn2d(endxb+alphab,kk,q14)
alocont_nn2d_tmp(endxb+alphab,kk+gammab,q23) = alocont_nn2d_tmp(endxb+alphab,kk+gammab,q23) + wall*alocont_o_nn2d(endxb+alphab,kk+gammab,q23)
!
!------------second step: calculate intensities at ghost point----------
!-----------one grid point further away from computational region-------
!---------------(point G1l for n_x>0 and point G2r for n_x < 0)---------
!--------------------------using SC method------------------------------
!
!G2l----------G1l--------x--------------x----------x-------G2r--------G1r
!
!set indices
ii=endxb
iim1=ii+alphab
iip1=ii-alphab
!
kkm1=kk+gammab
kkp1=kk-gammab
!
!calculate distance to intersection point with upwind x-grid-coordinate
dels_xu=(x(ii)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
dels_zu=(z(kk)-z(kkm1))/nn_z
dels_u=min(dels_xu,dels_zu)

dels_xd=(x(iip1)-x(ii))/nn_x
dels_zd=(z(kkp1)-z(kk))/nn_z
dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
opac_p=opac2d(ii,kk)
scont_p=scont2d(ii,kk)
!
!----------------------------upwind point-------------------------------
!
if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
   x_u = x(ii) - dels_u*nn_x
   z_u = z(kk) - dels_u*nn_z
 
   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(iim1,kk), &
                          scont2d(iim1,kkm1), scont2d(iim1,kk), &
                          int2d(iim1,kkm1), int2d(iim1,kk), &
                          z(kkm1), z(kk), z_u, &
                          c03_scontu, c05_scontu, &
                          c03_intu, c05_intu, &
                          opac_u, scont_u, int_u)
   c04_scontu=zero
   c04_intu=zero
!
elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
   x_u=x(ii)-dels_u*nn_x
   z_u=z(kk)-dels_u*nn_z
!
   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(ii,kkm1), &
                          scont2d(iim1,kkm1), scont2d(ii,kkm1), &
                          int2d(iim1,kkm1), int2d(ii,kkm1), &
                          x(iim1), x(ii), x_u, &
                          c03_scontu, c04_scontu, &
                          c03_intu, c04_intu, &
                          opac_u, scont_u, int_u)
   c05_scontu=zero
   c05_intu=zero
else
   write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
   stop 'error in set_boundaryp2d_lin: invalid dels_u'
endif
!
!---------------------------downwind point------------------------------
!
if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd_lin(opac2d(iip1,kk), opac2d(iip1,kkp1), &
                          scont2d(iip1,kk), scont2d(iip1,kkp1), &
                          z(kk), z(kkp1), z_d, &
                          c02_scontd, c05_scontd, &
                          opac_d, scont_d)
   c04_scontd=zero
elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd_lin(opac2d(ii,kkp1), opac2d(iip1,kkp1), &
                          scont2d(ii,kkp1), scont2d(iip1,kkp1), &
                          x(ii), x(iip1), x_d, &
                          c04_scontd, c05_scontd, &
                          opac_d, scont_d)
   c02_scontd=zero
else
   write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
   stop 'error in set_boundaryp2d_lin: invalid dels_d'
endif
!
!-----------------------radiative transfer------------------------------
!
call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)

!write(*,*) scont_p, scont_u, scont_d
!write(*,'(i5,13es16.8)') ii, x(ii), x_p, opac_u, opac_p, opac_d, dels_u, dels_d, alo_u, alo_p, dels_d, abs_sc
!write(*,*) startxb-alphab, ii
int2d(ii,kk) = int_sc
int2d(startxb-alphab,kk) = int_sc
!write(*,*) iim1, ii, kk, scont2d(iim1,kk), scont2d(1,kk), scont2d(10,kk)
intbound2d(bindxb,kk,oindx)=int2d(ii,kk)
alocont_o_nn2d(ii,kk,q14) = alo_p!*imask2d(ii,kk)
alocont_o_nn2d(ii,kk+gammab,q23) = (alo_u*c04_scontu + abs_sc*(c04_intu*alocont_o_nn2d(ii,kk+gammab,q14)))!*imask2d(ii,kk+gammab)
alocont_o_nn2d(ii+alphab,kk,q15) = (alo_u*c05_scontu + abs_sc*(c05_intu*alocont_o_nn2d(ii+alphab,kk,q14)))!*imask2d(ii+alphab,kk)
alocont_o_nn2d(ii+alphab,kk+gammab,q24) = alo_u*c03_scontu + abs_sc*(c03_intu*alocont_o_nn2d(ii+alphab,kk+gammab,q14) + &
                                                                     c04_intu*alocont_o_nn2d(ii+alphab,kk+gammab,q15) + &
                                                                     c05_intu*alocont_o_nn2d(ii+alphab,kk+gammab,q23))

!
!copy everyhing on right side
alocont_o_nn2d(startxb-alphab,kk,q14) = alocont_o_nn2d(ii,kk,q14)
alocont_o_nn2d(startxb-alphab,kk+gammab,q23) = alocont_o_nn2d(ii,kk+gammab,q23)
alocont_o_nn2d(startxb,kk,q15) = alocont_o_nn2d(ii+alphab,kk,q15)
alocont_o_nn2d(startxb,kk+gammab,q24) = alocont_o_nn2d(ii+alphab,kk+gammab,q24)
!alocont_o_nn2d(startxb,kk,q15) = alocont_o_nn2d(ii+alphab,kk,q15)

!perform angular integration
mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int2d(ii,kk)*wall
mint2d_tmp(startxb-alphab,kk) = mint2d_tmp(startxb-alphab,kk) + int2d(startxb-alphab,kk)*wall

!write(*,*) ii, kk, startxb-alphab, nn_z, int2d(ii, kk), abs_sc, scont_p, alo_p, contr_sc

normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
normalization2d_tmp(startxb-alphab,kk) = normalization2d_tmp(startxb-alphab,kk) + wall

alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alocont_o_nn2d(ii,kk,q14)
alocont_nn2d_tmp(ii,kk+gammab,q23) = alocont_nn2d_tmp(ii,kk+gammab,q23) + wall*alocont_o_nn2d(ii,kk+gammab,q23)
alocont_nn2d_tmp(ii+alphab,kk,q15) = alocont_nn2d_tmp(ii+alphab,kk,q15) + wall*alocont_o_nn2d(ii+alphab,kk,q15)
alocont_nn2d_tmp(ii+alphab,kk+gammab,q24) = alocont_nn2d_tmp(ii+alphab,kk+gammab,q24) + wall*alocont_o_nn2d(ii+alphab,kk+gammab,q24)

!copy everything on right side
alocont_nn2d_tmp(startxb-alphab,kk,q14) = alocont_nn2d_tmp(startxb-alphab,kk,q14) + wall*alocont_o_nn2d(startxb-alphab,kk,q14)
alocont_nn2d_tmp(startxb-alphab,kk+gammab,q23) = alocont_nn2d_tmp(startxb-alphab,kk+gammab,q23) + wall*alocont_o_nn2d(startxb-alphab,kk+gammab,q23)
alocont_nn2d_tmp(startxb,kk,q15) = alocont_nn2d_tmp(startxb,kk,q15) + wall*alocont_o_nn2d(startxb,kk,q15)
alocont_nn2d_tmp(startxb,kk+gammab,q24) = alocont_nn2d_tmp(startxb,kk+gammab,q24) + wall*alocont_o_nn2d(startxb,kk+gammab,q24)
!
!
end subroutine set_boundaryp2d_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine set_boundaryp2d(nn_x,nn_z,startxb,endxb,alphab,gammab,kk,bindxa,bindxb,oindx,wall,q4,q5,q6,q13,q14,q15,q22,q23,q24)
!
!           calculates intensity at left and right boundary
!         for periodic boundary conditions using the LC method  
!            upwind intensities interpolated by using Bezier
!
!only for nn_x > zero
!
use prog_type
use dime2d, only: int2d, opac2d, scont2d, x, z, mint2d_tmp, normalization2d_tmp, alocont_o_nn2d, alocont_nn2d_tmp, intbound2d, imask2d
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: startxb,endxb,alphab,gammab,kk,bindxa,bindxb,oindx
integer, intent(in) :: q4, q5, q6, q13, q14, q15, q22, q23, q24
real(dp), intent(in) :: nn_x, nn_z, wall
!
! ... local scalars
integer(i4b) :: i, j, ii, iim1, kkm1, iip1, kkp1, kkm2, nslab, im1
real(dp) :: x_u, x_p, x_d, z_u, z_p, z_d
real(dp) :: scont_u, scont_p, scont_d
real(dp) :: opac_u, opac_p, opac_d
real(dp) :: alo_u, alo_p, alo_d
real(dp) :: alo_uu, alo_pp, alo_pp2, abs_scc, cc03_scontu, cc04_scontu, cc05_scontu, cc03_intu, cc04_intu, cc05_intu
real(dp) :: int_u
real(dp) :: dels_xu, dels_zu, dels_u, dels_xd, dels_zd, dels_d
real(dp) :: abs_sc, contr_sc, int_sc
real(dp) :: abs_sct, contr_sct
real(dp) :: c01_intu, c02_intu, c03_intu, c04_intu, c05_intu, &
            c01_scontu, c02_scontu, c03_scontu, c04_scontu, c05_scontu, &
            c01_scontd, c02_scontd, c03_scontd, c04_scontd, c05_scontd
real(dp) :: iu2, su2, sp2, iu1, su1, sp1

real(dp) :: iu1_14, iu2_14, int_u14, alo_p14

real(dp) :: iu3_24, iu3_23, iu3_22, iu3_15, iu3_14, iu3_13, iu3_6, iu3_5, iu3_4
real(dp) :: a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3
real(dp) :: c03_scontu1, c04_scontu1, c05_scontu1, c04_scontu2
real(dp) :: c02_intu1, c03_intu1, c04_intu1, &
            c03_intu2, c04_intu2, &
            c04_intu3
!
! ... local logicals
logical :: lbound
!
! ... local characters
!
!if(nn_x.lt. zero) then
!   stop 'error in set_boundaryp2d: nn_x < 0 needs to be implemented'
!endif
!
!-------------first step: calculate intensities at ghost point----------
!-----------------directly adjacent to computational region-------------
!---------------(point G2r for n_x>0 and point G1l for n_x < 0)---------
!--------------------------using LC method------------------------------
!
!G2l----------G1l--------x--------------x----------x-------G2r--------G1r
!
!
!set index where boundary intensity shall be calculated
ii=startxb
iim1=ii+alphab
iip1=ii-alphab
!
!define total abosorption and source contribution
abs_sct=one
contr_sct=zero
!
!define logical to check if previous z-layer is hit
lbound=.false.
!
!calculate number of required calculations to sweep trough complete slab
nslab = floor((z(kk)-z(kk+gammab))/nn_z/abs(x(ii)-x(endxb+alphab)))+1
!
!set alo coefficients
iu3_24=zero
iu3_23=zero
iu3_22=zero
iu3_15=zero
iu3_14=zero
iu3_13=zero
iu3_6=zero
iu3_5=zero
iu3_4=zero
c03_scontu1 = zero
c04_scontu1 = zero
c05_scontu1 = zero
c03_intu1 = zero
c04_intu1 = zero
!
!define current z-coordinate and previous z-index
z_p = z(kk)
kkm1=kk+gammab
kkm2=kkm1+gammab
kkp1=kk-gammab
!
!--------------------current point--------------------------------------
!
opac_p=opac2d(ii,kk)
scont_p=scont2d(ii,kk)
!
!-----------------downwind point (required for dtau steps)--------------
!
dels_xd=(x(iip1)-x(ii))/nn_x
dels_zd=(z(kkp1)-z(kk))/nn_z
dels_d=min(dels_xd,dels_zd)
!
if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd(opac2d(iip1,kkm1), opac2d(iip1,kk), opac2d(iip1,kkp1), &
                      scont2d(iip1,kkm1), scont2d(iip1,kk), scont2d(iip1,kkp1), &
                      z(kkm1), z(kk), z(kkp1), z_d, &
                      c01_scontd, c02_scontd, c05_scontd, &
                      opac_d, scont_d)
   c03_scontd=zero
   c04_scontd=zero
!debug start
!   call coeff1d_contd_lin(opac2d(iip1,kk), opac2d(iip1,kkp1), &
!                          scont2d(iip1,kk), scont2d(iip1,kkp1), &
!                          z(kk), z(kkp1), z_d, &
!                          c02_scontd, c05_scontd, &
!                          opac_d, scont_d)
!   c01_scontd=zero
!   c03_scontd=zero
!   c04_scontd=zero
!debug end   
elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd(opac2d(iim1,kkp1), opac2d(ii,kkp1), opac2d(iip1,kkp1), &
                      scont2d(iim1,kkp1), scont2d(ii,kkp1), scont2d(iip1,kkp1), &
                      x(iim1), x(ii), x(iip1), x_d, &
                      c03_scontd, c04_scontd, c05_scontd, &
                     opac_d, scont_d)
   c01_scontd=zero
   c02_scontd=zero

!debug start
!   call coeff1d_contd_lin(opac2d(ii,kkp1), opac2d(iip1,kkp1), &
!                          scont2d(ii,kkp1), scont2d(iip1,kkp1), &
!                          x(ii), x(iip1), x_d, &
!                          c04_scontd, c05_scontd, &
!                          opac_d, scont_d)
!   c01_scontd=zero
!   c02_scontd=zero
!   c03_scontd=zero
!debug end   
else
   write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
   stop 'error in set_boundaryp2d: invalid dels_d'
endif
!
!-----------------------------------------------------------------------
!
!sweep trough the number of slabs
do j=1, nslab
!
!reset the x_p coordinate to beginnning of the slab
   x_p=x(ii)
!   write(*,*) j
!
!sweep through one slab
   do i=startxb, endxb, alphab
      
      im1=i+alphab
!
!calculate distance to intersection point with upwind x-grid-coordinate
      dels_xu=(x_p-x(im1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
      dels_zu=(z_p-z(kkm1))/nn_z
      dels_u=min(dels_xu,dels_zu)
!      write(*,*) nslab, i, dels_xu, dels_zu, dels_u
!      write(*,*) dels_xu, dels_zu, nn_z, opac_p
!
!----------------------------local point--------------------------------
!
!has been set at beginning of complete procedure, and will be updated
!using the upwind point
!
!----------------------------upwind point-------------------------------
!
      if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
         x_u = x_p - dels_u*nn_x
         z_u = z_p - dels_u*nn_z
!               
         call coeff1d_contu(opac2d(im1,kkm2), opac2d(im1,kkm1), opac2d(im1,kk), &
                            scont2d(im1,kkm2), scont2d(im1,kkm1), scont2d(im1,kk), &
                            int2d(im1,kkm2), int2d(im1,kkm1), int2d(im1,kk), &
                            z(kkm2), z(kkm1), z(kk), z_u, &
                            c01_scontu, c03_scontu, c05_scontu, &
                            c01_intu, c03_intu, c05_intu, &
                            opac_u, scont_u, int_u)
!debug start
!         call coeff1d_contu_lin(opac2d(im1,kkm1), opac2d(im1,kk), &
!                                scont2d(im1,kkm1), scont2d(im1,kk), &
!                                int2d(im1,kkm1), int2d(im1,kk), &
!                                z(kkm1), z(kk), z_u, &
!                                c03_scontu, c05_scontu, &
!                                c03_intu, c05_intu, &
!                                opac_u, scont_u, int_u)
!debug end
!
!store alo coefficients
         if(j.eq.1) then
            if(i.eq.startxb) then !step 1
               c03_scontu1 = c03_scontu
               c05_scontu1 = c05_scontu
            endif
         endif
!         
      elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
         x_u=x_p - dels_u*nn_x
         z_u=z_p - dels_u*nn_z
!
!only linear interpolations, because if left side is hit, then iim2 not given         
         call coeff1d_contu_lin(opac2d(im1,kkm1), opac2d(i,kkm1), &
                                scont2d(im1,kkm1), scont2d(i,kkm1), &
                                int2d(im1,kkm1), int2d(i,kkm1), &
                                x(im1), x(i), x_u, &
                                c03_scontu, c04_scontu, &
                                c03_intu, c04_intu, &
                                opac_u, scont_u, int_u)
!store alo coefficients
         if(j.eq.1) then
            if(i.eq.startxb) then !step 1
               c03_scontu1 = c03_scontu
               c04_scontu1 = c04_scontu
               iu3_24 = c03_intu*alocont_o_nn2d(iim1,kkm1,q14) + &
                            c04_intu*alocont_o_nn2d(iim1,kkm1,q15)
               iu3_23 = c03_intu*alocont_o_nn2d(ii,kkm1,q13) + &
                        c04_intu*alocont_o_nn2d(ii,kkm1,q14)
               iu3_22 = c04_intu*alocont_o_nn2d(iip1,kkm1,q13)
               iu3_15 = c03_intu*alocont_o_nn2d(iim1,kk,q5) + &
                        c04_intu*alocont_o_nn2d(iim1,kk,q6)
               iu3_14 = c03_intu*alocont_o_nn2d(ii,kk,q4) + &
                        c04_intu*alocont_o_nn2d(ii,kk,q5)
               iu3_13 = c04_intu*alocont_o_nn2d(iip1,kk,q4)
            elseif(i.eq.startxb+alphab) then !step 2
               c04_scontu2 = c04_scontu               
               iu3_24 = c03_intu*alocont_o_nn2d(iim1,kkm1,q13) + &
                        c04_intu*alocont_o_nn2d(iim1,kkm1,q14)
               iu3_23 = c04_intu*alocont_o_nn2d(ii,kkm1,q13)
               iu3_15 = c03_intu*alocont_o_nn2d(iim1,kk,q4) + &
                        c04_intu*alocont_o_nn2d(iim1,kk,q5)
               iu3_14 = c04_intu*alocont_o_nn2d(ii,kk,q4)
            elseif(i.eq.startxb+2*alphab) then !step 3
               iu3_24 = c04_intu*alocont_o_nn2d(iim1,kkm1,q13)
               iu3_15 = c04_intu*alocont_o_nn2d(iim1,kk,q4)
            endif
         endif
!         
!set logical and exit loop (note: outer loop is automatically at nslab)
         lbound=.true.
!
      else
         write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
         stop 'error in set_boundaryp2d: invalid dels_u'
      endif

!calculate absorption and source contribution
!      write(*,*) ii, kk, i, abs_sc
!      call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                     dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!      write(*,*) 'fsc_cont check'
      call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                        dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
      alo_d=zero
!add up source contribution and absorption (first source-contribution, then absorption!!!!)
      contr_sct = contr_sct + abs_sct*contr_sc
      abs_sct = abs_sct*abs_sc
!      if(kk.eq.15) then
!         write(*,*) contr_sct, contr_sc, scont_u, scont_p, scont_d
!      endif
!
!store alo coefficients
      if(j.eq.1) then
         if(i.eq.startxb) then !step 1
            a1 = abs_sc
            a2 = one
            a3 = one
            b1 = alo_u
            b2 = zero
            b3 = zero
            c1 = alo_p
            c2 = zero
            c3 = zero
            d1 = alo_d
            d2 = zero
            d3 = zero
         elseif(i.eq.startxb+alphab) then !step 2
            a2 = abs_sc
            b2 = alo_u
            c2 = alo_p
            d2 = alo_d
         elseif(i.eq.startxb+2*alphab) then !step 2
            a3 = abs_sc
            b3 = alo_u
            c3 = alo_p
            d3 = alo_d
         endif
      endif

!      if(i.eq.startxb) write(*,*) 'step 1', scont_u, c03_scontu1, contr_sc, b1*c03_scontu1
!      if(i.eq.startxb+alphab) write(*,*) 'step2', scont_u, c04_scontu2, contr_sc, b2*c04_scontu2+c2*c03_scontu1
!exit the loop if z-layer hit
      if(lbound) exit         
!
!update the downwind and current point for next step
         x_d=x_p
         z_d=z_p
         opac_d=opac_p
         dels_d=dels_u
         scont_d=scont_p
         x_p=x_u
         z_p=z_u
         opac_p=opac_u
         scont_p=scont_u
!      
   enddo
!
enddo
!
!set the intensity + alo
int2d(ii,kk)=int_u*abs_sct + contr_sct
int2d(endxb+alphab,kk)=int2d(ii,kk)
intbound2d(bindxa,kk,oindx)=int2d(ii,kk)
!
!if(kk.eq.20) write(*,*) int2d(ii,kk), int_u, abs_sct
!write(*,*) scont2d(iip1,kkm1), alocont_o_nn2d(ii,kkm1,q23), int2d(ii,kk), int2d(ii,kk)
!if(kk.eq.3) then
!   write(*,*) ii, kk, nn_z, int2d(ii,kk), int_u
!!   write(*,*) endxb+alphab, kk, nn_z, int2d(endxb+alphab,kk), int_u
!endif
!write(*,*) a1, a2, a3


!
alocont_o_nn2d(iim1,kkm1,q24) = a1*a2*a3*iu3_24 + (a1*a2*c3+a1*b2)*c04_scontu2 + (a1*a2*d3+a1*c2+b1)*c03_scontu1
alocont_o_nn2d(ii,kkm1,q23) = a1*a2*a3*iu3_23 + (a1*a2*d3+a1*c2+b1)*c04_scontu1
alocont_o_nn2d(iip1,kkm1,q22) = a1*a2*a3*iu3_22 + d1*c01_scontd
alocont_o_nn2d(iim1,kk,q15) = a1*a2*a3*iu3_15 + (a1*a2*d3+a1*c2+b1)*c05_scontu1
alocont_o_nn2d(ii,kk,q14) = a1*a2*a3*iu3_14 + a1*d2+c1
alocont_o_nn2d(iip1,kk,q13) = a1*a2*a3*iu3_13 + d1*c02_scontd
alocont_o_nn2d(iim1,kkp1,q6) = a1*a2*a3*iu3_6 + d1*c03_scontd
alocont_o_nn2d(ii,kkp1,q5) = a1*a2*a3*iu3_5 + d1*c04_scontd
alocont_o_nn2d(iip1,kkp1,q4) = a1*a2*a3*iu3_4 + d1*c05_scontd


!if(kk.eq.3) then
!   write(*,*) ii, kk, int2d(ii,kk), alocont_o_nn2d(iim1,kk,q15), contr_sct, c05_scontu1, (a1*a2*d3+a1*c2+b1)
!endif

!write(*,*) ii, kk, int2d(ii,kk), scont2d(ii,kk), alocont_o_nn2d(ii,kk,q14), contr_sct, c1, d2, a1

!
!copy everything on left side (not possible for q6, q15, q24, because index of previous grid point at left side not defined;
!This point will be explicitly set when setting up the coo-alo
alocont_o_nn2d(endxb+alphab,kkm1,q23) = alocont_o_nn2d(ii,kkm1,q23)
alocont_o_nn2d(endxb,kkm1,q22) = alocont_o_nn2d(iip1,kkm1,q22)
alocont_o_nn2d(endxb+alphab,kk,q14)=alocont_o_nn2d(ii,kk,q14) 
alocont_o_nn2d(endxb,kk,q13)=alocont_o_nn2d(iip1,kk,q13)
alocont_o_nn2d(endxb+alphab,kkp1,q5)=alocont_o_nn2d(ii,kkp1,q5)
alocont_o_nn2d(endxb,kkp1,q4)=alocont_o_nn2d(iip1,kkp1,q4) 
!
!
!perform angular integration
mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int2d(ii,kk)*wall
mint2d_tmp(endxb+alphab,kk) = mint2d_tmp(endxb+alphab,kk) + int2d(endxb+alphab,kk)*wall
normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
normalization2d_tmp(endxb+alphab,kk) = normalization2d_tmp(endxb+alphab,kk) + wall
!write(*,*) endxb+alphab, kk, normalization2d_tmp(endxb+alphab,kk), wall

alocont_nn2d_tmp(iim1,kkm1,q24) = alocont_nn2d_tmp(iim1,kkm1,q24) + wall*alocont_o_nn2d(iim1,kkm1,q24)
alocont_nn2d_tmp(ii,kkm1,q23) = alocont_nn2d_tmp(ii,kkm1,q23) + wall*alocont_o_nn2d(ii,kkm1,q23)
alocont_nn2d_tmp(iip1,kkm1,q22) = alocont_nn2d_tmp(iip1,kkm1,q23) + wall*alocont_o_nn2d(iip1,kkm1,q23)
alocont_nn2d_tmp(iim1,kk,q15) = alocont_nn2d_tmp(iim1,kk,q15) + wall*alocont_o_nn2d(iim1,kk,q15)
alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alocont_o_nn2d(ii,kk,q14)
alocont_nn2d_tmp(iip1,kk,q13) = alocont_nn2d_tmp(iip1,kk,q13) + wall*alocont_o_nn2d(iip1,kk,q13)
alocont_nn2d_tmp(iim1,kkp1,q6) = alocont_nn2d_tmp(iim1,kkp1,q6) + wall*alocont_o_nn2d(iim1,kkp1,q6)
alocont_nn2d_tmp(ii,kkp1,q5) = alocont_nn2d_tmp(ii,kkp1,q5) + wall*alocont_o_nn2d(ii,kkp1,q5)
alocont_nn2d_tmp(iip1,kkp1,q4) = alocont_nn2d_tmp(iip1,kkp1,q4) + wall*alocont_o_nn2d(iip1,kkp1,q4)

!again, left side for q6, q15 and q24 not possible
alocont_nn2d_tmp(endxb+alphab,kkm1,q23) = alocont_nn2d_tmp(endxb+alphab,kkm1,q23) + wall*alocont_o_nn2d(endxb+alphab,kkm1,q23)
alocont_nn2d_tmp(endxb,kkm1,q22) = alocont_nn2d_tmp(endxb,kkm1,q22) + wall*alocont_o_nn2d(endxb,kkm1,q22)
alocont_nn2d_tmp(endxb+alphab,kk,q14) = alocont_nn2d_tmp(endxb+alphab,kk,q14) + wall*alocont_o_nn2d(endxb+alphab,kk,q14)
alocont_nn2d_tmp(endxb,kk,q13) = alocont_nn2d_tmp(endxb,kk,q13) + wall*alocont_o_nn2d(endxb,kk,q13)
alocont_nn2d_tmp(endxb+alphab,kkp1,q5) = alocont_nn2d_tmp(endxb+alphab,kkp1,q5) + wall*alocont_o_nn2d(endxb+alphab,kkp1,q5)
alocont_nn2d_tmp(endxb,kkp1,q4) = alocont_nn2d_tmp(endxb,kkp1,q4) + wall*alocont_o_nn2d(endxb,kkp1,q4)


!write(*,*) scont2d(iim1,kkm1), alocont_o_nn2d(iim1,kkm1,q24), int2d(ii,kk)
!
!------------second step: calculate intensities at ghost point----------
!-----------one grid point further away from computational region-------
!---------------(point G1l for n_x>0 and point G2r for n_x < 0)---------
!--------------------------using SC method------------------------------
!
!G2l----------G1l--------x--------------x----------x-------G2r--------G1r
!
!set indices
ii=endxb
iim1=ii+alphab
iip1=ii-alphab
!
kkm1=kk+gammab
kkm2=kkm1+gammab
kkp1=kk-gammab
!
!calculate distance to intersection point with upwind x-grid-coordinate
dels_xu=(x(ii)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
dels_zu=(z(kk)-z(kkm1))/nn_z
dels_u=min(dels_xu,dels_zu)

dels_xd=(x(iip1)-x(ii))/nn_x
dels_zd=(z(kkp1)-z(kk))/nn_z
dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
opac_p=opac2d(ii,kk)
scont_p=scont2d(ii,kk)
!
!----------------------------upwind point-------------------------------
!
if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
   x_u = x(ii) - dels_u*nn_x
   z_u = z(kk) - dels_u*nn_z
 
   call coeff1d_contu(opac2d(iim1,kkm2), opac2d(iim1,kkm1), opac2d(iim1,kk), &
                      scont2d(iim1,kkm2), scont2d(iim1,kkm1), scont2d(iim1,kk), &
                      int2d(iim1,kkm2), int2d(iim1,kkm1), int2d(iim1,kk), &
                      z(kkm2), z(kkm1), z(kk), z_u, &
                      c01_scontu, c03_scontu, c05_scontu, &
                      c01_intu, c03_intu, c05_intu, &
                      opac_u, scont_u, int_u)
   c02_scontu=zero
   c04_scontu=zero
   c02_intu=zero
   c04_intu=zero

!debug start
!   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(iim1,kk), &
!                          scont2d(iim1,kkm1), scont2d(iim1,kk), &
!                          int2d(iim1,kkm1), int2d(iim1,kk), &
!                          z(kkm1), z(kk), z_u, &
!                          c03_scontu, c05_scontu, &
!                          c03_intu, c05_intu, &
!                          opac_u, scont_u, int_u)
!   c01_scontu=zero
!   c02_scontu=zero
!   c04_scontu=zero
!   c01_intu=zero
!   c02_intu=zero
!   c04_intu=zero
!debug end                  
!
elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
   x_u=x(ii)-dels_u*nn_x
   z_u=z(kk)-dels_u*nn_z
!
!only linear interpolations, because iim2 is not know (starts at right point again)   
   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(ii,kkm1), &
                          scont2d(iim1,kkm1), scont2d(ii,kkm1), &
                          int2d(iim1,kkm1), int2d(ii,kkm1), &
                          x(iim1), x(ii), x_u, &
                          c03_scontu, c04_scontu, &
                          c03_intu, c04_intu, &
                          opac_u, scont_u, int_u)
   c05_scontu=zero
   c05_intu=zero
else
   write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
   stop 'error in set_boundaryp2d: invalid dels_u'
endif
!
!---------------------------downwind point------------------------------
!
if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd(opac2d(iip1,kkm1), opac2d(iip1,kk), opac2d(iip1,kkp1), &
                      scont2d(iip1,kkm1), scont2d(iip1,kk), scont2d(iip1,kkp1), &
                      z(kkm1), z(kk), z(kkp1), z_d, &
                      c01_scontd, c02_scontd, c05_scontd, &
                      opac_d, scont_d)
   c03_scontd=zero
   c04_scontd=zero

!debug start
!   call coeff1d_contd_lin(opac2d(iip1,kk), opac2d(iip1,kkp1), &
!                          scont2d(iip1,kk), scont2d(iip1,kkp1), &
!                          z(kk), z(kkp1), z_d, &
!                          c02_scontd, c05_scontd, &
!                          opac_d, scont_d)
!   c01_scontd=zero
!   c03_scontd=zero
!   c04_scontd=zero
!debug end
!
elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd(opac2d(iim1,kkp1), opac2d(ii,kkp1), opac2d(iip1,kkp1), &
                      scont2d(iim1,kkp1), scont2d(i,kkp1), scont2d(iip1,kkp1), &
                      x(iim1), x(ii), x(iip1), x_d, &
                      c03_scontd, c04_scontd, c05_scontd, &
                      opac_d, scont_d)
   c01_scontd=zero
   c02_scontd=zero
!
!debug start   
!   call coeff1d_contd_lin(opac2d(ii,kkp1), opac2d(iip1,kkp1), &
!                          scont2d(ii,kkp1), scont2d(iip1,kkp1), &
!                          x(ii), x(iip1), x_d, &
!                          c04_scontd, c05_scontd, &
!                          opac_d, scont_d)
   c01_scontd=zero
   c02_scontd=zero
   c03_scontd=zero
!debug end
else
   write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
   stop 'error in set_boundaryp2d: invalid dels_d'
endif
!
!-----------------------radiative transfer------------------------------
!
!call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!              dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
alo_d=zero
!write(*,*) scont_p, scont_u, scont_d
!write(*,'(i5,13es16.8)') ii, x(ii), x_p, opac_u, opac_p, opac_d, dels_u, dels_d, alo_u, alo_p, dels_d, abs_sc
!write(*,*) startxb-alphab, ii
int2d(ii,kk) = int_sc
int2d(startxb-alphab,kk) = int_sc
intbound2d(bindxb,kk,oindx)=int2d(ii,kk)

!if(kk.eq.4) then
!   write(*,*) ii, kk, nn_z, int2d(ii,kk), int_u
!   write(*,*) startxb-alphab, kk, nn_z, int2d(startxb-alphab,kk), int_u
!endif


alocont_o_nn2d(iip1,kkp1,q4) = alo_d*c05_scontd
alocont_o_nn2d(ii,kkp1,q5) = (alo_d*c04_scontd + abs_sc*(c05_intu*alocont_o_nn2d(ii,kkp1,q4)))
alocont_o_nn2d(iim1,kkp1,q6) = (alo_d*c03_scontd + abs_sc*(c05_intu*alocont_o_nn2d(iim1,kkp1,q5)))
alocont_o_nn2d(iip1,kk,q13) = (alo_d*c02_scontd + abs_sc*(c04_intu*alocont_o_nn2d(iip1,kk,q4)))
alocont_o_nn2d(ii,kk,q14) = (alo_p + abs_sc*((c03_intu*alocont_o_nn2d(ii,kk,q4) + &
                                           c04_intu*alocont_o_nn2d(ii,kk,q5) + &
                                           c05_intu*alocont_o_nn2d(ii,kk,q13))))
alocont_o_nn2d(iim1,kk,q15) = (alo_u*c05_scontu + abs_sc*((c02_intu*alocont_o_nn2d(iim1,kk,q4) + &
                                                                     c03_intu*alocont_o_nn2d(iim1,kk,q5) + &
                                                                     c04_intu*alocont_o_nn2d(iim1,kk,q6) + &
                                                                     c05_intu*alocont_o_nn2d(iim1,kk,q14))))
alocont_o_nn2d(iip1,kkm1,q22) = (alo_d*c01_scontd + abs_sc*((c04_intu*alocont_o_nn2d(iip1,kkm1,q13))))
alocont_o_nn2d(ii,kkm1,q23) = (alo_u*c04_scontu + abs_sc*((c03_intu*alocont_o_nn2d(ii,kkm1,q13) + &
                                                         c04_intu*alocont_o_nn2d(ii,kkm1,q14) + &
                                                         c05_intu*alocont_o_nn2d(ii,kkm1,q22) + &
                                                         c01_intu*alocont_o_nn2d(ii,kkm1,q4))))
!write(*,*) scont2d(ii,kkm1), alocont_o_nn2d(ii,kkm1,q23), int2d(ii,kk), int2d(iim1,kk)
!write(*,*)
alocont_o_nn2d(iim1,kkm1,q24) = (alo_u*c03_scontu + abs_sc*((c02_intu*alocont_o_nn2d(iim1,kkm1,q13) + &
                                                            c03_intu*alocont_o_nn2d(iim1,kkm1,q14) + &
                                                            c04_intu*alocont_o_nn2d(iim1,kkm1,q15) + &
                                                            c05_intu*alocont_o_nn2d(iim1,kkm1,q23) + &
                                                            c01_intu*alocont_o_nn2d(iim1,kkm1,q5))))
!
!copy everyhing on right side (q4, q13, q22 cannot be stored because out of bounds)
alocont_o_nn2d(startxb-alphab,kkp1,q5) = alocont_o_nn2d(ii,kkp1,q5)
alocont_o_nn2d(startxb,kkp1,q6) = alocont_o_nn2d(iim1,kkp1,q6)
alocont_o_nn2d(startxb-alphab,kk,q14) = alocont_o_nn2d(ii,kk,q14)
alocont_o_nn2d(startxb,kk,q15) = alocont_o_nn2d(iim1,kk,q15)
alocont_o_nn2d(startxb-alphab,kkm1,q23) = alocont_o_nn2d(ii,kkm1,q23)
alocont_o_nn2d(startxb,kkm1,q24) = alocont_o_nn2d(iim1,kkm1,q24)
!
!perform angular integration
mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int2d(ii,kk)*wall
mint2d_tmp(startxb-alphab,kk) = mint2d_tmp(startxb-alphab,kk) + int2d(startxb-alphab,kk)*wall
!
normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
normalization2d_tmp(startxb-alphab,kk) = normalization2d_tmp(startxb-alphab,kk) + wall
!
alocont_nn2d_tmp(iip1,kkp1,q4) = alocont_nn2d_tmp(iip1,kkp1,q4) + wall*alocont_o_nn2d(iip1,kkp1,q4)
alocont_nn2d_tmp(ii,kkp1,q5) = alocont_nn2d_tmp(ii,kkp1,q5) + wall*alocont_o_nn2d(ii,kkp1,q5)
alocont_nn2d_tmp(iim1,kkp1,q6) = alocont_nn2d_tmp(iim1,kkp1,q6) + wall*alocont_o_nn2d(iim1,kkp1,q6)
alocont_nn2d_tmp(iip1,kk,q13) = alocont_nn2d_tmp(iip1,kk,q13) + wall*alocont_o_nn2d(iip1,kk,q13)
alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alocont_o_nn2d(ii,kk,q14)
alocont_nn2d_tmp(iim1,kk,q15) = alocont_nn2d_tmp(iim1,kk,q15) + wall*alocont_o_nn2d(iim1,kk,q15)
alocont_nn2d_tmp(iip1,kkm1,q22) = alocont_nn2d_tmp(iip1,kkm1,q22) + wall*alocont_o_nn2d(iip1,kkm1,q22)
alocont_nn2d_tmp(ii,kkm1,q23) = alocont_nn2d_tmp(ii,kkm1,q23) + wall*alocont_o_nn2d(ii,kkm1,q23)
alocont_nn2d_tmp(iim1,kkm1,q24) = alocont_nn2d_tmp(iim1,kkm1,q24) + wall*alocont_o_nn2d(iim1,kkm1,q24)
!
!copy everything on right side
alocont_nn2d_tmp(startxb-alphab,kkp1,q5) = alocont_nn2d_tmp(startxb-alphab,kkp1,q5) + wall*alocont_o_nn2d(startxb-alphab,kkp1,q5)
alocont_nn2d_tmp(startxb,kkp1,q6) = alocont_nn2d_tmp(startxb,kkp1,q6) + wall*alocont_o_nn2d(startxb,kkp1,q6)
alocont_nn2d_tmp(startxb-alphab,kk,q14) = alocont_nn2d_tmp(startxb-alphab,kk,q14) + wall*alocont_o_nn2d(startxb-alphab,kk,q14)
alocont_nn2d_tmp(startxb,kk,q15) = alocont_nn2d_tmp(startxb,kk,q15) + wall*alocont_o_nn2d(startxb,kk,q15)
alocont_nn2d_tmp(startxb-alphab,kkm1,q23) = alocont_nn2d_tmp(startxb-alphab,kkm1,q23) + wall*alocont_o_nn2d(startxb-alphab,kkm1,q23)
alocont_nn2d_tmp(startxb,kkm1,q24) = alocont_nn2d_tmp(startxb,kkm1,q24) + wall*alocont_o_nn2d(startxb,kkm1,q24)

!
!
end subroutine set_boundaryp2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine set_boundaryp2d_v02(nn_x,nn_z,startxb,endxb,alphab,gammab,kk,bindxa,bindxb,oindx,wall,q4,q5,q6,q13,q14,q15,q22,q23,q24)
!
!           calculates intensity at left and right boundary
!         for periodic boundary conditions using the LC method  
!            upwind intensities interpolated by using Bezier
!
!only for nn_x > zero
!
use prog_type
use dime2d, only: int2d, opac2d, scont2d, x, z, mint2d_tmp, normalization2d_tmp, alocont_o_nn2d, alocont_nn2d_tmp, intbound2d, imask2d
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: startxb,endxb,alphab,gammab,kk,bindxa,bindxb,oindx
integer, intent(in) :: q4, q5, q6, q13, q14, q15, q22, q23, q24
real(dp), intent(in) :: nn_x, nn_z, wall
!
! ... local scalars
integer(i4b) :: i, j, ii, iim1, kkm1, iip1, kkp1, kkm2, nslab, im1
real(dp) :: x_u, x_p, x_d, z_u, z_p, z_d
real(dp) :: scont_u, scont_p, scont_d
real(dp) :: opac_u, opac_p, opac_d
real(dp) :: alo_u, alo_p, alo_d
real(dp) :: alo_uu, alo_pp, alo_pp2, abs_scc, cc03_scontu, cc04_scontu, cc05_scontu, cc03_intu, cc04_intu, cc05_intu
real(dp) :: int_u
real(dp) :: dels_xu, dels_zu, dels_u, dels_xd, dels_zd, dels_d, dels_uu
real(dp) :: abs_sc, contr_sc, int_sc
real(dp) :: abs_sct, contr_sct
real(dp) :: c01_intu, c02_intu, c03_intu, c04_intu, c05_intu, &
            c01_scontu, c02_scontu, c03_scontu, c04_scontu, c05_scontu, &
            c01_scontd, c02_scontd, c03_scontd, c04_scontd, c05_scontd
real(dp) :: iu2, su2, sp2, iu1, su1, sp1

real(dp) :: iu1_14, iu2_14, int_u14, alo_p14

real(dp) :: iu3_24, iu3_23, iu3_22, iu3_15, iu3_14, iu3_13, iu3_6, iu3_5, iu3_4
real(dp) :: a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3
real(dp) :: c03_scontu1, c04_scontu1, c05_scontu1, c04_scontu2
real(dp) :: c02_intu1, c03_intu1, c04_intu1, &
            c03_intu2, c04_intu2, &
            c04_intu3
!
! ... local logicals
logical :: lbound
!
! ... local characters
!
!if(nn_x.lt. zero) then
!   stop 'error in set_boundaryp2d: nn_x < 0 needs to be implemented'
!endif
!
!-------------first step: calculate intensities at ghost point----------
!-----------------directly adjacent to computational region-------------
!---------------(point G2r for n_x>0 and point G1l for n_x < 0)---------
!--------------------------using LC method------------------------------
!
!G2l----------G1l--------x--------------x----------x-------G2r--------G1r
!
!
!set index where boundary intensity shall be calculated
ii=startxb
iim1=ii+alphab
iip1=ii-alphab
!
!define total abosorption and source contribution
abs_sct=one
contr_sct=zero
!
!define logical to check if previous z-layer is hit
lbound=.false.
!
!calculate number of required calculations to sweep trough complete slab
nslab = floor((z(kk)-z(kk+gammab))/nn_z/abs(x(ii)-x(endxb+alphab)))+1
!
!set alo coefficients
iu3_24=zero
iu3_23=zero
iu3_22=zero
iu3_15=zero
iu3_14=zero
iu3_13=zero
iu3_6=zero
iu3_5=zero
iu3_4=zero
c03_scontu1 = zero
c04_scontu1 = zero
c05_scontu1 = zero
c03_intu1 = zero
c04_intu1 = zero
!
!define current z-coordinate and previous z-index
z_p = z(kk)
kkm1=kk+gammab
kkm2=kkm1+gammab
kkp1=kk-gammab
!
!--------------------current point--------------------------------------
!
opac_p=opac2d(ii,kk)
scont_p=scont2d(ii,kk)
!
!-----------------downwind point (required for dtau steps)--------------
!
dels_xd=(x(iip1)-x(ii))/nn_x
dels_zd=(z(kkp1)-z(kk))/nn_z
dels_d=min(dels_xd,dels_zd)
!
if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd(opac2d(iip1,kkm1), opac2d(iip1,kk), opac2d(iip1,kkp1), &
                      scont2d(iip1,kkm1), scont2d(iip1,kk), scont2d(iip1,kkp1), &
                      z(kkm1), z(kk), z(kkp1), z_d, &
                      c01_scontd, c02_scontd, c05_scontd, &
                      opac_d, scont_d)
   c03_scontd=zero
   c04_scontd=zero
!debug start
!   call coeff1d_contd_lin(opac2d(iip1,kk), opac2d(iip1,kkp1), &
!                          scont2d(iip1,kk), scont2d(iip1,kkp1), &
!                          z(kk), z(kkp1), z_d, &
!                          c02_scontd, c05_scontd, &
!                          opac_d, scont_d)
!   c01_scontd=zero
!   c03_scontd=zero
!   c04_scontd=zero
!debug end   
elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd(opac2d(iim1,kkp1), opac2d(ii,kkp1), opac2d(iip1,kkp1), &
                      scont2d(iim1,kkp1), scont2d(ii,kkp1), scont2d(iip1,kkp1), &
                      x(iim1), x(ii), x(iip1), x_d, &
                      c03_scontd, c04_scontd, c05_scontd, &
                     opac_d, scont_d)
   c01_scontd=zero
   c02_scontd=zero

!debug start
!   call coeff1d_contd_lin(opac2d(ii,kkp1), opac2d(iip1,kkp1), &
!                          scont2d(ii,kkp1), scont2d(iip1,kkp1), &
!                          x(ii), x(iip1), x_d, &
!                          c04_scontd, c05_scontd, &
!                          opac_d, scont_d)
!   c01_scontd=zero
!   c02_scontd=zero
!   c03_scontd=zero
!debug end   
else
   write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
   stop 'error in set_boundaryp2d: invalid dels_d'
endif
!
!-----------------------------------------------------------------------
!
dels_uu=0.d0
!
!sweep trough the number of slabs
do j=1, nslab
!
!reset the x_p coordinate to beginnning of the slab
   x_p=x(ii)
!   write(*,*) j
!
!sweep through one slab
   do i=startxb, endxb, alphab
      
      im1=i+alphab
!
!calculate distance to intersection point with upwind x-grid-coordinate
      dels_xu=(x_p-x(im1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
      dels_zu=(z_p-z(kkm1))/nn_z
      dels_u=min(dels_xu,dels_zu)
!
      dels_uu=dels_uu+dels_u      
!
!----------------------------upwind point-------------------------------
!
      if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
         x_u = x_p - dels_u*nn_x
         z_u = z_p - dels_u*nn_z
!         
      elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
         x_u=x_p - dels_u*nn_x
         z_u=z_p - dels_u*nn_z
!
!only linear interpolations, because if left side is hit, then iim2 not given         
         call coeff1d_contu_lin(opac2d(im1,kkm1), opac2d(i,kkm1), &
                                scont2d(im1,kkm1), scont2d(i,kkm1), &
                                int2d(im1,kkm1), int2d(i,kkm1), &
                                x(im1), x(i), x_u, &
                                c03_scontu, c04_scontu, &
                                c03_intu, c04_intu, &
                                opac_u, scont_u, int_u)
!store alo coefficients
         if(i.eq.startxb) then !step 1
            c03_scontu1 = c03_scontu
            c04_scontu1 = c04_scontu
            iu3_24 = c03_intu*alocont_o_nn2d(iim1,kkm1,q14) + &
                     c04_intu*alocont_o_nn2d(iim1,kkm1,q15)
            iu3_23 = c03_intu*alocont_o_nn2d(ii,kkm1,q13) + &
                     c04_intu*alocont_o_nn2d(ii,kkm1,q14)
            iu3_22 = c04_intu*alocont_o_nn2d(iip1,kkm1,q13)
            iu3_15 = c03_intu*alocont_o_nn2d(iim1,kk,q5) + &
                     c04_intu*alocont_o_nn2d(iim1,kk,q6)
            iu3_14 = c03_intu*alocont_o_nn2d(ii,kk,q4) + &
                     c04_intu*alocont_o_nn2d(ii,kk,q5)
            iu3_13 = c04_intu*alocont_o_nn2d(iip1,kk,q4)
         elseif(i.eq.startxb+alphab) then !step 2
            c04_scontu2 = c04_scontu               
            iu3_24 = c03_intu*alocont_o_nn2d(iim1,kkm1,q13) + &
                     c04_intu*alocont_o_nn2d(iim1,kkm1,q14)
            iu3_23 = c04_intu*alocont_o_nn2d(ii,kkm1,q13)
            iu3_15 = c03_intu*alocont_o_nn2d(iim1,kk,q4) + &
                     c04_intu*alocont_o_nn2d(iim1,kk,q5)
            iu3_14 = c04_intu*alocont_o_nn2d(ii,kk,q4)
         endif
!         
!set logical and exit loop (note: outer loop is automatically at nslab)
         lbound=.true.
!
!calculate absorption and source contribution
!         call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                       dels_uu, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
         call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                           dels_uu, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
         alo_d=zero
!
         int2d(ii,kk)=int_u*abs_sc + contr_sc
         int2d(endxb+alphab,kk)=int2d(ii,kk)
         intbound2d(bindxa,kk,oindx)=int2d(ii,kk)
!
         alocont_o_nn2d(iim1,kkm1,q24) = abs_sc*iu3_24 + (c03_scontu1+c04_scontu2)*alo_u
         alocont_o_nn2d(ii,kkm1,q23) = abs_sc*iu3_23 + c04_scontu1*alo_u
         alocont_o_nn2d(iip1,kkm1,q22) = abs_sc*iu3_22 + alo_d*c01_scontd
         alocont_o_nn2d(iim1,kk,q15) = abs_sc*iu3_15
         alocont_o_nn2d(ii,kk,q14) = abs_sc*iu3_14 + alo_p
         alocont_o_nn2d(iip1,kk,q13) = abs_sc*iu3_13 + alo_d*c02_scontd
         alocont_o_nn2d(iim1,kkp1,q6) = abs_sc*iu3_6 + alo_d*c03_scontd
         alocont_o_nn2d(ii,kkp1,q5) = abs_sc*iu3_5 + alo_d*c04_scontd
         alocont_o_nn2d(iip1,kkp1,q4) = abs_sc*iu3_4 + alo_d*c05_scontd
!
!copy everything on left side (not possible for q6, q15, q24, because index of previous grid point at left side not defined;
!This point will be explicitly set when setting up the coo-alo
         alocont_o_nn2d(endxb+alphab,kkm1,q23) = alocont_o_nn2d(ii,kkm1,q23)
         alocont_o_nn2d(endxb,kkm1,q22) = alocont_o_nn2d(iip1,kkm1,q22)
         alocont_o_nn2d(endxb+alphab,kk,q14)=alocont_o_nn2d(ii,kk,q14) 
         alocont_o_nn2d(endxb,kk,q13)=alocont_o_nn2d(iip1,kk,q13)
         alocont_o_nn2d(endxb+alphab,kkp1,q5)=alocont_o_nn2d(ii,kkp1,q5)
         alocont_o_nn2d(endxb,kkp1,q4)=alocont_o_nn2d(iip1,kkp1,q4) 

!perform angular integration
         mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int2d(ii,kk)*wall
         mint2d_tmp(endxb+alphab,kk) = mint2d_tmp(endxb+alphab,kk) + int2d(endxb+alphab,kk)*wall
         normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
         normalization2d_tmp(endxb+alphab,kk) = normalization2d_tmp(endxb+alphab,kk) + wall

         alocont_nn2d_tmp(iim1,kkm1,q24) = alocont_nn2d_tmp(iim1,kkm1,q24) + wall*alocont_o_nn2d(iim1,kkm1,q24)
         alocont_nn2d_tmp(ii,kkm1,q23) = alocont_nn2d_tmp(ii,kkm1,q23) + wall*alocont_o_nn2d(ii,kkm1,q23)
         alocont_nn2d_tmp(iip1,kkm1,q22) = alocont_nn2d_tmp(iip1,kkm1,q23) + wall*alocont_o_nn2d(iip1,kkm1,q23)
         alocont_nn2d_tmp(iim1,kk,q15) = alocont_nn2d_tmp(iim1,kk,q15) + wall*alocont_o_nn2d(iim1,kk,q15)
         alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alocont_o_nn2d(ii,kk,q14)
         alocont_nn2d_tmp(iip1,kk,q13) = alocont_nn2d_tmp(iip1,kk,q13) + wall*alocont_o_nn2d(iip1,kk,q13)
         alocont_nn2d_tmp(iim1,kkp1,q6) = alocont_nn2d_tmp(iim1,kkp1,q6) + wall*alocont_o_nn2d(iim1,kkp1,q6)
         alocont_nn2d_tmp(ii,kkp1,q5) = alocont_nn2d_tmp(ii,kkp1,q5) + wall*alocont_o_nn2d(ii,kkp1,q5)
         alocont_nn2d_tmp(iip1,kkp1,q4) = alocont_nn2d_tmp(iip1,kkp1,q4) + wall*alocont_o_nn2d(iip1,kkp1,q4)

!again, left side for q6, q15 and q24 not possible
         alocont_nn2d_tmp(endxb+alphab,kkm1,q23) = alocont_nn2d_tmp(endxb+alphab,kkm1,q23) + wall*alocont_o_nn2d(endxb+alphab,kkm1,q23)
         alocont_nn2d_tmp(endxb,kkm1,q22) = alocont_nn2d_tmp(endxb,kkm1,q22) + wall*alocont_o_nn2d(endxb,kkm1,q22)
         alocont_nn2d_tmp(endxb+alphab,kk,q14) = alocont_nn2d_tmp(endxb+alphab,kk,q14) + wall*alocont_o_nn2d(endxb+alphab,kk,q14)
         alocont_nn2d_tmp(endxb,kk,q13) = alocont_nn2d_tmp(endxb,kk,q13) + wall*alocont_o_nn2d(endxb,kk,q13)
         alocont_nn2d_tmp(endxb+alphab,kkp1,q5) = alocont_nn2d_tmp(endxb+alphab,kkp1,q5) + wall*alocont_o_nn2d(endxb+alphab,kkp1,q5)
         alocont_nn2d_tmp(endxb,kkp1,q4) = alocont_nn2d_tmp(endxb,kkp1,q4) + wall*alocont_o_nn2d(endxb,kkp1,q4)
!
!exit the loop since everything is calculated

         exit
!
      else
         write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
         stop 'error in set_boundaryp2d: invalid dels_u'
      endif
!
!update the downwind and current point for next step
         x_d=x_p
         z_d=z_p
         x_p=x_u
         z_p=z_u
!      
   enddo
!
enddo
!
!------------second step: calculate intensities at ghost point----------
!-----------one grid point further away from computational region-------
!---------------(point G1l for n_x>0 and point G2r for n_x < 0)---------
!--------------------------using SC method------------------------------
!
!G2l----------G1l--------x--------------x----------x-------G2r--------G1r
!
!set indices
ii=endxb
iim1=ii+alphab
iip1=ii-alphab
!
kkm1=kk+gammab
kkm2=kkm1+gammab
kkp1=kk-gammab
!
!calculate distance to intersection point with upwind x-grid-coordinate
dels_xu=(x(ii)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
dels_zu=(z(kk)-z(kkm1))/nn_z
dels_u=min(dels_xu,dels_zu)

dels_xd=(x(iip1)-x(ii))/nn_x
dels_zd=(z(kkp1)-z(kk))/nn_z
dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
opac_p=opac2d(ii,kk)
scont_p=scont2d(ii,kk)
!
!----------------------------upwind point-------------------------------
!
if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
   x_u = x(ii) - dels_u*nn_x
   z_u = z(kk) - dels_u*nn_z
 
   call coeff1d_contu(opac2d(iim1,kkm2), opac2d(iim1,kkm1), opac2d(iim1,kk), &
                      scont2d(iim1,kkm2), scont2d(iim1,kkm1), scont2d(iim1,kk), &
                      int2d(iim1,kkm2), int2d(iim1,kkm1), int2d(iim1,kk), &
                      z(kkm2), z(kkm1), z(kk), z_u, &
                      c01_scontu, c03_scontu, c05_scontu, &
                      c01_intu, c03_intu, c05_intu, &
                      opac_u, scont_u, int_u)
   c02_scontu=zero
   c04_scontu=zero
   c02_intu=zero
   c04_intu=zero
!
!debug start
!   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(iim1,kk), &
!                          scont2d(iim1,kkm1), scont2d(iim1,kk), &
!                          int2d(iim1,kkm1), int2d(iim1,kk), &
!                          z(kkm1), z(kk), z_u, &
!                          c03_scontu, c05_scontu, &
!                          c03_intu, c05_intu, &
!                          opac_u, scont_u, int_u)
!   c01_scontu=zero
!   c02_scontu=zero
!   c04_scontu=zero
!   c01_intu=zero
!   c02_intu=zero
!   c04_intu=zero
!debug end                  
!
elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
   x_u=x(ii)-dels_u*nn_x
   z_u=z(kk)-dels_u*nn_z
!
!only linear interpolations, because iim2 is not know (starts at right point again)   
   call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(ii,kkm1), &
                          scont2d(iim1,kkm1), scont2d(ii,kkm1), &
                          int2d(iim1,kkm1), int2d(ii,kkm1), &
                          x(iim1), x(ii), x_u, &
                          c03_scontu, c04_scontu, &
                          c03_intu, c04_intu, &
                          opac_u, scont_u, int_u)
   c05_scontu=zero
   c05_intu=zero
else
   write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
   stop 'error in set_boundaryp2d: invalid dels_u'
endif
!
!---------------------------downwind point------------------------------
!
if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd(opac2d(iip1,kkm1), opac2d(iip1,kk), opac2d(iip1,kkp1), &
                      scont2d(iip1,kkm1), scont2d(iip1,kk), scont2d(iip1,kkp1), &
                      z(kkm1), z(kk), z(kkp1), z_d, &
                      c01_scontd, c02_scontd, c05_scontd, &
                      opac_d, scont_d)
   c03_scontd=zero
   c04_scontd=zero
!
!debug start
!   call coeff1d_contd_lin(opac2d(iip1,kk), opac2d(iip1,kkp1), &
!                          scont2d(iip1,kk), scont2d(iip1,kkp1), &
!                          z(kk), z(kkp1), z_d, &
!                          c02_scontd, c05_scontd, &
!                          opac_d, scont_d)
!   c01_scontd=zero
!   c03_scontd=zero
!   c04_scontd=zero
!debug end
!
elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
   x_d=x(ii)+dels_d*nn_x
   z_d=z(kk)+dels_d*nn_z
!
   call coeff1d_contd(opac2d(iim1,kkp1), opac2d(ii,kkp1), opac2d(iip1,kkp1), &
                      scont2d(iim1,kkp1), scont2d(i,kkp1), scont2d(iip1,kkp1), &
                      x(iim1), x(ii), x(iip1), x_d, &
                      c03_scontd, c04_scontd, c05_scontd, &
                      opac_d, scont_d)
   c01_scontd=zero
   c02_scontd=zero
!
!debug start   
!   call coeff1d_contd_lin(opac2d(ii,kkp1), opac2d(iip1,kkp1), &
!                          scont2d(ii,kkp1), scont2d(iip1,kkp1), &
!                          x(ii), x(iip1), x_d, &
!                          c04_scontd, c05_scontd, &
!                          opac_d, scont_d)
!   c01_scontd=zero
!   c02_scontd=zero
!   c03_scontd=zero
!debug end
else
   write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
   stop 'error in set_boundaryp2d: invalid dels_d'
endif
!
!-----------------------radiative transfer------------------------------
!
!call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!              dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                  dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
alo_d=zero
!write(*,*) scont_p, scont_u, scont_d
!write(*,'(i5,13es16.8)') ii, x(ii), x_p, opac_u, opac_p, opac_d, dels_u, dels_d, alo_u, alo_p, dels_d, abs_sc
!write(*,*) startxb-alphab, ii
int2d(ii,kk) = int_sc
int2d(startxb-alphab,kk) = int_sc
intbound2d(bindxb,kk,oindx)=int2d(ii,kk)

!if(kk.eq.4) then
!   write(*,*) ii, kk, nn_z, int2d(ii,kk), int_u
!   write(*,*) startxb-alphab, kk, nn_z, int2d(startxb-alphab,kk), int_u
!endif


alocont_o_nn2d(iip1,kkp1,q4) = alo_d*c05_scontd
alocont_o_nn2d(ii,kkp1,q5) = (alo_d*c04_scontd + abs_sc*(c05_intu*alocont_o_nn2d(ii,kkp1,q4)))
alocont_o_nn2d(iim1,kkp1,q6) = (alo_d*c03_scontd + abs_sc*(c05_intu*alocont_o_nn2d(iim1,kkp1,q5)))
alocont_o_nn2d(iip1,kk,q13) = (alo_d*c02_scontd + abs_sc*(c04_intu*alocont_o_nn2d(iip1,kk,q4)))
alocont_o_nn2d(ii,kk,q14) = (alo_p + abs_sc*((c03_intu*alocont_o_nn2d(ii,kk,q4) + &
                                           c04_intu*alocont_o_nn2d(ii,kk,q5) + &
                                           c05_intu*alocont_o_nn2d(ii,kk,q13))))
alocont_o_nn2d(iim1,kk,q15) = (alo_u*c05_scontu + abs_sc*((c02_intu*alocont_o_nn2d(iim1,kk,q4) + &
                                                                     c03_intu*alocont_o_nn2d(iim1,kk,q5) + &
                                                                     c04_intu*alocont_o_nn2d(iim1,kk,q6) + &
                                                                     c05_intu*alocont_o_nn2d(iim1,kk,q14))))
alocont_o_nn2d(iip1,kkm1,q22) = (alo_d*c01_scontd + abs_sc*((c04_intu*alocont_o_nn2d(iip1,kkm1,q13))))
alocont_o_nn2d(ii,kkm1,q23) = (alo_u*c04_scontu + abs_sc*((c03_intu*alocont_o_nn2d(ii,kkm1,q13) + &
                                                         c04_intu*alocont_o_nn2d(ii,kkm1,q14) + &
                                                         c05_intu*alocont_o_nn2d(ii,kkm1,q22) + &
                                                         c01_intu*alocont_o_nn2d(ii,kkm1,q4))))
!write(*,*) scont2d(ii,kkm1), alocont_o_nn2d(ii,kkm1,q23), int2d(ii,kk), int2d(iim1,kk)
!write(*,*)
alocont_o_nn2d(iim1,kkm1,q24) = (alo_u*c03_scontu + abs_sc*((c02_intu*alocont_o_nn2d(iim1,kkm1,q13) + &
                                                            c03_intu*alocont_o_nn2d(iim1,kkm1,q14) + &
                                                            c04_intu*alocont_o_nn2d(iim1,kkm1,q15) + &
                                                            c05_intu*alocont_o_nn2d(iim1,kkm1,q23) + &
                                                            c01_intu*alocont_o_nn2d(iim1,kkm1,q5))))
!
!copy everyhing on right side (q4, q13, q22 cannot be stored because out of bounds)
alocont_o_nn2d(startxb-alphab,kkp1,q5) = alocont_o_nn2d(ii,kkp1,q5)
alocont_o_nn2d(startxb,kkp1,q6) = alocont_o_nn2d(iim1,kkp1,q6)
alocont_o_nn2d(startxb-alphab,kk,q14) = alocont_o_nn2d(ii,kk,q14)
alocont_o_nn2d(startxb,kk,q15) = alocont_o_nn2d(iim1,kk,q15)
alocont_o_nn2d(startxb-alphab,kkm1,q23) = alocont_o_nn2d(ii,kkm1,q23)
alocont_o_nn2d(startxb,kkm1,q24) = alocont_o_nn2d(iim1,kkm1,q24)
!
!perform angular integration
mint2d_tmp(ii,kk) = mint2d_tmp(ii,kk) + int2d(ii,kk)*wall
mint2d_tmp(startxb-alphab,kk) = mint2d_tmp(startxb-alphab,kk) + int2d(startxb-alphab,kk)*wall
!
normalization2d_tmp(ii,kk) = normalization2d_tmp(ii,kk) + wall
normalization2d_tmp(startxb-alphab,kk) = normalization2d_tmp(startxb-alphab,kk) + wall
!
alocont_nn2d_tmp(iip1,kkp1,q4) = alocont_nn2d_tmp(iip1,kkp1,q4) + wall*alocont_o_nn2d(iip1,kkp1,q4)
alocont_nn2d_tmp(ii,kkp1,q5) = alocont_nn2d_tmp(ii,kkp1,q5) + wall*alocont_o_nn2d(ii,kkp1,q5)
alocont_nn2d_tmp(iim1,kkp1,q6) = alocont_nn2d_tmp(iim1,kkp1,q6) + wall*alocont_o_nn2d(iim1,kkp1,q6)
alocont_nn2d_tmp(iip1,kk,q13) = alocont_nn2d_tmp(iip1,kk,q13) + wall*alocont_o_nn2d(iip1,kk,q13)
alocont_nn2d_tmp(ii,kk,q14) = alocont_nn2d_tmp(ii,kk,q14) + wall*alocont_o_nn2d(ii,kk,q14)
alocont_nn2d_tmp(iim1,kk,q15) = alocont_nn2d_tmp(iim1,kk,q15) + wall*alocont_o_nn2d(iim1,kk,q15)
alocont_nn2d_tmp(iip1,kkm1,q22) = alocont_nn2d_tmp(iip1,kkm1,q22) + wall*alocont_o_nn2d(iip1,kkm1,q22)
alocont_nn2d_tmp(ii,kkm1,q23) = alocont_nn2d_tmp(ii,kkm1,q23) + wall*alocont_o_nn2d(ii,kkm1,q23)
alocont_nn2d_tmp(iim1,kkm1,q24) = alocont_nn2d_tmp(iim1,kkm1,q24) + wall*alocont_o_nn2d(iim1,kkm1,q24)
!
!copy everything on right side
alocont_nn2d_tmp(startxb-alphab,kkp1,q5) = alocont_nn2d_tmp(startxb-alphab,kkp1,q5) + wall*alocont_o_nn2d(startxb-alphab,kkp1,q5)
alocont_nn2d_tmp(startxb,kkp1,q6) = alocont_nn2d_tmp(startxb,kkp1,q6) + wall*alocont_o_nn2d(startxb,kkp1,q6)
alocont_nn2d_tmp(startxb-alphab,kk,q14) = alocont_nn2d_tmp(startxb-alphab,kk,q14) + wall*alocont_o_nn2d(startxb-alphab,kk,q14)
alocont_nn2d_tmp(startxb,kk,q15) = alocont_nn2d_tmp(startxb,kk,q15) + wall*alocont_o_nn2d(startxb,kk,q15)
alocont_nn2d_tmp(startxb-alphab,kkm1,q23) = alocont_nn2d_tmp(startxb-alphab,kkm1,q23) + wall*alocont_o_nn2d(startxb-alphab,kkm1,q23)
alocont_nn2d_tmp(startxb,kkm1,q24) = alocont_nn2d_tmp(startxb,kkm1,q24) + wall*alocont_o_nn2d(startxb,kkm1,q24)

!
!
end subroutine set_boundaryp2d_v02
!
!***********************************************************************
!***********************************************************************
!
!                 CONTINUUM TRASNFER ROUTINES
!
!***********************************************************************
!***********************************************************************
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_cont2d(oindx,nueindx)
!
!-----------------------------------------------------------------------
!------short characteristics for continuum radiative transfer in 2d-----
!---------------calculating intensties for given mu---------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime2d, only: nx, nz, x, z, int2d, opac2d, scont2d, alocont_o_nn2d, imask2d, imaskb2d, &
                  alocont_nn2d, alocont_o_nn2d, mint2d_tmp, alocont_nn2d_tmp, normalization2d_tmp
use angles, only: n_x, n_z, weight_omega, q_alo
use freq, only: nodes_nue
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: oindx,nueindx
!
! ... local scalars
integer(i4b) :: i, k
integer(i4b) :: alpha, gamma
integer(i4b) :: startx, startz, endx, endz
integer(i4b) :: startxb, endxb, alphab, gammab
integer(i4b) :: bindxa, bindxb
integer(i4b) :: iim2, iim1, ii, iip1, kkm2, kkm1, kk, kkp1
integer :: q4, q5, q6, q13, q14, q15, q22, q23, q24
real(dp) ::  xnue, nn_x, nn_z, wall
real(dp) :: opac_u, scont_u, dels_u, dels_xu, dels_zu, int_u, &
            opac_d, scont_d, dels_d, dels_xd, dels_zd, &
            opac_p, scont_p
real(dp) :: x_u, z_u, x_d, z_d
real(dp) :: abs_sc, int_sc, contr_sc, alo_p, alo_u, alo_d
real(dp) :: c01_scontu, c02_scontu, c03_scontu, c04_scontu, c05_scontu
real(dp) :: c01_intu, c02_intu, c03_intu, c04_intu, c05_intu
real(dp) :: c01_scontd, c02_scontd, c03_scontd, c04_scontd, c05_scontd
!
! ... for debugging
real(dp) :: fdum
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
! ... local logicals
!
!frequency
xnue=nodes_nue(nueindx)

!note: nn_x always positive because of symmetry
nn_x=n_x(oindx)
nn_z=n_z(oindx)
!
wall=weight_omega(oindx)
!
!indices for nearest neighbour alo
q4=q_alo(oindx,4)
q5=q_alo(oindx,5)
q6=q_alo(oindx,6)
q13=q_alo(oindx,13)
q14=q_alo(oindx,14)
q15=q_alo(oindx,15)
q22=q_alo(oindx,22)
q23=q_alo(oindx,23)
q24=q_alo(oindx,24)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
!index parameter:
!         if n_x >= 0                 if n_x < 0                
!                startx = 2                  startx = nx-1  
!                endx = nx-1                 endx = 2           
!                alpha=  1                   alpha=-1                  
!
!         if n_z >= 0                 if n_z < 0                
!                startz = 2                  startz = nz-1  
!                endz = nz-1                 endz = 2           
!                gamma = 1                   gamma = -1
!
!could make an own array for this!!!
if(nn_x.gt.zero) then
   startx = 3
   endx = nx-2
   alpha=  1
   startxb = nx-1
   endxb = 2
   alphab = -1
   bindxa = 1
   bindxb = 2
elseif(nn_x.lt.zero) then
   startx = nx-2
   endx = 3
   alpha=-1
   startxb = 2
   endxb = nx-1
   alphab = 1
   bindxa = 2
   bindxb = 1
else
   stop 'error in fsc_cont2d: n_x = 0 not allowed'
endif
!
if(nn_z.gt.zero) then
   startz = 3
   endz = nz-2
   gamma=  1
   gammab = -1
elseif(nn_z.lt.zero) then
   startz = nz-2
   endz = 3
   gamma=-1
   gammab = 1
else
   stop 'error in fsc_cont2d: n_z = 0 not allowed'
endif
!
!-----------------------reset the intensities---------------------------
!
call set_boundary2d(xnue,nn_x,nn_z)
!
!-----------------------------------------------------------------------
!
alocont_o_nn2d=zero
!
!-----------------------------------------------------------------------
!
do k=startz, endz, gamma
!
!   call set_boundaryp2d_v02(nn_x, nn_z, startxb,endxb,alphab,gammab,k,bindxa,bindxb,oindx,wall,q4,q5,q6,q13,q14,q15,q22,q23,q24)
   call set_boundaryp2d(nn_x, nn_z, startxb,endxb,alphab,gammab,k,bindxa,bindxb,oindx,wall,q4,q5,q6,q13,q14,q15,q22,q23,q24)
!   
   do i=startx, endx, alpha

      select case(imaskb2d(i,k))
!
!*************************standard rt procedure*************************
!
         case(9)

            iim1=i-alpha
            kkm1=k-gamma
            iip1=i+alpha
            kkp1=k+gamma
!
!calculate distance to intersection point with upwind x-grid-coordinate
            dels_xu=(x(i)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
            dels_zu=(z(k)-z(kkm1))/nn_z
            dels_u=min(dels_xu,dels_zu)
!
!calculate distance to intersection point with downwind x-grid-coordinate
            dels_xd=(x(iip1)-x(i))/nn_x
!calculate distance to instersection point with downwind z-grid-coordinate
            dels_zd=(z(kkp1)-z(k))/nn_z
            dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
            opac_p=opac2d(i,k)
            scont_p=scont2d(i,k)
!
!----------------------------upwind point-------------------------------
!
            if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
               x_u = x(i) - dels_u*nn_x
               z_u = z(k) - dels_u*nn_z
!
               kkm2=k-2*gamma
               call coeff1d_contu(opac2d(iim1,kkm2), opac2d(iim1,kkm1), opac2d(iim1,k), &
                                  scont2d(iim1,kkm2), scont2d(iim1,kkm1), scont2d(iim1,k), &
                                  int2d(iim1,kkm2), int2d(iim1,kkm1), int2d(iim1,k), &
                                  z(kkm2), z(kkm1), z(k), z_u, &
                                  c01_scontu, c03_scontu, c05_scontu, &
                                  c01_intu, c03_intu, c05_intu, &
                                  opac_u, scont_u, int_u)
               c02_scontu=zero
               c04_scontu=zero
               c02_intu=zero
               c04_intu=zero
!
!debug start
!               call coeff1d_contu_linb(opac2d(iim1,kkm1), opac2d(iim1,k), &
!                                      scont2d(iim1,kkm1), scont2d(iim1,k), &
!                                      int2d(iim1,kkm1), int2d(iim1,k), &
!                                      z(kkm1), z(k), z_u, &
!                                      c03_scontu, c05_scontu, &
!                                      c03_intu, c05_intu, &
!                                      opac_u, scont_u, int_u)
!               c01_scontu=zero
!               c02_scontu=zero
!               c04_scontu=zero
!               c01_intu=zero
!               c02_intu=zero
!               c04_intu=zero
!debug end               
            elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
               x_u=x(i)-dels_u*nn_x
               z_u=z(k)-dels_u*nn_z
!
               iim2=i-2*alpha

!               call coeff1d_contu(opac2d(iim2,kkm1), opac2d(iim1,kkm1), opac2d(i,kkm1), &
!                                  scont2d(iim2,kkm1), scont2d(iim1,kkm1), scont2d(i,kkm1), &
!                                  int2d(iim2,kkm1), int2d(iim1,kkm1), int2d(i,kkm1), &
!                                  x(iim2), x(iim1), x(i), x_u, &
!                                  c02_scontu, c03_scontu, c04_scontu, &
!                                  c02_intu, c03_intu, c04_intu, &
!                                  opac_u, scont_u, int_u)
!               c01_scontu=zero
!               c05_scontu=zero
!               c01_intu=zero
!               c05_intu=zero
!debug start
               call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(i,kkm1), &
                                      scont2d(iim1,kkm1), scont2d(i,kkm1), &
                                      int2d(iim1,kkm1), int2d(i,kkm1), &
                                      x(iim1), x(i), x_u, &
                                      c03_scontu, c04_scontu, &
                                      c03_intu, c04_intu, &
                                      opac_u, scont_u, int_u)
               c01_scontu=zero
               c02_scontu=zero
               c05_scontu=zero
               c01_intu=zero
               c02_intu=zero
               c05_intu=zero
!debug end
            else
               write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
               stop 'error in fsc_cont2d: invalid dels_u'
            endif
!
!--------------------------downwind point-------------------------------
!
            if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
               x_d=x(i)+dels_d*nn_x
               z_d=z(k)+dels_d*nn_z
!
               call coeff1d_contd(opac2d(iip1,kkm1), opac2d(iip1,k), opac2d(iip1,kkp1), &
                                  scont2d(iip1,kkm1), scont2d(iip1,k), scont2d(iip1,kkp1), &
                                  z(kkm1), z(k), z(kkp1), z_d, &
                                  c01_scontd, c02_scontd, c05_scontd, &
                                  opac_d, scont_d)
               c03_scontd=zero
               c04_scontd=zero
!debug start
!               call coeff1d_contd_lin(opac2d(iip1,k), opac2d(iip1,kkp1), &
!                                      scont2d(iip1,k), scont2d(iip1,kkp1), &
!                                      z(k), z(kkp1), z_d, &
!                                      c02_scontd, c05_scontd, &
!                                      opac_d, scont_d)
!               c01_scontd=zero
!               c03_scontd=zero
!               c04_scontd=zero
!debug end
            elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
               x_d=x(i)+dels_d*nn_x
               z_d=z(k)+dels_d*nn_z
!
!               call coeff1d_contd(opac2d(iim1,kkp1), opac2d(i,kkp1), opac2d(iip1,kkp1), &
!                                  scont2d(iim1,kkp1), scont2d(i,kkp1), scont2d(iip1,kkp1), &
!                                  x(iim1), x(i), x(iip1), x_d, &
!                                  c03_scontd, c04_scontd, c05_scontd, &
!                                  opac_d, scont_d)
!               c01_scontd=zero
!               c02_scontd=zero
!debug start
               call coeff1d_contd_lin(opac2d(i,kkp1), opac2d(iip1,kkp1), &
                                      scont2d(i,kkp1), scont2d(iip1,kkp1), &
                                      x(i), x(iip1), x_d, &
                                      c04_scontd, c05_scontd, &
                                      opac_d, scont_d)
               c01_scontd=zero
               c02_scontd=zero
               c03_scontd=zero
!debug end

            else
               write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
               stop 'error in fsc_cont2d: invalid dels_d'
            endif
!
!-----------------------radiative transfer------------------------------
!
!            write(*,*) i, k, dels_u, x(i), x(iim1), x(i)-x(iim1)
            call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                 dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!            write(*,*) alo_u, alo_p, alo_d, alo_u+alo_p+alo_d
!            call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                              dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!            alo_d=zero
!            int2d(i,k) = int_sc
!            write(*,*) i, k, x(i), z(k), alo_u+alo_p+alo_d+abs_sc
!             write(*,*)
            
            alocont_o_nn2d(iip1,kkp1,q4) = alo_d*c05_scontd
            alocont_o_nn2d(i,kkp1,q5) = (alo_d*c04_scontd + abs_sc*(c05_intu*alocont_o_nn2d(i,kkp1,q4)))
            alocont_o_nn2d(iim1,kkp1,q6) = (alo_d*c03_scontd + abs_sc*(c05_intu*alocont_o_nn2d(iim1,kkp1,q5)))
            alocont_o_nn2d(iip1,k,q13) = (alo_d*c02_scontd + abs_sc*(c04_intu*alocont_o_nn2d(iip1,k,q4)))
            alocont_o_nn2d(i,k,q14) = (alo_p + abs_sc*((c03_intu*alocont_o_nn2d(i,k,q4) + &
                                                       c04_intu*alocont_o_nn2d(i,k,q5) + &
                                                       c05_intu*alocont_o_nn2d(i,k,q13))))
            alocont_o_nn2d(iim1,k,q15) = (alo_u*c05_scontu + abs_sc*((c02_intu*alocont_o_nn2d(iim1,k,q4) + &
                                                                     c03_intu*alocont_o_nn2d(iim1,k,q5) + &
                                                                     c04_intu*alocont_o_nn2d(iim1,k,q6) + &
                                                                     c05_intu*alocont_o_nn2d(iim1,k,q14))))
            alocont_o_nn2d(iip1,kkm1,q22) = (alo_d*c01_scontd + abs_sc*((c04_intu*alocont_o_nn2d(iip1,kkm1,q13))))
            alocont_o_nn2d(i,kkm1,q23) = (alo_u*c04_scontu + abs_sc*((c03_intu*alocont_o_nn2d(i,kkm1,q13) + &
                                                                     c04_intu*alocont_o_nn2d(i,kkm1,q14) + &
                                                                     c05_intu*alocont_o_nn2d(i,kkm1,q22) + &
                                                                     c01_intu*alocont_o_nn2d(i,kkm1,q4))))
            alocont_o_nn2d(iim1,kkm1,q24) = (alo_u*c03_scontu + abs_sc*((c02_intu*alocont_o_nn2d(iim1,kkm1,q13) + &
                                                                        c03_intu*alocont_o_nn2d(iim1,kkm1,q14) + &
                                                                        c04_intu*alocont_o_nn2d(iim1,kkm1,q15) + &
                                                                        c05_intu*alocont_o_nn2d(iim1,kkm1,q23) + &
                                                                        c01_intu*alocont_o_nn2d(iim1,kkm1,q5))))


!            if(k.eq.18.and.i.eq.nx-2) then
!               write(*,*) i, k, int2d(i,k), int2d(iim1,k), int2d(iim1,kkm1), abs_sc*int_u, abs_sc*(c05_intu*alocont_o_nn2d(iim1,kkp1,q5))
!               write(*,*) i, k, int2d(i,k), scont2d(i+1,k), alocont_o_nn2d(iip1,k,q13), int_u, int2d(i,k)               
!            endif
            
!            write(*,*) i, k, scont2d(i,k), int2d(i,k), int_u, contr_sc, abs_sc!int2d(2,k), int2d(nx,k)
!            if(i.eq.5.and.k.eq.21) then
!               write(*,*) int2d(i,k), alocont_o_nn2d(iip1,kkm1,q22), int_u, c04_intu, &
!                    int2d(i,k), int2d(iim1,k), int2d(i,kkm1), int2d(2,k)
!               stop
!            endif
!            if(i.eq.16.and.k.eq.25) then
!               write(*,*) 'fsc', i, k, alocont_o_nn2d(i,k,q14), int2d(i,k)
!            endif
!            if(i.eq.17.and.k.eq.25) then
!               write(*,*) 'fsc', i, k, alocont_o_nn2d(iim1,k,q15), int2d(i,k)
!            endif
!            if(i.eq.15.and.k.eq.26) then
!               write(*,*) 'fsc', i, k, alocont_o_nn2d(iip1,kkm1,q22), int2d(i,k)
!               write(*,*) abs_sc, int_u, c05_intu, alocont_o_nn2d(iip1,kkm1,q13), c05_intu*alocont_o_nn2d(iip1,kkm1,q13), c05_intu*alocont_o_nn2d(iip1,kkm1,q13)*abs_sc, int2d(i,kkm!1)
!            endif
!            if(k.eq.4) then
!               write(*,*) i, k, nn_z, int2d(i,k), int_u, opac_d, opac2d(iip1,kkp1), (x_d-x(i)),(x(iip1)-x(i)), c04_scontd, c05_scontd
!            endif
!            if(k.eq.3.and.i.eq.16) then
!               write(*,*) i, k, alocont_o_nn2d(i,k,q14), alo_u, alo_p, alo_d, alo_u+alo_p+alo_d+abs_sc, alocont_o_nn2d(i,k,q4), alocont_o_nn2d(i,k,q13)
!               write(*,*) i, k, nn_z, opac_u, opac_p, opac_d, x_u, z_u, x_d, z_d
!               write(*,*)
!            endif
!
!            write(*,*) alocont_o_nn2d(iip1,kkp1,q4), alocont_o_nn2d(i,kkp1,q5), alocont_o_nn2d(iim1,kkp1,q6), alocont_o_nn2d(i,k,q14), alocont_o_nn2d(iip1,kkm1,q22), alocont_o_nn2d(i,kkm1,q23), alocont_o_nn2d(iim1,kkm1,q24)

            mint2d_tmp(i,k) = mint2d_tmp(i,k) + int2d(i,k)*wall
            normalization2d_tmp(i,k) = normalization2d_tmp(i,k) + wall
!        
            alocont_nn2d_tmp(iip1,kkp1,q4) = alocont_nn2d_tmp(iip1,kkp1,q4) + wall*alocont_o_nn2d(iip1,kkp1,q4)
            alocont_nn2d_tmp(i,kkp1,q5) = alocont_nn2d_tmp(i,kkp1,q5) + wall*alocont_o_nn2d(i,kkp1,q5)
            alocont_nn2d_tmp(iim1,kkp1,q6) = alocont_nn2d_tmp(iim1,kkp1,q6) + wall*alocont_o_nn2d(iim1,kkp1,q6)
            alocont_nn2d_tmp(iip1,k,q13) = alocont_nn2d_tmp(iip1,k,q13) + wall*alocont_o_nn2d(iip1,k,q13)
            alocont_nn2d_tmp(i,k,q14) = alocont_nn2d_tmp(i,k,q14) + wall*alocont_o_nn2d(i,k,q14)
            alocont_nn2d_tmp(iim1,k,q15) = alocont_nn2d_tmp(iim1,k,q15) + wall*alocont_o_nn2d(iim1,k,q15)
            alocont_nn2d_tmp(iip1,kkm1,q22) = alocont_nn2d_tmp(iip1,kkm1,q22) + wall*alocont_o_nn2d(iip1,kkm1,q22)
            alocont_nn2d_tmp(i,kkm1,q23) = alocont_nn2d_tmp(i,kkm1,q23) + wall*alocont_o_nn2d(i,kkm1,q23)
            alocont_nn2d_tmp(iim1,kkm1,q24) = alocont_nn2d_tmp(iim1,kkm1,q24) + wall*alocont_o_nn2d(iim1,kkm1,q24)
!
!*****************linear procedure near periodice boundary****************
!
         case(5,6)

            iim1=i-alpha
            kkm1=k-gamma
            iip1=i+alpha
            kkp1=k+gamma
!
!calculate distance to intersection point with upwind x-grid-coordinate
            dels_xu=(x(i)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
            dels_zu=(z(k)-z(kkm1))/nn_z
            dels_u=min(dels_xu,dels_zu)
!
!calculate distance to intersection point with downwind x-grid-coordinate
            dels_xd=(x(iip1)-x(i))/nn_x
!calculate distance to instersection point with downwind z-grid-coordinate
            dels_zd=(z(kkp1)-z(k))/nn_z
            dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
            opac_p=opac2d(i,k)
            scont_p=scont2d(i,k)
!
!----------------------------upwind point-------------------------------
!
            if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
               x_u = x(i) - dels_u*nn_x
               z_u = z(k) - dels_u*nn_z
!
               kkm2=k-2*gamma
               call coeff1d_contu(opac2d(iim1,kkm2), opac2d(iim1,kkm1), opac2d(iim1,k), &
                                  scont2d(iim1,kkm2), scont2d(iim1,kkm1), scont2d(iim1,k), &
                                  int2d(iim1,kkm2), int2d(iim1,kkm1), int2d(iim1,k), &
                                  z(kkm2), z(kkm1), z(k), z_u, &
                                  c01_scontu, c03_scontu, c05_scontu, &
                                  c01_intu, c03_intu, c05_intu, &
                                  opac_u, scont_u, int_u)
               c02_scontu=zero
               c04_scontu=zero
               c02_intu=zero
               c04_intu=zero
!
!debug start
!               call coeff1d_contu_linb(opac2d(iim1,kkm1), opac2d(iim1,k), &
!                                      scont2d(iim1,kkm1), scont2d(iim1,k), &
!                                      int2d(iim1,kkm1), int2d(iim1,k), &
!                                      z(kkm1), z(k), z_u, &
!                                      c03_scontu, c05_scontu, &
!                                      c03_intu, c05_intu, &
!                                      opac_u, scont_u, int_u)
!               c01_scontu=zero
!               c02_scontu=zero
!               c04_scontu=zero
!               c01_intu=zero
!               c02_intu=zero
!               c04_intu=zero
!debug end               
            elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
               x_u=x(i)-dels_u*nn_x
               z_u=z(k)-dels_u*nn_z
!
               iim2=i-2*alpha

               call coeff1d_contu(opac2d(iim2,kkm1), opac2d(iim1,kkm1), opac2d(i,kkm1), &
                                  scont2d(iim2,kkm1), scont2d(iim1,kkm1), scont2d(i,kkm1), &
                                  int2d(iim2,kkm1), int2d(iim1,kkm1), int2d(i,kkm1), &
                                  x(iim2), x(iim1), x(i), x_u, &
                                  c02_scontu, c03_scontu, c04_scontu, &
                                  c02_intu, c03_intu, c04_intu, &
                                  opac_u, scont_u, int_u)
               c01_scontu=zero
               c05_scontu=zero
               c01_intu=zero
               c05_intu=zero
!debug start
!               call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(i,kkm1), &
!                                      scont2d(iim1,kkm1), scont2d(i,kkm1), &
!                                      int2d(iim1,kkm1), int2d(i,kkm1), &
!                                      x(iim1), x(i), x_u, &
!                                      c03_scontu, c04_scontu, &
!                                      c03_intu, c04_intu, &
!                                      opac_u, scont_u, int_u)
!               c01_scontu=zero
!               c02_scontu=zero
!               c05_scontu=zero
!               c01_intu=zero
!               c02_intu=zero
!               c05_intu=zero
!debug end
            else
               write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
               stop 'error in fsc_cont2d: invalid dels_u'
            endif
!
!--------------------------downwind point-------------------------------
!
            if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
               x_d=x(i)+dels_d*nn_x
               z_d=z(k)+dels_d*nn_z
!
               call coeff1d_contd(opac2d(iip1,kkm1), opac2d(iip1,k), opac2d(iip1,kkp1), &
                                  scont2d(iip1,kkm1), scont2d(iip1,k), scont2d(iip1,kkp1), &
                                  z(kkm1), z(k), z(kkp1), z_d, &
                                  c01_scontd, c02_scontd, c05_scontd, &
                                  opac_d, scont_d)
               c03_scontd=zero
               c04_scontd=zero
!debug start
!               call coeff1d_contd_lin(opac2d(iip1,k), opac2d(iip1,kkp1), &
!                                      scont2d(iip1,k), scont2d(iip1,kkp1), &
!                                      z(k), z(kkp1), z_d, &
!                                      c02_scontd, c05_scontd, &
!                                      opac_d, scont_d)
!               c01_scontd=zero
!               c03_scontd=zero
!               c04_scontd=zero
!debug end
            elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
               x_d=x(i)+dels_d*nn_x
               z_d=z(k)+dels_d*nn_z
!
               call coeff1d_contd(opac2d(iim1,kkp1), opac2d(i,kkp1), opac2d(iip1,kkp1), &
                                  scont2d(iim1,kkp1), scont2d(i,kkp1), scont2d(iip1,kkp1), &
                                  x(iim1), x(i), x(iip1), x_d, &
                                  c03_scontd, c04_scontd, c05_scontd, &
                                  opac_d, scont_d)
               c01_scontd=zero
               c02_scontd=zero
!debug start
!               call coeff1d_contd_lin(opac2d(i,kkp1), opac2d(iip1,kkp1), &
!                                      scont2d(i,kkp1), scont2d(iip1,kkp1), &
!                                      x(i), x(iip1), x_d, &
!                                      c04_scontd, c05_scontd, &
!                                      opac_d, scont_d)
!               c01_scontd=zero
!               c02_scontd=zero
!               c03_scontd=zero
!debug end

            else
               write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
               stop 'error in fsc_cont2d: invalid dels_d'
            endif
!
!-----------------------radiative transfer------------------------------
!
!            write(*,*) i, k, dels_u, x(i), x(iim1), x(i)-x(iim1)
!            call fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
!                 dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!            write(*,*) alo_u, alo_p, alo_d, alo_u+alo_p+alo_d
            call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                              dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
            alo_d=zero
            int2d(i,k) = int_sc
!            write(*,*) i, k, x(i), z(k), alo_u+alo_p+alo_d+abs_sc
!             write(*,*)
            
            alocont_o_nn2d(iip1,kkp1,q4) = alo_d*c05_scontd
            alocont_o_nn2d(i,kkp1,q5) = (alo_d*c04_scontd + abs_sc*(c05_intu*alocont_o_nn2d(i,kkp1,q4)))
            alocont_o_nn2d(iim1,kkp1,q6) = (alo_d*c03_scontd + abs_sc*(c05_intu*alocont_o_nn2d(iim1,kkp1,q5)))
            alocont_o_nn2d(iip1,k,q13) = (alo_d*c02_scontd + abs_sc*(c04_intu*alocont_o_nn2d(iip1,k,q4)))
            alocont_o_nn2d(i,k,q14) = (alo_p + abs_sc*((c03_intu*alocont_o_nn2d(i,k,q4) + &
                                                       c04_intu*alocont_o_nn2d(i,k,q5) + &
                                                       c05_intu*alocont_o_nn2d(i,k,q13))))
            alocont_o_nn2d(iim1,k,q15) = (alo_u*c05_scontu + abs_sc*((c02_intu*alocont_o_nn2d(iim1,k,q4) + &
                                                                     c03_intu*alocont_o_nn2d(iim1,k,q5) + &
                                                                     c04_intu*alocont_o_nn2d(iim1,k,q6) + &
                                                                     c05_intu*alocont_o_nn2d(iim1,k,q14))))
            alocont_o_nn2d(iip1,kkm1,q22) = (alo_d*c01_scontd + abs_sc*((c04_intu*alocont_o_nn2d(iip1,kkm1,q13))))
            alocont_o_nn2d(i,kkm1,q23) = (alo_u*c04_scontu + abs_sc*((c03_intu*alocont_o_nn2d(i,kkm1,q13) + &
                                                                     c04_intu*alocont_o_nn2d(i,kkm1,q14) + &
                                                                     c05_intu*alocont_o_nn2d(i,kkm1,q22) + &
                                                                     c01_intu*alocont_o_nn2d(i,kkm1,q4))))
            alocont_o_nn2d(iim1,kkm1,q24) = (alo_u*c03_scontu + abs_sc*((c02_intu*alocont_o_nn2d(iim1,kkm1,q13) + &
                                                                        c03_intu*alocont_o_nn2d(iim1,kkm1,q14) + &
                                                                        c04_intu*alocont_o_nn2d(iim1,kkm1,q15) + &
                                                                        c05_intu*alocont_o_nn2d(iim1,kkm1,q23) + &
                                                                        c01_intu*alocont_o_nn2d(iim1,kkm1,q5))))


!            if(k.eq.18.and.i.eq.nx-2) then
!               write(*,*) i, k, int2d(i,k), int2d(iim1,k), int2d(iim1,kkm1), abs_sc*int_u, abs_sc*(c05_intu*alocont_o_nn2d(iim1,kkp1,q5))
!               write(*,*) i, k, int2d(i,k), scont2d(i+1,k), alocont_o_nn2d(iip1,k,q13), int_u, int2d(i,k)               
!            endif
            
!            write(*,*) i, k, scont2d(i,k), int2d(i,k), int_u, contr_sc, abs_sc!int2d(2,k), int2d(nx,k)
!            if(i.eq.5.and.k.eq.21) then
!               write(*,*) int2d(i,k), alocont_o_nn2d(iip1,kkm1,q22), int_u, c04_intu, &
!                    int2d(i,k), int2d(iim1,k), int2d(i,kkm1), int2d(2,k)
!               stop
!            endif
!            if(i.eq.16.and.k.eq.25) then
!               write(*,*) 'fsc', i, k, alocont_o_nn2d(i,k,q14), int2d(i,k)
!            endif
!            if(i.eq.17.and.k.eq.25) then
!               write(*,*) 'fsc', i, k, alocont_o_nn2d(iim1,k,q15), int2d(i,k)
!            endif
!            if(i.eq.15.and.k.eq.26) then
!               write(*,*) 'fsc', i, k, alocont_o_nn2d(iip1,kkm1,q22), int2d(i,k)
!               write(*,*) abs_sc, int_u, c05_intu, alocont_o_nn2d(iip1,kkm1,q13), c05_intu*alocont_o_nn2d(iip1,kkm1,q13), c05_intu*alocont_o_nn2d(iip1,kkm1,q13)*abs_sc, int2d(i,kkm!1)
!            endif
!            if(k.eq.4) then
!               write(*,*) i, k, nn_z, int2d(i,k), int_u, opac_d, opac2d(iip1,kkp1), (x_d-x(i)),(x(iip1)-x(i)), c04_scontd, c05_scontd
!            endif
!            if(k.eq.3.and.i.eq.16) then
!               write(*,*) i, k, alocont_o_nn2d(i,k,q14), alo_u, alo_p, alo_d, alo_u+alo_p+alo_d+abs_sc, alocont_o_nn2d(i,k,q4), alocont_o_nn2d(i,k,q13)
!               write(*,*) i, k, nn_z, opac_u, opac_p, opac_d, x_u, z_u, x_d, z_d
!               write(*,*)
!            endif
!
!            write(*,*) alocont_o_nn2d(iip1,kkp1,q4), alocont_o_nn2d(i,kkp1,q5), alocont_o_nn2d(iim1,kkp1,q6), alocont_o_nn2d(i,k,q14), alocont_o_nn2d(iip1,kkm1,q22), alocont_o_nn2d(i,kkm1,q23), alocont_o_nn2d(iim1,kkm1,q24)

            mint2d_tmp(i,k) = mint2d_tmp(i,k) + int2d(i,k)*wall
            normalization2d_tmp(i,k) = normalization2d_tmp(i,k) + wall
!        
            alocont_nn2d_tmp(iip1,kkp1,q4) = alocont_nn2d_tmp(iip1,kkp1,q4) + wall*alocont_o_nn2d(iip1,kkp1,q4)
            alocont_nn2d_tmp(i,kkp1,q5) = alocont_nn2d_tmp(i,kkp1,q5) + wall*alocont_o_nn2d(i,kkp1,q5)
            alocont_nn2d_tmp(iim1,kkp1,q6) = alocont_nn2d_tmp(iim1,kkp1,q6) + wall*alocont_o_nn2d(iim1,kkp1,q6)
            alocont_nn2d_tmp(iip1,k,q13) = alocont_nn2d_tmp(iip1,k,q13) + wall*alocont_o_nn2d(iip1,k,q13)
            alocont_nn2d_tmp(i,k,q14) = alocont_nn2d_tmp(i,k,q14) + wall*alocont_o_nn2d(i,k,q14)
            alocont_nn2d_tmp(iim1,k,q15) = alocont_nn2d_tmp(iim1,k,q15) + wall*alocont_o_nn2d(iim1,k,q15)
            alocont_nn2d_tmp(iip1,kkm1,q22) = alocont_nn2d_tmp(iip1,kkm1,q22) + wall*alocont_o_nn2d(iip1,kkm1,q22)
            alocont_nn2d_tmp(i,kkm1,q23) = alocont_nn2d_tmp(i,kkm1,q23) + wall*alocont_o_nn2d(i,kkm1,q23)
            alocont_nn2d_tmp(iim1,kkm1,q24) = alocont_nn2d_tmp(iim1,kkm1,q24) + wall*alocont_o_nn2d(iim1,kkm1,q24)


!
         case default
!
      end select
   enddo
enddo




!stop 'go on in fsc_cont2d'
!
end subroutine fsc_cont2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_cont2d_lin(oindx,nueindx)
!
!-----------------------------------------------------------------------
!------short characteristics for continuum radiative transfer in 2d-----
!---------------calculating intensties for given mu---------------------
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use dime2d, only: nx, nz, x, z, int2d, opac2d, scont2d, alocont_nn2d, imask2d, &
                  mint2d_tmp, normalization2d_tmp, alocont_o_nn2d, alocont_nn2d_tmp, intbound2d
use angles, only: n_x, n_z, weight_omega, q_alo
use freq, only: nodes_nue
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: oindx, nueindx
!
! ... local scalars
integer(i4b) :: i, k
integer(i4b) :: alpha, gamma
integer(i4b) :: startx, startz, endx, endz
integer(i4b) :: startxb, endxb, alphab, gammab
integer(i4b) :: bindxa, bindxb
integer(i4b) :: iim1, ii, iip1,kkm1, kk, kkp1
real(dp) ::  xnue, nn_x, nn_z, wall
real(dp) :: opac_u, scont_u, dels_u, dels_xu, dels_zu, int_u, &
            opac_d, scont_d, dels_d, dels_xd, dels_zd, &
            opac_p, scont_p, dels_r
real(dp) :: x_u, z_u, x_d, z_d
real(dp) :: abs_sc, int_sc, contr_sc, alo_p, alo_u, alo_d
real(dp) :: c03_scontu, c04_scontu, c05_scontu
real(dp) :: c03_intu, c04_intu, c05_intu
real(dp) :: c02_scontd, c04_scontd, c05_scontd
integer :: q14, q15, q23, q24
!
! ... for debugging
!
! ... local arrays
!
! ... local functions
!
! ... local characters
character(len=50) :: enter
!
! ... local logicals
!
!frequency
xnue=nodes_nue(nueindx)

!note: nn_x always positive because of symmetry
nn_x=n_x(oindx)
nn_z=n_z(oindx)
!
wall=weight_omega(oindx)
!
!indices for nearest neighbour alo
q14=q_alo(oindx,14)
q15=q_alo(oindx,15)
q23=q_alo(oindx,23)
q24=q_alo(oindx,24)
!
!-----------------------------------------------------------------------
!
!set directional index-parameter (phantom points are excluded from calculation)
!
if(nn_x.gt.zero) then
   startx = 3
   endx = nx-2
   alpha=  1
   startxb = nx-1
   endxb = 2
   alphab = -1
   bindxa = 1
   bindxb = 2
elseif(nn_x.lt.zero) then
   startx = nx-2
   endx = 3
   alpha=-1
   startxb = 2
   endxb = nx-1
   alphab = 1
   bindxa = 2
   bindxb = 1
else
   stop 'error in fsc_cont2d: n_x = 0 not allowed'
endif
!
if(nn_z.gt.zero) then
   startz = 3
   endz = nz-2
   gamma =  1
   gammab = -1
elseif(nn_z.lt.zero) then
   startz = nz-2
   endz = 3
   gamma = -1
   gammab = 1
else
   stop 'error in fsc_cont2d: n_z = 0 not allowed'
endif
!
!-----------------------reset the intensities---------------------------
!
call set_boundary2d(xnue,nn_x,nn_z)
!
!-----------------------------------------------------------------------
!
alocont_o_nn2d=zero
!
!scont2d=zero
!int2d=zero
!scont2d(nx-2,7)=one

!-----------------------------------------------------------------------
!
!do k=7, 7, gamma
do k=startz, endz, gamma
!
!   write(*,*) k, nz/2+1
   call set_boundaryp2d_lin(nn_x, nn_z, startxb,endxb,alphab,gammab,k,bindxa,bindxb,oindx,wall,q14,q15,q23,q24)
!   call set_boundaryp2d_linb(k,oindx,bindxa,bindxb,startx,alpha)
!   write(*,*)
!   
   do i=startx, endx, alpha

      select case(imask2d(i,k))
!
!*************************standard rt procedure*************************
!
         case(1)

            iim1=i-alpha
            kkm1=k-gamma
            iip1=i+alpha
            kkp1=k+gamma
!
!calculate distance to intersection point with upwind x-grid-coordinate
            dels_xu=(x(i)-x(iim1))/nn_x
!calculate distance to instersection point with upwind z-grid-coordinate
            dels_zu=(z(k)-z(kkm1))/nn_z
            dels_u=min(dels_xu,dels_zu)
!
!calculate distance to intersection point with downwind x-grid-coordinate
            dels_xd=(x(iip1)-x(i))/nn_x
!calculate distance to instersection point with downwind z-grid-coordinate
            dels_zd=(z(kkp1)-z(k))/nn_z
            dels_d=min(dels_xd,dels_zd)
!
!----------------------------local point--------------------------------
!
            opac_p=opac2d(i,k)
            scont_p=scont2d(i,k)
!
!----------------------------upwind point-------------------------------
!
            if(dels_xu.eq.dels_u) then
!intersection with z-line on x-level i-alpha
               x_u = x(i) - dels_u*nn_x
               z_u = z(k) - dels_u*nn_z
!               
               call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(iim1,k), &
                                      scont2d(iim1,kkm1), scont2d(iim1,k), &
                                      int2d(iim1,kkm1), int2d(iim1,k), &
                                      z(kkm1), z(k), z_u, &
                                      c03_scontu, c05_scontu, &
                                      c03_intu, c05_intu, &
                                      opac_u, scont_u, int_u)
               c04_scontu=zero
               c04_intu=zero
               
            elseif(dels_zu.eq.dels_u) then
!intersection with x-line on z-level k-gamma
               x_u=x(i)-dels_u*nn_x
               z_u=z(k)-dels_u*nn_z
!
               call coeff1d_contu_lin(opac2d(iim1,kkm1), opac2d(i,kkm1), &
                                      scont2d(iim1,kkm1), scont2d(i,kkm1), &
                                      int2d(iim1,kkm1), int2d(i,kkm1), &
                                      x(iim1), x(i), x_u, &
                                      c03_scontu, c04_scontu, &
                                      c03_intu, c04_intu, &
                                      opac_u, scont_u, int_u)
               c05_scontu=zero
               c05_intu=zero                
            else
               write(*,'(3es20.8,2l4)') dels_u, dels_xu, dels_zu, dels_u.eq.dels_xu, dels_u.eq.dels_zu
               stop 'error in fsc_cont2d: invalid dels_u'
            endif
!
!--------------------------downwind point-------------------------------
!
            if(dels_xd.eq.dels_d) then
!intersection with z-line on x-level i+alpha
               x_d=x(i)+dels_d*nn_x
               z_d=z(k)+dels_d*nn_z
!
               call coeff1d_contd_lin(opac2d(iip1,k), opac2d(iip1,kkp1), &
                                      scont2d(iip1,k), scont2d(iip1,kkp1), &
                                      z(k), z(kkp1), z_d, &
                                      c02_scontd, c05_scontd, &
                                      opac_d, scont_d)
               c04_scontd=zero
            elseif(dels_zd.eq.dels_d) then
!intersection with x-line on z-level k+gamma
               x_d=x(i)+dels_d*nn_x
               z_d=z(k)+dels_d*nn_z
!
               call coeff1d_contd_lin(opac2d(i,kkp1), opac2d(iip1,kkp1), &
                                      scont2d(i,kkp1), scont2d(iip1,kkp1), &
                                      x(i), x(iip1), x_d, &
                                      c04_scontd, c05_scontd, &
                                      opac_d, scont_d)
               c02_scontd=zero

            else
               write(*,'(3es20.8,2l4)') dels_d, dels_xd, dels_zd, dels_d.eq.dels_xd, dels_d.eq.dels_zd
               stop 'error in fsc_cont2d: invalid dels_d'
            endif
!
!-----------------------radiative transfer------------------------------
!
            call fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                                dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
            int2d(i,k) = int_sc
!            write(*,*) i, k, int_sc, int_u, scont_u, scont_p, scont_d

!            alocont_o_nn2d(i,k,q14) = imask2d(i,k)*alo_p
!            alocont_o_nn2d(iim1,k,q15) = imask2d(iim1,k) * &
!                                              (alo_u*c05_scontu + abs_sc*(c05_intu*alocont_o_nn2d(iim1,k,q14)))
!            alocont_o_nn2d(i,kkm1,q23) = imask2d(i,kkm1) * &
!                                              (alo_u*c04_scontu + abs_sc*(c04_intu*alocont_o_nn2d(i,kkm1,q14)))
!            alocont_o_nn2d(iim1,kkm1,q24) = imask2d(iim1,kkm1) * &
!                                                    (alo_u*c03_scontu + abs_sc*(c03_intu*alocont_o_nn2d(iim1,kkm1,q14) + &
!                                                                                c04_intu*alocont_o_nn2d(iim1,kkm1,q15) + &
!                                                                                c05_intu*alocont_o_nn2d(iim1,kkm1,q23)))
            alocont_o_nn2d(i,k,q14) = alo_p
            alocont_o_nn2d(iim1,k,q15) = alo_u*c05_scontu + abs_sc*(c05_intu*alocont_o_nn2d(iim1,k,q14))
            alocont_o_nn2d(i,kkm1,q23) = alo_u*c04_scontu + abs_sc*(c04_intu*alocont_o_nn2d(i,kkm1,q14))
            alocont_o_nn2d(iim1,kkm1,q24) = alo_u*c03_scontu + abs_sc*(c03_intu*alocont_o_nn2d(iim1,kkm1,q14) + &
                                                                       c04_intu*alocont_o_nn2d(iim1,kkm1,q15) + &
                                                                       c05_intu*alocont_o_nn2d(iim1,kkm1,q23))

            
!            write(*,'(2i5,13es16.8)') i, k, alocont_o_nn2d(i,k,q14), opac_u, opac_p, opac_d, dels_u, dels_d, alo_u, alo_p, abs_sc
!            write(*,*)
!            write(*,'(i5,13es16.8)') k, opac_u, opac_p, opac_d, alo_u, alo_p, alo_d, dels_u, dels_d
            
!perform angular integration
            mint2d_tmp(i,k) = mint2d_tmp(i,k) + int2d(i,k)*wall
            normalization2d_tmp(i,k) = normalization2d_tmp(i,k) + wall
!
            alocont_nn2d_tmp(i,k,q14) = alocont_nn2d_tmp(i,k,q14) + wall*alocont_o_nn2d(i,k,q14)
            alocont_nn2d_tmp(iim1,k,q15) = alocont_nn2d_tmp(iim1,k,q15) + wall*alocont_o_nn2d(iim1,k,q15)
            alocont_nn2d_tmp(i,kkm1,q23) = alocont_nn2d_tmp(i,kkm1,q23) + wall*alocont_o_nn2d(i,kkm1,q23)
            alocont_nn2d_tmp(iim1,kkm1,q24) = alocont_nn2d_tmp(iim1,kkm1,q24) + wall*alocont_o_nn2d(iim1,kkm1,q24)
!            write(*,'(i5, 13es20.8)') k, zero, alocont_nn2d_tmp(i,k,q14), alocont_nn2d_tmp(i,kkm1,q23)
!
         case default
!
      end select
   enddo
!
!calculate boundary intensities and set them to the beginning
!   call calc_boundaryp2d_linb(startx,endx,alpha,gamma,k,oindx,bindxa,bindxb,nn_x,nn_z,wall,q14)
enddo

!write(*,*)
!do k=1, nz
   !   write(*,*) k, alocont_o_nn2d(1,k,q14), alocont_o_nn2d(2,k,q14), alocont_o_nn2d(nx-1,k,q14), alocont_o_nn2d(nx,k,q14), alocont_o_nn2d(nx/2,k,q14), alocont_nn2d_tmp(1,k,q14), alocont_nn2d_tmp(nx/2+1,k,q14)
!   write(*,*) k, nn_z, int2d(nx/2+1,k), scont2d(nx/2+1,k)
!enddo
!stop 'go on in fsc_cont2d_lin'
!
end subroutine fsc_cont2d_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeffcr1d_mbez(fc, f_im1, f_i, at, bt, ct, at2, bt2, ct2)
!
!calculates new coefficients for 2d bezier interpolation to ensure monotonicity
!                 when interpolating in right interval
!
!on input: 
!   fc           - control point
!   f_im1, f_i   - function values that limit the interval
!   at, bt, ct   - coefficients that have been used for calculation of control point
!
!on output:
!   at2, bt2, ct2 - interpolation coefficients for control point sucht that monotonic interpolation can be used
!                   (e.g. at2=0.,bt2=0.,ct2=1. if f_c > f_i > f_im1)
!
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: fc, f_im1, f_i, at, bt, ct
real(dp), intent(out) :: at2, bt2, ct2
!
!--------------------------old version----------------------------------
!
!if(f_i.ge.f_im1) then
!   if(fc.gt.f_i) then
!      at2=zero
!      bt2=zero
!      ct2=one
!   elseif(fc.lt.f_im1) then
!      at2=zero
!      bt2=one
!      ct2=zero
!   else
!      at2=at
!      bt2=bt
!      ct2=ct
!   endif
!elseif(f_i.le.f_im1) then
!   if(fc.lt.f_i) then
!      at2=zero
!      bt2=zero
!      ct2=one
!   elseif(fc.gt.f_im1) then
!      at2=zero
!      bt2=one
!      ct2=zero
!   else
!      at2=at
!      bt2=bt
!      ct2=ct
!   endif
!endif
!
!--------------------------new version----------------------------------
!
if((f_i-f_im1)*(f_i-fc).le.zero) then
   at2=zero
   bt2=zero
   ct2=one
elseif((f_i-f_im1)*(fc-f_im1).le.zero) then
   at2=zero
   bt2=one
   ct2=zero
else
   at2=at
   bt2=bt
   ct2=ct
endif
!
end subroutine coeffcr1d_mbez
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff1d_contu(opac_im2, opac_im1, opac_i, &
                         scont_im2, scont_im1, scont_i, &
                         int_im2, int_im1, int_i, &
                         x_im2, x_im1, x_i, x_p, &
                         a_scont, b_scont, c_scont, &
                         a_inten, b_inten, c_inten, &
                         opac_p, scont_p, int_p)
!
!         interpolates opacity, continuum source function and intensity
!               values given on a 1d grid onto point x_p
!
!on input (f_* stands for opac_*, scont_* and int_*, respectivly):
!

!         f_im2-------f_im1----x-------f_i
!         x_im2-------x_im1---x_p------x_i
!                           
!        x_p coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_scont, b_scont, c_scont
!         a_inten, b_inten, c_inten
!
!      such that:
!         f_p = a*f_im2 + b*f_im1 + c*f_i
!
!   2. interpolated values at point p: opac_p, scont_p, int_p
!
use prog_type
use mod_interp2d, only: wp_interp2d
use fund_const
use dime2d, only: zmin, zmax
use mod_benchmark, only: tau_min, tau_max
!
implicit none
!
! ... arguments
real(dp), intent(in) :: opac_im2, opac_im1, opac_i, &
                        scont_im2, scont_im1, scont_i, &
                        int_im2, int_im1, int_i, &
                        x_im2, x_im1, x_i, x_p
real(dp), intent(out) :: a_scont, b_scont, c_scont, &
                         a_inten, b_inten, c_inten, &
                         opac_p, int_p, scont_p
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, tx, ax, bx, cx, &
            axt, bxt, cxt, axt2, bxt2, cxt2, &
            ax_opac, bx_opac, cx_opac, &
            axt_opac, bxt_opac, cxt_opac, &
            opacc, &
            axt_scont, bxt_scont, cxt_scont, &
            scontc, &
            axt_int, bxt_int, cxt_int, &
            intc, &
            dxp
real(dp) :: fac, fac2
!
! ... for debugging
real(dp) :: aopac, bopac
!
!
!
!define deltax, deltay
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
!
!define deltax-ratio
tx = (x_p-x_im1)/dxi
!
!20: standard quadratic interpolation
!50: control point from derivative (monotonic, with predefined weights for derivative)
goto 20
!
!-------------standard (non-monotonic) quadratic interpolation----------
!
!
20 continue
!
dxp = x_p-x_im1
!
a_scont = dxp*(dxp-dxi)/dxim1/dx
b_scont = (dxp+dxim1)*(dxi-dxp)/dxim1/dxi
c_scont = dxp*(dxp+dxim1)/dxi/dx
!
opac_p = a_scont*opac_im2 + b_scont*opac_im1 + c_scont*opac_i
scont_p = a_scont*scont_im2 + b_scont*scont_im1 + c_scont*scont_i
!
a_inten = a_scont
b_inten = b_scont
c_inten = c_scont
int_p = a_inten*int_im2 + b_inten*int_im1 + c_inten*int_i
!
!
!better: use linear interpolations for intensity and source function
a_scont=0.d0
b_scont=one-tx
c_scont=tx
scont_p = a_scont*scont_im2 + b_scont*scont_im1 + c_scont*scont_i
!
a_inten = a_scont
b_inten = b_scont
c_inten = c_scont
int_p = a_inten*int_im2 + b_inten*int_im1 + c_inten*int_i
!
opac_p = a_scont*opac_im2 + b_scont*opac_im1 + c_scont*opac_i
!
!analytic opacity law
!bopac=log(tau_max/tau_min)/(zmax-zmin)
!aopac=bopac*tau_max*exp(bopac*zmin)
!opac_p = aopac*exp(-bopac*x_p)

return
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!----------------and weights for derivative assigned--------------------
!
50 continue
!
ax = (one-tx)**2
bx = two*tx*(one-tx)
cx = tx**2
!
!calculate control points
!   scontc, scontc, scontc
!   intc, intc, intc
!   opacc, opacc, opacc
!
fac=max(wp_interp2d,dxim1/dx)
fac2=dxim1/dx
axt = (fac-one)*dxi/two/dxim1
bxt = ((two-fac)*dxim1 + (one-fac)*dxi)/two/dxim1
cxt = fac/two
axt2 = (fac2-one)*dxi/two/dxim1
bxt2 = ((two-fac2)*dxim1 + (one-fac2)*dxi)/two/dxim1
cxt2 = fac2/two
!
opacc   = opac_im2*axt2   + opac_im1*bxt2   + opac_i*cxt2
scontc  = scont_im2*axt   + scont_im1*bxt   + scont_i*cxt
intc    = int_im2*axt     + int_im1*bxt     + int_i*cxt
!
!ensure monotonicity
call coeffcr1d_mbez(opacc, opac_im1, opac_i, axt2, bxt2, cxt2, axt_opac, bxt_opac, cxt_opac)
call coeffcr1d_mbez(scontc, scont_im1, scont_i, axt, bxt, cxt, axt_scont, bxt_scont, cxt_scont)
call coeffcr1d_mbez(intc, int_im1, int_i, axt, bxt, cxt, axt_int, bxt_int, cxt_int)
ax_opac = axt_opac*bx
bx_opac = bxt_opac*bx + ax
cx_opac = cxt_opac*bx + cx
opac_p = ax_opac*opac_im2 + bx_opac*opac_im1 + cx_opac*opac_i
a_scont = axt_scont*bx
b_scont = bxt_scont*bx + ax
c_scont = cxt_scont*bx + cx
scont_p = a_scont*scont_im2 + b_scont*scont_im1 + c_scont*scont_i
a_inten = axt_int*bx
b_inten = bxt_int*bx + ax
c_inten = cxt_int*bx + cx
int_p = a_inten*int_im2 + b_inten*int_im1 + c_inten*int_i
!


!to test: use linear interpolations for source function and intensity
!a_scont=0.d0
!b_scont=one-tx
!c_scont=tx
!scont_p = a_scont*scont_im2 + b_scont*scont_im1 + c_scont*scont_i
!!!
!a_inten = a_scont
!b_inten = b_scont
!c_inten = c_scont
!int_p = a_inten*int_im2 + b_inten*int_im1 + c_inten*int_i

!write(*,*) x_im2, x_im1, (x_i+x_im1)/two, x_i
!write(*,*) int_im2, int_im1, intc, int_i
!write(*,*) int_im2, int_im1, (int_im1+int_i)/two, int_i
!write(*,*) opac_p
!write(*,*)
!
return
!
end subroutine coeff1d_contu
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff1d_contd(opac_im2, opac_im1, opac_i, &
                         scont_im2, scont_im1, scont_i, &
                         x_im2, x_im1, x_i, x_p, &
                         a_scont, b_scont, c_scont, &
                         opac_p, scont_p)
!
!            interpolates opacity and continuum source function
!               values given on a 1d grid onto point x_p
!
!on input (f_* stands for opac_* and scont_* respectivly):
!
!    f_im2-------f_im1-----x---------f_i
!    x_im2-------x_im1----x_p--------x_i
!                           
!        x_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function
!      (required for ALO calculations):
!         a_scont, b_scont, c_scont
!
!      such that:
!         f_p = a*f_im2 + b*f_im1 + c*f_i

!   2. interpolated values at point p: opac_p, scont_p
!
use prog_type
use mod_interp2d, only: wp_interp2d
use fund_const
use dime2d, only: zmin, zmax
use mod_benchmark, only: tau_min, tau_max
!
implicit none
!
! ... arguments
real(dp), intent(in) :: opac_im2, opac_im1, opac_i, &
                        scont_im2,   scont_im1,   scont_i, &
                        x_im2, x_im1, x_i, x_p
real(dp), intent(out) :: a_scont, b_scont, c_scont, opac_p, scont_p
!
! ... local scalars
real(dp) :: dxim1, dxi, dx, dxp, tx, ax, bx, cx, &
            axt, bxt, cxt, axt2, bxt2, cxt2, &
            ax_opac, bx_opac, cx_opac, &
            axt_opac, bxt_opac, cxt_opac, &
            opacc, &
            axt_scont, bxt_scont, cxt_scont, &
            scontc
real(dp) :: fac, fac2
!
! ... for debugging
real(dp) :: aopac, bopac

!
!define deltax
dxim1 = x_im1-x_im2
dxi = x_i-x_im1
dx = dxim1+dxi
!
!define deltax-ratio
tx = (x_p-x_im1)/dxi

!
!20: standard quadratic interpolation
!50: control point from derivative (monotonic, with predefined weights for derivative)
goto 20
!
!-------------standard (non-monotonic) quadratic interpolation----------
!
!
20 continue
!
dxp = x_p-x_im1
!
a_scont = dxp*(dxp-dxi)/dxim1/dx
b_scont = (dxp+dxim1)*(dxi-dxp)/dxim1/dxi
c_scont = dxp*(dxp+dxim1)/dxi/dx
!
opac_p = a_scont*opac_im2 + b_scont*opac_im1 + c_scont*opac_i
scont_p = a_scont*scont_im2 + b_scont*scont_im1 + c_scont*scont_i
!
!
!better: use linear interpolations for source function
a_scont=0.d0
b_scont=one-tx
c_scont=tx
scont_p = a_scont*scont_im2 + b_scont*scont_im1 + c_scont*scont_i
opac_p = a_scont*opac_im2 + b_scont*opac_im1 + c_scont*opac_i
!
!
!analytic opacity law
!bopac=log(tau_max/tau_min)/(zmax-zmin)
!aopac=bopac*tau_max*exp(bopac*zmin)
!opac_p = aopac*exp(-bopac*x_p)
!
!

return
!
!---------------monotonic quadratic bezier interpolation----------------
!--------(with derivatives for control point from weighted mean)--------
!----------------and weights for derivative assigned--------------------
!
50 continue
!
ax = (one-tx)**2
bx = two*tx*(one-tx)
cx = tx**2
!
!calculate control points on each j-level:
!   scontc_jm2, scontc_jm1, scontc_j
!   opacc_jm2, opacc_jm1, opacc_j
fac=max(wp_interp2d,dxim1/dx)
fac2=dxim1/dx
axt = (fac-one)*dxi/two/dxim1
bxt = ((two-fac)*dxim1 + (one-fac)*dxi)/two/dxim1
cxt = fac/two
axt2 = (fac2-one)*dxi/two/dxim1
bxt2 = ((two-fac2)*dxim1 + (one-fac2)*dxi)/two/dxim1
cxt2 = fac2/two
!
opacc   = opac_im2*axt2   + opac_im1*bxt2   + opac_i*cxt2
scontc  = scont_im2*axt   + scont_im1*bxt   + scont_i*cxt
!
!ensure monotonicity
call coeffcr1d_mbez(opacc, opac_im1, opac_i, axt2, bxt2, cxt2, axt_opac, bxt_opac, cxt_opac)
call coeffcr1d_mbez(scontc, scont_im1, scont_i, axt, bxt, cxt, axt_scont, bxt_scont, cxt_scont)
ax_opac = axt_opac*bx
bx_opac = bxt_opac*bx + ax
cx_opac = cxt_opac*bx + cx
opac_p = ax_opac*opac_im2 + bx_opac*opac_im1 + cx_opac*opac_i
a_scont = axt_scont*bx
b_scont = bxt_scont*bx + ax
c_scont = cxt_scont*bx + cx
scont_p = a_scont*scont_im2 + b_scont*scont_im1 + c_scont*scont_i
!
!
end subroutine coeff1d_contd
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff1d_contu_lin(opac_im1, opac_i, &
                             scont_im1, scont_i, &
                             int_im1, int_i, &
                             x_im1, x_i, x_p, &
                             a_scont, b_scont, &
                             a_inten, b_inten, &
                             opac_p, scont_p, int_p)
!
!         interpolates opacity, continuum source function and intensity
!               values given on a 1d grid onto point x_p
!
!on input (f_* stands for opac_*, scont_* and int_*, respectivly):
!
!    f_im1------x--------f_i
!    x_im1-----x_p-------x_i
!                           
!        x_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_scont, b_scont
!         a_inten, b_inten
!
!      such that:
!         f_p = a*f_im1 + b*f_i
!
!   2. interpolated values at point p: opac_p, scont_p, int_p
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: opac_im1, opac_i, scont_im1, scont_i, &
                        int_im1, int_i, x_im1, x_i, x_p
real(dp), intent(out) :: a_scont, b_scont, a_inten, b_inten, &
                         opac_p, int_p, scont_p
!
! ... local scalars
real(dp) :: dxi, tx
!
!
!define deltax
dxi = x_i-x_im1
!
!define deltax-ratio
tx = (x_p-x_im1)/dxi
!
!-------------------------linear interpolation------------------------
!
a_scont=one-tx
b_scont=tx
!
a_inten=a_scont
b_inten=b_scont
!
opac_p  = a_scont*opac_im1  + b_scont*opac_i
scont_p = a_scont*scont_im1 + b_scont*scont_i
int_p   = a_scont*int_im1   + b_scont*int_i
return
!
end subroutine coeff1d_contu_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff1d_contd_lin(opac_im1, opac_i, &
                             scont_im1, scont_i, &
                             x_im1, x_i, x_p, &
                             a_scont, b_scont, &
                             opac_p, scont_p)
!
!         interpolates opacity, continuum source function and intensity
!               values given on a 1d grid onto point x_p
!
!on input (f_* stands for opac_*, scont_* and int_*, respectivly):
!
!    f_im1------x--------f_i
!    x_im1-----x_p-------x_i
!                           
!        x_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_scont, b_scont
!         a_inten, b_inten
!
!      such that:
!         f_p = a*f_im1 + b*f_i
!
!   2. interpolated values at point p: opac_p, scont_p, int_p
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: opac_im1, opac_i, scont_im1, scont_i, &
                        x_im1, x_i, x_p
real(dp), intent(out) :: a_scont, b_scont, opac_p, scont_p
!
! ... local scalars
real(dp) :: dxi, tx
!
!
!define deltax
dxi = x_i-x_im1
!
!define deltax-ratio
tx = (x_p-x_im1)/dxi
!
!-------------------------linear interpolation------------------------
!
a_scont=one-tx
b_scont=tx
!
opac_p  = a_scont*opac_im1  + b_scont*opac_i
scont_p = a_scont*scont_im1 + b_scont*scont_i
!
return
!
end subroutine coeff1d_contd_lin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine coeff1d_contu_linb(opac_im1, opac_i, &
                             scont_im1, scont_i, &
                             int_im1, int_i, &
                             x_im1, x_i, x_p, &
                             a_scont, b_scont, &
                             a_inten, b_inten, &
                             opac_p, scont_p, int_p)
!
!         interpolates opacity (in log-scale), continuum source function and intensity
!               values given on a 1d grid onto point x_p
!
!on input (f_* stands for opac_*, scont_* and int_*, respectivly):
!
!    f_im1------x--------f_i
!    x_im1-----x_p-------x_i
!                           
!        x_p: coordinates of point onto which shall be interpolated
!
!on output:
!   1. interpolation coefficients for source function and intensity
!      (required for ALO calculations):
!         a_scont, b_scont
!         a_inten, b_inten
!
!      such that:
!         f_p = a*f_im1 + b*f_i
!
!   2. interpolated values at point p: opac_p, scont_p, int_p
!
use prog_type
use fund_const
!
implicit none
!
! ... arguments
real(dp), intent(in) :: opac_im1, opac_i, scont_im1, scont_i, &
                        int_im1, int_i, x_im1, x_i, x_p
real(dp), intent(out) :: a_scont, b_scont, a_inten, b_inten, &
                         opac_p, int_p, scont_p
!
! ... local scalars
real(dp) :: dxi, tx, x_im1t, x_it, x_pt, xshift, a_opac, b_opac
!
!
!define deltax
dxi = x_i-x_im1
!
!define deltax-ratio
tx = (x_p-x_im1)/dxi
!
!---------------linear interpolation of intensity and source function-------------
!
a_scont=one-tx
b_scont=tx
!
a_inten=a_scont
b_inten=b_scont
!
scont_p = a_scont*scont_im1 + b_scont*scont_i
int_p   = a_scont*int_im1   + b_scont*int_i
!
!-----------------------log-log interpolation of opacity---------------------------
!
!transform to a minimum x-value of 1 (to avoid negative x-values in log-log interpolation)
xshift=one-min(x_im1,x_i)
x_im1t = log10(x_im1+xshift)
x_it = log10(x_i+xshift)
x_pt = log10(x_p+xshift)
!
dxi = x_it-x_im1t
tx = (x_pt-x_im1t)/dxi
!
a_opac = one-tx
b_opac = tx
!
opac_p = opac_im1**a_opac * opac_i**b_opac
!
!write(*,*) a_scont*opac_im1+b_scont*opac_i, opac_p, a_opac, b_opac
!write(*,*) x_im1, x_i, x_p, xshift
!write(*,*) x_im1t, x_it, x_pt
!write(*,*) 
!stop

return
!
end subroutine coeff1d_contu_linb
