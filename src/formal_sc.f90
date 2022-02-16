!
!***********************************************************************
!***********************************************************************
!
!                     CONTINUUM ROUTINES
!
!***********************************************************************
!***********************************************************************
!
subroutine fsc_cont(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                      dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for continuum transport
!
!             u------------p----------d
!                 dels_u      dels_d
!
!      bezier interpolations/integrations
!
!   input: upwind intensity at u:       int_u
!          continuum opacity at u:      opac_u
!          continuum source-fct at u:   scont_u
!          continuum opacity at p:      opac_p
!          continuum source-fct at p:   scont_p
!          continuum opacity at d:      opac_d
!          continuum source-fct at d:   scont_d
!          path-length from u to p:     dels_u
!          path-length from p to d:     dels_d
!          coordinates of point u:      x_u, z_u
!
!   output: absorption part from u to p:       abs_sc
!           source contribution from u to p:   contr_sc
!           intensity at point p:              int_sc
!           alo_u, alo_p, alo_d:               alo-coefficients for points u, p, d
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use mod_integ1d, only: integ1d_tau_ud, coeff_source2
!
implicit none
!
! ... arguments
real(dp) :: int_u, opac_p, scont_p, &
                        opac_u, scont_u, dels_u, &
                        opac_d, scont_d, dels_d
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p, alo_d
!
! ... local scalars
integer(i4b) :: i, nd
integer(i8b) :: alo_u_int, alo_p_int, alo_d_int, abs_sc_int, norm_int
real(sp) :: alo_u_sp, alo_p_sp, alo_d_sp, abs_sc_sp
real(dp) :: delt_u, delt_d, delt_ud, rerr_contr
real(dp) :: s1, s2, s3, s4, s5, r1, r2, r3, r4, r5, scont2, scont3, scont4, &
            opac2, opac3, opac4, grad, dum, h, f1, f2, f3, f4, f5
real(dp) :: s_iim1, s_ii, opac_iim1, opac_ii, scont_iim1, scont_ii, int_iim1, int_ii, &
            ru, rp, r_ii
real(dp) :: delt_u2, delt_u3, int_sc3, int_sc2
real(dp) :: velr, rho, c1, c2, opac_u2
real(dp) :: tau_u, tau_p, tau_d, tau_iim1, tau_ii, tau_iip1, scont_iip1
real(dp) :: e0, e1, e2, a, b, c
real(dp) :: norm
!
! ... for debug: analytic opacities
real(dp) :: bvel, opac_thomson2
real(dp) :: bconst, xmloss_cgs, alpha_hamann, kappa_hamann, &
            x_p, z_p, x_d, z_d, &
            r_u, r_p, r_d, &
            vel_u, vel_p, vel_d, &
            rho_u, rho_p, rho_d

! ... local functions
real(dp) :: integral0, integral1, integral2, integral2b, integral3, integral3b, &
            integral_source0, integral_source1, integral_source2, integral_source3
real(dp) :: interpol_ypl, interpol_typ_quad3
!
!calculate delta-tau steps via monotonic bezier interpolation
call integ1d_tau_ud(abs(dels_u), abs(dels_d), opac_u, opac_p, opac_d, delt_u, delt_d)
!
abs_sc = exp(-delt_u)
if(delt_u.gt.20.) abs_sc=zero

!
!call coeff_source1(delt_u, alo_u, alo_p)
!call coeff_source4(delt_u, alo_u, alo_p)
!alo_d=zero
!
call coeff_source2(delt_u, delt_d, scont_u, scont_p, scont_d, alo_u, alo_p, alo_d)   !bezier interpolation
!call coeff_source3(delt_u, delt_d, scont_u, scont_p, scont_d, alo_u, alo_p, alo_d)   !warning: might cause osciallations
!call coeff_source3b(delt_u, delt_d, scont_u, scont_p, scont_d, alo_u, alo_p, alo_d)   !warning: might cause osciallations


alo_p=one-alo_u-alo_d-abs_sc
norm=alo_u+alo_p+alo_d+abs_sc


!   write(*,*) 'testa', norm, abs_sc, alo_u, alo_p, alo_d, delt_u, delt_d, dels_u, dels_d
!
if(norm.gt.one) then
!   alo_u_sp = alo_u
!   alo_p_sp = alo_p
!   alo_d_sp = alo_d
!   abs_sc_sp = abs_sc
!
!   alo_u_int = nint(alo_u_sp*1.d13, i8b)
!   alo_p_int = nint(alo_p_sp*1.d13, i8b)
!   alo_d_int = nint(alo_d_sp*1.d13, i8b)
!   abs_sc_int = nint(abs_sc_sp*1.d13, i8b)
!!
!   norm_int=alo_u_int+alo_p_int+alo_d_int+abs_sc_int
!!
!   alo_u = dble(alo_u_int)/dble(norm_int)
!   alo_p = dble(alo_p_int)/dble(norm_int)
!   alo_d = dble(alo_d_int)/dble(norm_int)
!   abs_sc = dble(abs_sc_int)/dble(norm_int)
!   norm=alo_u+alo_p+alo_d+abs_sc

!   write(*,*) 'testa', norm, abs_sc, alo_u, alo_p, alo_d, delt_u, delt_d
!   if(norm.gt.one) stop 'error in fsc_cont: alo coefficients too large'
endif
!
!write(*,*) 'test', opac_u, opac_p, opac_d, dels_u, dels_d, delt_u, delt_d
!
contr_sc = alo_u*scont_u + alo_p*scont_p + alo_d*scont_d
!
int_sc = int_u*abs_sc + contr_sc


!write(*,*) 'fsc_cont', delt_u, delt_d, opac_u, opac_p, opac_d, norm, contr_sc

return
!
end subroutine fsc_cont
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
subroutine fsc_cont_lin(int_u, opac_u, opac_p, opac_d, scont_u, scont_p, scont_d, &
                          dels_u, dels_d, abs_sc, contr_sc, int_sc, alo_u, alo_p)
!
!-----------------------------------------------------------------------
!
!     calculates formal solution along a line segment u-p, to obtain:
!           absorption-part     at point p
!           source-contribution at point p
!           intensity           at point p
!           for continuum transport
!
!     linear interpolations/integrations
!
!             u------------p---------d
!                 dels_u      dels_d
!
!   input: upwind intensity at u:       int_u
!          continuum opacity at u:      opac_u
!          continuum source-fct at u:   scont_u
!          continuum opacity at p:      opac_p
!          continuum source-fct at p:   scont_p
!          continuum opacity at d:      opac_d
!          continuum source-fct at d:   scont_d
!          path-length from p to d:     dels_d
!          path-length from u to p:     dels_u
!
!   output: absorption part from u to p:       abs_sc
!           source contribution from u to p:   contr_sc
!           intensity at point p:              int_sc
!           alo_u, alo_p, alo_d:               alo-coefficients for points u, p, d
!
!-----------------------------------------------------------------------
!
use prog_type
use fund_const
use mod_integ1d, only: integ1d_tau_ud, coeff_source1
!
implicit none
!
! ... arguments
real(dp), intent(in) :: int_u, opac_p, scont_p, &
                        opac_u, scont_u, dels_u, &
                        opac_d, scont_d, dels_d
real(dp), intent(out) :: abs_sc, contr_sc, int_sc, alo_u, alo_p
!
! ... local scalars
integer(i4b) :: i, nd
integer(i8b) :: alo_u_int, alo_p_int, abs_sc_int, norm_int
real(sp) :: alo_u_sp, alo_p_sp, abs_sc_sp
real(dp) :: delt_u, delt_d, norm
!
! ... local functions
real(dp) :: integral0, integral1, integral2, integral2b, integral3, integral3b, &
            integral_source0, integral_source1, integral_source2, integral_source3
!
!
!write(*,'(5es20.8)') dels_u, dels_d, opac_u, opac_p, opac_d
!calculate delta-tau steps via monotonic bezier interpolation
call integ1d_tau_ud(abs(dels_u), abs(dels_d), opac_u, opac_p, opac_d, delt_u, delt_d)
!
!write(*,*) delt_u, delt_d
!write(*,*)
abs_sc = exp(-delt_u)
if(delt_u.gt.20.) abs_sc=zero
!
call coeff_source1(delt_u, alo_u, alo_p)
!
!norm=alo_u+alo_p+abs_sc
!
!if(norm.gt.one) then
!   alo_u_sp = alo_u
!   alo_p_sp = alo_p
!   abs_sc_sp = abs_sc
!
!   alo_u_int = nint(alo_u_sp*1.d13, i8b)
!   alo_p_int = nint(alo_p_sp*1.d13, i8b)
!   abs_sc_int = nint(abs_sc_sp*1.d13, i8b)
!!
!   norm_int=alo_u_int+alo_p_int+abs_sc_int
!!
!   alo_u = dble(alo_u_int)/dble(norm_int)
!   alo_p = dble(alo_p_int)/dble(norm_int)
!   abs_sc = dble(abs_sc_int)/dble(norm_int)
!   norm=alo_u+alo_p+abs_sc
!   write(*,*) 'testa', alo_u+alo_p+alo_d+abs_sc, abs_sc, alo_u, alo_p, alo_d, delt_u, delt_d
!   if(norm.gt.one) stop 'error in fsc_cont: alo coefficients too large'!
!endif


contr_sc = alo_u*scont_u + alo_p*scont_p
int_sc = int_u*abs_sc + contr_sc
!
!
end subroutine fsc_cont_lin
