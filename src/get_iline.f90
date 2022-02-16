subroutine get_iline(iline, xnue0, na, gl, gu, flu)
!
!get xnue0 and na for a given line descirbed by integer iline
!note: whenever you change something in here,
!      you also need to change corresponding entries in
!         function get_opalbar
!         subroutines get_photprof, get_photprof_herrero, calc_photprof
!  
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: iline
integer(i4b), intent(out) :: na
real(dp), intent(out) :: xnue0, gl, gu, flu
!
write(*,*) '--------get_iline: getting xnue0, na, gl, gu, flu-------------'
write(*,'(a10,i20)') 'iline', iline
!
select case(iline)
   case(1)
!hydrogen 2->3
      write(*,*) 'Halpha recombination model'
      xnue0= 4.5680294d14
      gl = 8.d0
      gu = 18.d0
      flu = 6.4108d-1      
      na=1
   case(2)
!hydrogen 2->4
      write(*,*) 'Hbeta'
      xnue0= 6.1668776d14
      gl = 8.d0
      gu = 32.d0
      flu = 1.1938d-1
      na=1      
   case(10,11)
!CIV (1s2,2s) -> (1s2,2p,spin 3/2)
      write(*,*) 'CIV resonance line'
      xnue0 = 1.93798d15
      gl = 2.d0
      gu = 4.d0
      flu = 1.9d-1
      na=12
   case default
      stop 'error in get_iline: iline not properly specified'      
end select
!
write(*,'(a10,es20.8)') 'xnue0', xnue0
write(*,'(a10,es20.8)') 'gl', gl
write(*,'(a10,es20.8)') 'gu', gu
write(*,'(a10,es20.8)') 'flu', flu
write(*,'(a10,es20.8)') 'gf-value', gl*flu
write(*,'(a10,i20)') 'na', na
write(*,*)

!
end subroutine get_iline
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function get_opalbar(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)
!
!get opalbar for a given line transition described by iline
!opacity law for balmer lines by assuming completely ionized hydrogen
!
!input:  sr   stellar radius in cgs
!        yhe, hei:  helium abundance and number of free electrons per helium
!        temp:  temperature in K
!        vth_fiducial:  fiducial thermal velocity in cm/s
!        xnue0:   transition frequency
!        bl,bu:   NLTE-departure coefficients for lower and upper level
!        rho:    density in cgs
!        kline:   just a scaling factor to increase or decrease opacity (line strength)
!  
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: iline
real(dp), intent(in) :: sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
real(dp) :: get_opalbar
!
! ... local functions
real(dp) :: opalbar_halpha, opalbar_hbeta, opalbar_model_kline, opalbar_model_kline2
select case(iline)
   case(1)
      get_opalbar = opalbar_halpha(sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
   case(2)
      get_opalbar = opalbar_hbeta(sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
   case(10)
      get_opalbar = opalbar_model_kline(yhe, hei, rho, kline)*sr
   case(11)
      get_opalbar = opalbar_model_kline2(yhe, hei, rho, kline)*sr
   case default
      stop 'error in get_opalbar: iline not properly specified'
end select

!write(*,*) sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
!
end function get_opalbar
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function get_opalbar2(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)
!
!get opalbar for a given line transition described by iline
!opacity law for balmer lines by assuming neutral hydrogen in ground state  
!
!input:  sr   stellar radius in cgs
!        yhe, hei:  helium abundance and number of free electrons per helium
!        temp:  temperature in K
!        vth_fiducial:  fiducial thermal velocity in cm/s
!        xnue0:   transition frequency
!        bl,bu:   NLTE-departure coefficients for lower and upper level
!        rho:    density in cgs
!        kline:   just a scaling factor to increase or decrease opacity (line strength)
!  
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: iline
real(dp), intent(in) :: sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
real(dp) :: get_opalbar2
!
! ... local functions
real(dp) :: opalbar_halpha2, opalbar_hbeta2, opalbar_model_kline
!
!
select case(iline)
   case(1)
      get_opalbar2 = opalbar_halpha2(sr, yhe, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
!   case(2)
!      get_opalbar2 = opalbar_hbeta2(sr, yhe, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
   case(10)
      get_opalbar2 = opalbar_model_kline(yhe, hei, rho, kline)*sr
   case default
      stop 'error in get_opalbar2: iline not properly specified'
end select

!write(*,*) sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
!
end function get_opalbar2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
function get_opalbar3(iline, kline, sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)
!
!get opalbar for a given line transition described by iline
!opacity law for balmer lines by assuming LTE ionization balance for hydrogen (given hei)
!
!input:  sr   stellar radius in cgs
!        yhe, hei:  helium abundance and number of free electrons per helium
!        temp:  temperature in K
!        vth_fiducial:  fiducial thermal velocity in cm/s
!        xnue0:   transition frequency
!        bl,bu:   NLTE-departure coefficients for lower and upper level
!        rho:    density in cgs
!        kline:   just a scaling factor to increase or decrease opacity (line strength)
!  
use prog_type
use fund_const
!
implicit none
!
! ... arguments
integer(i4b), intent(in) :: iline
real(dp), intent(in) :: sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
real(dp) :: get_opalbar3
!
! ... local functions
real(dp) :: opalbar_halpha3, opalbar_hbeta3, opalbar_model_kline
!
!
select case(iline)
   case(1)
      get_opalbar3 = opalbar_halpha3(sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
!   case(2)
!      get_opalbar3 = opalbar_hbeta3(sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho)*kline
   case(10)
      get_opalbar3 = opalbar_model_kline(yhe, hei, rho, kline)*sr
   case default
      stop 'error in get_opalbar3: iline not properly specified'
end select

!write(*,*) sr, yhe, hei, temp, vth_fiducial, xnue0, bl, bu, rho, kline
!
end function get_opalbar3
