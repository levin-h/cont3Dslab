pro plot_sc1d, oname=oname, windx=windx
;
if(not keyword_set(windx)) then windx=0
;
dir_01='diff1d'
dir_02='sc1d'
;
lstr_01='M I'
lstr_01_diff='M II'
lstr_02='M III'
;
readcol, dir_01+'/output.dat', indx_01, z1d_01, opac1d_01, tau1d_01, scont1d_01, scont1d_diff_01, bnue1d_01
readcol, dir_02+'/output.dat', indx_02, z1d_02, opac1d_02, tau1d_02, scont1d_02, bnue1d_02
;
scont1d_01=scont1d_01/bnue1d_01
scont1d_diff_01=scont1d_diff_01/bnue1d_01
scont1d_02=scont1d_02/bnue1d_02
;
readcol, dir_01+'/convergence.dat', niter_01, epsmaxc_arr_01
readcol, dir_02+'/convergence.dat', niter_02, epsmaxc_arr_02
;
;---------------------------opacities-----------------------------------
;
window, windx
device, decomposed=0
windx=windx+1
;
xmax=max(tau1d_01)+0.1*mean(tau1d_01)
xmin=min(tau1d_01)
ymin=min(opac1d_01)-0.1*mean(opac1d_01)
ymax=max(opac1d_01)+0.1*mean(opac1d_01)
plot, tau1d_01, opac1d_01, $
      xrange=[xmax,xmin], /xlog, $
      yrange=[ymin,ymax], $
      xtitle=textoidl('\tau'), $
      ytitle=textoidl('\chi')
oplot, tau1d_02, opac1d_02, psym=1
;
legend, [lstr_01, lstr_02], $
   linestyle=[0,0], $
   psym=[0,1], $
   /right_legend, $
   charsize=2.    
;
;--------------------------planck function------------------------------
;
window, windx
device, decomposed=0
windx=windx+1
;
xmax=max(tau1d_01)+0.1*mean(tau1d_01)
xmin=min(tau1d_01)
ymin=min(bnue1d_01)-0.1*mean(bnue1d_01)
ymax=max(bnue1d_01)+0.1*mean(bnue1d_01)
plot, tau1d_01, bnue1d_01, $
      xrange=[xmax,xmin], /xlog, $
      yrange=[ymin,ymax], $
      xtitle=textoidl('\tau'), $
      ytitle=textoidl('B_\nue')
oplot, tau1d_02, bnue1d_02, psym=1

legend, [lstr_01, lstr_02], $
   linestyle=[0,0], $
   psym=[0,1], $
   /right_legend, $
   charsize=2.    
;
;--------------------------source function------------------------------
;
window, windx
device, decomposed=0
windx=windx+1
;
xmax=max(tau1d_01)+0.1*mean(tau1d_01)
xmin=min(tau1d_01)
ymin=min(scont1d_01)-0.1*mean(scont1d_01)
ymax=max(scont1d_01)+0.1*mean(scont1d_01)
plot, tau1d_01, scont1d_01, $
      xrange=[xmax,xmin], /xlog, $
      yrange=[ymin,ymax], $
      xtitle=textoidl('\tau'), $
      ytitle=textoidl('S/B')
oplot, tau1d_02, scont1d_02, psym=1

oplot, tau1d_01, scont1d_diff_01, line=2, thick=2.



legend, [lstr_01, lstr_01_diff, lstr_02], $
   linestyle=[0,2,0], $
   psym=[0,0,1], $
   thick=[2.,1.,1.], $
   /right_legend, $
   charsize=2.    
;
;------------------------convergence behaviour--------------------------
;
window, windx
device, decomposed=0
windx=windx+1


xmax=max([niter_01,niter_02])
xmin=min([niter_01,niter_02])
ymin=min([epsmaxc_arr_01,epsmaxc_arr_02])
yax=max([epsmaxc_arr_01,epsmaxc_arr_02])

plot, niter_01, epsmaxc_arr_01, $
      xrange=[xmin,xmax], $
      yrange=[ymin,ymax], $
      xtitle='# iterations', $
      ytitle='rel correction', $
      /ylog

print, niter_02
oplot, niter_02, epsmaxc_arr_02, $
       line=1

legend, [lstr_01, lstr_02], $
   linestyle=[0,1], $
   psym=[0,0], $
   /right_legend, $
   charsize=2.    
;
;------------------------different models-------------------------------
;
stop
dir='models/tau_max2d1/del_tau1d0/fvm1st'
read_scont, dir, tau1, siter1, sdiff1
dir='models/tau_max2d1/del_tau5d-1/fvm1st'
read_scont, dir, tau2, siter2, sdiff2
dir='models/tau_max2d1/del_tau1d-1/fvm1st'
read_scont, dir, tau3, siter3, sdiff3
dir='models/tau_max2d1/del_tau1d-2/fvm1st'
read_scont, dir, tau4, siter4, sdiff4
;
if(keyword_set(oname)) then begin
   print, 'writing output to: ', oname
   set_plot,'ps'
   device, file=oname, decomposed=0, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
loadct, 0
;
titlestr=textoidl('\tau = 20')
thick_bu=!p.thick
!p.thick=2.
;
plot, tau3, sdiff3, $
      xrange=[max(tau1),max([min(tau),1.d-1])], $
      /xlog, $
      xtitle=textoidl('\tau'), $
      ytitle='S/B', $
      title=titlestr
      line=0
!p.thick=thick_bu
oplot, tau1, siter1, line=1
oplot, tau2, siter2, line=2
oplot, tau3, siter3, line=3
oplot, tau4, siter4, line=4
;
lstr0=textoidl('diffusion equation')
lstr1=textoidl('FVM 1st order, \Delta \tau = 1')
lstr2=textoidl('FVM 1st order, \Delta \tau = 0.5')
lstr3=textoidl('FVM 1st order, \Delta \tau = 0.1')
lstr4=textoidl('FVM 1st order, \Delta \tau = 0.01')
;
legend, [lstr0, lstr1, lstr2, lstr3, lstr4], $
   linestyle=[0,1,2,3,4], $
   psym=[0,0,0,0,0], $
   thick=[2.,!p.thick,!p.thick,!p.thick,!p.thick], $
        /right_legend
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif

;
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
;
;-----------------------------------------------------------------------
;
pro read_convergence, dir, niter, scont_maxerr
;
readcol, dir+'/convergence.dat', niter, scont_maxerr
;
end
;
;-----------------------------------------------------------------------
;
pro read_searchlight, dir, tau, int1st, int2nd, intabs
;
readcol, dir+'/searchlight.dat', tau, int1st, int2nd, intabs
;
end
