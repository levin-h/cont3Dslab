pro plot_searchlight2d, windx=windx, oname=oname
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(not keyword_set(windx)) then windx=0
;
;-----------------------------------------------------------------------
;
;for bstar
fname = 'sc2d/searchlight2d.h5'
;
;-----------------------------------------------------------------------
;
read_searchlight2d, fname, ndxmax=ndxmax, ndzmax=ndzmax, $
                      xarr=x, zarr=z, int2d=int2d, inttheo2d=int2d_theo, $
                      mask2d=mask2d, opac2d=opac2d, nn_z=nn_z, epsmaxi_arr=epsmaxi_arr
norm=max(int2d_theo)
int2d=int2d/norm
int2d_theo=int2d_theo/norm

if(keyword_set(oname)) then oname='ps_files/contour_searchlight2d.ps'
ctitleStr=textoidl('I/I_c')
titleStr=textoidl('n_z=')+string(nn_z,format='(f5.2)')
contour_model2d, x, z, int2d, zlim=[-0.4,3.], xlim=[-1.5,1.5], windx=windx, oname=oname, titlestr=titlestr, ctitlestr=ctitlestr
;
;-----------------------------------------------------------------------
;
;plot at boundaries
int1da=int2d(0,*)
int1db=int2d(1,*)
int1d_theo=int2d_theo(0,*)

;
tau1d=fltarr(ndzmax)*0.d0
tau1dz=fltarr(ndzmax)*0.d0
tau1d(0)=1.d-6
tau1d(1)=1.d-6
tau1dz(0)=1.d-6
tau1dz(1)=1.d-6
for i=2, ndzmax-1 do begin
   dtau=0.5d0*(opac2d(ndxmax/2,i-1)+opac2d(ndxmax/2,i))*(z(i)-z(i-1))
   tau1d(i)=tau1d(i-1)+dtau/nn_z
   tau1dz(i)=tau1dz(i-1)+dtau
   print, tau1d(i), tau1dz(i)
endfor
indx=where(int1da gt 1.d-6)
xlim=[min(tau1d(indx)),max(tau1d(indx))]
ymin=min([int1da(indx),int1db(indx)])
ymax=max([int1da(indx),int1db(indx)])
ylim=[ymin,ymax]
;
loadct, 0
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
if(keyword_set(oname)) then oname='ps_files/plot_searchlight2d.ps'
if(keyword_set(oname)) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
titlestr=textoidl('n_z=')+string(nn_z,format='(f5.2)')
plot, tau1d, int1d_theo, $
      xtitle=textoidl('\tau_s'), $
      ytitle=textoidl('I'), $
      title=titlestr, $
      xrange=xlim, $
      yrange=ylim, /ylog, $
      charsize=1.2
oplot, tau1d, int1da, color=ci_blue
oplot, tau1d, int1db, color=ci_red


if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then oname='ps_files/convergence_searchlight2d.ps'
if(keyword_set(oname)) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
nit=n_elements(epsmaxi_arr)
niter=indgen(nit)
;
titlestr=textoidl('n_z=')+string(nn_z,format='(f5.2)')
plot, niter, abs(epsmaxi_arr), $
      ytitle='max relative correction', $
      xtitle='# iterations', $
      title=titlestr, $
      charsize=1.2, $
      /ylog


if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif


end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro contour_model2d, x, z, arr2d, xlim=xlim, zlim=zlim, clim=clim, ctitlestr=ctitlestr, titlestr=titlestr, windx=windx, oname=oname, isotropic=isotropic
;
;
pdefault=!p
;
;-----------------------------------------------------------------------
;
if(not keyword_set(windx)) then windx=0
if(not keyword_set(xlim)) then xlim=[-1.d0, 1.d0]
if(not keyword_set(zlim)) then zlim=[0.d0,5.d0]
if(not keyword_set(clim)) then begin
   cmin=min(arr2d)
   cmax=max(arr2d)
   clim=[cmin,cmax]
endif
;
;set titles
xtitleStr=textoidl('x')
ytitleStr=textoidl('z')
;
;get colors
get_contour_colors2, clim, ncolors, nlevels_iso, bottom, $
                    levels_final, levels_iso, c_colors, $
                     cb_indx, cb_tickmark_name
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then begin
   print, "writing output to: ", oname
   set_plot,'ps'
   device,file=oname, color=1, bits_per_pixel=8
endif else begin
   window, windx
   device, decomposed=0
   windx=windx+1
endelse
;
print, titlestr
;make contour plot
contourplots_single, arr2d, x, z, $
              ncolors, nlevels_iso, bottom, levels_final, levels_iso, $
              c_colors, $
              cb_indx, cb_tickmark_name, $ 
              titleStr=titleStr, xtitleStr=xtitleStr, ytitleStr=ytitleStr, $
              xlim=xlim, ylim=zlim, $
              ctable=13, ctitleStr=ctitleStr, /background2, isotropic=isotropic, /ogrid
;
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif

end
