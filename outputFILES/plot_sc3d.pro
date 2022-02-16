pro plot_sc3d, windx=windx, oname=oname
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
pdefault=!p
;
if(not keyword_set(windx)) then windx=0
;
;-----------------------------------------------------------------------
;
;for bstar
fname = 'sc3d/model3d.h5'
;
;-----------------------------------------------------------------------
;
read_sc3d, fname, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
                  xarr=x, yarr=y, zarr=z, opac3d=opac3d, scont3d=scont3d, t3d=t3d, $
                  nodes_nue=nodes_nue, nnue=nnue, epsmaxc_arr=epsmaxc_arr

xnue0=nodes_nue(0)
bnue3d=fltarr(ndxmax,ndymax,ndzmax)*0.d0
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         bnue3d(i,j,k) = bnue(xnue0,t3d(i,j,k))
      endfor
   endfor
endfor
scont3d=scont3d/bnue3d
;
;------------prepare the slice for which 2d shall be plotted------------
;
;xz-plane at y=yindx
scont2d=fltarr(ndxmax,ndzmax)
j=ndymax/2
for i=0, ndxmax-1 do begin
   for k=0, ndzmax-1 do begin
      scont2d(i,k)=scont3d(i,j,k)
   endfor
endfor
;
;contour_model2d, x, z, scont2d, zlim=[-.3d0,1.5d0], xlim=[-1.2,1.2], windx=windx
;
;-----------------------------------------------------------------------
;
plot_model3d, x, y, z, scont3d, opac3d, bnue3d, windx=windx
;
plot_convergence3d, epsmaxc_arr, windx=windx

!p=pdefault
end
;
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;-----------------------------------------------------------------------
;
pro contour_model2d, x, z, arr2d, xlim=xlim, zlim=zlim, clim=clim, ctitlestr=ctitlestr, windx=windx, oname=oname, isotropic=isotropic
;
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
;make contour plot
contourplots_single, arr2d, x, z, $
              ncolors, nlevels_iso, bottom, levels_final, levels_iso, $
              c_colors, $
              cb_indx, cb_tickmark_name, $ 
              titleStr=titleStr, xtitleStr=xtitleStr, ytitleStr=ytitleStr, $
              xlim=xlim, ylim=zlim, $
              ctable=13, ctitleStr=ctitleStr, /background2, isotropic=isotropic
;
;overplot the grid
nx=n_elements(x)
nz=n_elements(z)
loadct, 0
;for i=0, nx-1 do begin
;   oplot, [x(i),x(i)], $
;          [z(0),z(nz-1)]
;endfor
;for i=0, nz-1 do begin
;   oplot, [x(0),x(nx-1)], $
;          [z(i),z(i)]
;endfor
;
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
pro plot_model3d, x, y, z, scont3d, opac3d, bnue3d, xlim=xlim, windx=windx, oname=oname
;
;calculate 1d arrays along z
nx=n_elements(x)
ny=n_elements(y)
nz=n_elements(z)
;
ix=nx/2
iy=ny/2
;
opac1d=opac3d(ix,iy,*)
scont1d=scont3d(ix,iy,*)
bnue1d=bnue3d(ix,iy,*)
;
tau1d=fltarr(nz)*0.d0
tau1d(nz-1)=1.d-3
for i=nz-2, 0, -1 do begin
   dtau=0.5d0*(opac1d(i+1)+opac1d(i))*(z(i+1)-z(i))
   tau1d(i)=tau1d(i+1)+dtau
endfor
;
;------------------------------------------------------------------------
;
loadct, 0
if(not keyword_set(windx)) then windx=0
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
;---------------------------opacities-----------------------------------
;
;window, windx
;device, decomposed=0
;windx=windx+1
;;
;xmax=max(tau1d)+0.1*mean(tau1d)
;xmin=1.d-3
;ymin=min(opac1d)-0.1*mean(opac1d)
;ymax=max(opac1d)+0.1*mean(opac1d)
;plot, tau1d, opac1d, $
;      xrange=[xmax,xmin], /xlog, $
;      yrange=[ymin,ymax], $
;      xtitle=textoidl('\tau'), $
;      ytitle=textoidl('\chi')
;
;--------------------------planck function------------------------------
;
;window, windx
;device, decomposed=0
;windx=windx+1
;;
;xmax=max(tau1d)+0.1*mean(tau1d)
;xmin=1.d-3
;ymin=min(bnue1d)-0.1*mean(bnue1d)
;ymax=max(bnue1d)+0.1*mean(bnue1d)
;plot, tau1d, bnue1d, $
;      xrange=[xmax,xmin], /xlog, $
;      yrange=[ymin,ymax], $
;      xtitle=textoidl('\tau'), $
;      ytitle=textoidl('B_\nue')
;
;--------------------------source function------------------------------
;
window, windx
device, decomposed=0
windx=windx+1
;
xmax=max(tau1d)+0.1*mean(tau1d)
xmin=1.d-3
ymin=min(scont1d)-0.1*mean(scont1d)
ymax=max(scont1d)+0.1*mean(scont1d)
plot, tau1d, scont1d, $
      xrange=[xmax,xmin], /xlog, $
      yrange=[ymin,ymax], $
      xtitle=textoidl('\tau'), $
      ytitle=textoidl('S/B')

end
;
;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
;
pro plot_convergence3d, epsmaxc_arr, windx=windx
;
if(not keyword_set(windx)) then windx=0
;
;----------------------prepare convergence behaviour--------------------
;
itmaxc=n_elements(epsmaxc_arr)
;
epsmaxc_arr=abs(epsmaxc_arr)
iternrc_arr=indgen(itmaxc)
;
;----------------------convergence behaviour continuum------------------
;
xmin=0
xmax=max(where(epsmaxc_arr gt 0))
if(xmax eq -1) then begin
   ymin=1.d-5
endif else begin
   ymin=epsmaxc_arr(xmax)
endelse
ymax=max(epsmaxc_arr)
;
ylim=[ymin,ymax]
;
ytitleStr=textoidl('((S_c^{(k-1)} - S_c^{(k)}) /S_c^{(k)})_{max}')
xtitleStr='# iterations'
;
window, windx
device, decomposed=0
windx=windx+1

plot, iternrc_arr, epsmaxc_arr, $
      /ylog, $
      line=0, $
      yrange=ylim, $
      xrange=[xmin,xmax], $
      title=titlestr, $
      xtitle=xtitlestr, $
      ytitle=ytitlestr, $
      charsize=2.

end
