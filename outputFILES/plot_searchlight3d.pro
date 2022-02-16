pro plot_searchlight3d, windx=windx, oname=oname
;
;-------------------------OUTPUT IF HELP NEEDED-------------------------
;
if(not keyword_set(windx)) then windx=0
;
;-----------------------------------------------------------------------
;
;for bstar
fname = 'sc3d/searchlight3d.h5'
;
;-----------------------------------------------------------------------
;
read_searchlight3d, fname, ndxmax=ndxmax, ndymax=ndymax, ndzmax=ndzmax, $
                      xarr=x, yarr=y, zarr=z, int3d=int3d, inttheo3d=int3d_theo, $
                      mask3d=mask3d, opac3d=opac3d, nn_x=nn_x, nn_y=nn_y, nn_z=nn_z

norm=max(int3d_theo)
int3d_sc=int3d/norm
int3d_theo=int3d_theo/norm
;
;-----------------------------------------------------------------------
;
tau1d=fltarr(ndzmax)*0.d0
tau1dz=fltarr(ndzmax)*0.d0
tau1d(0)=1.d-6
tau1d(1)=1.d-6
tau1dz(0)=1.d-6
tau1dz(1)=1.d-6
for i=2, ndzmax-1 do begin
   dtau=0.5d0*(opac3d(ndxmax/2,ndymax/2,i-1)+opac3d(ndxmax/2,ndymax/2,i))*(z(i)-z(i-1))
   tau1d(i)=tau1d(i-1)+dtau/nn_z
   tau1dz(i)=tau1dz(i-1)+dtau
endfor
;
z3d=fltarr(ndxmax,ndymax,ndzmax)*0.d0
x3d=fltarr(ndxmax,ndymax,ndzmax)*0.d0
y3d=fltarr(ndxmax,ndymax,ndzmax)*0.d0
tau3d=fltarr(ndxmax,ndymax,ndzmax)*0.d0
tau3dz=fltarr(ndxmax,ndymax,ndzmax)*0.d0
for i=0, ndxmax-1 do begin
   for j=0, ndymax-1 do begin
      for k=0, ndzmax-1 do begin
         x3d(i,j,k) = x(i)
         y3d(i,j,k) = y(j)
         z3d(i,j,k)=z(k)
         tau3d(i,j,k)=tau1d(k)
         tau3dz(i,j,k)=tau1dz(k)
      endfor
   endfor
endfor
;
indx=where(int3d_sc gt 1.d-6)
xlim=[min(tau3d(indx)),max(tau3d(indx))]
ymin=min([int3d_sc(indx),int3d_theo(indx)])
ymax=max([int3d_sc(indx),int3d_theo(indx)])
ylim=[ymin,ymax]
;
loadct, 0
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
if(keyword_set(oname)) then oname='ps_files/plot_searchlight3d.ps'
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
;
plot, [xlim(0),xlim(0)], [ylim(0),ylim(0)], $
      xrange=xlim, /xs, $
      yrange=ylim, /ylog, $
      title=titlestr, $
      xtitle=textoidl('\tau_s'), $
      ytitle=textoidl('I(\tau_s)'), $
      charsize=1.2
oplot, tau1d, int3d_theo(ndxmax/2,ndymax/2,*), line=0
oplot, tau3d, int3d_sc, psym=2, color=ci_blue
;
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;-----------------------------------------------------------------------
;


window, windx
windx=windx+1
device, decomposed=0
;
loadct, 0
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
indx=where(int3d_sc gt 1.d-6)
xlim=[min(tau3d(indx)),max(tau3d(indx))]
ylim=[0.95,1.05]

plot, [xlim(0),xlim(0)], [ylim(0),ylim(0)], $
      xrange=xlim, /xs, $
      yrange=ylim, /ys, $
      xtitle=textoidl('\tau_s'), $
      ytitle=textoidl('I(\tau_s)/I_{theo}(\tau_s)')
oplot, tau3d(indx), int3d_sc(indx)/int3d_theo(indx), psym=2, color=ci_blue
;
;at left boundary
;indx=where(x3d eq x(2))
;oplot, tau3d(0,0,*), int3d_sc(0,0,*)/int3d_theo(0,0,*), psym=2, color=ci_red
;
;-----------------------------------------------------------------------
;
if(keyword_set(oname)) then oname='ps_files/contour_searchlight3d.ps'
ctitleStr=textoidl('I/I_c')
titleStr=textoidl('n_z=')+string(nn_z,format='(f5.2)')
;
;
iy=ndymax/2
int2d=fltarr(ndxmax,ndzmax)
for i=0, ndxmax-1 do begin
   for k=0, ndymax-1 do begin
      int2d(i,k)=int3d_theo(i,iy,k)
      int2d(i,k)=int3d_sc(i,iy,k)
   endfor
endfor
contour_model2d, x, z, int2d, zlim=[-0.4,11.], xlim=[-5.,5.], windx=windx, oname=oname, titlestr=titlestr, ctitlestr=ctitlestr
;
;
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
;
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
