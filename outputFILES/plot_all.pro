pro plot_all, oname=oname, windx=windx
;
if(not keyword_set(windx)) then windx=0
;
fname_01='diff1d/output.dat'
fname_02='sc1d/output.dat'
fname_03='sc2d/model2d.h5'
;
lstr_01='DIFF'
lstr_02='SC1D'
lstr_03='SC2D'
;
readcol, fname_01, indx_01, z1d_01, opac1d_01, tau1d_01, fdum, scont1d_01, bnue1d_01
readcol, fname_02, indx_02, z1d_02, opac1d_02, tau1d_02, scont1d_02, bnue1d_02
readcol, 'sc1d/convergence.dat', niter_02, epsmaxc_arr_02
epsmaxc_arr_02=abs(epsmaxc_arr_02)
;
read_sc2d, fname_03, ndxmax=ndxmax, ndzmax=ndzmax, xarr=x, zarr=z1d_03, opac2d=opac2d, scont2d=scont2d, t2d=t2d, nodes_nue=nodes_nue, nnue=nnue, epsmaxc_arr=epsmaxc_arr_03, nconv=nconv
;
nconv=floor(nconv)
epsmaxc_arr_03=abs(epsmaxc_arr_03(0:nconv-1))
niter_03=1.d0+findgen(nconv-1)
;
xnue0=nodes_nue(0)
bnue2d=fltarr(ndxmax,ndzmax)*0.d0
for i=0, ndxmax-1 do begin
   for j=0, ndzmax-1 do begin
;      bnue2d(i,j) = 1.d0;bnue(xnue0,t2d(i,j))
      bnue2d(i,j) = bnue(xnue0,t2d(i,j))      
   endfor
endfor
;
ix=ndxmax/2
;ix=10
opac1d_03=opac2d(ix,*)
scont1d_03=scont2d(ix,*)
bnue1d_03=bnue2d(ix,*)
tau1d_03=fltarr(ndzmax)*0.d0
tau1d_03(ndzmax-1)=tau1d_01(n_elements(tau1d_01)-1)
for i=ndzmax-2, 0, -1 do begin
   dtau=0.5d0*(opac1d_03(i+1)+opac1d_03(i))*(z1d_03(i+1)-z1d_03(i))
   tau1d_03(i)=tau1d_03(i+1)+dtau
endfor
;
;
;
scont1d_01=scont1d_01/bnue1d_01
scont1d_02=scont1d_02/bnue1d_02
scont1d_03=scont1d_03/bnue1d_03
;
;window, 0
;plot, z1d_02, opac1d_01
;;
;;test log-log inteprolation with shifts
;z1d_02=z1d_02-5.d0
;nd=1000
;z1d_test=min(z1d_02)+findgen(nd)*(max(z1d_02)-min(z1d_02))/(nd-1)
;;
;opac1d_test=fltarr(nd)*0.d0
;for i=0, nd-1 do begin
;;find the index
;   zp=z1d_test(i)
;   find_indx, zp,z1d_02,n_elements(z1d_02),iim1,ii
;   ziim1=z1d_02(iim1)
;   zii=z1d_02(ii)
;   const=1.d0-min([ziim1,zii])
;;
;   zp_tilde=alog10(zp+const)
;   ziim1_tilde=alog10(ziim1+const)
;   zii_tilde=alog10(zii+const)
;
;   fiim1=alog10(opac1d_01(iim1))
;   fii=alog10(opac1d_01(ii))
;
;   acoeff=1.d0-(zp_tilde-ziim1_tilde)/(zii_tilde-ziim1_tilde)
;   bcoeff=(zp_tilde-ziim1_tilde)/(zii_tilde-ziim1_tilde)
;   fp = acoeff*fiim1 + bcoeff*fii
;   opac1d_test(i)=opac1d_01(iim1)^acoeff * opac1d_01(ii)^bcoeff
;   
;   print, z1d_02(iim1), zp, z1d_02(ii)
;   print, z1d_test(i), z1d_02(ii), 10.d0^fiim1, 10.d0^fp, 10.d0^fii
;endfor;
;
;z1d_test=z1d_test+5.d0
;oplot, z1d_test, opac1d_test, psym=1
;stop
;stop
;
;---------------------------opacities-----------------------------------
;
loadct_lev, ci_black, ci_blue, ci_red, ci_cyan, ci_green, ci_magenta, ci_yellow, ci_white
;
window, windx
device, decomposed=0
windx=windx+1
;
!p.multi=[2,1,2]
;
xmax=max(tau1d_01)+0.1*mean(tau1d_01)
xmin=1.d-6
ymin=min(opac1d_01)-0.1*mean(opac1d_01)
ymax=max(opac1d_01)+0.1*mean(opac1d_01)
plot, tau1d_01, opac1d_01, $
      xrange=[xmax,xmin], /xlog, $
      yrange=[ymin,ymax], $
      xtitle=textoidl('\tau'), $
      ytitle=textoidl('\chi')
oplot, tau1d_02, opac1d_02, color=ci_blue
oplot, tau1d_03, opac1d_03, color=ci_red
;
legend, [lstr_01, lstr_02, lstr_03], $
   linestyle=[0,0,0], $
   psym=[0,0,0], $
   color=[0,ci_blue,ci_red], $  
   /right_legend, $
   charsize=2.    
;
;--------------------------planck function------------------------------
;
;window, windx
;device, decomposed=0
;windx=windx+1
;
!p.multi=[1,1,2]
;
xmax=max(tau1d_01)+0.1*mean(tau1d_01)
xmin=1.d-6
ymin=min(bnue1d_01)-0.1*mean(bnue1d_01)
ymax=max(bnue1d_01)+0.1*mean(bnue1d_01)
plot, tau1d_01, bnue1d_01, $
      xrange=[xmax,xmin], /xlog, $
      yrange=[ymin,ymax], $
      xtitle=textoidl('\tau'), $
      ytitle=textoidl('B_\nue')
oplot, tau1d_02, bnue1d_02, color=ci_blue
oplot, tau1d_03, bnue1d_03, color=ci_red

legend, [lstr_01, lstr_02, lstr_03], $
   linestyle=[0,0,0], $
   psym=[0,0,0], $
   color=[0,ci_blue,ci_red], $  
   /right_legend, $
   charsize=2.    

!p.multi=0
;
;--------------------------source function------------------------------
;
;
if(keyword_set(oname)) then oname='ps_files/plot_pp_scont.ps'
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
xmax=max(tau1d_01)+0.1*mean(tau1d_01)
xmin=min(tau1d_01)
ymin=min(scont1d_01)-0.1*mean(scont1d_01)
ymax=max(scont1d_01)+0.1*mean(scont1d_01)
plot, tau1d_01, scont1d_01, $
      xrange=[xmax,xmin], /xlog, $
      yrange=[ymin,ymax], $
      xtitle=textoidl('\tau'), $
      ytitle=textoidl('S/B')
oplot, tau1d_02, scont1d_02, color=ci_blue
oplot, tau1d_03, scont1d_03, color=ci_red
oplot, tau1d_03, scont1d_03, psym=1, color=ci_red
;
legend, [lstr_01, lstr_02, lstr_03], $
   linestyle=[0,0,0], $
   psym=[0,0,0], $
   thick=[2.,1.,1.], $
   color=[0,ci_blue,ci_red], $  
   /right_legend, $
        charsize=2.

if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;plot, z1d_01, scont1d_01
;oplot, z1d_02, scont1d_02, color=ci_blue
;oplot, z1d_03, scont1d_03, color=ci_red
;
;------------------------convergence behaviour--------------------------
;
if(keyword_set(oname)) then oname='ps_files/plot_pp_convergence.ps'
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
xmax=max([niter_02,niter_03])
xmin=min([niter_02,niter_03])
ymin=min([epsmaxc_arr_02,epsmaxc_arr_03])
ymax=max([epsmaxc_arr_02,epsmaxc_arr_03])

plot, [0.,0.], [0.,0.], $
      xrange=[xmin,xmax], $
      yrange=[ymin,ymax], $
      xtitle='# iterations', $
      ytitle='rel correction', $
      /ylog

oplot, niter_02, epsmaxc_arr_02, $
       color=ci_blue
oplot, niter_03, epsmaxc_arr_03, $
       color=ci_red

legend, [lstr_02, lstr_03], $
   linestyle=[0,0], $
   psym=[0,0], $
   color=[ci_blue,ci_red], $  
   /right_legend, $
   charsize=2.    
if keyword_set(oname) then begin
   device, /close
   set_plot,'x'
endif
;
;------------------------different models-------------------------------
;
print, tau1d_02
xmax=max(tau1d_01)+0.1*mean(tau1d_01)
xmin=min(tau1d_01)

window, windx
device, decomposed=0
windx=windx+1
plot, tau1d_02, scont1d_03/scont1d_02, $
      xrange=[xmax,xmin], /xlog, $
      xtitle=textoidl('\tau'), $
      ytitle=textoidl('S(SC2D)/S(SC1D)')
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
