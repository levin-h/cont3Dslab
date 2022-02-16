import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_models import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_interpol3d import *
from lev_misc import *
from lev_contour import *
#
def main(dir='../outputFILES/sc3d', oname3d='./ps_files/sc3d', ps=[0.,0.,0.], pe=[0.,0.,10.], vnorm=[0.,1.,0.], xlim=[-0.25,0.25], ylim=[0.9,6.], windx=0):
#
#plots contours of searchlight along a direction
#
#------------------read all information from hdf5-file------------------
#
   fname=dir+'/output_sc3d.h5'
   nx, ny, nz, kcont, eps_cont, x, y, z, mask3d, maskb3d, \
      rho3d, opac3d, velx3d, vely3d, velz3d, tgas3d, trad3d, eps_cont3d = read_scslab3d(fname, read='model')
   nx, ny, nz, x, y, z, scont3d, mint3d, fcontx3d, fconty3d, fcontz3d, \
      kcontxx3d, kcontxy3d, kcontxz3d, kcontyy3d, kcontyz3d, kcontzz3d = read_scslab3d(fname, read='solution')

   unit_length, unit_density, unit_temperature, unit_velocity = read_scslab3d(fname, read='units')
   unit_length_cgs = unit_length*cgs_rsu
#velocity in km/s
   velx3d=velx3d/1.e5
   vely3d=vely3d/1.e5
   velz3d=velz3d/1.e5
#
#eddington factor
   fedd_xx3d = kcontxx3d/mint3d
   fedd_yy3d = kcontyy3d/mint3d
   fedd_zz3d = kcontzz3d/mint3d

#opacity in cgs
   opac3d = opac3d/unit_length_cgs
#
#***********************************************************************
#
#                  plot along a line from ps to pe
#
#***********************************************************************
#
   xsize=18. #in cm
   xsize=xsize/1.54 #in inches
   aspect_ratio=1.5
   ysize=xsize/aspect_ratio

   fig0 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1   
   plt.ion()
   plt.show()
   ax0 = fig0.subplots(4,3)
#
#length of line and direction of line
   length=np.sqrt((pe[0]-ps[0])**2 + (pe[1]-ps[1])**2 + (pe[2]-ps[2])**2)
   nn_x = (pe[0]-ps[0])/length
   nn_y = (pe[1]-ps[1])/length
   nn_z = (pe[2]-ps[2])/length
   
   nray = 101
   zray = np.linspace(0.,length,nray)
   fedd_xx1d=np.zeros(nray)
   fedd_yy1d=np.zeros(nray)
   fedd_zz1d=np.zeros(nray)
   fedd_tot1d=np.zeros(nray)   
   fcontx1d=np.zeros(nray)
   fconty1d=np.zeros(nray)
   fcontz1d=np.zeros(nray)
   kcontxx1d=np.zeros(nray)
   kcontxy1d=np.zeros(nray)
   kcontxz1d=np.zeros(nray)
   kcontyy1d=np.zeros(nray)
   kcontyz1d=np.zeros(nray)
   kcontzz1d=np.zeros(nray)   
   scont1d=np.zeros(nray)
   mint1d=np.zeros(nray)
   kappa1d=np.zeros(nray)
   opac1d=np.zeros(nray)   
   velx1d=np.zeros(nray)
   vely1d=np.zeros(nray)
   velz1d=np.zeros(nray)
   rho1d=np.zeros(nray)
   tgas1d=np.zeros(nray)
   trad1d=np.zeros(nray)   

#
   for i in np.arange(0,nray):
#position in coordinate system
      px = ps[0]+zray[i]*nn_x
      py = ps[1]+zray[i]*nn_y
      pz = ps[2]+zray[i]*nn_z
      iim1, ii = find_indx(px,x,nx)
      jjm1, jj = find_indx(py,y,ny)
      kkm1, kk = find_indx(pz,z,nz)
#interpolation
      rho1d[i] = interpol3d_8p_lin(rho3d[kkm1,jjm1,iim1], rho3d[kkm1,jjm1,ii], rho3d[kkm1,jj,iim1], rho3d[kkm1,jj,ii],
                                   rho3d[kk,jjm1,iim1], rho3d[kk,jjm1,ii], rho3d[kk,jj,iim1], rho3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      tgas1d[i] = interpol3d_8p_lin(tgas3d[kkm1,jjm1,iim1], tgas3d[kkm1,jjm1,ii], tgas3d[kkm1,jj,iim1], tgas3d[kkm1,jj,ii],
                                 tgas3d[kk,jjm1,iim1], tgas3d[kk,jjm1,ii], tgas3d[kk,jj,iim1], tgas3d[kk,jj,ii],
                                 x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      trad1d[i] = interpol3d_8p_lin(trad3d[kkm1,jjm1,iim1], trad3d[kkm1,jjm1,ii], trad3d[kkm1,jj,iim1], trad3d[kkm1,jj,ii],
                                 trad3d[kk,jjm1,iim1], trad3d[kk,jjm1,ii], trad3d[kk,jj,iim1], trad3d[kk,jj,ii],
                                 x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)      
      
      velx1d[i] = interpol3d_8p_lin(velx3d[kkm1,jjm1,iim1], velx3d[kkm1,jjm1,ii], velx3d[kkm1,jj,iim1], velx3d[kkm1,jj,ii],
                                    velx3d[kk,jjm1,iim1], velx3d[kk,jjm1,ii], velx3d[kk,jj,iim1], velx3d[kk,jj,ii],
                                    x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      vely1d[i] = interpol3d_8p_lin(vely3d[kkm1,jjm1,iim1], vely3d[kkm1,jjm1,ii], vely3d[kkm1,jj,iim1], vely3d[kkm1,jj,ii],
                                    vely3d[kk,jjm1,iim1], vely3d[kk,jjm1,ii], vely3d[kk,jj,iim1], vely3d[kk,jj,ii],
                                    x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      velz1d[i] = interpol3d_8p_lin(velz3d[kkm1,jjm1,iim1], velz3d[kkm1,jjm1,ii], velz3d[kkm1,jj,iim1], velz3d[kkm1,jj,ii],
                                    velz3d[kk,jjm1,iim1], velz3d[kk,jjm1,ii], velz3d[kk,jj,iim1], velz3d[kk,jj,ii],
                                    x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      fedd_xx1d[i] = interpol3d_8p_lin(fedd_xx3d[kkm1,jjm1,iim1], fedd_xx3d[kkm1,jjm1,ii], fedd_xx3d[kkm1,jj,iim1], fedd_xx3d[kkm1,jj,ii],
                                   fedd_xx3d[kk,jjm1,iim1], fedd_xx3d[kk,jjm1,ii], fedd_xx3d[kk,jj,iim1], fedd_xx3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      fedd_yy1d[i] = interpol3d_8p_lin(fedd_yy3d[kkm1,jjm1,iim1], fedd_yy3d[kkm1,jjm1,ii], fedd_yy3d[kkm1,jj,iim1], fedd_yy3d[kkm1,jj,ii],
                                   fedd_yy3d[kk,jjm1,iim1], fedd_yy3d[kk,jjm1,ii], fedd_yy3d[kk,jj,iim1], fedd_yy3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      fedd_zz1d[i] = interpol3d_8p_lin(fedd_zz3d[kkm1,jjm1,iim1], fedd_zz3d[kkm1,jjm1,ii], fedd_zz3d[kkm1,jj,iim1], fedd_zz3d[kkm1,jj,ii],
                                   fedd_zz3d[kk,jjm1,iim1], fedd_zz3d[kk,jjm1,ii], fedd_zz3d[kk,jj,iim1], fedd_zz3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      mint1d[i] = interpol3d_8p_lin(mint3d[kkm1,jjm1,iim1], mint3d[kkm1,jjm1,ii], mint3d[kkm1,jj,iim1], mint3d[kkm1,jj,ii],
                                   mint3d[kk,jjm1,iim1], mint3d[kk,jjm1,ii], mint3d[kk,jj,iim1], mint3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      scont1d[i] = interpol3d_8p_lin(scont3d[kkm1,jjm1,iim1], scont3d[kkm1,jjm1,ii], scont3d[kkm1,jj,iim1], scont3d[kkm1,jj,ii],
                                   scont3d[kk,jjm1,iim1], scont3d[kk,jjm1,ii], scont3d[kk,jj,iim1], scont3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      fcontx1d[i] = interpol3d_8p_lin(fcontx3d[kkm1,jjm1,iim1], fcontx3d[kkm1,jjm1,ii], fcontx3d[kkm1,jj,iim1], fcontx3d[kkm1,jj,ii],
                                   fcontx3d[kk,jjm1,iim1], fcontx3d[kk,jjm1,ii], fcontx3d[kk,jj,iim1], fcontx3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      fconty1d[i] = interpol3d_8p_lin(fconty3d[kkm1,jjm1,iim1], fconty3d[kkm1,jjm1,ii], fconty3d[kkm1,jj,iim1], fconty3d[kkm1,jj,ii],
                                   fconty3d[kk,jjm1,iim1], fconty3d[kk,jjm1,ii], fconty3d[kk,jj,iim1], fconty3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      fcontz1d[i] = interpol3d_8p_lin(fcontz3d[kkm1,jjm1,iim1], fcontz3d[kkm1,jjm1,ii], fcontz3d[kkm1,jj,iim1], fcontz3d[kkm1,jj,ii],
                                   fcontz3d[kk,jjm1,iim1], fcontz3d[kk,jjm1,ii], fcontz3d[kk,jj,iim1], fcontz3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      kcontxx1d[i] = interpol3d_8p_lin(kcontxx3d[kkm1,jjm1,iim1], kcontxx3d[kkm1,jjm1,ii], kcontxx3d[kkm1,jj,iim1], kcontxx3d[kkm1,jj,ii],
                                   kcontxx3d[kk,jjm1,iim1], kcontxx3d[kk,jjm1,ii], kcontxx3d[kk,jj,iim1], kcontxx3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      kcontxy1d[i] = interpol3d_8p_lin(kcontxy3d[kkm1,jjm1,iim1], kcontxy3d[kkm1,jjm1,ii], kcontxy3d[kkm1,jj,iim1], kcontxy3d[kkm1,jj,ii],
                                   kcontxy3d[kk,jjm1,iim1], kcontxy3d[kk,jjm1,ii], kcontxy3d[kk,jj,iim1], kcontxy3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      kcontxz1d[i] = interpol3d_8p_lin(kcontxz3d[kkm1,jjm1,iim1], kcontxz3d[kkm1,jjm1,ii], kcontxz3d[kkm1,jj,iim1], kcontxz3d[kkm1,jj,ii],
                                   kcontxz3d[kk,jjm1,iim1], kcontxz3d[kk,jjm1,ii], kcontxz3d[kk,jj,iim1], kcontxz3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      kcontyy1d[i] = interpol3d_8p_lin(kcontyy3d[kkm1,jjm1,iim1], kcontyy3d[kkm1,jjm1,ii], kcontyy3d[kkm1,jj,iim1], kcontyy3d[kkm1,jj,ii],
                                   kcontyy3d[kk,jjm1,iim1], kcontyy3d[kk,jjm1,ii], kcontyy3d[kk,jj,iim1], kcontyy3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      kcontyz1d[i] = interpol3d_8p_lin(kcontyz3d[kkm1,jjm1,iim1], kcontyz3d[kkm1,jjm1,ii], kcontyz3d[kkm1,jj,iim1], kcontyz3d[kkm1,jj,ii],
                                   kcontyz3d[kk,jjm1,iim1], kcontyz3d[kk,jjm1,ii], kcontyz3d[kk,jj,iim1], kcontyz3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      kcontzz1d[i] = interpol3d_8p_lin(kcontzz3d[kkm1,jjm1,iim1], kcontzz3d[kkm1,jjm1,ii], kcontzz3d[kkm1,jj,iim1], kcontzz3d[kkm1,jj,ii],
                                   kcontzz3d[kk,jjm1,iim1], kcontzz3d[kk,jjm1,ii], kcontzz3d[kk,jj,iim1], kcontzz3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      opac1d[i] = interpol3d_8p_lin(opac3d[kkm1,jjm1,iim1], opac3d[kkm1,jjm1,ii], opac3d[kkm1,jj,iim1], opac3d[kkm1,jj,ii],
                                   opac3d[kk,jjm1,iim1], opac3d[kk,jjm1,ii], opac3d[kk,jj,iim1], opac3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      kappa1d[i] = opac1d[i]/rho1d[i]
      fedd_tot1d[i] = fedd_xx1d[i] + fedd_yy1d[i] + fedd_zz1d[i]

   tauz1d=np.zeros(nray)
   for i in np.arange(nray-2,-1,-1):
      dtau = 0.5*(opac1d[i+1]+opac1d[i])*np.abs(zray[i+1]-zray[i])*unit_length_cgs
      tauz1d[i]=tauz1d[i+1]+dtau      
#
#-----------------------------------------------------------------------
#
   titlestr=r'along $(x_s,y_s,z_s),(x_e,y_e,z_e) =$ ({xs:.2f},{ys:.2f},{zs:.2f}),({xe:.2f},{ye:.2f},{ze:.2f})'.format(xs=ps[0],ys=ps[1],zs=ps[2],xe=pe[0],ye=pe[1],ze=pe[2])
#

   fig0.suptitle(titlestr)
   
   ax0[0,0].set_xlabel(r'$s [unit\_length]$')
   ax0[0,0].set_ylabel(r'$\log(\rho)$')
   ax0[0,0].plot(zray,np.log10(rho1d))

   ax0[0,1].set_xlabel(r'$s [unit\_length]$')
   ax0[0,1].set_ylabel(r'$\log(T)$')
   ax0[0,1].plot(zray,np.log10(tgas1d), label=r'$T_{gas}$')
   ax0[0,1].plot(zray,np.log10(trad1d), label=r'$T_{gas}$')
   ax0[0,1].legend(ncol=4,fontsize='small',loc='upper left')   

   ax0[0,2].set_xlabel(r'$s [unit\_length]$')
   ax0[0,2].set_ylabel(r'$v_i [km/s]$')
   ax0[0,2].plot(zray,velx1d,color='black', label=r'$v_x$')
   ax0[0,2].plot(zray,vely1d,color='red', label=r'$v_y$')
   ax0[0,2].plot(zray,velz1d,color='blue', label=r'$v_z$')            
   ax0[0,2].legend(ncol=4,fontsize='small',loc='upper left')

   ax0[1,0].set_xlabel(r'$s [unit\_length]$')
   ax0[1,0].set_ylabel(r'$\log(\chi)$')
   ax0[1,0].plot(zray,np.log10(opac1d))

   ax0[1,1].set_xlabel(r'$s [unit\_length]$')
   ax0[1,1].set_ylabel(r'$\kappa$')
   ax0[1,1].plot(zray,kappa1d)      

   ax0[1,2].set_xlabel(r'$s [unit\_length]$')
   ax0[1,2].set_ylabel(r'$\log (S_c)$')
   ax0[1,2].plot(zray,np.log10(scont1d))

   ax0[2,0].set_xlabel(r'$s [unit\_length]$')
   ax0[2,0].set_ylabel(r'$\tau_{s}$')
   ax0[2,0].plot(zray,tauz1d)
   ax0[2,0].plot(zray,1.+np.zeros(nray),linestyle='dashed', color='black')   
   
   ax0[2,1].set_xlabel(r'$s [unit\_length]$')
   ax0[2,1].set_ylabel(r'$\log (J)$')
   ax0[2,1].plot(zray,np.log10(mint1d))

   ax0[2,2].set_xlabel(r'$s [unit\_length]$')
   ax0[2,2].set_ylabel(r'$F_i$')
   ax0[2,2].plot(zray,fcontx1d,color='black', label=r'$F_x$')
   ax0[2,2].plot(zray,fconty1d,color='red', label=r'$F_y$')
   ax0[2,2].plot(zray,fcontz1d,color='blue', label=r'$F_z$')            
   ax0[2,2].legend(ncol=4,fontsize='small',loc='upper right')

   ax0[3,0].set_xlabel(r'$s [unit\_length]$')
   ax0[3,0].set_ylabel(r'$K_{ii}$')
   ax0[3,0].plot(zray,np.log10(kcontxx1d),color='black', label=r'$K_{xx}$')
   ax0[3,0].plot(zray,np.log10(kcontyy1d),color='red', label=r'$K_{yy}$')
   ax0[3,0].plot(zray,np.log10(kcontzz1d),color='blue', label=r'$K_{zz}$')            
   ax0[3,0].legend(ncol=4,fontsize='small',loc='lower right')


   ax0[3,1].set_xlabel(r'$s [unit\_length]$')
   ax0[3,1].set_ylabel(r'$f_{Edd}$')
   ax0[3,1].set_ylim(0.,1.1)
   ax0[3,1].plot(zray,fedd_xx1d,color='black', label=r'$f_{xx}$')
   ax0[3,1].plot(zray,fedd_yy1d,color='red', label=r'$f_{yy}$')
   ax0[3,1].plot(zray,fedd_zz1d,color='blue', label=r'$f_{zz}$')
   ax0[3,1].plot(zray,fedd_tot1d,color='green', label=r'$\sum f_{ii}$')             
   ax0[3,1].legend(ncol=4,fontsize='small',loc='upper left')


   ax0[3,2].set_xlabel(r'$s [unit\_length]$')
   ax0[3,2].set_ylabel(r'$K_{ij}$')
   ax0[3,2].plot(zray,kcontxy1d,color='black', label=r'$K_{xy}$')
   ax0[3,2].plot(zray,kcontxz1d,color='red', label=r'$K_{xz}$')
   ax0[3,2].plot(zray,kcontyz1d,color='blue', label=r'$K_{yz}$')            
   ax0[3,2].legend(ncol=4,fontsize='small',loc='upper right')

#
#***********************************************************************
#
#          contour plots for a plane with normal vector vnorm
#
#for the moment, plot only xz-plane
#
#***********************************************************************
#
   xsize=18. #in cm
   xsize=xsize/1.54 #in inches
   aspect_ratio=.4
   ysize=xsize/aspect_ratio

   fig1 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   plt.ion()
   plt.show()

   ax1 = fig1.subplots(8,1)
#
#-----------------------------------------------------------------------
#
#   ax1[0].set_xlim(xlim)
#   ax1[0].set_ylim(ylim)
#   ax1[0].set_xlabel(r'$z$')
#   ax1[0].set_ylabel(r'$x$')   

#
   rho2d=np.zeros(shape=(nz,nx))
   velz2d=np.zeros(shape=(nz,nx))
   tgas2d=np.zeros(shape=(nz,nx))
   trad2d=np.zeros(shape=(nz,nx))   
   opac2d=np.zeros(shape=(nz,nx))
   kappa2d=np.zeros(shape=(nz,nx))
   scont2d=np.zeros(shape=(nz,nx))
   mint2d=np.zeros(shape=(nz,nx))
   fedd_xx2d=np.zeros(shape=(nz,nx))
   fedd_yy2d=np.zeros(shape=(nz,nx))
   fedd_zz2d=np.zeros(shape=(nz,nx))

#
   x_coord=np.zeros(shape=(nz,nx))
   y_coord=np.zeros(shape=(nz,nx))

   j=np.int(ny/2)-1
   for i in np.arange(0,nx):
      for k in np.arange(0,nz):
         x_coord[k][i] = x[i]
         y_coord[k][i] = z[k]
         rho2d[k][i] = rho3d[k][j][i]
         tgas2d[k][i] = tgas3d[k][j][i]
         trad2d[k][i] = trad3d[k][j][i]         
         velz2d[k][i] = velz3d[k][j][i]
         opac2d[k][i] = opac3d[k][j][i]
         scont2d[k][i] = scont3d[k][j][i]
         mint2d[k][i] = mint3d[k][j][i]
         fedd_xx2d[k][i] = fedd_xx3d[k][j][i]
         fedd_yy2d[k][i] = fedd_yy3d[k][j][i]
         fedd_zz2d[k][i] = fedd_zz3d[k][j][i]
         kappa2d[k][i] = opac2d[k][i]/rho2d[k][i]

#calculate total optical depth along z and along x
   tauz2d=np.zeros(shape=(nz,nx))
   taux2d=np.zeros(shape=(nz,nx))
   for i in np.arange(0,nx):   
      for k in np.arange(nz-2,-1,-1):
         dtau = 0.5*(opac2d[k+1][i]+opac2d[k][i])*np.abs(z[k+1]-z[k])*unit_length_cgs
         tauz2d[k][i] = tauz2d[k+1][i]+dtau   
#
#-------------------density---------------------------------------------
#
   lrho2d = get_logvals(rho2d)
   
   clevels, ticks = get_clevels(clim=[np.min(lrho2d),np.max(lrho2d)])

   
   contourplot = ax1[0].contourf(y_coord, x_coord, lrho2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[0],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(\rho)$')                 
#
#-------------------gas temperature-------------------------------------
#
   lt2d = get_logvals(tgas2d)
   clevels, ticks=get_clevels(clim=[np.min(lt2d),np.max(lt2d)])
#
   contourplot= ax1[1].contourf(y_coord, x_coord, lt2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(T_{gas})$')               
#
#-------------------radiation temperature-------------------------------
#
   lt2d = get_logvals(trad2d)
#
   contourplot= ax1[2].contourf(y_coord, x_coord, lt2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[2],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(T_{rad})$')
#
#------------------------kappa------------------------------------------
#
   clevels, ticks = get_clevels(clim=[np.min(kappa2d),np.max(kappa2d)])

   contourplot= ax1[3].contourf(y_coord, x_coord, kappa2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[3],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\kappa$')   
#
#------------------------tauz------------------------------------------
#
   ltauz2d = get_logvals(tauz2d)
   clevels, ticks = get_clevels(clim=[np.min(ltauz2d),np.max(ltauz2d)],nticks=5)

   contourplot= ax1[4].contourf(y_coord, x_coord, ltauz2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[4],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(\tau_z)$')   
#
#---------------------scont---------------------------------------------
#
   lscont2d = get_logvals(scont2d)
   clevels, ticks = get_clevels(clim=[np.min(lscont2d),np.max(lscont2d)])

   contourplot= ax1[5].contourf(y_coord, x_coord, lscont2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[5],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(S_c)$')   
#
#---------------------mean intensity-----------------------------------
#
   lmint2d = get_logvals(mint2d)
   clevels, ticks = get_clevels(clim=[np.min(lmint2d),np.max(lmint2d)])

   contourplot= ax1[6].contourf(y_coord, x_coord, lmint2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[6],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(J)$')
#
#---------------------mean intensity-----------------------------------
#
   clevels, ticks = get_clevels(clim=[0.333,1.])

   contourplot= ax1[7].contourf(y_coord, x_coord, fedd_zz2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[7],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$f_{zz}$')      
#
#------------------------output to files-------------------------------
#
   oname1 = oname3d+'_radial.png'
   oname2 = oname3d+'_radial.ps'
   fig0.savefig(oname1, bbox_inches='tight')
   fig0.savefig(oname2, bbox_inches='tight')

   oname1 = oname3d+'.png'
   oname2 = oname3d+'.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')
#
#
#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################
   
dir='../outputFILES/sc3d'
dir='/lhome/levin/Postdoc/hydro_sims/levin_wr2d/outputFILES/sc3d'
windx = 1
main(dir=dir,windx=windx,  ps=[0.,0.,1.], pe=[0.,0.,6.], )


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
