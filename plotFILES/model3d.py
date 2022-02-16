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
def main(dir='../inputFILES', oname3d='./ps_files/model3d', ps=[0.,0.,0.], pe=[0.,0.,10.], vnorm=[0.,1.,0.],
         clim_rho=[-20.,-12.], clim_temp=[10.e3,100.e3], clim_velz=[1.e2,4.e3], xlim=[-0.3,0.3], ylim=[0.9,6.], windx=0):
#
#plots contours of searchlight along a direction
#
#------------------read all information from hdf5-file------------------
#
   fname=dir+'/model3d.h5'

   nx, ny, nz, x, y, z, rho3d, velx3d, vely3d, velz3d, \
      tgas3d, trad3d, unit_length, unit_density, unit_velocity, unit_temperature = read_model_slab3d(fname)

#r in r_star
   radius=z[2]
#velocity in km/s
   velx3d=velx3d*unit_velocity/1.e5
   vely3d=vely3d*unit_velocity/1.e5
   velz3d=velz3d*unit_velocity/1.e5
#density and temperature in cgs
   rho3d = rho3d*unit_density
   tgas3d = tgas3d*unit_temperature
   trad3d = trad3d*unit_temperature   
#
#***********************************************************************
#
#                  plot along a line from ps to pe
#
#***********************************************************************
#
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=.4
   ysize=xsize/aspect_ratio

   fig0 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1   
   plt.ion()
   plt.show()
   ax0 = fig0.subplots(3,1)
#
#length of line and direction of line
   length=np.sqrt((pe[0]-ps[0])**2 + (pe[1]-ps[1])**2 + (pe[2]-ps[2])**2)
   nn_x = (pe[0]-ps[0])/length
   nn_y = (pe[1]-ps[1])/length
   nn_z = (pe[2]-ps[2])/length
   
   nray = 101
   zray = np.linspace(0.,length,nray)
   rho1d=np.zeros(nray)
   velx1d=np.zeros(nray)
   vely1d=np.zeros(nray)
   velz1d=np.zeros(nray)
   tgas1d=np.zeros(nray)
   trad1d=np.zeros(nray)   
#
   for i in range(0,nray):
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

#
#-----------------------------------------------------------------------
#
   titlestr=r'along $(x_s,y_s,z_s),(x_e,y_e,z_e) =$ ({xs:.2f},{ys:.2f},{zs:.2f}),({xe:.2f},{ye:.2f},{ze:.2f})'.format(xs=ps[0],ys=ps[1],zs=ps[2],xe=pe[0],ye=pe[1],ze=pe[2])
#
   ax0[0].set_xlabel(r'$s [unit\_length]$')
   ax0[0].set_ylabel(r'$\log(\rho)$')
   ax0[0].set_title(titlestr)
   ax0[0].plot(zray,np.log10(rho1d))

   ax0[1].set_xlabel(r'$s [unit\_length]$')
   ax0[1].set_ylabel(r'$\log(T)$')
   ax0[1].plot(zray,np.log10(tgas1d),label=r'$T_{gas}')
   ax0[1].plot(zray,np.log10(trad1d),label=r'$T_{rad}')
   ax0[1].legend(ncol=2,fontsize='small',loc='upper right')   

   ax0[2].set_xlabel(r'$s [unit\_length]$')
   ax0[2].set_ylabel(r'$v_i [km/s]$')
   ax0[2].plot(zray,velx1d,color='black', label=r'$v_x$')
   ax0[2].plot(zray,vely1d,color='red', label=r'$v_y$')
   ax0[2].plot(zray,velz1d,color='blue', label=r'$v_z$')            
   ax0[2].legend(fontsize=12)

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
   xsize=xsize/2.54 #in inches
   aspect_ratio=1.2
   ysize=xsize/aspect_ratio

   fig1 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   plt.ion()
   plt.show()

   ax1 = fig1.subplots(4,1)
#
#-----------------------------------------------------------------------
#
#   ax1[0].set_xlim(xlim)
#   ax1[0].set_ylim(ylim)
   ax1[0].set_xlabel(r'$z$')
   ax1[0].set_ylabel(r'$x$')

   ax1[1].set_xlabel(r'$z$')
   ax1[1].set_ylabel(r'$x$')
   
   ax1[2].set_xlabel(r'$z$')
   ax1[2].set_ylabel(r'$x$')
   
   ax1[3].set_xlabel(r'$z$')
   ax1[3].set_ylabel(r'$x$')            

#
   rho2d=np.zeros(shape=(nz,nx))
   velz2d=np.zeros(shape=(nz,nx))
   tgas2d=np.zeros(shape=(nz,nx))
   trad2d=np.zeros(shape=(nz,nx))   
#
   x_coord=np.zeros(shape=(nz,nx))
   y_coord=np.zeros(shape=(nz,nx))

   j=np.int((ny-1)/2)
   for i in range(0,nx):
      for k in range(0,nz):
         x_coord[k][i] = x[i]
         y_coord[k][i] = z[k]
         rho2d[k][i] = rho3d[k][j][i]
         tgas2d[k][i] = tgas3d[k][j][i]
         trad2d[k][i] = trad3d[k][j][i]         
         velz2d[k][i] = velz3d[k][j][i]
#
#-------------------density---------------------------------------------
#
   indx = np.where(rho2d > 0.)
   rho_min = np.min(rho2d[indx])
   lrho2d = rho2d*1.
   indx = np.where(lrho2d <= 0.)
   lrho2d[indx] = rho_min
   lrho2d = np.log10(lrho2d)
   clevels, ticks = get_clevels(clim=clim_rho)

   contourplot = ax1[0].contourf(y_coord, x_coord, lrho2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[0],orientation='vertical',aspect=20, ticks=ticks)
   cbar.set_label(r'$\log(\rho)$')                 
#
#-------------------gas temperature-------------------------------------
#
#  indx = np.where(rho2d > 0.)
   clevels, ticks = get_clevels(clim=clim_temp)

   contourplot= ax1[1].contourf(y_coord, x_coord, tgas2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[1],orientation='vertical',aspect=20, ticks=ticks)
   cbar.set_label(r'$T_{gas} [K]$')
#
#-------------------radiation temperature-------------------------------
#
#  indx = np.where(rho2d > 0.)
   clevels, ticks = get_clevels(clim=clim_temp)

   contourplot= ax1[2].contourf(y_coord, x_coord, trad2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[2],orientation='vertical',aspect=20, ticks=ticks)
   cbar.set_label(r'$T_{rad} [K]$')                      
#
#------------------------z velocity-------------------------------------
#
   clevels, ticks = get_clevels(clim=clim_velz)   

   contourplot= ax1[3].contourf(y_coord, x_coord, velz2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig1.colorbar(contourplot,ax=ax1[3],orientation='vertical',aspect=20, ticks=ticks)
   cbar.set_label(r'$v_r [km/s]$')   

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
   
dir='../inputFILES'
windx = 1
main(dir=dir,windx=windx,clim_rho=[-13.3,-7.1], clim_temp=[80000.,200000.],clim_velz=[0.,4500.],  ps=[0.,0.,1.], pe=[0.,0.,6.], )


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
