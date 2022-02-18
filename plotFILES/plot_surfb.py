import sys
sys.path.append('./levpy/')
sys.path.append('./levpy/version_sc3d')

import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_vbin import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_misc import *
from read_models import *
from lev_contour import *
import imageio
import os
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#from colorspacious import cspace_converter
#
def main(fname='../outputFILES/surfb_model00.h5', oname1='./ps_files/fluxem', oname2='./ps_files/surfb',
         oname3='./ps_files/surfbc',
         windx=1, create_animation=False):
#+
# NAME:
#	plot_surfb_vbin
#
# PURPOSE:
#	This procedure plots surface brightnes of spec_vbin.eo output
#
#-
#

   cdict = {'red':   [[0.0,  0.0, 0.0],
                      [0.4,  0.0, 0.0],
                      [1.0,  0.9, 0.9]],
            'green': [[0.0,  0.0, 0.0],
                      [0.4,  0.0, 0.0],
                      [1.0,  0.9, 0.9]],
            'blue':  [[0.0,  0.0, 0.0],
                      [0.5,  0.7, 0.7],
                      [1.0,  1.0, 1.0]]}

   cmap_lev = LinearSegmentedColormap('cmap_lev', segmentdata=cdict, N=256)   

   
#------------------read all information from hdf5-file------------------
#
   xic1, vth_fiducial, xnue0, nxobs_fs, \
      xobs_photprof, xic_photprof, xicc_photprof, \
      xobs_profile, ftot_profile, fcont_profile, femi_profile, fabs_profile, \
      nx_surfb, ny_surfb, nxobs_surfb, x_surfb, y_surfb, xobs_surfb, \
      iem3d_surfb, iemi3d_surfb, iabs3d_surfb, icont3d_surfb = read_surfb_slab3d(fname)
#
#-------------plot photospheric profile and emergent profile-------------
#
   xsize=18. #in cm
   xsize=xsize/2.54 #in inches
   aspect_ratio=1./.8
   ysize=xsize/aspect_ratio
   fig1 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
#
   plt.ion()
   plt.show()
#
   ax = fig1.subplots(2,1)
#
#photospheric profile
   ax[0].set_xlabel(r'$x_{obs} [v_{th}^\ast]$')
   ax[0].set_ylabel(r'$F_{core}/F_{c,core}$')
   ax[0].set_xlim(np.min(xobs_photprof),np.max(xobs_photprof))
   ax[0].plot(xobs_photprof, xic_photprof/xicc_photprof)

#emergent profile
   ax[1].set_xlabel(r'$x_{obs} [v_{th}^\ast]$')
   ax[1].set_ylabel(r'$F_{L}/F_{c}$')
   ax[1].set_xlim(np.max(xobs_profile),np.min(xobs_profile))   
   ax[1].plot(xobs_profile, ftot_profile/fcont_profile, linestyle='solid', color='black')
   ax[1].plot(xobs_profile, femi_profile/fcont_profile, linestyle='dotted', color='black')
   ax[1].plot(xobs_profile, fabs_profile/fcont_profile, linestyle='dashed', color='black')         
#
#-----------------------plotting----------------------------------------
#
   xsize=18. #in cm
   xsize=xsize/1.54 #in inches
   aspect_ratio=1./.8
   ysize=xsize/aspect_ratio
   fig2 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
#
   plt.ion()
   plt.show()
#
   ncol=3
   nrow=3
   ax2 = fig2.subplots(nrow,ncol)


   for i in np.arange(0,ncol*nrow):
      if i < nxobs_surfb:
         icont2d_surfb = trad(xnue0, icont3d_surfb[:,:][i])#icont3d_surfb[:,:][i]
         iem2d_surfb = trad(xnue0, iem3d_surfb[:,:][i])/icont2d_surfb
         iemi2d_surfb = iemi3d_surfb[:,:][i]/icont2d_surfb
         iabs2d_surfb = iabs3d_surfb[:,:][i]/icont2d_surfb
        
         icol = i%ncol
         irow = np.int((i - icol)/ncol)
         clevels, ticks = get_clevels(clim=[0.,1.000001])
         titlestr=r'$v [km/s]=${xobs:.2f}'.format(xobs=xobs_surfb[i]*100)
         ax2[irow,icol].set_title(titlestr)
         ax2[irow,icol].set_xlabel(r'$x_p$')
         ax2[irow,icol].set_ylabel(r'$y_p$')
         ax2[irow,icol].axis('equal') #isotropic
         ax2[irow,icol].set_xlim(np.min(x_surfb),np.max(x_surfb))
         ax2[irow,icol].set_ylim(np.min(y_surfb),np.max(y_surfb))                     
         contourplot = ax2[irow,icol].contourf(x_surfb, y_surfb, iem2d_surfb,
                                                 levels=clevels,
                                                 extend='both',
                                                 cmap=cmap_lev)
         contourplot.cmap.set_over('white')
         contourplot.cmap.set_under('black')
         contourplot.changed()
         cbar = fig2.colorbar(contourplot,ax=ax2[irow,icol],orientation='vertical',aspect=40,ticks=ticks)
         cbar.set_label(r'$T_{rad} [K]$')         
#
#-------------------------continuum surface brightness----------------------------
#
   xsize=18. #in cm
   xsize=xsize/1.54 #in inches
   aspect_ratio=1./.4
   ysize=xsize/aspect_ratio
   fig3 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
#
   plt.ion()
   plt.show()
#
   ax3 = fig3.subplots(1,2)

   i=8

   icont2d_surfb = trad(xnue0,icont3d_surfb[:,:][i])
   iem2d_surfb = trad(xnue0,iem3d_surfb[:,:][i])
   iemi2d_surfb = iemi3d_surfb[:,:][i]/icont2d_surfb
   iabs2d_surfb = iabs3d_surfb[:,:][i]/icont2d_surfb

#   indx=np.where(iem2d_surfb > 1.0000001)
#   print(indx)
#   exit()

#   clevels, ticks = get_clevels(clim=[np.min(icont2d_surfb),np.max(icont2d_surfb)])
   clevels, ticks = get_clevels(clim=[60.,120.])
#   titlestr=r'$v [km/s]=${xobs:.2f}'.format(xobs=xobs_surfb[i]*100.)
#   ax3[0].set_title(titlestr)
   ax3[0].set_xlabel(r'$x_p$')
   ax3[0].set_ylabel(r'$y_p$')
   ax3[0].axis('equal') #isotropic
   ax3[0].set_xlim(np.min(x_surfb),np.max(x_surfb))
   ax3[0].set_ylim(np.min(y_surfb),np.max(y_surfb))                     
   contourplot = ax3[0].contourf(x_surfb, y_surfb, icont2d_surfb/1.e3,
                                 levels=clevels,
                                 extend='both',
                                 cmap=cmap_lev)
   contourplot = ax3[0].contourf(x_surfb, y_surfb, icont2d_surfb/1.e3,
                                 levels=clevels,
                                 extend='both',
                                 cmap=cmap_lev)   
   contourplot.cmap.set_over('white')
   contourplot.cmap.set_under('black')
   contourplot.changed()
   cbar = fig3.colorbar(contourplot,ax=ax3[0],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$T_{rad} [kK]$')         


#   clevels, ticks = get_clevels(clim=[0.,1.000001])
   clevels, ticks = get_clevels(clim=[60.,120.])
#   titlestr=r'$x[v_t^\ast]=${xobs:.2f}'.format(xobs=xobs_surfb[i])
   titlestr=r'$v [km/s]=${xobs:.2f}'.format(xobs=xobs_surfb[i]*100.)
   ax3[1].set_title(titlestr)
   ax3[1].set_xlabel(r'$x_p$')
   ax3[1].set_ylabel(r'$y_p$')
   ax3[1].axis('equal') #isotropic
   ax3[1].set_xlim(np.min(x_surfb),np.max(x_surfb))
   ax3[1].set_ylim(np.min(y_surfb),np.max(y_surfb))                     
   contourplot = ax3[1].contourf(x_surfb, y_surfb, iem2d_surfb/1.e3,
                                 levels=clevels,
                                 extend='both',
                                 cmap=cmap_lev)                     
#   contourplot = ax3[1].contourf(x_surfb, y_surfb, iem2d_surfb,
#                                 levels=clevels,
#                                 extend='both',
#                                 cmap=cmap_lev)
   contourplot.cmap.set_over('white')
   contourplot.cmap.set_under('black')
   contourplot.changed()
   cbar = fig3.colorbar(contourplot,ax=ax3[1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$T_{rad} [kK]$')         
#
#-------------------------------create an animation-------------------------------
#
   if create_animation:
      animation_files = []
      animation_images = []
      for i in np.arange(0,nxobs_surfb):
         xsize=18. #in cm
         xsize=xsize/1.54 #in inches
         aspect_ratio=1./.8
         ysize=xsize/aspect_ratio
         fig4 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
#         plt.ion()
#         plt.show()
         ax4 = fig4.subplots()
#
         icont2d_surfb = icont3d_surfb[:,:][i]
         iem2d_surfb = iem3d_surfb[:,:][i]/icont2d_surfb
         iemi2d_surfb = iemi3d_surfb[:,:][i]/icont2d_surfb
         iabs2d_surfb = iabs3d_surfb[:,:][i]/icont2d_surfb

        
         clevels, ticks = get_clevels(clim=[0.,1.000001])
         titlestr=r'$x[v_t^\ast]=${xobs:.2f}'.format(xobs=xobs_surfb[i])
         ax4.set_title(titlestr)
         ax4.set_xlabel(r'$x_p$')
         ax4.set_ylabel(r'$y_p$')
         ax4.axis('equal') #isotropic
         ax4.set_xlim(np.min(x_surfb),np.max(x_surfb))
         ax4.set_ylim(np.min(y_surfb),np.max(y_surfb))                     
         contourplot = ax4.contourf(x_surfb, y_surfb, iem2d_surfb,
                                    levels=clevels,
                                    extend='both',
                                    cmap='Blues')
         contourplot.cmap.set_over('black')
         contourplot.cmap.set_under('grey')
         contourplot.changed()
         cbar = fig4.colorbar(contourplot,ax=ax4,orientation='vertical',aspect=40,ticks=ticks)
         cbar.set_label(r'$T_{rad} [K]$')         

         onamea = './animation_files/surfb{indx:05d}.png'.format(indx=i)
         fig4.savefig(onamea, bbox_inches='tight')
         fig4.clf()         


      animation_files = []
      animation_images = []
      for i in np.arange(0,nxobs_surfb):
         onamea = './animation_files/surfb{indx:05d}.png'.format(indx=i)
         animation_files.append(onamea)
         animation_images.append(imageio.imread(onamea))
#
      imageio.mimwrite('./animation_files/surfb.gif', animation_images, fps=2, loop=1)

#
#---------------------------------------------------------------------------------
#
#
   onamea = oname1+'.png'
   onameb = oname1+'.ps'
   fig1.savefig(onamea, bbox_inches='tight')
   fig1.savefig(onameb, bbox_inches='tight')

   onamea = oname2+'.png'
   onameb = oname2+'.ps'
   fig2.savefig(onamea, bbox_inches='tight')
   fig2.savefig(onameb, bbox_inches='tight')   

   onamea = oname3+'.png'
   onameb = oname3+'.ps'
   fig3.savefig(onamea, bbox_inches='tight')
   fig3.savefig(onameb, bbox_inches='tight')   

#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################


fname='../outputFILES/surfb_model00.h5'
#fname='/lhome/levin/Postdoc/for_jon/plots_wr3d/slab3d/surfb_model00.h5'


windx = 1
main(fname=fname,windx=windx,create_animation=False)


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
