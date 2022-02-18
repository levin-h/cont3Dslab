import sys
sys.path.append('./levpy/')
sys.path.append('./levpy/version_sc3d')

import matplotlib.pyplot as plt
import numpy as np
from const_cgs import *
from read_models import *
from lev_interpol1d import *
from lev_interpol2d import *
from lev_interpol3d import *
from lev_misc import *
from lev_contour import *
from amrvac_tools.datfiles.reading import amrvac_reader



def read_nico2d(file):
    unit_erad=250000000.0
    
    ds = amrvac_reader.load_file(file)
    ds.get_info()
#   print(ds.get_info())
    ad = ds.load_all_data()
    rho2d = ad['rho']
    velz2d = ad['v1']
    velx2d = ad['v2']
    vely2d = np.copy(velx2d)*0.
#required only to set the standard output
    trad2d = ad['Trad']
    tgas2d = ad['Tgas']
    egas2d = ad['e']
    opal2d = ad['OPAL']
    lambda2d = ad['lambda']
    edd2d = ad['Edd']
    gamma2d = ad['Gamma']
    kappa2d = ad['OPAL']
    erad2d = ad['r_e']
    fcontx2d = ad['F2']
    fcontz2d = ad['F1']        

    bounds_x, bounds_y = ds.get_bounds()

    z,x = ds.get_coordinate_arrays()

    return x, z, rho2d, velx2d, vely2d, velz2d, tgas2d, trad2d, lambda2d, edd2d, gamma2d, kappa2d, erad2d, fcontx2d, fcontz2d, unit_erad


def read_nico3d(file, yplane=0.):
    unit_erad=500000000.0
    ds = amrvac_reader.load_file(file)
    ds.get_info()
#   print(ds.get_info())
    ad = ds.load_all_data()
    rho3d = ad['rho']
    velz3d = ad['v1']
    vely3d = ad['v2']
    velx3d = ad['v3']    
#required only to set the standard output
    trad3d = ad['Trad']
    tgas3d = ad['Tgas']
    egas3d = ad['e']
    opal3d = ad['OPAL']
    lambda3d = ad['lambda']
    edd3d = ad['Edd']
    gamma3d = ad['Gamma']
    kappa3d = ad['OPAL']
    erad3d = ad['r_e']
    fcontx3d = ad['F3']
    fconty3d = ad['F2']    
    fcontz3d = ad['F1']        

    bounds_x, bounds_y, bounds_z = ds.get_bounds()

    z,y,x = ds.get_coordinate_arrays()

    nz=np.size(z)
    ny=np.size(y)
    nx=np.size(x)
    rho2d=np.zeros(shape=(nz,nx))
    velx2d=np.zeros(shape=(nz,nx))
    vely2d=np.zeros(shape=(nz,nx))
    velz2d=np.zeros(shape=(nz,nx))        
    tgas2d=np.zeros(shape=(nz,nx))
    trad2d=np.zeros(shape=(nz,nx))
    lambda2d=np.zeros(shape=(nz,nx))
    edd2d=np.zeros(shape=(nz,nx))
    gamma2d=np.zeros(shape=(nz,nx))    
    kappa2d=np.zeros(shape=(nz,nx))
    erad2d=np.zeros(shape=(nz,nx))
    fcontx2d=np.zeros(shape=(nz,nx))
    fcontz2d=np.zeros(shape=(nz,nx))

    #interpolate onto 2d array
    jjm1, jj = find_indx(yplane,y,ny)
    yjm1 = y[jjm1]
    yj = y[jj]
    j=int(ny/2)
    for i in np.arange(0,nx):
        for k in np.arange(0,nz):            
            rho2d[k][i] = interpol_yp(yjm1, yj, rho3d[k][jjm1][i], rho3d[k][jj][i], yplane)
            tgas2d[k][i] = interpol_yp(yjm1, yj, tgas3d[k][jjm1][i], tgas3d[k][jj][i], yplane)
            trad2d[k][i] = interpol_yp(yjm1, yj, trad3d[k][jjm1][i], trad3d[k][jj][i], yplane)
            velx2d[k][i] = interpol_yp(yjm1, yj, velx3d[k][jjm1][i], velx3d[k][jj][i], yplane)
            vely2d[k][i] = interpol_yp(yjm1, yj, vely3d[k][jjm1][i], vely3d[k][jj][i], yplane)
            lambda2d[k][i] = interpol_yp(yjm1, yj, lambda3d[k][jjm1][i], lambda3d[k][jj][i], yplane)
            edd2d[k][i] = interpol_yp(yjm1, yj, edd3d[k][jjm1][i], edd3d[k][jj][i], yplane)
            gamma2d[k][i] = interpol_yp(yjm1, yj, gamma3d[k][jjm1][i], gamma3d[k][jj][i], yplane)
            kappa2d[k][i] = interpol_yp(yjm1, yj, kappa3d[k][jjm1][i], kappa3d[k][jj][i], yplane)
            erad2d[k][i] = interpol_yp(yjm1, yj, erad3d[k][jjm1][i], erad3d[k][jj][i], yplane)
            fcontx2d[k][i] = interpol_yp(yjm1, yj, fcontx3d[k][jjm1][i], fcontx3d[k][jj][i], yplane)
            fcontz2d[k][i] = interpol_yp(yjm1, yj, fcontz3d[k][jjm1][i], fcontz3d[k][jj][i], yplane)
            
#            tgas2d[k][i] = tgas3d[k][j][i]
#            trad2d[k][i] = trad3d[k][j][i]
#            velx2d[k][i] = velx3d[k][j][i]
#            vely2d[k][i] = vely3d[k][j][i]
#            velz2d[k][i] = velz3d[k][j][i]
#            lambda2d[k][i] = lambda3d[k][j][i]
#            edd2d[k][i] = edd3d[k][j][i]
#            gamma2d[k][i] = gamma3d[k][j][i]
#            kappa2d[k][i] = kappa3d[k][j][i]
#            erad2d[k][i] = erad3d[k][j][i]
#            fcontx2d[k][i] = fcontx3d[k][j][i]
#            fcontz2d[k][i] = fcontz3d[k][j][i]

    return x, z, rho2d, velx2d, vely2d, velz2d, tgas2d, trad2d, lambda2d, edd2d, gamma2d, kappa2d, erad2d, fcontx2d, fcontz2d, unit_erad

########################################################################
#
def main(dir='../outputFILES/sc3d', oname3d='./ps_files/sc3d', yplane=0., ps=[0.,0.,0.], pe=[0.,0.,10.], vnorm=[0.,1.,0.], xlim=[-0.25,0.25], ylim=[0.9,6.], windx=0):
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

   print(nx,ny,nz)

   unit_length, unit_density, unit_temperature, unit_velocity = read_scslab3d(fname, read='units')
   unit_length_cgs = unit_length*cgs_rsu
#
#hydro variables
   mint3d = mint3d*4.*np.pi/cgs_clight
   fcontx3d = fcontx3d*4.*np.pi
   fconty3d = fconty3d*4.*np.pi
   fcontz3d = fcontz3d*4.*np.pi
   kcontxx3d = kcontxx3d*4.*np.pi/cgs_clight
   kcontxy3d = kcontxy3d*4.*np.pi/cgs_clight
   kcontxz3d = kcontxz3d*4.*np.pi/cgs_clight
   kcontyy3d = kcontyy3d*4.*np.pi/cgs_clight
   kcontyz3d = kcontyz3d*4.*np.pi/cgs_clight
   kcontzz3d = kcontzz3d*4.*np.pi/cgs_clight
#
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
#calculate div(E)
   dedx3d = np.zeros(shape=(nz,ny,nx))
   dedy3d = np.zeros(shape=(nz,ny,nx))
   dedz3d = np.zeros(shape=(nz,ny,nx))
   rcap3d = np.zeros(shape=(nz,ny,nx))
   lambda3d01 = np.zeros(shape=(nz,ny,nx))
   lambda3d02 = np.zeros(shape=(nz,ny,nx))         
   hx = (x[1]-x[0])*unit_length_cgs
   hy = (y[1]-y[0])*unit_length_cgs
   hz = (z[1]-z[0])*unit_length_cgs
   for i in np.arange(2,nx-2):
      dxi=x[i]-x[i-1]
      dxip1=x[i+1]-x[i]       
      for j in np.arange(2,ny-2):
         dyj=y[j]-y[j-1]
         dyjp1=y[j+1]-y[j]                 
         for k in np.arange(2,nz-2):
            dzk=z[k]-z[k-1]
            dzkp1=z[k+1]-z[k]                              
            dedx3d[k][j][i] = ((mint3d[k][j][i+1]-mint3d[k][j][i])*dxi+(mint3d[k][j][i]-mint3d[k][j][i-1])*dxip1)/2./dxi/dxip1/unit_length_cgs
            dedy3d[k][j][i] = ((mint3d[k][j+1][i]-mint3d[k][j][i])*dyj+(mint3d[k][j][i]-mint3d[k][j-1][i])*dyjp1)/2./dyj/dyjp1/unit_length_cgs
            dedz3d[k][j][i] = ((mint3d[k+1][j][i]-mint3d[k][j][i])*dzk+(mint3d[k][j][i]-mint3d[k-1][j][i])*dzkp1)/2./dzk/dzkp1/unit_length_cgs
#            dedx3d[k][j][i] = (mint3d[k][j][i-2]/12.-2.*mint3d[k][j][i-1]/3.+2.*mint3d[k][j][i+1]/3.-mint3d[k][j][i+2]/12.)/hx
#            dedy3d[k][j][i] = (mint3d[k][j-2][i]/12.-2.*mint3d[k][j-1][i]/3.+2.*mint3d[k][j+1][i]/3.-mint3d[k][j+2][i]/12.)/hy
#            dedz3d[k][j][i] = (mint3d[k-2][j][i]/12.-2.*mint3d[k-1][j][i]/3.+2.*mint3d[k+1][j][i]/3.-mint3d[k+2][j][i]/12.)/hz
            rcap3d[k][j][i] = np.sqrt(dedx3d[k][j][i]**2 + dedy3d[k][j][i]**2 + dedz3d[k][j][i]**2)/mint3d[k][j][i]/opac3d[k][j][i]
            lambda3d01[k][j][i] = (2.+rcap3d[k][j][i])/(6.+3.*rcap3d[k][j][i] + rcap3d[k][j][i]**2)
            lambda3d02[k][j][i] = -opac3d[k][j][i]*(fcontx3d[k][j][i]+fconty3d[k][j][i]+fcontz3d[k][j][i])/(dedx3d[k][j][i]+dedy3d[k][j][i]+dedz3d[k][j][i])/cgs_clight
#            if k>=nz-3:
#                print(rcap3d[k][j][i],opac3d[k][j][i],mint3d[k][j][i],np.sqrt(dedx3d[k][j][i]**2 + dedy3d[k][j][i]**2 + dedz3d[k][j][i]**2))
            
#
#-------------------read all information from nico----------------------
#
   l3d = False
   if dir == '../outputFILES/sc3d_nico_0001':
      fname=dir+'/WR_2D_alpha_0.66_0001.dat'
   elif dir == '../outputFILES/sc3d_nico_0014':
      fname=dir+'/WR_2D_alpha_0.66_0014.dat'
   elif dir == '../outputFILES/sc3d_nico_0013':
      fname='../models/nico_wr3d/WR_3D_alpha_LTE_0013.dat'
      l3d = True
   elif dir == '../outputFILES/sc3d_nico_0013_test':
      fname='../models/nico_wr3d/WR_3D_alpha_LTE_0013.dat'
      l3d = True      
   elif dir == '../outputFILES/sc3d_nico_0020':
      fname='../models/nico_wr3d/WR_3D_alpha_LTE_0020.dat'
      l3d = True
   else:
       print('nicos model not specified')
       exit()


   if l3d:
       x_nico, z_nico, rho2d_nico, velx2d_nico, vely2d_nico, velz2d_nico, \
           tgas2d_nico, trad2d_nico, lambda2d_nico, fedd2d_nico, \
           gamma2d_nico, kappa2d_nico, erad2d_nico, fcontx2d_nico, fcontz2d_nico, unit_erad = read_nico3d(fname, yplane=yplane)
   else:
       x_nico, z_nico, rho2d_nico, velx2d_nico, vely2d_nico, velz2d_nico, \
           tgas2d_nico, trad2d_nico, lambda2d_nico, fedd2d_nico, \
           gamma2d_nico, kappa2d_nico, erad2d_nico, fcontx2d_nico, fcontz2d_nico, unit_erad = read_nico2d(fname)

   nx_nico = np.size(x_nico)
   nz_nico = np.size(z_nico)   

   unit_opacity=0.2
   rho2d_nico=rho2d_nico*unit_density
   velx2d_nico=velx2d_nico*unit_velocity/1.e5
   vely2d_nico=vely2d_nico*unit_velocity/1.e5
   velz2d_nico=velz2d_nico*unit_velocity/1.e5
   kappa2d_nico=kappa2d_nico*unit_opacity
   unit_flux=10.0**(5.416)*cgs_lsu/4./np.pi/cgs_rsu**2
   erad2d_nico=erad2d_nico*unit_erad
   fcontx2d_nico=fcontx2d_nico*unit_flux
   fcontz2d_nico=fcontz2d_nico*unit_flux
   
   dedx2d_nico=np.zeros(shape=(nz_nico,nx_nico))
   dedz2d_nico=np.zeros(shape=(nz_nico,nx_nico))
   rcap2d_nico=np.zeros(shape=(nz_nico,nx_nico))
   lambda2d01_nico=np.zeros(shape=(nz_nico,nx_nico))
   lambda2d02_nico=np.zeros(shape=(nz_nico,nx_nico))
   hx = (x_nico[1]-x_nico[0])*unit_length_cgs
   hz = (z_nico[1]-z_nico[0])*unit_length_cgs
   for i in np.arange(2,nx_nico-2):
      for k in np.arange(2,nz_nico-2):
         dedx2d_nico[k][i] = (0.5*erad2d_nico[k][i+1]-0.5*erad2d_nico[k][i-1])/2./hx#(erad2d_nico[k][i-2]/12.-2.*erad2d_nico[k][i-1]/3.+2.*erad2d_nico[k][i+1]/3.-erad2d_nico[k][i+2]/12.)/hx
         dedz2d_nico[k][i] = (0.5*erad2d_nico[k+1][i]-0.5*erad2d_nico[k-1][i])/2./hz#(erad2d_nico[k-2][i]/12.-2.*erad2d_nico[k-1][i]/3.+2.*erad2d_nico[k+1][i]/3.-erad2d_nico[k+2][i]/12.)/hz
         rcap2d_nico[k][i] = np.sqrt(dedx2d_nico[k][i]**2 + dedz2d_nico[k][i]**2)/erad2d_nico[k][i]/rho2d_nico[k][i]/kappa2d_nico[k][i]
         lambda2d01_nico[k][i] = (2.+rcap2d_nico[k][i])/(6.+3.*rcap2d_nico[k][i] + rcap2d_nico[k][i]**2)
         lambda2d02_nico[k][i] = -kappa2d_nico[k][i]*rho2d_nico[k][i]*(fcontx2d_nico[k][i]+fcontz2d_nico[k][i])/(dedx2d_nico[k][i]+dedz2d_nico[k][i])/cgs_clight
#         print(lambda2d_nico[k][i]/lambda2d01_nico[k][i],lambda2d02_nico[k][i])
#   print(lambda2d01_nico/lambda2d_nico)
#   exit()
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

   plt.ion()
   plt.show()
#  
   fig0 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1   
   ax0 = fig0.subplots(3,1)

   fig1 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1   
   ax1 = fig1.subplots(3,1)
   
   fig2 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1   
   ax2 = fig2.subplots(3,1)      
#
#length of line and direction of line
   length=np.sqrt((pe[0]-ps[0])**2 + (pe[1]-ps[1])**2 + (pe[2]-ps[2])**2)
   nn_x = (pe[0]-ps[0])/length
   nn_y = (pe[1]-ps[1])/length
   nn_z = (pe[2]-ps[2])/length
   
   nray = 501
   zray = np.linspace(0.,length,nray)
   fedd_xx1d=np.zeros(nray)
   fedd_yy1d=np.zeros(nray)
   fedd_zz1d=np.zeros(nray)
   fedd_tot1d=np.zeros(nray)
   fedd1d=np.zeros(nray)
   fcontx1d=np.zeros(nray)
   fconty1d=np.zeros(nray)
   fcontz1d=np.zeros(nray)
   mint1d=np.zeros(nray)
   kappa1d=np.zeros(nray)
   velx1d=np.zeros(nray)
   vely1d=np.zeros(nray)
   velz1d=np.zeros(nray)
   rho1d=np.zeros(nray)
   opac1d=np.zeros(nray)   
   tgas1d=np.zeros(nray)
   trad1d=np.zeros(nray)   
   dedx1d=np.zeros(nray)
   dedy1d=np.zeros(nray)
   dedz1d=np.zeros(nray)
   rcap1d=np.zeros(nray)   
   lambda1d01=np.zeros(nray)
   lambda1d02=np.zeros(nray)      

   rho1d_nico = np.zeros(nray)
   kappa1d_nico = np.zeros(nray)
   tgas1d_nico = np.zeros(nray)
   trad1d_nico = np.zeros(nray)   
   fedd1d_nico = np.zeros(nray)
   erad1d_nico = np.zeros(nray)
   fcontx1d_nico = np.zeros(nray)
   fcontz1d_nico = np.zeros(nray)  
   dedx1d_nico=np.zeros(nray)
   dedz1d_nico=np.zeros(nray)
   rcap1d_nico=np.zeros(nray)      
   lambda1d_nico = np.zeros(nray)   
   lambda1d01_nico=np.zeros(nray)
   lambda1d02_nico=np.zeros(nray)                
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
      fcontx1d[i] = interpol3d_8p_lin(fcontx3d[kkm1,jjm1,iim1], fcontx3d[kkm1,jjm1,ii], fcontx3d[kkm1,jj,iim1], fcontx3d[kkm1,jj,ii],
                                   fcontx3d[kk,jjm1,iim1], fcontx3d[kk,jjm1,ii], fcontx3d[kk,jj,iim1], fcontx3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      fconty1d[i] = interpol3d_8p_lin(fconty3d[kkm1,jjm1,iim1], fconty3d[kkm1,jjm1,ii], fconty3d[kkm1,jj,iim1], fconty3d[kkm1,jj,ii],
                                   fconty3d[kk,jjm1,iim1], fconty3d[kk,jjm1,ii], fconty3d[kk,jj,iim1], fconty3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      fcontz1d[i] = interpol3d_8p_lin(fcontz3d[kkm1,jjm1,iim1], fcontz3d[kkm1,jjm1,ii], fcontz3d[kkm1,jj,iim1], fcontz3d[kkm1,jj,ii],
                                   fcontz3d[kk,jjm1,iim1], fcontz3d[kk,jjm1,ii], fcontz3d[kk,jj,iim1], fcontz3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      opac1d[i] = interpol3d_8p_lin(opac3d[kkm1,jjm1,iim1], opac3d[kkm1,jjm1,ii], opac3d[kkm1,jj,iim1], opac3d[kkm1,jj,ii],
                                   opac3d[kk,jjm1,iim1], opac3d[kk,jjm1,ii], opac3d[kk,jj,iim1], opac3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      dedx1d[i] = interpol3d_8p_lin(dedx3d[kkm1,jjm1,iim1], dedx3d[kkm1,jjm1,ii], dedx3d[kkm1,jj,iim1], dedx3d[kkm1,jj,ii],
                                 dedx3d[kk,jjm1,iim1], dedx3d[kk,jjm1,ii], dedx3d[kk,jj,iim1], dedx3d[kk,jj,ii],
                                 x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      dedy1d[i] = interpol3d_8p_lin(dedy3d[kkm1,jjm1,iim1], dedy3d[kkm1,jjm1,ii], dedy3d[kkm1,jj,iim1], dedy3d[kkm1,jj,ii],
                                 dedy3d[kk,jjm1,iim1], dedy3d[kk,jjm1,ii], dedy3d[kk,jj,iim1], dedy3d[kk,jj,ii],
                                 x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      dedz1d[i] = interpol3d_8p_lin(dedz3d[kkm1,jjm1,iim1], dedz3d[kkm1,jjm1,ii], dedz3d[kkm1,jj,iim1], dedz3d[kkm1,jj,ii],
                                 dedz3d[kk,jjm1,iim1], dedz3d[kk,jjm1,ii], dedz3d[kk,jj,iim1], dedz3d[kk,jj,ii],
                                 x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      lambda1d01[i] = interpol3d_8p_lin(lambda3d01[kkm1,jjm1,iim1], lambda3d01[kkm1,jjm1,ii], lambda3d01[kkm1,jj,iim1], lambda3d01[kkm1,jj,ii],
                                 lambda3d01[kk,jjm1,iim1], lambda3d01[kk,jjm1,ii], lambda3d01[kk,jj,iim1], lambda3d01[kk,jj,ii],
                                 x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      lambda1d02[i] = interpol3d_8p_lin(lambda3d02[kkm1,jjm1,iim1], lambda3d02[kkm1,jjm1,ii], lambda3d02[kkm1,jj,iim1], lambda3d02[kkm1,jj,ii],
                                 lambda3d02[kk,jjm1,iim1], lambda3d02[kk,jjm1,ii], lambda3d02[kk,jj,iim1], lambda3d02[kk,jj,ii],
                                 x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)
      rcap1d[i] = interpol3d_8p_lin(rcap3d[kkm1,jjm1,iim1], rcap3d[kkm1,jjm1,ii], rcap3d[kkm1,jj,iim1], rcap3d[kkm1,jj,ii],
                                   rcap3d[kk,jjm1,iim1], rcap3d[kk,jjm1,ii], rcap3d[kk,jj,iim1], rcap3d[kk,jj,ii],
                                   x[iim1], x[ii], y[jjm1], y[jj], z[kkm1], z[kk], px, py, pz)      

      kappa1d[i] = opac1d[i]/rho1d[i]
      fedd_tot1d[i] = fedd_xx1d[i] + fedd_yy1d[i] + fedd_zz1d[i]


      iim1, ii = find_indx(px,x_nico,nx_nico)
      kkm1, kk = find_indx(pz,z_nico,nz_nico)
      rho1d_nico[i] = interpol2d_4p_lin(rho2d_nico[kkm1][iim1], rho2d_nico[kkm1,ii], rho2d_nico[kk][iim1], rho2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      tgas1d_nico[i] = interpol2d_4p_lin(tgas2d_nico[kkm1][iim1], tgas2d_nico[kkm1,ii], tgas2d_nico[kk][iim1], tgas2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      trad1d_nico[i] = interpol2d_4p_lin(trad2d_nico[kkm1][iim1], trad2d_nico[kkm1,ii], trad2d_nico[kk][iim1], trad2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)      
      kappa1d_nico[i] = interpol2d_4p_lin(kappa2d_nico[kkm1][iim1], kappa2d_nico[kkm1,ii], kappa2d_nico[kk][iim1], kappa2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      erad1d_nico[i] = interpol2d_4p_lin(erad2d_nico[kkm1][iim1], erad2d_nico[kkm1,ii], erad2d_nico[kk][iim1], erad2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      fcontx1d_nico[i] = interpol2d_4p_lin(fcontx2d_nico[kkm1][iim1], fcontx2d_nico[kkm1,ii], fcontx2d_nico[kk][iim1], fcontx2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      fcontz1d_nico[i] = interpol2d_4p_lin(fcontz2d_nico[kkm1][iim1], fcontz2d_nico[kkm1,ii], fcontz2d_nico[kk][iim1], fcontz2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      lambda1d_nico[i] = interpol2d_4p_lin(lambda2d_nico[kkm1][iim1], lambda2d_nico[kkm1,ii], lambda2d_nico[kk][iim1], lambda2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      lambda1d01_nico[i] = interpol2d_4p_lin(lambda2d01_nico[kkm1][iim1], lambda2d01_nico[kkm1,ii], lambda2d01_nico[kk][iim1], lambda2d01_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      lambda1d02_nico[i] = interpol2d_4p_lin(lambda2d02_nico[kkm1][iim1], lambda2d02_nico[kkm1,ii], lambda2d02_nico[kk][iim1], lambda2d02_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)            
      fedd1d_nico[i] = interpol2d_4p_lin(fedd2d_nico[kkm1][iim1], fedd2d_nico[kkm1,ii], fedd2d_nico[kk][iim1], fedd2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      dedx1d_nico[i] = interpol2d_4p_lin(dedx2d_nico[kkm1][iim1], dedx2d_nico[kkm1,ii], dedx2d_nico[kk][iim1], dedx2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      dedz1d_nico[i] = interpol2d_4p_lin(dedz2d_nico[kkm1][iim1], dedz2d_nico[kkm1,ii], dedz2d_nico[kk][iim1], dedz2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)
      rcap1d_nico[i] = interpol2d_4p_lin(rcap2d_nico[kkm1][iim1], rcap2d_nico[kkm1,ii], rcap2d_nico[kk][iim1], rcap2d_nico[kk][ii],
                                        x_nico[iim1],x_nico[ii],z_nico[kkm1],z_nico[kk],px,pz)      

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
   
   ax0[0].set_xlabel(r'$s [unit\_length]$')
   ax0[0].set_ylabel(r'$\log(\rho)$')
   ax0[0].plot(zray,np.log10(rho1d),color='blue', label='Levin')
   ax0[0].plot(zray,np.log10(rho1d_nico), color='red', label='Nico')
   ax0[0].legend()

   ax0[1].set_xlabel(r'$s [unit\_length]$')
   ax0[1].set_ylabel(r'$\log(T)$')
   ax0[1].plot(zray,np.log10(tgas1d), color='blue', linestyle='solid', label=r'$T_{gas}$')
   ax0[1].plot(zray,np.log10(trad1d), color='blue', linestyle='dashed', label=r'$T_{rad}$')   
   ax0[1].plot(zray,np.log10(tgas1d_nico), linestyle='solid', color='red', label='Nico')
   ax0[1].plot(zray,np.log10(trad1d_nico), linestyle='dashed', color='red')

   ax0[2].set_xlabel(r'$s [unit\_length]$')
   ax0[2].set_ylabel(r'$\kappa$')
   ax0[2].plot(zray,kappa1d, color='blue')
   ax0[2].plot(zray,kappa1d_nico, color='red')

   ax1[0].set_xlabel(r'$s [unit\_length]$')
   ax1[0].set_ylabel(r'$\tau_{s}$')
   ax1[0].plot(zray,tauz1d)
   ax1[0].plot(zray,1.+np.zeros(nray),linestyle='dashed', color='black')   

   ax1[1].set_xlabel(r'$s [unit\_length]$')
   ax1[1].set_ylabel(r'$\log (E_{\rm rad})$')
   ax1[1].plot(zray,np.log10(mint1d), color='blue')
   ax1[1].plot(zray,np.log10(erad1d_nico), color='red')
   
   ax1[2].set_xlabel(r'$s [unit\_length]$')
   ax1[2].set_ylabel(r'$F_i$')
   ax1[2].plot(zray,fcontx1d,color='blue', label=r'$F_x$')
   ax1[2].plot(zray,fconty1d,color='blue', label=r'$F_y$', linestyle='dotted')
   ax1[2].plot(zray,fcontz1d,color='blue', label=r'$F_z$', linestyle='dashed')
   ax1[2].plot(zray,fcontx1d_nico,color='red', label=r'$F_x$')
   ax1[2].plot(zray,fcontz1d_nico,color='red', label=r'$F_z$', linestyle='dashed')               
   ax1[2].legend(ncol=5,fontsize='small',loc='upper right')

#   ax2[0].set_ylim(0.,3.)
   ax2[0].set_xlabel(r'$s [unit\_length]$')
   ax2[0].set_ylabel(r'$\lambda$')
   ax2[0].plot(zray,lambda1d01,color='blue',linestyle='solid')
   ax2[0].plot(zray,lambda1d02,color='blue',linestyle='dashed')      
   ax2[0].plot(zray,lambda1d01_nico,color='red',linestyle='solid')
   ax2[0].plot(zray,lambda1d02_nico,color='red',linestyle='dashed')
   ax2[0].plot(zray,lambda1d_nico,color='red',linestyle='dotted')   

   ax2[1].set_xlabel(r'$s [unit\_length]$')
   ax2[1].set_ylabel(r'$f_{Edd}$')
   ax2[1].plot(zray,fedd_xx1d,color='blue', linestyle='solid', label=r'$f_{xx}$')
   ax2[1].plot(zray,fedd_yy1d,color='blue', linestyle='dashed', label=r'$f_{yy}$')
   ax2[1].plot(zray,fedd_zz1d,color='blue', linestyle='dotted', label=r'$f_{zz}$')
   ax2[1].plot(zray,fedd_tot1d,color='green', linestyle='solid', label=r'$\sum f_{ii}$')   
   ax2[1].plot(zray,fedd1d_nico,color='red', linestyle='dashed')

   ax2[2].set_xlabel(r'$s [unit\_length]$')
   ax2[2].set_ylabel(r'$R$')
   ax2[2].plot(zray,rcap1d,color='blue', linestyle='solid')
   ax2[2].plot(zray,rcap1d_nico,color='red', linestyle='solid')

#   ax2[2].set_xlabel(r'$s [unit\_length]$')
#   ax2[2].set_ylabel(r'$gradE_{Edd}$')
#   ax2[2].plot(zray,dedx1d,color='blue', linestyle='solid', label=r'$\partial x$')
#   ax2[2].plot(zray,dedy1d,color='blue', linestyle='dashed', label=r'$\partial y$')
#   ax2[2].plot(zray,dedz1d,color='blue', linestyle='dotted', label=r'$\partial z$')
#   ax2[2].plot(zray,dedx1d_nico,color='red', linestyle='solid', label=r'$\partial x$')
#   ax2[2].plot(zray,dedz1d_nico,color='red', linestyle='dotted', label=r'$\partial z$')
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

   fig3 = plt.figure(windx,figsize=(xsize,ysize),constrained_layout=True)
   windx=windx+1
   plt.ion()
   plt.show()

   ax3 = fig3.subplots(6,2)
#
#-----------------------------------------------------------------------
#
   rho2d=np.zeros(shape=(nz,nx))
   tgas2d=np.zeros(shape=(nz,nx))
   trad2d=np.zeros(shape=(nz,nx))   
   opac2d=np.zeros(shape=(nz,nx))
   kappa2d=np.zeros(shape=(nz,nx))
   mint2d=np.zeros(shape=(nz,nx))
   lambda2d01=np.zeros(shape=(nz,nx))
   lambda2d02=np.zeros(shape=(nz,nx))
   rcap2d=np.zeros(shape=(nz,nx))      
#
   x_coord=np.zeros(shape=(nz,nx))
   y_coord=np.zeros(shape=(nz,nx))
   


   #interpolate onto 2d array
   jjm1, jj = find_indx(yplane,y,ny)
   yjm1 = y[jjm1]
   yj = y[jj]
   j=int(ny/2)
   for i in np.arange(0,nx):
       for k in np.arange(0,nz):
           x_coord[k][i] = x[i]
           y_coord[k][i] = z[k]
           
           rho2d[k][i] = interpol_yp(yjm1, yj, rho3d[k][jjm1][i], rho3d[k][jj][i], yplane)
           tgas2d[k][i] = interpol_yp(yjm1, yj, tgas3d[k][jjm1][i], tgas3d[k][jj][i], yplane)
           trad2d[k][i] = interpol_yp(yjm1, yj, trad3d[k][jjm1][i], trad3d[k][jj][i], yplane)
           opac2d[k][i] = interpol_yp(yjm1, yj, opac3d[k][jjm1][i], opac3d[k][jj][i], yplane)
           mint2d[k][i] = interpol_yp(yjm1, yj, mint3d[k][jjm1][i], mint3d[k][jj][i], yplane)
           lambda2d01[k][i] = interpol_yp(yjm1, yj, lambda3d01[k][jjm1][i], lambda3d01[k][jj][i], yplane)
           lambda2d02[k][i] = interpol_yp(yjm1, yj, lambda3d02[k][jjm1][i], lambda3d02[k][jj][i], yplane)           
           rcap2d[k][i] = interpol_yp(yjm1, yj, rcap3d[k][jjm1][i], rcap3d[k][jj][i], yplane)
           kappa2d[k][i] = opac2d[k][i]/rho2d[k][i]           

           #rho2d[k][i] = rho3d[k][j][i]
           #tgas2d[k][i] = tgas3d[k][j][i]
           #trad2d[k][i] = trad3d[k][j][i]
           #opac2d[k][i] = opac3d[k][j][i]
           #mint2d[k][i] = mint3d[k][j][i]
           #kappa2d[k][i] = opac2d[k][i]/rho2d[k][i]
           #lambda2d01[k][i] = lambda3d01[k][j][i]
           #lambda2d02[k][i] = lambda3d02[k][j][i]
           #rcap2d[k][i] = rcap3d[k][j][i]                  
#
#-------------------density---------------------------------------------
#
   lrho2d = get_logvals(rho2d)
   lrho2d_nico = get_logvals(rho2d_nico)   
   
   clevels, ticks = get_clevels(clim=[np.min(lrho2d),np.max(lrho2d)])

   
   contourplot = ax3[0,0].contourf(y_coord, x_coord, lrho2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')

   contourplot = ax3[0,1].contourf(z_nico, x_nico, np.transpose(lrho2d_nico),
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig3.colorbar(contourplot,ax=ax3[0,1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(\rho)$')                 
#
#-------------------gas temperature-------------------------------------
#
   lt2d = get_logvals(tgas2d)
   ltgas2d_nico = get_logvals(tgas2d_nico)    
   clevels, ticks=get_clevels(clim=[np.min(lt2d),np.max(lt2d)])
#
   contourplot= ax3[1,0].contourf(y_coord, x_coord, lt2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')

   contourplot = ax3[1,1].contourf(z_nico, x_nico, np.transpose(ltgas2d_nico),
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig3.colorbar(contourplot,ax=ax3[1,1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(T_{gas})$')
#   
#-------------------radiation temperature-------------------------------
#
   lt2d = get_logvals(trad2d)
   ltrad2d_nico = get_logvals(trad2d_nico)    
   clevels, ticks=get_clevels(clim=[np.min(lt2d),np.max(lt2d)])
#
   contourplot= ax3[2,0].contourf(y_coord, x_coord, lt2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')

   contourplot = ax3[2,1].contourf(z_nico, x_nico, np.transpose(ltrad2d_nico),
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   
   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()

   cbar = fig3.colorbar(contourplot,ax=ax3[2,1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\log(T_{rad})$')   
#
#------------------------kappa------------------------------------------
#
   clevels, ticks = get_clevels(clim=[np.min(kappa2d),np.max(kappa2d)])

   contourplot= ax3[3,0].contourf(y_coord, x_coord, kappa2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot= ax3[3,1].contourf(z_nico, x_nico, np.transpose(kappa2d_nico),
                  levels=clevels,
                  extend='both',
                  cmap='jet')   

   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()
   
   cbar = fig3.colorbar(contourplot,ax=ax3[3,1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\kappa$')   
#
#------------------------lambda----------------------------------------
#
   clevels, ticks = get_clevels(clim=[np.min(lambda2d01),np.max(lambda2d01)])

   contourplot= ax3[4,0].contourf(y_coord, x_coord, lambda2d01,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot= ax3[4,1].contourf(z_nico, x_nico, np.transpose(lambda2d01_nico),
                  levels=clevels,
                  extend='both',
                  cmap='jet')   

   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()
   
   cbar = fig3.colorbar(contourplot,ax=ax3[4,1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$\lambda$')

#
#------------------------R parameter-----------------------------------
#
#   clevels, ticks = get_clevels(clim=[np.min(rcap2d_nico),np.max(rcap2d_nico)])
   clevels, ticks = get_clevels(clim=[0.,100.])

   contourplot= ax3[5,0].contourf(y_coord, x_coord, rcap2d,
                  levels=clevels,
                  extend='both',
                  cmap='jet')
   contourplot= ax3[5,1].contourf(z_nico, x_nico, np.transpose(rcap2d_nico),
                  levels=clevels,
                  extend='both',
                  cmap='jet')   

   contourplot.cmap.set_over('black')
   contourplot.cmap.set_under('grey')
   contourplot.changed()
   
   cbar = fig3.colorbar(contourplot,ax=ax3[5,1],orientation='vertical',aspect=40,ticks=ticks)
   cbar.set_label(r'$R$')      
#
#
#---------------------feddington----------------------------------------
#
  
#
#------------------------output to files-------------------------------
#
   oname1 = oname3d+'_radial1.png'
   oname2 = oname3d+'_radial1.ps'
   fig0.savefig(oname1, bbox_inches='tight')
   fig0.savefig(oname2, bbox_inches='tight')

   oname1 = oname3d+'_radial2.png'
   oname2 = oname3d+'_radial2.ps'
   fig1.savefig(oname1, bbox_inches='tight')
   fig1.savefig(oname2, bbox_inches='tight')
   
   oname1 = oname3d+'_radial3.png'
   oname2 = oname3d+'_radial3.ps'
   fig2.savefig(oname1, bbox_inches='tight')
   fig2.savefig(oname2, bbox_inches='tight')   

   oname1 = oname3d+'.png'
   oname2 = oname3d+'.ps'
   fig3.savefig(oname1, bbox_inches='tight')
   fig3.savefig(oname2, bbox_inches='tight')
#
#
#############################################################################
#---------------------plot everything for a certain model--------------------
#############################################################################
   
dir='../outputFILES/sc3d_nico_0001'
dir='../outputFILES/sc3d_nico_0014'
dir='../outputFILES/sc3d_nico_0020'
dir='../outputFILES/sc3d_nico_0013'
dir='../outputFILES/sc3d_nico_0013_test'
windx = 1
main(dir=dir,windx=windx,  ps=[0.,0.,1.], pe=[0.,0.,6.], )


sdum = input("Press [q] to exit.")
if sdum == 'q':
    exit()
