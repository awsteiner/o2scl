"""
Plot data from ex_fptr
"""

# In order to find o2py.py
import sys
sys.path.append('../../src/other')

import numpy as np
import matplotlib.pyplot as plot
import o2py
import os
import math

def den(file,slice,title,prefix):
    fig=plot.figure(1,figsize=(6.4,6.4))
    fig.set_facecolor('white')
#
    axes=plot.axes([0.10,0.085,0.93,0.85])
    axes.minorticks_on()
    axes.tick_params('both',length=12,width=1,which='major')
    axes.tick_params('both',length=5,width=1,which='minor')
    plot.grid(False)
#
    dset=o2py.h5read_first_type(file,'table3d')
    plot.xlabel(r'$\log_{10}~\rho$',fontsize=16)
    plot.ylabel(r'$\mathrm{T}~(\mathrm{MeV})$',fontsize=16)
#
    sl=dset[slice].value
    sl=sl.transpose()
    xgrid=dset['xval'].value
    ygrid=dset['yval'].value
#
    for i in range(0,len(xgrid)):
        xgrid[i]=math.log(xgrid[i],10)
#
    img=plot.imshow(sl,interpolation='nearest',origin='lower',
                    extent=[xgrid[0],xgrid[len(xgrid)-1],ygrid[0],
                            ygrid[len(ygrid)-1]],aspect='auto')
    cbar=fig.colorbar(img)
#
    axes.text(0.5,1.03,title,transform=axes.transAxes,
              fontsize=16,verticalalignment='center',
              horizontalalignment='center')
#
    plot.savefig(prefix+'.png')
    plot.savefig(prefix+'.eps')
#    plot.show()
    plot.clf()

def con(file,slice,title,conlist,prefix):
    fig=plot.figure(1,figsize=(6.4,6.4))
    fig.set_facecolor('white')
#
    axes=plot.axes([0.10,0.085,0.85,0.85])
    axes.minorticks_on()
    axes.tick_params('both',length=12,width=1,which='major')
    axes.tick_params('both',length=5,width=1,which='minor')
    plot.grid(False)
#
    dset=o2py.h5read_first_type(file,'table3d')
    plot.xlabel(r'$\log_{10}~\rho$',fontsize=16)
    plot.ylabel(r'$\mathrm{T}~(\mathrm{MeV})$',fontsize=16)
#
    sl=dset[slice].value
    sl=sl.transpose()
    xgrid=dset['xval'].value
    ygrid=dset['yval'].value
#
    for i in range(0,len(xgrid)):
        xgrid[i]=math.log(xgrid[i],10)
#
    conts=plot.contour(xgrid,ygrid,sl,conlist)
    plot.clabel(conts,inline=1, fontsize=10)
#
    axes.text(0.5,1.03,title,transform=axes.transAxes,
              fontsize=16,verticalalignment='center',
              horizontalalignment='center')
#
    plot.savefig(prefix+'.png')
    plot.savefig(prefix+'.eps')
#    plot.show()
    plot.clf()

plot.rc('text',usetex=True)
plot.rc('font',family='serif')
plot.rcParams['lines.linewidth']=0.5
#
den('../oo_ls220_Atab.o2','data/A_Ye0.5',
    'Mass number for $Y_{e}$=0.5','oo_ls220_A0.5')
den('../oo_ls220_Atab.o2','data/A_Ye0.1',
    'Mass number for $Y_{e}$=0.1','oo_ls220_A0.1')
con('../oo_ls220_Atab.o2','data/mun_Ye0.5',
    '$\mu_{\mathrm{neutron}}$ for $Y_{e}$=0.5',
    [-5000,-2500,-1000,-500,-100,500,2000,5000],
    'oo_ls220_mun0.5')
con('../oo_ls220_Atab.o2','data/mun_Ye0.1',
    '$\mu_{\mathrm{neutron}}$ for $Y_{e}$=0.1',
    [-5000,-2500,-1000,-500,-100,500,2000,5000],
    'oo_ls220_mun0.1')
con('../oo_ls220_Atab.o2','data/mup_Ye0.5',
    '$\mu_{\mathrm{proton}}$ for $Y_{e}$=0.5',
    [-5000,-2500,-1000,-500,-100,500,2000,5000],
    'oo_ls220_mup0.5')
con('../oo_ls220_Atab.o2','data/mup_Ye0.1',
    '$\mu_{\mathrm{proton}}$ for $Y_{e}$=0.1',
    [-5000,-2500,-1000,-500,-100,500,2000,5000],
    'oo_ls220_mup0.1')
den('../oo_ls220_Atab.o2','data/Eint_Ye0.5',
    'Energy for $Y_{e}$=0.5','oo_ls220_Eint0.5')
den('../oo_ls220_Atab.o2','data/Eint_Ye0.1',
    'Energy for $Y_{e}$=0.1','oo_ls220_Eint0.1')
den('../oo_ls220_Atab.o2','data/Xn_Ye0.5',
    'Neutron fraction for $Y_{e}$=0.5','oo_ls220_Xn0.5')
den('../oo_ls220_Atab.o2','data/Xn_Ye0.1',
    'Neutron fraction for $Y_{e}$=0.1','oo_ls220_Xn0.1')
den('../oo_ls220_Atab.o2','data/Xp_Ye0.5',
    'Proton fraction for $Y_{e}$=0.5','oo_ls220_Xp0.5')
den('../oo_ls220_Atab.o2','data/Xp_Ye0.1',
    'Proton fraction for $Y_{e}$=0.1','oo_ls220_Xp0.1')
den('../oo_ls220_Atab.o2','data/Xalpha_Ye0.5',
    'Alpha fraction for $Y_{e}$=0.5','oo_ls220_Xalpha0.5')
den('../oo_ls220_Atab.o2','data/Xalpha_Ye0.1',
    'Alpha fraction for $Y_{e}$=0.1','oo_ls220_Xalpha0.1')
den('../oo_ls220_Atab.o2','data/Xnuclei_Ye0.5',
    'Heavy nucleus fraction for $Y_{e}$=0.5','oo_ls220_Xnuclei0.5')
den('../oo_ls220_Atab.o2','data/Xnuclei_Ye0.1',
    'Heavy nucleus fraction for $Y_{e}$=0.1','oo_ls220_Xnuclei0.1')




