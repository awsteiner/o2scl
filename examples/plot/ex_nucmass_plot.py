"""
Plot data from ex_nucmass
"""

# In order to find o2py.py
import sys
sys.path.append('../../src/other')

import numpy as np
import matplotlib.pyplot as plot
import o2py
import math

gc=o2py.plotter()
gc.read('../ex_nucmass_table.o2')
gc.canvas_flag=1

Zgrid=range(1,121)
Ngrid=range(1,181)
sl=np.zeros(shape=(120,180))

labels=['Semi-empirical',
       'Moller et al. (1995)',
       'HFB 14',
       'HFB 21',
       'HFB 27',
       'AME (2003)',
       'Duflo and Zuker (1996)',
       'Koura et al. (2005)',
       'Dieperink et al. (2009)',
       'Wang et al. (2010)',
       'Liu et al. (2011)']

names=['se','mnmsk','hfb14','hfb21','hfb27','ame03','dz96',
       'ktuy05','dvi','ws32','ws36']

for i in range(0,11):

    (fig,ax)=o2py.default_plot()

    # First initialize slice to zero
    for N in range(0,180):
        for Z in range(0,120):
            sl[Z,N]=0
    # Now fill with data
    print 'name:',names[i],'nlines:',gc.dset['nlines'][0]
    for row in range(0,gc.dset['nlines'][0]):
        if gc.dset['data/N'][row]>7:
            if gc.dset['data/Z'][row]>7:
                val=gc.dset['data/'+names[i]][row]
                sl[int(gc.dset['data/Z'][row]-0.99),
                   int(gc.dset['data/N'][row]-0.99)]=val
    # Now plot
    cax=plot.imshow(sl,interpolation='nearest',origin='lower',
                extent=[1,180,1,120],aspect='auto',cmap='PuOr')
    cbar=plot.colorbar(cax,orientation='vertical')
    ax.text(0.55,-0.08,'N',fontsize=16,va='center',ha='center',
            transform=ax.transAxes)
    ax.text(-0.1,0.55,'Z',fontsize=16,va='center',ha='center',
            transform=ax.transAxes)
    ax.text(0.55,0.92,labels[i],fontsize=16,va='center',ha='center',
            transform=ax.transAxes)
    plot.savefig('ex_nucmass_'+names[i]+'.png')
    plot.clf()


