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
o2py.default_plot()
gc.canvas_flag=1

Zgrid=range(1,121)
Ngrid=range(1,181)
slt=Ngrid
sl=[slt] * 120

names=['sm','mn','hfb','ame03','dz','ktuy']

for i in range(0,6):

    # First initialize slice to zero
    for N in range(0,180):
        for Z in range(0,120):
            print 'N=',N,'Z=',Z
            sl[Z][N]=0
    # Now fill with data
    for row in range(0,gc.dset['nlines'][0]):
        if gc.dset['data/N'][row]>7:
            if gc.dset['data/Z'][row]>7:
                val=abs(gc.dset['data/ame'][row]-
                        gc.dset['data/'+names[i]][row])
                sl[int(gc.dset['data/Z'][row]-0.99)][int(gc.dset['data/N'][row]-0.99)]=val
    # Now plot
    plot.imshow(sl,interpolation='nearest',origin='lower',
                extent=[1,180,1,120],aspect='auto')
    plot.savefig('ex_nucmass_'+names[i]+'.eps')


