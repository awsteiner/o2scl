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

os.system('cat ../ex_fptr.scr | grep -v t > exftemp1')
os.system('acol -generic ../ex_fptr.out -internal exftemp2.o2')

gc=o2py.plotter()
gc.read('exftemp2.o2')
gc.plot('x','y',color='red')
gc.text('x',0.5,-0.07)
gc.text('y',-0.1,0.5)

data=np.loadtxt('exftemp1')

plot.plot(data[0:8,0],data[0:8,1],marker='s')

plot.savefig('ex_fptr_plot.png')
plot.savefig('ex_fptr_plot.eps')

os.system('rm exftemp1')
os.system('rm exftemp2.o2')

# gc.show()


