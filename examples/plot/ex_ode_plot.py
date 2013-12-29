"""
Plot data from ex_ode
"""

# In order to find o2py.py
import sys
sys.path.append('../../src/other')

import numpy as np
import matplotlib.pyplot as plot
import o2py
import os

"""
Set up the plotter class
"""

gc=o2py.plotter()
gc.set('logy',1)

gc.axes=o2py.default_plot(0.14,0.12,0.05,0.08)
gc.canvas_flag=1

"""
Plot ex_ode_bessel
"""

# Cash-Karp data
gc.read_name('../ex_ode.o2','table_0')

calc=gc.dset['data/calc']
exact=gc.dset['data/exact']

diff=[]
for i in range(0,len(calc)):
    diff.append(abs(calc[i]-exact[i]))

gc.plot('x','err',color='black',ls='--')
plot.semilogy(gc.dset['data/x'],diff,color='red',ls='-.')

# Prince-Dormand data
gc.read_name('../ex_ode.o2','table_1')

calc=gc.dset['data/calc']
exact=gc.dset['data/exact']

diff=[]
for i in range(0,len(calc)):
    diff.append(abs(calc[i]-exact[i]))

gc.plot('x','err',color='purple',ls=':')
plot.semilogy(gc.dset['data/x'],diff,color='blue')

# Labels
gc.text('x',0.5,-0.07)
gc.text('y',-0.12,0.5)
gc.text('Bessel function with non-adaptive steppers',0.03,1.06)
gc.text('Prince-Dormand act. error',0.4,0.87)
gc.text('Cash-Karp act. error',0.50,0.67)
gc.text('Cash-Karp est. error',0.5,0.48)
gc.text('Prince-Dormand est. error',0.28,0.09)

# Save figure
plot.savefig('ex_ode_bessel.png')
plot.savefig('ex_ode_bessel.eps')

"""
Plot ex_ode_airy
"""

plot.clf()

gc.axes=o2py.default_plot(0.14,0.12,0.05,0.08)
plot.xlim([0.0,1.2])
plot.ylim([1.0e-17,3.0e-9])

gc.read_name('../ex_ode.o2','table_2')

calc=gc.dset['data/calc']
exact=gc.dset['data/exact']

diff=[]
for i in range(0,len(calc)):
    diff.append(abs(calc[i]-exact[i]))

gc.plot('x','err',color='black',ls='--')
plot.semilogy(gc.dset['data/x'],diff,color='red',ls='-.')

gc.read_name('../ex_ode.o2','table_3')

calc=gc.dset['data/calc']
exact=gc.dset['data/exact']
err=gc.dset['data/err']

diff=[]
for i in range(0,len(calc)):
    diff.append(abs(calc[i]-exact[i]))
abserr=[]
for i in range(0,len(calc)):
    abserr.append(abs(err[i]))

plot.semilogy(gc.dset['data/x'],abserr,color='purple',ls=':')
plot.semilogy(gc.dset['data/x'],diff,color='blue')

gc.text('x',0.5,-0.07)
gc.text('y',-0.12,0.5)
gc.text('Airy function with non-adaptive steppers',0.03,1.06)
gc.text('Prince-Dormand act. error',0.3,0.1)
gc.text('Cash-Karp act. error',0.1,0.95)
gc.text('Cash-Karp est. error',0.5,0.85)
gc.text('Prince-Dormand est. error',0.3,0.25)

plot.savefig('ex_ode_airy.png')
plot.savefig('ex_ode_airy.eps')

"""
Plot ex_ode_bessel2
"""

plot.clf()

gc.axes=o2py.default_plot(0.14,0.12,0.05,0.08)
plot.xlim([0.0,10])
plot.ylim([1.0e-10,3.0e-6])

gc.read_name('../ex_ode.o2','table_4')

calc=gc.dset['data/calc']
exact=gc.dset['data/exact']
err0=gc.dset['data/err0']

diff=[]
for i in range(0,len(calc)):
    diff.append(abs(calc[i]-exact[i]))
abserr=[]
for i in range(0,len(calc)):
    abserr.append(abs(err0[i]))

plot.semilogy(gc.dset['data/x'],abserr,color='black',ls='--')
plot.semilogy(gc.dset['data/x'],diff,color='black')

gc.read_name('../ex_ode.o2','table_5')

calc=gc.dset['data/calc']
exact=gc.dset['data/exact']
err0=gc.dset['data/err0']

diff=[]
for i in range(0,len(calc)):
    diff.append(abs(calc[i]-exact[i]))
abserr=[]
for i in range(0,len(calc)):
    abserr.append(abs(err0[i]))

plot.semilogy(gc.dset['data/x'],abserr,color='red',ls='--')
plot.semilogy(gc.dset['data/x'],diff,color='red')

gc.read_name('../ex_ode.o2','table_6')

calc=gc.dset['data/calc']
exact=gc.dset['data/exact']
err0=gc.dset['data/err0']

diff=[]
for i in range(0,len(calc)):
    diff.append(abs(calc[i]-exact[i]))
abserr=[]
for i in range(0,len(calc)):
    abserr.append(abs(err0[i]))

plot.semilogy(gc.dset['data/x'],abserr,color='blue',ls='--')
plot.semilogy(gc.dset['data/x'],diff,color='blue')

gc.text('x',0.5,-0.07)
gc.text('y',-0.12,0.5)
gc.text('Bessel function with adaptive steppers',0.11,1.06)
gc.text('Prince-Dormand act. error',0.04,0.52,color='blue')
gc.text('Prince-Dormand est. error',0.03,0.08,color='blue')
gc.text('Cash-Karp act. error',0.28,0.97)
gc.text('Cash-Karp est. error',0.02,0.91)
gc.text('Cash-Karp(2) act. error',0.47,0.59,color='red')
gc.text('Cash-Karp(2) est. error',0.1,0.3,color='red')

plot.savefig('ex_ode_bessel2.png')
plot.savefig('ex_ode_bessel2.eps')


