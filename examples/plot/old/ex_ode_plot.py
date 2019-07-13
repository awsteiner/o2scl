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
gc.ttext(0.5,-0.07,'x')
gc.ttext(-0.12,0.5,'y')
gc.ttext(0.03,1.06,'Bessel function with non-adaptive steppers')
gc.ttext(0.4,0.87,'Prince-Dormand act. error')
gc.ttext(0.50,0.67,'Cash-Karp act. error')
gc.ttext(0.5,0.48,'Cash-Karp est. error')
gc.ttext(0.28,0.09,'Prince-Dormand est. error')

# Save figure
plot.savefig('ex_ode_bessel.png')
#plot.savefig('ex_ode_bessel.eps')

"""
Plot ex_ode_airy
"""

plot.clf()

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

gc.ttext(0.5,-0.07,'x')
gc.ttext(-0.12,0.5,'y')
gc.ttext(0.03,1.06,'Airy function with non-adaptive steppers')
gc.ttext(0.3,0.1,'Prince-Dormand act. error')
gc.ttext(0.1,0.95,'Cash-Karp act. error')
gc.ttext(0.5,0.85,'Cash-Karp est. error')
gc.ttext(0.3,0.25,'Prince-Dormand est. error')

plot.savefig('ex_ode_airy.png')
#plot.savefig('ex_ode_airy.eps')

"""
Plot ex_ode_bessel2
"""

plot.clf()

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

gc.ttext(0.5,-0.07,'x')
gc.ttext(-0.12,0.5,'y')
gc.ttext(0.11,1.06,'Bessel function with adaptive steppers')
gc.ttext_color='blue'
gc.ttext(0.04,0.52,'Prince-Dormand act. error')
gc.ttext(0.03,0.08,'Prince-Dormand est. error')
gc.ttext_color='black'
gc.ttext(0.28,0.97,'Cash-Karp act. error')
gc.ttext(0.02,0.91,'Cash-Karp est. error')
gc.ttext_color='red'
gc.ttext(0.47,0.59,'Cash-Karp(2) act. error')
gc.ttext(0.1,0.3,'Cash-Karp(2) est. error')

plot.savefig('ex_ode_bessel2.png')
#plot.savefig('ex_ode_bessel2.eps')


