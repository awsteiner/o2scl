"""
Plot data from ex_eos_had
"""

def text_abs(axes,x,y,str,**kwargs):
    axes.text(x,y,str,transform=axes.transAxes,
              fontsize=16,va='center',ha='center',**kwargs)

def text(axes,x,y,str,**kwargs):
    axes.text(x,y,str,fontsize=16,va='center',ha='center',**kwargs)

def left_text(axes,x,y,str,**kwargs):
    axes.text(x,y,str,fontsize=16,va='center',ha='left',**kwargs)
              

# In order to find o2py.py
import sys
sys.path.append('../../src/other')

import numpy as np
import matplotlib.pyplot as plot
import o2py
import os

plot.rc('text',usetex=True)

"""
Plot EOS data
"""

if len(sys.argv)<2 or sys.argv[1]=='eos':
    
    axes=o2py.default_plot(0.14,0.12,0.04,0.04)
    
    plot.xlim([0.0,2000])
    plot.ylim([0.1,2000])
    
    # Read APR EOS data
    dset=o2py.h5read_first_type('../ex_eos_had_apr_nstar.o2','table')
    
    # Multiply both columns by hbar*c
    ed=dset['data/ed'].value
    pr=dset['data/pr'].value

    for i in range(0,len(ed)):
        ed[i]=ed[i]*197.33
        pr[i]=pr[i]*197.33

    # Plot
    plot.semilogy(ed,pr,color='black',ls='-')
    
    # Read SLy4 EOS data
    dset=o2py.h5read_first_type('../skyrme_SLy4_eos.o2','table')

    # Multiply both columns by hbar*c
    ed=dset['data/ed'].value
    pr=dset['data/pr'].value

    for i in range(0,len(ed)):
        ed[i]=ed[i]*197.33
        pr[i]=pr[i]*197.33

    # Plot
    plot.semilogy(ed,pr,color='black',ls=':')

    # Labels
    text_abs(axes,0.50,-0.09,r'$\varepsilon$ ($\mathrm{MeV}/\mathrm{fm}^{3}$)')
    text_abs(axes,-0.09,0.50,r'P ($\mathrm{MeV}/\mathrm{fm}^{3}$)',rotation=90)

    # Legend
    plot.plot([1000,1400],[2,2],color='black',ls='-')
    plot.plot([1000,1400],[1,1],color='black',ls=':')
    left_text(axes,1420,2,r'APR')
    left_text(axes,1420,1,r'SLy4')

    # Save 
    plot.savefig('ex_eos_had_eos.png',dpi=60)
    plot.show()

if len(sys.argv)>1 and sys.argv[1]=='gibbs':
    
    axes=o2py.default_plot(0.14,0.12,0.04,0.04)
    
    plot.xlim([0.19,0.25])
    plot.ylim([0.0,1])
    
    # Read APR EOS data
    dset=o2py.h5read_first_type('../ex_eos_had_apr_nstar.o2','table')
    
    # Multiply both columns by hbar*c
    nb=dset['data/nb'].value
    chi=dset['data/chi'].value
    nn=dset['data/nn'].value
    np=dset['data/np'].value
    nn2=dset['data/nn2'].value
    np2=dset['data/np2'].value

    for i in range(0,len(nb)):
        nn[i]=nn[i]/nb[i]*chi[i]
        np[i]=np[i]/nb[i]*chi[i]
        nn2[i]=nn2[i]/nb[i]*(1-chi[i])
        np2[i]=np2[i]/nb[i]*(1-chi[i])

    # Plot
    plot.plot(nb,nn,color='black',ls='-')
    plot.plot(nb,np,color='black',ls=':')
    plot.plot(nb,nn2,color='black',ls='--')
    plot.plot(nb,np2,color='black',ls='-.')
    
    # Labels
    text_abs(axes,0.50,-0.09,r'$n_B$ ($\mathrm{fm}^{-3}$)')
    text_abs(axes,-0.09,0.50,r'$x_i$',rotation=90)

    # Legend
    text(axes,0.205,0.8,r'$x_{nL}$')
    text(axes,0.195,0.1,r'$x_{pL}$')
    text(axes,0.235,0.8,r'$x_{nH}$')
    text(axes,0.245,0.1,r'$x_{pH}$')

    # Save 
    plot.savefig('ex_eos_had_gibbs.png',dpi=60)
    plot.show()

if len(sys.argv)>1 and sys.argv[1]=='mvsr':

    """
    Plot M-R data
    """

    axes=o2py.default_plot(0.14,0.12,0.04,0.04)

    plot.xlim([9,16])
    plot.ylim([0,2.3])

    # Read APR M-R data
    dset=o2py.h5read_first_type('../ex_eos_had_apr_mvsr.o2','table')

    gm=dset['data/gm'].value
    r=dset['data/r'].value

    for i in range(1,len(gm)):
        if gm[i]<gm[i-1]:
            gm[i]=gm[i-1]
            r[i]=r[i-1]

    plot.plot(r,gm,color='black',ls='-')

    # Read SLy4 M-R data
    dset=o2py.h5read_first_type('../skyrme_SLy4_mvsr.o2','table')
    
    gm=dset['data/gm'].value
    r=dset['data/r'].value
    
    for i in range(1,len(gm)):
        if gm[i]<gm[i-1]:
            gm[i]=gm[i-1]
            r[i]=r[i-1]

    plot.plot(r,gm,color='black',ls=':')

    # Labels
    text_abs(axes,0.50,-0.09,r'R (km)')
    text_abs(axes,-0.09,0.50,r'M ($\mathrm{M}_{\odot}$)',rotation=90)
    
    plot.plot([13,14.4],[0.9,0.9],color='black',ls='-')
    plot.plot([13,14.4],[0.7,0.7],color='black',ls=':')
    left_text(axes,14.6,0.9,r'APR')
    left_text(axes,14.6,0.7,r'SLy4')
    
    # Save figure
    plot.savefig('ex_eos_had_mvsr.png',dpi=60)
    plot.show()
    
