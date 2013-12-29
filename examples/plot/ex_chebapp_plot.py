"""
Plot data from ex_chebapp
"""

# In order to find o2mpl.py
import sys
sys.path.append('../../src/other')

import numpy as np
import matplotlib.pyplot as plot
import o2mpl

data=np.loadtxt('../ex_chebapp.out')

axes=o2mpl.default_plot()

plot.plot(data[:,0],data[:,1],color='black')
plot.plot([0.8,0.95],[-0.5,-0.5],color='black')
axes.text(0.78,-0.5,'Exact',
    fontsize=16,verticalalignment='center',horizontalalignment='right')

plot.plot(data[:,0],data[:,3],'-',color='red')
plot.plot([0.8,0.95],[-0.65,-0.65],'-',color='red')
axes.text(0.78,-0.65,'Approx. (n=50)',
    fontsize=16,verticalalignment='center',horizontalalignment='right',
    color='red')

plot.plot(data[:,0],data[:,4],':',color='blue')
plot.plot([0.8,0.95],[-0.8,-0.8],':',color='blue')
axes.text(0.78,-0.8,'Approx. (n=25)',
    fontsize=16,verticalalignment='center',horizontalalignment='right',
    color='blue')


plot.xlim([0,1])
plot.ylim([-1.1,1.1])

plot.xlabel('x',fontsize=16)

for label in axes.get_xticklabels():
    t=label.get_position()
    t2=t[0],t[1]-0.01
    label.set_position(t2)
    label.set_fontsize(16)

for label in axes.get_yticklabels():
    t=label.get_position()
    t2=t[0]-0.01,t[1]
    label.set_position(t2)
    label.set_fontsize(16)

plot.savefig('ex_chebapp_plot.png')
plot.savefig('ex_chebapp_plot.eps')
plot.show()


