"""
  -------------------------------------------------------------------
  
  Copyright (C) 2022-2023, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
"""

import numpy

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF

count=0

def fun(x):
    global count
    y=[]
    for i in range(0,len(x)):
        y.append(x[i]*numpy.pi)
        print('%d %7.6e %7.6e' % (count,x[i],y[i]))
    count=count+1
    return y

gpr=0

def train(x0):
    global gpr
    x=[]
    y=[]
    for i in range(0,100):
        x.append(numpy.sin(i));
        y.append(numpy.cos(x[i]))
    # Perform the GP fit
    kernel=RBF(1.0,(1.0e0,1.0e2))
    x2=numpy.asarray(x)
    gpr=GaussianProcessRegressor(kernel=kernel).fit(x2.reshape(-1,1),y)
    return x0

def feval(p):
    global gpr
    p2=gpr.predict([[p[0]]])
    return p2.tolist()

if __name__ == '__main__':
    fun([1,2])
    print(train([3,4]))
    print('%7.6e %7.6e' % (feval([0.5])[0],numpy.cos(0.5)))

    
    
