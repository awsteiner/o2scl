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

count=0

def fun(x):
    global count
    y=[]
    for i in range(0,len(x)):
        y.append(x[i]*numpy.pi)
        print('%d %7.6e %7.6e' % (count,x[i],y[i]))
    count=count+1
    return y

class mft:

    count=0
    
    def fun(self,x):
        y=[]
        for i in range(0,len(x)):
            y.append(x[i]*numpy.pi)
            print('%d %7.6e %7.6e' % (self.count,x[i],y[i]))
        self.count=self.count+1
        return y
            
if __name__ == '__main__':
    print(fun([1,2]))
    mft=mft()
    print(mft.fun([1,2]))

    
