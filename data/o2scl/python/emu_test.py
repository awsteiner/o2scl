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
#import sklearn
#from sklearn.gaussian_process import GaussianProcessRegressor
#from sklearn.gaussian_process.kernels import RBF
import gc

class emu_py:

    link=0
    gpr=0

    def __init__(self):
        print('init')
        gc.disable()
    
    def train(self,num_params,filename,ix_log_wgt,col_list,verbose):

        if verbose>0:
            print('num_params:',num_params)
            print('filename:',filename)
            print('ix_log_wgt:',ix_log_wgt)
            print('col_list:',col_list)
            print('verbose:',num_params)
        
        return

    def point(self,v):
        return v
    
    def __del__(self):
        print('del')
        gc.enable()
    
    
