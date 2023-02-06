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
import sklearn
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
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

        col_list2=col_list.split(',')
        lw_col=col_list2[0]
        print('col_list2,lw_col:',col_list2,lw_col)
            
        param_list=[]
        data_list=[]
            
        for i in range(0,len(col_list2)):
            if i<num_params:
                param_list.append(col_list2[i])
            else:
                data_list.append(col_list2[i])
        print('param_list:',param_list)
        print('data_list:',data_list)
            
        if True:
            
            import o2sclpy
            
            self.link=o2sclpy.linker()
            self.link.link_o2scl()
        
            tab=o2sclpy.table(self.link)
    
            # Open the file and read into tab2
            hf=o2sclpy.hdf_file(self.link)
            hf.open(filename,False,True)
            tab=o2sclpy.table(self.link)
            o2sclpy.hdf_input_table(self.link,hf,tab)
            hf.close()
    
            print('nlines',tab.get_nlines())
    
            x=numpy.zeros((tab.get_nlines(),num_params))
            for i in range(0,tab.get_nlines()):
                for j in range(0,num_params):
                    x[i,j]=tab.get(param_list[j],i)
            print('set x')
                    
            y=numpy.zeros((tab.get_nlines(),len(data_list)))
            for i in range(0,tab.get_nlines()):
                for j in range(0,len(data_list)):
                    y[i,j]=tab.get(data_list[j],i)
            print('set y')
            
        else:
            
            import h5py

            tfile=h5py.File(filename)
            print('nlines',tfile['tab/nlines'][0])
            
            x=numpy.zeros((tfile['tab/nlines'][0],num_params))
            for i in range(0,tfile['tab/nlines'][0]):
                for j in range(0,num_params):
                    x[i,j]=tfile['tab/data/'+param_list[j]][i]
            print('set x2')
                    
            y=numpy.zeros((tfile['tab/nlines'][0],len(data_list)))
            for i in range(0,tfile['tab/nlines'][0]):
                for j in range(0,len(data_list)):
                    y[i,j]=tfile['tab/data/'+data_list[j]][i]
            print('set y2')
            
        # Perform the GP fit
        kernel=RBF(1.0,(1.0e0,1.0e2))
        self.gpr=GaussianProcessRegressor(kernel=kernel).fit(x,y)
        
        return

    def point(self,v):
        yp=self.gpr.predict([v])
        return yp[0].tolist()

    def __del__(self):
        print('del')
        gc.enable()
        
if __name__ == '__main__':
    epy=emu_py();
    epy.train(2,'emu_data.o2',0,'z,x,y,d',1)
    epy.point([1,2])
    
    
