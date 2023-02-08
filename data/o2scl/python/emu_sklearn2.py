"""
  -------------------------------------------------------------------
  
  Copyright (C) 2022-2023, Satyajit Roy and Andrew W. Steiner
  
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

import o2sclpy
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from sklearn.gaussian_process.kernels import RBF, DotProduct
from sklearn.gaussian_process.kernels import RationalQuadratic
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import QuantileTransformer
from sklearn.preprocessing import MinMaxScaler

class emu_gpr:
    """
    Desc
    """

    def __init__(self):
        # Initialize class variables
        self.gpr=0
        self.kernel=1.0*RBF(1)
        self.train_table=0
        self.input_list=[]
        self.output_list=[]
        self.SS1=StandardScaler()

    def read_data(self,nd_in,hdf_file,ix_log_wgt,col_list,verbose):
        """
        Desc.
        """
        #ds_name, nd_in, param_list):

        link=o2sclpy.linker()
        link.link_o2scl()
        self.train_table=o2sclpy.table(link)
        col_list2=col_list.split(',')
        
        for i in range(0,nd_in):
            self.input_list.append(col_list2[i])

        for i in range(nd_in,len(col_list2)):
           self.output_list.append(col_list2[i])

        print("emu_gpr::read_data(): Filename:",hdf_file)
        hf=o2sclpy.hdf_file(link)
        hf.open(hdf_file)
        self.train_table.clear_table()
        o2sclpy.hdf_input_table(link,hf,self.train_table)
        hf.close()
        print("emu_gpr::read_data(): Table has",
              self.train_table.get_nlines(),"lines.")

        input_data=[]
        for i in self.input_list:
            input_data.append(np.array(self.train_table.__getitem__(i)))

        input_data=self.SS1.fit_transform(np.asarray(input_data).transpose())

        train_data =[]
        output_data=[]

        for i in self.input_list:
            train_data.append(np.array(self.train_table[i]
                                       [0:self.train_table.get_nlines()]))
        train_data=self.SS1.fit_transform(np.asarray(train_data).transpose())
         
        for i in self.output_list:
            output_data.append(np.array(self.train_table[i]
                                        [0:self.train_table.get_nlines()]))
        output_data=np.asarray(output_data).transpose()
        print("emu_gpr::read_data(): Input data:",
              np.shape(train_data))
        print("emu_gpr::read_data(): Output data:",
              np.shape(output_data))

        x_train, x_test, y_train, y_test=train_test_split(
            train_data, output_data, test_size=0.15)
        print("emu_gpr::read_data(): Training array:",
              x_train.shape)
        print("emu_gpr::read_data(): Target array:",
              y_train.shape)

        print("emu_gpr::read_data(): Training GPR model.")
        
        self.gpr=GPR(kernel=self.kernel,random_state=0,
                     n_restarts_optimizer=5, 
                     normalize_y=True).fit(x_train,y_train)
        
        print("emu_gpr::read_data(): Training done. Score:",
              self.gpr.score(x_test,y_test))
        
        return

    def predict(self, trial_input):
        """
        Desc.
        """

        print("emu_gpr::predict(): Starting.")
        trial_input=np.asarray(trial_input)
        trial_input=trial_input.reshape(1,len(self.input_list))
        trial_input=self.SS1.transform(trial_input)
        
        predicted, std_dev=self.gpr.predict(
            trial_input, return_std=True, return_cov=False)
        #predicted=self.SS2.inverse_transform(predicted)

        re_predicted=[]
        for i in range(0, len(self.output_list)):
            #print(self.output_list[i], predicted[0][i])
            re_predicted.append(predicted[0][i])
        re_predicted.append(std_dev[0][i])
             
        return re_predicted

if __name__ == '__main__':
    egpr=emu_gpr();
    egpr.read_data(2,'emu_data.o2',0,'z,x,y,d',1)
    print(egpr.predict([1,2]))
