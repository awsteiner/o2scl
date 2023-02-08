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
from sklearn.gaussian_process import GaussianProcessClassifier as GPC
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
from sklearn.gaussian_process.kernels import RBF, DotProduct
from sklearn.gaussian_process.kernels import RationalQuadratic
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import QuantileTransformer
from sklearn.preprocessing import MinMaxScaler
import tensorflow as tf

class emu_dnn:

    def __init__(self):
        # Class variable for DNN
        self.dnn=0
        self.train_table=0
        self.input_list=[]
        self.output_list=[]
        self.SS1=QuantileTransformer()
        self.SS2=QuantileTransformer()
        self.tracker=0

    def read_data(self,nd_in,hdf_file,ix_log_wgt,col_list,verbose):
        """
        Desc.
        """

        link=o2sclpy.linker()
        link.link_o2scl()
        
        new_data=o2sclpy.table(link)
        data=o2sclpy.table(link)
        train_table=o2sclpy.table(link)

        for i in range(0, nd_in):
            self.input_list.append(col_list[i])

        for i in range(nd_in, len(col_list)):
           self.output_list.append(col_list[i])

        print("emu_dnn::read_data(): Reading file:",hdf_file)
        hf=o2sclpy.hdf_file(link)
        hf.open(hdf_file)
        new_data.clear_table()
        o2sclpy.hdf_input_table(link,hf,new_data)
        hf.close()
        print("emu_dnn::read_data(): Done reading file.")
        print("emu_dnn::read_data(): Table in file", hdf_file,"has",
              new_data.get_nlines(),"lines.")

        train_table.add_table(data)
        self.train_table=train_table

        print("emu_dnn::read_data(): table now has",
              self.train_table.get_nlines(),
              "lines.")

        train_data=[]
        output_data=[]

        if (self.tracker!=0):
            train_data.clear()
            output_data.clear()
            print(np.shape(train_data),np.shape(output_data))
            print(len(self.input_list),len(self.output_list))
            train_data=[]
            output_data=[]
        
        for i in self.input_list:
            train_data.append(np.array(self.train_table[i]
                                       [0:self.train_table.get_nlines()]))
        #print("Pymodule: ",np.shape(np.array(data1.__getitem__("log10_Tcn"))))
        train_data=self.SS1.fit_transform(np.asarray(train_data).transpose())
         
        for i in self.output_list:
            output_data.append(np.array(self.train_table[i]
                                        [0:self.train_table.get_nlines()]))
        output_data=self.SS2.fit_transform(np.asarray(output_data).transpose())
        
        #self.output_data=self.SS2.fit_transform(output_data)
        #print("emu_dnn::read_data(): log_wgt: ", self.log_wgt)
        
        print("emu_dnn::read_data(): input data: ", np.shape(train_data))
        print("emu_dnn::read_data(): output data: ", np.shape(output_data))
        #print(output_data)

        x_train, x_test, y_train, y_test=train_test_split(
            train_data, output_data, test_size=0.15)
        print("emu_dnn::read_data(): Training array : ", x_train.shape)
        print("emu_dnn::read_data(): Target array : ", y_train.shape)

        def custom_loss(y_true, y_pred):
            # first output neuron
            mu=y_pred[:,:1]
            # second output neuron
            log_sig=y_pred[:,1:]
            # undo the log
            sig=tf.exp(log_sig)  

            return tf.reduce_mean(2*log_sig+((y_true-mu)/sig)**2)

        print("emu_dnn::read_data(): Training DNN model.")
        model=tf.keras.Sequential(
            [
                tf.keras.layers.Dense(
                    302,input_shape=(155,), activation='tanh'),
                tf.keras.layers.Dense(32, activation='tanh'),
                tf.keras.layers.Dense(21, activation='sigmoid')
            ])
        print(model.summary())

        model.compile(loss='mean_squared_error',
                      optimizer='adam', metrics=['accuracy'])
        model.fit(x_train,y_train,batch_size=128,epochs=150,
                  validation_data=(x_test,y_test),verbose=0)

        print("emu_dnn::read_data(): Training done.")
        print("emu_dnn::read_data(): Test Score: [loss, accuracy]: ",
              model.evaluate(x_test, y_test, verbose=0))
        self.dnn=model

        return
    
    def upTrain(self, vec):
        """
        Desc.
        """
        
        self.train_table.line_of_data(vec)
        self.tracker=self.tracker + 1
        print("emu_dnn::upTrain(): Tracker: ", self.tracker)
        if self.tracker == 300:
            self.train()
            self.tracker=0

        return
    
    def predict(self, trial_input):
        """
        Desc.
        """

        print("emu_dnn::predict(): starting predict")
        trial_input=np.asarray(trial_input)
        trial_input=trial_input.reshape(1, len(self.input_list))
        trial_input=self.SS1.transform(trial_input)
        
        # Predict gives the output and the std_dev in normal scale
        # no need to inverse transform anything
        
        re_predicted=[]

        predicted=self.dnn.predict(trial_input, verbose=1)
        trial_input.__del__()
        predicted=self.SS2.inverse_transform(predicted)
    
        for i in range(0, len(self.output_list)):
            #print(self.output_list[i], predicted[0][i])
            re_predicted.append(predicted[0][i])
            re_predicted.append(0)

        predicted.__del__()
        return re_predicted

if __name__ == '__main__':
    ednn=emu_dnn();
    ednn.read_data(2,'emu_data.o2',0,'z,x,y,d',1)
    print(ednn.predict([1,2]))
    
