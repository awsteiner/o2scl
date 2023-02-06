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
import h5py
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

class emu_gpr:

    def __init__(self):
        # Initialize class variables
        self.gpr = 0
        self.gpc = 0
        self.kernel = 1.0 * RBF(1)
        self.train_table = 0
        self.input_list = []
        self.output_list = []
        self.SS1 = StandardScaler()
        self.SS2 = StandardScaler()
        self.tracker = 0

    def read_data(self, file_list, ds_name, nd_in, param_list):
        """
        Desc.
        """

        link=o2sclpy.linker()
        link.link_o2scl()
        new_data = o2sclpy.table(link)
        data = o2sclpy.table(link)
        train_table = o2sclpy.table(link)
        p=o2sclpy.o2graph_plotter()

        for i in range(0, nd_in):
            self.input_list.append(param_list[i])

        for i in range(nd_in, len(param_list)):
           self.output_list.append(param_list[i])

        for hdf_file in file_list:
            print("PyModule : Reading File: ", hdf_file)
            hf=o2sclpy.hdf_file(link)
            hf.open(hdf_file)
            new_data.clear_table()
            o2sclpy.hdf_input_table(link,hf,new_data)
            hf.close()
            print("PyModule : Done Reading File: ")
            print("Pymodule : Table in file ", hdf_file, " has lines ",
                  new_data.get_nlines())

            # Delete any rows with a small log_wgt, which are emulated, or
            # have a value of "mult" which is zero
            print("Pymodule : table has ", data.get_nlines(), " lines.")

            if new_data.is_column("emulated"):
                if new_data.is_column("mult"):
                    new_data.delete_rows_func("emulated>0.5 || mult<0.5")
                else:
                    new_data.delete_rows_func("emulated>0.5")
            else:
                if new_data.is_column("mult"):
                    new_data.delete_rows_func("mult<0.5")
            
            data.add_table(new_data)
            print("Pymodule : table now has ", data.get_nlines(), " lines.")

        
        input_data = []
        for i in self.input_list:
            input_data.append(np.array(data.__getitem__(i)))
        input_data = self.SS1.fit_transform(np.asarray(input_data).transpose())
        # Only transform the input data, normalize_y is used in GPR itself
        # otherwise we lose data

        log_wgt = np.array(data.__getitem__("log_wgt")).transpose()
        log_wgt = log_wgt.reshape(len(log_wgt),1)
        for value in log_wgt:
            if value[0]<=(-80): value[0] = 0
            else: value[0] = 1

        print("PyModule 1: input data: ", np.shape(input_data))
        print("PyModule : log_wgt: ", np.shape(log_wgt))

        print("Pymodule : starting train")
        self.gpc = GPC(kernel=self.kernel).fit(input_data, log_wgt)

        data.delete_rows_func("log_wgt<(-80)")
        print("Pymodule : table now has ", data.get_nlines(), " lines.")

        train_table.add_table(data)
        self.train_table = train_table

        return

    def train(self):
        """
        Desc.
        """
        print("Pymodule : table now has ", self.train_table.get_nlines(),
              " lines.")

        train_data =[]
        output_data = []

        if (self.tracker !=0):
            train_data.clear()
            output_data.clear()
            print( np.shape(train_data),  np.shape(output_data))
            print(len(self.input_list), len(self.output_list))
            train_data = []
            output_data = []
        
        for i in self.input_list:
            train_data.append(np.array(self.train_table[i]
                                       [0:self.train_table.get_nlines()]))
        #print("Pymodule: ",np.shape(np.array(data1.__getitem__("log10_Tcn"))))
        train_data = self.SS1.fit_transform(np.asarray(train_data).transpose())
         
        for i in self.output_list:
            output_data.append(np.array(self.train_table[i]
                                        [0:self.train_table.get_nlines()]))
        output_data = np.asarray(output_data).transpose()
        #self.output_data = self.SS2.fit_transform(output_data)
        #print("PyModule : log_wgt: ", self.log_wgt)
        print("PyModule 2: input data: ", np.shape(train_data))
        print("PyModule 2: output data: ", np.shape(output_data))
        #print(output_data)

        x_train, x_test, y_train, y_test = train_test_split(
            train_data, output_data, test_size=0.15)
        print("PyModule : Training array : ", x_train.shape)
        print("PyModule : Target array : ", y_train.shape)

        print("PyModule : Training GPR model.")
        
        self.gpr = GPR(kernel=self.kernel, random_state=0,
                       n_restarts_optimizer= 5, 
                       normalize_y=True).fit(x_train, y_train)

        print("PyModule : Training done")
        print("PyModule: Test Score: ", self.gpr.score(x_test, y_test))
        
        return

    def upTrain(self, vec):
        """
        Desc.
        """

        self.train_table.line_of_data(vec)
        self.tracker = self.tracker + 1
        print("PyModule : Tracker: ", self.tracker)
        if self.tracker == 2:
            self.train()
            self.tracker = 0

        return

    def predict(self, trial_input):
        """
        Desc.
        """

        print("PyModule : starting predict")
        trial_input = np.asarray(trial_input)
        trial_input = trial_input.reshape(1, len(self.input_list))
        trial_input = self.SS1.transform(trial_input)
        # Predict gives the output and the std_dev in normal scale
        # no need to inverse transform anything
        re_predicted = []
        pred_lw = self.gpc.predict(trial_input)
        #print("pred_lw", pred_lw)
        if pred_lw>0:
            predicted, std_dev = self.gpr.predict(
                trial_input, return_std=True, return_cov=False)
            #predicted = self.SS2.inverse_transform(predicted)
            
            for i in range(0, len(self.output_list)):
            #print(self.output_list[i], predicted[0][i])
                re_predicted.append(predicted[0][i])
                re_predicted.append(std_dev[0][i])
        
        else:
            re_predicted.append(-800.0)
            for i in range(1, len(self.output_list)-1):
                re_predicted.append(0.0)
             
        return re_predicted


class emu_dnn:

    def __init__(self):
        # Class variable for DNN
        self.dnn = 0
        self.train_table = 0
        self.input_list = []
        self.output_list = []
        self.SS1 = QuantileTransformer()
        self.SS2 = QuantileTransformer()
        self.tracker = 0

    def read_data(self, file_list, ds_name, nd_in, param_list):
        """
        Desc.
        """

        link=o2sclpy.linker()
        link.link_o2scl()
        new_data = o2sclpy.table(link)
        data = o2sclpy.table(link)
        train_table = o2sclpy.table(link)
        p=o2sclpy.o2graph_plotter()

        for i in range(0, nd_in):
            self.input_list.append(param_list[i])

        for i in range(nd_in, len(param_list)):
           self.output_list.append(param_list[i])

        for hdf_file in file_list:
            print("PyModule : Reading File: ", hdf_file)
            hf=o2sclpy.hdf_file(link)
            hf.open(hdf_file)
            new_data.clear_table()
            o2sclpy.hdf_input_table(link,hf,new_data)
            hf.close()
            print("PyModule : Done Reading File: ")
            print("Pymodule : Table in file ", hdf_file, " has lines ",
                  new_data.get_nlines())

            # Delete any rows with a small log_wgt, which are emulated, or
            # have a value of "mult" which is zero
            print("Pymodule : table has ", data.get_nlines(), " lines.")
            """
            if new_data.is_column("emulated"):
                if new_data.is_column("mult"):
                    new_data.delete_rows_func("emulated>0.5 || mult<0.5")
                else:
                    new_data.delete_rows_func("emulated>0.5")
            else:
                if new_data.is_column("mult"):
                    new_data.delete_rows_func("mult<0.5")
            """
            data.add_table(new_data)
            print("Pymodule : table now has ", data.get_nlines(), " lines.")

        print("Pymodule : table now has ", data.get_nlines(), " lines.")

        train_table.add_table(data)
        self.train_table = train_table

        return
    
    def train(self):
        """
        Desc.
        """
        print("Pymodule : table now has ", self.train_table.get_nlines(),
              " lines.")

        train_data =[]
        output_data = []

        if (self.tracker !=0):
            train_data.clear()
            output_data.clear()
            print( np.shape(train_data),  np.shape(output_data))
            print(len(self.input_list), len(self.output_list))
            train_data = []
            output_data = []
        
        for i in self.input_list:
            train_data.append(np.array(self.train_table[i]
                                       [0:self.train_table.get_nlines()]))
        #print("Pymodule: ",np.shape(np.array(data1.__getitem__("log10_Tcn"))))
        train_data = self.SS1.fit_transform(np.asarray(train_data).transpose())
         
        for i in self.output_list:
            output_data.append(np.array(self.train_table[i]
                                        [0:self.train_table.get_nlines()]))
        output_data=self.SS2.fit_transform(np.asarray(output_data).transpose())
        #self.output_data = self.SS2.fit_transform(output_data)
        #print("PyModule : log_wgt: ", self.log_wgt)
        print("PyModule 2: input data: ", np.shape(train_data))
        print("PyModule 2: output data: ", np.shape(output_data))
        #print(output_data)

        x_train, x_test, y_train, y_test = train_test_split(
            train_data, output_data, test_size=0.15)
        print("PyModule : Training array : ", x_train.shape)
        print("PyModule : Target array : ", y_train.shape)

        def custom_loss(y_true, y_pred):
            mu = y_pred[:, :1]  # first output neuron
            log_sig = y_pred[:, 1:]  # second output neuron
            sig = tf.exp(log_sig)  # undo the log

            return tf.reduce_mean(2*log_sig + ((y_true-mu)/sig)**2)

        print("PyModule : Training DNN model.")
        model = tf.keras.Sequential(
            [
                tf.keras.layers.Dense(
                    302, input_shape=(155,), activation='tanh'),
                tf.keras.layers.Dense(32, activation='tanh'),
                tf.keras.layers.Dense(21, activation='sigmoid')
            ])
        print(model.summary())

        model.compile(loss='mean_squared_error',
                      optimizer='adam', metrics=['accuracy'])
        model.fit(x_train, y_train, batch_size=128, epochs=150,
                  validation_data=(x_test, y_test), verbose=0)

        print("PyModule : Training done")
        print("PyModule: Test Score: [loss, accuracy]: ",
              model.evaluate(x_test, y_test, verbose=0))
        self.dnn = model

        return
    
    def upTrain(self, vec):
        """
        Desc.
        """
        
        self.train_table.line_of_data(vec)
        self.tracker = self.tracker + 1
        print("PyModule : Tracker: ", self.tracker)
        if self.tracker == 300:
            self.train()
            self.tracker = 0

        return
    
    def predict(self, trial_input):
        """
        Desc.
        """

        print("PyModule : starting predict")
        trial_input = np.asarray(trial_input)
        trial_input = trial_input.reshape(1, len(self.input_list))
        trial_input = self.SS1.transform(trial_input)
        # Predict gives the output and the std_dev in normal scale
        # no need to inverse transform anything
        re_predicted = []

        predicted = self.dnn.predict(trial_input, verbose=1)
        trial_input.__del__()
        predicted = self.SS2.inverse_transform(predicted)
    
        for i in range(0, len(self.output_list)):
            #print(self.output_list[i], predicted[0][i])
            re_predicted.append(predicted[0][i])
            re_predicted.append(0)

        predicted.__del__()
        return re_predicted
