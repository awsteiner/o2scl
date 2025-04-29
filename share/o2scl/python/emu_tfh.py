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
        self.data=0
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
        
        self.data=o2sclpy.table(link)

        col_list_arr=col_list.split(',')

        for i in range(0, nd_in):
            self.input_list.append(col_list_arr[i])

        for i in range(nd_in, len(col_list_arr)):
           self.output_list.append(col_list_arr[i])
        nd_out=len(self.output_list)

        print("emu_dnn::read_data(): Reading file:",hdf_file)
        hf=o2sclpy.hdf_file(link)
        hf.open(hdf_file)
        o2sclpy.hdf_input_table(link,hf,self.data)
        hf.close()
        print("emu_dnn::read_data(): Done reading file.")
        print("emu_dnn::read_data(): Table in file", hdf_file,"has",
              self.data.get_nlines(),"lines.")

        print("emu_dnn::read_data(): table now has",
              self.data.get_nlines(),
              "lines.")

        train_data=[]
        output_data=[]

        for i in self.input_list:
            train_data.append(np.array(self.data[i]
                                       [0:self.data.get_nlines()]))
        train_data=self.SS1.fit_transform(np.asarray(train_data).transpose())
         
        for i in self.output_list:
            output_data.append(np.array(self.data[i]
                                        [0:self.data.get_nlines()]))
        output_data=self.SS2.fit_transform(np.asarray(output_data).transpose())
        
        print("emu_dnn::read_data(): Input data:",np.shape(train_data))
        print("emu_dnn::read_data(): Output data:",np.shape(output_data))

        x_train,x_test,y_train,y_test=train_test_split(
            train_data,output_data,test_size=0.15)
        print("emu_dnn::read_data(): Training array:",x_train.shape)
        print("emu_dnn::read_data(): Target array:",y_train.shape)

        def custom_loss(y_true, y_pred):
            # first output neuron
            mu=y_pred[:,:1]
            # second output neuron
            log_sig=y_pred[:,1:]
            # undo the log
            sig=tf.exp(log_sig)  

            return tf.reduce_mean(2*log_sig+((y_true-mu)/sig)**2)

        class RegressionHyperModel(HyperModel):
            def __init__(self, input_shape):
                self.input_shape = input_shape    
        
            def build(self, hp):
                model = keras.Sequential()
                # Tune the number of layers
                for i in range(hp.Int('num_layers', 1, 6)):
                    model.add(
                        layers.Dense(
                            # Tune number of units separately
                            units=hp.Int(f"units_{i}",
                                         min_value=50,
                                         max_value=200, step=50),
                            activation='relu'
                        )
                    )
                model.add(layers.Dense(1, activation='linear'))
                model.compile(loss='mean_squared_error', optimizer='adam')
                return model

        # Initialize the input shape
        input_shape = (x_tr.shape[1],)
        hypermodel = RegressionHyperModel(input_shape)

        tuner = kt.RandomSearch(
            # Pass the hypermodel object
            hypermodel,
            # Quantity to monitor during tuning
            objective='val_loss',
            # Set reproducibility of randomness
            seed=42,
            # Max number of trials with different hyperparameters
            max_trials=100,
            # Number of repeated trials with same hyperparameters
            executions_per_trial=1,
            # Set directory to store search results
            directory="random_search",
            # Set the subdirectory name
            project_name="np",
            # Choose if previous search results should be ignored
            overwrite=True             
        )
        
        # Set up callback for early stopping 
        stop_early=tf.keras.callbacks.EarlyStopping(monitor='loss',
                                                    min_delta=1.0e-10,
                                                    patience=10)
        
        # Print the summary of search space
        tuner.search_space_summary()

        tuner.search(x_tr, y_tr, batch_size=32, epochs=1000,
                     validation_data=(x_ts, y_ts),
                     callbacks=[stop_early], verbose=2)

        tuner.results_summary(num_trials=10)

        # Get the top model
        models = tuner.get_best_models(num_models=10)
        best_model = models[0]
        
        # Build the best model
        best_model.build(input_shape=(None, 91))
        
        # Show the best model
        best_model.summary()
        
        print("emu_dnn::read_data(): Training DNN model.")
        model=tf.keras.Sequential(
            [
                tf.keras.layers.Dense(
                    nd_in,input_shape=(nd_in,),activation='tanh'),
                tf.keras.layers.Dense(32,activation='tanh'),
                tf.keras.layers.Dense(32,activation='tanh'),
                tf.keras.layers.Dense(nd_out,activation='sigmoid')
            ])
        print(model.summary())

        model.compile(loss='mean_squared_error',
                      optimizer='adam',metrics=['accuracy'])
        model.fit(x_train,y_train,batch_size=128,epochs=150,
                  validation_data=(x_test,y_test),verbose=0)

        print("emu_dnn::read_data(): Training done.")
        print("emu_dnn::read_data(): Test Score: [loss, accuracy]: ",
              model.evaluate(x_test, y_test, verbose=0))
        self.dnn=model

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
        predicted=self.SS2.inverse_transform(predicted)
    
        return predicted

if __name__ == '__main__':
    ednn=emu_dnn();
    ednn.read_data(2,'emu_data.o2',0,'x,y,z,d',1)
    print(ednn.predict([1,2]))
    print(ednn.predict([-1.519097,-4.760777]))
    
