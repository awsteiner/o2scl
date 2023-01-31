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
    
    
