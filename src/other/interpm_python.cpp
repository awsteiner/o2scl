/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2022-2023, Mahamudul Hasan Anik, Satyajit Roy, and
  Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/interpm_python.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_PYTHON

interpm_python::interpm_python() {
  p_set_func=0;
  p_eval_func=0;
  p_set_args=0;
  p_eval_args=0;
  p_instance=0;
  p_class=0;
  p_module=0;
  p_name=0;
      
  n_params=0;
  n_outputs=0;
  n_points=0;
}
    
/** \brief Specify the Python module and function
 */
interpm_python::interpm_python(std::string module, std::string set_func,
                               std::string eval_func,
                               std::string eval_unc_func,
                               size_t n_pars, size_t n_dat, size_t n_out,
                               const o2scl::tensor<> &params,
                               const o2scl::tensor<> &outputs,
                               std::string options, 
                               std::string class_name, int v) {
                    
  verbose=v;

  if (o2scl_settings.py_initialized==false) {
    if (verbose>0) {
      std::cout << "Running py_init()." << std::endl;
    }
    o2scl_settings.py_init();
  }
  p_set_func=0;
  p_eval_func=0;
  p_eval_unc_func=0;
  p_set_args=0;
  p_eval_args=0;
  p_instance=0;
  p_class=0;
  p_module=0;
  p_name=0;
      
  n_params=0;
  n_outputs=0;
  n_points=0;
      
  if (module.length()>0) {
    set_function(module,set_func,eval_func,eval_unc_func,
                 n_pars,n_dat,n_out,
                 params,outputs,options,class_name,v);
  }
}      
    
void interpm_python::free() {
  if (verbose>0) {
    std::cout << "Starting interpm_python::free()." << std::endl;
  }
  if (p_set_func!=0) {
    if (verbose>0) {
      std::cout << "Decref set_func." << std::endl;
    }
    Py_DECREF(p_set_func);
  }
  if (p_eval_func!=0) {
    if (verbose>0) {
      std::cout << "Decref eval_func." << std::endl;
    }
    Py_DECREF(p_eval_func);
  }
  if (p_eval_unc_func!=0) {
    if (verbose>0) {
      std::cout << "Decref eval_unc_func." << std::endl;
    }
    Py_DECREF(p_eval_unc_func);
  }
  if (p_set_args!=0) {
    if (verbose>0) {
      std::cout << "Decref set_args." << std::endl;
    }
    Py_DECREF(p_set_args);
  }
  if (p_eval_args!=0) {
    if (verbose>0) {
      std::cout << "Decref eval_args." << std::endl;
    }
    Py_DECREF(p_eval_args);
  }
  if (p_instance!=0) {
    if (verbose>0) {
      std::cout << "Decref instance." << std::endl;
    }
    Py_DECREF(p_instance);
  }
  if (p_class!=0) {
    if (verbose>0) {
      std::cout << "Decref class." << std::endl;
    }
    Py_DECREF(p_class);
  }
  if (p_module!=0) {
    if (verbose>0) {
      std::cout << "Decref module." << std::endl;
    }
    Py_DECREF(p_module);
  }
  if (p_name!=0) {
    if (verbose>0) {
      std::cout << "Decref name." << std::endl;
    }
    Py_DECREF(p_name);
  }

  p_set_func=0;
  p_eval_func=0;
  p_eval_unc_func=0;
  p_eval_args=0;
  p_instance=0;
  p_class=0;
  p_module=0;
  p_name=0;
      
  n_params=0;
  n_outputs=0;
  n_points=0;
      
  if (verbose>0) {
    std::cout << "Done in interpm_python::free()." << std::endl;
  }
}      
    
interpm_python::~interpm_python() {
  free();
}      
  
int interpm_python::set_function(std::string module, std::string set_func,
                                 std::string eval_func,
                                 std::string eval_unc_func,
                                 size_t n_pars, size_t n_dat, size_t n_out,
                                 const o2scl::tensor<> &params,
                                 const o2scl::tensor<> &outputs,
                                 std::string options,
                                 std::string class_name, int v) {
  int ret;
  void *vp=set_function_internal(module,set_func,eval_func,eval_unc_func,
                                 n_pars,n_dat,n_out,params,outputs,options,
                                 class_name,v,ret);
  return ret;
}
      
void *interpm_python::set_function_internal
(std::string module, std::string set_func,
 std::string eval_func,
 std::string eval_unc_func,
 size_t n_pars, size_t n_dat, size_t n_out,
 const o2scl::tensor<> &params,
 const o2scl::tensor<> &outputs,
 std::string options,
 std::string class_name, int v, int &ret) {

  ret=0;
  
  if (params.get_rank()!=2 || outputs.get_rank()!=2) {
    O2SCL_ERR2("Invalid rank for input tensors in ",
               "interpm_python::set_function().",o2scl::exc_einval);
  }
      
  free();

  n_params=n_pars;
  n_points=n_dat;
  n_outputs=n_out;
      
  // Get the Unicode name of the user-specified module
  if (verbose>0) {
    std::cout << "Python version: "
              << o2scl_settings.py_version() << std::endl;
    std::cout << "Staring interpm_python::set_function()."
              << std::endl;
    std::cout << "  Getting unicode for module named "
              << module << std::endl;
  }
  p_name=PyUnicode_FromString(module.c_str());
  if (p_name==0) {
    O2SCL_ERR2("Create module name failed in ",
               "interpm_python::set_function().",
               o2scl::exc_efailed);
  }
      
  // Import the user-specified module
  if (verbose>0) {
    std::cout << "  Importing module." << std::endl;
  }
  p_module=PyImport_Import(p_name);
  if (p_module==0) {
    O2SCL_ERR2("Load module failed in ",
               "interpm_python::set_function().",
               o2scl::exc_efailed);
  }

  if (class_name.length()>0) {
    if (verbose>0) {
      std::cout << "  Obtaining python class." << std::endl;
    }
    p_class=PyObject_GetAttrString(p_module,class_name.c_str());
    if (p_class==0) {
      O2SCL_ERR2("Get class failed in ",
                 "interpm_python::set_function().",o2scl::exc_efailed);
    }
        
    // Create an instance of the class
    if (verbose>0) {
      std::cout << "  Loading python class." << std::endl;
    }
    if (PyCallable_Check(p_class)==false) {
      O2SCL_ERR2("Check class callable failed in ",
                 "funct_python_method::set_function().",
                 o2scl::exc_efailed);
    }
        
    if (verbose>0) {
      std::cout << "  Loading python class instance." << std::endl;
    }
    p_instance=PyObject_CallObject(p_class,0);
    if (p_instance==0) {
      O2SCL_ERR2("Instantiate class failed in ",
                 "funct_python_method::set_function().",
                 o2scl::exc_efailed);
    }
  }
      
  // Setup the arguments to the python function
  if (verbose>0) {
    std::cout << "  Making argument object for set function."
              << std::endl;
  }
  p_set_args=PyTuple_New(3);
  if (p_set_args==0) {
    O2SCL_ERR2("Create arg tuple failed in ",
               "interpm_python::set_function().",
               o2scl::exc_efailed);
  }

  // Setup the arguments to the python function
  if (verbose>0) {
    std::cout << "  Making argument object for eval function."
              << std::endl;
  }
  p_eval_args=PyTuple_New(1);
  if (p_eval_args==0) {
    O2SCL_ERR2("Create arg tuple failed in ",
               "interpm_python::set_function().",
               o2scl::exc_efailed);
  }

  if (class_name.length()>0) {

    // Load the python function
    if (verbose>0) {
      std::cout << "  Loading python member function eval: "
                << eval_func<< std::endl;
    }
    p_eval_func=PyObject_GetAttrString(p_instance,eval_func.c_str());
    if (p_eval_func==0) {
      O2SCL_ERR2("Get eval function failed in ",
                 "interpm_python::set_function().",
                 o2scl::exc_efailed);
    }
        
    // Load the python function
    if (verbose>0) {
      std::cout << "  Loading python member function eval_unc: "
                << eval_unc_func << std::endl;
    }
    p_eval_unc_func=PyObject_GetAttrString(p_instance,eval_unc_func.c_str());
    if (p_eval_unc_func==0) {
      O2SCL_ERR2("Get eval_unc function failed in ",
                 "interpm_python::set_function().",
                 o2scl::exc_efailed);
    }
        
    // Load the python function
    if (verbose>0) {
      std::cout << "  Loading python member function set: "
                << set_func << std::endl;
    }
    p_set_func=PyObject_GetAttrString(p_instance,set_func.c_str());
    if (p_set_func==0) {
      O2SCL_ERR2("Get set function failed in ",
                 "interpm_python::set_function().",
                 o2scl::exc_efailed);
    }

  } else {
    
    // Load the python function
    if (verbose>0) {
      std::cout << "  Loading python function set." << std::endl;
    }
    p_set_func=PyObject_GetAttrString(p_module,set_func.c_str());
    if (p_set_func==0) {
      O2SCL_ERR2("Get function failed in ",
                 "interpm_python::set_function().",
                 o2scl::exc_efailed);
    }

    // Load the python function
    if (verbose>0) {
      std::cout << "  Loading python function eval." << std::endl;
    }
    p_eval_func=PyObject_GetAttrString(p_module,eval_func.c_str());
    if (p_eval_func==0) {
      O2SCL_ERR2("Get function failed in ",
                 "interpm_python::set_function().",
                 o2scl::exc_efailed);
    }

    // Load the python function
    if (verbose>0) {
      std::cout << "  Loading python function eval_unc." << std::endl;
    }
    p_eval_unc_func=PyObject_GetAttrString(p_module,eval_unc_func.c_str());
    if (p_eval_unc_func==0) {
      O2SCL_ERR2("Get function failed in ",
                 "interpm_python::set_function().",
                 o2scl::exc_efailed);
    }
  }

  // AWS, 2/21/23: I'm not sure why it has to be done here and not in
  // a different function, but if I don't do it here I get a seg fault.
  //void *vp=o2scl_settings.py_import_array();
  import_array();

  if (params.get_size(0)!=n_points) {
    O2SCL_ERR("Input data does not have correct number of rows.",
              o2scl::exc_einval);
  }
  if (params.get_size(1)!=n_params) {
    O2SCL_ERR("Input data does not have correct number of columns.",
              o2scl::exc_einval);
  }
  if (outputs.get_size(0)!=n_points) {
    O2SCL_ERR("Output data does not have correct number of rows.",
              o2scl::exc_einval);
  }
  if (outputs.get_size(1)!=n_outputs) {
    O2SCL_ERR("Output data does not have correct number of columns.",
              o2scl::exc_einval);
  }
      
  npy_intp params_dims[]={(npy_intp)params.get_size(0),
    (npy_intp)params.get_size(1)};
  if (verbose>1) {
    std::cout << "interpm_python::operator():" << std::endl;
  }
  PyObject *array_in=PyArray_SimpleNewFromData
    (2,params_dims,NPY_DOUBLE,(void *)(&(params.get_data()[0])));
         
  int pret=PyTuple_SetItem(p_set_args,0,array_in);
  if (pret!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "mm_funct_python::operator().",o2scl::exc_efailed);
  }
      
  npy_intp outputs_dims[]={(npy_intp)outputs.get_size(0),
    (npy_intp)outputs.get_size(1)};
  if (verbose>1) {
    std::cout << "interpm_python::operator():" << std::endl;
  }
  PyObject *array_out=PyArray_SimpleNewFromData
    (2,outputs_dims,NPY_DOUBLE,(void *)(&(outputs.get_data()[0])));
      
  int ret2=PyTuple_SetItem(p_set_args,1,array_out);
  if (ret2!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "mm_funct_python::operator().",o2scl::exc_efailed);
  }

  if (verbose>0) {
    std::cout << "Creating python unicode for string: "
              << options.length() << " " << options << std::endl;
  }
  PyObject *p_options=PyUnicode_FromString(options.c_str());
  if (p_options==0) {
    O2SCL_ERR2("String creation failed in ",
               "emulator_python::set().",o2scl::exc_efailed);
  }
      
  int ret3=PyTuple_SetItem(p_set_args,2,p_options);
  if (ret3!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "mm_funct_python::operator().",o2scl::exc_efailed);
  }

  // Call the python function
  if (verbose>0) {
    std::cout << "  Calling python set function." << std::endl;
  }
  PyObject *result=PyObject_CallObject(p_set_func,p_set_args);
  if (result==0) {
    O2SCL_ERR2("Function call failed in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }

  if (verbose>0) {
    std::cout << "Done with interpm_python::set_function()."
              << std::endl;
  }

  return 0;
}
    
int interpm_python::eval(const std::vector<double> &x,
                         std::vector<double> &y) const {

  if (x.size()!=n_params) {
    O2SCL_ERR("Input vector does not have correct size.",
              o2scl::exc_einval);
  }
  if (y.size()!=n_outputs) {
    O2SCL_ERR("Output vector does not have correct size.",
              o2scl::exc_einval);
  }
  
  if (p_set_func==0 || p_eval_func==0 || p_eval_unc_func==0) {
    O2SCL_ERR2("No functions found in ",
               "interpm_python::operator().",
               o2scl::exc_efailed);
  }

  npy_intp x_dims[]={(npy_intp)x.size()};
  if (verbose>1) {
    std::cout << "interpm_python::operator():" << std::endl;
    std::cout << "  Array x: " << x.size() << std::endl;
  }
  PyObject *array_x=PyArray_SimpleNewFromData
    (1,x_dims,NPY_DOUBLE,(void *)(&(x[0])));
      
  int ret=PyTuple_SetItem(p_eval_args,0,array_x);
  if (ret!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "mm_funct_python::operator().",o2scl::exc_efailed);
  }
      
  // Call the python function
  if (verbose>1) {
    std::cout << "  Calling python eval function." << std::endl;
  }
  PyObject *result=PyObject_CallObject(p_eval_func,p_eval_args);
  if (result==0) {
    O2SCL_ERR2("Function call failed in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }

  if (PyArray_Check(result)==0) {
    O2SCL_ERR2("Function call did not return a numpy array in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }
      
  if (verbose>1) {
    std::cout << "  Obtaining output." << std::endl;
  }
  void *vp=PyArray_DATA((PyArrayObject *)result);
  double *dp=(double *)vp;
  for(size_t i=0;i<n_outputs;i++) {
    y[i]=dp[i];
    if (verbose>1) {
      std::cout << "  i,y[i]: " << i << " " << y[i] << std::endl;
    }
  }
      
  if (verbose>1) {
    std::cout << "  Decref result." << std::endl;
  }
  Py_DECREF(result);
  
  if (verbose>1) {
    std::cout << "Done in interpm_python::operator()."
              << std::endl;
  }

  return 0;
}      

int interpm_python::eval_unc(const std::vector<double> &x,
                             std::vector<double> &y,
                             std::vector<double> &yp) const {

  if (x.size()!=n_params) {
    O2SCL_ERR("Input vector does not have correct size.",
              o2scl::exc_einval);
  }
  if (y.size()!=n_outputs) {
    O2SCL_ERR("Output vector does not have correct size.",
              o2scl::exc_einval);
  }
  
  if (p_set_func==0 || p_eval_unc_func==0) {
    O2SCL_ERR2("No functions found in ",
               "interpm_python::operator().",
               o2scl::exc_efailed);
  }

  npy_intp x_dims[]={(npy_intp)x.size()};
  if (verbose>1) {
    std::cout << "interpm_python::operator():" << std::endl;
    std::cout << "  Array x: " << x.size() << std::endl;
  }
  PyObject *array_x=PyArray_SimpleNewFromData
    (1,x_dims,NPY_DOUBLE,(void *)(&(x[0])));
      
  int ret=PyTuple_SetItem(p_eval_args,0,array_x);
  if (ret!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "mm_funct_python::operator().",o2scl::exc_efailed);
  }
      
  // Call the python function
  if (verbose>1) {
    std::cout << "  Calling python eval_unc function." << std::endl;
  }
  PyObject *result=PyObject_CallObject(p_eval_unc_func,p_eval_args);
  if (result==0) {
    O2SCL_ERR2("Function call failed in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }

  if (PyTuple_Check(result)==0) {
    O2SCL_ERR2("Function call did not return a Python tuple in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }

  PyObject *p_y, *p_yp;
  p_y=PyTuple_GetItem(result,0);
  p_yp=PyTuple_GetItem(result,1);

  if (PyArray_Check(p_y)==0) {
    O2SCL_ERR2("Python tuple does not contain array in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }
  if (PyArray_Check(p_yp)==0) {
    O2SCL_ERR2("Python tuple does not contain array in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }
  
  if (verbose>1) {
    std::cout << "  Obtaining output." << std::endl;
  }
  
  void *vp=PyArray_DATA((PyArrayObject *)p_y);
  double *dp=(double *)vp;
  for(size_t i=0;i<n_outputs;i++) {
    y[i]=dp[i];
    if (verbose>1) {
      std::cout << "  i,y[i]: " << i << " " << y[i] << std::endl;
    }
  }
      
  void *vp2=PyArray_DATA((PyArrayObject *)p_yp);
  double *dp2=(double *)vp2;
  for(size_t i=0;i<n_outputs;i++) {
    yp[i]=dp2[i];
    if (verbose>1) {
      std::cout << "  i,yp[i]: " << i << " " << yp[i] << std::endl;
    }
  }
      
  if (verbose>1) {
    std::cout << "  Decref result." << std::endl;
  }
  Py_DECREF(result);
  
  if (verbose>1) {
    std::cout << "Done in interpm_python::operator()."
              << std::endl;
  }

  return 0;
}      

#endif
