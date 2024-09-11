/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2024, Andrew W. Steiner
  
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
#ifndef O2SCL_CLASSIFY_PYTHON_H
#define O2SCL_CLASSIFY_PYTHON_H

/** \file classify_python.h
    \brief File defining \ref o2scl::classify_python
*/

#include <iostream>
#include <string>
#include <cmath>

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>
#include <o2scl/tensor.h>
#include <o2scl/set_python.h>

#ifdef O2SCL_SET_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

namespace o2scl {

  /** \brief Multidimensional interpolation interface for python
   */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class vec_int_t=boost::numeric::ublas::vector<int>,
           class mat_x_t=o2scl::matrix_view_table<>,
           class mat_y_t=o2scl::matrix_view_table_transpose<> >
  class classify_python {
    
#if defined(O2SCL_SET_PYTHON) || defined(DOXYGEN)
    
  protected:
    
    /// Number of parameters
    size_t n_params;
    
    /// Number of outputs
    size_t n_outputs;
    
    /// Number of points
    size_t n_points;
    
  public:
    
    /// If true, throw exceptions on convergence errors (default true)
    bool err_nonconv;
    
    /** \brief Verbosity parameter (default 0)
     */
    int verbose;
    
  protected:

    /// Python module
    PyObject *p_module;
    
    /// Python class
    PyObject *p_class;

    /// An instance of the Python class
    PyObject *p_instance;

    /// Python function arguments for the set function
    PyObject *p_set_args;

    /// Python function arguments for the evaluation function
    PyObject *p_eval_args;

    /// Python set function
    PyObject *p_set_func;

    /// Python evaluation function
    PyObject *p_eval_func;

    /// Name of Python module 
    std::string c_module;

    /// Name of Python set function
    std::string c_set_func;

    /// Name of Python eval function
    std::string c_eval_func;

    /// Python class name
    std::string c_class_name;
    
    /// Python options
    std::string c_options;
    
  public:

    classify_python() {
      
      p_set_func=0;
      p_eval_func=0;
      p_set_args=0;
      p_eval_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      
      c_set_func="";
      c_eval_func="";
      c_class_name="";
      c_module="";
      c_options="";
      
      verbose=0;
      n_params=0;
      n_outputs=0;
      n_points=0;
      err_nonconv=true;
    }
    
    /** \brief Specify the Python module and function
     */
    classify_python(std::string class_name="", std::string options="",
                    int v=0, std::string module_name="o2sclpy",
                    std::string set_func="set_data_str",
                    std::string eval_func="eval") {
      
      verbose=0;
      n_params=0;
      n_outputs=0;
      n_points=0;
      err_nonconv=true;
      
      if (o2scl_settings.py_initialized==false) {
        if (this->verbose>0) {
          std::cout << "Running py_init()." << std::endl;
        }
        o2scl_settings.py_init();
      }
      
      p_set_func=0;
      p_eval_func=0;
      p_set_args=0;
      p_eval_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      
      set_functions(class_name,options,v,module_name,set_func,eval_func);
    }

    /** \brief Free memory associated with the Python objects and the
        interpolation

        \note This function does not finalize the Python interface.
    */
    void free() {
      if (this->verbose>0) {
        std::cout << "Starting classify_python::free()." << std::endl;
      }
      if (p_set_func!=0) {
        if (this->verbose>0) {
          std::cout << "Decref set_func." << std::endl;
        }
        Py_DECREF(p_set_func);
      }
      if (p_eval_func!=0) {
        if (this->verbose>0) {
          std::cout << "Decref eval_func." << std::endl;
        }
        Py_DECREF(p_eval_func);
      }
      if (p_set_args!=0) {
        if (this->verbose>0) {
          std::cout << "Decref set_args." << std::endl;
        }
        Py_DECREF(p_set_args);
      }
      if (p_eval_args!=0) {
        if (this->verbose>0) {
          std::cout << "Decref eval_args." << std::endl;
        }
        Py_DECREF(p_eval_args);
      }
      if (p_instance!=0) {
        if (this->verbose>0) {
          std::cout << "Decref instance." << std::endl;
        }
        Py_DECREF(p_instance);
      }
      if (p_class!=0) {
        if (this->verbose>0) {
          std::cout << "Decref class." << std::endl;
        }
        Py_DECREF(p_class);
      }
      if (p_module!=0) {
        if (this->verbose>0) {
          std::cout << "Decref module." << std::endl;
        }
        Py_DECREF(p_module);
      }
      
      p_set_func=0;
      p_eval_func=0;
      p_eval_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      
      this->n_params=0;
      this->n_outputs=0;
      this->n_points=0;
      
      if (this->verbose>0) {
        std::cout << "Done in classify_python::free()." << std::endl;
      }
    }      
    
    virtual ~classify_python() {
      free();
    }      
  
    /** \brief Specify the python module, class, functions, and options
    */
    int set_functions(std::string class_name="", std::string options="",
                      int v=0, std::string s_module="o2sclpy",
                      std::string set_func="set_data_str",
                      std::string eval_func="eval") {
      
      c_set_func=set_func;
      c_eval_func=eval_func;
      c_class_name=class_name;
      c_module=s_module;
      c_options=options;
      
      this->verbose=v;
      
      return 0;
    }
    
    /** \brief Set the data to be interpolated
     */
    virtual int set_data(size_t n_in, size_t n_out, size_t n_pts,
                         mat_x_t &user_x, mat_y_t &user_y) {

      tensor<> tin;
      tensor<int> tout;
      std::vector<size_t> in_size={n_pts,n_in}, out_size={n_pts,n_out};
      tin.resize(2,in_size);
      tout.resize(2,out_size);
      for(size_t j=0;j<n_pts;j++) {
        std::vector<size_t> ix;
        for(size_t k=0;k<n_in;k++) {
          ix={j,k};
          tin.get(ix)=user_x(j,k);
        }
        for(size_t k=0;k<n_out;k++) {
          ix={j,k};
          tout.get(ix)=user_y(j,k);
        }
      }

      return set_data_tensor(n_in,n_out,n_pts,tin,tout);
    }

    /** \brief Set the data to be interpolated (tensor form)
     */
    int set_data_tensor(size_t n_in, size_t n_out, size_t n_pts,
                        const o2scl::tensor<> &params,
                        const o2scl::tensor<int> &outputs) {
      int ret;
      set_data_internal(n_in,n_out,n_pts,params,outputs,ret);
      return ret;
    }

  protected:
    
    /** \brief The internal set data function
        
        \note This function returns a void pointer to be
        compatible with the Python import_array() macro. 
     */
    void *set_data_internal(size_t n_pars, size_t n_out, size_t n_dat, 
                            const o2scl::tensor<> &params,
                            const o2scl::tensor<int> &outputs,
                            int &ret) {
      
      ret=0;
      
      if (params.get_rank()!=2 || outputs.get_rank()!=2) {
        O2SCL_ERR2("Invalid rank for input tensors in ",
                   "classify_python::set_function().",o2scl::exc_einval);
      }
      
      free();
      
      this->n_params=n_pars;
      this->n_points=n_dat;
      this->n_outputs=n_out;
      
      p_module=o2scl_settings.py_import_module(c_module,this->verbose);
      
      if (c_class_name.length()>0) {
        if (this->verbose>0) {
          std::cout << "  Obtaining python class name "
                    << c_class_name << "." << std::endl;
        }
        p_class=PyObject_GetAttrString(p_module,c_class_name.c_str());
        if (p_class==0) {
          O2SCL_ERR2("Get class failed in ",
                     "classify_python::set_function().",o2scl::exc_efailed);
        }
        
        // Create an instance of the class
        if (this->verbose>0) {
          std::cout << "  Loading python class." << std::endl;
        }
        if (PyCallable_Check(p_class)==false) {
          O2SCL_ERR2("Check class callable failed in ",
                     "funct_python_method::set_function().",
                     o2scl::exc_efailed);
        }
        
        if (this->verbose>0) {
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
      if (this->verbose>0) {
        std::cout << "  Making argument object for set function."
                  << std::endl;
      }
      p_set_args=PyTuple_New(3);
      if (p_set_args==0) {
        O2SCL_ERR2("Create arg tuple failed in ",
                   "classify_python::set_function().",
                   o2scl::exc_efailed);
      }

      // Setup the arguments to the python function
      if (this->verbose>0) {
        std::cout << "  Making argument object for eval function."
                  << std::endl;
      }
      p_eval_args=PyTuple_New(1);
      if (p_eval_args==0) {
        O2SCL_ERR2("Create arg tuple failed in ",
                   "classify_python::set_function().",
                   o2scl::exc_efailed);
      }

      if (c_class_name.length()>0) {

        // Load the python function
        if (this->verbose>0) {
          std::cout << "  Loading python member function eval: "
                    << c_eval_func << std::endl;
        }
        p_eval_func=PyObject_GetAttrString(p_instance,c_eval_func.c_str());
        if (p_eval_func==0) {
          O2SCL_ERR2("Get eval function failed in ",
                     "classify_python::set_function().",
                     o2scl::exc_efailed);
        }
        
        // Load the python function
        if (this->verbose>0) {
          std::cout << "  Loading python member function set: "
                    << c_set_func << std::endl;
        }
        p_set_func=PyObject_GetAttrString(p_instance,c_set_func.c_str());
        if (p_set_func==0) {
          O2SCL_ERR2("Get set function failed in ",
                     "classify_python::set_function().",
                     o2scl::exc_efailed);
        }

      } else {
    
        // Load the python function
        if (this->verbose>0) {
          std::cout << "  Loading python function set." << std::endl;
        }
        p_set_func=PyObject_GetAttrString(p_module,c_set_func.c_str());
        if (p_set_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "classify_python::set_function().",
                     o2scl::exc_efailed);
        }

        // Load the python function
        if (this->verbose>0) {
          std::cout << "  Loading python function eval." << std::endl;
        }
        p_eval_func=PyObject_GetAttrString(p_module,c_eval_func.c_str());
        if (p_eval_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "classify_python::set_function().",
                     o2scl::exc_efailed);
        }

      }

      // AWS, 2/21/23: I'm not sure why it has to be done here and not in
      // a different function, but if I don't do it here I get a seg fault.
      //void *vp=o2scl_settings.py_import_array();
      import_array();

      if (params.get_size(0)!=this->n_points) {
        O2SCL_ERR("Input data does not have correct number of rows.",
                  o2scl::exc_einval);
      }
      if (params.get_size(1)!=this->n_params) {
        O2SCL_ERR("Input data does not have correct number of columns.",
                  o2scl::exc_einval);
      }
      if (outputs.get_size(0)!=this->n_points) {
        O2SCL_ERR("Output data does not have correct number of rows.",
                  o2scl::exc_einval);
      }
      if (outputs.get_size(1)!=this->n_outputs) {
        O2SCL_ERR("Output data does not have correct number of columns.",
                  o2scl::exc_einval);
      }
      
      npy_intp params_dims[]={(npy_intp)params.get_size(0),
                              (npy_intp)params.get_size(1)};
      if (this->verbose>1) {
        std::cout << "classify_python::operator():" << std::endl;
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
      if (this->verbose>1) {
        std::cout << "classify_python::operator():" << std::endl;
      }
      PyObject *array_out=PyArray_SimpleNewFromData
        (2,outputs_dims,NPY_INT,(void *)(&(outputs.get_data()[0])));
      
      int ret2=PyTuple_SetItem(p_set_args,1,array_out);
      if (ret2!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "mm_funct_python::operator().",o2scl::exc_efailed);
      }

      if (this->verbose>0) {
        std::cout << "Creating python unicode for string: "
                  << c_options.length() << " " << c_options << std::endl;
      }
      PyObject *p_options=PyUnicode_FromString(c_options.c_str());
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
      if (this->verbose>0) {
        std::cout << "  Calling python set function." << std::endl;
      }
      PyObject *result=PyObject_CallObject(p_set_func,p_set_args);
      if (result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "classify_python::operator().",o2scl::exc_efailed);
      }

      if (this->verbose>0) {
        std::cout << "Done with classify_python::set_function()."
                  << std::endl;
      }

      return 0;
    }      

  public:
    
    /** \brief Compute the function at point \c x and return the result
     */
    int eval_std_vec(const std::vector<double> &x,
                     std::vector<int> &y) const {
      
      if (x.size()!=this->n_params) {
        O2SCL_ERR("Input vector does not have correct size.",
                  o2scl::exc_einval);
      }
      if (y.size()!=this->n_outputs) {
        O2SCL_ERR("Output vector does not have correct size.",
                  o2scl::exc_einval);
      }
  
      if (p_set_func==0 || p_eval_func==0) {
        O2SCL_ERR2("No functions found in ",
                   "classify_python::operator().",
                   o2scl::exc_efailed);
      }

      npy_intp x_dims[]={(npy_intp)x.size()};
      if (this->verbose>1) {
        std::cout << "classify_python::operator():" << std::endl;
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
      if (this->verbose>1) {
        std::cout << "  Calling python eval function." << std::endl;
      }
      PyObject *result=PyObject_CallObject(p_eval_func,p_eval_args);
      if (result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "classify_python::operator().",o2scl::exc_efailed);
      }

      if (PyArray_Check(result)==0) {
        O2SCL_ERR2("Function call did not return a numpy array in ",
                   "classify_python::operator().",o2scl::exc_efailed);
      }
      
      if (this->verbose>1) {
        std::cout << "  Obtaining output." << std::endl;
      }
      void *vp=PyArray_DATA((PyArrayObject *)result);
      int *dp=(int *)vp;
      for(size_t i=0;i<this->n_outputs;i++) {
        y[i]=dp[i];
        if (this->verbose>1) {
          std::cout << "  i,y[i]: " << i << " " << y[i] << std::endl;
        }
      }
      
      if (this->verbose>1) {
        std::cout << "  Decref result." << std::endl;
      }
      Py_DECREF(result);
  
      if (this->verbose>1) {
        std::cout << "Done in classify_python::operator()."
                  << std::endl;
      }

      return 0;
    }

    /** \brief Evaluate the interpolation at point \c x,
        returning \c y
    */
    virtual int eval(const vec_t &x, vec_int_t &y) const {
      std::vector<double> x2(this->n_params);
      std::vector<int> y2(this->n_outputs);
      vector_copy(this->n_params,x,x2);
      int ret=eval_std_vec(x2,y2);
      vector_copy(this->n_outputs,y2,y);
      return ret;
    }
    
  private:

    classify_python(const classify_python &);
    classify_python& operator=(const classify_python&);

#endif

  };
  
}
    
#endif


