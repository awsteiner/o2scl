/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2023, Andrew W. Steiner
  
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
#ifndef O2SCL_KDE_PYTHON_H
#define O2SCL_KDE_PYTHON_H

/** \file kde_python.h
    \brief File defining \ref o2scl::kde_python
*/

#include <iostream>
#include <string>
#include <cmath>

#include <o2scl/err_hnd.h>
#include <o2scl/tensor.h>
#include <o2scl/exp_max.h>

#ifdef O2SCL_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

namespace o2scl {

#ifdef O2SCL_PYTHON

  /** \brief Multidimensional interpolation interface for python
   */
  template<class vec_t=std::vector<double> >
  class kde_python : public prob_dens_mdim<vec_t> {
    
  protected:

    /// Python unicode object containing function name
    PyObject *p_name;
    
    /// Python module containing function
    PyObject *p_module;
    
    /// The class
    PyObject *p_class;

    /// An instance of the class
    PyObject *p_instance;

    /// Function arguments
    PyObject *p_set_args;

    /// Function arguments
    PyObject *p_ld_args;

    /// Python function
    PyObject *p_set_func;

    /// Python function
    PyObject *p_sample_func;

    /// Python function
    PyObject *p_ld_func;

    /// Number of parameters
    size_t n_params;
    
    /// Number of points
    size_t n_points;

    /// Data
    o2scl::tensor<> data;
    
  public:

    /** \brief If true, then make a copy of the array before 
        calling log_pdf() (default true)

        \todo This should be done by a template specialization instead
        of a bool flag
    */
    bool array_copy;
    
    kde_python()  {
      
      if (o2scl_settings.py_initialized==false) {
        if (this->verbose>1) {
          std::cout << "Running py_init()." << std::endl;
        }
        o2scl_settings.py_init();
      }
      
      p_ld_func=0;
      p_sample_func=0;
      p_set_func=0;
      p_ld_args=0;
      p_set_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      p_name=0;
      
      n_params=0;
      n_points=0;
      array_copy=true;
    }
    
    /** \brief Specify the Python module and function
     */
    kde_python(std::string module, std::string set_func,
               std::string sample_func, std::string ld_func,
               size_t n_pars, size_t n_dat,
               o2scl::tensor<> &params,
               std::vector<double> array,
               std::string options="", 
               std::string class_name="", int v=0)  {
                    
      this->verbose=v;
      
      if (o2scl_settings.py_initialized==false) {
        if (this->verbose>1) {
          std::cout << "Running py_init()." << std::endl;
        }
        o2scl_settings.py_init();
      }
      p_ld_func=0;
      p_sample_func=0;
      p_set_func=0;
      p_ld_args=0;
      p_set_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      p_name=0;
      
      n_params=0;
      n_points=0;

      if (module.length()>0) {
        set_function(module,set_func,sample_func,ld_func,
                     n_pars,n_dat,params,array,options,class_name,v);
      }
    }      
    
    /** \brief Free memory associated with the KDE
     */
    virtual ~kde_python() {
      free();
    }      
    
    /** \brief Free the associated memory
     */
    void free() {
      if (this->verbose>1) {
        std::cout << "Starting kde_python::free()." << std::endl;
      }
      if (p_ld_func!=0) {
        if (this->verbose>1) {
          std::cout << "Decref ld_func." << std::endl;
        }
        Py_DECREF(p_ld_func);
      }
      if (p_sample_func!=0) {
        if (this->verbose>1) {
          std::cout << "Decref sample_func." << std::endl;
        }
        Py_DECREF(p_sample_func);
      }
      if (p_set_func!=0) {
        if (this->verbose>1) {
          std::cout << "Decref set_func." << std::endl;
        }
        Py_DECREF(p_set_func);
      }
      if (p_ld_args!=0) {
        if (this->verbose>1) {
          std::cout << "Decref ld_args." << std::endl;
        }
        Py_DECREF(p_ld_args);
      }
      if (p_set_args!=0) {
        if (this->verbose>1) {
          std::cout << "Decref set_args." << std::endl;
        }
        Py_DECREF(p_set_args);
      }
      if (p_instance!=0) {
        if (this->verbose>1) {
          std::cout << "Decref instance." << std::endl;
        }
        Py_DECREF(p_instance);
      }
      if (p_class!=0) {
        if (this->verbose>1) {
          std::cout << "Decref class." << std::endl;
        }
        Py_DECREF(p_class);
      }
      if (p_module!=0) {
        if (this->verbose>1) {
          std::cout << "Decref module." << std::endl;
        }
        Py_DECREF(p_module);
      }
      if (p_name!=0) {
        if (this->verbose>1) {
          std::cout << "Decref name." << std::endl;
        }
        Py_DECREF(p_name);
      }
      
      p_ld_func=0;
      p_sample_func=0;
      p_set_func=0;
      p_ld_args=0;
      p_set_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      p_name=0;
      
      n_params=0;
      n_points=0;
      
      if (this->verbose>1) {
        std::cout << "Done in kde_python::free()." << std::endl;
      }
    }      
    
    /** \brief Get the data
     */
    const o2scl::tensor<> &get_data() {
      return data;
    }
    
    /** \brief Specify the python and the parameters
        
        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string set_func,
                     std::string sample_func, std::string ld_func,
                     size_t n_pars, size_t n_dat, 
                     o2scl::tensor<> &params,
                     std::vector<double> array,
                     std::string options="",
                     std::string class_name="", int v=0) {
      int ret;
      void *vp=set_function_internal(module,set_func,sample_func,
                                     ld_func,n_pars,n_dat,params,
                                     array,ret,options,class_name,v);
      return ret;
    }

    /** \brief Internal version of set_function()
     */
    void *set_function_internal
      (std::string module, std::string set_func,
       std::string sample_func, std::string ld_func,
       size_t n_pars, size_t n_dat, 
       o2scl::tensor<> &params,
       std::vector<double> array, int &ret,
       std::string options="",
       std::string class_name="", int v=0) {

      ret=0;
      
      this->verbose=v;

      free();
      
      if (params.get_rank()!=2) {
        O2SCL_ERR2("Invalid rank for input tensors in ",
                   "kde_python().",o2scl::exc_einval);
      }
      
      n_params=n_pars;
      n_points=n_dat;
      if (n_pars==0) {
        O2SCL_ERR2("Invalid number of parameters in ",
                   "kde_python().",o2scl::exc_einval);
      }
      if (n_dat==0) {
        O2SCL_ERR2("Invalid number of data points in ",
                   "kde_python().",o2scl::exc_einval);
      }

      // AWS, 2/21/23: I'm not sure why it has to be done here and not in
      // a different function, but if I don't do it here I get a seg fault.
      //void *vp=o2scl_settings.py_import_array();
      import_array();

      // Get the Unicode name of the user-specified module
      if (this->verbose>1) {
        std::cout << "Python version: "
                  << o2scl_settings.py_version() << std::endl;
        std::cout << "Staring kde_python::set_function()."
                  << std::endl;
        std::cout << "  Getting unicode for module named "
                  << module << std::endl;
      }
      p_name=PyUnicode_FromString(module.c_str());
      if (p_name==0) {
        O2SCL_ERR2("Create module name failed in ",
                   "kde_python::set_function().",
                   o2scl::exc_efailed);
      }
      
      // Import the user-specified module
      if (this->verbose>1) {
        std::cout << "  Importing module " << module << std::endl;
      }
      p_module=PyImport_Import(p_name);
      if (p_module==0) {
        O2SCL_ERR2("Load module failed in ",
                   "kde_python::set_function().",
                   o2scl::exc_efailed);
      }

      if (class_name.length()>0) {
        if (this->verbose>1) {
          std::cout << "  Obtaining python class " << class_name
                    << std::endl;
        }
        p_class=PyObject_GetAttrString(p_module,class_name.c_str());
        if (p_class==0) {
          O2SCL_ERR2("Get class failed in ",
                     "kde_python::set().",o2scl::exc_efailed);
        }
        
        // Create an instance of the class
        if (this->verbose>1) {
          std::cout << "  Loading python class." << std::endl;
        }
        if (PyCallable_Check(p_class)==false) {
          O2SCL_ERR2("Check class callable failed in ",
                     "funct_python_method::set_function().",
                     o2scl::exc_efailed);
        }
        
        if (this->verbose>1) {
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
      if (this->verbose>1) {
        std::cout << "  Making argument object for set function."
                  << std::endl;
      }
      p_set_args=PyTuple_New(3);
      if (p_set_args==0) {
        O2SCL_ERR2("Create arg tuple failed in ",
                   "kde_python::set_function().",
                   o2scl::exc_efailed);
      }

      // Setup the arguments to the python function
      if (this->verbose>1) {
        std::cout << "  Making argument object for log_pdf function."
                  << std::endl;
      }
      p_ld_args=PyTuple_New(1);
      if (p_ld_args==0) {
        O2SCL_ERR2("Create arg tuple failed in ",
                   "kde_python::set_function().",
                   o2scl::exc_efailed);
      }

      if (class_name.length()>0) {

        // Load the python function
        if (this->verbose>1) {
          std::cout << "  Loading python member function sample: "
                    << sample_func<< std::endl;
        }
        p_sample_func=PyObject_GetAttrString(p_instance,sample_func.c_str());
        if (p_sample_func==0) {
          O2SCL_ERR2("Get sample function failed in ",
                     "kde_python::set_function().",
                     o2scl::exc_efailed);
        }
        
        // Load the python function
        if (this->verbose>1) {
          std::cout << "  Loading python member function get: "
                    << ld_func<< std::endl;
        }
        p_ld_func=PyObject_GetAttrString(p_instance,ld_func.c_str());
        if (p_ld_func==0) {
          O2SCL_ERR2("Get get function failed in ",
                     "kde_python::set_function().",
                     o2scl::exc_efailed);
        }
        
        // Load the python function
        if (this->verbose>1) {
          std::cout << "  Loading python member function set: "
                    << set_func << std::endl;
        }
        p_set_func=PyObject_GetAttrString(p_instance,set_func.c_str());
        if (p_set_func==0) {
          O2SCL_ERR2("Get set function failed in ",
                     "kde_python::set_function().",
                     o2scl::exc_efailed);
        }

      } else {
        
        // Load the python function
        if (this->verbose>1) {
          std::cout << "  Loading python function set." << std::endl;
        }
        p_set_func=PyObject_GetAttrString(p_module,set_func.c_str());
        if (p_set_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "kde_python::set_function().",
                     o2scl::exc_efailed);
        }

        // Load the python function
        if (this->verbose>1) {
          std::cout << "  Loading python function sample." << std::endl;
        }
        p_sample_func=PyObject_GetAttrString(p_module,sample_func.c_str());
        if (p_sample_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "kde_python::set_function().",
                     o2scl::exc_efailed);
        }

        // Load the python function
        if (this->verbose>1) {
          std::cout << "  Loading python function get." << std::endl;
        }
        p_ld_func=PyObject_GetAttrString(p_module,ld_func.c_str());
        if (p_ld_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "kde_python::set_function().",
                     o2scl::exc_efailed);
        }
      }

      // Swap the tensor data
      swap(data,params);
      
      if (data.get_size(0)!=n_points) {
        std::cout << "Sizes: " << data.get_size(0) << " "
                  << data.get_size(1) << " " 
                  << n_points << " " << n_params << std::endl;
        O2SCL_ERR("Input data does not have correct number of rows.",
                  o2scl::exc_einval);
      }
      if (data.get_size(1)!=n_params) {
        std::cout << "Sizes: " << data.get_size(0) << " "
                  << data.get_size(1) << " " 
                  << n_points << " " << n_params << std::endl;
        O2SCL_ERR("Input data does not have correct number of columns.",
                  o2scl::exc_einval);
      }

      // First argument to set function: the data
      
      npy_intp data_dims[]={(npy_intp)data.get_size(0),
        (npy_intp)data.get_size(1)};
      if (this->verbose>1) {
        std::cout << "kde_python::set_function():" << std::endl;
      }
      PyObject *array_in=PyArray_SimpleNewFromData
        (2,data_dims,NPY_DOUBLE,(void *)(&(data.get_data()[0])));
         
      int pret=PyTuple_SetItem(p_set_args,0,array_in);
      if (pret!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "kde_python::set_function().",o2scl::exc_efailed);
      }
      
      // Second argument to set function: the bandwidth array
      
      npy_intp extra_data_dims[]={(npy_intp)(array.size())};
      PyObject *extra_array_in=PyArray_SimpleNewFromData
        (1,extra_data_dims,NPY_DOUBLE,(void *)(&(array[0])));
      
      int ret2=PyTuple_SetItem(p_set_args,1,extra_array_in);
      if (ret2!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "kde_python::set_function().",o2scl::exc_efailed);
      }
      
      // Third argument to set function: the options string
      
      if (this->verbose>1) {
        std::cout << "Creating python unicode for string: "
                  << options.length() << " " << options << std::endl;
      }
      PyObject *p_options=PyUnicode_FromString(options.c_str());
      if (p_options==0) {
        O2SCL_ERR2("String creation failed in ",
                   "kde_python::set().",o2scl::exc_efailed);
      }
      
      int ret3=PyTuple_SetItem(p_set_args,2,p_options);
      if (ret3!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "kde_python::set_function().",o2scl::exc_efailed);
      }

      // Call the python set function
      
      if (this->verbose>1) {
        std::cout << "  Calling python set function." << std::endl;
      }
      PyObject *result=PyObject_CallObject(p_set_func,p_set_args);
      if (result==0) {
        O2SCL_ERR2("Function set call failed in ",
                   "kde_python::set_function().",o2scl::exc_efailed);
      }

      if (this->verbose>1) {
        std::cout << "Done with kde_python::set_function()."
                  << std::endl;
      }
      
      return 0;
    }

    
    /// Return the dimensionality
    virtual size_t dim() const {
      return n_params;
    }

    /// The normalized density 
    void *log_pdf_internal(const vec_t &x, double &dout) const {

      if (x.size()!=n_params) {
        O2SCL_ERR("Input vector does not have correct size.",
                  o2scl::exc_einval);
      }
  
      if (p_set_func==0 || p_ld_func==0) {
        O2SCL_ERR2("No functions found in ",
                   "kde_python::log_pdf().",
                   o2scl::exc_efailed);
      }

      std::vector<double> x2(x.size());
      if (array_copy) {
        o2scl::vector_copy(x.size(),x,x2);
      }

      import_array();

      npy_intp x_dims[]={(npy_intp)x.size()};
      if (this->verbose>1) {
        std::cout << "kde_python::log_pdf():" << std::endl;
        std::cout << "  Array x: " << x.size() << std::endl;
      }
      
      PyObject *array_x;
      if (array_copy) {
        array_x=PyArray_SimpleNewFromData
          (1,x_dims,NPY_DOUBLE,(void *)(&(x2[0])));
      } else {
        array_x=PyArray_SimpleNewFromData
          (1,x_dims,NPY_DOUBLE,(void *)(&(x[0])));
      }

      int ret=PyTuple_SetItem(p_ld_args,0,array_x);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "kde_python::log_pdf().",o2scl::exc_efailed);
      }
      
      // Call the python function
      if (this->verbose>1) {
        std::cout << "  Calling python ld function." << std::endl;
      }
      PyObject *result=PyObject_CallObject(p_ld_func,p_ld_args);
      if (result==0) {
        O2SCL_ERR2("Function ld call failed in ",
                   "kde_python::log_pdf().",o2scl::exc_efailed);
      }

      if (this->verbose>1) {
        std::cout << "  Obtaining output 1." << std::endl;
      }
      dout=PyFloat_AsDouble(result);

      if (this->verbose>1) {
        std::cout << "  Decref result." << std::endl;
      }
      Py_DECREF(result);
  
      if (this->verbose>1) {
        std::cout << "Done in kde_python::log_pdf()."
                  << std::endl;
      }

      return 0;
    }

    /// The normalized density 
    virtual double log_pdf(const vec_t &x) const {
      double dout;
      log_pdf_internal(x,dout);
      return dout;
    }
    
    /// Get the bandwidth
    virtual double get_bandwidth() const {

      // Load the python function
      if (this->verbose>1) {
        std::cout << "  Loading python member function get_bandwidth."
                  << std::endl;
      }
      PyObject *p_gb_func=PyObject_GetAttrString(p_instance,"get_bandwidth");
      if (p_gb_func==0) {
        O2SCL_ERR2("Get get function failed in ",
                   "kde_python::set_function().",
                   o2scl::exc_efailed);
      }
      
      // Call the python function
      if (this->verbose>1) {
        std::cout << "  Calling python gb function." << std::endl;
      }
      PyObject *result=PyObject_CallObject(p_gb_func,0);
      if (result==0) {
        O2SCL_ERR2("Function gb call failed in ",
                   "kde_python::get_bandwidth().",o2scl::exc_efailed);
      }

      if (this->verbose>1) {
        std::cout << "  Obtaining output 1." << std::endl;
      }
      double dret=PyFloat_AsDouble(result);

      if (this->verbose>1) {
        std::cout << "  Decref result." << std::endl;
      }
      Py_DECREF(result);
  
      if (this->verbose>1) {
        std::cout << "  Decref func." << std::endl;
      }
      Py_DECREF(p_gb_func);
  
      if (this->verbose>1) {
        std::cout << "Done in kde_python::get_bandwidth()."
                  << std::endl;
      }

      return dret;
    }
  
    /// The normalized density 
    virtual double pdf(const vec_t &x) const {
      double val=log_pdf(x);
      if (!std::isfinite(val)) {
        O2SCL_ERR2("Log PDF not finite or negative in ",
                   "prob_dens_mdim::log_pdf().",o2scl::exc_efailed);
      }
      double val2=exp(val);
      if (!std::isfinite(val2)) {
        std::cout << val << " " << val2 << std::endl;
        O2SCL_ERR2("PDF not finite in ",
                   "prob_dens_mdim::log_pdf().",o2scl::exc_efailed);
      }
      return val2;
    }
  
    /** \brief Sample the distribution
        
        This is a void * version of <tt>operator()</tt> which
        doesn't cause problems with the import_array() macro.
     */
    void* operator2(vec_t &x) const {

      import_array();

      if (p_sample_func==0) {
        O2SCL_ERR2("No functions found in ",
                   "kde_python::operator2().",
                   o2scl::exc_efailed);
      }
      
      // Call the python function
      if (this->verbose>1) {
        std::cout << "  Calling python sample function." << std::endl;
        std::cout << p_sample_func << std::endl;
        std::cout << data.get_rank() << " " << data.get_size(0) << " "
                  << data.get_size(1) << std::endl;
      }
      PyObject *result=PyObject_CallObject(p_sample_func,0);
      if (result==0) {
        O2SCL_ERR2("Function sample call failed in ",
                   "kde_python::operator().",o2scl::exc_efailed);
      }
      
      if (PyArray_Check(result)==0) {
        O2SCL_ERR2("Function call did not return a numpy array in ",
                   "kde_python::operator2().",o2scl::exc_efailed);
      }
      
      if (this->verbose>1) {
        std::cout << "  Obtaining output 1." << std::endl;
      }
      void *vp=PyArray_DATA((PyArrayObject *)result);
      double *dp=(double *)vp;
      for(size_t i=0;i<n_params;i++) {
        x[i]=dp[i];
        if (this->verbose>1) {
          std::cout << "  i,y[i]: " << i << " " << x[i] << std::endl;
        }
      }
      
      if (this->verbose>1) {
        std::cout << "  Decref result." << std::endl;
      }
      Py_DECREF(result);
      
      if (this->verbose>1) {
        std::cout << "Done in kde_python::operator2()."
                  << std::endl;
      }

      return 0;
    }

    /// Sample the distribution
    virtual void operator()(vec_t &x) const {
      void *vp=operator2(x);
      return;
    }
    
    
  private:

    kde_python(const kde_python &);
    kde_python& operator=(const kde_python&);

  };
  
#endif
  
}
    
#endif

