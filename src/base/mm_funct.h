/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
*/
#ifndef O2SCL_MM_FUNCT_H
#define O2SCL_MM_FUNCT_H

/** \file mm_funct.h
    \brief Function object classes for multi-dimensional functions
*/

#include <string>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/calc_utf8.h>
#include <o2scl/lib_settings.h>

#ifdef O2SCL_PYTHON
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

namespace o2scl {

  /// Array of multi-dimensional functions typedef in src/base/mm_funct.h
  typedef std::function<
    int(size_t,const boost::numeric::ublas::vector<double> &,
	boost::numeric::ublas::vector<double> &) > mm_funct;

  /** \brief Array of multi-dimensional functions in an array of strings
   */
  class mm_funct_strings {

  public:
      
    /** \brief Specify the strings
     */
    template<class vec_string_t=std::vector<std::string> >
      mm_funct_strings(int nv, vec_string_t &exprs,
			 vec_string_t &var_arr) {

      st_nv=nv;
      st_forms.resize(nv);
      st_vars.resize(nv);
      calc.resize(nv);
      for (int i=0;i<nv;i++) {
	calc[i].compile(exprs[i].c_str(),&vars);
	st_vars[i]=var_arr[i];
	st_forms[i]=exprs[i];
      }
    }
      
    virtual ~mm_funct_strings() {
    };
      
    /** \brief Set the values of the auxilliary parameters that were
	specified in 'parms' in the constructor
    */
    int set_parm(std::string name, double val) {
      vars[name]=val;
      return 0;
    }
      
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    template<class vec_t=boost::numeric::ublas::vector<double> >
      int operator()(size_t nv, const vec_t &x, vec_t &y) {

      for(size_t i=0;i<nv;i++) {
	vars[st_vars[i]]=x[i];
      }
      for(size_t i=0;i<nv;i++) {
	y[i]=calc[i].eval(&vars);
      }
      return 0;
    }
      
    /// Set the functions
    template<class vec_string_t=std::vector<std::string> >
      void set_function(int nv, vec_string_t &exprs,
			vec_string_t &var_arr) {

      st_nv=nv;
      st_forms.resize(nv);
      st_vars.resize(nv);
      calc.resize(nv);
      for (int i=0;i<nv;i++) {
	calc[i].compile(exprs[i],&vars);
	st_vars[i]=var_arr[i];
	st_forms[i]=exprs[i];
      }

      return;
    }
      
#ifndef DOXYGEN_INTERNAL
      
  protected:
      
    /// The function parsers
    std::vector<calc_utf8<> > calc;
      
    /// External variables to include in the function parsing
    std::map<std::string,double> vars;
      
    /// The expressions
    std::vector<std::string> st_forms;
      
    /// The variables
    std::vector<std::string> st_vars;
      
    /// The number of variables
    int st_nv;
      
    mm_funct_strings() {};
      
  private:
      
    mm_funct_strings(const mm_funct_strings &);
    mm_funct_strings& operator=(const mm_funct_strings&);
      
#endif
      
  };

#ifdef O2SCL_PYTHON
  
  /** \brief One-dimensional function from a Python function

      This class allows one to specify a Python function from a module
      (or optionally a class instide that module) and call that
      function from C++. The python function must have one positional
      argument which is a Python list (no keyword arguments will be
      passed) and it must return a Python list.

      If a class is specified in the constructor or with 
      \ref set_function(), then an instance of that Python class
      is automatically constructed.

      The Python library must be initialized before using
      this class (for example with \ref o2scl_settings::py_init())
      and finalized after this class calls its destructor (for
      example with \ref o2scl_settings::py_final()). 

      \note This class presumes that the input and output lists occupy
      different memory (i.e. they are different python objects). If
      this is not the case, the results will be unpredictable.

      \warning This function requires copying all of the elements to
      the C++ array to a Python object and vice versa. The 
      \ref o2scl::mm_funct_python_ndarray requires less copying.

      \future Find a way to transmit Python exception information
      back to this class and see if it can be handled with a C++
      try block. 
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class mm_funct_python {
    
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
    PyObject *p_args;

    /// Python function
    PyObject *p_func;

    /// Verbosity parameter
    int verbose;
    
  public:
    
    /** \brief Specify the Python module and function
     */
    mm_funct_python(std::string module="", std::string func="",
                    std::string class_name="", int v=0) {
                    
      verbose=v;

      if (o2scl_settings.py_initialized==false) {
        if (verbose>0) {
          std::cout << "Running py_init()." << std::endl;
        }
        o2scl_settings.py_init();
      }
      p_func=0;
      p_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      p_name=0;
      if (module.length()>0) {
        set_function(module,func,class_name);
      }
    }      
    
    void free() {
      if (p_func!=0) {
        if (verbose>0) {
          std::cout << "Decref func." << std::endl;
        }
        Py_DECREF(p_func);
      }
      if (p_args!=0) {
        if (verbose>0) {
          std::cout << "Decref args." << std::endl;
        }
        Py_DECREF(p_args);
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
      p_func=0;
      p_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      p_name=0;
      if (verbose>0) {
        std::cout << "Done in mm_funct_python::free()." << std::endl;
      }
    }      
    
    virtual ~mm_funct_python() {
      free();
    }      
  
    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string func,
                     std::string class_name, int v=0) {

      verbose=v;
      free();
      
      // Get the Unicode name of the user-specified module
      if (verbose>0) {
        std::cout << "Getting unicode for module name()." << std::endl;
      }
      p_name=PyUnicode_FromString(module.c_str());
      if (p_name==0) {
        O2SCL_ERR2("Create module name failed in ",
                   "mm_funct_python::set_function().",o2scl::exc_efailed);
      }
      
      // Import the user-specified module
      if (verbose>0) {
        std::cout << "Importing module." << std::endl;
      }
      p_module=PyImport_Import(p_name);
      if (p_module==0) {
        O2SCL_ERR2("Load module failed in ",
                   "mm_funct_python::set_function().",o2scl::exc_efailed);
      }

      if (class_name.length()>0) {
        if (verbose>0) {
          std::cout << "Obtaining python class." << std::endl;
        }
        p_class=PyObject_GetAttrString(p_module,class_name.c_str());
        if (p_class==0) {
          O2SCL_ERR2("Get class failed in ",
                     "emulator_python::set().",o2scl::exc_efailed);
        }
        
        // Create an instance of the class
        if (verbose>0) {
          std::cout << "Loading python class." << std::endl;
        }
        if (PyCallable_Check(p_class)==false) {
          O2SCL_ERR2("Check class callable failed in ",
                     "funct_python_method::set_function().",
                     o2scl::exc_efailed);
        }
        
        if (verbose>0) {
          std::cout << "Loading python class instance." << std::endl;
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
        std::cout << "Getting arguments for python function." << std::endl;
      }
      p_args=PyTuple_New(1);
      if (p_args==0) {
        O2SCL_ERR2("Create arg tuple failed in ",
                   "mm_funct_python::set_function().",o2scl::exc_efailed);
      }

      if (class_name.length()>0) {
        // Load the python function
        if (verbose>0) {
          std::cout << "Loading python function." << std::endl;
        }
        p_func=PyObject_GetAttrString(p_instance,func.c_str());
        if (p_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "mm_funct_python::set_function().",
                     o2scl::exc_efailed);
        }
      } else {
        // Load the python function
        if (verbose>0) {
          std::cout << "Loading python function." << std::endl;
        }
        p_func=PyObject_GetAttrString(p_module,func.c_str());
        if (p_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "mm_funct_python::set_function().",o2scl::exc_efailed);
        }
      }      
      return 0;
    }
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual int operator()(size_t n, const vec_t &v,
                           vec_t &y) const {
      
      if (p_func==0) {
        O2SCL_ERR2("No function found in ",
                   "mm_funct_python::operator().",o2scl::exc_efailed);
      }
      
      // Create the list object
      PyObject *pList_x=PyList_New(n);
      if (pList_x==0) {
        O2SCL_ERR2("List creation for x failed in ",
                   "mm_funct_python::operator().",o2scl::exc_efailed);
      }

      // Create a python object from the vector
      std::vector<PyObject *> p_values(n);
      if (verbose>0) {
        std::cout << "Creating python object from vector." << std::endl;
      }
      
      for(size_t i=0;i<n;i++) {
        p_values[i]=PyFloat_FromDouble(v[i]);
        if (p_values[i]==0) {
          O2SCL_ERR2("Value creation failed in ",
                     "mm_funct_python::operator().",o2scl::exc_efailed);
        }
        
        // AWS, 2/7/23: Set the python function arguments. Note that
        // after we use this function, the list "steals" the memory
        // for the individual elements, so no need to decref them
        // below. See https://pythonextensionpatterns.readthedocs.io/
        // en/latest/refcount.html for more.
        
        int iret=PyList_SetItem(pList_x,i,p_values[i]);
        if (iret!=0) {
          O2SCL_ERR2("Item set failed in ",
                     "mm_funct_python::operator().",o2scl::exc_efailed);
        }
      }
      
      int ret=PyTuple_SetItem(p_args,0,pList_x);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "mm_funct_python::operator().",o2scl::exc_efailed);
      }

      // Call the python function
      if (verbose>0) {
        std::cout << "Call python function." << std::endl;
      }
      PyObject *result=PyObject_CallObject(p_func,p_args);
      if (result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "mm_funct_python::operator().",o2scl::exc_efailed);
      }

      if (PyList_Check(result)==0) {
        O2SCL_ERR2("Function call did not return a list in ",
                   "mm_funct_python::operator().",o2scl::exc_efailed);
      }
      
      for(size_t i=0;i<n;i++) {
        PyObject *yval=PyList_GetItem(result,i);
        if (yval==0) {
          O2SCL_ERR2("Failed to get y list value in ",
                     "mm_funct_python::operator().",o2scl::exc_efailed);
        }
        y[i]=PyFloat_AsDouble(yval);
        /*
          AWS, 2/7/23: I believe the list maintains the ownership of 
          the memory so there's no need to decref the individual 
          elements
          
          if (verbose>0) {
          std::cout << "Decref yval " << i << " of " << p_values.size()
          << std::endl;
          }
          Py_DECREF(yval);
          std::cout << i << " " << yval << std::endl;
        */
      }

      /*
        AWS, 2/7/23: I believe the list maintains the ownership of 
        the memory so there's no need to decref the individual 
        elements

        for(size_t i=0;i<p_values.size();i++) {
        if (verbose>0) {
        std::cout << "Decref value " << i << " of " << p_values.size()
        << std::endl;
        }
        std::cout << i << " " << p_values[i] << std::endl;
        Py_DECREF(p_values[i]);
        }
      */
      
      if (verbose>0) {
        std::cout << "Decref list." << std::endl;
      }
      Py_DECREF(pList_x);
      
      if (verbose>0) {
        std::cout << "Decref result." << std::endl;
      }
      Py_DECREF(result);
  
      if (verbose>0) {
        std::cout << "Done in mm_funct_python::operator()." << std::endl;
      }

      return 0;
    }      

  protected:

    mm_funct_python() {};

  private:

    mm_funct_python(const mm_funct_python &);
    mm_funct_python& operator=(const mm_funct_python&);

  };

  /** \brief One-dimensional function from a Python function
      which uses numpy arrays

      This class allows one to specify a Python function from a module
      (or optionally a class instide that module) and call that
      function from C++. The python function must have one positional
      argument which is a numpy array (no keyword arguments will be
      passed) and it must return a numpy array.

      If a class is specified in the constructor or with 
      \ref set_function(), then an instance of that Python class
      is automatically constructed.

      The Python library must be initialized before using
      this class (for example with \ref o2scl_settings::py_init())
      and finalized after this class calls its destructor (for
      example with \ref o2scl_settings::py_final()). 

      \note This class presumes that the input and output lists occupy
      different memory (i.e. they are different python objects). If
      this is not the case, the results will be unpredictable.

      \warning This function requires copying all of the elements to
      the C++ array to a Python object and vice versa. The 
      \ref o2scl::mm_funct_python_ndarray requires less copying.

      \future Find a way to transmit Python exception information
      back to this class and see if it can be handled with a C++
      try block. 
  */
  template<class vec_t=std::vector<double> >
  class mm_funct_python_ndarray {
    
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
    PyObject *p_args;

    /// Python function
    PyObject *p_func;

    /// Verbosity parameter
    int verbose;
    
  public:
    
    /** \brief Specify the Python module and function
     */
    mm_funct_python_ndarray(std::string module="", std::string func="",
                    std::string class_name="", int v=0) {
                    
      verbose=v;

      if (o2scl_settings.py_initialized==false) {
        if (verbose>0) {
          std::cout << "Running py_init()." << std::endl;
        }
        o2scl_settings.py_init();
      }
      p_func=0;
      p_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      p_name=0;
      if (module.length()>0) {
        set_function(module,func,class_name);
      }
    }      
    
    void free() {
      if (verbose>0) {
        std::cout << "Starting mm_funct_python_ndarray::free()." << std::endl;
      }
      if (p_func!=0) {
        if (verbose>0) {
          std::cout << "Decref func." << std::endl;
        }
        Py_DECREF(p_func);
      }
      if (p_args!=0) {
        if (verbose>0) {
          std::cout << "Decref args." << std::endl;
        }
        Py_DECREF(p_args);
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
      p_func=0;
      p_args=0;
      p_instance=0;
      p_class=0;
      p_module=0;
      p_name=0;
      if (verbose>0) {
        std::cout << "Done in mm_funct_python_ndarray::free()." << std::endl;
      }
    }      
    
    virtual ~mm_funct_python_ndarray() {
      free();
    }      
  
    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string func,
                     std::string class_name) {

      free();
      
      // Get the Unicode name of the user-specified module
      if (verbose>0) {
        std::cout << "Staring mm_funct_python_ndarray::set_function()."
                  << std::endl;
        std::cout << "  Getting unicode for module name()." << std::endl;
      }
      p_name=PyUnicode_FromString(module.c_str());
      if (p_name==0) {
        O2SCL_ERR2("Create module name failed in ",
                   "mm_funct_python_ndarray::set_function().",
                   o2scl::exc_efailed);
      }
      
      // Import the user-specified module
      if (verbose>0) {
        std::cout << "  Importing module." << std::endl;
      }
      p_module=PyImport_Import(p_name);
      if (p_module==0) {
        O2SCL_ERR2("Load module failed in ",
                   "mm_funct_python_ndarray::set_function().",
                   o2scl::exc_efailed);
      }

      if (class_name.length()>0) {
        if (verbose>0) {
          std::cout << "  Obtaining python class." << std::endl;
        }
        p_class=PyObject_GetAttrString(p_module,class_name.c_str());
        if (p_class==0) {
          O2SCL_ERR2("Get class failed in ",
                     "emulator_python::set().",o2scl::exc_efailed);
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
        std::cout << "  Making argument object for function." << std::endl;
      }
      p_args=PyTuple_New(1);
      if (p_args==0) {
        O2SCL_ERR2("Create arg tuple failed in ",
                   "mm_funct_python_ndarray::set_function().",
                   o2scl::exc_efailed);
      }

      if (class_name.length()>0) {
        // Load the python function
        if (verbose>0) {
          std::cout << "  Loading python member function." << std::endl;
        }
        p_func=PyObject_GetAttrString(p_instance,func.c_str());
        if (p_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "mm_funct_python_ndarray::set_function().",
                     o2scl::exc_efailed);
        }
      } else {
        // Load the python function
        if (verbose>0) {
          std::cout << "  Loading python function." << std::endl;
        }
        p_func=PyObject_GetAttrString(p_module,func.c_str());
        if (p_func==0) {
          O2SCL_ERR2("Get function failed in ",
                     "mm_funct_python_ndarray::set_function().",
                     o2scl::exc_efailed);
        }
      }

      // I'm not sure why it has to be done here and not in
      // mm_funct_ts.cpp
      import_array();
      
      if (verbose>0) {
        std::cout << "Done with mm_funct_python_ndarray::set_function()."
                  << std::endl;
      }
      
      return 0;
    }
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual int operator()(size_t n, const vec_t &v,
                           vec_t &y) const {

      if (p_func==0) {
        O2SCL_ERR2("No function found in ",
                   "mm_funct_python_ndarray::operator().",
                   o2scl::exc_efailed);
      }

      /*
        AWS, 2/12/23, Check: 
        https://stackoverflow.com/questions/52731884/pyarray-simplenewfromdata
        if I have trouble with memory leaks
       */

      npy_intp dims[]={(npy_intp)n};
      if (verbose>0) {
        std::cout << "mm_funct_python_ndarray::operator():" << std::endl;
        std::cout << "  Input array: " << v[0] << " " << v[1] << std::endl;
      }
      PyObject *array_in=PyArray_SimpleNewFromData
        (1,dims,NPY_DOUBLE,(void *)(&(v[0])));
         
      int ret=PyTuple_SetItem(p_args,0,array_in);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "mm_funct_python::operator().",o2scl::exc_efailed);
      }
      
      // Call the python function
      if (verbose>0) {
        std::cout << "  Calling python function." << std::endl;
      }
      PyObject *result=PyObject_CallObject(p_func,p_args);
      if (result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "mm_funct_python_ndarray::operator().",o2scl::exc_efailed);
      }

      if (PyArray_Check(result)==0) {
        O2SCL_ERR2("Function call did not return a numpy array in ",
                   "mm_funct_python_ndarray::operator().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "  Obtaining output." << std::endl;
      }
      for(size_t i=0;i<n;i++) {
        void *vp=PyArray_GETPTR1(result,i);
        double *dp=(double *)vp;
        y[i]=*dp;
        std::cout << "  i,y[i]: " << i << " " << y[i] << std::endl;
      }
      
      if (verbose>0) {
        std::cout << "  Decref result." << std::endl;
      }
      Py_DECREF(result);
  
      if (verbose>0) {
        std::cout << "Done in mm_funct_python_ndarray::operator()."
                  << std::endl;
      }

      return 0;
    }      

  protected:

    mm_funct_python_ndarray() {};

  private:

    mm_funct_python_ndarray(const mm_funct_python_ndarray &);
    mm_funct_python_ndarray& operator=(const mm_funct_python_ndarray&);

  };

#endif
  
#ifdef O2SCL_NEVER_DEFINED
  /** \brief A wrapper to specify \ref o2scl::mm_funct-like objects 
      to GSL
  */
  template<class vec_t>
    class mm_funct_gsl : public gsl_multiroot_function {
    
  public:
    
    typedef std::function<int(size_t,const vec_t &, vec_t &)> func_t;

  protected:
    
    /// The function wrapper
    static int funct_wrap(const gsl_vector *x, void *params,
			  gsl_vector *f) {
      func_t *fp=(func_t *)params;
      vec_t x2(x->size), f2(x->size);
      o2scl::vector_copy<double *,vec_t>(x->size,x.data,x2);
      int ret=(*fp)(x->size,x2,f2);
      o2scl::vector_copy<vec_t,double *>(x->size,f2,f.data);
      return ret;
    }

  public:

    /// Create an object based on the specified function, \c f
    funct_gsl(func_t &f) {
      function=&funct_wrap;
      params=&f;
    }
    
  };
#endif

}

#endif
