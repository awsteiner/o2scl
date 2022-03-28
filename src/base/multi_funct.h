/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_MULTI_FUNCT_H
#define O2SCL_MULTI_FUNCT_H

/** \file multi_funct.h
    \brief Function object classes for a multi-dimensional function
*/

#include <string>
#include <functional>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/err_hnd.h>
#include <o2scl/calc_utf8.h>
#include <o2scl/lib_settings.h>

#ifdef O2SCL_PYTHON
#include <Python.h>
#endif

namespace o2scl {

  /// Multi-dimensional function typedef in src/base/multi_funct.h
  typedef std::function<
    double(size_t,const boost::numeric::ublas::vector<double> &)>
  multi_funct;
  
  /** \brief A multi-dimensional function from a string
   */
  class multi_funct_strings {
    
  public:
  
    /** \brief Specify the string and the parameters
     */
    template<class vec_string_t=std::vector<std::string> >
    multi_funct_strings(std::string expr, int nv,
                        vec_string_t &var_arr) {
    
      st_nv=nv;
      st_funct=expr;
      st_vars.resize(nv);
      for (int i=0;i<nv;i++) {
	calc.compile(expr.c_str(),&vars);
	st_vars[i]=var_arr[i];
      }
    }
  
    /** \brief Specify the string and the parameters
     */
    template<class vec_string_t=std::vector<std::string> >
    void set_function(std::string expr, int nv, vec_string_t &var_arr) {

      st_nv=nv;
      st_funct=expr;
      st_vars.resize(nv);
      for (int i=0;i<nv;i++) {
	calc.compile(expr.c_str(),&vars);
	st_vars[i]=var_arr[i];
      }
      return;
    }

    virtual ~multi_funct_strings() {
    };
  
    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms  in the constructor
    */
    int set_parm(std::string name, double val) {
      vars[name]=val;
      return 0;
    }

    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
    template<class vec_t=boost::numeric::ublas::vector<double> >
    double operator()(size_t nv, const vec_t &x) {

      for(size_t i=0;i<nv;i++) {
	vars[st_vars[i]]=x[i];
      }

      return calc.eval(&vars);
    }

  protected:

    /// The function parser
    calc_utf8<> calc;

    /// External variables to include in the function parsing
    std::map<std::string,double> vars;

    /// The number of variables
    int st_nv;

    /// The function string
    std::string st_funct;
    
    /// The variable string
    std::vector<std::string> st_vars;
  
    multi_funct_strings() {}
  
  private:

    multi_funct_strings(const multi_funct_strings &);
    multi_funct_strings& operator=(const multi_funct_strings&);

  };

  /** \brief One-dimensional function from a python function
   */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class multi_funct_python {
    
  protected:

    /// Python unicode object containing function name
    PyObject *pName;
    
    /// Python module containing function
    PyObject *pModule;
    
    /// Function arguments
    PyObject *pArgs;

    /// Python function
    PyObject *pFunc;

    /// Verbosity parameter
    int verbose;
    
  public:
    
    /** \brief Specify the python and the parameters
     */
    multi_funct_python(std::string module, std::string func, int v=0) {
      verbose=v;
      
      if (o2scl_settings.py_initialized==false) {
        if (verbose>0) {
          std::cout << "Running py_init()." << std::endl;
        }
        o2scl_settings.py_init();
      }
      set_function(module,func);
    }      
    
    virtual ~multi_funct_python() {
      if (verbose>0) {
        std::cout << "Decref func." << std::endl;
      }
      Py_DECREF(pFunc);
      if (verbose>0) {
        std::cout << "Decref args." << std::endl;
      }
      Py_DECREF(pArgs);
      if (verbose>0) {
        std::cout << "Decref module." << std::endl;
      }
      Py_DECREF(pModule);
      if (verbose>0) {
        std::cout << "Decref name." << std::endl;
      }
      Py_DECREF(pName);
      if (verbose>0) {
        std::cout << "Done in multi_funct_python destructor." << std::endl;
      }
    }      
  
    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string func) {
      
      // Get the Unicode name of the user-specified module
      if (verbose>0) {
        std::cout << "Getting unicode for module name()." << std::endl;
      }
      pName=PyUnicode_FromString(module.c_str());
      if (pName==0) {
        O2SCL_ERR2("Create module name failed in ",
                   "multi_funct_python::set_function().",o2scl::exc_efailed);
      }
      
      // Import the user-specified module
      if (verbose>0) {
        std::cout << "Importing module." << std::endl;
      }
      pModule=PyImport_Import(pName);
      if (pModule==0) {
        O2SCL_ERR2("Load module failed in ",
                   "multi_funct_python::set_function().",o2scl::exc_efailed);
      }

      // Setup the arguments to the python function
      if (verbose>0) {
        std::cout << "Getting arguments for python function." << std::endl;
      }
      pArgs=PyTuple_New(1);
      if (pArgs==0) {
        O2SCL_ERR2("Create arg tuple failed in ",
                   "multi_funct_python::set_function().",o2scl::exc_efailed);
      }
      
      // Load the python function
      if (verbose>0) {
        std::cout << "Loading python function." << std::endl;
      }
      pFunc=PyObject_GetAttrString(pModule,func.c_str());
      if (pFunc==0) {
        O2SCL_ERR2("Get function failed in ",
                   "multi_funct_python::set_function().",o2scl::exc_efailed);
      }
      
      return 0;
    }
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(size_t n, const vec_t &v) const {
      
      // Create the list object
      PyObject *pList=PyList_New(n);
      if (pList==0) {
        O2SCL_ERR2("List creation failed in ",
                   "multi_funct_python::operator().",o2scl::exc_efailed);
      }

      double y;
      
      // Create a python object from the vector
      std::vector<PyObject *> pValues(n);
      if (verbose>0) {
        std::cout << "Creating python object from vector." << std::endl;
      }
      
      for(size_t i=0;i<n;i++) {
        pValues[i]=PyFloat_FromDouble(v[i]);
        if (pValues[i]==0) {
          O2SCL_ERR2("Value creation failed in ",
                     "multi_funct_python::operator().",o2scl::exc_efailed);
        }
        
        // Set the python function arguments
        int iret=PyList_SetItem(pList,i,pValues[i]);
        if (iret!=0) {
          O2SCL_ERR2("Item set failed in ",
                     "multi_funct_python::operator().",o2scl::exc_efailed);
        }
      }
      
      int ret=PyTuple_SetItem(pArgs,0,pList);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "multi_funct_python::operator().",o2scl::exc_efailed);
      }

      // Call the python function
      if (verbose>0) {
        std::cout << "Call python function." << std::endl;
      }
      PyObject *result=PyObject_CallObject(pFunc,pArgs);
      if (result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "multi_funct_python::operator().",o2scl::exc_efailed);
      }

      y=PyFloat_AsDouble(result);

      for(size_t i=0;i<pValues.size();i++) {
        if (verbose>0) {
          std::cout << "Decref value " << i << " of " << pValues.size()
                    << std::endl;
        }
        Py_DECREF(pValues[i]);
      }
      if (verbose>0) {
        std::cout << "Decref list." << std::endl;
      }
      Py_DECREF(pList);
      
      if (verbose>0) {
        std::cout << "Decref result." << std::endl;
      }
      Py_DECREF(result);
  
      if (verbose>0) {
        std::cout << "Done in multi_funct_python::operator()." << std::endl;
      }

      return y;
    }      

  protected:

    multi_funct_python() {};

  private:

    multi_funct_python(const multi_funct_python &);
    multi_funct_python& operator=(const multi_funct_python&);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
