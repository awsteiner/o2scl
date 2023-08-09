/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_FUNCT_H
#define O2SCL_FUNCT_H

/** \file funct.h
    \brief Function object classes for one-dimensional functions
*/

#include <string>
#include <functional>

#include <gsl/gsl_math.h>

#include <o2scl/err_hnd.h>
#include <o2scl/set_python.h>

#ifdef O2SCL_SET_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#endif

namespace o2scl {

  // Experimental template typedef
  //template<class fp_t> using funct2 = std::function<fp_t(fp_t)>;
  
  /// One-dimensional function typedef in src/base/funct.h
  typedef std::function<double(double)> funct;

  /** \brief A wrapper to specify \ref o2scl::funct objects 
      to GSL
  */
  class funct_gsl : public gsl_function {

  protected:
    
    /// The function wrapper
    static double funct_wrap(double x, void *params) {
      funct *fp=(funct *)params;
      return (*fp)(x);
    }

  public:

    /// Create an object based on the specified function, \c f
    funct_gsl(funct &f) {
      function=&funct_wrap;
      params=&f;
    }
    
  };

#ifdef O2SCL_SET_PYTHON
  
  /** \brief One-dimensional function from a python function
   */
  class funct_python {

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
    funct_python(std::string module, std::string func, int v=0);
    
    virtual ~funct_python();
  
    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string func);
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const;

  protected:

    funct_python() {};

  private:

    funct_python(const funct_python &);
    funct_python& operator=(const funct_python&);

  };

  /** \brief One-dimensional function from a python function
   */
  class funct_python_method {

  protected:

    /// Python unicode object containing function name
    PyObject *pName;
    
    /// Python module containing function
    PyObject *pModule;
    
    /// Function arguments
    PyObject *pArgs;

    /// Python function
    PyObject *pFunc;

    /// Python class instance
    PyObject *pInstance;

    /// Python class
    PyObject *pClass;

    /// Verbosity parameter
    int verbose;
    
  public:
    
    /** \brief Specify the python and the parameters
     */
    funct_python_method(std::string module, std::string class_name,
                        std::string func, int v=0);
    
    virtual ~funct_python_method();
  
    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string class_name,
                     std::string func);
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const;

  protected:

    funct_python_method() {};

  private:

    funct_python_method(const funct_python_method &);
    funct_python_method& operator=(const funct_python_method&);

  };

#endif

  /// One-dimensional function typedef in src/base/funct.h
  typedef std::function<int(double,double &)> funct_ret;

}

#endif
