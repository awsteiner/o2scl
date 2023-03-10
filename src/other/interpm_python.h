/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
#ifndef O2SCL_INTERPM_PYTHON_H
#define O2SCL_INTERPM_PYTHON_H

/** \file interpm_python.h
    \brief File defining \ref o2scl::interpm_python
*/

#include <iostream>
#include <string>
#include <cmath>

#include <o2scl/err_hnd.h>
#include <o2scl/tensor.h>

#ifdef O2SCL_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

namespace o2scl {

#ifdef O2SCL_PYTHON

  /** \brief Multidimensional interpolation interface for python
   */
  class interpm_python {
    
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
    PyObject *p_eval_args;

    /// Python function
    PyObject *p_set_func;

    /// Python function
    PyObject *p_eval_func;

    /// Verbosity parameter
    int verbose;

    /// Number of parameters
    size_t n_params;
    
    /// Number of outputs
    size_t n_outputs;
    
    /// Number of points
    size_t n_points;
    
  public:

    interpm_python();
    
    /** \brief Specify the Python module and function
     */
    interpm_python(std::string module, std::string set_func,
                   std::string eval_func, 
                   size_t n_pars, size_t n_dat, size_t n_out,
                   const o2scl::tensor<> &params,
                   const o2scl::tensor<> &outputs,
                   std::string options="", 
                   std::string class_name="", int v=0);
    
    void free();
    
    virtual ~interpm_python();
  
    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string set_func,
                     std::string eval_func, 
                     size_t n_pars, size_t n_dat, size_t n_out,
                     const o2scl::tensor<> &params,
                     const o2scl::tensor<> &outputs,
                     std::string options="",
                     std::string class_name="", int v=0);
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual int eval(const std::vector<double> &x,
                     std::vector<double> &y) const;
    
  private:

    interpm_python(const interpm_python &);
    interpm_python& operator=(const interpm_python&);

  };
  
#endif
  
}
    
#endif

