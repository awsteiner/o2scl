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
#ifndef O2SCL_GMM_PYTHON_H
#define O2SCL_GMM_PYTHON_H

/** \file gmm_python.h
    \brief File defining \ref o2scl::gmm_python
*/

#include <iostream>
#include <string>
#include <cmath>

#include <o2scl/err_hnd.h>
#include <o2scl/tensor.h>
#include <o2scl/exp_max.h>
#include <o2scl/kde_python.h>

#ifdef O2SCL_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

namespace o2scl {

#ifdef O2SCL_PYTHON

  /** \brief Multidimensional interpolation interface for python
   */
  class gmm_python {
    
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
    PyObject *p_components_args;

    /// Function arguments
    PyObject *p_get_args;

    /// Python function
    PyObject *p_set_func;

    /// Python function
    PyObject *p_components_func;

    /// Python function
    PyObject *p_get_func;

    /// Number of parameters
    size_t n_params;
    
    /// Number of points
    size_t n_points;

    /// Number of components
    size_t n_components;

    /// Underlying probability distribution
    prob_dens_mdim_gmm<> pdm_gmm;

    /** \brief An internal version of \ref get_python which returns a 
        void * for the import_array() macro
     */
    void *get_python_internal(int &ret);
    
    /** \brief Internal version of \ref set_function()
        which returns a void * for the import_array() macro

        Note that the \c params object is only used for the
        duration of this function and can be destroyed afterwards
    */
    void *set_function_internal
    (std::string module, size_t n_comp, const o2scl::tensor<> &params,
     std::string options, std::string class_name, int v, int &ret);
      
  public:

    gmm_python();
    
    virtual ~gmm_python();

    /** \brief Specify the Python module and function
     */
    gmm_python(std::string module, size_t n_comp,
               const o2scl::tensor<> &params, std::string options="",
               std::string class_name="", int v=0);

    /// The name of the Python set function (default "set_data_str")
    std::string set_func;
    
    /// The name of the Python get function (default "get_data")
    std::string get_func;
    
    /// The name of the Python components function (default "components")
    std::string components_func;
    
    /// The name of the Python components function (default "log_pdf")
    std::string log_pdf_func;
    
    /** \brief Free the associated memory
     */
    void free();
    
    /// Verbosity parameter
    int verbose;

    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, size_t n_comp,
                     const o2scl::tensor<> &params,
                     std::string options="",
                     std::string class_name="", int v=0);

    /** \brief At point \c x, compute the components for 
        each Gaussian in the mixture and place the result in \c y
    */
    virtual int components(const std::vector<double> &x,
                           std::vector<double> &y) const;
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual int get_python();

    /** \brief Get the underlying Gaussian mixture probability
        density
     */
    const prob_dens_mdim_gmm<> &get_gmm() {
      return pdm_gmm;
    }
    
  private:

    gmm_python(const gmm_python &);
    gmm_python& operator=(const gmm_python&);

  };

#endif
  
}
    
#endif

