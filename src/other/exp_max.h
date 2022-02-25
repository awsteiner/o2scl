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
#ifndef O2SCL_EXP_MAX_H
#define O2SCL_EXP_MAX_H

/** \file exp_max.h
    \brief File defining \ref o2scl::exp_max
*/

#include <iostream>
#include <string>
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_combination.h>

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>
#include <o2scl/vec_stats.h>
#include <o2scl/linear_solver.h>
#include <o2scl/columnify.h>
#include <o2scl/table.h>
#include <o2scl/tensor.h>

namespace o2scl {

  /** \brief Expectation maximization
   */
  template<class mat_t=const_matrix_view_table<> >
  class exp_mat {
    
  protected:
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    
    exp_max() {
      data_set=false;
    }
    
    /** \brief Verbosity parameter (default 0)
     */
    int verbose;
    
    /** \brief Initialize the data
        The object \c vecs should be a matrix with a
        first index of size <tt>n_in+n_out</tt> and a second 
        index of size <tt>n_points</tt>. It may have be
        any type which allows the use of <tt>operator(,)</tt>
        and <tt>std::swap</tt>.
    */
    void set_data(size_t n_in, size_t n_points,
                  mat_t &dat) {
      
      if (n_points<points) {
        O2SCL_ERR2("Not enough points provided in ",
                   "exp_max::set_data()",exc_efailed);
      }
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
                   "exp_max::set_data()",exc_efailed);
      }
      if (n_out<1) {
        O2SCL_ERR2("Must provide at least one output column in ",
                   "exp_max::set_data()",exc_efailed);
      }
      np=n_points;
      nd_in=n_in;
      std::swap(data,dat);
      data_set=true;
      
      return;
    }

    /** \brief Get the data used for interpolation
     */
    void get_data(size_t &n_in, size_t &n_points,
                  mat_t &dat) {
      n_points=np;
      n_in=nd_in;
      std::swap(data,dat);
      data_set=false;
      n_points=0;
      n_in=0;
      return;
    }
    //@}

    /// \name Evaluate interpolation
    //@{
    /** \brief Perform the interpolation over the first function
     */
    template<class mat_t> int compute(size_t n_gauss,
                                      const mat_t &mean_init) {

      if (n_gauss==0) {
        O2SCL_ERR("Cannot select zero gaussians in compute().",
                  o2scl::exc_einval);
      }
      
      weights.resize(n_gauss);
      resps.resize(n_gauss,np);
      means.resize(n_gauss,nd_in);
      covars.resize(n_gauss,nd_in,nd_in);
      pdmg.resize(n_gauss);
      
      // ------------------------------------------------------------
      // Initialize by setting the responsibilities to 1 for the
      // closest Gaussian and 0 for the other Gaussians. Additionally,
      // set the weights equal to the number of points closest to each
      // mean. This code also works for the trivial single Gaussian
      // case.
      
      vector_set_all(n_gauss,weights,0.0);

      for(size_t i=0;i<np;i++) {
        double dist_min;
        size_t j_min=0;
        for(size_t j=0;j<n_gauss;j++) {
          double dist=0.0;
          for(size_t k=0;k<nd_in;k++) {
            dist+=pow(mean_init[j][k]-data[i][k],2.0);
          }
          if (j==0 || dist<dist_min) {
            dist_min=dist;
            j_min=j;
          }
        }
        for(size_t j=0;j<n_gauss;j++) {
          if (j==j_min) resp(j,i)=1.0;
          else resp(j,i)=0.0;
        }
        weights[j_min]+=1;
      }

      vals.resize(np);
      vector_set_all(np,vals,1.0);
      
      pdmg.set(nd_in,np,data,vals,
      
      return 0;
    }
    //@}
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /// The number of points
    size_t np;
    /// The number of dimensions of the inputs
    size_t nd_in;
    /// The copy of the data
    mat_t data;
    /// True if the data has been specified
    bool data_set;

    ubvector weights;
    ubmatrix resps;
    ubmatrix means;
    tensor3 covars;

    std::vector<prob_dens_mdim_gaussian> pdmg;
    
#endif
    
  };
    
}
    
#endif

