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
#include <o2scl/prob_dens_func.h>

namespace o2scl {

  /** \brief Expectation maximization for a Gaussian mixture model
   */
  template<class mat_t=const_matrix_view_table<> > class exp_max_gmm {
    
  protected:
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    
    exp_max_gmm() {
      data_set=false;
      verbose=1;
    }
    
    /** \brief Verbosity parameter (default 0)
     */
    int verbose;
    
    /** \brief Initialize the data
        The object \c vecs should be a matrix with a
        first index of size <tt>n_in</tt> and a second 
        index of size <tt>n_points</tt>. It may have be
        any type which allows the use of <tt>operator(,)</tt>
        and <tt>std::swap</tt>.
    */
    void set_data(size_t n_in, size_t n_points, mat_t &dat) {
      
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
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
      return;
    }
    //@}

    /// \name Evaluate interpolation
    //@{
    /** \brief Perform the interpolation over the first function
     */
    int compute_auto(size_t n_gauss=1) {
      
      if (n_gauss==0) {
        O2SCL_ERR("Cannot select zero gaussians in compute().",
                  o2scl::exc_einval);
      }
      
      weights.resize(n_gauss);
      resps.resize(np,n_gauss);
      pdmg.resize(n_gauss);
      
      // Randomly initialize the Gaussians and the weights
      
      ubmatrix tcovar(nd_in,nd_in);
      ubvector tmean(nd_in);
      pdmg[0].set_ret(nd_in,np,data,tmean,tcovar);
      
      if (verbose>0) {
        std::cout << "tmean: " << std::endl;
        vector_out(std::cout,tmean,true);
        std::cout << "tcovar: " << std::endl;
        matrix_out(std::cout,nd_in,nd_in,tcovar);
      }
      
      for(int k=n_gauss-1;k>=0;k--) {
        pdmg[0](tmean);
        if (verbose>0) {
          std::cout << "tmean: " << k << std::endl;
          vector_out(std::cout,tmean,true);
        }
        pdmg[k].set_covar(nd_in,tmean,tcovar);
      }
      for(size_t k=0;k<n_gauss;k++) {
        weights[k]=0.5;
      }

      compute(n_gauss);
      
      return 0;
    }

    /** \brief Desc
     */
    template<class ten_t, class mat2_t, class vec_t>
    int compute_guess(size_t n_gauss, vec_t &wgts,
                      mat2_t &means, ten_t &covars) {
      
      if (n_gauss==0) {
        O2SCL_ERR("Cannot select zero gaussians in compute().",
                  o2scl::exc_einval);
      }
      
      weights.resize(n_gauss);
      resps.resize(np,n_gauss);
      pdmg.resize(n_gauss);

      for(size_t k=0;k<n_gauss;k++) {
        weights[k]=wgts[k];
        ubvector mean(nd_in);
        ubmatrix covar(nd_in,nd_in);
        for(size_t j=0;j<nd_in;j++) {
          mean[j]=means(k,j);
        }
        for(size_t i=0;i<nd_in;i++) {
          for(size_t j=0;j<nd_in;j++) {
            covar(i,j)=covars(k,i,j);
          }
        }
      }

      compute(n_gauss);
      
      return 0;
    }

    /** \brief Desc
     */
    int compute(size_t n_gauss=1) {
      
      for(size_t it=0;it<20;it++) {
        
        // Compute responsibilities (expectation step)
        
        double total=0.0;
        for(size_t i=0;i<np;i++) {
          double total=0.0;
          ubvector data_i(nd_in);
          for(size_t j=0;j<nd_in;j++) {
            data_i[j]=data(i,j);
          }
          for(size_t k=0;k<n_gauss;k++) {
            total+=pdmg[k].pdf(data_i)*weights[k];
          }
          for(size_t k=0;k<n_gauss;k++) {
            resps(i,k)=pdmg[k].pdf(data_i)*weights[k]/total;
            if (verbose>1) {
              std::cout << "resp: " << i << " " << k << " "
                        << resps(i,k) << std::endl;
            }
          }
        }
        
        // Maximization step

        // Compute weights
        for(size_t k=0;k<n_gauss;k++) {
          weights[k]=0.0;
          for(size_t i=0;i<nd_in;i++) {
            weights[k]+=resps(i,k);
          }
          weights[k]/=np;
          if (verbose>0) {
            std::cout << "weights: " << k << " " << weights[k] << std::endl;
          }
        }
        
        // Compute means and covariances
        for(size_t k=0;k<n_gauss;k++) {
          matrix_column_gen<ubmatrix> resp_k(resps,k);
          pdmg[k].set_wgts(nd_in,np,data,resp_k);
          if (verbose>0) {
            std::cout << "Means: " << k << std::endl;
            const ubvector &peak=pdmg[k].get_peak();
            vector_out(std::cout,peak,true);
          }
        }

        std::cout << "Iteration " << it << " of 20" << std::endl;
        char ch;
        std::cin >> ch;
      }
      
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

    /// Weights 
    ubvector weights;
    /// Responsibilities
    ubmatrix resps;
    /// The gaussians
    std::vector<o2scl::prob_dens_mdim_gaussian<>> pdmg;
    
#endif
    
  };
    
}
    
#endif

