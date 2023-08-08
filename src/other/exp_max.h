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
#ifndef O2SCL_EXP_MAX_H
#define O2SCL_EXP_MAX_H

/** \file exp_max.h
    \brief File defining \ref o2scl::exp_max_gmm
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

      (very experimental and not working yet)
   */
  template<class data_mat_t=const_matrix_view_table<>,
           class gauss_vec_t=boost::numeric::ublas::vector<double>,
           class gauss_mat_t=boost::numeric::ublas::matrix<double>>
  class exp_max_gmm {

  public:

    /** \brief If true, call the error handler if the calculation
        fails
    */
    bool err_nonconv;

    /** \brief Get the underlying Gaussian mixture probability
        density
     */
    const prob_dens_mdim_gmm<> &get_gmm() {
      return pdmg;
    }
    
  protected:
    
    /// The underlying Gaussian mixture probability density
    prob_dens_mdim_gmm<gauss_vec_t,gauss_mat_t> pdmg;
    
    typedef boost::numeric::ublas::vector<double> internal_vec_t;
    typedef boost::numeric::ublas::matrix<double> internal_mat_t;
    
    /** \brief Use the expectation-maximization algorithm to 
        optimize the Gaussian mixture
     */
    int calc_internal(size_t n_gauss=1) {
      
      std::cout << "exp_max_gmm::calc_internal()" << std::endl;

      internal_mat_t last_means(n_gauss,nd_in);
      for(size_t k=0;k<n_gauss;k++) {
        const gauss_vec_t &peak=pdmg.pdmg[k].get_peak();
        if (verbose>0) {
          std::cout << "  Mean " << k+1 << " of " << n_gauss
                    << " before: " << k << std::endl;
          std::cout << "  ";
          vector_out(std::cout,peak,true);
        }
        for(size_t i=0;i<nd_in;i++) {
          last_means(k,i)=peak[i];
        }
      }
      
      bool done=false;
      for(int it=0;it<ntrial && done==false;it++) {

        // -----------------------------------------------------
        // Compute responsibilities (expectation step)
        
        for(size_t i=0;i<np;i++) {
          double total=0.0;
          gauss_vec_t data_i(nd_in);
          for(size_t j=0;j<nd_in;j++) {
            data_i[j]=data(i,j);
          }
          for(size_t k=0;k<n_gauss;k++) {
            total+=pdmg.pdmg[k].pdf(data_i)*pdmg.weights[k];
          }
          for(size_t k=0;k<n_gauss;k++) {
            resps(i,k)=pdmg.pdmg[k].pdf(data_i)*pdmg.weights[k]/total;
            if (verbose>2) {
              std::cout << "  resp: " << i << " " << k << " "
                        << resps(i,k) << std::endl;
            }
          }
        }
        
        // -----------------------------------------------------
        // Maximization step

        // Compute pdmg.weights
        for(size_t k=0;k<n_gauss;k++) {
          pdmg.weights[k]=0.0;
          for(size_t i=0;i<np;i++) {
            pdmg.weights[k]+=resps(i,k);
          }
          pdmg.weights[k]/=np;
          if (verbose>0) {
            std::cout << "  weight: " << k+1 << " of "
                      << n_gauss << ": " << pdmg.weights[k] << std::endl;
          }
        }
        
        // Compute means and covariances
        for(size_t k=0;k<n_gauss;k++) {
          matrix_column_gen<internal_mat_t> resp_k(resps,k);
          pdmg.pdmg[k].verbose=1;
          pdmg.pdmg[k].set_wgts(nd_in,np,data,resp_k);
        }
        
        // Compare means with previous values to test for convergence
        done=true;
        for(size_t k=0;k<n_gauss;k++) {
          const gauss_vec_t &peak=pdmg.pdmg[k].get_peak();
          if (verbose>0) {
            std::cout << "  Mean: " << k+1 << " of " << n_gauss << ":"
                      << std::endl;
            std::cout << "  ";
            vector_out(std::cout,peak,true);
          }
          for(size_t i=0;i<nd_in;i++) {
            if (done==true && (last_means(k,i)!=0.0 || peak[i]!=0.0)) {
              double rel_diff;
              if (last_means(k,i)==0.0) {
                rel_diff=abs(peak[i]-last_means(k,i))/abs(peak[i]);
              } else {
                rel_diff=abs(peak[i]-last_means(k,i))/abs(last_means(k,i));
              }
              if (tol_abs>0.0 && peak[i]>tol_abs) {
                if (verbose>2) {
                  std::cout << "Setting done to false because peak(k,i) "
                            << "at k,i=" << k << "," << i << " is "
                            << peak[i] << " and tol_abs is "
                            << tol_abs << std::endl;
                }
                done=false;
              }
              if (rel_diff>tol_rel) {
                if (verbose>2) {
                  std::cout << "Setting done to false because peak(k,i) "
                            << "at k,i=" << k << "," << i << " is "
                            << peak[i] << " last_mean is "
                            << last_means(k,i) << " and tol_rel is "
                            << tol_rel << std::endl;
                }
                done=false;
              }
            }
            last_means(k,i)=peak[i];
          }
        }

        if (verbose>0) {
          std::cout << "Iteration " << it+1 << " of " << ntrial << std::endl;
          if (verbose>1) {
            char ch;
            std::cin >> ch;
          }
        }
      }

      if (done==false) {
        O2SCL_CONV_RET("Did not converge in "
                       "exp_max_gmm::calc_internal().",
                       o2scl::exc_einval,err_nonconv);
      }
      
      return 0;
    }
    
  public:
    
    exp_max_gmm() {
      data_set=false;
      verbose=0;
      ntrial=100;
      tol_rel=1.0e-8;
      tol_abs=0.0;
      nd_in=0;
      np=0;
      err_nonconv=true;
    }
    
    /** \brief Verbosity parameter (default 0)
     */
    int verbose;
    
    /// Maximum number of iterations
    int ntrial;

    /// Relative tolerance (default \f$ 10^[-6} \f$)
    double tol_rel;

    /// Absolute tolerance (default \f$ 0 \f$)
    double tol_abs;
    
    /** \brief Base random number generator

        This is used to automatically generate initial means for the
        Gaussians if a user-specified guess is not provided.
     */
    rng<> r2;
    
    /** \brief Initialize the data
        The object \c vecs should be a matrix with a
        first index of size <tt>n_in</tt> and a second 
        index of size <tt>n_points</tt>. It may have be
        any type which allows the use of <tt>operator(,)</tt>
        and <tt>std::swap</tt>.
    */
    void set_data(size_t n_in, size_t n_points, data_mat_t &dat) {
      
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
                   "exp_max_gmm::set_data()",exc_efailed);
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
                  data_mat_t &dat) {
      n_points=np;
      n_in=nd_in;
      std::swap(data,dat);
      data_set=false;
      return;
    }
    //@}

    /// \name Compute the Gaussian mixture model
    //@{
    /** \brief Compute the Gaussian mixture model using random 
        guesses based on the data
     */
    int calc_auto(size_t n_gauss) {

      if (np==0 || nd_in==0) {
        O2SCL_ERR2("Data not set in ","exp_max_gmm::calc_auto().",
                   o2scl::exc_einval);
      }
      
      if (n_gauss==0) {
        O2SCL_ERR2("Cannot select zero gaussians in ",
                   "exp_max_gmm::calc_auto().",
                   o2scl::exc_einval);
      }
      
      pdmg.weights.resize(n_gauss);
      resps.resize(np,n_gauss);
      pdmg.pdmg.resize(n_gauss);
      
      // Randomly initialize the Gaussians and the pdmg.weights
      
      gauss_mat_t tcovar(nd_in,nd_in);
      gauss_vec_t tmean(nd_in);
      pdmg.pdmg[0].set_ret(nd_in,np,data,tmean,tcovar);
      
      if (verbose>0) {
        std::cout << "exp_max_gmm::calc_auto(): initial mean: " << std::endl;
        std::cout << "  ";
        vector_out(std::cout,tmean,true);
        std::cout << "exp_max_gmm::calc_auto(): initial covar: " << std::endl;
        matrix_out(std::cout,nd_in,nd_in,tcovar,"  ");
      }
      
      for(int k=n_gauss-1;k>=0;k--) {
        pdmg.pdmg[0](tmean);
        if (verbose>0) {
          std::cout << "exp_max_gmm::calc_auto(): random mean: "
                    << k << std::endl;
          std::cout << "  ";
          vector_out(std::cout,tmean,true);
        }
        pdmg.pdmg[k].set_covar(nd_in,tmean,tcovar);
      }
      for(size_t k=0;k<n_gauss;k++) {
        pdmg.weights[k]=0.5;
      }

      return calc_internal(n_gauss);
    }

    /** \brief Compute the Gaussian mixture model using the 
        specified initial pdmg.weights, means, and covariance matrices
    */
    template<class vec_t, class mat2_t, class ten_t>
    int calc_guess(size_t n_gauss, vec_t &wgts,
                      mat2_t &means, ten_t &covars) {
      
      if (np==0 || nd_in==0) {
        O2SCL_ERR2("Data not set in ","exp_max_gmm::calc_guess().",
                   o2scl::exc_einval);
      }
      
      if (n_gauss==0) {
        O2SCL_ERR2("Cannot select zero gaussians in ",
                   "exp_max_gmm::calc_guess().",
                   o2scl::exc_einval);
      }
      
      pdmg.weights.resize(n_gauss);
      resps.resize(np,n_gauss);
      pdmg.pdmg.resize(n_gauss);

      for(size_t k=0;k<n_gauss;k++) {
        pdmg.weights[k]=wgts[k];
        if (verbose>0) {
          std::cout << "User-specified pdmg.weights: "
                    << pdmg.weights[k] << std::endl;
        }
        internal_vec_t mean(nd_in);
        internal_mat_t covar(nd_in,nd_in);
        for(size_t j=0;j<nd_in;j++) {
          mean[j]=means(k,j);
        }
        if (verbose>0) {
          std::cout << "User-specified mean: ";
          vector_out(std::cout,mean,true);
        }
        for(size_t i=0;i<nd_in;i++) {
          for(size_t j=0;j<nd_in;j++) {
            covar(i,j)=covars.get(k,i,j);
          }
        }
        if (verbose>0) {
          std::cout << "User-specified covariance matrix: ";
          matrix_out(std::cout,nd_in,nd_in,covar);
        }
        pdmg.pdmg[k].set_covar(2,mean,covar);
      }

      return calc_internal(n_gauss);
    }
    //@}
    
  protected:
    
    /// The number of points
    size_t np;
    /// The number of dimensions of the inputs
    size_t nd_in;
    /// The copy of the data
    data_mat_t data;
    /// True if the data has been specified
    bool data_set;

    /// The responsibilities
    internal_mat_t resps;
    
  };
    
}
    
#endif

