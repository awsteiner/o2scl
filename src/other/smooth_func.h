/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2022, Andrew W. Steiner
  
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
/** \file smooth_func.h
    \brief File for ..
*/
#ifndef O2SCL_SMOOTH_FUNC_H
#define O2SCL_SMOOTH_FUNC_H

#include <iostream>
#include <vector>

#include <gsl/gsl_qrng.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_vector.h>

#include <o2scl/funct.h>
#include <o2scl/vector.h>

namespace o2scl {

  /** \brief Smooth a function by averaging in a neighborhood
      of points defined by a Sobol sequence
      
      \warning The function \ref o2scl::smooth_func::set_func() stores a
      pointer to the function specified by the user, so the user must
      ensure that this pointer is still valid when 
      \ref o2scl::smooth_func::operator()() is called.
      
      \future Move memory allocation outside of \ref
      o2scl::smooth_func::operator()() .
  */
  template<class vec_t, class func_t> class smooth_func {
    
  protected:
  
    /** \brief The pointer to the original function
     */
    func_t *f;
  
    /** \brief Step size defining the neighborhood (default 0.01)
     */
    std::vector<double> step;

    /** \brief Number of points in the Sobol sequence (default 40)
     */
    size_t N;
  
  public:

    smooth_func() {
      N=40;
      step.resize(1);
      step[0]=1.0e-2;
      f=0;
    }

    /** \brief Set the base function
     */
    void set_func(func_t &func) {
      f=&func;
      return;
    }

    /** \brief Set the number of points to use in the average

	If \c n_new is zero then the error handler will be called.
    */
    void set_n(size_t n_new) {
      if (n_new==0) {
	O2SCL_ERR2("Cannot call set_n() with argument 0 in ",
		   "smooth_func::set_n().",o2scl::exc_einval);
      }
      N=n_new;
      return;
    }
  
    /** \brief Set the stepsize
     */
    template<class vec2_t> void set_step(vec2_t &v) {
      if (v.size()==0) {
	O2SCL_ERR2("Sent an empty vector to ",
		   "smooth_func::set_step().",o2scl::exc_einval);
      }
      step.resize(v.size());
      o2scl::vector_copy(v.size(),v,step);
      return;
    }

    /** \brief Evaluate the smoothed function

	If the user-specified function returns a non-zero value for any
	point, then that contribution to the average is ignored. This
	function will return a non-zero value if the user-specified
	function returns a non-zero value for all of the points.
    */
    int operator()(size_t nv, const vec_t &x, vec_t &y) {

      if (f==0) {
	O2SCL_ERR2("Function not set in ",
		   "smooth_func::operator().",o2scl::exc_einval);
      }
    
      // Allocate the Sobol object
      gsl_qrng *gq=gsl_qrng_alloc(gsl_qrng_sobol,nv);
    
      std::vector<double> v(nv);
      vec_t x2(nv), y2(nv);

      for(size_t j=0;j<nv;j++) {
	y(j)=0.0;
      }

      int count=0;
      for(size_t i=0;i<N;i++) {

	// Create the new x point
	gsl_qrng_get(gq,&(v[0]));
	for(size_t j=0;j<nv;j++) {
	  x2[j]=(1.0+(v[j]*2.0-1.0)*step[j%step.size()])*x[j];
	}

	// Evaluate the function
	int fret=(*f)(nv,x2,y2);

	// Add the y value to the running total
	if (fret==0) {
	  for(size_t j=0;j<nv;j++) {
	    y[j]+=y2[j];
	  }
	  count++;
	}
      }

      if (count==0) {
	return o2scl::exc_efailed;
      }
    
      // Compute the average from the running total
      for(size_t j=0;j<nv;j++) {
	y[j]/=((double)count);
      }
    
      gsl_qrng_free(gq);

      return 0;
    }
  
  };

  /** \brief Apply a Gaussian filter to a function
   */
  template<class func_t=funct> class gauss_filter {
    
  protected:
  
    /** \brief The pointer to the original function
     */
    func_t *f;
  
    /** \brief Number of points (default 5)
     */
    size_t K;

    /** \brief GSL workspace
     */
    gsl_filter_gaussian_workspace *w;

    /** \brief Windowing parameter
     */
    double alpha;

    /** \brief Vector storage for kernel
     */
    gsl_vector *v;
    
  public:

    gauss_filter() {
      K=5;
      f=0;
      w=gsl_filter_gaussian_alloc(K);
      v=gsl_vector_alloc(K);
      alpha=1.0;
      h_rel=1.0e-14;
      h_abs=0.0;
      lower_limit=0.0;
      upper_limit=0.0;
    }

    /// Relative stepsize (default \f$ 10^{-14} \f$
    double h_rel;
    
    /// Absolute stepsize (default \f$ 0 \f$)
    double h_abs;

    /// Lower limit
    double lower_limit;
    
    /// Upper limit
    double upper_limit;
    
    virtual ~gauss_filter() {
      gsl_filter_gaussian_free(w);
      gsl_vector_free(v);
    }
    
    /** \brief Set the base function
     */
    void set_func(func_t &func) {
      f=&func;
      return;
    }

    /** \brief Set the window parameter
     */
    void set_alpha(double alpha_new) {
      alpha=alpha_new;
      return;
    }
    
    /** \brief Set the number of points to use in the average

	If \c n_new is zero then the error handler will be called.
    */
    void set_K(size_t K_new) {
      if (K_new==0) {
	O2SCL_ERR2("Cannot call set_n() with argument 0 in ",
		   "gauss_filter::set_n().",o2scl::exc_einval);
      }
      gsl_filter_gaussian_free(w);
      gsl_vector_free(v);
      K=K_new;
      w=gsl_filter_gaussian_alloc(K);
      v=gsl_vector_alloc(K);
      return;
    }
  
    /** \brief Evaluate the function
    */
    double operator()(double x) {
      
      if (f==0) {
	O2SCL_ERR2("Function not set in ",
		   "gauss_filter::operator().",o2scl::exc_einval);
      }

      double h=0.0;
      if (h_rel>0.0) {
        if (h_abs>0.0) {
          h=x*h_rel+h_abs;
        } else {
          h=x*h_rel;
        }
      } else {
        if (h_abs>0.0) {
          h=h_abs;
        } else {
          O2SCL_ERR("Step not set in gauss_filter::operator().",
                    o2scl::exc_einval);
        }
      }
      
      gsl_filter_gaussian_kernel(alpha,0,1,v);
      
      double y=0.0;

      // If the smallest or largest value is beyond the lower and
      // upper limits, then adjust the sum accordingly
      
      if (lower_limit<upper_limit && (x-h<lower_limit || x+h>upper_limit)) {
        
        double w_sum=0.0;
        for(int j=0;j<((int)K);j++) {
          double x2=x+h*((double)(j-((int)K)/2))/((double)(K/2));
          if (x2>=lower_limit && x2<=upper_limit) {
            w_sum+=gsl_vector_get(v,j);
            y+=(*f)(x2)*gsl_vector_get(v,j);
          }
        }
        y/=w_sum;
        
      } else {
        
        for(int j=0;j<((int)K);j++) {
          double x2=x+h*((double)(j-((int)K)/2))/((double)(K/2));
          y+=(*f)(x2)*gsl_vector_get(v,j);
        }
        
      }

      return y;
    }
  
  };
  
}

#endif
