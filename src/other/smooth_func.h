/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2018, Andrew W. Steiner
  
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
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

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

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
