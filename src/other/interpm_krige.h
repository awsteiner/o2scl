/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017, Andrew W. Steiner
  
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
#ifndef O2SCL_INTERPM_KRIGE_H
#define O2SCL_INTERPM_KRIGE_H

/** \file interpm_krige.h
    \brief File defining \ref o2scl::interpm_krige
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

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multi-dimensional interpolation by kriging
  */
  template<class vec_t> class interpm_krige {

  protected:

    /** \brief Inverse covariance matrix times function vector
     */
    std::vector<ubvector> Kinvf;

    /** \brief Pointer to user-specified covariance function
     */
    covar_func_t *f;
    
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    
    interpm_krige() {
      data_set=false;
      verbose=0;
    }

    /** \brief Verbosity parameter (default 0)
     */
    int verbose;

    /** \brief Initialize the data for the interpolation

	The object \c vecs should be a vector (of size <tt>n_in+n_out</tt>)
	of vectors (all of size <tt>n_points</tt>). It may have be
	any time which allows the use of <tt>std::swap</tt> for
	each vector in the list. 
    */
    template<class vec_vec_t, class vec_vec2_t, class vec_vec3_t>
      void set_data_noise(size_t n_in, size_t n_out, size_t n_points,
			  vec_vec_t &x, vec_vec2_t &y, vec_vec3_t &vars) {

      if (n_points<3) {
	O2SCL_ERR2("Must provide at least three points in ",
		   "interpm_krige::set_data()",exc_efailed);
      }
      if (n_in<1) {
	O2SCL_ERR2("Must provide at least one input column in ",
		   "interpm_krige::set_data()",exc_efailed);
      }
      if (n_out<1) {
	O2SCL_ERR2("Must provide at least one output column in ",
		   "interpm_krige::set_data()",exc_efailed);
      }
      np=n_points;
      nd_in=n_in;
      nd_out=n_out;
      ptrs_x.resize(n_in);
      for(size_t i=0;i<n_in;i++) {
	std::swap(ptrs_x[i],x[i]);
      }
      data_set=true;

      Kinvf.resize(n_out);

      // Store pointer to covariance function
      f=&fcovar;

      // Loop over all output functions
      for(size_t iout=0;iout<n_out;iout++) {
	
	// Construct the KXX matrix
	ubmatrix KXX(n_points,n_points);
	for(size_t irow=0;irow<n_points;irow++) {
	  for(size_t icol=0;icol<n_points;icol++) {
	    if (irow>icol) {
	      KXX(irow,icol)=KXX(icol,irow);
	    } else if (irow==icol) {
	      KXX(irow,icol)=fcovar(x[irow],x[icol])+var[iout][irow];
	    } else {
	      KXX(irow,icol)=fcovar(x[irow],x[icol]);
	    }
	  }
	}
	
	// Construct the inverse of KXX
	ubmatrix inv_KXX(n_points,n_points);
	o2scl::permutation p;
	int signum;
	o2scl_linalg::LU_decomp(n_points,KXX,p,signum);
	o2scl_linalg::LU_invert<ubmatrix,ubmatrix,ubmatrix_column>
	  (n_points,KXX,p,inv_KXX);

	// Inverse covariance matrix times function vector
	Kinvf[iout].resize(n_dim);
	boost::numeric::ublas::axpy_prod(inv_KXX,y[iout],Kinvf[iout],true);
	
      }
      
      return;
    }

    /** \brief Perform the interpolation over the first function
     */
    template<class vec2_t, class vec3_t>
      double eval(const vec2_t &x, vec3_t &y) const {
    
      if (data_set==false) {
	O2SCL_ERR("Data not set in interpm_krige::eval().",
		  exc_einval);
      }

      y.resize(nd_out);
      for(size_t iout=0;iout<nd_out;iout++) {
	y[iout]=0.0;
	for(size_t iin=0;iin<nd_in;iin++) {
	  y[iout]+=(*f)(x,ptrs_x[iin])*Kinvf[iout][iin];
	}
      }

      return ret;
      
    }
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /// The number of points
    size_t np;
    /// The number of dimensions of the inputs
    size_t nd_in;
    /// The number of dimensions of the outputs
    size_t nd_out;
    /// A vector of pointers holding the data
    std::vector<vec_t> ptrs_x;
    /// True if the data has been specified
    bool data_set;
    /// Number of points to include in each interpolation (default 3)
    size_t order;
    
#endif
    
  };
    
#ifndef DOXYGEN_NO_O2NS
}
#endif
    
#endif



