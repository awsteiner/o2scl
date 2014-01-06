/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
/* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Gerard
 * Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 * 
 */
#ifndef O2SCL_SERIES_ACC_H
#define O2SCL_SERIES_ACC_H

/** \file series_acc.h
    \brief File defining \ref o2scl::series_acc
*/

#include <iostream>
#include <cmath>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_machine.h>
#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Series acceleration by Levin u-transform (GSL)

      Given an array of terms in a sum, this attempts to evaluate the
      entire sum with an estimate of the error.

      \future Move the workspaces to classes?
      \future Create an example
  */
  class series_acc {

  public:
    
    /** \brief \c size is the number of terms in the series
     */
    series_acc(size_t size=0);

    virtual ~series_acc();

    /** \brief Return the accelerated sum of the series with
	a simple error estimate

	The input vector \c x should be an array with \c n values
	from <tt>x[0]</tt> to <tt>x[n-1]</tt>.
    */
    template<class vec_t>
      double series_accel(size_t na, vec_t &array, double &abserr_trunc) {

      if (na!=size) set_size(na);

      double sum_accel;
      size_t min_terms=0;
      size_t max_terms=size-1;

      if (size == 0) {

	sum_accel=0.0;
	abserr_trunc=0.0;
	wt->sum_plain=0.0;
	wt->terms_used=0;
	return sum_accel;

      } else if (size == 1) {

	sum_accel=array[0];
	abserr_trunc=GSL_POSINF;
	wt->sum_plain=array[0];
	wt->terms_used=1;
	return sum_accel;

      } else {

	const double SMALL=0.01;
	const size_t nmax=GSL_MAX(max_terms,size)-1;
	double trunc_n=0.0, trunc_nm1=0.0;
	double actual_trunc_n=0.0, actual_trunc_nm1=0.0;
	double result_n=0.0, result_nm1=0.0;
	size_t n;
	int better=0;
	int before=0;
	int converging=0;
	double least_trunc=GSL_DBL_MAX;
	double result_least_trunc;
	
	/* Calculate specified minimum number of terms. No convergence
	   tests are made, and no truncation information is stored. */
	
	for (n=0; n < min_terms; n++) {
          const double t=array[n];
	  
          result_nm1=result_n;
          levin_utrunc_step(t,n,result_n);
        }

	/* Assume the result after the minimum calculation is the best. */

	result_least_trunc=result_n;

	/* Calculate up to maximum number of terms. Check truncation
	   condition. */

	for (; n <= nmax; n++) {
          const double t=array[n];

          result_nm1=result_n;
          levin_utrunc_step(t,n,result_n);

          /* Compute the truncation error directly */

          actual_trunc_nm1=actual_trunc_n;
          actual_trunc_n=fabs(result_n-result_nm1);

          /* Average results to make a more reliable estimate of the
             real truncation error */

          trunc_nm1=trunc_n;
          trunc_n=0.5*(actual_trunc_n+actual_trunc_nm1);

          /* Determine if we are in the convergence region. */

          better=(trunc_n < trunc_nm1 || trunc_n < SMALL*fabs(result_n));
          converging=converging || (better && before);
          before=better;

          if (converging) {
	    if (trunc_n < least_trunc) {
	      /* Found a low truncation point in the convergence
		 region. Save it. */

	      least_trunc=trunc_n;
	      result_least_trunc=result_n;
	    }

	    if (fabs(trunc_n/result_n) < 10.0*GSL_MACH_EPS)
	      break;
	  }
        }

	if (converging) {

          /* Stopped in the convergence region. Return result and
             error estimate. */
	  sum_accel=result_least_trunc;
	  abserr_trunc=least_trunc;
          wt->terms_used=n;
          return sum_accel;

	} else {

          /* Never reached the convergence region. Use the last
             calculated values. */
	  sum_accel=result_n;
	  abserr_trunc=trunc_n;
          wt->terms_used=n;
          return sum_accel;
        }
      }
    

      return sum_accel;
    }

    /** \brief Return the accelerated sum of the series with
	an accurate error estimate

	The input vector \c x should be an array with \c n values
	from <tt>x[0]</tt> to <tt>x[n-1]</tt>.
    */
    template<class vec_t>
      double series_accel_err(size_t na, vec_t &array, double &abserr) {

      if (na!=size) set_size(na);

      double sum_accel;
      size_t min_terms=0;
      size_t max_terms=size-1;

      /* Ignore any trailing zeros in the array */
      size_t size2=size;
      
      while (size2 > 0 && array[size2-1]==0) {
	size2--;
      }
      
      if (size2==0) {

	sum_accel=0.0;
	abserr=0.0;
	w->sum_plain=0.0;
	w->terms_used=0;
	return sum_accel;

      } else if (size2==1) {

	sum_accel=array[0];
	abserr=0.0;
	w->sum_plain=array[0];
	w->terms_used=1;
	return sum_accel;

      } else {

	const double SMALL=0.01;
	const size_t nmax=GSL_MAX(max_terms,size)-1;
	double noise_n=0.0,noise_nm1=0.0;
	double trunc_n=0.0, trunc_nm1=0.0;
	double actual_trunc_n=0.0, actual_trunc_nm1=0.0;
	double result_n=0.0, result_nm1=0.0;
	double variance=0;
	size_t n;
	unsigned int i;
	int better=0;
	int before=0;
	int converging=0;
	double least_trunc=GSL_DBL_MAX;
	double least_trunc_noise=GSL_DBL_MAX;
	double least_trunc_result;

	/* Calculate specified minimum number of terms.  No convergence
	   tests are made, and no truncation information is stored.  */

	for (n=0; n < min_terms; n++) {
	  const double t=array[n];
	  result_nm1=result_n;
	  levin_u_step(t,n,nmax,result_n);
	}

	least_trunc_result=result_n;

	variance=0;
	for (i=0; i < n; i++) {
	  double dn=w->dsum[i]*GSL_MACH_EPS*array[i];
	  variance += dn*dn;
	}
	noise_n=std::sqrt(variance);

	/* Calculate up to maximum number of terms.  Check truncation
	   condition.  */

	for (; n <= nmax; n++) {
	  const double t=array[n];

	  result_nm1=result_n;
	  levin_u_step(t,n,nmax,result_n);

	  /* Compute the truncation error directly */

	  actual_trunc_nm1=actual_trunc_n;
	  actual_trunc_n=fabs(result_n-result_nm1);

	  /* Average results to make a more reliable estimate of the
	     real truncation error */

	  trunc_nm1=trunc_n;
	  trunc_n=0.5*(actual_trunc_n+actual_trunc_nm1);

	  noise_nm1=noise_n;
	  variance=0;

	  for (i=0; i <= n; i++) {
	    double dn=w->dsum[i]*GSL_MACH_EPS*array[i];
	    variance += dn*dn;
	  }
	  
	  noise_n=std::sqrt(variance);

	  /* Determine if we are in the convergence region.  */

	  better=(trunc_n < trunc_nm1 || 
		  trunc_n < SMALL*fabs(result_n));
	  converging=converging || (better && before);
	  before=better;

	  if (converging) {
	    if (trunc_n < least_trunc) {
	      /* Found a low truncation point in the convergence
		 region. Save it. */

	      least_trunc_result=result_n;
	      least_trunc=trunc_n;
	      least_trunc_noise=noise_n;
	    }

	    if (noise_n > trunc_n/3.0) {
	      break;
	    }

	    if (trunc_n < 10.0*GSL_MACH_EPS*fabs(result_n)) {
	      break;
	    }
	  }

	}

	if (converging) {

	  /* Stopped in the convergence region.  Return result and
	     error estimate.  */
	  sum_accel=least_trunc_result;
	  abserr=GSL_MAX_DBL(least_trunc,least_trunc_noise);
	  w->terms_used=n;
	  return sum_accel;

	} else {

	  /* Never reached the convergence region.  Use the last
	     calculated values.  */
	  sum_accel=result_n;
	  abserr=GSL_MAX_DBL(trunc_n,noise_n);
	  w->terms_used=n;
	  return sum_accel;
	}
      }

      return sum_accel;
    }

    /** \brief Set the number of terms */
    void set_size(size_t new_size);

#ifndef DOXYGEN_NO_O2NS

  protected:
    
    /** \brief An internal function reducing two matrix indices, i and j,
	to index of a single array
    */
    size_t series_index(size_t i, size_t j, size_t nmax) {
      return i*(nmax+1)+j;
    }

    /** \brief Perform a step
     */
    int levin_u_step(const double term, const size_t n, 
		     const size_t nmax, double &sum_accel);    

    /** \brief Perform a step
     */
    int levin_utrunc_step(const double term, const size_t n, 
			  double &sum_accel);
    
    /// The GSL workspace
    gsl_sum_levin_u_workspace *w;

    /// The GSL workspace
    gsl_sum_levin_utrunc_workspace *wt;

    /// The workspace size
    size_t size;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
