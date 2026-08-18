/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef O2SCL_INTERP_SMART_H
#define O2SCL_INTERP_SMART_H

/*
  AWS 12/14/12: Completely deprecating this for now until it's fixed.
*/
#ifdef O2SCL_NEVER_DEFINED

/** \file interp_smart.h
    \brief File for "smart" interpolation routines

    In addition to the smart interpolation routines, this file
    contains the template functions \ref vector_find_level(), \ref
    vector_integ_linear() and \ref vector_invert_enclosed_sum().
*/

#include <o2scl/interp.h>
#include <o2scl/vector.h>
#include <o2scl/ovector_rev_tlate.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief Smart interpolation class with pre-specified vectors
      
      \todo This class needs updating with the new base interpolation 
      classes

      This class can semi-intelligently handle arrays which are
      not-well formed outside the interpolation region. In particular,
      if an initial interpolation or derivative calculation fails, the
      arrays are searched for the largest neighborhood around the
      point \c x for which an interpolation or differentiation will
      likely produce a finite result.

      \future Properly implement handling of non-monotonic regions in
      the derivative functions as well as the interpolation function.
  */
  template<class vec_t, class svec_t, class alloc_vec_t, class alloc_t>
    class interp_smart {

#ifndef DOXYGEN_INTERNAL
    
  protected:

    /// If true, then \ref sx and \ref sy have been allocated
    bool sxalloc;
    /// Storage for internally created subvector
    svec_t *sx;
    /// Storage for internally created subvector
    svec_t *sy;

    /// Pointer to base interpolation object
    interp_base<vec_t> *rit1;
    /// Pointer to base interpolation object
    interp_base<svec_t> *rit2;
    
    /// Memory allocator for objects of type \c alloc_vec_t
    alloc_t ao;

    /// True if the user-specified x vector is increasing
    bool inc;

    /// Pointer to user-specified vector
    const vec_t *lx;
    /// Pointer to user-specified vector
    const vec_t *ly;

    /// Size of user-specifed vector
    size_t ln;

    /// Interpolation type
    size_t itype;
    
#endif

  public:

    /// Create with base interpolation objects \c it and \c rit
    interp_smart(size_t n, const vec_t &x, const vec_t &y,
		     size_t interp_type=1) {

      sx=0;
      sy=0;
      sxalloc=false;

      if (interp_type==itp_linear) {
	rit1=new interp_linear<vec_t>;
	rit2=new interp_linear<svec_t>;
      } else if (interp_type==itp_cspline) {
	rit1=new interp_cspline<vec_t>;
	rit2=new interp_cspline<svec_t>;
      } else if (interp_type==itp_cspline_peri) {
	rit1=new interp_cspline_peri<vec_t>;
	rit2=new interp_cspline_peri<svec_t>;
      } else if (interp_type==itp_akima) {
	rit1=new interp_akima<vec_t>;
	rit2=new interp_akima<svec_t>;
      } else if (interp_type==itp_akima_peri) {
	rit1=new interp_akima_peri<vec_t>;
	rit2=new interp_akima_peri<svec_t>;
      } else {
	O2SCL_ERR2("Invalid interpolation type in ",
		   "interp_smart::interp_smart().",gsl_einval);
      }

      ln=0;
	
      if (n<rit1->min_size || n<rit2->min_size) {
	O2SCL_ERR("Vector size too small in o2scl_interp_vec().",gsl_edom);
      } else {

	rit1->allocate(n);
	rit1->init(x,y,n);
	ln=n;
	lx=&x;
	ly=&y;
	inc=true;

      } 
	
    }

    virtual ~interp_smart() {
      if (sxalloc) {
	delete sx;
	delete sy;
      }
      
      if (ln>0) {
	rit1->free();
	rit2->free();
      }
      
      /*
	We don't want to use "bim1->free_interp(rit1)" here instead
	of "delete rit1" because the pointer "bim1" may refer to 
	a user-defined interpolation manager, and it may have gone
	out of scope earlier.
      */
      delete rit1;
      delete rit2;
    }

    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double operator()(const double x0) {
      return interp(x0);
    }

    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double interp(const double x0) {
      
      double ret=0.0;
      int err;
      if (ln>0) {
	
	rit1->interp(*lx,*ly,ln,x0,ret);
	size_t nn=0;
	int retfs=find_inc_subset(x0,ln,*lx,*ly,nn);
	
	if (retfs==0 && nn>1 && nn>=rit2->min_size) {
	  rit2->allocate(nn);
	  rit2->init(*sx,*sy,nn);
	  rit2->interp(*sx,*sy,nn,x0,ret);
	  rit2->free();
	} else {
	  O2SCL_ERR2("Interpolation failed in ",
		     "interp_smart::interp().",gsl_efailed);
	  return 0.0;
	}
	
      }
      
      return ret;
    }
    
    /// Give the value of the derivative \f$ y^{prime}(x=x_0) \f$ .
    virtual double deriv(const double x0) {
      double ret=0.0;
      if (ln>0) {
	rit1->deriv(*lx,*ly,ln,x0,ret);
      }
      return ret;
    }		      
    
    /** \brief Give the value of the second derivative  
	\f$ y^{prime \prime}(x=x_0) \f$ .
    */
    virtual double deriv2(const double x0) {
      double ret=0.0;
      if (ln>0) {
	rit1->deriv2(*lx,*ly,ln,x0,ret);
      }
      return ret;
    }		      
    
    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(const double x1, const double x2) {
      double ret=0.0;
      if (ln>0) {
	rit1->integ(*lx,*ly,ln,x1,x2,ret);
      }
      return ret;
    }		      

#ifndef DOXYGEN_INTERNAL

  protected:

    /// A lookup function for generic vectors
    size_t local_lookup(size_t n, const vec_t &x, double x0) {
      size_t row=0, i=0;
      while(!gsl_finite(x[i]) && i<n-1) i++;
      if (i==n-1) {
	return 0;
      }
      double best=x[i], bdiff=fabs(x[i]-x0);
      for(;i<n;i++) {
	if (gsl_finite(x[i]) && fabs(x[i]-x0)<bdiff) {
	  row=i;
	  best=x[i];
	  bdiff=fabs(x[i]-x0);
	}
      }
      return row;
    }

    /** \brief Try to find the largest monotonically increasing and
	finite region around the desired location

	This function looks through the vector \c x near the element
	closest to \c x0 to find the largest possible monotonic
	region. If it succeeds, it returns \ref gsl_success, and if it
	fails, it returns \ref gsl_efailed. It does not call the error
	handler.

	\todo Variables \c row and \c row appear to be unused? (2/17/11)
    */
    int find_inc_subset(const double x0, size_t sz, const vec_t &x, 
			const vec_t &y, size_t &nsz) {

      size_t row=local_lookup(sz,x,x0), row2=row++;
      size_t left=row, right=row;
      
      if (row2>sz) row2=row-1;

      if (!gsl_finite(x[row]) || !gsl_finite(y[row])) {
	return gsl_efailed;
      }
  
      // Increasing case
      
      while(left>0 && x[left-1]<x[left] && 
	    gsl_finite(x[left-1]) && gsl_finite(y[left-1])) left--;
      while(right<sz-1 && x[right]<x[right+1] &&
	    gsl_finite(x[right+1]) && gsl_finite(y[right+1])) right++;
      
      nsz=right-left+1;
      
      if (sxalloc) {
	delete sx;
	delete sy;
      } 
      sx=new svec_t(x,left,nsz);
      sy=new svec_t(y,left,nsz);
      sxalloc=true;
      
      return 0;
    }

#endif

  };

  // sm_interp_vec typedef
  typedef interp_smart<ovector_const_view,ovector_const_subvector,ovector,
    ovector_alloc> interp_sm;
  
  /** \brief A specialization of interp_smart for C-style double arrays
   */
  template<class arr_t> class interp_sma : 
  public interp_smart<arr_t,array_const_subvector,arr_t,
    array_alloc<arr_t> >
      {
	
      public:
    
    
	/** \brief Create an interpolation object with user-specified 
	    interpolation types
	*/
      interp_sma(size_t n, const arr_t &x, const arr_t &y,
		     size_t interp_type) :
	interp_smart<arr_t,array_const_subvector,arr_t,
	  array_alloc<arr_t> >(n,x,y,interp_type) {
	}
	
	/// Create an interpolation object with the default interpolation types
      interp_sma(size_t n, const arr_t &x, const arr_t &y) :    
	interp_smart<arr_t,array_const_subvector,arr_t,
	  array_alloc<arr_t> >(n,x,y) {
	}
	
      };

#ifndef DOXYGENP
}
#endif

#endif

#endif
