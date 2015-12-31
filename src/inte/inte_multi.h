/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#ifndef O2SCL_MULTI_INTE_H
#define O2SCL_MULTI_INTE_H

/** \file inte_multi.h
    \brief File defining \ref o2scl::inte_multi
*/
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/inte.h>
#include <o2scl/multi_funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Multi-dimensional integration over a hypercube 
      [abstract base]
 
      Multi-dimensional integration over a region defined by constant
      limits. For more general regions of integration, use children of
      the class \ref inte_gen.
  */
  template<class func_t=multi_funct11, 
    class vec_t=boost::numeric::ublas::vector<double> > class inte_multi {
      
      public:
      
      inte_multi() {
	tol_rel=1.0e-8;
	verbose=0;
	interror=0.0;
	err_nonconv=true;
      }
      
      virtual ~inte_multi() {}
      
      /// If true, call the error handler if the routine does not "converge"
      bool err_nonconv;
      
      /** \brief Verbosity
       */
      int verbose;
  
      /** \brief The maximum "uncertainty" in the value of the integral
	  (default \f$ 10^{-8} \f$).
       */
      double tol_rel;
  
      /** \brief Integrate function \c func over the hypercube from
	  \f$ x_i=a_i \f$ to \f$ x_i=b_i \f$ for
          \f$ 0<i< \f$ ndim-1
	  
	  \comment 
	  (I had some problems with Doxygen 
	  using \mathrm in the formula above.)
	  \endcomment
      */
      virtual double minteg(func_t &func, size_t ndim, const vec_t &a, 
			    const vec_t &b) {
	double err, res;
	minteg_err(func,ndim,a,b,res,err);
	return res;
      }
      
      /** \brief Integrate function \c func over the hypercube from
	  \f$ x_i=a_i \f$ to \f$ x_i=b_i \f$ for
          \f$ 0<i< \f$ ndim-1
      */
      virtual int minteg_err(func_t &func, size_t ndim, const vec_t &a, 
			     const vec_t &b, double &res, double &err)=0;
      
      /** \brief Return the error in the result from the last call to 
	  minteg() or minteg_err()
	  
	  This will quietly return zero if no integrations have been
	  performed.
      */
      double get_error() { return interror; }
  
      /// Return string denoting type ("inte_multi")
      const char *type() { return "inte_multi"; }
  
#ifndef DOXYGEN_INTERNAL
  
      protected:
  
      /// The uncertainty for the last integration computation
      double interror;
  
#endif
  
    };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif





  
