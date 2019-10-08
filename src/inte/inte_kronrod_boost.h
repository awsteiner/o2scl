/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019, Andrew W. Steiner
  
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
#ifndef O2SCL_INTE_KRONROD_BOOST_H
#define O2SCL_INTE_KRONROD_BOOST_H

/** \file inte_kronrod_boost.h
    \brief File defining \ref o2scl::inte_kronrod_boost
*/

#include <cmath>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <o2scl/inte.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Gauss-Kronrod integration class (Boost)

      The rule parameter should be either 15, 31, 41, 51, or 61. 
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .

      \future Figure out what to do with L1norm. The boost
      documentation claims that "the error estimates provided by the
      routine are woefully pessimistic" and the integral appears to be
      correct, but the boost documentation also says "if there is a
      significant difference between this [the L1 norm] and the
      returned value, then the result is likely to be
      ill-conditioned". It would be nice to test L1 norm in some
      reasonable way.
   */
  template<class func_t=funct, size_t rule=15, class fp_t=double>
  class inte_kronrod_boost :
    public inte<func_t,fp_t> {
    
  protected:

    /// Maximum depth
    size_t max_depth;

  /// L1 norm
  fp_t L1norm;
  
    public:

    inte_kronrod_boost() {
      max_depth=15;
    }
  
    virtual ~inte_kronrod_boost() {
    }

  void set_max_depth(size_t md) {
    max_depth=md;
    return;
  }
    
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			  fp_t &res, fp_t &err) {
      res=boost::math::quadrature::gauss_kronrod<fp_t,rule>::integrate
	(func,a,b,max_depth,this->tol_rel,&err,&L1norm);
      if (err>this->tol_rel) {
	O2SCL_ERR2("Failed to achieve tolerance in ",
		   "inte_kronrod_boost::integ_err().",o2scl::exc_efailed);
      }
      return 0;
    }
  
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
