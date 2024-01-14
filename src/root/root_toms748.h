/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2024, Andrew W. Steiner
  
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
/** \file root_toms748.h
    \brief File defining \ref o2scl::root_toms748
*/
#ifndef O2SCL_ROOT_TOMS748_H
#define O2SCL_ROOT_TOMS748_H

#include <string>

#include <boost/math/tools/roots.hpp>

#include <o2scl/root.h>

namespace o2scl {

  /** \brief Convergence test similar to
      <tt>gsl_root_test_interval()</tt> for \ref root_toms748
   */
  template <class fp_t=double> class gsl_tolerance {
    
  public:
    
    gsl_tolerance(fp_t tol_abs, fp_t tol_rel) {
      epsabs=tol_abs;
      epsrel=tol_rel;
    }
    
    bool operator()(const fp_t &a, const fp_t &b) {
      fp_t abs_a=fabs(a);
      fp_t abs_b=fabs(b);
      fp_t min_abs;
      if ((a > 0.0 && b > 0.0) || (a < 0.0 && b < 0.0)) {
	if (abs_a<abs_b) min_abs=abs_a;
	else min_abs=abs_b;
      } else {
	min_abs=0.0;
      }
      
      fp_t tolerance = epsabs+epsrel*min_abs;
      
      if (fabs(b-a) < tolerance) {
	return true;
      }
      return false;
    }
    
  private:
    
    fp_t epsabs;
    fp_t epsrel;
    
  };
  
  /** \brief Bracketing solver based the Boost implementation of TOMS
      748
      
      This class currently uses \ref o2scl::gsl_tolerance as a test,
      since this works even when the root is zero.
   */
  template<class func_t=funct, class fp_t=double> class root_toms748 : 
    public root_bkt<func_t,func_t,fp_t> {
    
  public:
    
  root_toms748() {
  }
  
  virtual ~root_toms748() {}
  
  /// Return the type, \c "root_toms748".
  virtual const char *type() { return "root_toms748"; }
  
  /// Solve \c func using \c x as an initial guess, returning \c x.
  virtual int solve_bkt(fp_t &x1, fp_t x2, func_t &func) {
    std::pair<fp_t,fp_t> res;
    size_t digits;
    if (this->tol_rel>1.0) digits=1;
    else if (this->tol_rel<=0.0) digits=18;
    else digits=((size_t)(-log10(this->tol_rel)));
    gsl_tolerance<fp_t> tol(this->tol_abs,this->tol_rel);
    //boost::math::tools::eps_tolerance<fp_t> tol(digits);
    size_t niter=((size_t)this->ntrial);
    res=boost::math::tools::toms748_solve(func,x1,x2,tol,niter);
    this->last_ntrial=niter;
    x1=res.first;
    return 0;
  }
      
  };
  
}
  
#endif
