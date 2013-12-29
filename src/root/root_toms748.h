/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Andrew W. Steiner
  
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

  /** \brief Wrapper for Boost 1-d solver
      
      This class does not work with \ref funct objects, only with
      functions of type <tt>std::function<double(double)></tt>.
   */
#ifdef O2SCL_CPP11
  template<class func_t=funct11> class root_toms748 : public root_bkt<func_t>
#else
    template<class func_t> class root_toms748 : public root_bkt<func_t>
#endif
  {
    
  public:
    
  root_toms748() {
  }
  
  virtual ~root_toms748() {}
  
  /// Return the type, \c "root_toms748".
  virtual const char *type() { return "root_toms748"; }
  
  /// Solve \c func using \c x as an initial guess, returning \c x.
  virtual int solve_bkt(double &x1, double x2, func_t &func) {
    std::pair<double,double> res;
    size_t digits;
    if (this->tol_rel>1.0) digits=1;
    else if (this->tol_rel<=0.0) digits=18;
    else digits=((size_t)(-log10(this->tol_rel)));
    boost::math::tools::eps_tolerance<double> tol(digits);
    size_t niter=((size_t)this->ntrial);
    res=boost::math::tools::toms748_solve(func,x1,x2,tol,niter);
    this->last_ntrial=niter;
    x1=res.first;
    return 0;
  }
      
  };
  
}
  
#endif
