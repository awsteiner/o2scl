/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2017, Andrew W. Steiner
  
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
#ifndef O2SCL_MIN_BRENT_BOOST_H
#define O2SCL_MIN_BRENT_BOOST_H

/** \file min_brent_boost.h
    \brief File defining \ref o2scl::min_brent_boost
*/

#include <boost/math/tools/minima.hpp>

#include <o2scl/min.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
   
  /** \brief One-dimensional minimization using Brent's method (GSL)

      This is a wrapper for
      <tt>boost::math::tools::brent_find_minima()</tt> .

      The value of fourth argument, <tt>digits</tt> is 
      taken to be \f$ - \log_{10}(t) \f$ where \f$ t \f$ is
      the value of \ref o2scl::min_base::tol_rel .
  */
  template<class func_t=funct11> class min_brent_boost : 
  public min_bkt_base<func_t> {
    
  public:
 
  min_brent_boost() {
  }
      
  virtual ~min_brent_boost() {}
  
  /** \brief Calculate the minimum \c fmin of \c func 
      with \c x2 bracketed between \c x1 and \c x3.

      The initial value of \c x2 is ignored.
  */
  virtual int min_bkt(double &x2, double x1, double x3, double &fmin,
		      func_t &func) {

    std::pair<double,double> res;
    size_t digits;
    if (this->tol_rel>1.0) digits=1;
    else if (this->tol_rel<=0.0) digits=18;
    else digits=((size_t)(-log10(this->tol_rel)));
    res=boost::math::tools::brent_find_minima(func,x1,x3,digits);
    x2=res.first;
    fmin=res.second;
    return success;
  }
  
  /// Return string denoting type ("min_brent_boost")
  virtual const char *type() { return "min_brent_boost"; }
  
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
