/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner and Jerry Gagelman
  
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
#ifndef O2SCL_INTE_H
#define O2SCL_INTE_H

/** \file inte.h
    \brief File defining \ref o2scl::inte
*/

#include <cmath>
#include <o2scl/err_hnd.h>
#include <o2scl/funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Base integration class [abstract base]

      \future It might be useful to have all of the integration
      classes report the number of function evaluations used
      in addition to the number of iterations which were taken
  */
  template<class func_t=funct> class inte {
    
  public:
  
  inte() {
    tol_rel=1.0e-8;
    tol_abs=1.0e-8;
    verbose=0;
    interror=0.0;
    err_nonconv=true;
  }

  virtual ~inte() {}

  /** \brief Verbosity
   */
  int verbose;

  /// The most recent number of iterations taken
  size_t last_iter;
  
  /** \brief The maximum relative uncertainty 
      in the value of the integral (default \f$ 10^{-8} \f$)
  */
  double tol_rel;
    
  /** \brief The maximum absolute uncertainty 
      in the value of the integral (default \f$ 10^{-8} \f$)
  */
  double tol_abs;
  
  /** \brief If true, call the error handler if the routine does not 
      converge or reach the desired tolerance (default true)

      If this is false, the function proceeds normally and 
      may provide convergence information in the integer return value.
  */
  bool err_nonconv;
  
  /** \brief Integrate function \c func from \c a to \c b.
   */
  virtual double integ(func_t &func, double a, double b) {
    double res;
    int ret=integ_err(func,a,b,res,this->interror);
    if (ret!=0) {
      O2SCL_ERR2("Integration failed in inte::integ(), ",
		 "but cannot return int.",o2scl::exc_efailed);
    }
    return res;
  }

  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, double a, double b, 
			double &res, double &err)=0;
  
  /** \brief Return the numerically estimated error in the result from
      the last call to integ()
      
      This will quietly return zero if no integrations have been
      performed or if the integrator does not estimate the error.
  */
  double get_error() { return interror; }
  
  /// Return string denoting type ("inte")
  virtual const char *type() { return "inte"; }
  
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




