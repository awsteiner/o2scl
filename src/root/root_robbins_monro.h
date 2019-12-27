/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2020, Andrew W. Steiner
  
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
#ifndef O2SCL_ROOT_ROBBINS_MONRO_H
#define O2SCL_ROOT_ROBBINS_MONRO_H

/** \file root_robbins_monro.h
    \brief File defining \ref o2scl::root_robbins_monro 
*/

#include <limits>
#include <o2scl/funct.h>
#include <o2scl/root.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief One-dimensional root-finding for noisy functions
   */
  template<class func_t=funct, class dfunc_t=func_t>
    class root_robbins_monro : public root<func_t,dfunc_t> {
    
  public:
  
  root_robbins_monro() {
    coeff=1.0;
    gamma=1.0;
  }

  /// If true, apply Polyak-Ruppert averaging
  bool pr_averaging;

  /// Coefficient for sequence (default 1.0)
  double coeff;

  /// Exponent for sequence (default 1.0)
  double gamma;
  
  /// Return the type, \c "root_robbins_monro".
  virtual const char *type() { return "root_robbins_monro"; }

  /** \brief Solve \c func using \c x as an initial guess
   */
  virtual int solve(double &x, func_t &func) {
    double x2;

    if (pr_averaging) {
      
      double avg=0.0;
      for(size_t i=0;i<this->ntrial/2;i++) {
	double fx=func(x);
	if (fabs(fx)<this->tol_rel) {
	  return o2scl::success;
	}
	x2=x-fx*coeff/pow(((double)(i+1)),gamma);
	avg+=(x2*((double)i)+avg)/((double)(i+1));
	double favg=func(avg);
	if (fabs(favg)<this->tol_rel) {
	  x=avg;
	  return o2scl::success;
	}
	x=x2;
      }
      O2SCL_CONV2_RET("Function root_robbins_monro::solve() exceeded ",
		      "maximum number of iterations.",o2scl::exc_emaxiter,
		      this->err_nonconv);
    } 

    for(size_t i=0;i<this->ntrial;i++) {
      double fx=func(x);
      if (fabs(fx)<this->tol_rel) {
	return o2scl::success;
      }
      x2=x-f(x)*coeff/pow(((double)(i+1)),gamma);
      x=x2;
    }
    fx=func(x);
    if (fabs(fx)<this->tol_rel) {
      return o2scl::success;
    }
    O2SCL_CONV2_RET("Function root_robbins_monro::solve() exceeded ",
		    "maximum number of iterations.",o2scl::exc_emaxiter,
		    this->err_nonconv);
    return o2scl::exc_emaxiter;
  }
  
  /// Solve \c func in region \f$ x_1<x<x_2 \f$ returning \f$ x_1 \f$.
  virtual int solve_bkt(double &x1, double x2, func_t &f) {
    x1=(x1+x2)/2.0;
    return solve(x1,f);
  }

  /** \brief Solve \c func using \c x as an initial
      guess using derivatives \c df.
  */
  virtual int solve_de(double &x, func_t &func, dfunc_t &df) {
    double x2;

    if (pr_averaging) {
      
      double avg=0.0;
      for(size_t i=0;i<this->ntrial/2;i++) {
	double fx=func(x);
	double fpx=df(x);
	if (fabs(fx)<this->tol_rel) {
	  return o2scl::success;
	}
	x2=x-fx/df*coeff/pow(((double)(i+1)),gamma);
	avg+=(x2*((double)i)+avg)/((double)(i+1));
	double favg=func(avg);
	if (fabs(favg)<this->tol_rel) {
	  x=avg;
	  return o2scl::success;
	}
	x=x2;
      }
      O2SCL_CONV2_RET("Function root_robbins_monro::solve() exceeded ",
		      "maximum number of iterations.",o2scl::exc_emaxiter,
		      this->err_nonconv);
    } 

    for(size_t i=0;i<this->ntrial;i++) {
      double fx=func(x);
      double fpx=df(x);
      if (fabs(fx)<this->tol_rel) {
	return o2scl::success;
      }
      x2=x-fx/df*coeff/pow(((double)(i+1)),gamma);
      x=x2;
    }
    fx=func(x);
    if (fabs(fx)<this->tol_rel) {
      return o2scl::success;
    }
    O2SCL_CONV2_RET("Function root_robbins_monro::solve() exceeded ",
		    "maximum number of iterations.",o2scl::exc_emaxiter,
		    this->err_nonconv);
    return o2scl::exc_emaxiter;
  }
  
  };
   
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
