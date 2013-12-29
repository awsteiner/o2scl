/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Jerry Gagelman
  and Andrew W. Steiner
  
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
#ifndef O2SCL_GSL_INTE_QAGS_H
#define O2SCL_GSL_INTE_QAGS_H

#include <o2scl/inte.h>
#include <o2scl/funct.h>
#include <o2scl/inte_singular_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Integrate a function with a singularity (GSL)

      If a function is unbounded but has a finite integral, using the
      adaptive algorithm described for \ref inte_qag_gsl to compute
      that integral (up to specified tolerance) will converge very
      slowly. The integration routine of this class remedies this by
      combining the adaptive algorithm with a series-acceleration
      method.
      
      See \ref gslinte_subsect in the User's guide for general
      information about the GSL integration classes.

      \comment
      Note that it's important that this is separate from
      inte_singular_gsl::qags(), since this class uses set_rule(2)
      while other children of inte_singular_gsl do not.
      \endcomment
  */
  template<class func_t=funct> class inte_qags_gsl : 
  public inte_singular_gsl<func_t> {
    
  public:

  inte_qags_gsl() {
    this->set_rule(2);
  }
      
  virtual ~inte_qags_gsl() {
  }

  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, double a, double b, 
			double &res, double &err) {
    return this->qags(func,a,b,this->tol_abs,this->tol_rel,&res,&err);
  }
           
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
