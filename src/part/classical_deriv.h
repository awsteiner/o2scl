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
#ifndef O2SCL_SN_CLASSICAL_H
#define O2SCL_SN_CLASSICAL_H

/** \file classical_deriv.h
    \brief File defining \ref o2scl::classical_deriv
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include <o2scl/constants.h>
#include <o2scl/part_deriv.h>
#include <o2scl/classical.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a classical particle with derivatives
   */
  class classical_deriv {

  protected:
    
    /// For computing non-derivative properties
    classical cl;

  public:

    classical_deriv();

    virtual ~classical_deriv();
    
    /** \brief Compute the properties of particle \c p at temperature 
	\c temper from its chemical potential
    */
    virtual void calc_mu(part_deriv &p, double temper);
    
    /** \brief Compute the properties of particle \c p at temperature 
	\c temper from its density
    */
    virtual void calc_density(part_deriv &p, double temper);

    /// Return string denoting type ("classical_deriv")
    virtual const char *type() { return "classical_deriv"; };

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
