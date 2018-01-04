/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_CLASSICAL_H
#define O2SCL_CLASSICAL_H

/** \file classical.h
    \brief File defining \ref o2scl::classical
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/part.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Classical particle class
      
      \future Write a calc_density_zerot() function for 
      completeness?
  */
  class classical {

  public:

    /** \brief Create a classical particle with mass \c m  
	and degeneracy \c g 

	\note This class attempts to handle zero temperature limit
	somewhat gracefully, even though the classical limit doesn't
	necessarily make physical sense there.
    */
    classical();

    virtual ~classical() {
    }

    /** \brief Calculate properties as function of chemical potential
	
	If the temperature is less than zero, the error
	handler will be called. 

	\future Handle the case \f$ \mu/T>308 \f$ properly.
     */
    virtual void calc_mu(part &p, double temper);

    /** \brief Calculate properties as function of density

	If the density or the temperature is less than zero, the error
	handler will be called. In the case of zero density, the
	chemical potential is set to the mass and the energy density,
	pressure, and entropy are set to zero. 
     */
    virtual void calc_density(part &p, double temper);

    /// Return string denoting type ("classical")
    virtual const char *type() { return "classical"; }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
