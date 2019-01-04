/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
/** \file eos_quark_bag.h
    \brief File defining \ref o2scl::eos_quark_bag
*/
#ifndef O2SCL_BAG_EOS_H
#define O2SCL_BAG_EOS_H

#include <o2scl/constants.h>
#include <o2scl/part.h>
#include <o2scl/eos_quark.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Simple bag model

      An equation of state with \f$ P=-B+P_{u,FG}+P_{d,FG}+P_{s,FG}
      \f$ where \f$ P_{i,FG} \f$ is the Fermi gas contribution from
      particle \f$ i \f$ and \f$ B \f$ is a density- and
      temperature-independent bag constant.

      The finite temperature functions run the zero temperature code 
      if the temperature is less than or equal to 0.
  */

  class eos_quark_bag : public eos_quark {

  public:

    eos_quark_bag();

    virtual ~eos_quark_bag() {};

    virtual int calc_p(quark &u, quark &d, quark &s, thermo &th);

    virtual int calc_e(quark &u, quark &d, quark &s, thermo &th);

    /** \brief Calculate equation of state as a function of 
	the chemical potentials
	
	This function returns zero (success) unless the 
	call to quark::pair_mu() fails.
    */
    virtual int calc_temp_p(quark &u, quark &d, quark &s, 
			    double temper, thermo &th);

    /** \brief Calculate equation of state as a function of 
	the densities
	
	This function returns zero (success) unless the 
	call to quark::pair_density() fails.
    */
    virtual int calc_temp_e(quark &u, quark &d, quark &s, 
			    double temper, thermo &th);

    /// The bag constant in \f$ fm^{-4}\f$ (default \f$200/(\hbar c) \f$).
    double bag_constant;

    /// Return string denoting type ("eos_quark_bag")
    virtual const char *type() { return "eos_quark_bag"; }
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
