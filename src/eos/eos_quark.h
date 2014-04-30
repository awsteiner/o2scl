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
#ifndef O2SCL_QUARK_EOS_H
#define O2SCL_QUARK_EOS_H

#include <iostream>
#include <o2scl/eos_base.h>
#include <o2scl/quark.h>
#include <o2scl/fermion_eff.h>
#include <o2scl/deriv.h>
#include <o2scl/mroot.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Quark matter equation of state base
  */
  class eos_quark : public eos_base {

  public:

    eos_quark();

    virtual ~eos_quark() {};

    /// Calculate equation of state as a function of chemical potentials
    virtual int calc_p(quark &u, quark &d, quark &s, thermo &th);

    /// Calculate equation of state as a function of density
    virtual int calc_e(quark &u, quark &d, quark &s, thermo &th);

    /// Calculate equation of state as a function of chemical potentials
    virtual int calc_temp_p(quark &u, quark &d, quark &s, 
			    double temper, thermo &th);
  
    /// Calculate equation of state as a function of density
    virtual int calc_temp_e(quark &u, quark &d, quark &s, 
			    double temper, thermo &th);
  
    /// Return string denoting type ("eos_quark")
    virtual const char *type() { return "eos_quark"; }

    /// Object for computing fermion thermodynamics
    fermion_eval_thermo *fet;

    /// Default fermion thermodynamics
    fermion_eff def_fet;
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



