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
/** \file eos_had_ddc.h
    \brief File defining \ref o2scl::eos_had_ddc
*/
#ifndef O2SCL_DDC_EOS_H
#define O2SCL_DDC_EOS_H

#include <string>
#include <cmath>
#include <o2scl/lib_settings.h>
#include <o2scl/constants.h>
#include <o2scl/part.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/fermion.h>
#include <o2scl/mm_funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Relativistic mean field EOS with density dependent couplings
      
      Based on \ref Typel99.

      \future Implement the finite temperature EOS properly.
  */
  class eos_had_ddc : public eos_had_eden_base {
  public:

    /// \name Masses
    //@{
    /// nucleon mass
    double mnuc;
    /// \f$ \phi \f$ mass (in \f$ \mathrm{fm}^{-1} \f$ )
    double ms;
    /// \f$ A_{\omega} \f$ mass (in \f$ \mathrm{fm}^{-1} \f$ )
    double mw;
    /// \f$ A_{\rho} \f$ mass (in \f$ \mathrm{fm}^{-1} \f$ )
    double mr;
    //@}

    /// \name Parameters for couplings
    //@{
    /// The coupling \f$ \Gamma_{\sigma}(\rho_{\mathrm{sat}}) \f$
    double Gs;
    /// The coupling \f$ \Gamma_{\omega}(\rho_{\mathrm{sat}}) \f$
    double Gw;
    /// The coupling \f$ \Gamma_{\rho}(\rho_{\mathrm{sat}}) \f$
    double Gr;
    /// \f$ a_{\sigma} \f$
    double as;
    /// \f$ a_{\omega} \f$
    double aw;
    /// \f$ a_{\rho} \f$
    double ar;
    /// \f$ b_{\sigma} \f$
    double bs;
    /// \f$ b_{\omega} \f$
    double bw;
    /// \f$ c_{\sigma} \f$
    double cs;
    /// \f$ c_{\omega} \f$
    double cw;
    /// \f$ d_{\sigma} \f$
    double ds;
    /// \f$ d_{\omega} \f$
    double dw;
    //@}

    // The saturation density
    double rho0;

    eos_had_ddc();

    /// Equation of state as a function of the densities
    virtual int calc_e(fermion &n, fermion &p, thermo &th) {
      return exc_eunimpl;
    }
    
    /** \brief Equation of state and meson field equations
	as a function of the density
      
	This calculates the pressure and energy density as a function
	of \f$ \mu_n, \mu_p, \phi, A_{\omega}, A_{\rho} \f$ . When the
	field equations have been solved, \c f1, \c f2, and \c f3 are
	all zero.
	
	\todo Is the thermodynamic identity is satisfied even when the
	field equations are not solved? Check this.
    */
    virtual int calc_eq_e(fermion &neu, fermion &p, double sig, 
			  double ome, double rho, double &f1, 
			  double &f2, double &f3, thermo &th);
    

    /// Return string denoting type ("eos_had_ddc")
    virtual const char *type() { return "eos_had_ddc"; }

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Zero-temperature fermion thermodynamics
    fermion_zerot fzt;
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
