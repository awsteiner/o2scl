/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
/** \file eos_quark.h
    \brief File defining \ref o2scl::eos_quark
*/
#ifndef O2SCL_QUARK_EOS_H
#define O2SCL_QUARK_EOS_H

#include <iostream>
#include <o2scl/eos_base.h>
#include <o2scl/quark.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/deriv.h>
#include <o2scl/mroot.h>

namespace o2scl {

  /** \brief Quark matter equation of state base
  */
  class eos_quark : public eos_base {

  public:

    /** \name The default quark objects
	
	The masses are automatically set in the constructor to
	\c up_default_mass, \c down_default_mass, and 
	\c strange_default_mass.c
    */
    //@{
    quark def_up; 
    quark def_down;
    quark def_strange;
    //@}
    
    eos_quark();

    virtual ~eos_quark() {};

    /// Calculate equation of state as a function of chemical potentials
    virtual int calc_p(quark &u, quark &d, quark &s, thermo &th);

    /// Calculate equation of state as a function of density
    virtual int calc_e(quark &u, quark &d, quark &s, thermo &th);

    /** \brief Calculate equation of state as a function of chemical 
	potentials at finite temperature
    */
    virtual int calc_temp_p(quark &u, quark &d, quark &s, 
			    double temper, thermo &th);
  
    /** \brief Calculate equation of state as a function of density
	at finite temperature
    */
    virtual int calc_temp_e(quark &u, quark &d, quark &s, 
			    double temper, thermo &th);
  
    /** \brief Equation of state as a function of baryon, charge,
        and strangeness density at finite temperature
    */
    virtual int calc_temp_f_gen(double nB, double nQ, double nS, double T,
                                thermo &th) {
      
      def_up.n=0.75*(2.0*nB+nQ);
      def_down.n=0.25*(6.0*nB-3.0*nQ-4.0*nS);
      def_strange.n=nS;
      return calc_temp_e(def_up,def_down,def_strange,T,th);
    }
    
    /// Return string denoting type ("eos_quark")
    virtual const char *type() { return "eos_quark"; }

    /// Object for computing fermion thermodynamics
    o2scl::fermion_rel fet;

  };

}

#endif


