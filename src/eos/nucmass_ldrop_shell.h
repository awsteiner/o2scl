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
/** \file nucmass_ldrop_shell.h
    \brief File defining \ref o2scl::nucmass_ldrop_shell
*/
#ifndef LDROP_SHELL_H
#define LDROP_SHELL_H

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>
#include <o2scl/constants.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/mmin_conp.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/nucmass_ldrop.h>
#include <o2scl/nucmass_frdm.h>

namespace o2scl {

  /** \brief Liquid drop model with shell effects
  */
  class nucmass_ldrop_shell : public nucmass_ldrop_pair, 
    public nucmass_ibm_shell {

  public:

    nucmass_ldrop_shell();

    virtual ~nucmass_ldrop_shell() {}

    /// If true, include shell effects (default true)
    bool inc_shell;

    /** \brief Return the free binding energy of a \nucleus in a many-body 
	environment
    */
    virtual double drip_binding_energy_d
      (double Z, double N, double npout, double nnout,
       double chi, double T);

    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);

  };

  /** \brief Mass formula adding simple shell effects to the 
      finite-range liquid droplet model
   */
  class nucmass_frdm_shell : public nucmass_frdm, public nucmass_ibm_shell {
    
  public:

    nucmass_frdm_shell();

    virtual ~nucmass_frdm_shell() {}

    /// Compute the mass excess
    double mass_excess_d(double Z, double N);

    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);

  };

}

#endif
