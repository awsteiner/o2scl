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
/** \file eos_base.h
    \brief File defining \ref o2scl::eos_base
*/
#ifndef O2SCL_EOS_BASE_H
#define O2SCL_EOS_BASE_H

#include <o2scl/part.h>
#include <o2scl/fermion.h>
#include <o2scl/boson.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/part_deriv.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/eos_leptons.h>

namespace o2scl {

  /** \brief Equation of state base class [abstract virtual base]
    
      A base class for an equation of state
  */
  class eos_base {

  public:

    eos_base();

    virtual ~eos_base() {};

    /// Set class thermo object
    virtual void set_thermo(thermo &th);

    /// The default thermo object
    thermo def_thermo;

    /// Return string denoting type ("eos_base")
    virtual const char *type() { return "eos_base"; }
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    
    /** \brief Compute the EOS in beta-equilibrium at 
	zero temperature

        This solves the function \ref solve_beta_eq_T0(). This
        function is different from \ref nstar_cold because it is a
        generic interface which works for non-hadronic EOSs.
    */
    virtual int beta_eq_T0(ubvector &nB_grid, ubvector &guess,
                           eos_leptons &elep,
			   std::shared_ptr<table_units<> > results)=0;

    /** \brief Equation of state as a function of baryon, charge,
        and strangeness density at finite temperature
    */
    virtual int calc_temp_f_gen(double nB, double nQ, double nS,
                                double T, thermo &th)=0;
    
  protected:

    /// A pointer to the thermo object
    thermo *eos_thermo;

  };

  
}

#endif
