/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#ifndef O2SCL_BOSON_REL_H
#define O2SCL_BOSON_REL_H

/** \file boson_rel.h
    \brief File defining \ref o2scl::boson_rel
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/root.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/root_cern.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>

#include <o2scl/boson.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a relativistic boson
      
      \todo Testing not completely finished.
  */
  class boson_rel {

  public:

    /// Create a boson with mass \c m and degeneracy \c g
    boson_rel();

    virtual ~boson_rel();
    
    /** \brief Calculate properties as function of chemical potential
     */
    virtual void calc_mu(boson &b, double temper);
    
    /** \brief Calculate properties as function of density
     */
    virtual void calc_density(boson &b, double temper);
    
    /** \brief Calculate properties with antiparticles as function of
        chemical potential
    */
    virtual void pair_mu(boson &b, double temper);

    /** \brief Calculate properties with antiparticles as function of
	density
    */
    virtual void pair_density(boson &b, double temper) {
      O2SCL_ERR("Function boson_rel::pair_density() unimplemented.",
		exc_eunimpl);
    }

    /// Calculate effective chemical potential from density
    virtual void nu_from_n(boson &b, double temper);

    /// Set degenerate and nondegenerate integrators
    void set_inte(inte<> &l_nit, inte<> &l_dit);

    /** \brief Set the solver for use in calculating the chemical
	potential from the density */
    void set_density_root(root<> &rp) {
      density_root=&rp;
      return;
    }

    /// The default solver for calc_density().
    root_cern<> def_density_root;

    /// Default nondegenerate integrator
    inte_qagiu_gsl<> def_nit;

    /// Default degenerate integrator
    inte_qag_gsl<> def_dit;

    /// Return string denoting type ("boson_rel")
    virtual const char *type() { return "boson_rel"; }

  protected:

#ifndef DOXYGEN_NO_O2NS

    /// The non-degenerate integrator
    inte<> *nit;
    /// The degenerate integrator
    inte<> *dit;
    /// The solver for calc_density()
    root<> *density_root;

    /// Non-degenerate density integral
    double density_fun(double u, boson &b, double T);
    /// Non-degenerate energy density integral
    double energy_fun(double u, boson &b, double T);
    /// Non-degenerate entropy integral
    double entropy_fun(double u, boson &b, double T);
    /// Degenerate density integral
    double deg_density_fun(double u, boson &b, double T);
    /// Degenerate energy density integral
    double deg_energy_fun(double u, boson &b, double T);
    /// Degenerate entropy integral
    double deg_entropy_fun(double u, boson &b, double T);
    /// Solve for the density in calc_density()
    double solve_fun(double x, boson &b, double T);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
