/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/mroot_hybrids.h>
#include <o2scl/inte.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>

#include <o2scl/boson.h>

namespace o2scl {

  /** \brief Equation of state for a relativistic boson
      
      \verbatim embed:rst

      .. todo:: 

         - In class boson_rel: Testing not completely finished.
         
      \endverbatim
  */
  class boson_rel : public boson_thermo {

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
    
    /** \brief Calculate the maximum density as a function of temperature
     */
    virtual void calc_max_density(boson &b, double temper);
    
    /** \brief Calculate properties with antiparticles as function of
        chemical potential
    */
    virtual void pair_mu(boson &b, double temper);

    /** \brief Calculate properties with antiparticles as function of
	density
    */
    virtual void pair_density(boson &b, double temper);

    /// Calculate effective chemical potential from density
    virtual void nu_from_n(boson &b, double temper);

    /// Set degenerate and nondegenerate integrators
    void set_inte(inte<> &l_nit, inte<> &l_dit);

    /** \brief Set the solver for use in calculating the chemical
	potential from the density */
    void set_density_mroot(mroot<> &rp) {
      density_mroot=&rp;
      return;
    }

    /// The default solver for calc_density().
    mroot_hybrids<> def_density_mroot;

    /// Default nondegenerate integrator
    inte_qagiu_gsl<> def_nit;

    /// Default degenerate integrator
    inte_qag_gsl<> def_dit;

    /// Return string denoting type ("boson_rel")
    virtual const char *type() { return "boson_rel"; }

    /** \brief If true, verify the thermodynamic identity
     */
    bool verify_ti;

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief Verbosity parameter
     */
    bool use_expansions;

    double deg_limit;

    double upper_limit_fac;
    
  protected:

    /// The non-degenerate integrator
    inte<> *nit;
    /// The degenerate integrator
    inte<> *dit;
    /// The solver for calc_density()
    mroot<> *density_mroot;

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
    int solve_fun(size_t nv, const ubvector &x, ubvector &y,
                  double density, boson &b, double T);
    /// Desc
    int pair_density_fun(size_t nv, const ubvector &x, ubvector &y,
                         double density, boson &b, double T);

  };

}

#endif
