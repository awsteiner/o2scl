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
/** \file eos_base.h
    \brief File defining \ref o2scl::eos_base
*/
#ifndef O2SCL_EOS_H
#define O2SCL_EOS_H

#include <o2scl/part.h>
#include <o2scl/fermion.h>
#include <o2scl/boson.h>
#include <o2scl/fermion_rel.h>

namespace o2scl {

  /** \brief Lepton and photon EOS
   */
  class eos_leptons {

  public:

    // Vector type
    typedef boost::numeric::ublas::vector<double> ubvector;
    
  protected:

    /** \brief Electron thermodynamics from the electron density

        \note This internal function presumes that 
        \ref fermion_rel::err_nonconv is false .
     */
    int electron_density(double T);

    /** \brief Function to solve for \ref pair_density_eq()

        \note This internal function presumes that include_muons is
        true (otherwise there's nothing to solve) and that \ref
        fermion_rel::err_nonconv is false .
     */
    int pair_density_eq_fun(size_t nv, const ubvector &x,
                            ubvector &y, double T, double nq);

  public:
    
    /// Solver for \ref pair_density_eq() 
    mroot_hybrids<> mh;
    
    /** \brief Relativistic fermion thermodyanmics
     */
    fermion_rel frel;
    
    /** \brief Electron
     */
    fermion e;
    
    /** \brief Muon
     */
    fermion mu;

    /** \brief Photon
     */
    boson ph;
    
    /** \brief If true, call the error handler if msolve()
        does not converge (default true)
    */
    bool err_nonconv;
    
    /** \brief If true, include muons (default true)
     */
    bool include_muons;

    /** \brief If true, include photons (default false)
     */
    bool include_photons;

    /** \brief Thermodynamic quantities
     */
    thermo th;

    eos_leptons();

    /** \brief Thermodynamics from the electron and muon 
        chemical potentials
     */
    int pair_mu(double T);

    /** \brief Thermodynamics from the electron chemical potential
        in weak equilibrium
    */
    int pair_mu_eq(double T);

    /** \brief Thermodynamics from the electron and muon densities
     */
    int pair_density(double T);
    
    /** \brief Thermodynamics from the charge density in 
        weak equilibrium

        The first argument \c nq, is the total negative charge density
        including electrons (and muons if \ref include_muons is true)
        and \c T is the temperature.
    */
    int pair_density_eq(double nq, double T);
    
  };
  
  /** \brief Equation of state base class
    
      A base class for an equation of state
  */
  class eos_base {

  public:

    eos_base();

    virtual ~eos_base() {};

    /// Set class thermo object
    virtual void set_thermo(thermo &th);

    /// Get class thermo object
    virtual const thermo &get_thermo();

    /// The default thermo object
    thermo def_thermo;

    /// Return string denoting type ("eos_base")
    virtual const char *type() { return "eos_base"; }
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    
    /** \brief Compute the EOS in beta-equilibrium at 
	zero temperature
    */
    virtual int beta_eq_T0(ubvector &nB_grid, ubvector &guess,
			   fermion &e, bool include_muons,
			   fermion &mu, fermion_rel &frel,
			   std::shared_ptr<table_units<> > results);
    
  protected:

    /// A pointer to the thermo object
    thermo *eos_thermo;

  };

  
}

#endif
