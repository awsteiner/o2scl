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

  /** \brief
   */
  class eos_leptons {
    
  public:
    
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
    
    /** \brief If true, include muons (default true)
     */
    bool include_muons;

    /** \brief If true, include photons (default false)
     */
    bool include_photons;

    /** \brief Thermodynamic quantities
     */
    thermo th;

    eos_leptons() {
      include_muons=true;
      include_photons=false;
      
      convert_units<double> &cu=o2scl_settings.get_convert_units();
      e.init(cu.convert("kg","1/fm",o2scl_mks::mass_electron),2.0);
      mu.init(cu.convert("kg","1/fm",o2scl_mks::mass_muon),2.0);
                         
      ph.init(0.0,2.0);
    }

    /** \brief Thermodynamics from the electron and muon 
        chemical potentials
     */
    int pair_mu(double T) {
      frel.pair_mu(e,T);
      th.ed=e.ed;
      th.pr=e.pr;
      th.en=e.en;
      if (include_muons) {
        frel.pair_mu(mu,T);
        th.ed+=mu.ed;
        th.pr+=mu.pr;
        th.en+=mu.en;
      }
      if (include_photons) {
        ph.massless_calc(T);
        th.ed+=ph.ed;
        th.pr+=ph.pr;
        th.en+=ph.en;
      }
      return 0;
    }

    /** \brief Thermodynamics from the electron chemical potential
        in weak equilibrium
    */
    int pair_mu_eq(double T) {
      if (include_muons) {
        mu.mu=e.mu;
      }
      pair_mu(T);
      return 0;
    }

    /** \brief Thermodynamics from the electron and muon densities
     */
    int pair_density(double T) {
      
      bool fr_en=frel.err_nonconv;
      frel.err_nonconv=false;

      int retx;
      if (e.inc_rest_mass) {
        e.inc_rest_mass=false;
        e.mu-=e.m;
        retx=frel.pair_density(e,T);
        e.inc_rest_mass=true;
        e.mu+=e.m;
      } else {
        retx=frel.pair_density(e,T);
      }
      
      // Sometimes the solver fails, but we can recover by adjusting
      // the upper limit for degenerate fermions and tightening the
      // integration tolerances
      if (retx!=0) {
        
        frel.upper_limit_fac=40.0;
        frel.def_dit.tol_rel/=1.0e2;
        frel.def_dit.tol_abs/=1.0e2;
        frel.def_nit.tol_rel/=1.0e2;
        frel.def_nit.tol_abs/=1.0e2;
        
        if (e.inc_rest_mass) {
          e.inc_rest_mass=false;
          e.mu-=e.m;
          retx=frel.pair_density(e,T);
          e.inc_rest_mass=true;
          e.mu+=e.m;
        } else {
          retx=frel.pair_density(e,T);
        }

        if (retx!=0) {
          O2SCL_ERR2("Function pair_density() for electrons failed in ",
                     "class eos_leptons().",o2scl::exc_efailed);
        }
        
        frel.upper_limit_fac=20.0;
        frel.def_dit.tol_rel*=1.0e2;
        frel.def_dit.tol_abs*=1.0e2;
        frel.def_nit.tol_rel*=1.0e2;
        frel.def_nit.tol_abs*=1.0e2;
        
      }
      
      th.ed=e.ed;
      th.pr=e.pr;
      th.en=e.en;

      if (include_muons) {

        if (mu.inc_rest_mass) {
          mu.inc_rest_mass=false;
          mu.mu-=mu.m;
          retx=frel.pair_density(mu,T);
          mu.inc_rest_mass=true;
          mu.mu+=mu.m;
        } else {
          retx=frel.pair_density(mu,T);
        }
        
        // Sometimes the solver fails, but we can recover by adjusting
        // the upper limit for degenerate fermions and tightening the
        // integration tolerances
        if (retx!=0) {
          
          frel.upper_limit_fac=40.0;
          frel.def_dit.tol_rel/=1.0e2;
          frel.def_dit.tol_abs/=1.0e2;
          frel.def_nit.tol_rel/=1.0e2;
          frel.def_nit.tol_abs/=1.0e2;
          
          if (mu.inc_rest_mass) {
            mu.inc_rest_mass=false;
            mu.mu-=mu.m;
            retx=frel.pair_density(mu,T);
            mu.inc_rest_mass=true;
            mu.mu+=mu.m;
          } else {
            retx=frel.pair_density(mu,T);
          }
          
          if (retx!=0) {
            O2SCL_ERR2("Function pair_density() for muons failed in ",
                       "class eos_leptons().",o2scl::exc_efailed);
          }
        
          frel.upper_limit_fac=20.0;
          frel.def_dit.tol_rel*=1.0e2;
          frel.def_dit.tol_abs*=1.0e2;
          frel.def_nit.tol_rel*=1.0e2;
          frel.def_nit.tol_abs*=1.0e2;
          
        }
      
        th.ed+=mu.ed;
        th.pr+=mu.pr;
        th.en+=mu.en;
      }
      
      if (include_photons) {
        ph.massless_calc(T);
        th.ed+=ph.ed;
        th.pr+=ph.pr;
        th.en+=ph.en;
      }

      frel.err_nonconv=fr_en;
      
      return 0;
    }
    
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
