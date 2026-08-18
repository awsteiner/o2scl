/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2026, Andrew W. Steiner
  
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
#ifndef O2SCL_EOS_LEPTONS_MULTIP_H
#define O2SCL_EOS_LEPTONS_MULTIP_H

#include <o2scl/set_multip.h>
#include <o2scl/fermion_multip.h>
#include <o2scl/eos_leptons.h>

namespace o2scl {

#if defined (O2SCL_SET_MULTIP) || defined (DOXYGEN)

  /** \brief Multiprecision version of \ref eos_leptons

      This class provides the additional functionality on top of \ref
      eos_leptons to compute the lepton and photon EOS with
      multiprecision.
   */
  class eos_leptons_multip : public eos_leptons {

  public:

    typedef boost::numeric::ublas::vector<long double> ubvector_ld;
    
  protected:
    
    /** \brief Photon (long double)
     */
    boson_ld ph_ld;
    
    /** \brief Photon (25-digit version)
     */
    boson_fp25 ph_fp25;

    /** \brief Photon derivatives (long double version)
     */
    part_deriv_press_ld phd_ld;

    /** \brief Photon derivatives (25-digit version)
     */
    part_deriv_press_fp25 phd_fp25;

    /** \brief Muon derivatives (long double version)
     */
    part_deriv_press_ld mud_ld;
    
    /** \brief Muon derivatives (25-digit version)
     */
    part_deriv_press_fp25 mud_fp25;
    
    /** \brief Electron derivatives (long double version)
     */
    part_deriv_press_ld ed_ld;
    
    /** \brief Electron derivatives (25-digit version)
     */
    part_deriv_press_fp25 ed_fp25;
    
    /** \brief Relativistic fermion thermodynamics with derivatives
        (long double)
    */
    fermion_deriv_rel_ld fdrel_ld;
    
    /** \brief Relativistic fermion thermodynamics with derivatives
        (25-digit version)
     */
    fermion_deriv_rel_fp25 fdrel_fp25;
    
    /// \name Unit conversion objects to set the lepton masses
    //@{
    /// Long double precision unit conversion object
    convert_units<long double> cu_ld;

    /// 25-digit precision unit conversion object
    convert_units<o2fp_25> cu_fp25;
    //@}

    /** \brief Compute electron thermodynamics from the electron 
        density
    */
    virtual int electron_density(double T);
    
    /** \brief Compute electron thermodynamics from the electron 
        density (long double version)
    */
    virtual int electron_density_ld(long double T);

    /** \brief Compute electron thermodynamics from the electron 
        density (25-digit version)
    */
    virtual int electron_density_fp25(o2fp_25 T);

    /** \brief Compute particle thermodynamics from the density
     */
    template <class part_t, class part_thermo_t, class fp_t>
    int particle_density_tl(part_t &pa, part_thermo_t &pt, fp_t T) {
      
      int retx;
      
      // I find that the calculation without the rest mass is a bit more
      // stable, so we use that method and add the rest mass back in
      // later if necessary.
      bool inc_rest_mass=false;
      if (pa.inc_rest_mass) {
        
        inc_rest_mass=true;
        pa.inc_rest_mass=false;
        pa.mu-=pa.m;
      }

      if (verbose>1) {
        std::cout << "eos_leptons_multip::particle_density_tl(): "
                  << "n : ";
        std::cout.setf(std::ios::showpos);
        std::cout << dtos(pa.n,0) << std::endl;
        std::cout.unsetf(std::ios::showpos);
        pt.verbose=3;
      }
      retx=pt.pair_density(pa,T);
      if (verbose>1) {
        std::cout << "eos_leptons_multip::particle_density_tl(): "
                  << "mu: ";
        std::cout.setf(std::ios::showpos);
        std::cout << dtos(pa.mu,0) << std::endl;
        std::cout.unsetf(std::ios::showpos);
        pt.verbose=0;
      }
      
      if (inc_rest_mass) {
        pa.inc_rest_mass=true;
        pa.mu+=pa.m;
        pa.ed+=pa.m*e.n;
      }
      
      return retx;
    }

    /** \brief Template version of function to solve
     */
    template <class part_t, class part_thermo_t, class fp_t>
    int pair_density_eq_fun_tl(size_t nv,
                               const boost::numeric::ublas::vector<fp_t> &x,
                               boost::numeric::ublas::vector<fp_t> &y,
                               fp_t T, fp_t nq, part_t &pe, 
                               part_t &pmu, part_thermo_t &pt) {

      if (pde_from_density) {

        pe.n=x[0]*nq;
        int retx=part_density_tl(pe,pt,T);
        if (retx!=0) return retx;
        
      } else {
        
        pe.mu=x[0];
        
        bool inc_rest_mass=false;
        if (pe.inc_rest_mass) {
          inc_rest_mass=true;
          pe.inc_rest_mass=false;
          pe.mu-=pe.m;
        }

        pt.pair_mu(pe,T);
        
        if (inc_rest_mass) {
          pe.inc_rest_mass=true;
          pe.mu+=pe.m;
          pe.ed+=pe.n*pe.m;
        }
      }
      
      if (pe.inc_rest_mass) {
        if (pmu.inc_rest_mass) {
          pmu.mu=pe.mu;
        } else {
          pmu.mu=pe.mu-pmu.m;
        }
      } else {
        if (pmu.inc_rest_mass) {
          pmu.mu=pe.mu+pe.m;
        } else {
          pmu.mu=pe.mu+pe.m-pmu.m;
        }
      }
      
      if (pmu.inc_rest_mass) {
        pmu.inc_rest_mass=false;
        pmu.mu-=pmu.m;
        pt.pair_mu(pmu,T);
        pmu.inc_rest_mass=true;
        pmu.mu+=pmu.m;
        pmu.ed+=pmu.m*pmu.n;
      } else {
        pt.pair_mu(pmu,T);
      }
      
      y[0]=(pe.n+pmu.n-nq)/fabs(nq);
      
      return 0;
    }
    
    /** \brief Function to solve for \ref pair_density_eq()
     */
    virtual int pair_density_eq_fun(size_t nv, const ubvector &x,
                                    ubvector &y, double T, double nq);
    
    /** \brief Function to solve for \ref pair_density_eq_ld()
     */
    virtual long double pair_density_eq_ld_fun
    (long double x, long double T, long double nq);
    
    /** \brief Function to solve for \ref pair_density_eq_fp25()
     */
    virtual o2fp_25 pair_density_eq_fp25_fun
    (o2fp_25 x, o2fp_25 T, o2fp_25 nq);
    
  public:

    eos_leptons_multip();
    
    virtual ~eos_leptons_multip() {
    }
    
    /** \brief Electron in long double precision
     */
    fermion_ld eld;
    
    /** \brief Electron in 25-digit precision
     */
    fermion_fp25 efp25;

    /** \brief Muon in long double precision
     */
    fermion_ld muld;
    
    /** \brief Muon in 25-digit precision
     */
    fermion_fp25 mufp25;

    /** \brief Tau in long double precision
     */
    fermion_ld tauld;
    
    /** \brief Tau in 25-digit precision
     */
    fermion_fp25 taufp25;
    
    /** \brief Thermodynamic quantities for the full EOS (long double)
     */
    thermo_ld th_ld;

    /** \brief Thermodynamic quantities for the full EOS (25-digit)
     */
    thermo_fp25 th_fp25;

    /** \brief Photon derivatives
     */
    part_deriv_press_ld thd_ld;
    
    /** \brief Photon derivatives
     */
    part_deriv_press_fp25 thd_fp25;
    
    /** \brief Thermodynamics from the charge density in 
        weak equilibrium

        The first argument \c nq, is the total negative charge density
        including electrons (and muons if \ref include_muons is true)
        and \c T is the temperature.

        When \ref include_muons is false, this function is essentially
        equivalent to \ref pair_density() using \c nq for the electron
        density.

        The charge density should be in units of \f$
        1/\mathrm{fm}^{-3} \f$ and the temperature should be in units
        of \f$ 1/\mathrm{fm} \f$.

        The current values of the electron chemical potential
        potentials is used as initial guess. If \ref
        pde_from_density is true, then the current value
        of the electron density is also used as an initial guess.
    */
    virtual int pair_density_eq(double nq, double T);

    /** \brief Thermodynamics from the charge density in 
        weak equilibrium (long double version)
    */
    virtual int pair_density_eq_ld(long double nq, long double T);

    /** \brief Thermodynamics from the charge density in 
        weak equilibrium (25-digit version)
    */
    virtual int pair_density_eq_fp25(o2fp_25 nq, o2fp_25 T);

    /** \brief Relativistic fermion thermodynamics in long double precision
     */
    fermion_rel_ld frel_ld;
    
    /** \brief Relativistic fermion thermodynamics in 25 digit precision
     */
    fermion_rel_fp25 frel_fp25;

    /// One-dimensional solvers with different floating-point types
    //@{
    root_cern<funct,double> rc;
    root_cern<funct_ld,long double> rc_ld;
    root_cern<funct_fp25,o2fp_25> rc_fp25;
    //@}
    
  };

#endif  
}

#endif
