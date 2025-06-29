/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
#ifndef O2SCL_EOS_LEPTONS_H
#define O2SCL_EOS_LEPTONS_H

#include <o2scl/part.h>
#include <o2scl/fermion.h>
#include <o2scl/boson.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/part_deriv.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/set_multip.h>

namespace o2scl {

  /** \brief Lepton and photon EOS

      The function \ref pair_density_eq() computes the thermodynamic
      functions of the electrons and muons in in weak equilibrium,
      i.e. when \f$ \mu_e = \mu_{\mu} \f$. The inputs to this function
      are the temperature and the total charge density, i.e. the sum
      of the electron, positron, muon, and antimuon densities, There
      are two methods to solve for weak equilibrium. The first is to
      use the electron density as a variable, use the electron density
      to compute the electron chemical potential, and then use the
      muon chemical potential to compare the muon density and solve
      the equation \f$ n_e + n_{\mu} = n_q \f$. This method is used if
      \ref pde_from_density is true (the default). The second is to
      use the electron chemical potential as a variable and then
      compute the electron density from that. It is unclear which
      method is more stable or accurate. 
   */
  class eos_leptons {

  public:

    // Vector type
    typedef boost::numeric::ublas::vector<double> ubvector;
    
  protected:

    /// \name Internal functions [protected]
    //@{
    /** \brief Compute electron thermodynamics from the electron 
        density

        \note This internal function is designed to be used when the
        flag \ref fermion_rel::err_nonconv for the \ref frel object is
        set to false . If this function fails, it returns a non-zero
        value rather than calling the error handler. This allows it to
        be used inside \ref pair_density_eq_fun(). Thus this function
        is protected instead of public.
        
        Because this function uses \ref fermion_rel::pair_density(),
        the current electron chemical potential (from \ref e) is used
        as an initial guess.
        
        The temperature should be in units of \f$ 1/\mathrm{fm} \f$ .
    */
    virtual int electron_density(double T);

    /** \brief Function to solve for \ref pair_density_eq()

        \note This internal function presumes that include_muons is
        true (otherwise there's nothing to solve) and that the \ref
        fermion_rel::err_nonconv flag in the \ref frel object is set
        to false .

        The temperature should be in units of \f$ 1/\mathrm{fm} \f$ .

        The value \c x is used for the electron density or chemical
        potential, depending on the value of \ref pde_from_density .

        The value 
        \f[
        \frac{\left(n_e + n_{\mu} - n_q\right)}{|n_q|}
        \f]
        is stored in \c y.
    */
    virtual int pair_density_eq_fun(size_t nv, const ubvector &x,
                                    ubvector &y, double T, double nq);
    //@}

    /// \name Internal particle objects [protected]
    //@{
    /** \brief Relativistic fermion thermodynamics with derivatives
     */
    fermion_deriv_rel fdrel;
    //@}
    
    /// \name Accuracy control
    //@{
    /// Accuracy flag
    int accuracy;
    /// Default accuracy
    static const int acc_default=0;
    /// Improved accuracy
    static const int acc_improved=1;
    /// Use long double internally
    static const int acc_ld=2;
    /// Use 25-digit floating point numbers internally
    static const int acc_fp_25=3;
    //@}

  public:

    /// \name Constructor
    //@{
    eos_leptons();
    //@}

    virtual ~eos_leptons() {
    }
    
    /// \name Particle objects
    //@{
    /** \brief Electron
     */
    fermion e;

    /** \brief Electron derivatives
     */
    part_deriv_press ed;

    /** \brief Muon
     */
    fermion mu;

    /** \brief Tau
     */
    fermion tau;

    /** \brief Electron neutrino 
     */
    fermion nu_e;

    /** \brief Muon neutrino 
     */
    fermion nu_mu;
    
    /** \brief Tau neutrino 
     */
    fermion nu_tau;

    /** \brief Muon derivatives
     */
    part_deriv_press mud;
    
    /** \brief Photon
     */
    boson ph;

    /** \brief Photon derivatives
     */
    part_deriv_press phd;

    /** \brief Thermodynamic quantities for the full EOS
     */
    thermo th;

    /** \brief Photon derivatives
     */
    part_deriv_press thd;
    //@}

    /// \name Settings
    //@{
    /** \brief If true, use the electron density for 
        \ref pair_density_eq_fun() (default true)
     */
    bool pde_from_density;
    
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

    /** \brief If true, include derivatives (default false)
     */
    bool include_deriv;

    /// Verbosity parameter (default 0)
    int verbose;
    
    /// High accuracy
    void improved_acc() {
      frel.upper_limit_fac=40.0;
      frel.dit.tol_abs=1.0e-13;
      frel.dit.tol_rel=1.0e-13;
      frel.nit.tol_abs=1.0e-13;
      frel.nit.tol_rel=1.0e-13;
      frel.density_root.tol_rel=1.0e-10;

      fdrel.upper_limit_fac=40.0;
      fdrel.dit.tol_abs=1.0e-13;
      fdrel.dit.tol_rel=1.0e-13;
      fdrel.nit.tol_abs=1.0e-13;
      fdrel.nit.tol_rel=1.0e-13;

      accuracy=acc_improved;
      return;
    }

    /// Default accuracy
    void default_acc() {
      frel.upper_limit_fac=20.0;
      frel.dit.tol_abs=1.0e-8;
      frel.dit.tol_rel=1.0e-8;
      frel.nit.tol_abs=1.0e-8;
      frel.nit.tol_rel=1.0e-8;
      frel.dit.tol_abs=1.0e-8;
      frel.dit.tol_rel=1.0e-8;
      frel.nit.tol_abs=1.0e-8;
      frel.nit.tol_rel=1.0e-8;
      frel.density_root.tol_rel=4.0e-7;

      fdrel.upper_limit_fac=20.0;
      fdrel.dit.tol_abs=1.0e-8;
      fdrel.dit.tol_rel=1.0e-8;
      fdrel.nit.tol_abs=1.0e-8;
      fdrel.nit.tol_rel=1.0e-8;
      fdrel.dit.tol_abs=1.0e-8;
      fdrel.dit.tol_rel=1.0e-8;
      fdrel.nit.tol_abs=1.0e-8;
      fdrel.nit.tol_rel=1.0e-8;

      accuracy=acc_default;
      return;
    }

    /// Use long double arithmetic to improve precision
    void ld_acc() {
      accuracy=acc_ld;
      return;
    }
    
    /// Use 25-digit precision to improve precision
    void fp_25_acc() {
      accuracy=acc_fp_25;
      return;
    }
    //@}

    /// \name Main methods
    //@{
    /** \brief Thermodynamics from the electron and muon 
        chemical potentials

        The temperature should be in units of \f$ 1/\mathrm{fm} \f$ .
     */
    int pair_mu(double T);

    /** \brief Thermodynamics from the electron chemical potential
        in weak equilibrium

        When \ref include_muons is false, this function is essentially
        equivalent to \ref pair_mu().

        The temperature should be in units of \f$ 1/\mathrm{fm} \f$ .
    */
    int pair_mu_eq(double T);

    /** \brief Thermodynamics from the electron and muon densities

        The temperature should be in units of \f$ 1/\mathrm{fm} \f$ .
        The electron and muon densities should be stored in the \ref
        part::n field of \ref e and \ref mu (in units of \f$
        1/\mathrm{fm}^3 \f$ before calling this function.

        If \ref include_muons is false, then muons are ignored.

        The current values of the electron and muon chemical
        potentials are used as initial guesses to \ref
        fermion_rel::pair_density() for both the electrons and (if
        included) muons.

        The final total energy density, pressure, and entropy
        are stored in \ref th. 
     */
    int pair_density(double T);
    
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
    //@}

    /// \name Other objects
    //@{
    /// Solver for \ref pair_density_eq() 
    mroot_hybrids<> mh;

    /** \brief Relativistic fermion thermodynamics
     */
    fermion_rel frel;

    //@}
    
  };
  
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
    boson_cdf25 ph_cdf25;

    /** \brief Photon derivatives (long double version)
     */
    part_deriv_press_ld phd_ld;

    /** \brief Photon derivatives (25-digit version)
     */
    part_deriv_press_cdf25 phd_cdf25;

    /** \brief Muon derivatives (long double version)
     */
    part_deriv_press_ld mud_ld;
    
    /** \brief Muon derivatives (25-digit version)
     */
    part_deriv_press_cdf25 mud_cdf25;
    
    /** \brief Electron derivatives (long double version)
     */
    part_deriv_press_ld ed_ld;
    
    /** \brief Electron derivatives (25-digit version)
     */
    part_deriv_press_cdf25 ed_cdf25;
    
    /** \brief Relativistic fermion thermodynamics with derivatives
        (long double)
    */
    fermion_deriv_rel_ld fdrel_ld;
    
    /** \brief Relativistic fermion thermodynamics with derivatives
        (25-digit version)
     */
    fermion_deriv_rel_cdf25 fdrel_cdf25;
    
    /// \name Unit conversion objects to set the lepton masses
    //@{
    /// Long double precision unit conversion object
    convert_units<long double> cu_ld;

    /// 25-digit precision unit conversion object
    convert_units<cpp_dec_float_25> cu_cdf25;
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
    virtual int electron_density_cdf25(cpp_dec_float_25 T);
    
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

    /** \brief Function to solve for \ref pair_density_eq_cdf25()
     */
    virtual cpp_dec_float_25 pair_density_eq_cdf25_fun
    (cpp_dec_float_25 x, cpp_dec_float_25 T, cpp_dec_float_25 nq);
    
  public:

    eos_leptons_multip();
    
    virtual ~eos_leptons_multip() {
    }
    
    /** \brief Electron in long double precision
     */
    fermion_ld eld;
    
    /** \brief Electron in 25-digit precision
     */
    fermion_cdf25 ecdf25;

    /** \brief Muon in long double precision
     */
    fermion_ld muld;
    
    /** \brief Muon in 25-digit precision
     */
    fermion_cdf25 mucdf25;

    /** \brief Tau in long double precision
     */
    fermion_ld tauld;
    
    /** \brief Tau in 25-digit precision
     */
    fermion_cdf25 taucdf25;
    
    /** \brief Thermodynamic quantities for the full EOS (long double)
     */
    thermo_ld th_ld;

    /** \brief Thermodynamic quantities for the full EOS (25-digit)
     */
    thermo_cdf25 th_cdf25;

    /** \brief Photon derivatives
     */
    part_deriv_press_ld thd_ld;
    
    /** \brief Photon derivatives
     */
    part_deriv_press_cdf25 thd_cdf25;
    
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
    virtual int pair_density_eq_cdf25(cpp_dec_float_25 nq,
                                      cpp_dec_float_25 T);

    /** \brief Relativistic fermion thermodynamics in long double precision
     */
    fermion_rel_ld frel_ld;
    
    /** \brief Relativistic fermion thermodynamics in 25 digit precision
     */
    fermion_rel_cdf25 frel_cdf25;

    /// One-dimensional solvers with different floating-point types
    //@{
    root_cern<funct,double> rc;
    root_cern<funct_ld,long double> rc_ld;
    root_cern<funct_cdf25,cpp_dec_float_25> rc_cdf25;
    //@}
    
  };

#endif  
}

#endif
