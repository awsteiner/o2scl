/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_EOS_H
#define O2SCL_EOS_H

#include <o2scl/part.h>
#include <o2scl/fermion.h>
#include <o2scl/boson.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/part_deriv.h>
#include <o2scl/fermion_deriv_rel.h>

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

    /** \brief Electron thermodynamics from the electron density

        \note This internal function is designed to be used when \ref
        fermion_rel::err_nonconv is false . If this function fails,
        returns a non-zero value rather than calling the error
        handler. This allows it to be used inside \ref
        pair_density_eq_fun(). This is also the reason why this
        function is protected instead of public.
        
        Because this function uses \ref fermion_rel::pair_density(),
        the current electron chemical potential is used as an initial
        guess.

        The temperature should be in units of \f$ 1/\mathrm{fm} \f$ .
    */
    int electron_density(double T);

    /** \brief Function to solve for \ref pair_density_eq()

        \note This internal function presumes that include_muons is
        true (otherwise there's nothing to solve) and that \ref
        fermion_rel::err_nonconv is false .

        The temperature should be in units of \f$ 1/\mathrm{fm} \f$ .
     */
    int pair_density_eq_fun(size_t nv, const ubvector &x,
                            ubvector &y, double T, double nq);

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

    /** \brief Relativistic fermion thermodynamics with derivatives
     */
    fermion_deriv_rel fdrel;
    
    /// \name Accuracy control
    //@{
    int accuracy;
    static const int acc_default=0;
    static const int acc_improved=1;
    static const int acc_ld=2;
    static const int acc_fp_25=3;
    //@}

    /// Long double precision unit conversion object
    convert_units<long double> cu_ld;

    /// 25-digit precision unit conversion object
    convert_units<cpp_dec_float_25> cu_cdf25;
    
  public:

    /// \name Constructor
    //@{
    eos_leptons();
    //@}

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
      frel.fri.dit.tol_abs=1.0e-13;
      frel.fri.dit.tol_rel=1.0e-13;
      frel.fri.nit.tol_abs=1.0e-13;
      frel.fri.nit.tol_rel=1.0e-13;
      frel.fri.dit.tol_abs=1.0e-13;
      frel.fri.dit.tol_rel=1.0e-13;
      frel.fri.nit.tol_abs=1.0e-13;
      frel.fri.nit.tol_rel=1.0e-13;
      frel.density_root->tol_rel=1.0e-10;

      fdrel.upper_limit_fac=40.0;
      fdrel.def_dit.tol_abs=1.0e-13;
      fdrel.def_dit.tol_rel=1.0e-13;
      fdrel.def_nit.tol_abs=1.0e-13;
      fdrel.def_nit.tol_rel=1.0e-13;
      fdrel.def_dit.tol_abs=1.0e-13;
      fdrel.def_dit.tol_rel=1.0e-13;
      fdrel.def_nit.tol_abs=1.0e-13;
      fdrel.def_nit.tol_rel=1.0e-13;
      fdrel.def_density_root.tol_rel=1.0e-10;

      accuracy=acc_improved;
      return;
    }

    /// Default accuracy
    void default_acc() {
      frel.upper_limit_fac=20.0;
      frel.fri.dit.tol_abs=1.0e-8;
      frel.fri.dit.tol_rel=1.0e-8;
      frel.fri.nit.tol_abs=1.0e-8;
      frel.fri.nit.tol_rel=1.0e-8;
      frel.fri.dit.tol_abs=1.0e-8;
      frel.fri.dit.tol_rel=1.0e-8;
      frel.fri.nit.tol_abs=1.0e-8;
      frel.fri.nit.tol_rel=1.0e-8;
      frel.density_root->tol_rel=4.0e-7;

      fdrel.upper_limit_fac=20.0;
      fdrel.def_dit.tol_abs=1.0e-8;
      fdrel.def_dit.tol_rel=1.0e-8;
      fdrel.def_nit.tol_abs=1.0e-8;
      fdrel.def_nit.tol_rel=1.0e-8;
      fdrel.def_dit.tol_abs=1.0e-8;
      fdrel.def_dit.tol_rel=1.0e-8;
      fdrel.def_nit.tol_abs=1.0e-8;
      fdrel.def_nit.tol_rel=1.0e-8;
      fdrel.def_density_root.tol_rel=4.0e-7;

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
    int pair_density_eq(double nq, double T);
    //@}

    /// \name Other objects
    //@{
    /// Solver for \ref pair_density_eq() 
    mroot_hybrids<> mh;

    /** \brief Relativistic fermion thermodynamics
     */
    fermion_rel frel;

    /** \brief Relativistic fermion thermodynamics in long double precision
     */
    fermion_rel_ld frel_ld;
    
    /** \brief Relativistic fermion thermodynamics in 25 digit precision
     */
    fermion_rel_cdf25 frel_cdf25;
    //@}
    
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
                           eos_leptons &elep,
			   //fermion &e, bool include_muons,
			   //fermion &mu, fermion_rel &frel,
			   std::shared_ptr<table_units<> > results);
    
  protected:

    /// A pointer to the thermo object
    thermo *eos_thermo;

  };

  
}

#endif
