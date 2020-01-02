/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
/** \file nstar_cold.h
    \brief File defining \ref o2scl::nstar_cold
*/
#ifndef O2SCL_COLD_NSTAR_H
#define O2SCL_COLD_NSTAR_H

#include <o2scl/eos_had_base.h>
#include <o2scl/tov_solve.h>
#include <o2scl/tov_solve.h>
#include <o2scl/table.h>
#include <o2scl/fermion.h>
#include <o2scl/root_cern.h>
#include <o2scl/mroot_cern.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/eos_tov.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Naive static cold neutron star

      This uses eos_had_base::calc_e() to compute the equation of 
      state of zero-temperature beta-equilibrated neutron star
      matter and tov_solve::mvsr() to compute the mass versus
      radius curve. By default the crust EOS is given by
      that in \ref o2scl::eos_tov_interp::default_low_dens_eos() .

      The neutron, proton, electron and muon are given masses
      according to their values in \ref o2scl_mks
      after a conversion to units of \f$ 1/\mathrm{fm} \f$.

      There is an example for the usage of this class given
      in the \ref ex_nstar_cold_sect.

      If \ref err_nonconv is true and the solver fails,
      the error handler is called. 
      
      \hline
      \b EOS \b Output

      The function calc_eos() generates an object of type
      \ref table_units, which contains the following columns
      - \c ed in units of \f$ 1/\mathrm{fm}^4 \f$, the total energy 
      density of neutron star matter, \f$ \varepsilon \f$
      - \c pr in units of \f$ 1/\mathrm{fm}^4 \f$, the total pressure
      of neutron star matter, \f$ P \f$
      - \c nb in units of \f$ 1/\mathrm{fm}^3 \f$, the baryon 
      number density, \f$ n_B \f$
      - \c mun in units of \f$ 1/\mathrm{fm} \f$, the neutron 
      chemical potential, \f$ \mu_n \f$
      - \c mup in units of \f$ 1/\mathrm{fm} \f$, the proton chemical
      potential, \f$ \mu_p \f$
      - \c mue in units of \f$ 1/\mathrm{fm} \f$, the electron 
      chemical potential, \f$ \mu_e \f$
      - \c nn in units of \f$ 1/\mathrm{fm}^3 \f$, the neutron number
      density, \f$ n_n \f$
      - \c np in units of \f$ 1/\mathrm{fm}^3 \f$, the proton number
      density, \f$ n_p \f$
      - \c ne in units of \f$ 1/\mathrm{fm}^3 \f$, the electron 
      number density, \f$ n_e \f$
      - \c kfn in units of \f$ 1/\mathrm{fm} \f$, the neutron Fermi
      momentum, \f$ k_{F,n} \f$
      - \c kfp in units of \f$ 1/\mathrm{fm} \f$, the proton Fermi
      momentum, \f$ k_{F,p} \f$
      - \c kfe in units of \f$ 1/\mathrm{fm} \f$, the electron 
      Fermi momentum, \f$ k_{F,e} \f$
      - \c dednb_Ye in units of \f$ 1/\mathrm{fm} \f$, 
      \f$ ( d \varepsilon / d n_B )_{Y_e} \f$ where
      \f$ Y_e = n_e / n_B \f$ is the electron fraction
      (computed using \ref eos_had_base::const_pf_derivs() )
      - \c dPdnb_Ye in units of \f$ 1/\mathrm{fm} \f$, 
      \f$ ( d P / d n_B )_{Y_e} \f$ .
      - \c fcs2, the squared speed of sound at fixed electron
      fraction, the ratio of the previous two quantities

      If \ref include_muons is true, the table has 
      additional columns
      - \c mumu in units of \f$ 1/\mathrm{fm} \f$, the muon chemical
      potential, \f$ \mu_{\mu} \f$
      - \c nmu in units of \f$ 1/\mathrm{fm}^3 \f$, the muon number
      density, \f$ n_{\mu} \f$
      - \c kfmu in units of \f$ 1/\mathrm{fm} \f$, the muon Fermi 
      momentum, \f$ k_{F,\mu} \f$

      If the energy density is always positive and increasing, and the
      pressure is always positive and increasing, then the EOS is
      well-formed and \ref well_formed is \c true. The variable \ref
      pressure_dec records the lowest baryon density where the
      pressure decreases with increasing density. 
      If \ref err_nonconv is true and the EOS is not well
      formed, the error handler is called, and the remaining
      columns below are not computed. 
      
      After computing the equation of state, \ref calc_eos()
      also adds the following columns
      - \c cs2 (unitless), the squared speed of sound divided by \f$ c^2 \f$
      - \c logp, the natural logarithm of the pressure stored in \c pr
      - \c loge, the natural logarithm of the energy density 
      stored in \c ed
      - \c s in units of \f$ 1/\mathrm{fm} \f$, 
      the semi-perimeter of the Urca triangle
      - \c urca in units of \f$ 1/\mathrm{fm}^4 \f$, 
      the squared area of the Urca triangle
      - \c ad_index, the adiabatic index, \f$ \Gamma \f$
      If the eos is not well-formed and \ref well_formed is <tt>false</tt>,
      then the columns <tt>cs2</tt>, <tt>logp</tt>, and <tt>loge</tt> 
      are set to zero. The columns \c cs2 and \c ad_indes are computing
      from derivatives using the current table interpolation type. 

      The condition for the direct Urca process is the area of the
      triangle formed by the neutron, proton, and electron Fermi
      momenta. Using the definition of the semi-perimeter, 
      \f[
      s \equiv \left( k_{F,n}+k_{F,p}+k_{F,e} \right)/2
      \f]
      Heron's formula gives the triangle area as
      \f[
      a=\sqrt{s(s-k_{F,n})(s-k_{F,p})(s-k_{F,e})} \, .
      \f]
      The column in the eos \ref table labeled \c urca is \f$ a^2 \f$
      . If this quantity is positive, then direct Urca is allowed. The
      variable \ref allow_urca is the smallest density for which the
      direct Urca process turns on, and \ref deny_urca is the smallest
      density for which the direct Urca process turns off.
      
      The squared speed of sound (in units of  \f$ c \f$ )
      is calculated by 
      \f[
      c_s^2 = \frac{ d P }{d \varepsilon}
      \f]
      and this is placed in the column labeled \c cs2. If the
      EOS is not well-formed, then this column is set to zero. If 
      \c cs2 is larger than 1, the EOS is said to be "acausal". The
      variables \ref acausal, \ref acausal_ed, and \ref acausal_pr
      record the baryon density, energy density, and pressure where
      the EOS becomes acausal. The adabatic index is calculated by
      \f[
      \Gamma = \frac{ d \ln P} { d \ln \varepsilon}
      \f]
      Note that \f$ \Gamma \f$ must be greater than \f$ 4/3 \f$ 
      at the center of the neutron star for stability. (This
      is a necessary, but not sufficient condition.)  If
      the EOS is not well-formed then this column is set to zero.

      \hline
      \b TOV \b Output

      The TOV table contains all the columns typically 
      generated for mass versus radius tables in \ref tov_solve,
      as well as columns containing the central values of 
      al the densities and chemical potentials, and all the 
      other columns computed for the EOS above. 

      \hline

      \future Warn if the EOS becomes pure neutron matter.
      \future Some of the auxillary quantities can be computed
      directly without using the table methods and the 
      EOS calculation would be a bit faster.
  */

  class nstar_cold {

  public:

    nstar_cold();

    /// \name Basic operation
    //@{
    /** \brief Set the equation of state
	
        This should be set before calling calc_eos().
    */
    void set_eos(eos_had_base &he) {
      hep=&he;
      eos_set=true;
      return;
    }

    /** \brief Calculate the given equation of state
     */
    int calc_eos(double np_0=0.0);
    
    /** \brief Compute the density at which the direct Urca process is allowe

	This is faster than using calc_eos() since it does nothing
	other than computes the critical density. It does not store
	the equation of state.
     */
    double calc_urca(double np_0=0.0);

    /** \brief Calculate the M vs. R curve
     */
    int calc_nstar();

    /** \brief Calculate the profile for a fixed gravitational mass
     */
    int fixed(double target_mass);
    //@}

    /// \name Output 
    //@{
    /** \brief If true, the energy density and pressure 
	of the EOS is monotonically 
	increasing and the pressure is always positive
    */
    bool well_formed;

    /** \brief The smallest baryon density where the pressure starts 
	to decrease

	If this is zero after calling calc_eos(), then 
	the pressure does not decrease in the specified range
	of baryon density
     */
    double pressure_dec;

    /** \brief The smallest density where Urca becomes allowed

        If this is zero after calling calc_eos(), then direct
	Urca is never allowed.
     */
    double allow_urca;

    /** \brief The smallest density where Urca becomes disallowed
	
        If this is zero after calling calc_eos(), then direct
	Urca is not disallowed at a higher density than 
	it becomes allowed.
     */
    double deny_urca;

    /** \brief The density at which the EOS becomes acausal

        If this is zero, then the EOS is causal at all baryon densities
	in the specified range
     */
    double acausal;

    /** \brief The pressure at which the EOS becomes acausal

        If this is zero, then the EOS is causal at all baryon densities
	in the specified range
     */
    double acausal_pr;

    /** \brief The energy density at which the EOS becomes acausal

        If this is zero, then the EOS is causal at all baryon densities
	in the specified range
     */
    double acausal_ed;

    /** \brief Solver tolerance (default \f$ 10^{-4} \f$)
    */
    double solver_tol;

    /// Verbosity parameter (default 0)
    int verbose;

    /** \brief Set the EOS table

	In order for the \ref calc_nstar() function to use this 
	table, it must
	contain at least the columns <tt>ed, pr</tt>, and
	<tt>nB</tt> which store the energy density, pressure,
	and baryon density. 
     */
    void set_eos_table(std::shared_ptr<table_units<> > t) {
      eost=t;
      return;
    }
    
    /** \brief Get the eos table (after having called calc_eos())
     */
    std::shared_ptr<table_units<> > get_eos_results() {
      return eost;
    }
    
    /** \brief Get the results from the TOV (after having called calc_nstar())
     */
    std::shared_ptr<table_units<> > get_tov_results() {
      return tp->get_results();
    }
    //@}

    /** \name Configuration
     */
    //@{
    /** \brief The starting baryon density (default 0.05)
     */
    double nb_start;

    /** \brief The final baryon density (default 2.0)
     */
    double nb_end;

    /** \brief The baryon density stepsize (default 0.01)
     */
    double dnb;

    /** \brief If true, include muons (default false)
     */
    bool include_muons;

    /** \brief If true, throw an exception if the solver fails
	or if the EOS is not well-formed (default true)
    */
    bool err_nonconv;

    /** \brief Set the equation solver for the EOS
     */
    int set_root(root<> &rf) {
      rp=&rf;
      return 0;
    }

    /** \brief Specify the object for solving the TOV equations
	
        The default uses the low-density equation of state with
	tov::verbose=0. In calc_nstar(), the units are set by calling
	tov_solve::set_units().
    */
    int set_tov(tov_solve &ts) {
      tp=&ts;
      return 0;
    }
    //@}

    /** \name Default objects */
    //@{
    /// The default neutron
    fermion np;

    /// The default proton
    fermion pp;

    /// Zero-temperature fermion thermodynamics
    fermion_zerot fzt;

    /** \brief The default TOV equation solver
     */
    tov_solve def_tov;

    /** \brief The default equation solver for the EOS
    */
    root_cern<> def_root;

    /// Default EOS object for the TOV solver
    eos_tov_interp def_eos_tov;
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Solve to ensure zero charge in \f$ \beta \f$-equilibrium
    double solve_fun(double x, thermo &hb);

    /// True if equation of state has been set
    bool eos_set;

    /// The electron
    fermion e;

    /// The muon
    fermion mu;

    /// A pointer to the equation of state
    eos_had_base *hep;

    /// A pointer to the TOV object
    tov_solve *tp;
    
    /// A pointer to the solver
    root<> *rp;

    /// Storage for the EOS table
    std::shared_ptr<table_units<> > eost;

    /// The baryon density
    double barn;

#endif

  };


#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
