/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#ifndef O2SCL_NSTAR_COLD_H
#define O2SCL_NSTAR_COLD_H

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

  /** \brief Static neutron star at zero temperature

      This class provides a simplified interface for computing the
      equation of state (EOS) of zero-temperature beta-equilibrated
      neutron star matter and computing the mass-radius curve.

      The beta-equilibrium EOS can be computed with the function \ref
      calc_eos() and the mass-radius curve with the function \ref
      calc_nstar(). The function \ref fixed() computes the profile of
      a fixed gravitational mass star. By default the crust EOS is
      given by that in \ref
      o2scl::eos_tov_interp::default_low_dens_eos() .

      There is an example for the usage of this class given
      in the :ref:`Cold Neutron Star Structure` section
      of the User's Guide.

      The neutron, proton, electron and muon are given masses
      according to their values in \ref o2scl_mks
      after a conversion to units of \f$ 1/\mathrm{fm} \f$.

      \note The function \ref set_eos() stores a pointer to the EOS
      object so the user must take care that the pointer is valid.

      If \ref err_nonconv is true (the default) and the solver fails,
      the error handler is called. 

      The baryon density range is specified by \ref nb_start, \ref nb_end
      and \ref dnb. Given the value of \f$ f \f$ defined by
      \f[
      f = \frac{\mathrm{nb\_start} - \mathrm{nb\_end}}/{\mathrm{dnb}}
      \f]
      if \f$ f<1 \f$ or \f$ f>10^8 \f$, then the error handler is 
      called.

      Note that not all equations of state lead to physical neutron
      stars, and the function \ref calc_eos() will not warn you about
      this. However, the function \ref calc_nstar() will call the
      error handler if the speed of sound becomes negative or becomes
      larger than the speed of light at a density smaller than the
      central density of the maximum mass star.
      
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

      \comment
      - \c dednb_Ye in units of \f$ 1/\mathrm{fm} \f$, 
      \f$ ( d \varepsilon / d n_B )_{Y_e} \f$ where
      \f$ Y_e = n_e / n_B \f$ is the electron fraction
      (computed using \ref eos_had_base::const_pf_derivs() )
      - \c dPdnb_Ye in units of \f$ 1/\mathrm{fm} \f$, 
      \f$ ( d P / d n_B )_{Y_e} \f$ .
      - \c fcs2, the squared speed of sound at fixed electron
      fraction, the ratio of the previous two quantities
      \endcomment

      If \ref include_muons is true, the table has 
      additional columns
      - \c mumu in units of \f$ 1/\mathrm{fm} \f$, the muon chemical
      potential, \f$ \mu_{\mu} \f$
      - \c nmu in units of \f$ 1/\mathrm{fm}^3 \f$, the muon number
      density, \f$ n_{\mu} \f$
      - \c kfmu in units of \f$ 1/\mathrm{fm} \f$, the muon Fermi 
      momentum, \f$ k_{F,\mu} \f$
      
      After computing the equation of state,
      then \ref calc_eos() also adds the following columns
      - \c cs2 (unitless), the squared speed of sound divided by \f$ c^2 \f$
      - \c logp, the natural logarithm of the pressure stored in \c pr
      - \c loge, the natural logarithm of the energy density 
      stored in \c ed
      - \c s in units of \f$ 1/\mathrm{fm} \f$, 
      the semi-perimeter of the Urca triangle
      - \c urca in units of \f$ 1/\mathrm{fm}^4 \f$, 
      the squared area of the Urca triangle
      - \c ad_index, the adiabatic index, \f$ \Gamma \f$

      The columns \c cs2 and \c ad_index are computing
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
      variable \ref allow_urca_nb is the smallest density for which the
      direct Urca process turns on, and \ref deny_urca_nb is the smallest
      density for which the direct Urca process turns off.
      
      The squared speed of sound (in units of \f$ c \f$ )
      is calculated by 
      \f[
      c_s^2 = \frac{ d P }{d \varepsilon}
      \f]
      and this is placed in the column labeled \c cs2. If 
      \c cs2 is larger than 1, the EOS is said to be "acausal". The
      variables \ref acausal_nb, \ref acausal_ed, and \ref acausal_pr
      record the baryon density, energy density, and pressure where
      the EOS becomes acausal. The adabatic index is calculated by
      \f[
      \Gamma = \frac{ d \ln P} { d \ln \varepsilon}
      \f]
      Note that \f$ \Gamma \f$ must be greater than \f$ 4/3 \f$ 
      at the center of the neutron star for stability. (This
      is a necessary, but not sufficient condition.)  

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

    typedef boost::numeric::ublas::vector<double> ubvector;
    
    /// \name Basic operation
    //@{
    /** \brief Set the equation of state
	
        This function stores a pointer to the EOS object and
        should be called before calling calc_eos().
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
    /** \brief If true, the pressure or energy density became negative
        at some point
     */
    bool eos_neg;

    /** \brief The smallest baryon density where the pressure starts 
	to decrease

	If this is zero after calling calc_eos(), then 
	the pressure does not decrease in the specified range
	of baryon density
     */
    double pressure_dec_nb;

    /** \brief The smallest density where Urca becomes allowed

        If this is zero after calling calc_eos(), then direct
	Urca is never allowed.
     */
    double allow_urca_nb;

    /** \brief The smallest density where Urca becomes disallowed
	
        If this is zero after calling calc_eos(), then direct
	Urca is not disallowed at a higher density than 
	it becomes allowed.
     */
    double deny_urca_nb;

    /** \brief The density at which the EOS becomes acausal

        If this is zero, then the EOS is causal at all baryon densities
	in the specified range
     */
    double acausal_nb;

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

    /// The table row which contains the maximum gravitational mass
    size_t max_row;
    
    /** \brief Set the EOS table

	In order for the \ref calc_nstar() function to use this table,
	it must contain at least the columns <tt>ed, pr</tt>, and
	<tt>nB</tt> which store the energy density, pressure, and
	baryon density.
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
    /// Verbosity parameter (default 0)
    int verbose;

    /** \brief If true, remove rows beyond the maximum mass (default true)

        Note that if the M-R curve has multiple branches, this class will
        not remove all unstable configurations.
     */
    bool remove_rows;
    
    /** \brief The starting baryon density (default 0.05)
     */
    double nb_start;

    /** \brief The final baryon density (default 2.0)
     */
    double nb_end;

    /** \brief The baryon density stepsize (default 0.01)
     */
    double dnb;

    /** \brief If true, include muons (default true)
     */
    bool include_muons;

    /** \brief If true, throw an exception if the solver fails
	or if the EOS is unphysical (default true)
    */
    bool err_nonconv;

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
    fermion neut;

    /// The default proton
    fermion prot;

    /// Zero-temperature fermion thermodynamics (for the leptons)
    fermion_zerot fzt;

    /** \brief The default TOV equation solver
     */
    tov_solve def_tov;

    /// Default EOS object for the TOV solver
    eos_tov_interp def_eos_tov;
    
    /** \brief The EOS solver

        The value of \ref mroot_hybrids::err_nonconv is set to false
        in the constructor so that this class can manage solver
        convergence errors.
    */
    mroot_hybrids<> mh;
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Solve to ensure zero charge in \f$ \beta \f$-equilibrium

        This function takes as input (in <tt>x[0]</tt>) the proton
        density and outputs (in <tt>y[0]</tt>) returns the net charge
        density in beta-equilibrium, i.e.
        \f[
        n_p - n_e - n_{\mu}
        \f]

        This function returns a non-zero value if the densities are
        negative or if the EOS fails.
    */
    int solve_fun(size_t nv, const ubvector &x, ubvector &y,
                  thermo &hb, double n_B);
    
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
    
    /// Storage for the EOS table
    std::shared_ptr<table_units<> > eost;

#endif

  };

  /** \brief Neutron stars at finite temperature and entropy
   */
  class nstar_hot : public nstar_cold {

  protected:
    
    typedef boost::numeric::ublas::vector<double> ubvector;

    /// A pointer to the finite temperature equation of state
    eos_had_temp_base *hepT;

    /** \brief Fermion thermodynamics (for the leptons)
     */
    fermion_rel ft;

    /** \brief Solve for beta equilibrium at finite temperature
     */
    int solve_fun_T(size_t nv, const ubvector &x,
                     ubvector &y, thermo &hb, double T,
                     double n_B);
    
    /** \brief Solve for beta equilibrium at finite temperature
     */
    int solve_fun_s(size_t nv, const ubvector &x, ubvector &y,
                    thermo &hb, double s, double n_B);

    /** \brief Solve for beta equilibrium at finite temperature
     */
    int solve_fun_s_YLe(size_t nv, const ubvector &x, ubvector &y,
                        thermo &hb, double s, double YLe, double n_B);

    /** \brief If true, then the hadronic EOS has been set
     */
    bool eos_T_set;
    
  public:

    nstar_hot() {
      eos_T_set=false;
    }
    
    /** \brief Set the finite-temperature hadronic EOS
     */
    void set_eos_T(eos_had_temp_base &he) {
      hepT=&he;
      eos_T_set=true;
      return;
    }
    
    /** \brief Compute the EOS in betq-equilibrium at finite temperature
     */
    int calc_eos_T(double T, double np_0=0.0);
    
    /** \brief Compute the EOS in betq-equilibrium at fixed entropy per baryon
     */
    int calc_eos_s(double s, double np_0=0.0);
    
    /** \brief Compute the EOS in betq-equilibrium at fixed entropy 
	per baryon at a fixed number of electron-type leptons per
	baryon
    */
    int calc_eos_s_YLe(double s, double YLe, double np_0=0.0);
    
  };
  

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
