/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_REL_FERMION_H
#define O2SCL_REL_FERMION_H

/** \file fermion_rel.h
    \brief File defining \ref o2scl::fermion_rel
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/shared_ptr.h>
#include <o2scl/fermion.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a relativistic fermion

      This class computes the thermodynamics of a relativistic fermion
      either as a function of the density or the chemical potential.
      It employs direct integration, using two different integrators
      for the degenerate and non-degenerate regimes. The default
      integrators are inte_qag_gsl (for degenerate fermions) and
      inte_qagiu_gsl (for non-degenerate fermions). For the functions
      calc_mu() and calc_density(), if the temperature argument is
      less than or equal to zero, the functions \ref
      fermion_zerot::calc_mu_zerot() and \ref
      fermion_zerot::calc_density_zerot() will be used to compute the
      result.

      \hline 
      <b>Degeneracy parameter:</b>

      Define the degeneracy parameter 
      \f[
      \psi=(\nu-m^{*})/T 
      \f] 
      where \f$ \nu \f$ is the effective chemical potential (including
      the rest mass) and \f$
      m^{*} \f$ is the effective mass. For \f$ \psi \f$ smaller than
      \ref min_psi, the non-degenerate expansion in \ref
      fermion_eval_thermo::calc_mu_ndeg() is attempted first. If that
      fails, then integration is used. For \f$ \psi \f$ greater than
      \ref deg_limit (degenerate regime), a finite interval integrator
      is used and for \f$ \psi \f$ less than \ref deg_limit
      (non-degenerate regime), an integrator over the interval from
      \f$ [0,\infty) \f$ is used. In the case where \ref
      part::inc_rest_mass is false, the degeneracy parameter is
      \f[
      \psi=(\nu+m-m^{*})/T 
      \f] 

      <b>Integration limits:</b>

      The upper limit on the degenerate integration is given by
      \f[
      \mathrm{upper~limit} = \sqrt{{\cal L}^2-m^{*,2}}
      \f]
      where \f$ {\cal L}\equiv u T+\nu \f$ and \f$ u \f$ is \ref
      fermion_rel::upper_limit_fac . In the case where \ref
      part::inc_rest_mass is false, the result is
      \f[
      \mathrm{upper~limit} = \sqrt{(m+{\cal L})^2-m^{*2}}
      \f]
      
      The entropy is only significant at the Fermi surface, thus
      in the degenerate case, the lower limit of the entropy
      integral can be given be determined by the value of \f$ k \f$ 
      which solves
      \f[
      - u = \frac{\sqrt{k^2+m^{* 2}}-\nu}{T}
      \f]
      The solution is 
      \f[
      \mathrm{lower~limit} = \sqrt{(-u T+{\nu})^2-m^{*,2}}
      \f]
      but this solution is only valid if \f$ (m^{*}-\nu)/T < -u \f$.
      In the case where part::inc_rest_mass is false, the result is
      \f[
      \mathrm{lower~limit} = \sqrt{(-u T + m +\nu)^2-m^{*,2}}
      \f]
      which is valid if \f$ (m^{*}-\nu - m)/T < -u \f$.

      <b>Entropy integrand:</b>

      In the degenerate regime, the entropy integrand
      \f[
      - k^2 \left[ f \log f + \left(1-f\right) \log 
      \left(1-f \right) \right]
      \f]
      where \f$ f \f$ is the fermionic distribution function can lose
      precision when \f$ (E^{*} - \nu)/T \f$ is negative and
      sufficiently large in absolute magnitude. Thus when \f$ (E^{*} -
      \nu)/T < S \f$ where \f$ S \f$ is stored in \ref deg_entropy_fac
      (default is -30), the integrand is written as
      \f[
      -k^2 \left( E/T-\nu/T \right) e^{E/T-\nu/T} \, .
      \f]
      If \f$ (E - \nu)/T < S \f$ is less than -1 times \ref exp_limit
      (e.g. less than -200), then the entropy integrand is assumed 
      to be zero.
      
      <b>Non-degenerate integrands:</b>
      
      \comment
      It's not at all clear that this dimensionless form is more
      accurate than other potential alternatives. On the other hand,
      it seems that the uncertainties in the integrations are larger
      than the errors made by the integrand at present.
      \endcomment
      The integrands in the non-degenerate regime are written
      in a dimensionless form, by defining \f$ u \f$ with
      the relation
      \f$ p = \sqrt{\left(T u + m^{*}\right)^2-m^{* 2}} \f$,
      \f$ y \equiv \nu/ T \f$, and 
      \f$ \mathrm{mx} \equiv m^{*}/T \f$. 
      The density integrand is 
      \f[
      \left(\mathrm{mx}+u\right) \sqrt{u^2+2 (\mathrm{mx}) u}
      \left(\frac{e^{y}}{e^{\mathrm{mx}+u}+e^{y}}\right) \, , 
      \f]
      the energy integrand is 
      \f[
      \left(\mathrm{mx}+u\right)^2 \sqrt{u^2+2 (\mathrm{mx}) u}
      \left(\frac{e^{y}}{e^{\mathrm{mx}+u}+e^{y}}\right) \, ,
      \f]
      and the entropy integrand is 
      \f[
      \left(\mathrm{mx}+u\right) \sqrt{u^2+2 (\mathrm{mx}) u} 
      \left(t_1+t_2\right) \, ,
      \f]
      where 
      \f{eqnarray*}
      t_1 &=& \log \left(1+e^{y-\mathrm{mx}-u}\right)/
      \left(1+e^{y-\mathrm{mx}-u}\right) \nonumber \\
      t_2 &=& \log \left(1+e^{\mathrm{mx}+u-y}\right)/
      \left(1+e^{\mathrm{mx}+u-y}\right) \, .
      \f}

      \hline 
      <b>Accuracy:</b>

      The default settings for for this class give an accuracy of at
      least 1 part in \f$ 10^6 \f$ (and frequently better than this).

      When the integrators provide numerical uncertainties, these
      uncertainties are stored in \ref unc. In the case of
      calc_density() and pair_density(), the uncertainty from the
      numerical accuracy of the solver is not included. (There is also
      a relatively small inaccuracy due to the mathematical evaluation
      of the integrands which is not included in \ref unc.)
     
      One can improve the accuracy to within 1 part in \f$ 10^{10} \f$ 
      using
      \code
      fermion_rel rf(1.0,2.0);
      rf.upper_limit_fac=40.0;
      rf.dit->tol_abs=1.0e-13;
      rf.dit->tol_rel=1.0e-13;
      rf.nit->tol_abs=1.0e-13;
      rf.nit->tol_rel=1.0e-13;
      rf.density_root->tol_rel=1.0e-10;
      \endcode
      which decreases the both the relative and absolute tolerances
      for both the degenerate and non-degenerate integrators and
      improves the accuracy of the solver which determines the
      chemical potential from the density. Of course if these
      tolerances are too small, the calculation may fail.

      \hline 
      <b>Todos:</b>

      \future The expressions which appear in in the integrand
      functions density_fun(), etc. could likely be improved,
      especially in the case where \ref o2scl::part::inc_rest_mass is
      <tt>false</tt>. There should not be a need to check if
      <tt>ret</tt> is finite.

      \future It appears this class doesn't compute the uncertainty in
      the chemical potential or density with calc_density(). This
      could be fixed.

      \future I'd like to change the lower limit on the entropy 
      integration, but the value in the code at the moment (stored
      in <tt>ll</tt>) makes bm_part2.cpp worse.

      \future The function pair_mu() should set the antiparticle
      integrators as done in fermion_deriv_rel.
  */
  class fermion_rel : public fermion_eval_thermo {

  public:

    /// \name Numerical parameters
    //@{
    /** \brief If true, call the error handler when convergence 
	fails (default true)
    */
    bool err_nonconv;

    /** \brief The smallest value of \f$ (\mu-m)/T \f$ for which 
	integration is used
     */
    double min_psi;

    /** \brief The critical degeneracy at which to switch integration 
	techniques (default 2)
    */
    double deg_limit;
    
    /** \brief The limit for exponentials to ensure integrals are finite 
	(default 200)
    */
    double exp_limit;

    /// The factor for the degenerate upper limits (default 20)
    double upper_limit_fac;

    /// A factor for the degenerate entropy integration (default 30)
    double deg_entropy_fac;
    //@}

    /// Storage for the uncertainty
    fermion unc;

    /// If true, use expansions for extreme conditions (default true)
    bool use_expansions;

    /// Create a fermion with mass \c m and degeneracy \c g
    fermion_rel();

    virtual ~fermion_rel();
    
    /** \brief Calculate properties as function of chemical potential
    */
    virtual void calc_mu(fermion &f, double temper);

    /** \brief Calculate properties as function of density

        This function uses the current value of \c nu (or \c mu if the
	particle is non interacting) for an initial guess to solve for
	the chemical potential. If this guess is too small, then this
	function may fail.
     */
    virtual int calc_density(fermion &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	chemical potential
    */
    virtual void pair_mu(fermion &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	density
     */
    virtual int pair_density(fermion &f, double temper);

    /** \brief Calculate effective chemical potential from density

	\future This function might be improved by generating a
	bracket for a bracketing solver, rather than \ref
	o2scl::root_cern which is the default for \ref
	o2scl::fermion_rel::density_root.
     */
    virtual int nu_from_n(fermion &f, double temper);
    
    /// The non-degenerate integrator
    o2_shared_ptr<inte<> >::type nit;

    /// The degenerate integrator
    o2_shared_ptr<inte<> >::type dit;

    /// The solver for calc_density()
    o2_shared_ptr<root<> >::type density_root;

    /// The backup solver for calc_density()
    o2_shared_ptr<root<> >::type density_root2;
    
    /// Return string denoting type ("fermion_rel")
    virtual const char *type() { return "fermion_rel"; }

  protected:
    
#ifndef DOXYGEN_INTERNAL
    
    /// The integrand for the density for non-degenerate fermions
    double density_fun(double u, fermion &f, double T);

    /// The integrand for the energy density for non-degenerate fermions
    double energy_fun(double u, fermion &f, double T);

    /// The integrand for the entropy density for non-degenerate fermions
    double entropy_fun(double u, fermion &f, double T);

    /// The integrand for the density for degenerate fermions
    double deg_density_fun(double u, fermion &f, double T);

    /// The integrand for the energy density for degenerate fermions
    double deg_energy_fun(double u, fermion &f, double T);

    /// The integrand for the entropy density for degenerate fermions
    double deg_entropy_fun(double u, fermion &f, double T);

    /// Solve for the chemical potential given the density
    double solve_fun(double x, fermion &f, double T);

    /** \brief Solve for the chemical potential given the density 
	with antiparticles
	
	\future Particles and antiparticles have different degeneracy
	factors, so we separately use the expansions one at a time. It
	is probably better to separately generate a new expansion
	function which automatically handles the sum of particles and
	antiparticles.
    */
    double pair_fun(double x, fermion &f, double T, bool log_mode);
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
