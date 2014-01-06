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
#ifndef O2SCL_SN_FERMION_H
#define O2SCL_SN_FERMION_H

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/root_cern.h>
#include <o2scl/inte.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>

#include <o2scl/part_deriv.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a relativistic fermion

      \note This class does not work with inc_rest_mass=true. For
      example, the integration limits in calc_mu() need to be reworked
      for this case.

      This implements an equation of state for a relativistic fermion
      using direct integration. After subtracting the rest mass from
      the chemical potentials, the distribution function is
      \f[
      \left\{1+\exp[(\sqrt{k^2+m^{* 2}}-m-\nu)/T]\right\}^{-1}
      \f]
      where \f$ k \f$ is the momentum, \f$ \nu \f$ is the effective
      chemical potential, \f$ m \f$ is the rest mass, and \f$ m^{*}
      \f$ is the effective mass.  For later use, we define \f$ E^{*} =
      \sqrt{k^2 + m^{*2}} \f$ . The degeneracy parameter is
      \f[
      \psi=(\nu+(m-m^{*}))/T 
      \f] 
      For \f$ \psi \f$ greater than \ref deg_limit (degenerate
      regime), a finite interval integrator is used and for \f$ \psi
      \f$ less than \ref deg_limit (non-degenerate regime), an
      integrator over the interval from \f$ [0,\infty) \f$ is
      used. The upper limit on the degenerate integration is given by
      the value of the momentum \f$ k \f$ which is the solution of
      \f[
      (\sqrt{k^2+m^{*,2}}-m-\nu)/T=\mathrm{f{l}imit}
      \f]
      which is
      \f[
      \sqrt{(m+{\cal L})^2-m^{*2}}
      \f]
      where \f$ {\cal L}\equiv\mathrm{f{l}imit}\times T+\nu \f$ .

      For the entropy integration, we set the lower limit
      to
      \f[
      2 \sqrt{\nu^2+2 \nu m} - \mathrm{upper~limit}
      \f]
      since the only contribution to the entropy is at the Fermi surface.
      \comment
      I'm not sure, but I think this is an expression determined
      from a small T taylor expansion of the argument of the 
      exponential.
      \endcomment

      In the non-degenerate regime, we make the substitution
      \f$ u=k/T \f$ to ensure that the variable of integration
      scales properly.

      Uncertainties are given in \ref unc.

      \todo This needs to be corrected to calculate \f$ \sqrt{k^2+m^{*
      2}}-m \f$ gracefully when \f$ m^{*}\approx m << k \f$ .
      \todo Call error handler if inc_rest_mass is true or update
      to properly treat the case when inc_rest_mass is true.
      
      \b Evaluation \b of \b the \b derivatives

      The relevant
      derivatives of the distribution function are
      \f[
      \frac{\partial f}{\partial T}=
      f(1-f)\frac{E^{*}-m-\nu}{T^2}
      \f]
      \f[
      \frac{\partial f}{\partial \nu}=
      f(1-f)\frac{1}{T}
      \f]
      \f[
      \frac{\partial f}{\partial k}=
      -f(1-f)\frac{k}{E^{*} T}
      \f]
      \f[
      \frac{\partial f}{\partial m^{*}}=
      -f(1-f)\frac{m^{*}}{E^{*} T}
      \f]

      We also need the derivative of the entropy integrand w.r.t. the 
      distribution function, which is
      \f[
      {\cal S}\equiv f \ln f +(1-f) \ln (1-f) \qquad
      \frac{\partial {\cal S}}{\partial f} = \ln 
      \left(\frac{f}{1-f}\right) = 
      \left(\frac{\nu-E^{*}+m}{T}\right)
      \f]
      where the entropy density is
      \f[
      s = - \frac{g}{2 \pi^2} \int_0^{\infty} {\cal S} k^2 d k
      \f]

      The derivatives can be integrated directly (\ref method = \ref
      direct) or they may be converted to integrals over the
      distribution function through an integration by parts (\ref
      method = \ref byparts)
      \f[
      \int_a^b f(k) \frac{d g(k)}{dk} dk = \left.f(k) g(k)\right|_{k=a}^{k=b}
      - \int_a^b g(k) \frac{d f(k)}{dk} dk 
      \f]
      using the distribution function for \f$ f(k) \f$ and 0 and 
      \f$ \infty \f$ as the limits, we have
      \f[
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{d g(k)}{dk} f dk =
      \frac{g}{2 \pi^2} \int_0^{\infty} g(k) f (1-f) \frac{k}{E^{*} T} dk 
      \f]
      as long as \f$ g(k) \f$ vanishes at \f$ k=0 \f$ .
      Rewriting,
      \f[
      \frac{g}{2 \pi^2} \int_0^{\infty} h(k) f (1-f) dk =
      \frac{g}{2 \pi^2} \int_0^{\infty} f \frac{T}{k} 
      \left[ h^{\prime} E^{*}-\frac{h E^{*}}{k}+\frac{h k}{E^{*}} \right] dk
      \f]
      as long as \f$ h(k)/k \f$ vanishes at \f$ k=0 \f$ .

      \b Explicit \b forms

      1) The derivative of the density wrt the chemical potential
      \f[
      \left(\frac{d n}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2}{T} f (1-f) dk
      \f]
      Using \f$ h(k)=k^2/T \f$ we get
      \f[
      \left(\frac{d n}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} 
      \left(\frac{k^2+E^{*2}}{E^{*}}\right) f dk
      \f]

      2) The derivative of the density wrt the temperature
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2(E^{*}-m-\nu)}{T^2} 
      f (1-f) dk
      \f]
      Using \f$ h(k)=k^2(E^{*}-\nu)/T^2 \f$ we get
      \f[
      \left(\frac{d n}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f}{T} 
      \left[2 k^2+E^{*2}-E^{*}\left(\nu+m\right)-
      k^2 \left(\frac{\nu+m}{E^{*}}\right)\right] dk
      \f]

      3) The derivative of the entropy wrt the chemical potential
      \f[
      \left(\frac{d s}{d \mu}\right)_T = 
      \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
      \frac{(E^{*}-m-\nu)}{T^2} dk
      \f]
      This verifies the Maxwell relation
      \f[
      \left(\frac{d s}{d \mu}\right)_T =
      \left(\frac{d n}{d T}\right)_{\mu}
      \f]

      4) The derivative of the entropy wrt the temperature
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
      \frac{(E^{*}-m-\nu)^2}{T^3} dk
      \f]
      Using \f$ h(k)=k^2 (E^{*}-\nu)^2/T^3 \f$ 
      \f[
      \left(\frac{d s}{d T}\right)_{\mu} = 
      \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f(E^{*}-m-\nu)}{E^{*}T^2} 
      \left[E^{* 3}+3 E^{*} k^2- (E^{* 2}+k^2)(\nu+m)\right] d k
      \f]

      5) The derivative of the density wrt the effective mass
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      -\frac{g}{2 \pi^2} \int_0^{\infty} 
      \frac{k^2 m^{*}}{E^{*} T} f (1-f) dk
      \f]
      Using \f$ h(k)=-(k^2 m^{*})/(E^{*} T) \f$ we get
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
      -\frac{g}{2 \pi^2} \int_0^{\infty} 
      m^{*} f dk
      \f]

      \note The dsdT integration may fail if the system is
      very degenerate. When method is byparts, the integral involves a
      large cancellation between the regions from \f$ k \in (0,
      \mathrm{ulimit/2}) \f$ and \f$ k \in (\mathrm{ulimit/2},
      \mathrm{ulimit}) \f$. Switching to method=direct and setting the
      lower limit to \f$ \mathrm{llimit} \f$, may help, but recent
      testing on this gave negative values for dsdT. For very
      degenerate systems, an expansion may be better than trying
      to perform the integration. The value of the integrand
      at k=0 also looks like it might be causing difficulties.

      \future It might be worth coding up direct differentiation, or
      differentiating the eff results, as these may succeed more
      generally.

      \future This class will have difficulty with extremely degenerate
      or extremely non-degnerate systems. Fix this.

      \future Create a more intelligent method for dealing with bad 
      initial guesses for the chemical potential in calc_density().
  */
  class sn_fermion : public fermion_deriv_thermo {
    
  public:

    /// Create a fermion with mass \c m and degeneracy \c g
    sn_fermion();

    virtual ~sn_fermion();
    
    /** \brief Limit of arguments of exponentials for Fermi functions 
	(default 200.0)
    */
    double exp_limit;

    /** \brief The critical degeneracy at which to switch integration 
	techniques (default 2.0)
    */
    double deg_limit;
    
    /** \brief The limit for the Fermi functions (default 20.0)
	
	sn_fermion will ignore corrections smaller than about
	\f$ \exp(-\mathrm{f{l}imit}) \f$ . 
    */
    double upper_limit_fac;
    
    /// Storage for the most recently calculated uncertainties 
    fermion_deriv unc;
    
    /** \name Method of computing derivatives
     */
    //@{
    /// Method (default is \ref byparts)
    int method;
    /// In the form containing \f$ f(1-f) \f$ .
    static const int direct=1;
    /// Integrate by parts
    static const int byparts=2;
    //@}

    /** \brief Calculate properties as function of chemical potential
    */
    virtual void calc_mu(fermion_deriv &f, double temper);

    /** \brief Calculate properties as function of density
     */
    virtual void calc_density(fermion_deriv &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	chemical potential
    */
    virtual void pair_mu(fermion_deriv &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	density
     */
    virtual void pair_density(fermion_deriv &f, double temper);

    /// Calculate effective chemical potential from density
    virtual void nu_from_n(fermion_deriv &f, double temper);

    /** \brief Set inte objects
	
	The first integrator is used for non-degenerate integration
	and should integrate from 0 to \f$ \infty \f$ (like \ref
	inte_qagiu_gsl). The second integrator is for the degenerate
	case, and should integrate between two finite values.
    */
    void set_inte(inte<funct> &unit, inte<funct> &udit);

    /** \brief Set the solver for use in calculating the chemical
	potential from the density */
    void set_density_root(root<funct> &rp) {
      density_root=&rp;
      return;
    }

    /// The default integrator for the non-degenerate regime
    inte_qagiu_gsl<funct> def_nit;

    /// The default integrator for the degenerate regime
    inte_qag_gsl<funct> def_dit;

    /// The default solver for npen_density() and pair_density()
    root_cern<funct> def_density_root;
    
    /// Return string denoting type ("sn_fermion")
    virtual const char *type() { return "sn_fermion"; };

    /// Calibrate with more accurate tabulated results
    double deriv_calibrate(fermion_deriv &f, int verbose, 
			   std::string fname="");

  protected:

#ifndef DOXYGEN_NO_O2NS

    /// Pointer to the data object
    fermion_deriv *fp;

    /// Temperature
    double T;

    /// The integrator for non-degenerate fermions
    inte<funct> *nit;

    /// The integrator for degenerate fermions
    inte<funct> *dit;

    /// The solver for calc_density() and pair_density()
    root<funct> *density_root;

    /** \name The integrands, as a function of \f$ u=k/T \f$, for 
	non-degenerate integrals
    */
    //@{
    double density_fun(double u);
    double energy_fun(double u);
    double entropy_fun(double u);
    double density_T_fun(double k);
    double density_mu_fun(double k);
    double entropy_T_fun(double k);
    double density_ms_fun(double k);
    //@}

    /** \name The integrands, as a function of momentum, for the 
	degenerate integrals
    */
    //@{
    double deg_density_fun(double u);
    double deg_energy_fun(double u);
    double deg_entropy_fun(double u);
    double deg_density_T_fun(double k);
    double deg_density_mu_fun(double k);
    double deg_entropy_T_fun(double k);
    double deg_density_ms_fun(double k);
    //@}

    /** \brief Solve for the chemical potential from the density 
	for calc_density()
    */
    double solve_fun(double x);

    /** \brief Solve for the chemical potential from the density 
	for pair_density()
    */
    double pair_fun(double x);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
