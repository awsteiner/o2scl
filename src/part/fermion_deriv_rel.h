/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#ifndef O2SCL_FERMION_DERIV_REL_H
#define O2SCL_FERMION_DERIV_REL_H

/** \file fermion_deriv_rel.h
    \brief File defining \ref o2scl::fermion_deriv_rel
*/

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
#include <o2scl/fermion_rel.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a relativistic fermion

      \note This class only has preliminary support for
      inc_rest_mass=true (more testing should be done, particularly
      for the "pair" functions)

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

      In the non-degenerate regime, we make the substitution \f$ u=k/T
      \f$ to help ensure that the variable of integration scales
      properly.

      Uncertainties are given in \ref unc.

      \future This class may need more corrections to ensure
      quantities like \f$ \sqrt{k^2+m^{*2}}-m \f$ are computed
      accurately when \f$ m^{*}\approx m \ll k \f$ .
      
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
      method = \ref by_parts)
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
      \comment
      This derivative may be written in terms of the 
      others
      \f[
      \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = \frac{3 n}{m^{*}}
      - \frac{T}{m^{*}}\left[ \left(\frac{d n}{d T}\right)_{\mu}
      +\frac{\mu}{T} \left(\frac{d n}{d \mu}\right)_{T}
      \right] - \left(\frac{d n}{d \mu}\right)_{T}
      \f]
      \endcomment

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

      \future The option err_nonconv=false is not really implemented
      yet.

      \future It might be worth coding up direct differentiation, or
      differentiating the eff results, as these may succeed more
      generally.

      \future This class will have difficulty with extremely degenerate
      or extremely non-degnerate systems. Fix this by using the
      expansions similar to the method used in \ref o2scl::fermion_rel .

      \future Create a more intelligent method for dealing with bad 
      initial guesses for the chemical potential in calc_density().
  */
  class fermion_deriv_rel : public fermion_deriv_thermo {
    
  public:

    /// Create a fermion with mass \c m and degeneracy \c g
    fermion_deriv_rel();

    virtual ~fermion_deriv_rel();
    
    /** \brief Limit of arguments of exponentials for Fermi functions 
	(default 200.0)
    */
    double exp_limit;

    /** \brief The critical degeneracy at which to switch integration 
	techniques (default 2.0)
    */
    double deg_limit;
    
    /** \brief The limit for the Fermi functions (default 20.0)
	
	fermion_deriv_rel will ignore corrections smaller than about
	\f$ \exp(-\mathrm{f{l}imit}) \f$ . 
    */
    double upper_limit_fac;
    
    /// Storage for the most recently calculated uncertainties 
    fermion_deriv unc;

    /// Object for computing non-derivative quantities
    fermion_rel fr;

    /** \name Method of computing derivatives
     */
    //@{
    /// Method (default is \ref automatic)
    int method;
    /// Automatically choose method
    static const int automatic=0;
    /// In the form containing \f$ f(1-f) \f$ .
    static const int direct=1;
    /// Integrate by parts
    static const int by_parts=2;
    //@}

    /** \brief If true, call the error handler when convergence 
	fails (default true)
    */
    bool err_nonconv;

    /** \brief Calculate properties as function of chemical potential
    */
    virtual int calc_mu(fermion_deriv &f, double temper);

    /** \brief Calculate properties as function of density
     */
    virtual int calc_density(fermion_deriv &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	chemical potential
    */
    virtual int pair_mu(fermion_deriv &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	density
     */
    virtual int pair_density(fermion_deriv &f, double temper);

    /// Calculate effective chemical potential from density
    virtual int nu_from_n(fermion_deriv &f, double temper);

    /** \brief Set inte objects
	
	The first integrator is used for non-degenerate integration
	and should integrate from 0 to \f$ \infty \f$ (like \ref
	o2scl::inte_qagiu_gsl). The second integrator is for the
	degenerate case, and should integrate between two finite
	values.
    */
    void set_inte(inte<> &unit, inte<> &udit);
    
    /** \brief Set the solver for use in calculating the chemical
	potential from the density */
    void set_density_root(root<> &rp) {
      density_root=&rp;
      return;
    }

    /// The default integrator for the non-degenerate regime
    inte_qagiu_gsl<> def_nit;

    /// The default integrator for the degenerate regime
    inte_qag_gsl<> def_dit;

    /// The default solver for npen_density() and pair_density()
    root_cern<> def_density_root;
    
    /// Return string denoting type ("fermion_deriv_rel")
    virtual const char *type() { return "fermion_deriv_rel"; };

    /** \brief Calibrate with more accurate tabulated results

	This compares the approximation to the exact results over a
	grid with \f$ T = \left\{10^{-2},1,10^{2}\right\} \f$, \f$
	\log_{10} (m/T) = \left\{ -3,-2,-1,0,1,2,3\right\} \f$, and
	\f$ \log_{10} \psi = \left\{ -3,-2,-1,0,1\right\} \f$, where
	\f$ \psi \equiv \left(\mu-m\right)/T \f$ using
	calc_density() and calc_mu(), with both inc_rest_mass
	taking both values <tt>true</tt> and <tt>false</tt>.
	
	The <tt>verbose</tt> parameter controls the amount of
	output, and \c fname is the filename for the file
	<tt>fermion_cal.o2</tt>.
     */
    double deriv_calibrate(fermion_deriv &f, int verbose, 
			   std::string fname="");

  protected:

#ifndef DOXYGEN_NO_O2NS

    /// The internal integration method
    int intl_method;

    /// The integrator for non-degenerate fermions
    inte<> *nit;

    /// The integrator for degenerate fermions
    inte<> *dit;

    /// The solver for calc_density() and pair_density()
    root<> *density_root;

    /** \name The integrands, as a function of \f$ u=k/T \f$, for 
	non-degenerate integrals
    */
    //@{
    double density_fun(double u, fermion_deriv &f, double T);
    double energy_fun(double u, fermion_deriv &f, double T);
    double entropy_fun(double u, fermion_deriv &f, double T);
    double density_T_fun(double k, fermion_deriv &f, double T);
    double density_mu_fun(double k, fermion_deriv &f, double T);
    double entropy_T_fun(double k, fermion_deriv &f, double T);
    double density_ms_fun(double k, fermion_deriv &f, double T);
    //@}

    /** \name The integrands, as a function of momentum, for the 
	degenerate integrals
    */
    //@{
    double deg_density_fun(double u, fermion_deriv &f, double T);
    double deg_energy_fun(double u, fermion_deriv &f, double T);
    double deg_entropy_fun(double u, fermion_deriv &f, double T);
    double deg_density_T_fun(double k, fermion_deriv &f, double T);
    double deg_density_mu_fun(double k, fermion_deriv &f, double T);
    double deg_entropy_T_fun(double k, fermion_deriv &f, double T);
    double deg_density_ms_fun(double k, fermion_deriv &f, double T);
    //@}

    /** \brief Solve for the chemical potential from the density 
	for calc_density()
    */
    double solve_fun(double x, fermion_deriv &f, double T);

    /** \brief Solve for the chemical potential from the density 
	for pair_density()
    */
    double pair_fun(double x, fermion_deriv &f, double T);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
