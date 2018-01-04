/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_FERMION_H
#define O2SCL_FERMION_H

/** \file fermion.h
    \brief File defining \ref o2scl::fermion
*/
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

// For gsl_sf_fermi_dirac_int()
#include <gsl/gsl_specfunc.h>

#include <o2scl/constants.h>
#include <o2scl/funct.h>
#include <o2scl/root.h>
#include <o2scl/root_cern.h>
#include <o2scl/part.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Fermion class

      This class adds two member data variables, \ref kf and \ref
      del, for the Fermi momentum and the gap, respectively.
  */
  class fermion : public part {

  public:

    /// Fermi momentum
    double kf;
    /// Gap
    double del;

    /// Create a fermion with mass \c mass and degeneracy \c dof.
    fermion(double mass=0, double dof=0);

    virtual ~fermion() {
    }

    /// Return string denoting type ("fermion")
    virtual const char *type() { return "fermion"; }
    
    /// Copy constructor
    fermion(const fermion &f) {
      g=f.g;
      m=f.m;
      ms=f.ms;
      n=f.n;
      ed=f.ed;
      pr=f.pr;
      mu=f.mu;
      en=f.en;
      nu=f.nu;
      inc_rest_mass=f.inc_rest_mass;
      non_interacting=f.non_interacting;
      kf=f.kf;
      del=f.del;
    }

    /// Copy construction with operator=()
    fermion &operator=(const fermion &f) {
      if (this!=&f) {
	g=f.g;
	m=f.m;
	ms=f.ms;
	n=f.n;
	ed=f.ed;
	pr=f.pr;
	mu=f.mu;
	en=f.en;
	nu=f.nu;
	inc_rest_mass=f.inc_rest_mass;
	non_interacting=f.non_interacting;
	kf=f.kf;
	del=f.del;
      }
      return *this;
    }

  };

  /** \brief Fermion properties at zero temperature

      This is a base class for the computation of fermionic statistics
      at zero temperature. The more general case of finite temperature
      is taken care of by \ref fermion_eval_thermo objects. The
      primary functions are calc_mu_zerot() and calc_density_zerot()
      which compute all the thermodynamic quantities as a function of
      the chemical potential, or the density, respectively.
      
      \future Use hypot() and other more accurate functions for the
      analytic expressions for the zero temperature integrals. [Progress
      has been made, but there are probably other functions which may
      break down for small but finite masses and temperatures]
  */
  class fermion_zerot {

  public:

    fermion_zerot() {
    };

    virtual ~fermion_zerot() {
    }

    /// \name Zero-temperature fermions 
    //@{
    /** \brief Calculate the Fermi momentum from the density
	
	Uses the relation \f$ k_F = ( 6 \pi^2 n /g )^{1/3} \f$
    */
    void kf_from_density(fermion &f); 

    /** \brief Energy density at T=0 from \ref fermion::kf and \ref part::ms

	Calculates the integral 
	\f[
	\varepsilon = \frac{g}{2 \pi^2} \int_0^{k_F} k^2 
	\sqrt{k^2+m^{* 2}} d k
	\f]
    */
    void energy_density_zerot(fermion &f); 

    /** \brief Pressure at T=0 from \ref fermion::kf and \ref part::ms

	Calculates the integral 
	\f[
	P=\frac{g}{6 \pi^2} \int_0^{k_F} \frac{k^4}{\sqrt{k^2+m^{* 2}}} d k
	\f]
    */
    void pressure_zerot(fermion &f); 

    /** \brief Zero temperature fermions from \ref part::mu or 
	\ref part::nu and \ref part::ms
    */
    virtual void calc_mu_zerot(fermion &f);
    
    /** \brief Zero temperature fermions from \ref part::n and \ref part::ms
     */
    virtual void calc_density_zerot(fermion &f);
    //@}

  };
  
  /** \brief Fermion with finite-temperature thermodynamics
      [abstract base]

      This is an abstract base for the computation of
      finite-temperature fermionic statistics. Different children
      (e.g. \ref fermion_eff and \ref fermion_rel) use different
      techniques to computing the momentum integrations.

      Because massless fermions at finite temperature are much
      simpler, there are separate member functions included in this
      class to handle them. The functions massless_calc_density() and
      massless_calc_mu() compute the thermodynamics of massless
      fermions at finite temperature given the density or the chemical
      potentials. The functions massless_pair_density() and
      massless_pair_mu() perform the same task, but automatically
      include antiparticles.

      The function massless_calc_density() uses a \ref root object to
      solve for the chemical potential as a function of the density.
      The default is an object of type root_cern. The function
      massless_pair_density() does not need to use the \ref root
      object because of the simplification afforded by the inclusion
      of antiparticles.

      \future Create a Chebyshev approximation for inverting the 
      the Fermi functions for massless_calc_density() functions?
  */
  class fermion_eval_thermo : public fermion_zerot {

  public:
    
    fermion_eval_thermo();

    virtual ~fermion_eval_thermo() {
    }

    /** \brief Non-degenerate expansion for fermions

	Attempts to evaulate thermodynamics of a non-degenerate
	fermion. If the result is accurate to within the requested
	precision, this function returns <tt>true</tt>, and otherwise
	this function returns <tt>false</tt> and the values in stored
	in the <tt>pr</tt>, <tt>n</tt>, <tt>en</tt>, and <tt>ed</tt>
	field are meaningless.

	If \f$ \mu \f$ is negative and sufficiently far from zero,
	then the thermodynamic quantities are smaller than the smallest
	representable double-precision number. In this case,
	this function will return <tt>true</tt> and report all
	quantities as zero.
      
	Defining \f$ \psi \equiv (\mu-m)/T \f$, \f$ t \equiv T/m \f$,
	and \f$ d \equiv g~m^4/(2 \pi^2) \f$ the pressure 
	in the non-degenerate limit (\f$ \psi \rightarrow - \infty \f$)
	is (\ref Johns96)
	\f[
	P = d \sum_{n=1}^{\infty} P_n
	\f]
	where 
	\f[
	P_n \equiv \left(-1\right)^{n+1} \left(\frac{t^2}{n^2}\right)
	e^{n \left(\psi+1/t\right)} K_2 \left( \frac{n}{t} \right)
	\f]
	The density is then
	\f[
	n = d \sum_{n=1}^{\infty} \frac{n P_n}{T}
	\f]
	and the entropy density is
	\f[
	s = \frac{d}{m} \sum_{n=1}^{\infty} \left\{ \frac{2 P_n}{t}
	-\frac{n P_n}{t^2}+ 
	\frac{\left(-1\right)^{n+1}}{2 n} 
	e^{n \left(\psi+1/t\right)} \left[ K_1 \left( \frac{n}{t} 
	\right)+K_3 \left( \frac{n}{t} \right) \right]
	\right\}
	\f]

	This function is accurate over a wide range of conditions
	when \f$ \psi < -4 \f$.

	The ratio of the nth term to the first term in the pressure
	series is
	\f[
	R_n \equiv \frac{P_{n}}{P_{1}} = \frac{(-1)^{n+1} 
	e^{(n-1)(\psi+1/t)} K_2(n/t) }{n^2 K_2(1/t)}
	\f]
	This function currently uses 20 terms in the series and
	immediately returns <tt>false</tt> if \f$ |R_{20}| \f$
	is greater than <tt>prec</tt>

	In the nondegenerate and nonrelativistic (\f$ t \rightarrow 0
	\f$) limit, the argument to the Bessel functions and the
	exponential becomes too large. In this case, it's better to
	use the expansions, e.g. for \f$ x \equiv n/t \rightarrow
	\infty \f$,
	\f[
	\sqrt{\frac{2 x}{\pi}} e^{x} K_2(x) \approx
	1 + \frac{3}{8 x} - \frac{15}{128 x^2} + ...
	\f]
	The current code currently goes up to \f$ x^{-12} \f$ in the
	expansion, which is enough for the default precision of \f$
	10^{-18} \f$ since \f$ (20/700)^{12} \sim 10^{-19} \f$.
    */
    virtual bool calc_mu_ndeg(fermion &f, double temper, 
			      double prec=1.0e-18, bool inc_antip=false);

    /** \brief Degenerate expansion for fermions

	Attempts to evaulate thermodynamics of a degenerate fermion.
	If the result is accurate to within the requested precision,
	this function returns <tt>true</tt>, and otherwise this
	function returns <tt>false</tt> and the values in stored in
	the <tt>pr</tt>, <tt>n</tt>, <tt>en</tt>, and <tt>ed</tt>
	field are meaningless.

	The pressure, density, and energy density, should be accurate
	to the requested precision, but the first term in the series
	expansion for the entropy is zero, so the entropy is one order
	lower in accuracy.

	\future Make a function like this for dndm, dsdT, etc. 
	for fermion_deriv .
    */
    virtual bool calc_mu_deg(fermion &f, double temper, 
			     double prec=1.0e-18);

    /** \brief Calculate properties as function of chemical potential
     */
    virtual void calc_mu(fermion &f, double temper)=0;

    /** \brief Calculate properties as function of density

	\note This function returns an integer value, in contrast to
	\ref calc_mu(), because of the potential for non-convergence.
     */
    virtual int calc_density(fermion &f, double temper)=0;

    /** \brief Calculate properties with antiparticles as function of
	chemical potential
    */
    virtual void pair_mu(fermion &f, double temper)=0;

    /** \brief Calculate properties with antiparticles as function of
	density

	\note This function returns an integer value, in contrast to
	\ref pair_mu(), because of the potential for non-convergence.
    */
    virtual int pair_density(fermion &f, double temper)=0;

    /// \name Massless fermions
    //@{
    /// Finite temperature massless fermions
    virtual void massless_calc_mu(fermion &f, double temper);

    /// Finite temperature massless fermions
    virtual void massless_calc_density(fermion &f, double temper);

    /** \brief Finite temperature massless fermions and antifermions 
     */
    virtual void massless_pair_mu(fermion &f, double temper);

    /** \brief Finite temperature massless fermions and antifermions 

	In the cases \f$ n^3 \gg T \f$ and \f$ T \gg n^3 \f$ ,
	expansions are used instead of the exact formulas to avoid
	loss of precision.

	In particular, using the parameter
	\f[
	\alpha = \frac{g^2 \pi^2 T^6}{243 n^2}
	\f]
	and defining the expression
	\f[
	\mathrm{cbt} = \alpha^{-1/6} \left( -1 + \sqrt{1+\alpha}\right)^{1/3}
	\f]
	we can write the chemical potential as
	\f[
	\mu = \frac{\pi T}{\sqrt{3}} \left(\frac{1}{\mathrm{cbt}} -
	\mathrm{cbt} \right)
	\f]
	
	These expressions, however, do not work well when \f$ \alpha
	\f$ is very large or very small, so series expansions are
	used whenever \f$ \alpha > 10^{4} \f$ or 
	\f$ \alpha < 3 \times 10^{-4} \f$. For small \f$ \alpha \f$, 
	\f[
	\left(\frac{1}{\mathrm{cbt}} -
	\mathrm{cbt} \right) \approx
	\frac{2^{1/3}}{\alpha^{1/6}} - 
	\frac{\alpha^{1/6}}{2^{1/3}} +
	\frac{\alpha^{5/6}}{6{\cdot}2^{2/3}} +
	\frac{\alpha^{7/6}}{12{\cdot}2^{1/3}} -
	\frac{\alpha^{11/6}}{18{\cdot}2^{2/3}} -
	\frac{5 \alpha^{13/6}}{144{\cdot}2^{1/3}} +
	\frac{77 \alpha^{17/6}}{2592{\cdot}2^{2/3}}
	\f]
	and for large \f$ \alpha \f$, 
	\f[
	\left(\frac{1}{\mathrm{cbt}} -
	\mathrm{cbt} \right) \approx
	\frac{2}{3} \sqrt{\frac{1}{\alpha}} - 
	\frac{8}{81} \left(\frac{1}{\alpha}\right)^{3/2} +
	\frac{32}{729} \left(\frac{1}{\alpha}\right)^{5/2}
	\f]

	This approach works to within about 1 \part in \f$ 10^{12} \f$,
	and is tested in <tt>fermion_ts.cpp</tt>.
	
	\future This could be improved by including more terms
	in the expansions.
    */
    virtual void massless_pair_density(fermion &f, double temper);
    //@}
    
    /** \brief Set the solver for use in massless_calc_density() */ 
    void set_massless_root(root<> &rp) {
      massless_root=&rp;
      return;
    }

    /** \brief The default solver for massless_calc_density()
	
	We default to a solver of type root_cern here since we
	don't have a bracket or a derivative.
    */
    root_cern<> def_massless_root;

    /// Return string denoting type ("fermion_eval_thermo")
    virtual const char *type() { return "fermion_eval_thermo"; }

    /** \brief Test the thermodynamics of calc_density() and 
	calc_mu()

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

	\future Also calibrate massless fermions?
	\future Convert into separate class?
    */
    virtual double calibrate(fermion &f, int verbose=0, bool test_pair=false,
			     std::string fname="");

#ifndef DOXYGEN_NO_O2NS

  protected:
    
    /// A pointer to the solver for massless fermions
    root<> *massless_root;

    /// Solve for the chemical potential for massless fermions
    double massless_solve_fun(double x, fermion &f, double temper);

#endif
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
