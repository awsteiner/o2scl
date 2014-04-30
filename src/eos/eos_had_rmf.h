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

#ifndef O2SCL_RMF_EOS_H
#define O2SCL_RMF_EOS_H

#include <string>
#include <cmath>
#include <o2scl/lib_settings.h>
#include <o2scl/constants.h>
#include <o2scl/mm_funct.h>

#include <o2scl/part.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/fermion.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Relativistic mean field theory EOS

      This class computes the properties of nucleonic matter using a
      mean-field approximation to a field-theoretical model.
      
      Before sending neutrons and protons to these member functions,
      the masses should be set to and the degeneracy factor should be
      set to 2. Some models which can be loaded using
      <tt>o2scl_hdf::rmf_load()</tt> expect that the neutron and
      proton masses are set to the value stored in \ref mnuc.
    
      \note Since this EOS uses the effective masses and chemical
      potentials in the \ref o2scl::part class, the values of
      <tt>o2scl::part::non_interacting</tt> for neutrons and protons
      are set to false in many of the functions.

      \note Matter at two different densities can have the same
      chemical potentials, so the behavior of the function \ref
      o2scl::eos_had_rmf::calc_temp_p() is ambiguous. This arises because
      the field equations have more than one solution for a specified
      chemical potential. Internally, \ref
      o2scl::eos_had_rmf::calc_temp_p() either uses the initial guess
      specified by a call to \ref o2scl::eos_had_rmf::set_fields(), or
      uses hard-coded initial guess values typical for saturation
      densities. In order to ensure that the user gets the desired
      solution to the field equations, it may be necessary to specify
      a sufficiently accurate initial guess. There is no ambiguity in
      the behavior of \ref o2scl::eos_had_rmf::calc_eq_temp_p(), however.

      \note This class can fail to solve the meson field equations or
      fail to solve for the nucleon densities. By default the error
      handler is called when this happens. If \ref err_nonconv is
      false, then functions which don't converge (which also return
      <tt>int</tt>) will return a non-zero value. Note that the
      solvers (in \ref def_sat_mroot and \ref
      o2scl::eos_had_base::def_mroot) also has its own data member
      indicating how to handle nonconvergence \ref
      o2scl::mroot::err_nonconv which is separate.

      \comment
      AWS, 11/17/13: It is not clear that this is entirely necessary
      as almost all the CONV_ERR calls in eos_had_rmf.cpp are due to calls
      to solvers. It could be that then err_nonconv can be removed and
      all the eos_had_rmf functions just always directly return any
      nonzero values they get from solvers. One nice thing about the
      explicit CONV_ERR calls in eos_had_rmf.cpp is that it makes the code
      easier to read. In any case err_nonconv should probably be
      pushed up to eos_had_base.
      \endcomment

      \hline
      \b Background

      The full Lagragian can be written as a sum of several terms
      \f[
      {\cal L} = {\cal L}_{\mathrm{Dirac}} + {\cal L}_{\sigma} + 
      {\cal L}_{\omega} + {\cal L}_{\rho} + {\cal L}_{\mathrm{int}} \, .
      \f]

      The part for the nucleon fields is
      \f[
      {\cal L}_{\mathrm{Dirac}} = 
      \bar{\Psi}_i \left[ i {{\partial}\!\!\!{\slash}} - 
      g_{\omega} {{\omega}\!\!\!{\slash}} - \frac{g_{\rho}}{2} 
      {{\vec{\rho}}\!\!\!{\slash}}
      \vec{\tau} - M_i + g_{\sigma} \sigma - 
      e q_i A\!\!\!{\slash} \right] \Psi_i 
      \f]
      where \f$ \Psi \f$ is the nucleon field and \f$ \sigma, \omega
      \f$ and \f$ \rho \f$ are the meson fields. The meson masses
      are \f$ m_{\sigma}, m_{\omega} \f$ and \f$ m_{\rho}
      \f$ and meson-nucleon
      couplings are \f$ g_{\sigma}, g_{\omega} \f$ and \f$ g_{\rho}
      \f$ . The couplings \c cs, \c cw, and \c cr are related to \f$
      g_{\sigma}, g_{\omega} \f$ and \f$ g_{\rho} \f$ by
      \f[
      c_{\sigma} = g_{\sigma}/m_{\sigma} \quad
      c_{\omega} = g_{\omega}/m_{\omega} \quad \mathrm{and} \quad
      c_{\rho} = g_{\rho}/m_{\rho}
      \f]
      The nucleon masses are in \f$ M_i \f$ and stored in
      <tt>part::m</tt> and \f$ q_i \f$ just represents the charge (1
      for protons and 0 for neutrons). The Coulomb field, \f$ A_{\mu}
      \f$, is ignored in this class, but used in \ref
      o2scl::nucleus_rmf.

      The part for the \f$ \sigma \f$ field is
      \f[
      {\cal L}_{\sigma} =
      {\textstyle \frac{1}{2}} \left( \partial_{\mu} \sigma \right)^2 
      - {\textstyle \frac{1}{2}} m^2_{\sigma} \sigma^2 
      - \frac{b M}{3} \left( g_{\sigma} \sigma\right)^3 
      - \frac{c}{4} \left( g_{\sigma} \sigma\right)^4 \, .
      \f]
      where \f$ m_{\sigma} \f$ is the meson mass, 
      \f$ b \f$ and \f$ c \f$ are unitless couplings and
      \f$ M \f$ is a dimensionful scale, ususally taken to be
      939 MeV (which need not be equal to \f$ M_i \f$ above). 
      The coefficients \f$ b \f$ and \f$ c \f$ are related to the somewhat
      standard \f$ \kappa \f$ and \f$ \lambda \f$ by:
      \f[
      \kappa=2 M b \quad \lambda=6 c;
      \f]

      The part for the \f$ \omega \f$ field is
      \f[
      {\cal L}_{\omega} =
      - {\textstyle \frac{1}{4}} f_{\mu \nu} f^{\mu \nu} 
      + {\textstyle \frac{1}{2}} m^2_{\omega}\omega^{\mu}\omega_{\mu} 
      + \frac{\zeta}{24} g_{\omega}^4 \left(\omega^\mu \omega_\mu\right)^2
      \f]
      where \f$ m_{\omega} \f$ is the meson mass.
      
      The part for the \f$ \rho \f$ field is
      \f[
      {\cal L}_{\rho} = 
      - {\textstyle \frac{1}{4}} \vec{B}_{\mu \nu} \cdot \vec{B}^{\mu \nu}
      + {\textstyle \frac{1}{2}} m^2_{\rho} \vec{\rho}^{~\mu} \cdot 
      \vec{\rho}_{~\mu} 
      + \frac{\xi}{24} g_{\rho}^4 \left(\vec{\rho}^{~\mu}\right) \cdot 
      \vec{\rho}_{~\mu} 
      \f]

      Finally, additional meson interactions are
      \f[
      {\cal L}_{\mathrm{int}} = 
      g_{\rho}^2 f (\sigma, \omega) \vec{\rho}^{~\mu} \cdot 
      \vec{\rho}_{~\mu} \nonumber \\
      \f]
      The function \f$ f \f$ is the coefficient of \f$ g_r^2 \rho^2 \f$ 
      \f$ f(\sigma,\omega) = b_1 \omega^2 + b_2 \omega^4 + b_3 \omega^6 +
      a_1 \sigma + a_2 \sigma^2 + a_3 \sigma^3 + a_4 \sigma^4 +
      a_5 \sigma^5 + a_6 \sigma^6 \f$ 
      where the notation from \ref Horowitz01 is:
      \f$ f(\sigma,\omega) = \lambda_4 g_s^2 \sigma^2 + 
      \lambda_v g_w^2 \omega^2 \f$ 
      This implies \f$ b_1=\lambda_v g_w^2 \f$ and 
      \f$ a_2=\lambda_4 g_s^2 \f$ 

      The couplings, \c cs, \c cw, and \c cr all have units of \f$
      \mathrm{fm} \f$, and the couplings \c b, \c c, \c zeta and \c xi are
      unitless. The additional couplings from \ref Steiner05b, \f$ a_i
      \f$ have units of \f$ \mathrm{fm}^{(i-2)} \f$ and the couplings
      \f$ b_j \f$ have units of \f$ \mathrm{fm}^{(2j-2)} \f$ .

      When the variable \ref zm_mode is true, the effective mass is
      fixed using the approach of \ref Zimanyi90 .
    
      The expressions for the energy densities are often simplified in
      the literature using the field equations. These expressions are
      not used in this code since they are only applicable in infinite
      matter where the field equations hold, and are not suitable for
      use in applications (such as to finite nuclei in \ref
      o2scl::nucleus_rmf) where the spatial derivatives of the fields
      are non-zero. Notice that in the proper expressions for the
      energy density the similarity between terms in the pressure up
      to a sign. This procedure allows one to verify the thermodynamic
      identity even if the field equations are not solved and allows
      the user to add gradient terms to the energy density and
      pressure.

      See also \ref Muller96 .

      \hline
      \b Field \b equations

      The field equations are:
      \f[
      0 = m_{\sigma}^2 \sigma - g_{\sigma} \left( n_{s n} + n_{s p} \right)
      + b M g_{\sigma}^3 \sigma^2 + c g_{\sigma}^4 \sigma^3 -
      g_{\rho}^2 \rho^2 \frac{\partial f}{\partial \sigma}
      \f]
      \f[
      0 = m_{\omega}^2 \omega - g_{\omega} \left(n_n+n_p\right)
      + \frac{\zeta}{6} g_{\omega}^4 \omega^3 + g_{\rho}^2 \rho^2 
      \frac{\partial f}{\partial \omega}
      \f]
      \f[
      0 = m_{\rho}^2 \rho + \frac{1}{2} g_{\rho} \left(n_n-n_p\right)
      + 2 g_{\rho}^2 \rho f + \frac{\xi}{6} g_{\rho}^4 \rho^3
      \f]

      \hline
      \b Saturation \b properties

      Defining
      \f[
      U(\sigma)=\frac{1}{2} m_\sigma^2\sigma^2+\frac{b M}{3}(g_\sigma\sigma)^3
      +\frac{c}{4}(g_\sigma\sigma)^4\;, 
      \f]
      the binding energy per particle in symmetric matter at equilibrium
      is given by
      \f[
      \frac{E}{A} = \frac{1}{n_0} \left[U(\sigma_0)+
      \frac{1}{2} m_\omega\omega_0^2+
      \frac{\zeta}{8}(g_\omega\omega_0)^4+\frac{2}{\pi^2}
      \int\limits_0^{k_F} dk k^2\sqrt{k^2+M^{*2}} \right] 
      \f]
      where the Dirac
      effective mass is  \f$ M^{*}_i = M_i - g_{\sigma}\sigma_0 \f$ .
      The compressibility is given by
      \f[
      K=9\frac{g_\omega^2}{m_\omega^2}n_0+3\frac{k_F^2}{E_F^*}
      -9n_0\frac{M^{*2}}{E_F^{*2}}\left[\left(\frac{1}{g_\sigma^2}
      \frac{\partial^2}{\partial\sigma_0^2}+\frac{3}{g_\sigma M^*}
      \frac{\partial}{\partial\sigma_0}\right)
      U(\sigma_0)-3\frac{n_0}{E_F^*}\right]^{-1}\;.
      \f]
      The symmetry energy of bulk matter is given by
      \f[
      E_{sym} = \frac{k_F^2}{6 E_F^{*}} + \frac{ n }
      {8 \left(g_{\rho}^2/m_{\rho}^2 + 2 f (\sigma_0, \omega_0) 
      \right)} \, .
      \f]
    
      In the above equations, the subscipt \f$ 0 \f$ denotes the mean
      field values of \f$ \sigma \f$ and \f$ \omega \f$ .  For the case
      \f$ f=0 \f$ , the symmetry energy varies linearly with the density at
      large densities. The function \f$ f \f$ permits variations in the
      density dependence of the symmetry energy above nuclear matter
      density.
    
      \hline

      \todo 
      - The functions fcomp_fields(), fkprime_fields(), and fesym_fields()
      are not quite correct if the neutron and proton masses are different.
      For this reason, they are currently unused by saturation().
      - The fix_saturation() and calc_cr() functions use mnuc, and should
      be modified to allow different neutron and proton masses.
      - Check the formulas in the "Background" section
      - Make sure that this class properly handles particles for which 
      inc_rest_mass is true/false
      - The error handler is called sometimes when calc_e() is used
      to compute pure neutron matter. This should be fixed.

      \future
      - Finish putting the err_nonconv system into calc_p(),
      calc_temp_e() and fix_saturation(), etc.
      - It might be nice to remove explicit reference to the meson
      masses in functions which only compute nuclear matter since they
      are unnecessary. This might, however, demand redefining some of
      the couplings.
      - Fix calc_p() to be better at guessing
      - The number of couplings is getting large, maybe new
      organization is required.
      - Overload eos_had_base::fcomp() with an exact version
      - It would be nice to analytically compute the Jacobian
      of the field equations for the solver

  */
  class eos_had_rmf : public eos_had_base_temp_pres {

  public:

    /// \name Other data members
    //@{
    /** \brief The number of separate calls to the solver 
	that the <tt>calc_e</tt> functions take (default 20)

	Values larger than about \f$ 10^4 \f$ are probably
	not useful. 
    */
    size_t calc_e_steps;

    /** \brief If true, solve for relative densities rather than 
	absolute densities (default false)

	Setting this to true makes \ref calc_temp_e() and \ref
	calc_e() more accurate at low densities. 
    */
    bool calc_e_relative;

    /// Modifies method of calculating effective masses (default false)
    bool zm_mode;
    
    /** \brief Verbosity parameter

	If this is greater than zero, then some functions report
	on their progress.
	- The function \ref saturation() reports progress towards
	computing the properties of nuclear matter near saturation.
	- The functions \ref calc_e() and \ref calc_temp_e() report
	progress on solving for matter at a fixed density.
     */
    int verbose;

    /** \brief If true, throw exceptions when the function calc_e()
	does not converge (default true)
    */
    bool err_nonconv;
    //@}
    
    /// \name Masses
    //@{
    /** \brief The scale \f$ M \f$
	
	This need not be exactly equal to the neutron or proton mass, 
	but provides the scale for the coupling \c b.
    */
    double mnuc;

    /// \f$ \sigma \f$ mass (in \f$ \mathrm{fm}^{-1} \f$ )
    double ms;

    /// \f$ \omega \f$ mass (in \f$ \mathrm{fm}^{-1} \f$ )
    double mw;

    /// \f$ \rho \f$ mass (in \f$ \mathrm{fm}^{-1} \f$ )
    double mr;

    //@}

    /// \name Standard couplings (including nonlinear sigma terms)
    //@{
    double cs, cw, cr, b, c;
    //@}

    /// \name Quartic terms for omega and rho.
    //@{
    double zeta, xi;
    //@}
  
    /// \name Additional isovector couplings
    //@{
    double a1, a2, a3, a4, a5, a6, b1, b2, b3;
    //@}

    eos_had_rmf();

    /* \brief Load parameters for model named 'model'
	
	Presently accepted values from file rmfdata/model_list:
	\include rmfdata/model_list
	
	In these files, the nucleon and meson masses are by default
	specified in MeV, and cs, cw, and cr are given in fm. The
	parameters b and c are both unitless. If the bool 'oakstyle' is
	true, then load() assumes that gs, gw, and gr have been given
	where gs and gw are as usual, but gr is a factor of two smaller
	than usual, and g2 and g3 have been given where g2 = -b M gs^3
	and g3 = c gs^4. If tokistyle is true, then it is additionally
	assumed that c3 is given where c3=zeta/6*gw^4.
	
	If \c external is true, then model is the filename (relative
	to the current directory) of the file containing the model
	parameters. Otherwise, the model is assumed to be present in
	the \o2 library data directory.
    */
    //int load(std::string model, bool external=false);

    /// \name Compute EOS 
    //@{
    /** \brief Equation of state as a function of density

	Initial guesses for the chemical potentials are taken
	from the user-given values. Initial guesses for the fields
	can be set by set_fields(), or default values will be used.
	After the call to calc_e(), the final values of the fields
	can be accessed through get_fields(). 

	This is a little more robust than the standard version
	in the parent \ref eos_had_base.
	
	\future Improve the operation of this function when the
	proton density is zero.
       
    */
    virtual int calc_e(fermion &ne, fermion &pr, thermo &lth);
  
    /**  \brief Equation of state as a function of chemical potential
	 
	 Solves for the field equations automatically.
	 
	 \future It may be possible to make the solver for the
	 field equations more robust
    */
    virtual int calc_p(fermion &ne, fermion &pr, thermo &lth);

    /** \brief Equation of state and meson field equations 
	as a function of chemical potentials
      
	This calculates the pressure and energy density as a function of
	\f$ \mu_n,\mu_p,\sigma,\omega,rho \f$ . When the field equations
	have been solved, \c f1, \c f2, and \c f3 are all zero. 
	
	The thermodynamic identity is satisfied even when the field
	equations are not solved.

	\future Probably best to have f1, f2, and f3 scaled
	in some sensible way, i.e. scaled to the fields?
    */
    virtual int calc_eq_p(fermion &neu, fermion &p, double sig, 
			  double ome, double rho, double &f1, 
			  double &f2, double &f3, thermo &th);

    /** \brief Equation of state and meson field equations as a 
	function of chemical potentials at finite temperature

	Analogous to \ref calc_eq_p() except at finite temperature.
    */
    virtual int calc_eq_temp_p(fermion &ne, fermion &pr, double temper, 
			       double sig, double ome, double rho, double &f1, 
			       double &f2, double &f3, thermo &th);
    
    /** \brief Equation of state as a function of chemical potential
	
	Solves for the field equations automatically.
    */
    virtual int calc_temp_p(fermion &ne, fermion &pr, double T,
			    thermo &lth);

    /** \brief Equation of state as a function of densities at 
	finite temperature
    */
    int calc_temp_e(fermion &ne, fermion &pr, double T, 
		    thermo &lth);
    //@}

    /// \name Saturation properties
    //@{
    /** \brief Calculate cs, cw, cr, b, and c from the saturation 
	properties

	Note that the meson masses and \ref mnuc must be specified
	before calling this function.

	This function does not give correct results when bool zm_mode 
	is true. 
	
	\c guess_cs, \c guess_cw, \c guess_b, and \c guess_c are
	initial guesses for \c cs, \c cw, \c b, and \c c respectively.

	\todo 
	- Fix this for zm_mode=true
	- Ensure solver is more robust
	
    */
    int fix_saturation(double guess_cs=4.0, double guess_cw=3.0, 
		       double guess_b=0.001, double guess_c=-0.001);
    
    /** \brief Calculate properties of nuclear matter at the
	saturation density

	This function first constructs an initial guess, increasing
	the chemical potentials if required to ensure the neutron and
	proton densities are finite, and then uses \ref
	eos_had_rmf::sat_mroot to solve the field equations and ensure
	that the neutron and proton densities are equal and the
	pressure is zero. The quantities \ref eos_had_base::n0, \ref
	eos_had_base::eoa, and \ref eos_had_base::msom can be computed
	directly, and the compressibility, the skewness, and the
	symmetry energy are computed using the functions
	fkprime_fields() and fesym_fields(). This function overrides
	the generic version in \ref eos_had_base.

	If \ref verbose is greater than zero, then then this function
	reports details on the initial iterations to get the initial
	guess for the solver.
    */
    virtual void saturation();
  
    /** \brief Calculate symmetry energy assuming the field
	equations have already been solved
	
	This may only work at saturation density and may assume
	equal neutron and proton masses.
    */
    double fesym_fields(double sig, double ome, double nb);

    /** \brief Calculate the compressibility assuming the field
	equations have already been solved
	
	This may only work at saturation density and may assume
	equal neutron and proton masses.
    */
    double fcomp_fields(double sig, double ome, double nb);

    /** \brief Calculate compressibilty and \c kprime assuming the field
	equations have already been solved

	This may only work at saturation density and may assume
	equal neutron and proton masses.
	
	\todo This function, \ref o2scl::eos_had_rmf::fkprime_fields() is
	currently untested.
    */
    void fkprime_fields(double sig, double ome, double nb,
		       double &k, double &kprime);
    //@}

    /// \name Fields and field equations
    //@{
    /** \brief A function for solving the field equations

	The values <tt>x[0], x[1]</tt>, and <tt>x[2]</tt> should be
	set to \f$ \sigma, \omega \f$ , and \f$ \rho \f$ on input (in
	\f$ \mathrm{fm}^{-1} \f$ ) and on exit, <tt>y[0], y[1]</tt>
	and <tt>y[2]</tt> contain the field equations and are zero
	when the field equations have been solved.
    */
    int field_eqs(size_t nv, const ubvector &x, ubvector &y);

    /** \brief A function for solving the field equations at finite 
	temperature
	
	The values <tt>x[0], x[1]</tt>, and <tt>x[2]</tt> should be
	set to \f$ \sigma, \omega \f$ , and \f$ \rho \f$ on input (in
	\f$ \mathrm{fm}^{-1} \f$ ) and on exit, <tt>y[0], y[1]</tt>
	and <tt>y[2]</tt> contain the field equations and are zero
	when the field equations have been solved.
    */
    int field_eqsT(size_t nv, const ubvector &x, ubvector &y);

    /** \brief Set a guess for the fields for the next call to calc_e(), 
	calc_p(), or saturation()
    */
    virtual int set_fields(double sig, double ome, double lrho) {
      sigma=sig;
      omega=ome;
      rho=lrho;
      guess_set=true;
      return 0;
    }

    /** \brief Return the most recent values of the meson fields 
	
	This returns the most recent values of the meson fields set by
	a call to \ref saturation(), \ref calc_e(), or 
	\ref calc_p(fermion &, fermion &, thermo &).
    */
    int get_fields(double &sig, double &ome, double &lrho) {
      sig=sigma;
      ome=omega;
      lrho=rho;
      return 0;
    }
    //@}

    /// Return string denoting type ("eos_had_rmf")
    virtual const char *type() { return "eos_had_rmf"; }

    /// \name Solver
    //@{
    /** \brief Set class mroot object for use calculating saturation density
     */
    virtual int set_sat_mroot(mroot<mm_funct<>,ubvector,
			      jac_funct<> > &mrx) {
      sat_mroot=&mrx;
      return 0;
    }

    /** \brief The default solver for calculating the saturation 
	density
	
	Used by fn0() (which is called by saturation()) to solve
	saturation_matter_e() (1 variable).
    */
    mroot_hybrids<mm_funct<>,ubvector,ubmatrix,
      jac_funct<> > def_sat_mroot;
    //@}

    /// \name Functions dealing with naturalness
    //@{
        /** \brief Set the coefficients of a eos_had_rmf object to their 
	limits from naturalness

	As given in \ref Muller96 .

	The definition of the vector-isovector field and coupling
	matches what is done here. Compare the Lagrangian above
	with Eq. 10 from the reference.

	The following couplings should all be of the same
	size:
	\f[
	\frac{1}{2 c_s^2 M^2}, \frac{1}{2 c_v^2 M^2} 
	\frac{1}{8 c_{\rho}^2 M^2},~\mathrm{and}~\frac{
	\bar{a}_{ijk} M^{i+2 j+2 k-4}}{2^{2 k}}
	\f]
	which are equivalent to 
	\f[
	\frac{m_s^2}{2 g_s^2 M^2}, \frac{m_v^2}{2 g_v^2 M^2} 
	\frac{m_{\rho}^2}{8 g_{\rho}^2 M^2},~\mathrm{and}~\frac{
	a_{ijk} M^{i+2 j+2 k-4}}{g_s^i g_v^{2 j} 
	g_{\rho}^{2 k} 2^{2 k}}
	\f]
	
	The connection the \f$ a_{ijk} \f$ 's and the coefficients 
	that are used here is 
	\f{eqnarray*}
	\frac{b M}{3} g_{\sigma}^3 \sigma^3 &=& a_{300}~\sigma^3
	\nonumber \\
	\frac{c}{4} g_{\sigma}^4 \sigma^4 &=& a_{400}~\sigma^4
	\nonumber \\
	\frac{\zeta}{24} g_{\omega}^4 \omega^4 &=& a_{020}~\omega^4
	\nonumber \\
	\frac{\xi}{24} g_{\rho}^4 \rho^4 &=& a_{002}~\rho^4
	\nonumber \\
	b_1 g_{\rho}^2 \omega^2 \rho^2 &=& a_{011}~\omega^2 \rho^2 
	\nonumber \\
	b_2 g_{\rho}^2 \omega^4 \rho^2 &=& a_{021}~\omega^4 \rho^2 
	\nonumber \\
	b_3 g_{\rho}^2 \omega^6 \rho^2 &=& a_{031}~\omega^6 \rho^2 
	\nonumber \\
	a_1 g_{\rho}^2 \sigma^1 \rho^2 &=& a_{101}~\sigma^1 \rho^2 
	\nonumber \\
	a_2 g_{\rho}^2 \sigma^2 \rho^2 &=& a_{201}~\sigma^2 \rho^2 
	\nonumber \\
	a_3 g_{\rho}^2 \sigma^3 \rho^2 &=& a_{301}~\sigma^3 \rho^2 
	\nonumber \\
	a_4 g_{\rho}^2 \sigma^4 \rho^2 &=& a_{401}~\sigma^4 \rho^2 
	\nonumber \\
	a_5 g_{\rho}^2 \sigma^5 \rho^2 &=& a_{501}~\sigma^5 \rho^2 
	\nonumber \\
	a_6 g_{\rho}^2 \sigma^6 \rho^2 &=& a_{601}~\sigma^6 \rho^2 
	\nonumber
	\f}

	Note that Muller and Serot use the notation 
	\f[
	\frac{\bar{\kappa} g_s^3 }{2} = \frac{\kappa}{2} = b M 
	g_s^3 \qquad \mathrm{and} \qquad
	\frac{\bar{\lambda} g_s^4}{6} = \frac{\lambda}{6}
	= c g_s^4
	\f]
	which differs slightly from the "standard" notation above.

	We need to compare the values of
	\f{eqnarray*}
	&\frac{m_s^2}{2 g_s^2 M^2}, \frac{m_v^2}{2 g_v^2 M^2} 
	\frac{m_{\rho}^2}{8 g_{\rho}^2 M^2},\frac{b}{3},
	\frac{c}{4}
	&
	\nonumber \\
	&\frac{\zeta}{24}, \frac{\xi}{384},
	\frac{b_1}{4 g_{\omega}^2},
	\frac{b_2 M^2}{4 g_{\omega}^4},
	\frac{b_3 M^4}{4 g_{\omega}^6},
	\frac{a_1}{4 g_{\sigma} M},&
	\nonumber \\
	&\frac{a_2}{4 g_{\sigma}^2},
	\frac{a_3 M}{4 g_{\sigma}^3},
	\frac{a_4 M^2}{4 g_{\sigma}^4},
	\frac{a_5 M^3}{4 g_{\sigma}^5},~\mathrm{and}~\frac{a_6 M^4}
	{4 g_{\sigma}^6}\, .&
	\f}

	These values are stored in the variables cs, cw, cr, b, c,
	zeta, xi, b1, etc. in the specified \ref eos_had_rmf object. All
	of the numbers should be around 0.001 or 0.002.

	For the scale \f$ M \f$, \ref mnuc is used.

	\todo I may have ignored some signs in the above, which are
	unimportant for this application, but it would be good to fix
	them for posterity.

    */
    void check_naturalness(eos_had_rmf &re) {
      
      double gs=cs*ms;
      double gw=cw*mw;
      double gr=cr*mr;

      re.cs=0.5/cs/cs/mnuc/mnuc;
      re.cw=0.5/cw/cw/mnuc/mnuc;
      re.cr=0.125/cr/cr/mnuc/mnuc;
      re.b=b/3.0;
      re.c=c/4.0;

      re.zeta=zeta/24.0;
      re.xi=xi/384.0;

      re.b1=b1/gw/gw/4.0;
      re.b2=b2/pow(gw,4.0)/4.0*mnuc*mnuc;
      re.b3=b3/pow(gw,6.0)/4.0*pow(mnuc,4.0);

      re.a1=a1/gs/4.0/mnuc;
      re.a2=a2/pow(gs,2.0)/4.0;
      re.a3=a3/pow(gs,3.0)/4.0*mnuc;
      re.a4=a4/pow(gs,4.0)/4.0*mnuc*mnuc;
      re.a5=a5/pow(gs,5.0)/4.0*pow(mnuc,3.0);
      re.a6=a6/pow(gs,6.0)/4.0*pow(mnuc,4.0);

      return;
    }
    
    /** \brief Provide the maximum values of the couplings assuming
	a limit on naturalness

	The limits for the couplings are function of the nucleon and
	meson masses, except for the limits on \c b, \c c, \c zeta,
	and \c xi which are independent of the masses because of the
	way that these four couplings are defined.
    */
    void naturalness_limits(double value, eos_had_rmf &re) {
      
      double gs=cs*ms;
      double gw=cw*mw;
      double gr=cr*mr;

      re.cs=value*2.0*mnuc*mnuc;
      re.cw=value*2.0*mnuc*mnuc;
      re.cr=value*8.0*mnuc*mnuc;
      re.b=value*3.0;
      re.c=value*4.0;

      re.zeta=value*24.0;
      re.xi=value*384.0;
      
      re.b1=value*gw*gw*4.0;
      re.b2=value*pow(gw,4.0)*4.0/mnuc/mnuc;
      re.b3=value*pow(gw,6.0)*4.0/pow(mnuc,4.0);

      re.a1=value*gs*4.0*mnuc;
      re.a2=value*pow(gs,2.0)*4.0;
      re.a3=value*pow(gs,3.0)*4.0/mnuc;
      re.a4=value*pow(gs,4.0)*4.0/mnuc/mnuc;
      re.a5=value*pow(gs,5.0)*4.0/pow(mnuc,3.0);
      re.a6=value*pow(gs,6.0)*4.0/pow(mnuc,4.0);
      
      return;
    }
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Temporary baryon density
     */
    double n_baryon;
    
    /** \brief Temporary charge density 

	\future Should use eos_had_base::proton_frac instead?
    */
    double n_charge;

    /// \name The meson fields
    //@{
    double sigma, omega, rho;
    //@}

    /// Temperature for solving field equations at finite temperature
    double fe_temp;

    /// For calc_e(), if true, then solve for neutron matter
    bool ce_neut_matter;
    
    /// For calc_e(), if true, then solve for proton matter
    bool ce_prot_matter;

    /// True if a guess for the fields has been given
    bool guess_set;

    /// The solver to compute saturation properties
    mroot<mm_funct<>,ubvector,jac_funct<> > *sat_mroot;
    
    /// The function for fix_saturation()
    int fix_saturation_fun(size_t nv, const ubvector &x, ubvector &y);
    
    /// Compute matter at zero pressure (for saturation())
    virtual int zero_pressure(size_t nv, const ubvector &ex, 
			      ubvector &ey);

    /// The function for calc_e()
    virtual int calc_e_solve_fun(size_t nv, const ubvector &ex, 
				 ubvector &ey);

    /// The function for calc_temp_e()
    virtual int calc_temp_e_solve_fun(size_t nv, const ubvector &ex, 
				      ubvector &ey);

    /** \brief Calculate the \c cr coupling given \c sig and \c ome 
	at the density 'nb'.
	
	Used by fix_saturation().
    */
    int calc_cr(double sig, double ome, double nb);

    /// Temperature storage for calc_temp_e()
    double ce_temp;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
