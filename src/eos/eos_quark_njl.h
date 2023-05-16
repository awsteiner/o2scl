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
/** \file eos_quark_njl.h
    \brief File defining \ref o2scl::eos_quark_njl
*/
#ifndef O2SCL_EOS_QUARK_NJL_H
#define O2SCL_EOS_QUARK_NJL_H

#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/eos_quark.h>
#include <o2scl/quark.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/inte.h>
#include <o2scl/inte_qag_gsl.h>

namespace o2scl {

  /** \brief Nambu Jona-Lasinio EOS

      \verbatim embed:rst
      This class is based on [Buballa99]_, but the notation is a bit
      different. The effective mass is denoted :math:`m^{*}` rather
      than :math:`m`, and the bare mass is denoted :math:`m` rather
      than :math:`m_0`.
      \endverbatim

      The quantities \ref G, and \ref K are the coupling constants and
      \ref L (i.e. \f$ \Lambda \f$) is the three-momentum cutoff. In
      order to use the EOS, the user should either (i) set the bag
      constant, \ref B0 directly, or (ii) use \ref set_parameters() to
      modify the parameters (and then the \ref set_parameters()
      function also automatically computes the bag constant).

      This class can compute the EOS from the quark condensates
      (stored in \ref o2scl::quark::qq) by setting \ref from_qq to
      <tt>true</tt> (this is the default) or from the dynamical masses
      (denoted \f$ m_i^{*} \f$ and 
      stored in \ref o2scl::part_tl::ms) by setting \ref from_qq to
      <tt>false</tt>.

      The Fermi gas-like contribution to the pressure due
      plus the the contribution from the bag pressure is stored in
      \ref o2scl::part_tl::pr . For the \f$ T=0 \f$ EOS, the 
      energy density for each quark is set so that
      \f$ \varepsilon + P = \mu n \f$ . 

      \note This code does not currently check to see if \f$
      k_F>\Lambda \f$ and thus gives unphysical results at higher
      densities.

      The Lagrangian is
      \f[
      {\cal L} = \bar q ( i \partial{\hskip-2.0mm}/ - {\hat m_0}) q \;+\;
      G \sum_{k=0}^8 [\,({\bar q}\lambda_k q)^2 + ({\bar q}
      i\gamma_5\lambda_k q)^2\,] + {\cal L}_6
      \f]
      \f[
      {\cal L}_6 = - K \,[ \,{\rm det}_f ({\bar
      q}(1+\gamma_5) q) + {\rm det}_f ({\bar q}(1-\gamma_5) q) \,] \, .
      \f]

      <b>Finite T documentation</b>

      This implementation includes contributions from antiquarks.
      
      <b>Zero temperature</b>

      In the mean-field approximation, 
      the corresponding thermodynamic potential at \f$ T=0 \f$ is
      \f[
      \Omega = \Omega_{FG} + \Omega_{\mathrm{int}}
      \f]
      where \f$\Omega_{\mathrm{FG}}\f$ is the Fermi gas contribution and
      \f$\Omega_{\mathrm{int}}\f$ is the contribution from interactions.
      \f$\Omega_{\mathrm{int}}\f$ is
      \f[
      \frac{\Omega_{\mathrm{int}}}{V} = - 2 N_c \sum_{i=u,d,s}
      \int \frac {d^3p}{(2\pi)^3} \sqrt{m_i^{*2} + p^2} +
      \frac{\Omega_{V}}{V} \, .
      \f]
      The lst term is the vacuum contribution, 
      \f$\Omega_{V}\f$:
      \f[
      \frac{\Omega_{V}}{V} = 
      \sum_{i=u,d,s} 2 G \langle\bar{q}_i q_i \rangle^2
      - 4 K \langle \bar{q}_u q_u \rangle \langle \bar{q}_d q_d \rangle
      \langle \bar{q}_s q_s \rangle + B_0\,.
      \f]
      where \f$B_0\f$ is a constant defined to ensure that the 
      energy density and the pressure of the vacuum is zero.

      \verbatim embed:rst
      Unlike [Buballa99]_, the bag constant,
      :math:`\Omega_{\mathrm{int}}/V` is defined without the term
      \endverbatim

      \f[
      \sum_{i=u,d,s} 2 N_C \int_0^{\Lambda} 
      \frac{d^3 p}{(2 \pi)^3} \sqrt{ m_{i}^2+p^2 } ~dp
      \f]
      where \f$ m_i \f$ is the bare mass, since this allows an easier
      comparison to the finite temperature EOS. The constant \f$B_0\f$
      in this case is therefore significantly larger, but the energy
      density and pressure are still zero in the vacuum.

      <b>Finite Temperature</b>

      The thermodynamic potential (including the contribution from
      antiquarks) in the mean-field approximation is
      \f{eqnarray*}
      \Omega &=& 6 \sum_{i=u,d,s}
      \int_0^{\Lambda} \frac {d^3p}{(2\pi)^3} \left\{ E_i^{*} +
      T \log \left[ 1 + \exp (-E^{*}_i+\mu_i)/T \right]
      + T \log \left[ 1 + \exp (-E^{*}_i-\mu_i)/T \right]
      \right\} \nonumber \\
      && -
      \sum_{i=u,d,s} 2 G \langle\bar{q}_i q_i \rangle^2
      + 4 K \langle \bar{q}_u q_u \rangle \langle \bar{q}_d q_d \rangle
      \langle \bar{q}_s q_s \rangle + B_0
      \f}
      where \f$ E_i^{*} \equiv \sqrt{m_i^{* 2}+p^2} \f$ .
      The quark densities are
      \f[
      n_i = 6 
      \int_0^{\Lambda} \frac{d^3p}{(2\pi)^3} 
      \left\{
      \left[ \frac{1}{1 + \exp (E_i^{*}-\mu_i)/T} \right]
      -
      \left[ \frac{1}{1 + \exp (E_i^{*}+\mu_i)/T} \right]
      \right\} \, ,
      \f]
      and the quark condensates are
      \f[
      \langle \bar{q}_i q_i \rangle = - 6
      \int_0^{\Lambda} \frac {d^3p}{(2\pi)^3} 
      \frac{m_i^{*} k^2}{E_i^{*}}
      \left[ 1 - 
      \frac{1}{1 + \exp (E^{*}_i-\mu_i)/T} -
      \frac{1}{1 + \exp (E^{*}_i+\mu_i)/T} \right] \, .
      \f]
      \verbatim embed:rst
      See, e.g. [Hatsuda94]_ or [Alaverdyan21]_.
      \endverbatim

      \verbatim embed:rst
      The Feynman-Hellman theorem [Bernard88]_, gives
      \endverbatim
      \f[
      \left< \bar{q} q \right> = \frac{\partial m^{*}}{\partial m}
      \f]

      The functions calc_e() and calc_p() never return a value other
      than zero, but will give nonsensical results for nonsensical
      inputs.

      <b>References</b>

      \verbatim embed:rst
      This class was created for [Steiner00]_. 
      See also [Bernard88]_, [Hatsuda94]_,
      [Buballa99]_, and [Alaverdyan21]_.

      .. todo::

         - In class eos_quark_njl: better documentation.
         - future Remove the stored quark pointers if they are
           unnecessary?

      \endverbatim

  */
  class eos_quark_njl : public eos_quark {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Class constructor
    eos_quark_njl();

    /// \name Settings
    //@{
    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief If true, call the error handler when the solver fails
        to converge (default true)
    */
    bool err_on_fail;
    
    /** \brief Accuracy limit for Fermi integrals for finite temperature
        (default 20)

	\ref limit is used for the finite temperature integrals to
	ensure that no numbers larger than <tt>exp(limit)</tt> or
	smaller than <tt>exp(-limit)</tt> are avoided.
    */
    double limit;

    /** \brief Determine variables for solver (default true)
        
	If this is false, then the internal solvers use the effective
	masses as variables, otherwise, the quark condensates are used
	as variables. The two forms are equivalent, but having both
        versions can be useful for testing.
     */
    bool from_qq;
    //@}

    /// \name NJL model parameters
    //@{
    /// The momentum cutoff, \f$ \Lambda \f$ (in \f$ \mathrm{fm}^{-1} \f$)
    double L;
    
    /// The four-fermion coupling (in \f$ \mathrm{fm}^{2} \f$)
    double G;
    
    /** \brief The 't Hooft six-fermion interaction coupling
	(in \f$ \mathrm{fm}^{5} \f$)
    */
    double K;
    
    /// The bag constant (in \f$ \mathrm{fm}^{-4} \f$)
    double B0;
    
    /** \brief Set the parameters and the bag constant \c B0 
      
	This function allows the user to specify the momentum cutoff,
	\c lambda, the four-fermion coupling \c fourferm 
	and the six-fermion coupling from the 't Hooft interaction
	\c sixferm . If 0.0 is given for any of the values, then
	the default is used (\f$ \Lambda=602.3/(\hbar c),
	G=1.835/\Lambda^2, K=12.36/\Lambda^5 \f$).

	The value of the shift in the bag constant \c B0 is
	automatically calculated to ensure that the energy density and
	the pressure of the vacuum are zero. The functions
	set_quarks() and set_thermo() can be used before hand to
	specify the \ref quark and \ref thermo objects.
    */
    virtual int set_parameters(double lambda=0.0, double fourferm=0.0, 
			       double sixferm=0.0);
    //@}

    /// \name Zero-temperature EOS
    //@{
    /** \brief Equation of state as a function of chemical potentials
	
	This function automatically solves the gap equations.
    */
    virtual int calc_p(quark &u, quark &d, quark &s, thermo &lth);

    /** \brief Equation of state and gap equations as a function of 
	chemical potential
    */
    virtual int calc_eq_p(quark &u, quark &d, quark &s, double &gap1, 
			  double &gap2, double &gap3, thermo &lth);
    
    /** \brief Equation of state and gap equations as a function of
	the densities
    */
    virtual int calc_eq_e(quark &u, quark &d, quark &s, double &gap1, 
			  double &gap2, double &gap3, thermo &lth);

    /** \brief Compute the thermodynamic potential at \f$ T=0 \f$
     */
    virtual double f_therm_pot(double qqu, double qqd, double qqs,
                               double msu, double msd, double mss,
                               bool vac_terms=true);
    //@}

    /// \name Finite-temperature EOS
    //@{
    /** \brief Equation of state as a function of chemical potentials at 
	finite temperature

	This function automatically solves the gap equations.
    */
    virtual int calc_temp_p(quark &u, quark &d, quark &s, 
			    double T, thermo &th);

    /** \brief Equation of state and gap equations
	as a function of chemical potentials at finite temperature
    */
    int calc_eq_temp_p(quark &tu, quark &td, quark &ts,
		       double &gap1, double &gap2, double &gap3,
		       thermo &qb, double temper);

    /** \brief Compute the thermodynamic potential at \f$ T>0 \f$
     */
    virtual double f_therm_pot_T(double qqu, double qqd, double qqs,
                                 double msu, double msd, double mss,
                                 double T, bool vac_terms=true);
    //@}

    /// \name Gap equations
    //@{
    /** \brief Calculates gap equations in \c y as a function of the 
	constituent masses in \c x
	
	The function utilizes the \ref quark objects which can
	be specified in set_quarks() and the \ref thermo object
	which can be specified in eos::set_thermo().
    */
    int gap_func_ms(size_t nv, const ubvector &x, ubvector &y);
    
    /** \brief Calculates gap equations in \c y as a function of the 
	quark condensates in \c x

	The function utilizes the \ref quark objects which can
	be specified in set_quarks() and the \ref thermo object
	which can be specified in eos::set_thermo().
    */
    int gap_func_qq(size_t nv, const ubvector &x, ubvector &y);
    
    /** \brief Calculates gap equations in \c y as a function of the 
	constituent masses in \c x

	The function utilizes the \ref quark objects which can
	be specified in set_quarks() and the \ref thermo object
	which can be specified in eos::set_thermo().
    */
    int gap_func_ms_T(size_t nv, const ubvector &x, ubvector &y, double T);
    
    /** \brief Calculates gap equations in \c y as a function of the 
	quark condensates in \c x

	The function utilizes the \ref quark objects which can
	be specified in set_quarks() and the \ref thermo object
	which can be specified in eos::set_thermo().
    */
    int gap_func_qq_T(size_t nv, const ubvector &x, ubvector &y, double T);
    //@}

    /** \name The default quark masses (in inverse femtometers)

        \verbatim embed:rst
        These default masses are taken from [Buballa99]_, where they
        were used to fix the pion and kaon decay constants, and the
        pion, kaon, and eta prime masses.
        \endverbatim

	They are set in the constructor. The default values are 5.5
	MeV for the up and down quark and 140.7 MeV for the strange
	quark (then divided by \ref o2scl_const::hc_mev_fm for the
	conversion).
    */
    //@{
    double up_default_mass;
    double down_default_mass;
    double strange_default_mass;
    //@}

    /** \brief Set the quark objects to use

	The quark objects are used in gap_func_ms(), gap_func_qq(), 
	gap_func_ms_T(), gap_func_qq_T(), and B0_func().
    */
    int set_quarks(quark &u, quark &d, quark &s);

    /** \name The default quark objects
	
	The masses are automatically set in the constructor to
	\c up_default_mass, \c down_default_mass, and 
	\c strange_default_mass.c
    */
    //@{
    quark def_up; 
    quark def_down;
    quark def_strange;
    //@}

    /// Return string denoting type ("eos_quark_njl")
    virtual const char *type() { return "eos_quark_njl"; }

    /// \name Numerical methods
    //@{
    /// Set solver to use in set_parameters()
    virtual int set_solver
      (mroot<mm_funct,boost::numeric::ublas::vector<double>,
       jac_funct> &s) {
      solver=&s;
      return 0;
    }

    /// Set integration object
    virtual int set_inte(inte<> &i) {
      it=&i;
      return 0;
    }

    /// The default solver
    mroot_hybrids<mm_funct,boost::numeric::ublas::vector<double>, 
      boost::numeric::ublas::matrix<double>,jac_funct> def_solver;
    
    /// The default integrator
    inte_qag_gsl<> def_it;
    //@}
    
  protected:

#ifndef DOXYGEN_INTERNAL

    /// \name Protected functions and variables
    //@{
    /// The integrator for finite temperature integrals
    inte<> *it;
    
    /// The solver to use for set_parameters()
    mroot<mm_funct,boost::numeric::ublas::vector<double>, 
      jac_funct> *solver;
    
    /// Used by calc_B0() to compute the bag constant
    int B0_func(size_t nv, const ubvector &x, ubvector &y);
    
    /** \brief Calculates the contribution to the bag constant from 
        quark \c q 

        This function evaluates the integral
        \f[
        I = 6 \int_0^{\Lambda} \frac{d^3 p}{(2 \pi)^3} \sqrt{m^{*2}+p^2}~dp
        \f]

        \verbatim embed:rst
        This is part of the bag pressure (at zero temperature) as
        given in Eq. 5 of [Buballa99]_.
        \endverbatim

        This function is used in \ref calc_eq_p() and \ref calc_eq_e() .
     */
    void njl_bag(quark &q);

    /// The up quark
    quark *up;
    /// The down quark
    quark *down;
    /// The strange quark
    quark *strange;

    /// The integrand for the quark condensate
    double integ_qq(double x, double T, double mu, double m,
                    double ms);
    /// The integrand for the density
    double integ_density(double x, double T, double mu, double m,
                         double ms);
    /// The integrand for the energy density
    double integ_edensity(double x, double T, double mu, double m,
                          double ms);
    /// The integrand for the pressure
    double integ_pressure(double x, double T, double mu, double m,
                          double ms);
    //@}
  
#endif
    
  };

  /** \brief The Nambu-Jona-Lasinio model with vector interactions

      \verbatim embed:rst
      This class is based on [Alaverdyan21]_, but see also
      older work e.g. [Klimt90]_. We use the \ref part_tl::nu
      field to store what is referred to in [Alaverdyan21]_
      as :math:`{\tilde{\mu}}`. 
      \endverbatim
   */
  class eos_quark_njl_vec : public eos_quark_njl {
    
  public:

    /// The vector coupling constant
    double GV;

    /// \name The zero temperature EOS
    //@{
    /** \brief Compute the equation of state as a function of the
        chemical potentials
	
	This function automatically solves the gap equations.
    */
    virtual int calc_p(quark &u, quark &d, quark &s, thermo &lth);

    /** \brief Compute the gap equations and the equation of state
        as a function of the chemical potentials
     */
    virtual int calc_eq_p_vec(quark &tu, quark &td, quark &ts,
                              double &gap1, double &gap2, double &gap3,
                              double &vec1, double &vec2, double &vec3,
                              thermo &th);

    /** \brief Compute the thermodynamic potential at 
        zero temperature
     */
    virtual double f_therm_pot_vec(double qqu, double qqd, double qqs,
                               double msu, double msd, double mss,
                               double nuu, double nud, double nus,
                               bool vac_terms=true);

    //@}

    /// \name The finite temperature EOS
    //@{
    /** \brief Equation of state as a function of chemical potentials at 
	finite temperature

	This function automatically solves the gap equations.
    */
    virtual int calc_temp_p(quark &u, quark &d, quark &s, 
			    double T, thermo &th);

    /** \brief Compute the gap equations and the equation of state at
        finite temperature as a function of the chemical potentials
     */
    virtual int calc_eq_temp_p(quark &tu, quark &td, quark &ts,
                               double &gap1, double &gap2,
                               double &gap3, double &vec1, double &vec2,
                               double &vec3, thermo &th, double temper);
    
    /** \brief Compute the thermodynamic potential at \f$ T>0 \f$
     */
    virtual double f_therm_pot_T_vec(double qqu, double qqd, double qqs,
                                     double msu, double msd, double mss,
                                     double nuu, double nud, double nus,
                                     double T, bool vac_terms=true);
    //@}
    
    /// \name The gap equations
    //@{
    /** \brief The gap equations starting from the effective masses
     */
    virtual int gap_func_ms_vec(size_t nv, const ubvector &x, ubvector &y);

    /** \brief The gap equations starting from the effective masses
     */
    virtual int gap_func_qq_vec(size_t nv, const ubvector &x, ubvector &y);

    /** \brief Calculates gap equations in \c y as a function of the 
	quark condensates in \c x

	The function utilizes the \ref quark objects which can
	be specified in set_quarks() and the \ref thermo object
	which can be specified in eos::set_thermo().
    */
    int gap_func_qq_T(size_t nv, const ubvector &x, ubvector &y, double T);

    /** \brief Calculates gap equations in \c y as a function of the 
	constituent masses in \c x

	The function utilizes the \ref quark objects which can
	be specified in set_quarks() and the \ref thermo object
	which can be specified in eos::set_thermo().
    */
    int gap_func_ms_T(size_t nv, const ubvector &x, ubvector &y, double T);
    //@}
    
  };
  
}

#endif
