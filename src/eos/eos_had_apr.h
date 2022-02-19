/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
/** \file eos_had_apr.h
    \brief File defining \ref o2scl::eos_had_apr
*/
#ifndef O2SCL_APR_EOS_H
#define O2SCL_APR_EOS_H

#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/fermion.h>
#include <o2scl/fermion_nonrel.h>

namespace o2scl {

  /** \brief EOS from Akmal, Pandharipande, and Ravenhall

      \verbatim embed:rst
      The EOS of Akmal, Pandharipande, and Ravenhall, from 
      [Akmal98eo]_ (APR).
      \endverbatim

      The Hamiltonian is:
      \f[
      {\cal H}_{APR} = {\cal H}_{kin} + {\cal H}_{pot}
      \f]

      \f[
      {\cal H}_{kin} = \left\{ \frac{\hbar^2}{2 m} + 
      \left[ p_3 + \left( 1 - x \right) 
      p_5 \right] n e^{-p_4 n} \right\} \tau_n + 
      \left\{ \frac{\hbar^2}{2 m} + \left[ p_3 + x 
      p_5 \right] n e^{-p_4 n} \right\} \tau_p
      \f]

      \f[
      {\cal H}_{pot} =
      g_1 \left[ 1 - \left( 1 - 2 x \right)^2 \right] +
      g_2 \left( 1 - 2 x \right)^2
      \f]

      The following are definitions for \f$ g_i \f$ in the low-density
      phase (LDP) or the high-density phase (HDP):

      \f[
      g_{1,LDP} = -n^2 \left[ p_1 + p_2 n + p_6 n^2 + 
      \left( p_{10} + p_{11} n \right) e^{-p_9^2 n^2} \right] 
      \f]

      \f[
      g_{2,LDP} = -n^2 \left[ \frac{p_{12}}{n} + p_7 + p_8 n + 
      p_{13} e^{-p_9^2 n^2} \right]
      \f]

      \f[
      g_{1,HDP} = g_{1,LDP} -n^2 \left[ p_{17} \left( n - p_{19} \right)
      + p_{21} \left( n - p_{19} \right)^2 e^{p_{18} 
      \left( n - p_{19} \right) } 
      \right]
      \f]

      \f[
      g_{2,HDP} = g_{2,LDP} -n^2 \left[ p_{15} \left( n - p_{20} \right)
      + p_{14} \left( n - p_{20} \right)^2 e^{p_{16} 
      \left( n - p_{20} \right)} \right]
      \f] 

      \note APR seems to have been designed to be used with
      non-relativistic neutrons and protons with equal masses of 939
      MeV. This gives a saturation density very close to 0.16.
  
      The variables \f$ \nu_n\f$ and \f$ \nu_p\f$ contain the
      expressions \f$ (-\mu_n+V_n)/T \f$ and \f$ (-\mu_p+V_p)/T \f$
      respectively, where \f$ V \f$ is the potential part of the
      single particle energy for particle i (i.e. the derivative of
      the Hamiltonian w.r.t. density while energy density held
      constant). Equivalently, \f$ \nu_n\f$ is just \f$ -k_{F_n}^2/ 2
      m^{*} \f$.
      
      The selection between the LDP and HDP is controlled by the value
      of \ref pion. The default is to use the LDP at densities below
      0.16 \f$ \mathrm{fm}^{-3} \f$, and for larger densities to just
      use whichever minimizes the energy.

      \verbatim embed:rst
      The finite temperature approximations from [Prakash97ca]_ 
      are used in testing.
      \endverbatim

      \note Since this EOS uses the effective masses and chemical
      potentials in the fermion class, the values of
      part::non_interacting for neutrons and protons are set to false
      in many of the functions.

      \note The parameter array is unit indexed, so that 
      <tt>par[0]</tt> is unused. This choice makes the connection
      between the code and the paper a bit more transparent.
      
      \future There might be room to improve the testing
      of the finite temperature \part a bit.
      \future There is some repetition between calc_e() and calc_temp_e() 
      that possibly could be removed.
      \future Include the analytic relations from Constantinou et al. 
  */
  class eos_had_apr : public eos_had_temp_eden_base {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Non-relativistic fermion thermodyanmics
    fermion_nonrel nrf;
    
    /// The variable indicating which parameter set is to be used
    int choice;
    
    /// Storage for the parameters
    double par[22];

    /// An integer to indicate which phase was used in calc_e()
    int lp;

#endif

  public:

    /** \brief Create an EOS object with the default parameter 
	set (\f$ A18 + \mathrm{UIX}^{*}+\delta v \f$).
    */
    eos_had_apr();

    virtual ~eos_had_apr();

    /// \name Choice of phase
    //@{
    /** \brief use LDP for densities less than 0.16 and for higher
	densities, use the phase which minimizes energy (default)
    */
    static const int best=0;
    /// LDP (no pion condensation)
    static const int ldp=1;
    /// HDP (pion condensation)
    static const int hdp=2;
    /// Choice of phase (default \ref best)
    int pion;
    /** \brief Return the phase of the most recent call to calc_e()

	This function always returns either \ref eos_had_apr::ldp or
	\ref eos_had_apr::hdp .
     */
    int last_phase() { return lp; }
    //@}
    
    /// \name Basic usage
    //@{
    /** \brief Equation of state as a function of density
     */
    virtual int calc_e(fermion &n, fermion &p, thermo &th);
    
    /// Equation of state as a function of densities
    virtual int calc_temp_e(fermion &n, fermion &pr, double temper, 
			    thermo &th);
    
    /** \brief Compute the compressibility in nuclear
	(isospin-symmetric matter)

	See general notes at eos_had_base::fcomp(). This computes the
	compressibility (at fixed proton fraction = 0.5) exactly,
	unless \ref parent_method is true in which case the derivative
	is taken numerically in eos_had_base::fcomp().
     */
    double fcomp_nuc(double nb);

    /** \brief Calculate symmetry energy of matter as energy of 
	neutron matter minus the energy of nuclear matter

	This function returns the energy per baryon of neutron matter
	minus the energy per baryon of nuclear matter. This will
	deviate significantly from the results from fesym() only if
	the dependence of the symmetry energy on \f$ \delta \f$ is not
	quadratic.
    */
    double fesym_diff(double nb);
    //@}

    /// \name Model selection
    //@{
    /** \brief Select model

	Valid values for \c model_index are: \n
	1 - A18+UIX*+deltav (preferred by Akmal, et. al. - this is 
	the default) \n
	2 - A18+UIX* \n
	3 - A18+deltav \n
	4 - A18 \n

	If any other integer is given, the error handler is called.

	\note This cannot be virtual because it is called by the
	constructor
    */
    void select(int model_index=1);
    /// With three body forces and relativistic corrections
    static const int a18_uix_deltav=1;
    /// With three body forces
    static const int a18_uix=2;
    /// With relativistic corrections
    static const int a18_deltav=3;
    /// No three body forces or relativistic corrections
    static const int a18=4;
    /** \brief Get the value of one of the model parameters
     */
    double get_par(size_t n) {
      if (n>=23 || n==0) {
	O2SCL_ERR("Index out of range in eos_had_apr::set_par().",
		  exc_einval);
      }
      return par[n];
    }

    /** \brief Set the value of one of the model parameters
     */
    void set_par(size_t n, double x) {
      if (n>=23 || n==0) {
	O2SCL_ERR("Index out of range in eos_had_apr::set_par().",
		  exc_einval);
      }
      par[n]=x;
      return;
    }
    //@}

    /// \name Other functions
    //@{
    /** \brief Calculate Q's for semi-infinite nuclear matter
    
	For general discussion, see the documentation to eos_had_base::qs().

	For APR, we set \f$ x_1=x_2=0 \f$ so that \f$ Q_i=P_i/2 \f$ and then
	\f{eqnarray*}
	P_1 &=& \left(\frac{1}{2} p_3-p_5 \right) e^{-p_4 n}
	\nonumber \\
	P_2 &=& \left(\frac{1}{2} p_3+p_5 \right) e^{-p_4 n}
	\f}

	This gives
	\f{eqnarray*}
	Q_{nn}&=&\frac{1}{4} e^{-p_4 \rho}
	\left[ -6 p_5 - p_4 (p_3 - 2 p_5) (n_n + 2 n_p) \right]
	\nonumber \\
	Q_{np}&=&\frac{1}{8} e^{-p_4 \rho}
	\left[ 4 (p_3 - 4 p_5) - 3 p_4 (p_3 - 2 p_5) (n_n + n_p)\right]
	\nonumber \\
	Q_{pp}&=&\frac{1}{4} e^{-p_4 \rho}
	\left[ -6 p_5 - p_4 (p_3 - 2 p_5) (n_p + 2 n_n) \right]
	\f}
    */
    int gradient_qij2(double nn, double np, 
		      double &qnn, double &qnp, double &qpp, 
		      double &dqnndnn, double &dqnndnp,
		      double &dqnpdnn, double &dqnpdnp,
		      double &dqppdnn, double &dqppdnp);

    /// Return string denoting type ("eos_had_apr")
    virtual const char *type() { return "eos_had_apr"; }
    //@}

    /** \brief If true, use the methods from eos_had_base for \ref
	fcomp() and \ref fesym_diff() (default true)

	This can be set to true to check the difference in the
	compressibility and symmety energy between the exact
	expressions and the numerical values from class \ref
	o2scl::eos_had_base.

	\future This variable is probably unnecessary, as the
	syntax
	\code
	eos_had_apr apr;
	ccout << apr.eos_had_base::fcomp(0.16) << endl;
	\endcode
	works just as well.
    */
    bool parent_method;
    
  };

}

#endif
