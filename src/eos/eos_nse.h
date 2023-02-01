/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
/** \file eos_nse.h
    \brief File defining \ref o2scl::eos_nse
*/
#ifndef O2SCL_NSE_EOS_H
#define O2SCL_NSE_EOS_H 

#include <o2scl/classical.h>
#include <o2scl/constants.h>
#include <o2scl/nucdist.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/mmin_simp2.h>

namespace o2scl {

  /** \brief Equation of state for nuclei in statistical equilibrium

      This class computes the composition of matter in nuclear
      statistical equilibrium. The chemical potential of a nucleus X
      with proton number \f$ Z_X \f$ and neutron number \f$ N_X \f$ is
      given by
      \f[
      \mu_X = N \mu_n + Z \mu_p - E_{\mathrm{bind},X}
      \f]
      where \f$ \mu_n \f$ and \f$ \mu_p \f$ are the neutron and proton
      chemical potentials and \f$ E_{\mathrm{bind},X} \f$ is the
      binding energy of the nucleus. The chemical potentials are
      assumed to be in units of \f$ \mathrm{fm}^{-1} \f$.

      The baryon number density and electron fraction are then given 
      by
      \f[
      n_B = \sum_X n_{X} (N_X + Z_X) \qquad Y_e n_B = \sum_X n_X Z_X
      \f]
      where \f$ n_X \f$ is the number density which is determined from
      the chemical potential above. 
 
      The nuclei in specified in the parameter named \c nd, must have
      their proton number, neutron number, mass number, binding
      energy, and spin degeracy already specified. This class
      implicitly assumes that the nuclei are non-interacting and that
      the values of <tt>part::inc_rest_mass</tt> are false. The
      chemical potential arguments also do not include the rest mass.
      The nuclear rest mass is presumed to be \f$ Z_X m_p + N_X m_n
      \f$. 

      The function \ref calc_density() attempts to solve for the
      neutron and proton chemical potentials given the neutron and
      proton densities. However, this is relatively difficult. At low
      enough temperatures, \f$ n(\mu) \f$ is a staircase-like function
      with alernating regions which are very flat and or nearly
      vertical. For this reason, derivative-based methods often fail
      without extremely good guesses. The current method of solution
      combines \ref make_guess(), \ref density_min() and \ref
      direct_solve() in order to obtain the solution.

      Note also that \ref calc_density() will fail if there are
      no nuclei in the distribution which equal, or surround the
      requested value of \f$ Y_e=n_p/(n_n+n_p) \f$ determined from \c nn and 
      \c np .
  */
  class eos_nse {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /// Function to solve to match neutron and proton densities
    int solve_fun(size_t nv, const ubvector &x, ubvector &y, 
		  double nn, double np, double T,
		  std::vector<o2scl::nucleus> &nd);

    /// Function to minimize to match neutron and proton densities
    double minimize_fun(size_t nv, const ubvector &x, double T,
			double nn, double np, o2scl::thermo &th,
			std::vector<o2scl::nucleus> &nd);
    
    /// Solver
    mroot<> *mroot_ptr;

    /// Minimizer
    mmin_base<> *mmin_ptr;

    /// Compute particle properties assuming classical thermodynamics
    classical_thermo cla;

#endif

  public:

    eos_nse();

    /// \name Basic usage
    //@{
    /// Verbosity parameter (default 1)
    int verbose;
    
    /** \brief If true, call the error handler if calc_density() does
	not converge (default true)
    */
    bool err_nonconv;

    /** \brief Calculate the equation of state as a function of the
	chemical potentials

	Given \c mun, \c mup and \c T, this computes the composition
	(the individual densities are stored in the distribution \c
	nd), the neutron number density \c nn, and the proton number
	density \c np. Note that the densities can be infinite if
	the chemical potentials are sufficiently large.

	This function does not use the solver or the minimizer.
    */
    void calc_mu(double mun, double mup, double T, double &nn, 
		 double &np, thermo &th, std::vector<nucleus> &nd);

  /** \brief Calculate the equation of state as a function of the densities

	Given the neutron number density \c nn in \f$ \mathrm{fm}^{-3}
	\f$, the proton number density \c np and the temperature \c T
	in \f$ \mathrm{fm}^{-1} \f$, this computes the composition
	(the individual densities are stored in the distribution \c
	nd) and the chemical potentials are given in \c mun and \c mup
	. The nuclei in \c nd must have their proton number, neutron
	number, atomic number, binding energy, and spin degeracy
	already specified.
	
	This function uses \ref make_guess(), \ref direct_solve(),
	and \ref density_min(), to self-consistently compute the
	chemical potentials.
    */
    int calc_density(double nn, double np, double T, double &mun, 
		     double &mup, thermo &th, std::vector<nucleus> &nd);
    //@}

    /// \name Tools for fixing chemical potentials from the densities
    //@{
    /// The maximum number of iterations for \ref make_guess() (default 60)
    size_t make_guess_iters;

    /** \brief The initial stepsize for the chemical potentials relative
	to the temperature (default \f$ 10^5 \f$ )
    */
    double make_guess_init_step;
    
    /** \brief Find values for the chemical potentials which ensure
	that the densities are within a fixed range

	This function improves initial guesses for the chemical
	potentials in order to ensure the densities are within a
	specified range. It can sometimes even succeed when the
	chemical potentials are so far off as to make the densities
	infinite or zero. This function is used by \ref calc_density()
	to improve the initial guesses for the chemical potentials if
	necessary.

	The algorithm can fail in several different ways. This is
	particularly likely if the density range specified by \c
	nn_min, \c nn_max, \c np_min, and \c np_max is small. This
	function ignores the value of \ref err_nonconv, and throws an
	exception on failure only if \c err_on_fail is true (which is
	the default).
    */
    int make_guess(double &mun, double &mup, double T,
		   thermo &th, std::vector<nucleus> &nd,
		   double nn_min=1.0e-20, double nn_max=1.0e8,
		   double np_min=1.0e-20, double np_max=1.0e8,
		   bool err_on_fail=true);

    /** \brief Obtain chemical potentials from densities directly
	using a solver

	This function often requires extremely good guesses for the
	chemical potentials, especially at low temperatures.
    */
    int direct_solve(double nn, double np, double T, 
		     double &mun, double &mup, thermo &th, 
		     std::vector<nucleus> &nd);
    
    /** \brief Obtain chemical potentials from densities 
	using a minimizer
	
	This function often requires extremely good guesses for the
	chemical potentials, especially at low temperatures. By
	default, this calls the minimizer five times, as this seems to
	improve convergence using the default minimizer. By default,
	the value of \ref o2scl::mmin_base::err_nonconv is set to
	false for \ref def_mmin .
    */
    int density_min(double nn, double np, double T, 
		    double &mun, double &mup, thermo &th, 
		    std::vector<nucleus> &nd);
    //@}

    /// \name Numerical methods
    //@{
    /// Default solver 
    mroot_hybrids<> def_mroot;
    
    /// Default minimizer
    mmin_simp2<> def_mmin;
    
    /** \brief Set the solver for use in \ref direct_solve()
     */
    void set_mroot(mroot<> &rp) {
      mroot_ptr=&rp;
      return;
    }

    /** \brief Set the minimizer for use in \ref density_min()
     */
    void set_mmin(mmin_base<> &mp) {
      mmin_ptr=&mp;
      return;
    }
    //@}
    
  };

}

#endif
