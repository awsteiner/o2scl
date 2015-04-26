/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

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

      The function \ref calc_density() can, for low enough
      temperatures, require a very good guess in order to successfully
      solve for the chemical potentials. This is particularly a
      problem also when the typical \f$ Z/A \f$ of the nuclei is not
      close to the desired \f$ Y_e \f$ or the nuclear distribution has
      only a few nuclei.
  */
  class eos_nse {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /// Function to solve for baryon and charge conservation
    int solve_fun(size_t nv, const ubvector &x, ubvector &y, 
		  double nn, double np, double T,
		  std::vector<o2scl::nucleus> &nd);

    /// Desc
    double minimize_fun(size_t nv, const ubvector &x, double T,
			double nn, double np, o2scl::thermo &th,
			std::vector<o2scl::nucleus> &nd);
    
    /// Solver
    mroot<> *mroot_ptr;

    /// Minimizer
    mmin_base<> *mmin_ptr;

    /// Compute particle properties assuming classical thermodynamics
    classical cla;

#endif

  public:

    eos_nse();

    /// Verbosity parameter (default 1)
    int verbose;
    
    /// The maximum number of iterations for \ref make_guess() (default 40)
    size_t make_guess_iters;

    /// Desc
    double make_guess_init_step;
    
    /** \brief If true, call the error handler if calc_density() does
	not converge (default true)
    */
    bool err_nonconv;

    /** \brief Calculate the equation of state as a function of the
	chemical potentials

	Given \c mun, \c mup and \c T, this computes the composition
	(the individual densities are stored in the distribution \c
	nd), the neutron number density \c nn, and the proton number
	density \c np. 

	This function does not use the solver.
    */
    void calc_mu(double mun, double mup, double T, double &nn, 
		 double &np, thermo &th, std::vector<nucleus> &nd);

    /** \brief Find values for the chemical potentials which ensure
	that the densities are within a fixed range

	This function is used by \ref calc_density() to improve
	the initial guesses for the chemical potentials if
	necessary.

	\note This function can fail, for example if the density range
	specified is too small or if the specified distribution
	consists of nuclei with different values of Ye from that
	specified by the density range.
    */
    int make_guess(double &mun, double &mup, double T,
		   thermo &th, std::vector<nucleus> &nd,
		   double nn_min=1.0e-20, double nn_max=1.0e8,
		   double np_min=1.0e-20, double np_max=1.0e8,
		   bool err_on_fail=true);

  /** \brief Calculate the equation of state as a function of the densities

	Given the neutron number density \c nn in \f$ \mathrm{fm}^{-3}
	\f$, the proton number density \c np and the temperature \c T
	in \f$ \mathrm{fm}^{-1} \f$, this computes the composition
	(the individual densities are stored in the distribution \c
	nd) and the chemical potentials are given in \c mun and \c mup
	. The nuclei in \c nd must have their proton number, neutron
	number, atomic number, binding energy, and spin degeracy
	already specified.

	This function uses the solver to self-consistently compute
	the chemical potentials. 
     */
    int calc_density(double nn, double np, double T, double &mun, 
		     double &mup, thermo &th, std::vector<nucleus> &nd);
    
    int direct_solve(double nn, double np, double T, 
		     double &mun, double &mup, thermo &th, 
		     std::vector<nucleus> &nd, bool err_on_fail=true);
    
    int density_min(double nn, double np, double T, 
		    double &mun, double &mup, thermo &th, 
		    std::vector<nucleus> &nd);
    
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
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
