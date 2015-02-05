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

      \future Is it better to solve for the log(mu/T) instead of for
      mu/T?
  */
  class eos_nse {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /// Function to solve for baryon and charge conservation
    int solve_fun(size_t nv, const ubvector &x, ubvector &y, 
		  double nB, double Ye, double T,
		  std::vector<o2scl::nucleus> &nd);
    
    /// Solver
    mroot<> *root;

    /// Compute particle properties assuming classical thermodynamics
    classical cla;

#endif

  public:

    eos_nse();
    
    /** \brief If true, call the error handler if calc_density() does
	not converge (default true)
    */
    bool err_nonconv;

    /** \brief Calculate the equation of state as a function of the
	chemical potentials

	Given \c mun, \c mup and \c T, this computes the composition
	(the individual densities are stored in the distribution \c
	nd) the baryon number density \c nB, and the electron fraction
	\c Ye. 

	This function does not use the solver.
    */
    void calc_mu(double mun, double mup, double T, double &nB, 
		 double &Ye, thermo &th, std::vector<nucleus> &nd);
    
    /** \brief Calculate the equation of state as a function of the densities

	Given the baryon number density \c nB in \f$ \mathrm{fm}^{-3}
	\f$, the electron fraction \c Ye and the temperature \c T in
	\f$ \mathrm{fm}^{-1} \f$, this computes the composition (the
	individual densities are stored in the distribution \c nd) and
	the chemical potentials are given in \c mun and \c mup . The
	nuclei in \c nd must have their proton number, neutron number,
	atomic number, binding energy, and spin degeracy already
	specified.

	This function uses the solver to self-consistently compute
	the chemical potentials. 
     */
    int calc_density(double nB, double Ye, double T, double &mun, 
		     double &mup, thermo &th, std::vector<nucleus> &nd);
    
    /// Default solver 
    mroot_hybrids<> def_root;

    /** \brief Set the solver for use in computing the chemical potentials
     */
    void set_mroot(mroot<> &rp) {
      root=&rp;
      return;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
