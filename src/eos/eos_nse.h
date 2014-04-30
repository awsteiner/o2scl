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

      This class computes the composition of matter in nuclear statistical
      equilibrium. The chemical potential of a nucleus X with proton number
      \f$ Z_X \f$ and neutron number \f$ N_X \f$ is given by
      \f[
      \mu_X = N \mu_n + Z \mu_p - E_{\mathrm{bind},X}
      \f]
      where \f$ \mu_n \f$ and \f$ \mu_p \f$ are the neutron and proton
      chemical potentials and \f$ E_{\mathrm{bind},X} \f$ is the binding
      energy of the nucleus. 

      The baryon number density and electron fraction are then given 
      by
      \f[
      n_B = n_{X} (N_X + Z_X) \qquad Y_e n_B = n_X Z_X
      \f]
      where \f$ n_X \f$ is the number density which is determined from
      the chemical potential above. 
 
      This implicitly assumes that the nuclei are non-interacting.

      \future Right now calc_density() needs a very good guess. This 
      could be fixed, probably by solving for the log(mu/T) instead 
      of mu. 
  */
  class nse_eos {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// Parameter structure
    typedef struct {
      double nb, Ye, T;
      nucdist *ndp;
    } solve_parms;
    
    /// Function to solve for baryon and charge conservation
    int solve_fun(size_t nv, const ubvector &x, ubvector &y, 
		  solve_parms &pa);
    
    /// Solver
    mroot<mm_funct<>,ubvector,jac_funct<> > *root;

    /// Compute particle properties assuming classical thermodynamics
    classical cla;

#endif

  public:

    nse_eos();
    
    /** \brief If true, call the error handler if calc_density() does
	not converge (default true)
    */
    bool err_nonconv;

    /** \brief Calculate the equation of state as a function of the
	chemical potentials

	Given \c mun, \c mup and \c T, this computes the composition
	(the individual densities are stored in the distribution \c
	nd) the baryon number density \c nb, and the electron fraction
	\c Ye. 

	This function does not use the solver.
    */
    void calc_mu(double mun, double mup, double T,
		 double &nb, double &Ye, thermo &th, nucdist &nd);
    
    /** \brief Calculate the equation of state as a function of the densities

	Given the baryon number density \c nb, and the electron
	fraction \c Ye and the temperature \c T, this computes the
	composition (the individual densities are stored in the
	distribution \c nd) and the chemical potentials are given in
	\c mun and \c mup .

	This function uses the solver to self-consistently compute
	the chemical potentials. 
     */
    int calc_density(double nb, double Ye, double T, 
		     double &mun, double &mup, thermo &th, nucdist &nd);
    
    /// Default solver 
    mroot_hybrids<mm_funct<>,ubvector,ubmatrix,jac_funct<> > def_root;

    /** \brief Set the solver for use in computing the chemical potentials
     */
    void set_mroot(mroot<mm_funct<>,ubvector,jac_funct<> > &rp) {
      root=&rp;
      return;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
