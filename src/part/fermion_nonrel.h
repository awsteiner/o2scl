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
#ifndef O2SCL_NONREL_FERMION_H
#define O2SCL_NONREL_FERMION_H

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/root_cern.h>
#include <o2scl/inte_qagiu_gsl.h>

#include <o2scl/fermion.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Nonrelativistic fermion class

      The rest mass energy density is given by n*m not n*ms. Note that
      the effective mass here is the Landau mass, not the Dirac mass.
      
      Pressure is computed with
      \f[
      P = 2 \varepsilon/3
      \f]
      and entropy density with
      \f[
      s = \frac{5 \varepsilon}{3 T} - \frac{n \mu}{T}
      \f]
      These relations can be verified with an integration by
      parts. See, e.g. \ref Callen pg. 403 or \ref Landau pg. 164.
      
      The functions pair_density() and pair_mu() have not 
      been implemented. 

      \todo Check behaviour of calc_density() at zero density, and
      compare with that from \ref o2scl::fermion_eff, \ref
      o2scl::fermion_rel, and \ref o2scl::classical.

      \todo Implement pair_density() and pair_mu().

      \todo Make sure to test with non-interacting equal to 
      true or false, and document whether or not it works
      with both inc_rest_mass equal to true or false

      \future This could be improved by performing a Chebyshev
      approximation (for example) to invert the density integral so
      that we don't need to use a solver.
  */
  class fermion_nonrel : public fermion_eval_thermo {

  public:

    /// Create a nonrelativistic fermion with mass 'm' and degeneracy 'g'
    fermion_nonrel();
    
    virtual ~fermion_nonrel();
    
    /** \brief Zero temperature fermions
    */
    virtual void calc_mu_zerot(fermion &f);

    /** \brief Zero temperature fermions
    */
    virtual void calc_density_zerot(fermion &f);
    
    /** \brief Calculate properties as function of chemical potential
     */
    virtual void calc_mu(fermion &f, double temper);

    /** \brief Calculate properties as function of density

	If the density is zero, this function just sets part::mu,
	part::nu, part::ed, part::pr, and part::en to zero and returns
	without calling the error handler (even though at 
	zero density and finite	temperature, the chemical potentials
	formally are equal to \f$ -\infty \f$). 
     */
    virtual int calc_density(fermion &f, double temper);

    virtual void pair_mu(fermion &f, double temper) {
      O2SCL_ERR2("Function fermion_nonrel::pair_mu() not ",
		 "implemented.",exc_eunimpl);
      return;
    }

    virtual void pair_density(fermion &f, double temper) {
      O2SCL_ERR2("Function fermion_nonrel::pair_density() not ",
		 "implemented.",exc_eunimpl);
      return;
    }
    
    /// Calculate effective chemical potential from density
    virtual void nu_from_n(fermion &f, double temper);

    /** \brief Set the solver for use in calculating the chemical
	potential from the density 
    */
    int set_density_root(root<funct11> &rp) {
      density_root=&rp;
      return 0;
    }

    /// The default solver for calc_density().
    root_cern<funct11> def_density_root;

    /// Return string denoting type ("fermion_nonrel")
    virtual const char *type() { return "fermion_nonrel"; }

  protected:

#ifndef DOXYGEN_NO_O2NS

    /// Solver to compute chemical potential from density
    root<funct11> *density_root;
    
    /** \brief Function to compute chemical potential from density

	Variable \c nog is the target baryon density divided by
	the spin degeneracy, and \c msT is the effective mass
	times the temperature.
     */
    double solve_fun(double x, double nog, double msT);

  private:

    fermion_nonrel(const fermion_nonrel &);
    fermion_nonrel& operator=(const fermion_nonrel&);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
