/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014-2016, Andrew W. Steiner
  
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
/** \file eos_had_hlps.h
    \brief File defining \ref o2scl::eos_had_hlps
*/
#ifndef O2SCL_HLPS_EOS_H
#define O2SCL_HLPS_EOS_H

#include <iostream>

#include <o2scl/poly.h>
#include <o2scl/fermion.h>
#include <o2scl/eos_had_base.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Schematic EOS from Hebeler et al.

      The energy per baryon is
      \f[
      E/A = \left( 3 \pi^2 n_0/2 \right)^{2/3} \frac{1}{2 M}
      \left\{ \frac{3}{5} \left[ x^{5/3} + (1-x)^{5/3} \right] (2 u)^{2/3}-
      [(2 \alpha - 4 \alpha_L) x (1-x)+\alpha_L] u +
      \left[ (2 \eta - 4 \eta_L) x (1-x) + \eta_L \right] u^{\gamma} \right\}
      \f]
      where \f$ u = n/n_0 \f$ .

      One can fix the values of \f$ \alpha, \eta, \f$ and \f$ \gamma \f$
      by the requirement that the pressure is zero at saturation and
      by fixing the binding energy and incompressibility.

      Note that the original reference has a typo in the pressure in
      Eq. 3. The \f$ 2/5 \f$ factor in front should be \f$ 1/5 \f$ .

      See Ref. \ref Hebeler13eo .
  */
  class eos_had_hlps : public eos_had_eden_base {

  protected:

    /// To solve quadratic equation for 'gamma'
    quadratic_real_coeff_gsl quad;
    
  public:

    /// \name Constants (all unitless)
    //@{
    double gamma;
    double alpha;
    double eta;
    double alphaL;
    double etaL;
    //@}

    eos_had_hlps();

    virtual ~eos_had_hlps() {};

    /** \brief Fix 'alpha', 'eta' and 'gamma' from saturation
	properties

	All inputs must be in \f$ \mathrm{fm}^{-1} \f$. This employs a
	simple iterative method that may not always converge.
    */
    void fix_coeffs(double M, double B, double K);

    /** \brief Fix 'alphaL' and 'etaL' from neutron matter EOS and its
	derivative

	The parameters \c M and \c Eneut must be in
	\f$ \mathrm{fm}^{-1} \f$ and \c dEneut must be in 
	\f$ \mathrm{fm}^{2} \f$
    */
    void fix_neutron_matter(double M, double Eneut, double dEneut);

    /** \brief Equation of state as a function of density
    */
    virtual int calc_e(fermion &ln, fermion &lp, 
		       thermo &lth);

    /// Return string denoting type ("eos_had_hlps")
    virtual const char *type() { return "eos_had_hlps"; }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
