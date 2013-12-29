/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef O2SCL_SCHEMATIC_EOS_H
#define O2SCL_SCHEMATIC_EOS_H

#include <iostream>
#include <cmath>
#include <o2scl/hadronic_eos.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Schematic hadronic equation of state

      A schematic equation of state defined by the energy density:
      \f[
      \epsilon = n_n m_n + n_p m_p + 
      n \left\{ eoa+\frac{comp}{18}(n/n0-1)^2+
      \frac{kprime}{162}(n/n0-1)^3+
      \frac{kpp}{1944}(n/n0-1)^4+(1- 2 x)^2 
      \left[a \left(\frac{n}{n0}\right)^{2/3}+
      b \left(\frac{n}{n0}\right)^{\gamma} \right] \right\}
      \f]

      Symmetry energy at nuclear matter density is a+b. 

      Note that it doesn't really matter what kind of particle
      object is used, since the calc_e() function doesn't use
      any of the particle thermodynamics functions. 
  */
  class schematic_eos : public hadronic_eos_eden {

  public:

    /** \brief The kinetic energy symmetry coefficient in inverse fm 
	(default \f$ 17~mathrm{MeV}~/(\hbar c) \f$)

	The default value corresponds to an effective mass of about
	0.7. 
    */
    double a;
    
    /** \brief The potential energy symmetry coefficient in inverse 
	fm (default \f$ 13~mathrm{MeV}~/(\hbar c) \f$) 
    */
    double b;

    /** \brief The coefficient of a density to the fourth term in 
	inverse fm (default 0)
    */
    double kpp;
    
    /** \brief The exponent of the high-density symmetry energy 
	(unitless, default 1.0)
    */
    double gamma;
    
    schematic_eos();

    virtual ~schematic_eos() {};

    /** \brief Equation of state as a function of density
    */
    virtual int calc_e(fermion &ln, fermion &lp, thermo &lth);

    /** \brief Set kprime so that the energy per baryon of zero-density 
	matter is zero
    */
    virtual int set_kprime_zeroden() {
      kprime=162.0*eoa+9.0*comp;
      return 0;
    }
    
    /** \brief Set kpp so that the energy per baryon of zero-density 
	matter is zero
    */
    virtual int set_kpp_zeroden() {
      kpp=12.0*kprime-108.0*comp-1944.0*eoa;
      return 0;
    }

    /** \brief Fix the kinetic energy symmetry coefficient from 
	the reduced nucleon effective mass and the saturation density

	This assumes the nucleons are non-relativistic and that the
	neutrons and protons have equal mass. The relativistic
	corrections are around 1 part in \f$ 10^{6} \f$.

	\todo This was computed in schematic_sym.nb, which might be 
	added to the documentation?
    */
    virtual int set_a_from_mstar(double u_msom, double mnuc) {
      a=cbrt(n0*n0*o2scl_const::pi2*o2scl_const::pi2/4.0/3.0)/
	(2.0*u_msom*mnuc);
      return 0;
    }

    /** \brief Return the energy per baryon of matter at zero density

	This is inaccessible from calc_e() so is available separately
	here. Using set_kprime_zeroden() or set_kpp_zeroden() will 
	fix kprime or kpp (respectively) to ensure that this is zero.

	The result provided here does not include the nucleon mass and
	is given in \f$ \mathrm{fm}^{-1} \f$.
    */
    virtual double eoa_zeroden() {
      return eoa+comp/18.0-kprime/162.0+kpp/1944.0;
    }

    /// Return string denoting type ("schematic_eos")
    virtual const char *type() { return "schematic_eos"; }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
