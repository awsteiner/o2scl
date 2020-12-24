/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
/** \file eos_had_schematic.h
    \brief File defining \ref o2scl::eos_had_schematic
*/
#ifndef O2SCL_SCHEMATIC_EOS_H
#define O2SCL_SCHEMATIC_EOS_H

#include <iostream>
#include <cmath>
#include <o2scl/eos_had_base.h>

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

      Symmetry energy at nuclear matter density is \f$ a+b \f$. 
  */
  class eos_had_schematic : public eos_had_eden_base {

  public:

    /** \brief The kinetic energy symmetry coefficient in inverse fm 
        (default \f$ 17~\mathrm{MeV}~/(\hbar c) \f$)

        The default value corresponds to an effective mass of about
        0.7. 
    */
    double a;
    
    /** \brief The potential energy symmetry coefficient in inverse 
        fm (default \f$ 13~\mathrm{MeV}~/(\hbar c) \f$) 
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
    
    eos_had_schematic();

    virtual ~eos_had_schematic() {};

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

        \verbatim embed:rst
        .. todo:: 

           - In eos_had_schematic::set_a_from_mstar(): 
             This was computed in schematic_sym.nb, which might 
             be added to the documentation?

        \endverbatim
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
    
    /** \brief Return the baryon number susceptibility, \f$ \partial \mu_B /
        \partial n_B \f$ in \f$ \mathrm{fm}^{2} \f$. 

        \verbatim embed:rst
        .. todo:: 

           - This function, eos_had_schematic::baryon_suscep() 
             is untested.

        \endverbatim
    */
    virtual double baryon_suscep(double n, double x) {
      double alpha=n*(1.0-2.0*x);
      double ret=(kpp*n*(5.0-2.0*n0)*(n-n0)*(n-n0)+
                  18.0*n0*(3.0*comp*n*(3.0*n-2.0*n0)*n0)+
                  kprime*n*(2.0*n*n-3.0*n*n0+n0*n0)+3.0*pow(n0,3.0)*alpha*
                  (-1.0*a*pow(n/n0,2.0/3.0)*(-6.0+alpha)+
                   9*b*pow(n/n0,gamma)*(2.0+alpha*(-1.0+gamma))))/
        (486.0*n*pow(n0,4.0));
      return ret;
    }

    /// Return string denoting type ("eos_had_schematic")
    virtual const char *type() { return "eos_had_schematic"; }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
