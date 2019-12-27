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
#ifndef O2SCL_CLASSICAL_H
#define O2SCL_CLASSICAL_H

/** \file classical.h
    \brief File defining \ref o2scl::classical_thermo_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/part.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Classical particle class

      \note Note that it makes no sense to include
      \f$ T=0 \f$ functions here
  */
  template<class fp_t=double> class classical_thermo_tl {

  protected:

  /// Desc
  fp_t pi;
    
  public:

  /** \brief Create a classical particle with mass \c m  
      and degeneracy \c g 

      \note This class attempts to handle zero temperature limit
      somewhat gracefully, even though the classical limit doesn't
      necessarily make physical sense there.
  */
  classical_thermo_tl() {
    pi=boost::math::constants::pi<fp_t>();
  }

  virtual ~classical_thermo_tl() {
  }

  /** \brief Calculate properties as function of chemical potential
	
      If the temperature is less than zero, the error handler will
      be called.
  */
  virtual void calc_mu(part_tl<fp_t> &p, fp_t temper) {

    if (temper<0.0) {
      O2SCL_ERR2("Temperature less than zero in ",
		 "classical_thermo::calc_mu().",exc_einval);
    }

    if (p.non_interacting==true) { p.nu=p.mu; p.ms=p.m; }

    // Handle zero temperature case
    if (temper==0.0) {
      if (p.inc_rest_mass) {
	p.n=0.0;
	p.ed=p.n*p.m;
      } else {
	p.n=0.0;
	p.ed=0.0;
      }
      p.pr=0.0;
      p.en=0.0;
      return;
    }

    if (p.inc_rest_mass) {
      if ((p.nu-p.m)/temper<std::numeric_limits<fp_t>::min_exponent10) {
	p.n=0.0;
      } else {
	p.n=exp((p.nu-p.m)/temper)*p.g*pow(p.ms*temper/o2scl_const::pi/2.0,1.5);
      }
      p.ed=1.5*temper*p.n+p.n*p.m;
    } else {
      if (p.nu/temper<std::numeric_limits<fp_t>::min_exponent10) {
	p.n=0.0;
      } else {
	p.n=exp(p.nu/temper)*p.g*pow(p.ms*temper/o2scl_const::pi/2.0,1.5);
      }
      p.ed=1.5*temper*p.n;
    }
    p.pr=p.n*temper;
    p.en=(p.ed+p.pr-p.n*p.nu)/temper;
    return;
  }


  /** \brief Calculate properties as function of density

      If the density or the temperature is less than zero, the error
      handler will be called. In the case of zero density, the
      chemical potential is set to the mass and the energy density,
      pressure, and entropy are set to zero. 
  */
  virtual void calc_density(part_tl<fp_t> &p, fp_t temper) {

    if (p.n<0.0 || temper<0.0) {
      O2SCL_ERR2("Density or temperature less than zero in ",
		 "classical_thermo::calc_density().",exc_einval);
    }
    
    if (p.non_interacting==true) { p.ms=p.m; }

    // Handle zero density first
    if (p.n==0.0) {
      if (p.inc_rest_mass) {
	p.nu=p.m;
      } else {
	p.nu=0.0;
      }
      if (p.non_interacting==true) { p.mu=p.nu; }
      p.ed=0.0;
      p.pr=0.0;
      p.en=0.0;
      return;
    }

    // Handle the zero temperature case
    if (temper==0.0) {
      if (p.inc_rest_mass) {
	p.nu=p.m;
	p.ed=p.n*p.m;
      } else {
	p.nu=0.0;
	p.ed=0.0;
      }
      p.pr=0.0;
      p.en=0.0;
      return;
    }

    if (p.inc_rest_mass) {
      p.nu=p.m+temper*log(p.n/p.g*pow(2.0*o2scl_const::pi/p.ms/temper,1.5));
      p.ed=1.5*temper*p.n+p.n*p.m;
    } else {
      p.nu=temper*log(p.n/p.g*pow(2.0*o2scl_const::pi/p.ms/temper,1.5));
      p.ed=1.5*temper*p.n;
    }
  
    if (p.non_interacting==true) { p.mu=p.nu; }

    p.pr=p.n*temper;
    p.en=(p.ed+p.pr-p.n*p.nu)/temper;

    return;
  }


  /// Return string denoting type ("classical_thermo")
  virtual const char *type() { return "classical_thermo"; }
    
  };

  /** \brief Double-precision version of \ref o2scl::classical_thermo_tl 
   */
  typedef classical_thermo_tl<double> classical_thermo;

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
