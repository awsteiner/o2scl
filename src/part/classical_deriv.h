/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#ifndef O2SCL_CLASSICAL_DERIV_H
#define O2SCL_CLASSICAL_DERIV_H

/** \file classical_deriv.h
    \brief File defining \ref o2scl::classical_deriv_thermo_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include <o2scl/constants.h>
#include <o2scl/part_deriv.h>
#include <o2scl/classical.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a classical particle with derivatives
   */
  template<class fp_t=double>
    class classical_deriv_thermo_tl : public deriv_thermo_base_tl<fp_t> {

  protected:
    
  /// For computing non-derivative properties
  classical_thermo_tl<fp_t> cl;

  public:

  classical_deriv_thermo_tl() {
  }
  
  virtual ~classical_deriv_thermo_tl() {
  }
    
  /** \brief Compute the properties of particle \c p at temperature 
      \c temper from its chemical potential
  */
  virtual void calc_mu(part_deriv &p, fp_t temper) {
      
    cl.calc_mu(p,temper);

    if (temper==0.0) {
      p.dndT=0.0;
      p.dndmu=0.0;
      p.dsdT=0.0;
      return;
    }

    p.dndT=-p.nu/temper/temper*p.n+1.5*p.n/temper;
    p.dndmu=p.n/temper;
    p.dsdT=2.5*p.dndT-p.nu*p.dndT/temper+p.n*p.nu/temper/temper;

    return;
  }


  /** \brief Compute the properties of particle \c p at temperature 
      \c temper from its density
  */
  virtual void calc_density(part_deriv &p, fp_t temper) {

    cl.calc_density(p,temper);

    // Handle zero density first
    if (p.n==0.0 || temper==0.0) {
      p.dndT=0.0;
      p.dndmu=0.0;
      p.dsdT=0.0;
      return;
    }

    p.dndT=-p.nu/temper/temper*p.n+1.5*p.n/temper;
    p.dndmu=p.n/temper;
    p.dsdT=2.5*p.dndT-p.nu*p.dndT/temper+p.n*p.nu/temper/temper;
  
    return;
  }

  /// Return string denoting type ("classical_deriv_thermo")
  virtual const char *type() { return "classical_deriv_thermo"; };

  };

  /** \brief Double-precision version of 
      \ref o2scl::classical_deriv_thermo_tl 
   */
  typedef classical_deriv_thermo_tl<double> classical_deriv_thermo;

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
