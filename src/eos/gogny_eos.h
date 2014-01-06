/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_GOGNY_EOS_H
#define O2SCL_GOGNY_EOS_H

#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/hadronic_eos.h>
#include <o2scl/table3d.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Gogny EOS
      
      Gogny EOS from \ref Chappert08 with data kindly supplied
      by Michel Girod. 

      \warning The table might not have a sufficiently dense grid to
      accurately compute derivatives, including compressibility and
      symmetry energy.
  */
  class gogny_eos : public hadronic_eos_eden {
    
  public:

    /// The original EOS data 
    table3d t3d;
    
    gogny_eos() {
    }

    virtual ~gogny_eos() {
    }

    /** \brief Equation of state as a function of density
    */
    virtual int calc_e(fermion &ne, fermion &pr, thermo &th) {
      
      if (t3d.get_nx()==0) {
	O2SCL_ERR_RET("No data loaded in gogny_eos::calc_e().",
		      exc_einval);
      }

      double rho=ne.n+pr.n;
      double asym=(ne.n-pr.n)/(ne.n+pr.n);
      double hc=o2scl_const::hc_mev_fm;

      // Arbitrarily set effective masses to bare masses
      ne.ms=ne.m;
      pr.ms=pr.m;

      // The energy density
      double E=t3d.interp(rho,asym,"enneut")/hc;
      th.ed=E*rho+ne.n*ne.m+pr.n*pr.m;
      
      // The derivatives of the energy densities wrt the baryon density
      double dEdrho=t3d.deriv_x(rho,asym,"enneut")/hc;
      double dEda=t3d.deriv_y(rho,asym,"enneut")/hc;
      ne.mu=ne.m+E+ne.n*(dEdrho+dEda*2.0*pr.n/rho/rho);
      pr.mu=pr.m+E+pr.n*(dEdrho-dEda*2.0*ne.n/rho/rho);

      // The pressure
      th.pr=-th.ed+ne.n*ne.mu+pr.n*pr.mu;

      return 0;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
