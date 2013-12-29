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
/* 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifndef O2SCL_INTE_GSL_H
#define O2SCL_INTE_GSL_H

#include <limits>

#include <gsl/gsl_machine.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief GSL integration base
      
      This base class does not perform any actual integration, but
      just provides functions to be used in the integration
      classes based on GSL.
  */
  class inte_gsl {
    
  public:

    inte_gsl() {
    }
    
  protected:

    /** \brief QUADPACK's nonlinear rescaling of the absolute-error 
	estimate.
	  
        The values \f$ \rho_{\mathrm{abs}} \f$ (stored in
	<tt>result_abs</tt>) and \f$ \rho_{\mathrm{abs}} \f$ (stored
	in <tt>result_asc</tt>) are assumed to be
	\f{eqnarray*}{
	\rho_{\mathrm{abs}} &=& \int_a^b |f|\,dx, \\
	\rho_{\mathrm{asc}} &=& \int_a^b |f - \mu(f)|\, dx, \qquad
	\mu(f) = \frac{1}{b-a}\int_a^b f\, dx,
	\f}
	all of which are computed from the best (i.e., finest-grid)
	approximation of the integrals.  The rescaled error, \f$
	\sigma_\mathrm{err}, \f$ is computed from the raw error, \c
	err, by
	\f[
	\sigma_\mathrm{err} =
	\rho_\mathrm{asc} \cdot 
	\min \left\{1, \; 
	\left(\frac{200 |\mathrm{err}|}{\rho_\mathrm{asc}} \right)^{3/2}
	\right\},
	\f]
	or
	\f[
	\sigma_\mathrm{err} = 
	50\cdot \epsilon_\mathrm{mach} \cdot \rho_\mathrm{abs},
	\f]
	whichever of the two is greater. The value \f$
	\epsilon_\mathrm{mach} \f$ denotes "machine epsilon." (In the
	case that the second value underflows, the first value is
	automatically accepted.)

	This function is used in \ref inte_qng_gsl and \ref
	inte_kronrod_gsl::gauss_kronrod_base().
    */
    double rescale_error(double err, const double result_abs, 
			 const double result_asc) {

      err=fabs(err);
      
      if (result_asc != 0 && err != 0) {
	
	double scale=pow((200*err/result_asc),1.5);
	
	if (scale < 1) {
	  err=result_asc*scale;
	} else {
	  err=result_asc;
	}
      }

#ifdef O2SCL_CPP11
      double dbl_eps=std::numeric_limits<double>::epsilon();
      double dbl_min=std::numeric_limits<double>::min();
#else 
      double dbl_eps=GSL_DBL_EPSILON;
      double dbl_min=GSL_DBL_MIN;
#endif

      if (result_abs > dbl_min/(50*dbl_eps)) {

	double min_err=50*dbl_eps*result_abs;
	
	if (min_err > err) {
	  err=min_err;
	}
      }
      
      return err;
    }
    
  };
  
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
