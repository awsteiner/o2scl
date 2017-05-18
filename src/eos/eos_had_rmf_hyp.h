/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2012, Andrew W. Steiner
  
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

#ifndef O2SCL_EOS_HAD_RMF_H
#define O2SCL_EOS_HAD_RMF_H

#include <string>
#include <cmath>
#include <o2scl/lib_settings.h>
#include <o2scl/constants.h>
#include <o2scl/mm_funct.h>

#include <o2scl/part.h>
#include <o2scl/hadronic_eos.h>
#include <o2scl/fermion.h>

#ifndef DOXYGENP
namespace o2scl {
#endif
  
  /** \brief Relativistic mean field theory EOS with hyperons
      
      Based on \ref Glendenning91, but generalized for higher-order
      couplings as in \ref eos_had_rmf .
      
      \anchor Glendenning91 Glendenning91:
      \htmlonly
      <a href="http://dx.doi.org/">
      N.K. Glendenni and S.A. Moszkowski</a>,
      \endhtmlonly
      \latexonly
      \href{http://dx.doi.org/}{
      N.K. Glendenni and S.A. Moszkowski},
      \endlatexonly
      Phys. Rev. Lett. \b 67, 1805 (1991).
      
   */
  class eos_had_rmf_hyp : public eos_had_rmf {

  public:

    /// \name Hyperon-meson couplings
    //@{
    double xs;
    double xw;
    double xr;
    //@}

    /// If true, include cascade hyperons (default true)
    bool inc_cascade;
    
    eos_had_rmf_hyp();
    
    /** \brief Equation of state and meson field equations 
	as a function of chemical potentials
    */
    int calc_eq_p
      (fermion &ne, fermion &pr, fermion &lam, fermion &sigp, fermion &sigz, 
       fermion &sigm, fermion &casz, fermion &casm, double sig, double ome, 
       double lrho, double &f1, double &f2, double &f3, thermo &lth);

    /** \brief Compute xs assuming a fixed value of the \f$ \Lambda \f$
	binding energy in nuclear matter in \f$ \mathrm{fm}^{-1} \f$
     */
    void calc_xs(double lam_be);

    /** \brief Compute xs assuming a fixed value of the \f$ \Lambda \f$
	binding energy in nuclear matter in \f$ \mathrm{fm}^{-1} \f$
     */
    void calc_xw(double lam_be);

#ifndef DOXYGEN_INTERNAL

  protected:

#endif

  };

#ifndef DOXYGENP
}
#endif

#endif
