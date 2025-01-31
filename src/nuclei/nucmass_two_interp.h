/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2025, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifndef NUCMASS_TWO_INTERP_H
#define NUCMASS_TWO_INTERP_H

/** \file nucmass_two_interp.h
*/

#include <cmath>

#include <o2scl/nucleus.h>
#include <o2scl/nucmass_dz.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/nucmass_fit.h>

namespace o2scl {
  
  /** \brief Desc
  */
  class nucmass_two_interp : public nucmass_fit_base {

  public:

    nucmass_two_interp();

    virtual ~nucmass_two_interp();
    
    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x) {
      return nfb->fit_fun(nv,x);
    }

    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x) {
      return nfb->guess_fun(nv,x);
    }

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess_d(double Z, double N);

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N) {
      return mass_excess_d(Z,N);
    }

    /// Set the nuclear mass formula
    void set_fit_base(nucmass_fit_base &nfb_user) {
      nfb=&nfb_user;
      nfit=nfb->nfit;
      return;
    }

    /// Set the default interpolator
    void set_interpm_base(interpm_base<> &ib_user) {
      ib=&ib_user;
      return;
    }
    
    /// Default nuclear mass formula
    nucmass_dz_fit def_fit;

    /// Default interpolation object
    interpm_idw<> def_ib;

    /// Desc
    void set_default();
    
  protected:

    /// Desc
    nucmass_fit_base *nfb;

    /// Desc
    interpm_base<> *ib;

  };

}

#endif

