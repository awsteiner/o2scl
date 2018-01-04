/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCLEUS_H
#define O2SCL_NUCLEUS_H

/** \file nucleus.h
    \brief File defining \ref o2scl::nucleus
*/

#include <o2scl/part.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A simple nucleus class
      
      The variable part::m is typically used for the mass of the
      nucleus with no electrons.

      The binding energy of the nucleus (\ref be) is typically defined
      as the mass of the nucleus (without the electrons or their
      binding energy) minus Z times the mass of the proton minus N
      times the mass of the neutron.

      The mass excess (\ref mex) is defined as the mass of the nucleus
      including the electron contribution (but not including the
      electron binding energy) minus A times the mass of the atomic
      mass unit.

      The variable \ref part::inc_rest_mass is set to \c false by
      default, to insure that energies and chemical potentials do not
      include the rest mass. This is typically appropriate for nuclei.
  */
  class nucleus : public part {

  public:

    nucleus();

    /// Proton number
    int Z;

    /// Neutron number
    int N;

    /// Baryon number
    int A;

    /// Mass excess (typically in \f$ \mathrm{fm}^{-1} \f$)
    double mex;

    /** \brief Binding energy with a minus sign for bound nuclei
	(typically in \f$ \mathrm{fm}^{-1} \f$)
    */
    double be;

    /// Copy constructor
    nucleus(const nucleus &p) {
      g=p.g;
      m=p.m;
      ms=p.ms;
      n=p.n;
      ed=p.ed;
      pr=p.pr;
      mu=p.mu;
      en=p.en;
      nu=p.nu;
      N=p.N;
      Z=p.Z;
      A=p.A;
      mex=p.mex;
      be=p.be;
      inc_rest_mass=p.inc_rest_mass;
      non_interacting=p.non_interacting;
    }

    /// Copy construction with operator=()
    nucleus &operator=(const nucleus &p) {
      if (this!=&p) {
	g=p.g;
	m=p.m;
	ms=p.ms;
	n=p.n;
	ed=p.ed;
	pr=p.pr;
	mu=p.mu;
	en=p.en;
	nu=p.nu;
	N=p.N;
	Z=p.Z;
	A=p.A;
	mex=p.mex;
	be=p.be;
	inc_rest_mass=p.inc_rest_mass;
	non_interacting=p.non_interacting;
      }
      return *this;
    }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
