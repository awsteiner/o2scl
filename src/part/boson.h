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
#ifndef O2SCL_BOSON_H
#define O2SCL_BOSON_H

/** \file boson.h
    \brief File defining \ref o2scl::boson
*/

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>

#include <o2scl/part.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Boson class
   */
  class boson : public part {
    
  public:

    /// Create a boson with mass \c mass and degeneracy \c dof 
    boson(double mass=0.0, double dof=0.0);

    /** \brief The condensate
	
	The condensate variable is provided principally for
	user storage and is mostly ignored by \o2p classes. 
    */
    double co;

    /** \brief Calculate properties of massless bosons 
	
	The expressions used are exact. The chemical potentials are
	ignored.
    */
    virtual void massless_calc(double temper);
    
    /// Return string denoting type ("boson")
    virtual const char *type() { return "boson"; }

  };

  /** \brief Compute the thermodynamic properties of a boson 
      [abstract base]
   */
  class boson_thermo {
    
  public:
    
    /** \brief Calculate thermodynamic properties as function of
	chemical potential
    */
    virtual void calc_mu(boson &b, double temper)=0;

    /** \brief Calculate thermodynamic properties as function of
	density
    */
    virtual void calc_density(boson &b, double temper)=0;

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of chemical potential
    */
    virtual void pair_mu(boson &b, double temper)=0;

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of density
    */
    virtual void pair_density(boson &b, double temper)=0;

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
