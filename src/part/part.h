/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_PART_H
#define O2SCL_PART_H

#include <string>
#include <iostream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/inte.h>
#include <o2scl/funct.h>
#include <o2scl/mroot.h>

/** \file part.h
    \brief File for definitions for \ref o2scl::thermo and \ref o2scl::part 
*/

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A class holding some useful thermodynamical variables (energy
      density, pressure, entropy density)
  */
  class thermo {

  public:

    /// pressure
    double pr;
    /// energy density
    double ed;
    /// entropy density
    double en;

    /// Return string denoting type ("thermo")
    const char *type() { return "thermo"; }

  };

  /** \brief Addition operator
  */
  extern thermo operator+(const thermo &left, const thermo &right);

  /** \brief Subtraction operator
  */
  extern thermo operator-(const thermo &left, const thermo &right);

  /** \brief Particle base class 
   */
  class part {
    
  public:

    /// Degeneracy (e.g. spin and color if applicable)
    double g;
    /// Mass
    double m;
    /// Number density
    double n;
    /// Energy density
    double ed;
    /// Pressure
    double pr;
    /// Chemical potential
    double mu;
    /// Entropy density
    double en;
    /// Effective mass (Dirac unless otherwise specified)
    double ms;
    /// Effective chemical potential
    double nu;
    /** \brief If true, include the mass in the energy 
	density and chemical potential (default true) 
    */
    bool inc_rest_mass;
    /// True if the particle is non-interacting (default true)
    bool non_interacting;
    
    /// Make a particle of mass \c mass and degeneracy \c dof.
    part(double mass=0.0, double dof=0.0);

    virtual ~part();
  
    /// Set the mass \c mass and degeneracy \c dof.
    virtual void init(double mass, double dof);

    /** \brief Make \c ap an anti-particle with the same mass
	and degeneracy

	This sets the \ref m, \ref g, \ref ms, \ref inc_rest_mass
	and \ref non_interacting fields of \c ap equal to that
	of the current object. If \ref inc_rest_mass is true,
	then it sets 
	\f[
	\mu_{\mathrm{anti}} = - \mu
	\qquad\mathrm{and}\qquad
	\nu_{\mathrm{anti}} = - \nu
	\f]
	and if \ref inc_rest_mass is false, it sets
	\f[
	\mu_{\mathrm{anti}} = - \mu - 2 m
	\qquad\mathrm{and}\qquad
	\nu_{\mathrm{anti}} = - \nu - 2 m
	\f]
    */
    virtual void anti(part &ap);

    /// Return string denoting type ("part")
    virtual const char *type() { return "part"; }
    
  };

  /** \brief Addition operator
  */
  extern thermo operator+(const thermo &left, const part &right);

  /** \brief Subtraction operator
  */
  extern thermo operator-(const thermo &left, const part &right);

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
