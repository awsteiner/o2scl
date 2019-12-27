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
#ifndef O2SCL_BOSON_EFF_H
#define O2SCL_BOSON_EFF_H

/** \file boson_eff.h
    \brief File defining \ref o2scl::boson_eff
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/constants.h>
#include <o2scl/funct.h>
#include <o2scl/mm_funct.h>
#include <o2scl/root.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/root_cern.h>
#include <o2scl/mroot_hybrids.h>

#include <o2scl/boson.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Boson class from fitting method
     
      Based on the fitting method of \ref Johns96 which is an update
      of the method from \ref Eggleton73 . This method is approximate,
      but very fast. For a more accurate (but slower) method, use
      \ref o2scl::boson_rel.
      
      Given the chemical potential and the temperature the functions
      \ref calc_mu() and \ref pair_mu() work by solving the equation
      (c.f. Eq. 26 in \ref Johns96)
      \f[
      \psi = \frac{h}{(h+\sqrt{a})} - \ln \left( 
      \frac{h+\sqrt{a}}{\sqrt{a}}\right)
      \f]
      for \f$ h \f$ given \f$ \psi=(\mu-m)/T \f$. The pressure, energy
      density, and entropy, are determined as polynomials in \f$ h \f$
      with a set of precomputed coefficients as done in \ref Johns96 .

      When the density and temperature is given instead (\ref
      calc_density() and \ref pair_density()), then there are two ways
      to proceed:
      - use the density to solve for \f$ f \f$ , or
      - use the density to solve for the chemical potential.

      Because the density is a complicated polynomial in \f$ f \f$,
      the former procedure does not work very well (the polynomial
      produces spurious solutions) even though it might be less time
      consuming. In this class, the density is solved for the
      effective chemical potential instead. The initial guess is just
      taken from the present value of part::mu or, if 
      part::non_interacting is false, from part::nu .
  */
  class boson_eff : public boson_thermo {

  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    /// Create a boson with mass \c m and degeneracy \c g 
    boson_eff();
    
    virtual ~boson_eff();
  
    /** \brief Load coefficients for finite-temperature approximation
	
	Presently acceptable values of \c fn are: 
	\c bosejel21, \c bosejel22, \c bosejel34, and
	\c bosejel34cons from \ref Johns96.
    */
    int load_coefficients(int ctype);
    /// A set of coefficients from \ref Johns96
    static const int cf_bosejel21=2;
    /// A set of coefficients from \ref Johns96
    static const int cf_bosejel22=3;
    /// A set of coefficients from \ref Johns96 (default)
    static const int cf_bosejel34=4;
    /** \brief The set of coefficients from \ref Johns96 which retains 
	better thermodynamic consistency
     */
    static const int cf_bosejel34cons=5;
  
    /** \brief Calculate thermodynamic 
	properties as function of chemical potential
    */
    virtual void calc_mu(boson &b, double temper);

    /** \brief Calculate thermodynamic 
	properties as function of density
    */
    virtual void calc_density(boson &b, double temper);

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of chemical potential
    */
    virtual void pair_mu(boson &b, double temper);

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of density
    */
    virtual void pair_density(boson &b, double temper);

    /** \brief Set the solver for use in calculating \f$ \psi \f$ 
     */
    void set_psi_root(root<> &rp) {
      psi_root=&rp;
      return;
    }

    /** \brief Set the solver for use in calculating the chemical
	potential from the density 
    */
    void set_density_mroot(mroot<> &rp) {
      density_mroot=&rp;
      return;
    }

    /** \brief The default solver for calc_density() and pair_density()
     */
    mroot_hybrids<> def_density_mroot;

    /** \brief The default solver for \f$ \psi \f$
     */
    root_cern<> def_psi_root;

    /// Return string denoting type ("boson_eff")
    virtual const char *type() { return "boson_eff"; }

#ifndef DOXYGEN_INTERNAL

  protected:
  
    /// The coefficients
    ubmatrix Pmnb;
    /// The number of coefficient rows
    int sizem;
    /// The number of coefficient columns
    int sizen;
    /// The parameter, \f$ a \f$
    double parma;
    /// Temporary storage
    double fix_density;

    /// The solver for calc_density()
    mroot<> *density_mroot;

    /// The solver to compute \f$ h \f$ from \f$ \psi \f$.
    root<> *psi_root;

    /// The function which solves for \f$ h \f$ from \f$ \psi \f$.
    double solve_fun(double x, double psi);

    /// Fix density for \ref calc_density()
    int density_fun(size_t nv, const ubvector &x, ubvector &y,
		    boson &b, double T);

    /// Fix density for \ref pair_density()
    int pair_density_fun(size_t nv, const ubvector &x, ubvector &y,
			 boson &b, double T);

#endif
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
