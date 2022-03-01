/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_EFF_FERMION_H
#define O2SCL_EFF_FERMION_H

/** \file fermion_eff.h
    \brief File defining \ref o2scl::fermion_eff
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
#include <o2scl/misc.h>

#include <o2scl/fermion.h>

namespace o2scl {

  /** \brief Fermion class from fitting method

      \verbatim embed:rst
      Based on the fitting method of [Johns96]_ which is an update
      of the method from [Eggleton73]_. 
      \endverbatim
      
      For a more accurate (but slower) method, use
      \ref o2scl::fermion_rel_tl.

      If the temperature is less than or equal to \ref tlimit (which 
      defaults to zero), the zero-temperature expressions
      from the parent class \ref fermion_zerot_tl are used. 
      
      Given the chemical potential and the temperature the functions
      calc_mu() and pair_mu() work by solving the equation (c.f. Eq. 15
      in Johns et al. (1996))
      \f[
      \psi= 2 \sqrt{1+f/a}+\log\left(\frac{\sqrt{1+f/a}-1}
      {\sqrt{1+f/a}+1}\right)
      \f]
      for \f$ f \f$ given \f$ \psi=(\mu-m)/T \f$.
      If \f$ f/a<10^{-10} \f$, then the alternative expression
      \f[
      \psi= 2 \left[1+f/(2 a)\right]+\log\left\{\frac{f/(2 a)}
      {\left[1+f/(2 a)\right]}\right\}
      \f]
      is used. The pressure, energy density, and entropy, are
      determined as polynomials in \f$ f \f$ with a
      set of precomputed coefficients as done in Johns et al. (1996) .

      If \f$ \psi \f$ is less than \ref min_psi (which defaults to -4)
      then the non-dengenerate approximation from \ref
      fermion_thermo::calc_mu_ndeg() is used. The value of \ref
      min_psi can be decreased to ensure that the expansion is not
      used, but values of \f$ \psi \f$ less than about -200 can cause
      the Johns et al. (1996) procedure outlined above to fail. Values of
      \ref min_psi larger than -4 are not useful.

      When the density and temperature is given instead
      (calc_density() and pair_density()), then there are two ways to
      proceed.
      - Use the density to solve for \f$ f \f$ .
      - Use the density to solve for the chemical potential.

      Because the density is a complicated polynomial in \f$ f \f$,
      the former procedure does not work very well even though it
      might be less time consuming. In this class, the density is
      solved for the effective chemical potential instead. The initial
      guess is just taken from the present value of part::nu .

      \verbatim embed:rst

      .. todo::

         In class fermion_eff:

         - Future: Use bracketing to speed up one-dimensional root finding.

         - Future: Implement the same method of allowing the user to
           turn on or off the degenerate and nondegenerate expansions as
           in \ref o2scl::fermion_rel_tl .

      \endverbatim
  */
  class fermion_eff : public fermion_thermo {

  protected:

    /** \brief The function which solves for the chemical potential
	given the density
     */
    double density_fun(double x, fermion &f, double temper);

    /** \brief The function which solves for the chemical potential
	given the density (including antiparticles)
     */
    double pair_density_fun(double x, fermion &f, double temper);
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    /// Create a fermion with mass \c mass and degeneracy \c dof 
    fermion_eff();
    
    virtual ~fermion_eff();
    
    /** \name Coefficients for finite-temperature approximation
     */
    //@{
    /** \brief Load coefficients
	
	The argument \c ctype Should be one of the constants below.
    */
    void load_coefficients(int ctype);

    /// A set of coefficients from Jim Lattimer 
    static const int cf_fermilat3=1;
    /// The smaller set of coefficients from Johns et al. (1996)
    static const int cf_fermijel2=2;
    /// The larger set of coefficients from Johns et al. (1996)
    static const int cf_fermijel3=3;
    /** \brief The set of coefficients from Johns et al. (1996) which retains 
	better thermodynamic consistency
    */
    static const int cf_fermijel3cons=4;
    //@}

    /** \brief Calculate thermodynamic 
	properties as function of chemical potential
	
	If the quantity \f$ (\mu-m)/T \f$ (or \f$ (\nu-m^{*})/T \f$ in
	the case of interacting particles) is less than -200, then
	this quietly sets the density, the scalar density, the energy
	density, the pressure and the entropy to zero and exits.
    */
    virtual int calc_mu(fermion &f, double temper);

    /** \brief Calculate thermodynamic properties as function of
	density

	\warning This function needs a guess for the chemical
	potential, and will fail if that guess is not sufficiently
	accurate.
    */
    virtual int calc_density(fermion &f, double temper);

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of chemical potential

	\warning This function needs a guess for the chemical
	potential, and will fail if that guess is not sufficiently
	accurate.
    */
    virtual void pair_mu(fermion &f, double temper);

    /** \brief Calculate thermodynamic properties with antiparticles
	as function of density
    */
    virtual int pair_density(fermion &f, double temper);

    /** \brief Set the solver for use in calculating \f$ \psi \f$ 
     */
    int set_psi_root(root<> &rp) {
      psi_root=&rp;
      return 0;
    }

    /** \brief Set the solver for use in calculating the chemical
	potential from the density
    */
    int set_density_root(root<> &rp) {
      density_root=&rp;
      return 0;
    }
    
    /** \brief If the temperature is less than \c tlimit then the
	zero-temperature functions are used (default 0).
    */
    double tlimit;
    
    /** \brief If true, call the error handler when convergence 
	fails (default true)
    */
    bool err_nonconv;

    /** \brief The default solver for \f$ \psi \f$
     */
    root_cern<> def_psi_root;
    
    /** \brief The default solver for calc_density() and pair_density()
     */
    root_cern<> def_density_root;
    
    /// Return string denoting type ("fermion_eff")
    virtual const char *type() { return "fermion_eff"; }

    /// The minimum value of \f$ \psi \f$ (default -200)
    double min_psi;

    /// If true, verify the thermodynamic identity
    bool verify_ti;
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /// The matrix of coefficients
    ubmatrix Pmnf;
    /// The parameter \f$ a \f$
    double parma;
    /// The array row size
    int sizem;
    /// The array column size
    int sizen;
    
    /// The solver for \f$ \psi \f$
    root<> *psi_root;
    /// The other solver for calc_density()
    root<> *density_root;
    
    /// The function which solves for \f$ f \f$ from \f$ \psi \f$.
    double solve_fun(double x, double psi);
    
#endif

  };

}

#endif
