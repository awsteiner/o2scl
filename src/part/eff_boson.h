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
#ifndef O2SCL_EFF_BOSON_H
#define O2SCL_EFF_BOSON_H

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
     
     \todo Better documentation (see eff_fermion)
     \todo Remove the 'meth2' stuff
     \todo Remove static variables fix_density and stat_temper
     \todo Remove exit() calls
  */
  class eff_boson {

  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    /// Create a boson with mass \c m and degeneracy \c g 
    eff_boson();
    
    virtual ~eff_boson();
  
    /** \brief Load coefficients for finite-temperature approximation
	
	Presently acceptable values of \c fn are: \c boselat3 from
	Lattimer's notes \c bosejel21, \c bosejel22, \c bosejel34, and
	\c bosejel34cons from \ref Johns96.
    */
    int load_coefficients(int ctype);
    /// A set of coefficients from Jim Lattimer 
    static const int cf_boselat3=1;
    /// A set of coefficients from \ref Johns96
    static const int cf_bosejel21=2;
    /// A set of coefficients from \ref Johns96
    static const int cf_bosejel22=3;
    /// A set of coefficients from \ref Johns96
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
    void set_psi_root(root<funct> &rp) {
      psi_root=&rp;
      return;
    }

    /** \brief Set the solver for use in calculating the chemical
	potential from the density 
    */
    void set_density_mroot(mroot<mm_funct<>,
			   boost::numeric::ublas::vector<double>, 
			   jac_funct<> > &rp) {
      density_mroot=&rp;
      return;
    }

    /** \brief Set the solver for use in calculating the chemical
	potential from the density (meth2=true) 
    */
    void set_meth2_root(root<funct> &rp) {
      meth2_root=&rp;
      return;
    }

    /** \brief The default solver for calc_density() and pair_density()
     */
    mroot_hybrids<mm_funct<>,
      boost::numeric::ublas::vector<double>, 
      boost::numeric::ublas::matrix<double>,
      jac_funct<> > def_density_mroot;

    /** \brief The default solver for \f$ \psi \f$
     */
    root_cern<funct> def_psi_root;

    /** \brief The default solver for calc_density() and pair_density()
     */
    root_cern<funct> def_meth2_root;

    virtual const char *type() { return "eff_boson"; }

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

    /// Desc
    boson *bp;

    /// Desc
    double T;
    
    /// The solver for calc_density()
    mroot<mm_funct<>,boost::numeric::ublas::vector<double>, 
      jac_funct<> > *density_mroot;

    /// The solver to compute \f$ h \f$ from \f$ \psi \f$.
    root<funct> *psi_root;

    /// The solver for calc_density()
    root<funct> *meth2_root;

    /// The function which solves for \f$ h \f$ from \f$ \psi \f$.
    double solve_fun(double x, double &psi);

    /// Fix density for calc_density()
    int density_fun(size_t nv, const ubvector &x, ubvector &y);

    /// Fix density for pair_density()
    int pair_density_fun(size_t nv, const ubvector &x, ubvector &y);

#endif
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
