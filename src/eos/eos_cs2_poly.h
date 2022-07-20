/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2022, Andrew W. Steiner
  
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
/** \file eos_cs2_poly.h
    \brief File defining \ref o2scl::eos_cs2_poly
*/
#ifndef O2SCL_EOS_CS2_POLY_H
#define O2SCL_EOS_CS2_POLY_H

#include <gsl/gsl_sf_hyperg.h>

#include <o2scl/constants.h>
#include <o2scl/err_hnd.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/eos_tov.h>

namespace o2scl {

  /** \brief An EOS based on a polynomial speed of sound

      \verbatim embed:rst
      Based on [Constantinou17]_.
      \endverbatim
      
      This class constructs an EOS based on a speed of sound of the form
      \f[
      c_s^2 = a_1 + \frac{a_2 n_B^{a_3}}{1+a_4 n_B^{a_3}}
      \f]
      where \f$ n_B \f$ is the baryon number density .
      
      The EOS requires a hypergeometric function which only converges
      under specific conditions on the parameters.
  */
  class eos_cs2_poly : public eos_tov {

  protected:
    
    /** \brief Solver
     */
    root_brent_gsl<> rbg;
    
    /** \brief First speed of sound parameter
     */
    double a1i;
    /** \brief Second speed of sound parameter
     */
    double a2i;
    /** \brief Third speed of sound parameter
     */
    double a3i;
    /** \brief Fourth speed of sound parameter
     */
    double a4i;
    /** \brief Chemical potential integration constant
     */
    double C1;
    /** \brief Energy density integration constant
     */
    double C2;

    /// If true, then a guess for the baryon density has been given
    bool nb_guess_set;

    /// An initial guess when solving for the baryon density
    double nb_guess;
    
  public:

    eos_cs2_poly() {
      C1=0.0;
      C2=0.0;
      nb_guess_set=false;
      nb_guess=0.0;
    }

    void set_nb_guess(double nb) {
      nb_guess=nb;
      nb_guess_set=true;
      return;
    }
    
    /** \brief Fix \f$ a_1 \f$ and \f$ a_2 \f$ based on fitting
	to the sound speed at two different densities
    */
    void fix_params(double nb0, double cs20, double nb1, double cs21,
		    double a3, double a4) {
      a3i=a3;
      a4i=a4;
      double nb0a3=pow(nb0,a3);
      double nb1a3=pow(nb1,a3);
      a1i=(cs21*nb0a3-cs20*nb1a3-a4*cs20*nb0a3*nb1a3+a4*cs21*nb0a3*nb1a3)/
	(nb0a3-nb1a3);
      a2i=(cs20-cs21)*(1.0+a4*nb0a3)*(1.0+a4*nb1a3)/(nb0a3-nb1a3);
      baryon_column=true;
      return;
    }
  
    /** \brief Fix the integration constants by specifying the
	chemical potential at some baryon density and the energy density
	at another baryon density
    */
    void fix_integ_consts(double nb1, double mu1, double nb2, double ed2) {
      C1=mu1*pow(nb1,-a1i)*pow(1.0+a4i*pow(nb1,a3i),-a2i/a3i/a4i);
      double na3=pow(nb2,a3i);
      if (fabs(-a4i*na3)>1.0) {
	O2SCL_ERR2("Fourth argument of hyperg_2F1 greater than 1 ",
		   "in eos_cs2_poly::fix_coeffs().",o2scl::exc_einval);
      }
      C2=((1.0+a1i)*ed2-pow(nb2,1.0+a1i)*C1*
	  gsl_sf_hyperg_2F1((1.0+a1i)/a3i,-a2i/a3i/a4i,(1.0+a1i+a3i)/a3i,
			    -a4i*na3))/(1.0+a1i);
      return;
    }

    /** \brief Return the squared sound speed given the baryon density 
	in \f$ \mathrm{fm}^{-3} \f$
    */
    double cs2_from_nb(double nb) {
      double nba3=pow(nb,a3i);
      return a1i+a2i*nba3/(1.0+a4i*nba3);
    }

    /** \brief Return the chemical potential in \f$ \mathrm{fm}^{-1}
	\f$, including the rest mass, given the baryon density in \f$
	\mathrm{fm}^{-3} \f$
    */
    double mu_from_nb(double nb) {
      return pow(nb,a1i)*pow(1.0+a4i*pow(nb,a3i),a2i/a3i/a4i)*C1;
    }

    /** \brief Return the energy density in \f$ \mathrm{fm}^{-4} \f$,
	including the rest mass energy density, given the baryon density
	in \f$ \mathrm{fm}^{-3} \f$
    */
    double ed_from_nb(double nb) {
      double na3=pow(nb,a3i);
      if (fabs(-a4i*na3)>1.0) {
	O2SCL_ERR2("Fourth argument of hyperg_2F1 greater than 1 ",
		   "in eos_cs2_poly::ed_from_nb().",o2scl::exc_einval);
      }
      return C2+pow(nb,1.0+a1i)*C1*
	gsl_sf_hyperg_2F1((1.0+a1i)/a3i,-a2i/a3i/a4i,1.0+(1.0+a1i)/a3i,
			  -a4i*na3)/(1.0+a1i);
    }

    /** \brief Return the pressure in \f$ \mathrm{fm}^{-4} \f$ 
	given the baryon density 
	in \f$ \mathrm{fm}^{-3} \f$
    */
    double pr_from_nb(double nb) {
      return -ed_from_nb(nb)+nb*mu_from_nb(nb);
    }

    /// Convert ed_from_nb() into a function for the solver
    double ed_from_nb_function(double ed0, double nb) {
      return (ed_from_nb(nb)-ed0)/ed0;
    }
    
    /// Convert pr_from_nb() into a function for the solver
    double pr_from_nb_function(double pr0, double nb) {
      return (pr_from_nb(nb)-pr0)/pr0;
    }
    
    /** \brief Compute the pressure from the energy density using the
        specified guess for the baryon density
    */
    double pr_from_ed_guess(double ed, double nb_guess_loc) {

      // Solve for r prime
      funct fx=std::bind
        (std::mem_fn<double(double,double)>
         (&eos_cs2_poly::ed_from_nb_function),this,ed,std::placeholders::_1);
      rbg.solve(nb_guess_loc,fx);

      // Now compute pressure
      return pr_from_nb(nb_guess_loc);
    }
    
    /** \brief Compute the baryon density from the energy density
        using the specified guess for the baryon density
    */
    double nb_from_ed_guess(double ed, double nb_guess_loc) {

      // Solve for r prime
      funct fx=std::bind
        (std::mem_fn<double(double,double)>
         (&eos_cs2_poly::ed_from_nb_function),this,ed,std::placeholders::_1);
      rbg.solve(nb_guess_loc,fx);

      // Now compute pressure
      return nb_guess_loc;
    }

    /** \brief Compute the energy density from the pressure using the
        specified guess for the baryon density
    */
    double ed_from_pr_guess(double pr, double nb_guess_loc) {

      // Solve for r prime
      funct fx=std::bind
        (std::mem_fn<double(double,double)>
         (&eos_cs2_poly::pr_from_nb_function),this,pr,std::placeholders::_1);
      rbg.solve(nb_guess_loc,fx);

      // Now compute pressure
      return ed_from_nb(nb_guess_loc);
    }

    /** \brief Compute the baryon density from the pressure using the
        specified guess for the baryon density
    */
    double nb_from_pr_guess(double pr, double nb_guess_loc) {

      // Solve for r prime
      funct fx=std::bind
        (std::mem_fn<double(double,double)>
         (&eos_cs2_poly::pr_from_nb_function),this,pr,std::placeholders::_1);
      rbg.solve(nb_guess_loc,fx);

      // Now compute pressure
      return nb_guess_loc;
    }

    /// Compute the energy density from the pressure
    double pr_from_ed(double ed) {
      nb_guess_set=false;
      if (nb_guess_set) {
        double nb=pr_from_ed_guess(ed,nb_guess);
      }
      return pr_from_ed_guess(ed,1.0);
    }
    
    /// Compute the baryon density from the energy density
    double nb_from_ed(double ed) {
      nb_guess_set=false;
      if (nb_guess_set) {
        return nb_from_ed_guess(ed,nb_guess);
      }
      return nb_from_ed_guess(ed,1.0);
    }
    
    /// Compute the energy density from the pressure
    double ed_from_pr(double pr) {
      nb_guess_set=false;
      if (nb_guess_set) {
        return ed_from_pr_guess(pr,nb_guess);
      }
      return ed_from_pr_guess(pr,1.0);
    }
    
    /// Compute the baryon density from the pressure
    double nb_from_pr(double pr) {
      nb_guess_set=false;
      if (nb_guess_set) {
        return nb_from_pr_guess(pr,nb_guess);
      }
      return nb_from_pr_guess(pr,1.0);
    }
    
    /// Compute the energy and baryon densities from the pressure
    virtual void ed_nb_from_pr(double pr, double &ed, double &nb) {
      nb_from_pr(pr);
      ed_from_nb(nb);
      return;
    }
    
  };

}
 
#endif
