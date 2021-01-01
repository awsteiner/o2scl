/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2021, Andrew W. Steiner
  
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
  class eos_cs2_poly {

  protected:
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
  
  public:

    eos_cs2_poly() {
      C1=0.0;
      C2=0.0;
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

  };

  /** \brief EOS with a constant speed of sound
   */
  class eos_cs2_const {
    
  protected:

    /** \brief Chemical potential integration constant
     */
    double C1;
    /** \brief Energy density integration constant
     */
    double C2;
  
  public:

    /** \brief Speed of sound
     */
    double cs2;
    
    eos_cs2_const() {
      cs2=1.0;
      C1=0.0;
      C2=0.0;
    }
  
    /** \brief Fix the integration constants by specifying the
	pressure, baryon chemical potential, and energy density
    */
    void fix_integ_consts(double mub1, double ed1, double pr1) {
      double nb1=(ed1+pr1)/mub1;
      C1=mub1*pow(nb1,-cs2);
      C2=(ed1*cs2-pr1)/(1.0+cs2);
      return;
    }
    
    /** \brief Return the chemical potential in \f$ \mathrm{fm}^{-1}
	\f$, including the rest mass, given the baryon density in \f$
	\mathrm{fm}^{-3} \f$
    */
    double mu_from_nb(double nb) {
      return C1*pow(nb,cs2);
    }

    /** \brief Return the energy density in \f$ \mathrm{fm}^{-4} \f$,
	including the rest mass energy density, given the baryon density
	in \f$ \mathrm{fm}^{-3} \f$
    */
    double ed_from_nb(double nb) {
      return C1*pow(nb,cs2+1.0)/(1.0+cs2)+C2;
    }

    /** \brief Return the energy density in \f$ \mathrm{fm}^{-4} \f$,
	including the rest mass energy density, given the baryon density
	in \f$ \mathrm{fm}^{-3} \f$
    */
    double nb_from_ed(double ed) {
      return pow((ed-C2)*(1.0+cs2)/C1,1.0/(1.0+cs2));
    }

    /** \brief Return the pressure in \f$ \mathrm{fm}^{-4} \f$ 
	given the baryon density 
	in \f$ \mathrm{fm}^{-3} \f$
    */
    double pr_from_nb(double nb) {
      return C1*cs2*pow(nb,cs2+1.0)/(1.0+cs2)-C2;
    }

  };

}
 
#endif
