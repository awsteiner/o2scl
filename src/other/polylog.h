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
/** \file polylog.h
    \brief File defining \ref o2scl::polylog
*/
#ifndef O2SCL_POLYLOG_H
#define O2SCL_POLYLOG_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <gsl/gsl_specfunc.h>

#include <o2scl/constants.h>
#include <o2scl/err_hnd.h>
#include <o2scl/lib_settings.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/inte_double_exp_boost.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Fermi-Dirac integral by integration

      This class performs direct computation of the
      Fermi-Dirac integral
      \f[
      F_{a}(\mu) = \int_0^{\infty} \frac{x^a}{1+\exp^{x-\mu}} \, .
      \f]
   */
  template<class inte_t, class fp_t=double> class fermi_dirac_integ_tl {

  protected:
  
  /// Internal function type
  typedef std::function<fp_t(fp_t)> func_t;
  
  /** \brief The Fermi-Dirac function
   */
  fp_t obj_func(fp_t x, fp_t a, fp_t mu) {
    fp_t res;
    if (x==0.0) res=0;
    else res=pow(x,a)/(1.0+exp(x-mu));
    if (!std::isfinite(res)) {
      std::cout << x << " " << a << " " << mu << " " << x-mu << " "
		<< res << std::endl;
    }
    return res;
  }

  public:
  
  /** \brief The integrator
   */
  inte_t iiu;
  
  /** \brief Compute the integral, storing the result in 
      \c res and the error in \c err
  */
  void calc_err(fp_t a, fp_t mu, fp_t &res, fp_t &err) {
    func_t f=
    std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
	      (&fermi_dirac_integ_tl::obj_func),
	      this,std::placeholders::_1,a,mu);
    iiu.integ_iu_err(f,0.0,res,err);
    return;
  }
  
  };
  
  /** \brief Compute the fermion integrals for a non-relativistic
      particle using the GSL functions

      This class is used in \o2p in \ref o2scl::fermion_nonrel to
      compute the Fermi-Dirac integrals for non-relativistic
      fermions.
   */
  class fermion_nr_integ_gsl {

  public:
    
    /** \brief Fermi-Dirac integral of order \f$ 1/2 \f$
     */
    double calc_1o2(double y) {
      return gsl_sf_fermi_dirac_half(y)*sqrt(o2scl_const::pi)/2.0;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ -1/2 \f$
     */
    double calc_m1o2(double y) {
      return gsl_sf_fermi_dirac_mhalf(y)*sqrt(o2scl_const::pi);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3/2 \f$
     */
    double calc_3o2(double y) {
      return gsl_sf_fermi_dirac_3half(y)*sqrt(o2scl_const::pi)*0.75;
    }

  };

#if defined(O2SCL_LD_TYPES) || defined(DOXYGEN)

#ifdef O2SCL_NEW_BOOST_INTEGRATION

  /** \brief Compute the fermion integrals for a non-relativistic
      particle by directly integrating in long double precision
   */
  class fermion_nr_integ_direct {

  protected:
    
    /** \brief The integrator
     */
    fermi_dirac_integ_tl<o2scl::inte_exp_sinh_boost
      <funct_ld,15,long double>,long double> it;
    
  public:

    fermion_nr_integ_direct() {
      it.iiu.tol_rel=1.0e-18;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 1/2 \f$
     */
    double calc_1o2(double y) {
      long double y2=y, res, err;
      it.calc_err(0.5L,y2,res,err);
      return ((double)res);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ -1/2 \f$
     */
    double calc_m1o2(double y) {
      long double y2=y, res, err;
      it.calc_err(-0.5L,y2,res,err);
      return ((double)res);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3/2 \f$
     */
    double calc_3o2(double y) {
      long double y2=y, res, err;
      it.calc_err(1.5L,y2,res,err);
      return ((double)res);
    }
    
    /** \brief Polylogarithm function

	\note This currently only works for negative y.

	The relationship between the polylogarithm and the 
	Fermi-Dirac distribution is:
	\f[
	\mathrm{Li}_{1+s}(-e^{\mu}) = - \frac{1}{\Gamma(s+1)} 
	\int_0^{\infty} \frac{k^{s}}{e^{k-\mu}+1}
	\f]
	or 
	\f[
	\mathrm{Li}_{s}(z) = - \frac{1}{\Gamma(s)} 
	\int_0^{\infty} \frac{k^{s-1}}{e^{k-\ln(-z)}+1}
	\f]
     */
    double calc_polylog(double s, double y) {
      long double a=s-1, mu=log(-y), res, err;
      it.calc_err(a,mu,res,err);
      return -((double)res/boost::math::tgamma(s));
    }

  };

#endif
  
#endif
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
