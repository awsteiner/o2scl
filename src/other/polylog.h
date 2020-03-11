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

      Note that the GSL definition of the Fermi-Dirac integral
      includes an additional factor of the Gamma function which is not
      included here.
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
    else if (x>std::numeric_limits<fp_t>::max_exponent) res=0;
    else res=pow(x,a)/(1.0+exp(x-mu));
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

  /** \brief Bose-Einstein integral by integration

      This class performs direct computation of the
      Bose-Einstein integral
      \f[
      F_{a}(\mu) = \int_0^{\infty} \frac{x^a}{\exp^{x-\mu}-1} \, .
      \f]
   */
  template<class inte_t, class fp_t=double> class bose_einstein_integ_tl {
    
  protected:
  
  /// Internal function type
  typedef std::function<fp_t(fp_t)> func_t;
  
  /** \brief The Bose-Einstein function
   */
  fp_t obj_func(fp_t x, fp_t a, fp_t mu) {
    fp_t res;
    if (x==0.0) res=0;
    else if (x>std::numeric_limits<fp_t>::max_exponent) res=0;
    else res=pow(x,a)/(exp(x-mu)-1.0);
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
	      (&bose_einstein_integ_tl::obj_func),
	      this,std::placeholders::_1,a,mu);
    iiu.integ_iu_err(f,0.0,res,err);
    return;
  }
  
  };

  /** \brief Exponentially scaled modified Bessel function of the
      second kind by integration
      
      \warning Not working yet. 

      This class uses an integral representation of the exponentially
      scaled modified Bessel function of the second kind
      \f[
      K_n(z) \exp^{z} = \frac{\sqrt{\pi} z^{n}}{2^{n} \Gamma(n+1/2)}
      \int_1^{\infty} e^{z(1-t)} 
      \left(t^2-1\right)^{n-1/2}~dt
      \f]
      (see
      http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/07/01/01/)
      
  */
  template<class inte_t, class fp_t=double> class bessel_K_exp_integ_tl {

  protected:
  
  /// Internal function type
  typedef std::function<fp_t(fp_t)> func_t;
  
  /** \brief The exponentially scaled modified Bessel integrand
   */
  fp_t obj_func(fp_t t, size_t n, fp_t z) {
    fp_t res;
    fp_t arg=(1-t)*z;
    if (arg<std::numeric_limits<fp_t>::min_exponent) return 0;
    res=exp(arg)*pow(t*t-1,n-0.5);
    if (!std::isfinite(res)) {
      std::cout << t << " x " << n << " " << z << " "
                << res << " " << -t*z+z << " " << t*t-1 << std::endl;
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
    void calc_err(size_t n, fp_t z, fp_t &res, fp_t &err) {
      func_t f=std::bind(std::mem_fn<fp_t(fp_t,size_t,fp_t)>
			 (&bessel_K_exp_integ_tl::obj_func),
			 this,std::placeholders::_1,n,z);
      iiu.integ_iu_err(f,1.0,res,err);
      fp_t fact=o2scl_const::root_pi*pow(z/2,n)/boost::math::tgamma(n+0.5);
      res*=fact;
      err*=fact;
      return;
    }
    
  };
  
  /** \brief Compute the fermion integrals for a non-relativistic
      particle using the GSL functions

      This class is used in \o2p in \ref o2scl::fermion_nonrel to
      compute the Fermi-Dirac integrals for non-relativistic
      fermions.
   */
  class fermi_dirac_integ_gsl {

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

    /** \brief Fermi-Dirac integral of order \f$ 2 \f$
     */
    double calc_2(double y) {
      return gsl_sf_fermi_dirac_int(2,y)*2.0;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    double calc_3(double y) {
      return gsl_sf_fermi_dirac_int(3,y)*6.0;
    }
    
  };

  /** \brief Desc
   */
  class bessel_K_exp_integ_gsl {
    
  public:

    double K1exp(double x) {
      return gsl_sf_bessel_Kn_scaled(1.0,x);
    }
    
    double K2exp(double x) {
      return gsl_sf_bessel_Kn_scaled(2.0,x);
    }
    
    double K3exp(double x) {
      return gsl_sf_bessel_Kn_scaled(3.0,x);
    }
    
  };
  
#if defined(O2SCL_LD_TYPES) || defined(DOXYGEN)

#ifdef O2SCL_NEW_BOOST_INTEGRATION

  /** \brief Compute the fermion integrals for a non-relativistic
      particle by directly integrating in long double precision
   */
  class fermi_dirac_integ_direct {

  protected:
    
    /** \brief The integrator
     */
    fermi_dirac_integ_tl<o2scl::inte_exp_sinh_boost
      <funct_ld,15,long double>,long double> it;
    
  public:

    fermi_dirac_integ_direct() {
      it.iiu.tol_rel=1.0e-17;
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
    
    /** \brief Fermi-Dirac integral of order \f$ 2 \f$
     */
    double calc_2(double y) {
      long double y2=y, res, err;
      it.calc_err(2.0L,y2,res,err);
      return ((double)res);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    double calc_3(double y) {
      long double y2=y, res, err;
      it.calc_err(3.0L,y2,res,err);
      return ((double)res);
    }
    
    /** \brief Polylogarithm function

	\note This currently only works for negative y, even though
	the polylogarithm functions are well defined for \f$ y \in
	[0,1] \f$ .

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

  /** \brief Desc
   */
  class bessel_K_exp_integ_direct {
    
  protected:
    
    /** \brief The integrator
     */
    bessel_K_exp_integ_tl<o2scl::inte_exp_sinh_boost
      <funct_ld,15,long double>,long double> it;
    
  public:

    double K1exp(double x) {
      long double x2=x, res, err;
      it.calc_err(1.0L,x2,res,err);
      return ((double)res);
    }
    
    double K2exp(double x) {
      long double x2=x, res, err;
      it.calc_err(2.0L,x2,res,err);
      return ((double)res);
    }
    
    double K3exp(double x) {
      long double x2=x, res, err;
      it.calc_err(3.0L,x2,res,err);
      return ((double)res);
    }
    
  };
  
  /** \brief Polylogarithm function

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
      The relationship between the polylogarithm and the 
      Bose-Einstein distribution is:
      \f[
      \mathrm{Li}_{1+s}(e^{\mu}) = \frac{1}{\Gamma(s+1)} 
      \int_0^{\infty} \frac{k^{s}}{e^{k-\mu}-1}
      \f]
      or 
      \f[
      \mathrm{Li}_{s}(z) = \frac{1}{\Gamma(s)} 
      \int_0^{\infty} \frac{k^{s-1}}{e^{k-\ln(z)}-1}
      \f]
  */
  class polylog {

  protected:
    
    /** \brief The integrator for negative arguments
     */
    fermi_dirac_integ_tl<o2scl::inte_exp_sinh_boost
      <funct_ld,15,long double>,long double> it;
    
    /** \brief The integrator for positive arguments
     */
    bose_einstein_integ_tl<o2scl::inte_exp_sinh_boost
      <funct_ld,15,long double>,long double> it2;
    
  public:
    
    polylog() {
      it.iiu.tol_rel=1.0e-17;
      it2.iiu.tol_rel=1.0e-17;
    }
    
    /** \brief Polylogarithm function
     */
    double calc(double s, double y) {
      if (y<0.0) {
	// Fermi-Dirac integral representation
	long double a=s-1, mu=log(-y), res, err;
	it.calc_err(a,mu,res,err);
	return -((double)res/boost::math::tgamma(s));
      } else if (y==0.0) {
	return 0.0;
      }
      // Bose-Einstein integral representation
      long double a=s-1, mu=log(y), res, err;
      it2.calc_err(a,mu,res,err);
      return ((double)res/boost::math::tgamma(s));
    }

  };
  
#endif
  
#endif
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
