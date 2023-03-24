/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
    \brief File defining various integrals and polylogs
*/
#ifndef O2SCL_POLYLOG_H
#define O2SCL_POLYLOG_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <gsl/gsl_specfunc.h>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <o2scl/constants.h>
#include <o2scl/err_hnd.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/inte_kronrod_boost.h>
#include <o2scl/inte_double_exp_boost.h>
#include <o2scl/exception.h>

namespace o2scl {

  /** \brief Compute a Fermi-Dirac integral by direct integration

      This class performs direct computation of the
      Fermi-Dirac integral
      \f[
      F_{a}(\mu) = \int_0^{\infty} \frac{x^a}{1+e^{x-\mu}} \, .
      \f]
      using an integrator which is specified as a class
      template parameter. This class is used in 
      \ref o2scl::fermi_dirac_integ_direct and \ref o2scl::polylog .

      Note that the GSL definition of the Fermi-Dirac integral
      includes an additional factor of \f$ 1/\Gamma(a+1) \f$
      which is not included here.

      \verbatim embed:rst
      
      .. todo::

         In fermi_dirac_integ_tl, better testing of accuracy.

      \endverbatim
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
      else if (x-mu>std::numeric_limits<fp_t>::max_exponent) res=0;
      else res=pow(x,a)/(1+exp(x-mu));
      return res;
    }
    
  public:
    
    /** \brief The integrator
     */
    inte_t iiu;
    
    /** \brief Compute the integral, storing the result in 
        \c res and the error in \c err
    */
    int calc_err(fp_t a, fp_t mu, fp_t &res, fp_t &err) {
      func_t f=
        std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
                  (&fermi_dirac_integ_tl::obj_func),
                  this,std::placeholders::_1,a,mu);
      int iret=iiu.integ_iu_err(f,0.0,res,err);
      return iret;
    }
    
  };

  /** \brief Compute a Bose-Einstein integral by direct integration

      This class performs direct computation of the
      Bose-Einstein integral
      \f[
      F_{a}(\mu) = \int_0^{\infty} \frac{x^a}{e^{x-\mu}-1} \, .
      \f]
      This class is used in \ref o2scl::polylog .

      \verbatim embed:rst

      .. todo::

         In bose_einstein_integ_tl, better testing of accuracy.

      \endverbatim
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
    int calc_err(fp_t a, fp_t mu, fp_t &res, fp_t &err) {
      func_t f=
        std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
                  (&bose_einstein_integ_tl::obj_func),
                  this,std::placeholders::_1,a,mu);
      int iret=iiu.integ_iu_err(f,0.0,res,err);
      return iret;
    }
    
  };

  /** \brief Compute an exponentially scaled modified Bessel function
      of the second kind by direct integration
      
      This class uses an integral representation of the exponentially
      scaled modified Bessel function of the second kind
      \f[
      K_n(z) e^{z} = \frac{\sqrt{\pi} z^{n}}{2^{n} \Gamma(n+1/2)}
      \int_1^{\infty} e^{z(1-t)} 
      \left(t^2-1\right)^{n-1/2}~dt
      \f]
      (see
      http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/07/01/01/)
      
      \verbatim embed:rst

      .. todo::

         In bessel_K_exp_integ_tl, better testing of accuracy.

      \endverbatim
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
    // 0.4343 is approximately log10(exp(1)) so this value
    // is approximately -706 in double precision
    if (arg<std::numeric_limits<fp_t>::min_exponent10/0.4343) {
      return 0;
    }
    res=exp(arg)*pow(t*t-1,n-0.5);
    if (!isfinite(res)) {
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
    int calc_err(size_t n, fp_t z, fp_t &res, fp_t &err) {
      func_t f=std::bind(std::mem_fn<fp_t(fp_t,size_t,fp_t)>
			 (&bessel_K_exp_integ_tl::obj_func),
			 this,std::placeholders::_1,n,z);
      int iret=iiu.integ_iu_err(f,1.0,res,err);
      fp_t fact=o2scl_const::root_pi*pow(z/2,n)/boost::math::tgamma(n+0.5);
      res*=fact;
      err*=fact;
      return iret;
    }
    
  };
  
  /** \brief Compute several Fermi-Dirac integrals
      useful for non-relativistic fermions using GSL

      This class performs direct computation of the
      Fermi-Dirac integral
      \f[
      F_{a}(\mu) = \int_0^{\infty} \frac{x^a}{1+e^{x-\mu}} \, .
      \f]
      using the functions from GSL where 
      \f$ a \in [-1/2,1/2,3/2,2,3] \f$ .

      Note that the GSL definition of the Fermi-Dirac integral
      includes an additional factor of \f$ 1/\Gamma(a+1) \f$
      which is not included here.

      The GSL functions internally throw overflow exceptions if the
      argument is too large. In this class, the value <tt>inf</tt>
      is returned in those cases instead. 

      This class is used in \o2p in <tt>o2scl::fermion_thermo_tl</tt>,
      <tt>o2scl::fermion_tl</tt> and <tt>o2scl::fermion_nonrel_tl</tt>
      to compute the Fermi-Dirac integrals for non-relativistic
      fermions.
   */
  class fermi_dirac_integ_gsl {

  public:
    
    /** \brief Fermi-Dirac integral of order \f$ 1/2 \f$
     */
    double calc_1o2(double y) {
      double ret;
      try {
	ret=gsl_sf_fermi_dirac_half(y)*sqrt(o2scl_const::pi)/2.0;
      } catch (const o2scl::exc_overflow_error &e) {
	return std::numeric_limits<double>::infinity();
      }
      return ret;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ -1/2 \f$
     */
    double calc_m1o2(double y) {
      double ret;
      try {
	ret=gsl_sf_fermi_dirac_mhalf(y)*sqrt(o2scl_const::pi);
      } catch (const o2scl::exc_overflow_error &e) {
	return std::numeric_limits<double>::infinity();
      }
      return ret;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3/2 \f$
     */
    double calc_3o2(double y) {
      double ret;
      try {
	ret=gsl_sf_fermi_dirac_3half(y)*sqrt(o2scl_const::pi)*0.75;
      } catch (const o2scl::exc_overflow_error &e) {
	return std::numeric_limits<double>::infinity();
      }
      return ret;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 2 \f$
     */
    double calc_2(double y) {
      double ret;
      try {
	ret=gsl_sf_fermi_dirac_int(2,y)*2.0;
      } catch (const o2scl::exc_overflow_error &e) {
	return std::numeric_limits<double>::infinity();
      }
      return ret;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    double calc_3(double y) {
      double ret;
      try {
	ret=gsl_sf_fermi_dirac_int(3,y)*6.0;
      } catch (const o2scl::exc_overflow_error &e) {
	return std::numeric_limits<double>::infinity();
      }
      return ret;
    }
    
  };

  /** \brief Compute exponentially scaled modified bessel function of
      the second kind using the GSL functions

      This class computes \f$ K_n(z) e^z\f$ for \f$ n=1,2,3 \f$.

      \note This class is used in \o2p in
      <tt>o2scl::fermion_thermo_tl</tt>,
      <tt>o2scl::fermion_nonrel_tl</tt>, and
      <tt>o2scl::fermion_rel_tl</tt>
   */
  class bessel_K_exp_integ_gsl {
    
  public:

    /** \brief Compute \f$ K_1(x) e^x \f$
     */
    double K1exp(double x) {
      return gsl_sf_bessel_Kn_scaled(1.0,x);
    }
    
    /** \brief Compute \f$ K_2(x) e^x \f$
     */
    double K2exp(double x) {
      return gsl_sf_bessel_Kn_scaled(2.0,x);
    }
    
    /** \brief Compute \f$ K_3(x) e^x \f$
     */
    double K3exp(double x) {
      return gsl_sf_bessel_Kn_scaled(3.0,x);
    }
    
  };
  
  /** \brief Compute exponentially scaled modified bessel function of
      the second kind using Boost

      This class computes \f$ K_n(z) e^z\f$ for \f$ n=1,2,3 \f$.
   */
  template<class fp_t, class internal_fp_t=fp_t>
  class bessel_K_exp_integ_boost {
    
  public:

    /** \brief 
        
        (This is needed for fermion_rel)
     */
    void set_tol(double t) {
      return;
    }
    
    /** \brief Compute \f$ K_1(x) e^x \f$
     */
    fp_t K1exp(fp_t x) {
      internal_fp_t x2=static_cast<internal_fp_t>(x), ret;
      ret=exp(x2)*boost::math::cyl_bessel_k(1,x2);
      fp_t ret2=static_cast<fp_t>(ret);
      return ret2;
    }
    
    /** \brief Compute \f$ K_2(x) e^x \f$
    */
    fp_t K2exp(fp_t x) {
      internal_fp_t x2=static_cast<internal_fp_t>(x), ret;
      ret=exp(x2)*boost::math::cyl_bessel_k(2,x2);
      fp_t ret2=static_cast<fp_t>(ret);
      return ret2;
    }
    
    /** \brief Compute \f$ K_3(x) e^x \f$
     */
    fp_t K3exp(fp_t x) {
      internal_fp_t x2=static_cast<internal_fp_t>(x), ret;
      ret=exp(x2)*boost::math::cyl_bessel_k(3,x2);
      fp_t ret2=static_cast<fp_t>(ret);
      return ret2;
    }
    
  };
  

  /** \brief Compute several Fermi-Dirac integrals useful for
      non-relativistic fermions by integrating with a higher
      precision type

      This class performs direct computation of the
      Fermi-Dirac integral
      \f[
      F_{a}(\mu) = \int_0^{\infty} \frac{x^a}{1+e^{x-\mu}} \, .
      \f]
      where \f$ a \in [-1/2,1/2,3/2,2,3] \f$ .
      See also \ref o2scl::fermi_dirac_integ_gsl .
      The integration is handled by an object of 
      type \ref o2scl::fermi_dirac_integ_tl using an
      integrator of type \ref o2scl::inte_exp_sinh_boost .

      Note that the GSL definition of the Fermi-Dirac integral
      includes an additional factor of \f$ 1/\Gamma(a+1) \f$
      which is not included here. 

      \verbatim embed:rst
      
      .. todo::

         Future: In fermi_dirac_integ_direct, create a new function
         allowing arbitrary values of 'a' in the equation above.

      \endverbatim
   */
  template <class fp_t=double, class func_t=funct_ld,
	    size_t max_refine=30, 
	    class internal_fp_t=long double>
  class fermi_dirac_integ_direct {

  protected:
    
    internal_fp_t half;
    internal_fp_t three_half;
    internal_fp_t three;
    internal_fp_t two;
    
  public:

    /** \brief The integrator
     */
    fermi_dirac_integ_tl<o2scl::inte_exp_sinh_boost
                         <func_t,max_refine,internal_fp_t>,
                         internal_fp_t> it;

    fermi_dirac_integ_direct() {
      // AWS 8/14/21: changed from 1.0e-17 to 1.0e-14 because it
      // appeard to more frequently converge (see polylog_ts) without
      // sacrificing accuracy.
      it.iiu.tol_rel=1.0e-14;

      // Set up the arguments in the internal precision
      three=3;
      two=2;
      half=1;
      half/=two;
      three_half=three;
      three_half/=two;
    }

    /** \brief Set tolerance
     */
    void set_tol(const internal_fp_t &tol) {
      it.iiu.tol_rel=tol;
      return;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 1/2 \f$
     */
    fp_t calc_1o2(fp_t y) {
      internal_fp_t y2=y, res, err;
      it.calc_err(half,y2,res,err);
      return ((fp_t)res);
    }

    /** \brief Fermi-Dirac integral of order \f$ 1/2 \f$
     */
    int calc_1o2_ret(fp_t y, fp_t &res, fp_t &err) {
      internal_fp_t y2=static_cast<internal_fp_t>(y), res2, err2;
      int iret=it.calc_err(half,y2,res2,err2);
      res=(fp_t)res2;
      err=(fp_t)err2;
      return iret;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ -1/2 \f$
     */
    fp_t calc_m1o2(fp_t y) {
      internal_fp_t y2=y, res, err;
      it.calc_err(-half,y2,res,err);
      return ((fp_t)res);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ -1/2 \f$
     */
    int calc_m1o2_ret(fp_t y, fp_t &res, fp_t &err) {
      internal_fp_t y2=static_cast<internal_fp_t>(y), res2, err2;
      int iret=it.calc_err(-half,y2,res2,err2);
      res=(fp_t)res2;
      err=(fp_t)err2;
      return iret;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3/2 \f$
     */
    fp_t calc_3o2(fp_t y) {
      internal_fp_t y2=y, res, err;
      it.calc_err(three_half,y2,res,err);
      return ((fp_t)res);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3/2 \f$
     */
    int calc_3o2_ret(fp_t y, fp_t &res, fp_t &err) {
      internal_fp_t y2=static_cast<internal_fp_t>(y), res2, err2;
      int iret=it.calc_err(three_half,y2,res2,err2);
      res=(fp_t)res2;
      err=(fp_t)err2;
      return iret;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 2 \f$
     */
    fp_t calc_2(fp_t y) {
      internal_fp_t y2=y, res, err;
      it.calc_err(two,y2,res,err);
      return ((fp_t)res);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 2 \f$
     */
    int calc_2_ret(fp_t y, fp_t &res, fp_t &err) {
      internal_fp_t y2=static_cast<internal_fp_t>(y), res2, err2;
      int iret=it.calc_err(two,y2,res2,err2);
      res=(fp_t)res2;
      err=(fp_t)err2;
      return iret;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    fp_t calc_3(fp_t y) {
      internal_fp_t y2=y, res, err;
      it.calc_err(three,y2,res,err);
      return ((fp_t)res);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    int calc_3_ret(fp_t y, fp_t &res, fp_t &err) {
      internal_fp_t y2=static_cast<internal_fp_t>(y), res2, err2;
      int iret=it.calc_err(three,y2,res2,err2);
      res=(fp_t)res2;
      err=(fp_t)err2;
      return iret;
    }
    
  };

  /** \brief Fermi-Dirac integral using multiprecison

      \verbatim embed:rst

      .. todo::
      
         - In class fermi_dirac_multip: implement degenerate
           and nondegenerate expansions.

      \endverbatim
   */
  class fermi_dirac_multip {
    
  protected:
    
    /// Tolerance
    double tol;

    /** \brief The Fermi-Dirac function
     */
    template<class fp_t, class fp2_t> fp_t obj_func(fp_t x, fp2_t a2,
                                                    fp2_t mu2) {
      fp_t res;
      fp_t a=static_cast<fp_t>(a2);
      fp_t mu=static_cast<fp_t>(mu2);
      if (x==0.0) res=0;
      else if (x-mu>std::numeric_limits<fp_t>::max_exponent) res=0;
      else res=pow(x,a)/(1+exp(x-mu));
      return res;
    }
    
  public:

    inte_kronrod_boost<61> ikb;
    inte_multip_double_exp_boost ideb;
    inte_adapt_cern iac;
    
    fermi_dirac_multip() {
      tol=1.0e-16;
      err_nonconv=true;
      ikb.err_nonconv=false;
      ideb.err_nonconv=false;
      iac.err_nonconv=false;

      // Set up the arguments in the internal precision
    }

    /** \brief If true, then convergence failures call the error 
        handler (default true)
    */
    bool err_nonconv;

    /** \brief Set tolerance
     */
    void set_tol(const double &tol_) {
      tol=tol_;
      return;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 1/2 \f$
     */
    template<class fp_t>
    int calc_1o2_ret_full(fp_t y, fp_t &res, fp_t &err, int &method) {

      fp_t two=2;
      fp_t one=1;
      fp_t half=one/two;
      fp_t zero=0;
      
      int ret1=ikb.integ_iu_err_multip([this,half,y](auto &&x) mutable 
      { return this->obj_func(x,half,y); },zero,res,err,tol);
      if (ret1==0) {
        method=1;
        return 0;
      }

      int ret2=ideb.integ_iu_err_multip([this,half,y](auto &&x) mutable 
      { return this->obj_func(x,half,y); },zero,res,err,tol);
      if (ret2==0) {
        method=2;
        return 0;
      }
      
      int ret3=iac.integ_iu_err_multip([this,half,y](auto &&x) mutable 
      { return this->obj_func(x,half,y); },zero,res,err,tol);
      if (ret3==0) {
        method=3;
        return 0;
      }
      
      O2SCL_CONV_RET("Function calc_1o2 failed.",o2scl::exc_efailed,
                     err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Fermi-Dirac integral of order \f$ 1/2 \f$
     */
    template<class fp_t>
    int calc_1o2_ret(fp_t y, fp_t &res, fp_t &err) {
      int method;
      return calc_1o2_ret_full(y,res,err,method);
      return 0;
    }

    /** \brief Fermi-Dirac integral of order \f$ 1/2 \f$
     */
    template<class fp_t>
    fp_t calc_1o2(fp_t y) {
      fp_t res, err;
      calc_1o2_ret(y,res,err);
      return res;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ -1/2 \f$
     */
    template<class fp_t>
    int calc_m1o2_ret_full(fp_t y, fp_t &res, fp_t &err, int &method) {

      fp_t two=2;
      fp_t one=1;
      fp_t mhalf=-one/two;
      fp_t zero=0;

      int ret1=ikb.integ_iu_err_multip([this,mhalf,y](auto &&x) mutable 
      { return this->obj_func(x,mhalf,y); },zero,res,err,tol);
      if (ret1==0) {
        method=1;
        return 0;
      }

      int ret2=ideb.integ_iu_err_multip([this,mhalf,y](auto &&x) mutable 
      { return this->obj_func(x,mhalf,y); },zero,res,err,tol);
      if (ret2==0) {
        method=2;
        return 0;
      }
      
      int ret3=iac.integ_iu_err_multip([this,mhalf,y](auto &&x) mutable 
      { return this->obj_func(x,mhalf,y); },zero,res,err,tol);
      if (ret3==0) {
        method=3;
        return 0;
      }
      
      O2SCL_CONV_RET("Function calc_m1o2_ret_full failed.",
                     o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Fermi-Dirac integral of order \f$ -1/2 \f$
     */
    template<class fp_t>
    int calc_m1o2_ret(fp_t y, fp_t &res, fp_t &err) {
      int method;
      return calc_m1o2_ret_full(y,res,err,method);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ -1/2 \f$
     */
    template<class fp_t>
    fp_t calc_m1o2(fp_t y) {
      fp_t res, err;
      calc_m1o2_ret(y,res,err);
      return res;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3/2 \f$
     */
    template<class fp_t>
    int calc_3o2_ret_full(fp_t y, fp_t &res, fp_t &err, int &method) {

      fp_t three=3;
      fp_t two=2;
      fp_t three_half=-three/two;
      fp_t zero=0;
      
      int ret1=ikb.integ_iu_err_multip([this,three_half,y](auto &&x) mutable 
      { return this->obj_func(x,three_half,y); },zero,res,err,tol);
      if (ret1==0) {
        method=1;
        return 0;
      }

      int ret2=ideb.integ_iu_err_multip([this,three_half,y](auto &&x) mutable 
      { return this->obj_func(x,three_half,y); },zero,res,err,tol);
      if (ret2==0) {
        method=2;
        return 0;
      }
      
      int ret3=iac.integ_iu_err_multip([this,three_half,y](auto &&x) mutable 
      { return this->obj_func(x,three_half,y); },zero,res,err,tol);
      if (ret3==0) {
        method=3;
        return 0;
      }
      
      O2SCL_CONV_RET("Function calc_3o2_ret_full failed.",
                     o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Fermi-Dirac integral of order \f$ 3/2 \f$
     */
    template<class fp_t>
    int calc_3o2_ret(fp_t y, fp_t &res, fp_t &err) {
      int method;
      return calc_3o2_ret_full(y,res,err,method);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3/2 \f$
     */
    template<class fp_t>
    fp_t calc_3o2(fp_t y) {
      fp_t res, err;
      calc_3o2_ret(y,res,err);
      return res;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 2 \f$
     */
    template<class fp_t>
    int calc_2_ret_full(fp_t y, fp_t &res, fp_t &err, int &method) {

      fp_t two=2;
      fp_t zero=0;
      
      int ret1=ikb.integ_iu_err_multip([this,two,y](auto &&x) mutable 
      { return this->obj_func(x,two,y); },zero,res,err,tol);
      if (ret1==0) {
        method=1;
        return 0;
      }

      int ret2=ideb.integ_iu_err_multip([this,two,y](auto &&x) mutable 
      { return this->obj_func(x,two,y); },zero,res,err,tol);
      if (ret2==0) {
        method=2;
        return 0;
      }
      
      int ret3=iac.integ_iu_err_multip([this,two,y](auto &&x) mutable 
      { return this->obj_func(x,two,y); },zero,res,err,tol);
      if (ret3==0) {
        method=3;
        return 0;
      }
      
      O2SCL_CONV_RET("Function calc_2_ret_full failed.",
                     o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Fermi-Dirac integral of order \f$ 2 \f$
     */
    template<class fp_t>
    int calc_2_ret(fp_t y, fp_t &res, fp_t &err) {
      int method;
      return calc_2_ret_full(y,res,err,method);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 2 \f$
     */
    template<class fp_t>
    fp_t calc_2(fp_t y) {
      fp_t res, err;
      calc_2_ret(y,res,err);
      return res;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    template<class fp_t>
    int calc_3_ret_full(fp_t y, fp_t &res, fp_t &err, int &method) {

      fp_t three=3;
      fp_t zero=0;
      
      int ret1=ikb.integ_iu_err_multip([this,three,y](auto &&x) mutable 
      { return this->obj_func(x,three,y); },zero,res,err,tol);
      if (ret1==0) {
        method=1;
        return 0;
      }

      int ret2=ideb.integ_iu_err_multip([this,three,y](auto &&x) mutable 
      { return this->obj_func(x,three,y); },zero,res,err,tol);
      if (ret2==0) {
        method=2;
        return 0;
      }
      
      int ret3=iac.integ_iu_err_multip([this,three,y](auto &&x) mutable 
      { return this->obj_func(x,three,y); },zero,res,err,tol);
      if (ret3==0) {
        method=3;
        return 0;
      }
      
      O2SCL_CONV_RET("Function calc_3_ret_full failed.",
                     o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    template<class fp_t>
    int calc_3_ret(fp_t y, fp_t &res, fp_t &err) {
      int method;
      return calc_3_ret_full(y,res,err,method);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    template<class fp_t>
    fp_t calc_3(fp_t y) {
      fp_t res, err;
      calc_3_ret(y,res,err);
      return res;
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    template<class fp_t>
    int calc_err_full(fp_t n, fp_t y, fp_t &res, fp_t &err, int &method) {

      fp_t zero=0;
      
      int ret1=ikb.integ_iu_err_multip([this,n,y](auto &&x) mutable 
      { return this->obj_func(x,n,y); },zero,res,err,tol);
      if (ret1==0) {
        method=1;
        return 0;
      }

      int ret2=ideb.integ_iu_err_multip([this,n,y](auto &&x) mutable 
      { return this->obj_func(x,n,y); },zero,res,err,tol);
      if (ret2==0) {
        method=2;
        return 0;
      }
      
      int ret3=iac.integ_iu_err_multip([this,n,y](auto &&x) mutable 
      { return this->obj_func(x,n,y); },zero,res,err,tol);
      if (ret3==0) {
        method=3;
        return 0;
      }
      
      O2SCL_CONV_RET("Function calc_err_full failed.",
                     o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    template<class fp_t>
    int calc_err(fp_t n, fp_t y, fp_t &res, fp_t &err) {
      int method;
      return calc_err_full(n,y,res,err,method);
    }
    
    /** \brief Fermi-Dirac integral of order \f$ 3 \f$
     */
    template<class fp_t>
    fp_t calc(fp_t n, fp_t y) {
      fp_t res, err;
      calc_err(n,y,res,err);
      return res;
    }
    
  };

  /** \brief Desc
   */
  class bose_einstein_multip {
    
  protected:
    
    /// Tolerance
    double tol;

    /** \brief The Fermi-Dirac function
     */
    template<class fp_t, class fp2_t> fp_t obj_func(fp_t x, fp2_t a2,
                                                    fp2_t mu2) {
      fp_t res;
      fp_t a=static_cast<fp_t>(a2);
      fp_t mu=static_cast<fp_t>(mu2);
      if (x==0.0) res=0;
      else if (x-mu>std::numeric_limits<fp_t>::max_exponent) res=0;
      else res=pow(x,a)/(exp(x-mu)-1);
      return res;
    }
    
  public:

    inte_kronrod_boost<61> ikb;
    inte_multip_double_exp_boost ideb;
    inte_adapt_cern iac;
    
    bose_einstein_multip() {
      tol=1.0e-16;
      err_nonconv=true;
      ikb.err_nonconv=false;
      ideb.err_nonconv=false;
      iac.err_nonconv=false;

      // Set up the arguments in the internal precision
    }

    /** \brief If true, then convergence failures call the error 
        handler (default true)
    */
    bool err_nonconv;

    /** \brief Set tolerance
     */
    void set_tol(const double &tol_) {
      tol=tol_;
      return;
    }
    
    /** \brief Bose-Einstein integral
     */
    template<class fp_t>
    int calc_err_full(fp_t a, fp_t y, fp_t &res, fp_t &err, int &method) {

      fp_t zero=0;
      
      int ret1=ikb.integ_iu_err_multip([this,a,y](auto &&x) mutable 
      { return this->obj_func(x,a,y); },zero,res,err,tol);
      if (ret1==0) {
        method=1;
        return 0;
      }

      int ret2=ideb.integ_iu_err_multip([this,a,y](auto &&x) mutable 
      { return this->obj_func(x,a,y); },zero,res,err,tol);
      if (ret2==0) {
        method=2;
        return 0;
      }
      
      int ret3=iac.integ_iu_err_multip([this,a,y](auto &&x) mutable 
      { return this->obj_func(x,a,y); },zero,res,err,tol);
      if (ret3==0) {
        method=3;
        return 0;
      }
      
      O2SCL_CONV_RET("Function calc_err_full failed.",
                     o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Bose-Einstein integral
     */
    template<class fp_t>
    int calc_err(fp_t a, fp_t y, fp_t &res, fp_t &err) {
      int method;
      return calc_err_full(a,y,res,err,method);
    }
    
    /** \brief Bose-Einstein integral
     */
    template<class fp_t> fp_t calc(fp_t a, fp_t y) {
      fp_t res, err;
      calc_err(a,y,res,err);
      return res;
    }
    
  };

  /** \brief Compute exponentially scaled modified Bessel function of
      the second kind by direct integration

      This class computes \f$ K_n(z) e^z\f$ for \f$ n=1,2,3 \f$
      by directly integrating. It integrates the representation
      \f[
      K_n(z) e^{z} = \frac{\sqrt{\pi} z^{n}}{2^{n} \Gamma(n+1/2)}
      \int_1^{\infty} e^{z(1-t)} 
      \left(t^2-1\right)^{n-1/2}~dt
      \f]
      (see
      http://functions.wolfram.com/Bessel-TypeFunctions/BesselK/07/01/01/)
      by applying an integrator (of
      type \ref o2scl::bessel_K_exp_integ_tl) with a larger floating
      point type and then casting the result back to \c fp t. This
      should work with boost multiprecision types but is only
      currently tested with <tt>internal_fp_t=long double</tt>.

      With the default types, this class should give almost identical
      results to \ref o2scl::bessel_K_exp_integ_gsl .
  */
  template <class fp_t=double, class func_t=funct_ld, size_t max_refine=15,
	    class internal_fp_t=long double>
  class bessel_K_exp_integ_direct {
    
  public:

    /** \brief The integrator
     */
    bessel_K_exp_integ_tl<o2scl::inte_exp_sinh_boost
                          <func_t,max_refine,internal_fp_t>,internal_fp_t> it;
    
    bessel_K_exp_integ_direct() {
      it.iiu.tol_rel=1.0e-12;
    }

    /** \brief Set tolerance
     */
    void set_tol(const fp_t &tol) {
      it.iiu.tol_rel=tol;
      return;
    }
    
    /** \brief Compute \f$ K_1(x) e^x \f$
     */
    fp_t K1exp(fp_t x) {
      internal_fp_t x2=x, res, err;
      it.calc_err(1,x2,res,err);
      return ((fp_t)res);
    }

    /** \brief Compute \f$ K_1(x) e^x \f$
     */
    int K1exp_ret(fp_t x, fp_t &res, fp_t &err) {
      internal_fp_t x2=x, res2, err2;
      int iret=it.calc_err(1,x2,res2,err2);
      res=(fp_t)res2;
      err=(fp_t)err2;
      return iret;
    }
    
    /** \brief Compute \f$ K_2(x) e^x \f$
     */
    fp_t K2exp(fp_t x) {
      internal_fp_t x2=x, res, err;
      it.calc_err(2,x2,res,err);
      return ((fp_t)res);
    }
    
    /** \brief Compute \f$ K_2(x) e^x \f$
     */
    int K2exp_ret(fp_t x, fp_t &res, fp_t &err) {
      internal_fp_t x2=x, res2, err2;
      int iret=it.calc_err(2,x2,res2,err2);
      res=(fp_t)res2;
      err=(fp_t)err2;
      return iret;
    }
    
    /** \brief Compute \f$ K_3(x) e^x \f$
     */
    fp_t K3exp(fp_t x) {
      internal_fp_t x2=x, res, err;
      it.calc_err(3,x2,res,err);
      return ((fp_t)res);
    }
    
    /** \brief Compute \f$ K_3(x) e^x \f$
     */
    int K3exp_ret(fp_t x, fp_t &res, fp_t &err) {
      internal_fp_t x2=x, res2, err2;
      int iret=it.calc_err(3,x2,res2,err2);
      res=(fp_t)res2;
      err=(fp_t)err2;
      return iret;
    }
    
  };

  /** \brief Bessel K times exponential by brute force
   */
  template<class fp_t, size_t max1, size_t max2, size_t max3,
           class fp1_t, class fp2_t, class fp3_t>
  class bessel_K_exp_integ_bf {
    
  protected:
    
    /// Tolerance
    fp_t tol;

  public:
    
    /// Lowest precision integrator
    bessel_K_exp_integ_direct<fp_t,std::function<fp1_t(fp1_t)>,max1,
                              fp1_t> bke1;
    
    /// Medium precision integrator
    bessel_K_exp_integ_direct<fp_t,std::function<fp2_t(fp2_t)>,max2,
                              fp2_t> bke2;
    
    /// Highest precision integrator
    bessel_K_exp_integ_direct<fp_t,std::function<fp3_t(fp3_t)>,max3,
                              fp3_t> bke3;

    bessel_K_exp_integ_bf() {
      tol=1.0e-17;
      bke1.it.iiu.err_nonconv=false;
      bke2.it.iiu.err_nonconv=false;
      bke3.it.iiu.err_nonconv=false;
      err_nonconv=true;
    }

    /** \brief If true, then convergene failures call the error 
        handler (default true)
    */
    bool err_nonconv;

    /** \brief Set tolerance
     */
    void set_tol(const fp_t &tol_) {
      tol=tol_;
      return;
    }
    
    /** \brief Compute \f$ K_1(x) e^x \f$
     */
    int K1exp_ret_full(fp_t x, fp_t &res, fp_t &err, int &method) {
      bke1.set_tol(tol);
      int ret1=bke1.K1exp_ret(x,res,err);
      if (ret1==0) {
        method=1;
        return 0;
      }
      bke2.set_tol(tol);
      int ret2=bke2.K1exp_ret(x,res,err);
      if (ret2==0) {
        method=2;
        return 0;
      }
      bke3.set_tol(tol);
      int ret3=bke3.K1exp_ret(x,res,err);
      if (ret3==0) {
        method=3;
      } else {
        method=0;
      }
      return ret3;
    }
    
    /** \brief Compute \f$ K_1(x) e^x \f$
     */
    int K1exp_ret(fp_t x, fp_t &res, fp_t &err) {
      int method;
      int iret=K1exp_ret_full(2,res,err,method);
      if (iret!=0) {
        O2SCL_CONV_RET("Function K1exp failed.",o2scl::exc_efailed,
                       err_nonconv);
      }
      return iret;
    }
    
    /** \brief Compute \f$ K_1(x) e^x \f$
     */
    fp_t K1exp(fp_t x) {
      fp_t res, err;
      K1exp_ret(x,res,err);
      return res;
    }

    /** \brief Compute \f$ K_2(x) e^x \f$
     */
    int K2exp_ret_full(fp_t x, fp_t &res, fp_t &err, int &method) {
      bke1.set_tol(tol);
      int ret1=bke1.K2exp_ret(x,res,err);
      if (ret1==0) {
        method=1;
        return 0;
      }
      bke2.set_tol(tol);
      int ret2=bke2.K2exp_ret(x,res,err);
      if (ret2==0) {
        method=2;
        return 0;
      }
      bke3.set_tol(tol);
      int ret3=bke3.K2exp_ret(x,res,err);
      if (ret3==0) {
        method=3;
      } else {
        method=0;
      }
      return ret3;
    }
    
    /** \brief Compute \f$ K_2(x) e^x \f$
     */
    int K2exp_ret(fp_t x, fp_t &res, fp_t &err) {
      int method;
      int iret=K2exp_ret_full(2,res,err,method);
      if (iret!=0) {
        O2SCL_CONV_RET("Function K2exp failed.",o2scl::exc_efailed,
                       err_nonconv);
      }
      return iret;
    }
    
    /** \brief Compute \f$ K_2(x) e^x \f$
     */
    fp_t K2exp(fp_t x) {
      fp_t res, err;
      K2exp_ret(x,res,err);
      return res;
    }

    /** \brief Compute \f$ K_3(x) e^x \f$
     */
    int K3exp_ret_full(fp_t x, fp_t &res, fp_t &err, int &method) {
      bke1.set_tol(tol);
      int ret1=bke1.K3exp_ret(x,res,err);
      if (ret1==0) {
        method=1;
        return 0;
      }
      bke2.set_tol(tol);
      int ret2=bke2.K3exp_ret(x,res,err);
      if (ret2==0) {
        method=2;
        return 0;
      }
      bke3.set_tol(tol);
      int ret3=bke3.K3exp_ret(x,res,err);
      if (ret3==0) {
        method=3;
      } else {
        method=0;
      }
      return ret3;
    }
    
    /** \brief Compute \f$ K_3(x) e^x \f$
     */
    int K3exp_ret(fp_t x, fp_t &res, fp_t &err) {
      int method;
      int iret=K3exp_ret_full(2,res,err,method);
      if (iret!=0) {
        O2SCL_CONV_RET("Function K3exp failed.",o2scl::exc_efailed,
                       err_nonconv);
      }
      return iret;
    }
    
    /** \brief Compute \f$ K_3(x) e^x \f$
     */
    fp_t K3exp(fp_t x) {
      fp_t res, err;
      K3exp_ret(x,res,err);
      return res;
    }
    
    };
  
  /** \brief Class to compute the polylogarithm function

      \note experimental

      This class uses integral representations based on the
      Fermi-Dirac or Bose-Einstein functions to compute the polylog
      functions.

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
      this representation works for negative values of \f$ z \f$.

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
      this representation works for positive values of \f$ z \f$.

      \verbatim embed:rst
      .. todo::
      
         In class polylog, test with higher accuracy floating point
         types.

      A classic reference for the polylogarithm function is [Lewin81]_.
      \endverbatim
  */
  template<class fp_t=double, class func_t=funct_ld, size_t max_refine=15,
           class internal_fp_t=long double> class polylog {

  protected:
    
    /** \brief The integrator for negative arguments
     */
    fermi_dirac_integ_tl<o2scl::inte_exp_sinh_boost
                         <func_t,max_refine,internal_fp_t>,
                         internal_fp_t> it_fd;
    
    /** \brief The integrator for positive arguments
     */
    bose_einstein_integ_tl<o2scl::inte_exp_sinh_boost
      <func_t,max_refine,internal_fp_t>,internal_fp_t> it_be;
    
  public:
    
    polylog() {
      it_fd.iiu.tol_rel=1.0e-17;
      it_be.iiu.tol_rel=1.0e-17;
    }
    
    /** \brief Set tolerance
     */
    void set_tol(fp_t tol) {
      it_fd.iiu.tol_rel=tol;
      it_be.iiu.tol_rel=tol;
      return;
    }
    
    /** \brief Polylogarithm function
     */
    fp_t calc(fp_t s, fp_t y) {
      if (y<0.0) {
	// Fermi-Dirac integral representation
	internal_fp_t a=s-1, mu=log(-y), res, err;
	it_fd.calc_err(a,mu,res,err);
	return -((fp_t)res/boost::math::tgamma(s));
      } else if (y==0.0) {
	return 0.0;
      }
      // Bose-Einstein integral representation
      internal_fp_t a=s-1, mu=log(y), res, err;
      it_be.calc_err(a,mu,res,err);
      return ((fp_t)res/boost::math::tgamma(s));
    }

  };

  /** \brief Polylogarithms with multiprecision integration
   */
  template<class fp_t, class internal_fp_t> class polylog_multip {

  protected:

    /** \brief The integrator for negative arguments
     */
    fermi_dirac_multip fdm;
    
    /** \brief The integrator for positive arguments
     */
    bose_einstein_multip bem;
    
  public:

    polylog_multip() {
      fdm.set_tol(std::numeric_limits<fp_t>::epsilon()/2.0);
      bem.set_tol(std::numeric_limits<fp_t>::epsilon()/2.0);
    }
    
    /** \brief Set tolerance
     */
    void set_tol(fp_t tol) {
      fdm.set_tol(tol);
      bem.set_tol(tol);
      return;
    }
    
    /** \brief Polylogarithm function
     */
    fp_t calc(fp_t s, fp_t y) {
      if (y<0.0) {
	// Fermi-Dirac integral representation
	internal_fp_t a=s-1, mu=log(-y), res, err;
	fdm.calc_err(a,mu,res,err);
	return -((fp_t)res/boost::math::tgamma(s));
      } else if (y==0.0) {
	return 0.0;
      }
      // Bose-Einstein integral representation
      internal_fp_t a=s-1, mu=log(y), res, err;
      bem.calc_err(a,mu,res,err);
      return ((fp_t)res/boost::math::tgamma(s));
    }

  };
  
}

#endif
