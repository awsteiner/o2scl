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
//#include <gsl/gsl_sf_dilog.h>
#include <o2scl/inte_adapt_cern.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Polylogarithms (approximate) \f$ Li_n(x)\f$ 
    
      This class is experimental.

      This gives an approximation to the polylogarithm functions.

      Only works at present for \f$n=0,1,...,6\f$. Uses GSL library
      for n=2.

      Uses linear interpolation for \f$-1<x<0\f$
      and a series expansion for \f$x<-1\f$

      \future
      - Give error estimate? 
      - Improve accuracy?
      - Use more sophisticated interpolation?
      - Add the series \f$Li(n,x)=x+2^{-n} x^2+3^{-n} x^3+...\f$ 
      for \f$ x \rightarrow 0\f$?
      - Implement for positive arguments < 1.0
      - Make another polylog class which implements series acceleration?

      For reference, there are exact relations
      \f[
      \mathrm{Li}_2 \left(\frac{1}{2}\right) =
      \frac{1}{12}\left[\pi^2-6\left(\ln 2\right)^2\right]
      \f]
      \f[
      \mathrm{Li}_3 \left(\frac{1}{2}\right) =
      \frac{1}{24}\left[ 4\left(\ln 2\right)^3 - 2 \pi^2 \ln 2 +
      21 \zeta (3) \right]
      \f]
      \f[
      \mathrm{Li}_{-1} (x) = \frac{x}{\left(1-x\right)^2}
      \f]
      \f[
      \mathrm{Li}_{-2} (x) = \frac{x\left(x+1\right)}{\left(1-x\right)^3}
      \f]

  */
  class polylog {

  public:

    /// 0-th order polylogarithm = \f$ x/(1-x)\f$
    double li0(double x);

    /// 1-st order polylogarithm = \f$ -\ln(1-x) \f$
    double li1(double x);

    /// 2-nd order polylogarithm
    double li2(double x);

    /// 3-rd order polylogarithm
    double li3(double x);

    /// 4-th order polylogarithm
    double li4(double x);

    /// 5-th order polylogarithm
    double li5(double x);

    /// 6-th order polylogarithm
    double li6(double x);

    polylog();
    
    ~polylog();

  protected:

#ifndef DOXYGEN_NO_O2NS

    double *arg;
    
    double *two;
    
    double *three;
    
    double *four;
    
    double *five;
    
    double *six;
    
    double li2neg1;
    
    double li4neg1;
    
    double li6neg1;

#endif

  };

  /** \brief Fermi-Dirac integral

      This class performs direct computation of the
      Fermi-Dirac integral
      \f[
      F_{a}(\mu) = \int_0^{\infty} \frac{x^a}{1+\exp^{x-\mu}}
      \f]
      using \ref o2scl::inte_adapt_cern . 
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
    iiu.integ_err(f,0.0,0.0,res,err);
    return;
  }
  
  };
  
  /** \brief Double-precision version of \ref o2scl::fermi_dirac_integ_tl
   */
  typedef fermi_dirac_integ_tl
    <inte_iu<funct,inte_adapt_cern
    <funct,inte_gauss56_cern<funct,double,
    inte_gauss56_coeffs_double>,
    100,double>,double>,double> fermi_dirac_integ_double;

#if defined(O2SCL_LD_TYPES) || defined(DOXYGEN)
  
  /** \brief Long double version of \ref o2scl::fermi_dirac_integ_tl
   */
  typedef fermi_dirac_integ_tl
    <inte_iu<funct_ld,inte_adapt_cern
    <funct_ld,inte_gauss56_cern<funct_ld,long double,
    inte_gauss56_coeffs_long_double>,
    1000,long double>,long double>,long double>
    fermi_dirac_integ_long_double;
  
#endif

  /** \brief Compute the fermion integrals for a non-relativistic
      particle using the GSL functions
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
  
  /** \brief Compute the fermion integrals for a non-relativistic
      particle by directly integrating in long double precision
   */
  class fermion_nr_integ_direct {

  protected:
    
    /** \brief The integrator
     */
    fermi_dirac_integ_long_double it;
    
  public:

    fermion_nr_integ_direct() {
      // it.iiu is the inte_iu object
      // it.iiu.def_inte is the inte_adapt_cern object
      it.iiu.def_inte.tol_rel=1.0e-16;
      it.iiu.def_inte.tol_abs=1.0e-16;
      it.iiu.def_inte.verbose=1;
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

  };
  
#endif
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
