/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_INTE_KRONROD_BOOST_H
#define O2SCL_INTE_KRONROD_BOOST_H

/** \file inte_kronrod_boost.h
    \brief File defining \ref o2scl::inte_kronrod_boost
*/

#include <cmath>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <o2scl/inte.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Gauss-Kronrod integration class (Boost)

      The rule parameter should be either 15, 31, 41, 51, or 61. 
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .

      \future Figure out what to do with L1norm. The boost
      documentation claims that "the error estimates provided by the
      routine are woefully pessimistic" and the integral appears to be
      correct, but the boost documentation also says "if there is a
      significant difference between this [the L1 norm] and the
      returned value, then the result is likely to be
      ill-conditioned". It would be nice to test L1 norm in some
      reasonable way.
  */
  template<class func_t=funct, size_t rule=15, class fp_t=double>
    class inte_kronrod_boost : public inte<func_t,fp_t> {
    
  protected:

  /// Maximum depth
  size_t max_depth;

  public:

  inte_kronrod_boost() {
    max_depth=15;
  }
  
  virtual ~inte_kronrod_boost() {
  }

  /** \brief Set the maximum number of interval splittings
   */
  void set_max_depth(size_t md) {
    max_depth=md;
    return;
  }
    
  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			fp_t &res, fp_t &err) {
    res=boost::math::quadrature::gauss_kronrod<fp_t,rule>::integrate
      (func,a,b,max_depth,this->tol_rel,&err,&L1norm);
    if (err>this->tol_rel) {
      if (this->verbose>0) {
        std::cout << "Function inte_kronrod_boost::integ_err() failed."
                  << std::endl;
        std::cout << "Values err,tol_rel,L1norm,max: "
                  << err << " " << this->tol_rel << " "
                  << L1norm << " " << max_depth
                  << std::endl;
      }
      O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                      "inte_kronrod_boost::integ_err().",o2scl::exc_efailed,
                      this->err_nonconv);
    }
    return 0;
  }
  
    /// L1 norm
    fp_t L1norm;
    
    /// Number of refinement levels in last integral computed
    size_t levels;
    
  };

  template<class func_t=funct_multip<>>
  class inte_multip_kronrod_boost {
    
  protected:
    
    typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<50>> cpp_dec_float_50;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

    /// \name The derivative objects for varying levels of precision
    //@{
    inte_kronrod_boost<func_t,15,double> ikb_d;
    inte_kronrod_boost<func_t,15,long double> ikb_ld;
    inte_kronrod_boost<func_t,15,cpp_dec_float_25> ikb_cdf25;
    inte_kronrod_boost<func_t,15,cpp_dec_float_35> ikb_cdf35;
    inte_kronrod_boost<func_t,15,cpp_dec_float_50> ikb_cdf50;
    inte_kronrod_boost<func_t,15,cpp_dec_float_100> ikb_cdf100;
    //@}
    
  public:

    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief Power for tolerance of function evaluations 
        (default 1.33)
     */
    double pow_tol_func;

    /** \brief Verbosity parameter
     */
    int verbose;

    inte_multip_kronrod_boost() {
      tol_rel=-1.0;
      verbose=0;
      pow_tol_func=1.33;
      ikb_d.err_nonconv=false;
      ikb_ld.err_nonconv=false;
      ikb_cdf25.err_nonconv=false;
      ikb_cdf35.err_nonconv=false;
      ikb_cdf50.err_nonconv=false;
      ikb_cdf100.err_nonconv=false;
    }

    /** \brief Calculate the first derivative of \c func  w.r.t. x and 
	uncertainty
    */
    template<class fp_t>
    int integ_err(func_t &func, fp_t a, fp_t b, 
                  fp_t &res, fp_t &err, double tol_loc=-1.0) {
      
      if (tol_loc<=0.0) {
        if (tol_rel<=0.0) {
          tol_loc=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          tol_loc=tol_rel;
        }
      } 

      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err(): set "
                  << "tolerance to: " << tol_loc << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      func.tol_rel=pow(tol_loc,pow_tol_func);
      
      int ret;
      
      if (tol_loc>pow(10.0,-std::numeric_limits<double>::digits10+3)) {
        double a_d=static_cast<double>(a);
        double b_d=static_cast<double>(b);
        double res_d, err_d;
        
        ikb_d.tol_rel=tol_loc;
        ret=ikb_d.integ_err(func,a_d,b_d,res_d,err_d);
        
        if (ret==0 && err_d<tol_loc) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        }
      }

      if (tol_loc>pow(10.0,-std::numeric_limits<long double>::digits10+3)) {
        long double a_ld=static_cast<long double>(a);
        long double b_ld=static_cast<long double>(b);
        long double res_ld, err_ld;
        
        ikb_ld.tol_rel=tol_loc;
        ret=ikb_ld.integ_err(func,a_ld,b_ld,res_ld,err_ld);
        
        if (ret==0 && err_ld<tol_loc) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        }
      }

      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_25>::digits10+3)) {
        cpp_dec_float_25 a_cdf25=static_cast<cpp_dec_float_25>(a);
        cpp_dec_float_25 b_cdf25=static_cast<cpp_dec_float_25>(b);
        cpp_dec_float_25 res_cdf25, err_cdf25;
        
        ikb_cdf25.tol_rel=tol_loc;
        ret=ikb_cdf25.integ_err(func,a_cdf25,b_cdf25,res_cdf25,err_cdf25);
        
        if (ret==0 && err_cdf25<tol_loc) {
          res=static_cast<fp_t>(res_cdf25);
          err=static_cast<fp_t>(err_cdf25);
          return 0;
        }
      }

      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_35>::digits10+3)) {
        cpp_dec_float_35 a_cdf35=static_cast<cpp_dec_float_35>(a);
        cpp_dec_float_35 b_cdf35=static_cast<cpp_dec_float_35>(b);
        cpp_dec_float_35 res_cdf35, err_cdf35;
        
        ikb_cdf35.tol_rel=tol_loc;
        ret=ikb_cdf35.integ_err(func,a_cdf35,b_cdf35,res_cdf35,err_cdf35);
        
        if (ret==0 && err_cdf35<tol_loc) {
          res=static_cast<fp_t>(res_cdf35);
          err=static_cast<fp_t>(err_cdf35);
          return 0;
        }
      }

      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_50>::digits10+3)) {
        cpp_dec_float_50 a_cdf50=static_cast<cpp_dec_float_50>(a);
        cpp_dec_float_50 b_cdf50=static_cast<cpp_dec_float_50>(b);
        cpp_dec_float_50 res_cdf50, err_cdf50;
        
        ikb_cdf50.tol_rel=tol_loc;
        ret=ikb_cdf50.integ_err(func,a_cdf50,b_cdf50,res_cdf50,err_cdf50);
        
        if (ret==0 && err_cdf50<tol_loc) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        }
      }

      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_100>::digits10+3)) {
        cpp_dec_float_100 a_cdf100=static_cast<cpp_dec_float_100>(a);
        cpp_dec_float_100 b_cdf100=static_cast<cpp_dec_float_100>(b);
        cpp_dec_float_100 res_cdf100, err_cdf100;
        
        ikb_cdf100.tol_rel=tol_loc;
        ret=ikb_cdf100.integ_err(func,a_cdf100,b_cdf100,res_cdf100,
                                 err_cdf100);
        
        if (ret==0 && err_cdf100<tol_loc) {
          res=static_cast<fp_t>(res_cdf100);
          err=static_cast<fp_t>(err_cdf100);
          return 0;
        }
      }

      if (verbose>0) {
        std::cout << "Function root_multip_brent_gsl::deriv_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << tol_loc << std::endl;
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in root_multip_brent_gsl::deriv_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

  };
  

  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
