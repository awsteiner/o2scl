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

namespace o2scl {

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
    
  };

  /** \brief Gauss-Kronrod multiprecision integration class (Boost)

      \note The uncertainties reported by this class depend on those
      returned by the boost integration functions and are occasionally
      be underestimated. 
      
   */
  template<size_t rule=15> class inte_multip_kronrod_boost :
    public inte<funct,double>{
    
  protected:
    
    /// Maximum depth
    size_t max_depth;
    
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template <class fp_t>
    int integ_err_funct(funct &f, fp_t a, fp_t b, 
                       fp_t &res, fp_t &err, fp_t &L1norm_loc,
                       double target_tol, double integ_tol) {
      int ret=0;
      
      res=boost::math::quadrature::gauss_kronrod<fp_t,rule>::integrate
        (f,a,b,max_depth,target_tol,&err,&L1norm_loc);
      this->L1norm=static_cast<double>(L1norm_loc);
      
      if (verbose>1) {
        std::cout << "inte_multip_kronrod_boost::integ_err() "
                  << "tols(target,integ),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << err << std::endl;
      }

      if (err>integ_tol) {
        return 1;
      }
      
      return ret;
    }
    
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err

        There are three tolerances:
        - \c target_tol is the target tolerance which is sent to
        the boost integration function. The error value returned 
        by the boost integration function is often larger than this
        - \c integ_tol is the desired final tolerance of the integration.
        This function regards the integration as a failure if the 
        error value is larger than \c integ_tol
        - \c func_tol is the tolerance for evaluations of the 
        integrand. This value is passed to \ref o2scl::funct_multip.

    */
    template <typename func_t, class fp_t>
    int integ_err_int(func_t &&func, fp_t a, fp_t b, 
                      fp_t &res, fp_t &err, fp_t &L1norm_loc,
                      double target_tol, double integ_tol,
                      double func_tol) {
      
      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;

      std::function<fp_t(fp_t)> fx=[fm2,func](fp_t x) mutable -> fp_t
      { return fm2(func,x); };
      
      integ_err_funct(fx,a,b,res,err,L1norm_loc,target_tol,
                      integ_tol);

      return 0;
    }

    /// \name Typedefs for multiprecision
    //@{
    typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<50>> cpp_dec_float_50;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;
    //@}

  public:

    /** \brief Set the maximum number of interval splittings
     */
    void set_max_depth(size_t md) {
      max_depth=md;
      return;
    }
    
    /// Number of refinement levels in last integral computed
    size_t levels;
    
    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief Relative tolerance
     */
    double tol_abs;

    /** \brief Relative tolerance for multiprecision integrations
     */
    double tol_rel_multip;

    /** \brief Power for tolerance of function evaluations 
        (default 1.33)
     */
    double pow_tol_func;

    /// L1 norm
    double L1norm;
    
    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief If true, call the error handler if the integration
        does not succeed (default true)
    */
    bool err_nonconv;
    
    inte_multip_kronrod_boost() {
      tol_rel=-1.0;
      verbose=0;
      pow_tol_func=1.33;
      max_depth=15;
      err_nonconv=true;
      tol_rel=1.0e-8;
      tol_abs=1.0e-8;
    }

    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<class fp_t>
    int integ_err(funct &func, fp_t a, fp_t b, 
                  fp_t &res, fp_t &err) {
      
      fp_t L1norm_loc;
      int ret=integ_err_int2(func,a,b,res,err,L1norm_loc,
                             this->tol_rel,this->tol_rel/10.0);
      
      if (ret!=0) {
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
      
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    virtual int integ_err(funct &func, double a, double b, 
                          double &res, double &err) {
      return integ_err<double>(func,a,b,res,err);
    }
    
    /** \brief Calculate the first derivative of \c func  w.r.t. x and 
	uncertainty
    */
    template <typename func_t, class fp_t>
    int integ_err_multip(func_t &&func, fp_t a, fp_t b, 
                         fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (integ_tol<=0.0) {
        if (tol_rel_multip<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel_multip;
        }
      } 

      if (verbose>0) {
        std::cout << "int_multip_kronrod_boost::integ_err(): set "
                  << "tolerance to: " << integ_tol << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      double func_tol=pow(integ_tol,pow_tol_func);

      double target_tol=integ_tol/10.0;
      
      int ret;

      if (integ_tol>pow(10.0,-std::numeric_limits<double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_kronrod_boost::integ_err(): "
            << integ_tol << " > "
            << pow(10.0,-std::numeric_limits<double>::digits10+3)
            << "\n  for double integration." << std::endl;
        }
        double a_d=static_cast<double>(a);
        double b_d=static_cast<double>(b);
        double res_d, err_d, L1norm_d;
        
        ret=integ_err_int(func,a_d,b_d,res_d,err_d,L1norm_d,
                          target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_kronrod_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double a_ld=static_cast<long double>(a);
        long double b_ld=static_cast<long double>(b);
        long double res_ld, err_ld, L1norm_ld;
        
        ret=integ_err_int(func,a_ld,b_ld,res_ld,err_ld,L1norm_ld,
                          target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_ld<integ_tol) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_25>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_kronrod_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 a_cdf25=static_cast<cpp_dec_float_25>(a);
        cpp_dec_float_25 b_cdf25=static_cast<cpp_dec_float_25>(b);
        cpp_dec_float_25 res_cdf25, err_cdf25, L1norm_cdf25;
        
        ret=integ_err_int(func,a_cdf25,b_cdf25,res_cdf25,
                          err_cdf25,L1norm_cdf25,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf25<integ_tol) {
          res=static_cast<fp_t>(res_cdf25);
          err=static_cast<fp_t>(err_cdf25);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_35>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_kronrod_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 a_cdf35=static_cast<cpp_dec_float_35>(a);
        cpp_dec_float_35 b_cdf35=static_cast<cpp_dec_float_35>(b);
        cpp_dec_float_35 res_cdf35, err_cdf35, L1norm_cdf35;
        
        ret=integ_err_int(func,a_cdf35,b_cdf35,res_cdf35,
                          err_cdf35,L1norm_cdf35,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf35<integ_tol) {
          res=static_cast<fp_t>(res_cdf35);
          err=static_cast<fp_t>(err_cdf35);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_50>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_kronrod_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 a_cdf50=static_cast<cpp_dec_float_50>(a);
        cpp_dec_float_50 b_cdf50=static_cast<cpp_dec_float_50>(b);
        cpp_dec_float_50 res_cdf50, err_cdf50, L1norm_cdf50;
        
        ret=integ_err_int(func,a_cdf50,b_cdf50,res_cdf50,
                          err_cdf50,L1norm_cdf50,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf50<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_100>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_kronrod_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_100>::digits10+3)
                    << "\n  for cpp_dec_float_100 integration." << std::endl;
        }
        cpp_dec_float_100 a_cdf100=static_cast<cpp_dec_float_100>(a);
        cpp_dec_float_100 b_cdf100=static_cast<cpp_dec_float_100>(b);
        cpp_dec_float_100 res_cdf100, err_cdf100, L1norm_cdf100;
        
        ret=integ_err_int(func,a_cdf100,b_cdf100,res_cdf100,
                          err_cdf100,L1norm_cdf100,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf100<integ_tol) {
          res=static_cast<fp_t>(res_cdf100);
          err=static_cast<fp_t>(err_cdf100);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_multip_kronrod_boost::integ_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_multip_kronrod_boost::integ_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

  };
  
}

#endif
