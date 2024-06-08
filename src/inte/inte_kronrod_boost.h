/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2019-2024, Andrew W. Steiner
  
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
  
  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_INTE_KRONROD_BOOST_H
#define O2SCL_INTE_KRONROD_BOOST_H

/** \file inte_kronrod_boost.h
    \brief File defining \ref o2scl::inte_kronrod_boost
*/

#include <cmath>

#define BOOST_DISABLE_ASSERTS
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <o2scl/inte.h>
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
#include <o2scl/funct_multip.h>
#endif

namespace o2scl {

  /** \brief Gauss-Kronrod integration class with multiprecision
      (Boost)

      If the default value of \ref tol_rel is used, then this class
      uses the square root of \c numeric_limits::epsilon for the
      relative tolerance. For double precision numbers, this tolerance
      is usually about \f$ 10^{-8} \f$. If the final uncertainty
      exceeds this value, then the error handler is called, unless
      \ref err_nonconv is false. Internally, the boost integration
      function is called with a tolerance which is a factor of 10
      smaller, because this is often necessary to ensure convergence.

      The multiprecision integration functions require a template
      function input, and their default tolerance is given by 
      \f$ 10^{-d} \f$ where \f$ d \f$ is \c numeric_limits::digits10 .

      \note The uncertainties reported by this class depend on those
      returned by the boost integration functions and are occasionally
      underestimated.

      \note The default maximum depth may be insufficient, especially
      for high-precision types or multiprecision integration, and can
      be changed with \ref set_max_depth().

      \warning For sufficiently difficult integrands, the 
      multiprecision 
      functions may take a very long time to complete.
      
  */
  template<size_t rule=15,
           class fp_25_t=o2fp_25, class fp_35_t=o2fp_35,
           class fp_50_t=o2fp_50, class fp_100_t=o2fp_100>
  class inte_kronrod_boost {
    
  protected:
    
    /// \name Internal methods and data [protected]
    //@{
    /// Maximum depth (default 15)
    size_t max_depth;
    
    /** \brief Internal integration wrapper of the boost function
        which stores the L1 norm and tests if the uncertainty is
        sufficiently small

        This function is used by both \ref integ_err() and \ref
        integ_err_int() .
    */
    template<typename func_t, class fp_t>
    int integ_err_funct(func_t &f, fp_t a, fp_t b, 
                        fp_t &res, fp_t &err, fp_t &L1norm_loc,
                        double target_tol, double integ_tol) {
      int ret=0;

      if (verbose>2) {
        std::cout << "  inte_kronrod_boost::integ_err_funct(): "
                  << "Calling boost::integrate." << std::endl;
      }
      
      // AWS, 9/25/23: I don't think that the boost kronrod
      // integrate() function can accept an rvalue reference (even
      // though the double_exp integration function can), so this
      // function must be a bit different than integ_err_int() below.
      res=boost::math::quadrature::gauss_kronrod<fp_t,rule>::integrate
        (f,a,b,max_depth,target_tol,&err,&L1norm_loc);
      this->L1norm=static_cast<double>(L1norm_loc);
      
      if (verbose>1) {
        std::cout << "  inte_kronrod_boost::integ_err_funct() "
                  << "tols(target,integ),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << err << std::endl;
      }
      
      if (err/abs(res)>integ_tol) {
        if (verbose>0) {
          std::cout << "  inte_kronrod_boost::integ_err_funct() failed "
                    << "because " << err/abs(res) << " > "
                    << integ_tol << std::endl;
        }
        return 1;
      }
      
      return ret;
    }

    /** \brief Integrate function \c func, the internal wrapper
        which uses a \ref funct_multip object

        There are three tolerances:
        - \c target_tol is the target tolerance which is sent to
        the boost integration function. The error value returned 
        by the boost integration function is often larger than this
        - \c integ_tol is the desired final tolerance of the integration.
        This function regards the integration as a failure if the 
        error value is larger than \c integ_tol
        - \c func_tol is the tolerance for evaluations of the 
        integrand. This value is passed to \ref o2scl::funct_multip.

        This function is used by \ref integ_err_multip()
    */
    template <typename func_t, class fp_t>
    int integ_err_int(func_t &&func, fp_t a, fp_t b, 
                      fp_t &res, fp_t &err, fp_t &L1norm_loc,
                      double target_tol, double integ_tol, double func_tol) {
      
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      
      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;
      fm2.verbose=fm_verbose;

      std::function<fp_t(fp_t)> fx=[fm2,func](fp_t x) mutable -> fp_t
      { return fm2(func,x); };

      if (verbose>2) {
        std::cout << "inte_kronrod_boost::integ_err_int(): "
                  << "Going to integ_err_funct()." << std::endl;
      }
      
      integ_err_funct(fx,a,b,res,err,L1norm_loc,target_tol,
                      integ_tol);

      if (verbose>1) {
        std::cout << "inte_kronrod_boost::integ_err_int() "
                  << "tols(target,integ,func),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << std::endl;
      }

#else
      res=0;
      err=std::numeric_limits<fp_t>::infinity();
#endif
      
      if (err/abs(res)>integ_tol) {
        if (verbose>0) {
          std::cout << "  inte_kronrod_boost::integ_err_int() failed "
                    << "because " << err/abs(res) << " > "
                    << integ_tol << std::endl;
        }
        return 1;
      }
      return 0;
    }
    //@}
    
  public:

    /// \name Constructor
    //@{
    inte_kronrod_boost() {
      verbose=0;
      pow_tol_func=1.33;
      max_depth=15;
      err_nonconv=true;
      tol_rel=-1.0;
      tol_abs=-1.0;
      fm_verbose=0;
    }
    //@}
    
    /// \name Integration settings
    //@{
    /** \brief Power for tolerance of function evaluations in
        multiprecision integrations (default 1.33)
    */
    double pow_tol_func;

    /** \brief The maximum relative uncertainty 
	in the value of the integral (default \f$ -1 \f$)
    */
    double tol_rel;

    /** \brief The maximum absolute uncertainty 
	in the value of the integral (default \f$ -1 \f$)
        
        \note This value is unused by this integrator, but this
        is included for compatibility with the other integrators. 
    */
    double tol_abs;

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief If true, call the error handler if the integration
        does not succeed (default true)
    */
    bool err_nonconv;

    /** \brief Verbosity parameter for the internal
        \ref funct_multip object
    */
    int fm_verbose;
    
    /** \brief Set the maximum number of interval splittings
        (default 15)
    */
    void set_max_depth(size_t md) {
      max_depth=md;
      return;
    }
    //@}

    /// \name Integration output quantities
    //@{
    /** \brief \f$ L_1 \f$ norm from the last integration
     */
    double L1norm;
    //@}

    /// \name Main integration functions
    //@{
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_err(func_t &func, fp_t a, fp_t b, fp_t &res, fp_t &err) {

      double tol_rel_loc;
      if (tol_rel<=0.0) {
        tol_rel_loc=sqrt(static_cast<double>
                         (std::numeric_limits<fp_t>::epsilon()));
      } else {
        tol_rel_loc=tol_rel;
      }
      
      fp_t L1norm_loc;
      int ret=integ_err_funct(func,a,b,res,err,L1norm_loc,
                              tol_rel_loc/10.0,tol_rel_loc);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "  Values err,tol_rel,L1norm,max: "
                    << err << " " << tol_rel_loc << " "
                    << L1norm << " " << max_depth
                    << std::endl;
        }
        O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_kronrod_boost::integ_err().",o2scl::exc_efailed,
                        this->err_nonconv);
      }
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \c b.
     */
    template<typename func_t, class fp_t>
    fp_t integ(func_t &func, fp_t a, fp_t b) {
      fp_t res, interror;
      int ret=integ_err(func,a,b,res,interror);
      if (ret!=0) {
	O2SCL_ERR2("Integration failed in inte_kronrod_",
                   "boost::integ().",o2scl::exc_efailed);
      }
      return res;
    }
    //@}

    /// \name Multiprecision integration functions
    //@{
    /** \brief Integrate function \c func from \c a to \c b using
        multipreicsion, placing the result in \c res and the error in
        \c err
    */
    template <typename func_t, class fp_t>
    int integ_err_multip(func_t &&func, fp_t a, fp_t b, 
                         fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (integ_tol<=0.0) {
        if (tol_rel<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel;
        }
      } 

      if (verbose>0) {
        std::cout << "inte_kronrod_boost::integ_err_multip(): set "
                  << "tolerance to: " << integ_tol << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      double func_tol=pow(integ_tol,pow_tol_func);

      // We set the target tolerance an order of magnitude smaller
      // than the desired tolerance to make sure we achieve the
      // requested tolerance
      double target_tol=integ_tol/10.0;
      
      int ret;

      // We require that there are 3 more digits in the floating point
      // type than the required integration tolerance
      if (integ_tol>pow(10.0,-std::numeric_limits<double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "  " << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << " for double integration." << std::endl;
        }
        double a_d=static_cast<double>(a);
        double b_d=static_cast<double>(b);
        double res_d, err_d, L1norm_d;

        ret=integ_err_int(func,a_d,b_d,res_d,err_d,L1norm_d,
                          target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d/abs(res_d)<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          if (verbose>0) {
            if (ret!=0) {
              std::cout << "  Failed. Returned non-zero value."
                        << std::endl;
            } else {
              std::cout << "  Failed. Relative error "
                        << err_d/abs(res_d) << " >= " << integ_tol
                        << std::endl;
            }
          }
          // AWS 3/23/23 I'm not sure why this is necessary
          //target_tol/=10;
        }
      }

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      
      if (integ_tol>pow(10.0,
                        -std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "  " << integ_tol << " > "
                    << pow(10.0,
                           -std::numeric_limits<long double>::digits10+3)
                    << " for long double integration." << std::endl;
        }
        long double a_ld=static_cast<long double>(a);
        long double b_ld=static_cast<long double>(b);
        long double res_ld, err_ld, L1norm_ld;
        
        ret=integ_err_int(func,a_ld,b_ld,res_ld,err_ld,L1norm_ld,
                          target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_ld/abs(res_ld)<integ_tol) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        } else {
          if (verbose>0) {
            if (ret!=0) {
              std::cout << "  Failed. Returned non-zero value."
                        << std::endl;
            } else {
              std::cout << "  Failed. Relative error "
                        << err_ld/abs(res_ld) << " >= " << integ_tol
                        << std::endl;
            }
          }
          // AWS 3/23/23 I'm not sure why this is necessary
          //target_tol/=10;
        }
      }
      
      if (integ_tol>pow(10.0,-std::numeric_limits
                        <fp_25_t>::digits10+3)) {
        if (verbose>0) {
          std::cout << "  " << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <fp_25_t>::digits10+3)
                    << " for fp_25_t integration." << std::endl;
        }
        fp_25_t a_cdf25=static_cast<fp_25_t>(a);
        fp_25_t b_cdf25=static_cast<fp_25_t>(b);
        fp_25_t res_cdf25, err_cdf25, L1norm_cdf25;

        ret=integ_err_int(func,a_cdf25,b_cdf25,res_cdf25,
                          err_cdf25,L1norm_cdf25,target_tol,
                          integ_tol,func_tol);

        if (verbose>1) {
          std::cout << "ret,res,err,tol: " << ret << " "
                    << res_cdf25 << " " << err_cdf25 << " "
                    << integ_tol << std::endl;
        }
        if (ret==0 && err_cdf25/abs(res_cdf25)<integ_tol) {
          res=static_cast<fp_t>(res_cdf25);
          err=static_cast<fp_t>(err_cdf25);
          return 0;
        } else {
          if (verbose>0) {
            if (ret!=0) {
              std::cout << "  Failed. Returned non-zero value."
                        << std::endl;
            } else {
              std::cout << "  Failed. Relative error "
                        << err_cdf25/abs(res_cdf25) << " >= " << integ_tol
                        << std::endl;
            }
          }
          // AWS 3/23/23 I'm not sure why this is necessary
          //target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <fp_35_t>::digits10+3)) {
        if (verbose>0) {
          std::cout << "  " << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <fp_35_t>::digits10+3)
                    << " for fp_35_t integration." << std::endl;
        }
        fp_35_t a_cdf35=static_cast<fp_35_t>(a);
        fp_35_t b_cdf35=static_cast<fp_35_t>(b);
        fp_35_t res_cdf35, err_cdf35, L1norm_cdf35;
        
        ret=integ_err_int(func,a_cdf35,b_cdf35,res_cdf35,
                          err_cdf35,L1norm_cdf35,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf35/abs(res_cdf35)<integ_tol) {
          res=static_cast<fp_t>(res_cdf35);
          err=static_cast<fp_t>(err_cdf35);
          return 0;
        } else {
          if (verbose>0) {
            if (ret!=0) {
              std::cout << "  Failed. Returned non-zero value."
                        << std::endl;
            } else {
              std::cout << "  Failed. Relative error "
                        << err_cdf35/abs(res_cdf35) << " >= " << integ_tol
                        << std::endl;
            }
          }
          // AWS 3/23/23 I'm not sure why this is necessary
          //target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <fp_50_t>::digits10+3)) {
        if (verbose>0) {
          std::cout << "  " << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <fp_50_t>::digits10+3)
                    << " for fp_50_t integration." << std::endl;
        }
        fp_50_t a_cdf50=static_cast<fp_50_t>(a);
        fp_50_t b_cdf50=static_cast<fp_50_t>(b);
        fp_50_t res_cdf50, err_cdf50, L1norm_cdf50;
        
        ret=integ_err_int(func,a_cdf50,b_cdf50,res_cdf50,
                          err_cdf50,L1norm_cdf50,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf50/abs(res_cdf50)<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          if (verbose>0) {
            if (ret!=0) {
              std::cout << "  Failed. Returned non-zero value."
                        << std::endl;
            } else {
              std::cout << "  Failed. Relative error "
                        << err_cdf50/abs(res_cdf50) << " >= " << integ_tol
                        << std::endl;
            }
          }
          // AWS 3/23/23 I'm not sure why this is necessary
          //target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <fp_100_t>::digits10+3)) {
        if (verbose>0) {
          std::cout << "  " << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <fp_100_t>::digits10+3)
                    << " for fp_100_t integration." << std::endl;
        }
        fp_100_t a_cdf100=static_cast<fp_100_t>(a);
        fp_100_t b_cdf100=static_cast<fp_100_t>(b);
        fp_100_t res_cdf100, err_cdf100, L1norm_cdf100;
        
        ret=integ_err_int(func,a_cdf100,b_cdf100,res_cdf100,
                          err_cdf100,L1norm_cdf100,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf100/abs(res_cdf100)<integ_tol) {
          res=static_cast<fp_t>(res_cdf100);
          err=static_cast<fp_t>(err_cdf100);
          return 0;
        } else {
          if (verbose>0) {
            if (ret!=0) {
              std::cout << "  Failed. Returned non-zero value."
                        << std::endl;
            } else {
              std::cout << "  Failed. Relative error "
                        << err_cdf100/abs(res_cdf100) << " >= " << integ_tol
                        << std::endl;
            }
          }
          // AWS 3/23/23 I'm not sure why this is necessary
          //target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_kronrod_boost::integ_err_multip() "
                  << "failed after fp_100_t.\n  "
                  << "Original tolerance: " << integ_tol << std::endl;
      }

#endif
    
      O2SCL_CONV2_RET("Failed to compute with requested accuracy ",
                      "in inte_kronrod_boost::integ_err_multip().",
                      o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Integrate function \c func from \c a to \c b.
     */
    template<typename func_t, class fp_t>
    fp_t integ_multip(func_t &&func, fp_t a, fp_t b) {
      fp_t res, interror;
      int ret=integ_err_multip(func,a,b,res,interror);
      if (ret!=0) {
	O2SCL_ERR2("Integration failed in inte_kronrod_",
                   "boost::integ_multip().",o2scl::exc_efailed);
      }
      return res;
    }

    /** \brief Integrate function \c func from \c a to \f$ \infty \f$ using
        multipreicsion, placing the result in \c res and the error in
        \c err
    */
    template <typename func_t, class fp_t>
    int integ_iu_err_multip(func_t &&func, fp_t a, 
                            fp_t &res, fp_t &err, double integ_tol=-1.0) {
      return integ_err_multip(func,a,
                              std::numeric_limits<fp_t>::infinity(),
                              res,err,integ_tol);
    }
    
    /** \brief Integrate function \c func from \c a to \c b.
     */
    template<typename func_t, class fp_t>
    fp_t integ_iu_multip(func_t &&func, fp_t a, fp_t b) {
      fp_t res, interror;
      int ret=integ_iu_err_multip(func,a,b,res,interror);
      if (ret!=0) {
	O2SCL_ERR2("Integration failed in inte_kronrod_",
                   "boost::integ_iu_multip().",o2scl::exc_efailed);
      }
      return res;
    }

    /** \brief Integrate function \c func from \f$ -\infty \f$ to \c b using
        multipreicsion, placing the result in \c res and the error in
        \c err
    */
    template <typename func_t, class fp_t>
    int integ_il_err_multip(func_t &&func, fp_t b, 
                            fp_t &res, fp_t &err, double integ_tol=-1.0) {
      return integ_err_multip(func,
                              -std::numeric_limits<fp_t>::infinity(),
                              b,res,err,integ_tol);
    }
    
    /** \brief Integrate function \c func from \c a to \c b.
     */
    template<typename func_t, class fp_t>
    fp_t integ_il_multip(func_t &&func, fp_t a, fp_t b) {
      fp_t res, interror;
      int ret=integ_il_err_multip(func,a,b,res,interror);
      if (ret!=0) {
	O2SCL_ERR2("Integration failed in inte_kronrod_",
                   "boost::integ_il_multip().",o2scl::exc_efailed);
      }
      return res;
    }

    /** \brief Integrate function \c func from \f$ -\infty \f$ to \f$
        \infty \f$ using multipreicsion, placing the result in \c res
        and the error in \c err
    */
    template <typename func_t, class fp_t>
    int integ_i_err_multip(func_t &&func, 
                           fp_t &res, fp_t &err, double integ_tol=-1.0) {
      return integ_err_multip(func,
                              -std::numeric_limits<fp_t>::infinity(),
                              std::numeric_limits<fp_t>::infinity(),
                              res,err,integ_tol);
    }
    
    /** \brief Integrate function \c func from \c a to \c b.
     */
    template<typename func_t, class fp_t>
    fp_t integ_i_multip(func_t &&func, fp_t a, fp_t b) {
      fp_t res, interror;
      int ret=integ_i_err_multip(func,a,b,res,interror);
      if (ret!=0) {
	O2SCL_ERR2("Integration failed in inte_kronrod_",
                   "boost::integ_i_multip().",o2scl::exc_efailed);
      }
      return res;
    }
    //@}

  };
  
}

#endif
