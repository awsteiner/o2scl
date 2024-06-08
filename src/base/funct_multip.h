/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
#ifndef O2SCL_FUNCT_MULTIP_H
#define O2SCL_FUNCT_MULTIP_H

/** \file funct.h
    \brief Multiprecisions extension to Function object classes
*/

// For o2scl::dtos()
#include <o2scl/string_conv.h>

// for typeid()
#include <typeinfo>

#ifndef O2SCL_NO_BOOST_MULTIPRECISION

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <o2scl/set_mpfr.h>
#ifdef O2SCL_SET_MPFR
#include <boost/multiprecision/mpfr.hpp>
#endif

#endif

namespace o2scl {

#ifndef O2SCL_NO_BOOST_MULTIPRECISION

  /// \name Floating point typedefs in src/base/funct_multip.h
  //@{
#ifdef O2SCL_SET_MPFR  
  typedef boost::multiprecision::number<
    boost::multiprecision::mpfr_float_backend<25>> mpfr_25;
  typedef boost::multiprecision::number<
    boost::multiprecision::mpfr_float_backend<35>> mpfr_35;
  typedef boost::multiprecision::number<
    boost::multiprecision::mpfr_float_backend<50>> mpfr_50;
  typedef boost::multiprecision::number<
    boost::multiprecision::mpfr_float_backend<100>> mpfr_100;
#endif

  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<50>> cpp_dec_float_50;
  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

  // Choose the best floating point type depending on what is
  // available. 7/28/23: I'm currently having problems with mpfr types
  // so they're commented out, however, in the future, this should be
  // fixed as I think the MPFR types are faster.
  
  //#ifdef O2SCL_SET_MPFR
  //  typedef mpfr_25 o2fp_25;
  //typedef mpfr_35 o2fp_35;
  //typedef mpfr_50 o2fp_50;
  //typedef mpfr_100 o2fp_100;
  //#else
  typedef cpp_dec_float_25 o2fp_25;
  typedef cpp_dec_float_35 o2fp_35;
  typedef cpp_dec_float_50 o2fp_50;
  typedef cpp_dec_float_100 o2fp_100;
  //#endif
  //@}

#else

  typedef long double o2fp_25;
  typedef long double o2fp_35;
  typedef long double o2fp_50;
  typedef long double o2fp_100;
  
  // end of #ifndef O2SCL_NO_BOOST_MULTIPRECISION
#endif
  
  /// \name One-dimensional function typedefs in src/base/funct_multip.h
  //@{
  /** \brief One-dimensional long double function in 
      src/base/funct_multip.h
  */
  typedef std::function<long double(long double)> funct_ld;

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  
  /** \brief One-dimensional Boost 25-digit function in 
      src/base/funct_multip.h
  */
  typedef std::function<cpp_dec_float_25(cpp_dec_float_25)>
  funct_cdf25;
  
  /** \brief One-dimensional Boost 35-digit function in 
      src/base/funct_multip.h
  */
  typedef std::function<cpp_dec_float_35(cpp_dec_float_35)>
  funct_cdf35;
  
  /** \brief One-dimensional Boost 50-digit function in 
      src/base/funct_multip.h
  */
  typedef std::function<cpp_dec_float_50(cpp_dec_float_50)>
  funct_cdf50;
  
  /** \brief One-dimensional Boost 100-digit function in 
      src/base/funct_multip.h
  */
  typedef std::function<cpp_dec_float_100(cpp_dec_float_100)>
  funct_cdf100;

  /** \brief One-dimensional long double function with integer
      return in src/base/funct_multip.h
  */
  typedef std::function<int(long double,long double &)> funct_ret_ld;
  
  /** \brief One-dimensional Boost 25-digit function with integer
      return in src/base/funct_multip.h
  */
  typedef std::function<int(cpp_dec_float_25,cpp_dec_float_25 &)>
  funct_ret_cdf25;
  
  /** \brief One-dimensional Boost 35-digit function with integer
      return in src/base/funct_multip.h
  */
  typedef std::function<int(cpp_dec_float_35,cpp_dec_float_35 &)>
  funct_ret_cdf35;
  
  /** \brief One-dimensional Boost 50-digit function with integer
      return in src/base/funct_multip.h
  */
  typedef std::function<int(cpp_dec_float_50,cpp_dec_float_50 &)>
  funct_ret_cdf50;

  /** \brief One-dimensional Boost 100-digit function with integer
      return in src/base/funct_multip.h
  */
  typedef std::function<int(cpp_dec_float_100,cpp_dec_float_100 &)>
  funct_ret_cdf100;

  // end of #ifndef O2SCL_NO_BOOST_MULTIPRECISION
#endif
  
  //@}

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
#if defined (O2SCL_SET_MPFR) || defined (DOXYGEN)
  
  /** \brief One-dimensional function typedef in src/base/funct.h
      
      This typedef is defined only if --enable-mpfr is selected when
      when O2scl is configured.
  */
  typedef std::function<mpfr_25(mpfr_25)> funct_mpfr25;

  /** \brief One-dimensional function typedef in src/base/funct.h
      
      This typedef is defined only if --enable-mpfr is selected when
      when O2scl is configured.
  */
  typedef std::function<mpfr_35(mpfr_35)> funct_mpfr35;

  /** \brief One-dimensional function typedef in src/base/funct.h
      
      This typedef is defined only if --enable-mpfr is selected when
      when O2scl is configured.
  */
  typedef std::function<mpfr_50(mpfr_50)> funct_mpfr50;

  /** \brief One-dimensional function typedef in src/base/funct.h
      
      This typedef is defined only if --enable-mpfr is selected when
      when O2scl is configured.
  */
  typedef std::function<mpfr_100(mpfr_100)> funct_mpfr100;

  // end of #if defined (O2SCL_SET_MPFR) || defined (DOXYGEN)
#endif
  // end of #ifndef O2SCL_NO_BOOST_MULTIPRECISION
#endif
  
#ifndef O2SCL_NO_BOOST_MULTIPRECISION

  /** \brief Use multiprecision to automatically evaluate a function to
      a specified level of precision

      \note Experimental.

      The function must be specified as a template, i.e. it must be of
      the form <tt>template<class fp_t> fp_t function(fp_t x)</tt>.
      
      By default, this class uses a relative tolerance equal to \f$
      10^{-d+1} \f$, where \f$ d \f$ is the value returned by
      <tt>numeric_limits<fp_t>::digits10</tt>. This choice ensures
      that most exact or nearly exact functions will require only two
      function evaluations, one at \c double and one at \c long \c
      double precision. However, only simple functions can be
      evaluated to within this accuracy without more precise inputs.
      Preferably, the user should choose this tolerance carefully. If
      the tolerance is sufficiently small, no low-precision function
      evaluations will be performed.

      This class will fail to evalate a function with the requested
      precision if:
      - the user-specified input and result data type does not have enough
      precision to compute or store the result (i.e. the tolerance
      is less than \f$ 10^{-m} \f$ where \f$ m \f$ is the value,
      returned by <tt>numeric_limits<fp_t>::max_digits10</tt>).
      - the requested precision is near to or smaller than 1.0e-50, or
      - the function is noisy, non-deterministic, or is 
      discontinuous in the local neighborhood.

      If \ref verbose is 0, no output will be generated. If it is 1,
      then the tolerance and result will be output. If it is 2, then
      more diagnostics will be output.

      \note The algorithm attempts not to be wasteful, but is not
      necessarily optimized for speed. One way to improve it would be
      to more intelligently choose the number of digits used in the
      boost multiprecision numbers based on the tolerance which was
      specified. Another way to improve this class would be to use
      other multiprecision types beyond boost.
  */
  template<class fp_25_t, class fp_35_t, class fp_50_t, class fp_100_t>
  class funct_multip_tl {

  public:

    funct_multip_tl() {
      verbose=0;
      tol_rel=-1.0;
      err_nonconv=true;
    }

    ~funct_multip_tl() {
    }
    
    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief If true, call the error handler if the function
        evaluation fails
    */
    bool err_nonconv;
    
    /** \brief Evaluate the function and return the error estimate
        with the specified tolerance

        If \c tol_loc is positive and non-zero, then this is the
        relative tolerance used. If \c tol_loc is zero or negative 
        and \ref tol_rel is positive and non-zero, then \ref tol_rel
        is used for the relative tolerance. Otherwise, if both of
        these values is negative, then the default relative
        tolerance is used.
    */
    template<typename func_t, class fp_t>
    int eval_tol_err(func_t &&f, const fp_t &x, fp_t &val,
                     fp_t &err, double tol_loc=-1) const {

      // Tolerance choice and verification logic
      
      if (tol_loc<=0.0 && tol_rel<=0.0) {
        // Add one to the value returned by digits10 to get a
        // reasonable precision goal for the user-specified type.
        // This choice means that most exact or nearly exact functions
        // will require only two function evaluations.
        tol_loc=pow(10.0,-std::numeric_limits<fp_t>::digits10+1);
        if (verbose>0) {
          std::cout << "funct_multip_tl::eval_tol_err(): "
                    << "Set tolerance from data type to: "
                    << tol_loc << std::endl;
        }
      } else if (tol_loc<=0.0) {
        // If the data type is not sufficiently accurate to hold
        // the requested tolerance, then call the error handler
        if (tol_rel<pow(10.0,-std::numeric_limits<fp_t>::max_digits10)) {
          std::cerr << "Class data member tol_rel is " << tol_rel
                    << " but data type " << typeid(fp_t).name()
                    << " only stores "
                    << std::numeric_limits<fp_t>::digits10
                    << " digits." << std::endl;
          O2SCL_ERR("Cannot compute to required precision",
                    o2scl::exc_einval);
        }
        // If the user specified tol_rel, then use that value
        tol_loc=tol_rel;
        if (verbose>0) {
          std::cout << "funct_multip_tl::eval_tol_err(): "
                    << "Set tolerance from value of tol_rel to: "
                    << tol_loc << std::endl;
        }
      } else {
        // If the data type is not sufficiently accurate to hold
        // the requested tolerance, then call the error handler
        if (tol_loc<pow(10.0,-std::numeric_limits<fp_t>::max_digits10)) {
          std::cerr << "Caller requested tolerance " << tol_loc
                    << " but data type only stores "
                    << std::numeric_limits<fp_t>::digits10
                    << " digits." << std::endl;
          O2SCL_ERR("Cannot compute to required precision",
                    o2scl::exc_einval);
        } else if (verbose>0) {
          std::cout << "funct_multip_tl::eval_tol_err(): "
                    << "Set tolerance from user-specified value to: "
                    << tol_loc << std::endl;
        }
        // Use the value of tol_loc
      }

      /// First pass, compare double and long double

      bool d_eval=false;
      double x_d=0, y_d=0;
      if (tol_loc>pow(10.0,-std::numeric_limits<double>::digits10)) {
        x_d=static_cast<double>(x);
        y_d=f(x_d);
        d_eval=true;
      }
      
      bool ld_eval=false;
      long double x_ld=0, y_ld=0;
      if (tol_loc>pow(10.0,-std::numeric_limits<long double>::digits10)) {
        x_ld=static_cast<long double>(x);
        y_ld=f(x_ld);
        ld_eval=true;
      }

      if (d_eval && ld_eval) {
        if (y_ld==0 && y_d==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_tl::eval_tol_err() "
                      << "double and long double both got zero."
                      << std::endl;
          }
          return 0;
        }
      
        if (y_ld!=0) {
          err=static_cast<fp_t>(abs(y_ld-y_d)/abs(y_ld));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_ld);
            if (verbose>0) {
              std::cout << "funct_multip_tl::eval_tol_err() "
                        << "succeeded with double and long double:\n  "
                        << "arg,val,err: " << x << " "
                        << val << " " << err << std::endl;
              if (verbose>2) {
                std::cout << "funct_multip_tl::eval_tol_err() waiting "
                          << "for a key: " << std::flush;
                char ch;
                std::cin >> ch;
              }
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (d_eval && ld_eval) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed first round: " << dtos(y_ld,0) << " "
                    << dtos(y_d,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (ld_eval && verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed first round (d_eval is false): "
                    << dtos(y_ld,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed first round "
                    << "(d_eval and ld_eval are both false)." << std::endl;
        }
        if (verbose>2) {
          std::cout << "funct_multip_tl::eval_tol_err() waiting "
                    << "for a key: " << std::flush;
          char ch;
          std::cin >> ch;
        }
      }
    
      /// Second pass, compare long double and 25-digit precision
      
      bool cdf25_eval=false;
      fp_25_t x_cdf25=0, y_cdf25=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <fp_25_t>::digits10)) {
        x_cdf25=static_cast<fp_25_t>(x);
        y_cdf25=f(x_cdf25);
        cdf25_eval=true;
      }

      if (ld_eval && cdf25_eval) {
        if (y_cdf25==0 && y_ld==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_tl::eval_tol_err() "
                      << "long double and 25-digit both got zero."
                      << std::endl;
            if (verbose>2) {
              std::cout << "funct_multip_tl::eval_tol_err() waiting "
                        << "for a key: " << std::flush;
              char ch;
              std::cin >> ch;
            }
          }
          return 0;
        }
        if (y_cdf25!=0) {
	  fp_25_t temp1=abs(y_cdf25-static_cast<fp_25_t>(y_ld));
	  fp_25_t temp2=abs(y_cdf25);
          err=static_cast<fp_t>(temp1/temp2);
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf25);
            if (verbose>0) {
              std::cout << "funct_multip_tl::eval_tol_err() "
                        << "succeeded with long double and 25-digit:\n  "
                        << "arg,val,err: " << x << " "
                        << val << " " << err << std::endl;
              if (verbose>2) {
                std::cout << "funct_multip_tl::eval_tol_err() waiting "
                          << "for a key: " << std::flush;
                char ch;
                std::cin >> ch;
              }
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (ld_eval && cdf25_eval) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed second round: " << dtos(y_cdf25,0) << " "
                    << dtos(y_ld,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf25_eval && verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed second round (ld_eval is false): "
                    << dtos(y_cdf25,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed second round (ld_eval and "
                    << "cdf25_eval are both false)." << std::endl;
        }
        if (verbose>2) {
          std::cout << "funct_multip_tl::eval_tol_err() waiting "
                    << "for a key: " << std::flush;
          char ch;
          std::cin >> ch;
        }
      }
    
      /// Third pass, compare 25- and 35-digit precision

      bool cdf35_eval=false;
      fp_35_t x_cdf35=0, y_cdf35=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <fp_35_t>::digits10)) {
        x_cdf35=static_cast<fp_35_t>(x);
        y_cdf35=f(x_cdf35);
        cdf35_eval=true;
      }

      if (cdf25_eval && cdf35_eval) {
        if (y_cdf35==0 && y_cdf25==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_tl::eval_tol_err() "
                      << "25-digit and 35-digit both got zero."
                      << std::endl;
            if (verbose>2) {
              std::cout << "funct_multip_tl::eval_tol_err() waiting "
                        << "for a key: " << std::flush;
              char ch;
              std::cin >> ch;
            }
          }
          return 0;
        }
        if (y_cdf35!=0) {
	  fp_35_t temp1=abs(y_cdf35-static_cast<fp_35_t>(y_cdf25));
	  fp_35_t temp2=abs(y_cdf35);
          err=static_cast<fp_t>(temp1/temp2);
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf35);
            if (verbose>0) {
              std::cout << "funct_multip_tl::eval_tol_err() "
                        << "succeeded with 25-digit and 35-digit:\n  "
                        << "arg,val,err: " << x << " "
                        << val << " " << err << std::endl;
              if (verbose>2) {
                std::cout << "funct_multip_tl::eval_tol_err() waiting "
                          << "for a key: " << std::flush;
                char ch;
                std::cin >> ch;
              }
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf25_eval && cdf35_eval) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed third round: " << dtos(y_cdf35,0) << " "
                    << dtos(y_cdf25,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf35_eval && verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval is false): "
                    << dtos(y_cdf35,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval and "
                    << "cdf35_eval are both false)." << std::endl;
        }
        if (verbose>2) {
          std::cout << "funct_multip_tl::eval_tol_err() waiting "
                    << "for a key: " << std::flush;
          char ch;
          std::cin >> ch;
        }
      }
    
      /// Fourth pass, compare 35- and 50-digit precision
      
      bool cdf50_eval=false;
      fp_50_t x_cdf50=0, y_cdf50=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <fp_50_t>::digits10)) {
        x_cdf50=static_cast<fp_50_t>(x);
        y_cdf50=f(x_cdf50);
        cdf50_eval=true;
      }

      if (cdf35_eval && cdf50_eval) {
        if (y_cdf50==0 && y_cdf35==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_tl::eval_tol_err() "
                      << "35-digit and 50-digit both got zero."
                      << std::endl;
          }
          return 0;
        }
        if (y_cdf50!=0) {
	  fp_50_t temp1=abs(y_cdf50-static_cast<fp_50_t>(y_cdf35));
	  fp_50_t temp2=abs(y_cdf50);
          err=static_cast<fp_t>(temp1/temp2);
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf50);
            if (verbose>0) {
              std::cout << "funct_multip_tl::eval_tol_err() "
                        << "succeeded with 35-digit and 50-digit:\n  "
                        << "arg,val,err: " << x << " "
                        << val << " " << err << std::endl;
              if (verbose>2) {
                std::cout << "funct_multip_tl::eval_tol_err() waiting "
                          << "for a key: " << std::flush;
                char ch;
                std::cin >> ch;
              }
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf35_eval && cdf50_eval) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed fourth round: " << dtos(y_cdf50,0) << " "
                    << dtos(y_cdf35,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf50_eval && verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval is false): "
                    << dtos(y_cdf50,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval and "
                    << "cdf50_eval are both false)." << std::endl;
        }
        if (verbose>2) {
          std::cout << "funct_multip_tl::eval_tol_err() waiting "
                    << "for a key: " << std::flush;
          char ch;
          std::cin >> ch;
        }
      }
    
      /// Final pass, compare 50- and 100-digit precision
      
      bool cdf100_eval=false;
      fp_100_t x_cdf100=0, y_cdf100=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <fp_100_t>::digits10)) {
        x_cdf100=static_cast<fp_100_t>(x);
        y_cdf100=f(x_cdf100);
        cdf100_eval=true;
      }

      if (cdf100_eval && cdf50_eval) {
        if (y_cdf100==0 && y_cdf50==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_tl::eval_tol_err() "
                      << "50-digit and 100-digit both got zero."
                      << std::endl;
            if (verbose>2) {
              std::cout << "funct_multip_tl::eval_tol_err() waiting "
                        << "for a key: " << std::flush;
              char ch;
              std::cin >> ch;
            }
          }
          return 0;
        }
        if (y_cdf100!=0) {
	  fp_100_t temp1=abs(y_cdf100-static_cast<fp_100_t>(y_cdf50));
	  fp_100_t temp2=abs(y_cdf100);
          err=static_cast<fp_t>(temp1/temp2);
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf100);
            if (verbose>0) {
              std::cout << "funct_multip_tl::eval_tol_err() "
                        << "succeeded with 50-digit and 100-digit:\n  "
                        << "arg,val,err: " << x << " "
                        << val << " " << err << std::endl;
              if (verbose>2) {
                std::cout << "funct_multip_tl::eval_tol_err() waiting "
                          << "for a key: " << std::flush;
                char ch;
                std::cin >> ch;
              }
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf50_eval && cdf100_eval) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed last round: " << dtos(y_cdf100,0) << " "
                    << dtos(y_cdf50,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf100_eval && verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval is false): "
                    << dtos(y_cdf100,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_tl::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval and "
                    << "cdf100_eval are both false)." << std::endl;
          if (verbose>2) {
            std::cout << "funct_multip_tl::eval_tol_err() waiting "
                      << "for a key: " << std::flush;
            char ch;
            std::cin >> ch;
          }
        }
      }
    
      /// Algorithm failed
      
      O2SCL_CONV2("Failed to compute with requested accuracy ",
                  "in funct_multip_tl::eval_tol_err().",
                  o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Evaluate the function and return the error estimate
        with the default tolerance for the specified type
    */
    template<typename func_t, class fp_t>
    int eval_err(func_t &&f, const fp_t &x,
                 fp_t &val, fp_t &err) const {
      return eval_tol_err(f,x,val,err);
    }
  
    /** \brief Evalulate the function without an error estimate
     */
    template<typename func_t, class fp_t>
    fp_t operator()(func_t &&f, const fp_t &x) const {
      fp_t val;
      fp_t err;
      
      eval_err(f,x,val,err);
    
      return val;
    }
      
  };

  /** \brief The default multiprecision function object
   */
  typedef funct_multip_tl<o2fp_25,o2fp_35,o2fp_50,o2fp_100>
  funct_multip;

  /** \brief The multiprecision function object with \t cpp_dec_float
      types
  */
  typedef funct_multip_tl<cpp_dec_float_25,cpp_dec_float_35,
                          cpp_dec_float_50,cpp_dec_float_100>
  funct_multip_cdf;

#ifdef O2SCL_SET_MPFR
  
  /** \brief The multiprecision function object using \t mpfr 
      types
  */
  typedef funct_multip_tl<mpfr_25,mpfr_35,mpfr_50,mpfr_100>
  funct_multip_mpfr;
  
#endif
  
  /** \brief A multiprecision function evaluation class with 
      transformations useful for integrals

      This class is used in \ref inte_adapt_cern .
  */
  template<class lim_fp_t, class fp_25_t, class fp_35_t, class fp_50_t,
           class fp_100_t> class funct_multip_transform_tl {

  public:

    funct_multip_transform_tl() {
      verbose=0;
      tol_rel=-1.0;
      err_nonconv=true;
    }

    ~funct_multip_transform_tl() {
    }
    
    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief If true, call the error handler if the function
        evaluation fails
    */
    bool err_nonconv;

    /// The upper limit (when finite)
    lim_fp_t upper_lim;
    
    /// The lower limit (when finite)
    lim_fp_t lower_lim;
    
    /** \brief Evaluate the function and return the error estimate
        with the specified tolerance

        If \c tol_loc is positive and non-zero, then this is the
        relative tolerance used. If \c tol_loc is zero or negative 
        and \ref tol_rel is positive and non-zero, then \ref tol_rel
        is used for the relative tolerance. Otherwise, if both of
        these values is negative, then the default relative
        tolerance is used.
    */
    template<typename func_t, class fp_t>
    int eval_tol_err(char mode, func_t &&f, const fp_t &t, fp_t &val,
                     fp_t &err, double tol_loc=-1) const {

      // Tolerance choice and verification logic
      
      if (tol_loc<=0.0 && tol_rel<=0.0) {
        // Add one to the value returned by digits10 to get a
        // reasonable precision goal for the user-specified type.
        // This choice means that most exact or nearly exact functions
        // will require only two function evaluations.
        tol_loc=pow(10.0,-std::numeric_limits<fp_t>::digits10+1);
        if (verbose>0) {
          std::cout << "funct_multip_transform_tl::eval_tol_err(): "
                    << "Set tolerance from data type to: "
                    << tol_loc << std::endl;
        }
      } else if (tol_loc<=0.0) {
        // If the data type is not sufficiently accurate to hold
        // the requested tolerance, then call the error handler
        if (tol_rel<pow(10.0,-std::numeric_limits<fp_t>::max_digits10)) {
          std::cerr << "Class data member tol_rel is " << tol_rel
                    << " but data type only stores "
                    << std::numeric_limits<fp_t>::digits10
                    << " digits." << std::endl;
          O2SCL_ERR("Cannot compute to required precision",
                    o2scl::exc_einval);
        }
        // If the user specified tol_rel, then use that value
        tol_loc=tol_rel;
        if (verbose>0) {
          std::cout << "funct_multip_transform_tl::eval_tol_err(): "
                    << "Set tolerance from value of tol_rel to: "
                    << tol_loc << std::endl;
        }
      } else {
        // If the data type is not sufficiently accurate to hold
        // the requested tolerance, then call the error handler
        if (tol_loc<pow(10.0,-std::numeric_limits<fp_t>::max_digits10)) {
          std::cerr << "Caller requested tolerance " << tol_loc
                    << " but data type only stores "
                    << std::numeric_limits<fp_t>::digits10
                    << " digits." << std::endl;
          O2SCL_ERR("Cannot compute to required precision",
                    o2scl::exc_einval);
        } else if (verbose>0) {
          std::cout << "funct_multip_transform_tl::eval_tol_err(): "
                    << "Set tolerance from user-specified value to: "
                    << tol_loc << std::endl;
        }
        // Use the value of tol_loc
      }

      /// First pass, compare double and long double

      bool d_eval=false;
      double t_d=0, y_d=0;
      
      if (tol_loc>pow(10.0,-std::numeric_limits<double>::digits10)) {
        t_d=static_cast<double>(t);
        double x_d;
        if (mode=='u') {
          x_d=static_cast<double>(lower_lim)+(1-t_d)/t_d;
          y_d=f(x_d)/t_d/t_d;
        } else if (mode=='l') {
          x_d=static_cast<double>(upper_lim)-(1-t_d)/t_d;
          y_d=f(x_d)/t_d/t_d;
        } else {
          x_d=(1-t_d)/t_d;
          y_d=(f(x_d)+f(-x_d))/t_d/t_d;
        }
        d_eval=true;
      }
      
      bool ld_eval=false;
      long double x_ld=0, y_ld=0;
      if (tol_loc>pow(10.0,-std::numeric_limits<long double>::digits10)) {
        long double t_ld=static_cast<long double>(t);
        if (mode=='u') {
          x_ld=static_cast<long double>(lower_lim)+(1-t_ld)/t_ld;
          y_ld=f(x_ld)/t_ld/t_ld;
        } else if (mode=='l') {
          x_ld=static_cast<long double>(upper_lim)-(1-t_ld)/t_ld;
          y_ld=f(x_ld)/t_ld/t_ld;
        } else {
          x_ld=(1-t_ld)/t_ld;
          y_ld=(f(x_ld)+f(-x_ld))/t_ld/t_ld;
        }
        ld_eval=true;
      }

      if (d_eval && ld_eval) {
        if (y_ld==0 && y_d==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform_tl::eval_tol_err() "
                      << "double and long double both got zero."
                      << std::endl;
          }
          return 0;
        }
      
        if (y_ld!=0) {
          err=static_cast<fp_t>(abs(y_ld-y_d)/abs(y_ld));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_ld);
            if (verbose>0) {
              std::cout << "funct_multip_transform_tl::eval_tol_err() "
                        << "succeeded with double and long double:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (d_eval && ld_eval) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed first round: " << dtos(y_ld,0) << " "
                    << dtos(y_d,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (ld_eval && verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed first round (d_eval is false): "
                    << dtos(y_ld,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed first round "
                    << "(d_eval and ld_eval are both false)." << std::endl;
        }
      }
    
      /// Second pass, compare long double and 25-digit precision
      
      bool cdf25_eval=false;
      fp_25_t x_cdf25=0, y_cdf25=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <fp_25_t>::digits10)) {
        fp_25_t t_cdf25=static_cast<fp_25_t>(t);
        if (mode=='u') {
          x_cdf25=static_cast<fp_25_t>(lower_lim)+
            (1-t_cdf25)/t_cdf25;
          y_cdf25=f(x_cdf25)/t_cdf25/t_cdf25;
        } else if (mode=='l') {
          x_cdf25=static_cast<fp_25_t>(upper_lim)-
            (1-t_cdf25)/t_cdf25;
          y_cdf25=f(x_cdf25)/t_cdf25/t_cdf25;
        } else {
          x_cdf25=(1-t_cdf25)/t_cdf25;
          fp_25_t res_p=f(x_cdf25);
          fp_25_t x2=-x_cdf25;
          fp_25_t res_m=f(x2);
          y_cdf25=(res_p+res_m)/t_cdf25/t_cdf25;
        }
        cdf25_eval=true;
      }

      if (ld_eval && cdf25_eval) {
        if (y_cdf25==0 && y_ld==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform_tl::eval_tol_err() "
                      << "long double and 25-digit both got zero."
                      << std::endl;
          }
          return 0;
        }
        if (y_cdf25!=0) {
          err=static_cast<fp_t>(abs(y_cdf25-y_ld)/abs(y_cdf25));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf25);
            if (verbose>0) {
              std::cout << "funct_multip_transform_tl::eval_tol_err() "
                        << "succeeded with long double and 25-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (ld_eval && cdf25_eval) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed second round: " << dtos(y_cdf25,0) << " "
                    << dtos(y_ld,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf25_eval && verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed second round (ld_eval is false): "
                    << dtos(y_cdf25,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed second round (ld_eval and "
                    << "cdf25_eval are both false)." << std::endl;
        }
      }
    
      /// Third pass, compare 25- and 35-digit precision

      bool cdf35_eval=false;
      fp_35_t x_cdf35=0, y_cdf35=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <fp_35_t>::digits10)) {
        fp_35_t t_cdf35=static_cast<fp_35_t>(t);
        if (mode=='u') {
          x_cdf35=static_cast<fp_35_t>(lower_lim)+
            (1-t_cdf35)/t_cdf35;
          y_cdf35=f(x_cdf35)/t_cdf35/t_cdf35;
        } else if (mode=='l') {
          x_cdf35=static_cast<fp_35_t>(upper_lim)-
            (1-t_cdf35)/t_cdf35;
          y_cdf35=f(x_cdf35)/t_cdf35/t_cdf35;
        } else {
          x_cdf35=(1-t_cdf35)/t_cdf35;
          fp_35_t res_p=f(x_cdf35);
          fp_35_t x2=-x_cdf35;
          fp_35_t res_m=f(x2);
          y_cdf35=(res_p+res_m)/t_cdf35/t_cdf35;
        }
        cdf35_eval=true;
      }

      if (cdf25_eval && cdf35_eval) {
        if (y_cdf35==0 && y_cdf25==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform_tl::eval_tol_err() "
                      << "25-digit and 35-digit both got zero."
                      << std::endl;
          }
          return 0;
        }
        if (y_cdf35!=0) {
          err=static_cast<fp_t>(abs(y_cdf35-y_cdf25)/abs(y_cdf35));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf35);
            if (verbose>0) {
              std::cout << "funct_multip_transform_tl::eval_tol_err() "
                        << "succeeded with 25-digit and 35-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf25_eval && cdf35_eval) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed third round: " << dtos(y_cdf35,0) << " "
                    << dtos(y_cdf25,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf35_eval && verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval is false): "
                    << dtos(y_cdf35,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval and "
                    << "cdf35_eval are both false)." << std::endl;
        }
      }
    
      /// Fourth pass, compare 35- and 50-digit precision
      
      bool cdf50_eval=false;
      fp_50_t x_cdf50=0, y_cdf50=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <fp_50_t>::digits10)) {
        fp_50_t t_cdf50=static_cast<fp_50_t>(t);
        if (mode=='u') {
          x_cdf50=static_cast<fp_50_t>(lower_lim)+
            (1-t_cdf50)/t_cdf50;
          y_cdf50=f(x_cdf50)/t_cdf50/t_cdf50;
        } else if (mode=='l') {
          x_cdf50=static_cast<fp_50_t>(upper_lim)-
            (1-t_cdf50)/t_cdf50;
          y_cdf50=f(x_cdf50)/t_cdf50/t_cdf50;
        } else {
          x_cdf50=(1-t_cdf50)/t_cdf50;
          fp_50_t res_p=f(x_cdf50);
          fp_50_t x2=-x_cdf50;
          fp_50_t res_m=f(x2);
          y_cdf50=(res_p+res_m)/t_cdf50/t_cdf50;
        }
        cdf50_eval=true;
      }

      if (cdf35_eval && cdf50_eval) {
        if (y_cdf50==0 && y_cdf35==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform_tl::eval_tol_err() "
                      << "35-digit and 50-digit both got zero."
                      << std::endl;
          }
          return 0;
        }
        if (y_cdf50!=0) {
          err=static_cast<fp_t>(abs(y_cdf50-y_cdf35)/abs(y_cdf50));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf50);
            if (verbose>0) {
              std::cout << "funct_multip_transform_tl::eval_tol_err() "
                        << "succeeded with 35-digit and 50-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf35_eval && cdf50_eval) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed fourth round: " << dtos(y_cdf50,0) << " "
                    << dtos(y_cdf35,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf50_eval && verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval is false): "
                    << dtos(y_cdf50,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval and "
                    << "cdf50_eval are both false)." << std::endl;
        }
      }
    
      /// Final pass, compare 50- and 100-digit precision
      
      bool cdf100_eval=false;
      fp_100_t x_cdf100=0, y_cdf100=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <fp_100_t>::digits10)) {
        fp_100_t t_cdf100=static_cast<fp_100_t>(t);
        if (mode=='u') {
          x_cdf100=static_cast<fp_100_t>(lower_lim)+
            (1-t_cdf100)/t_cdf100;
          y_cdf100=f(x_cdf100)/t_cdf100/t_cdf100;
        } else if (mode=='l') {
          x_cdf100=static_cast<fp_100_t>(upper_lim)-
            (1-t_cdf100)/t_cdf100;
          y_cdf100=f(x_cdf100)/t_cdf100/t_cdf100;
        } else {
          x_cdf100=(1-t_cdf100)/t_cdf100;
          fp_100_t res_p=f(x_cdf100);
          fp_100_t x2=-x_cdf100;
          fp_100_t res_m=f(x2);
          y_cdf100=(res_p+res_m)/t_cdf100/t_cdf100;
        }
        cdf100_eval=true;
      }

      if (cdf100_eval && cdf50_eval) {
        if (y_cdf100==0 && y_cdf50==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform_tl::eval_tol_err() "
                      << "50-digit and 100-digit both got zero."
                      << std::endl;
          }
          return 0;
        }
        if (y_cdf100!=0) {
          err=static_cast<fp_t>(abs(y_cdf100-y_cdf50)/abs(y_cdf100));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf100);
            if (verbose>0) {
              std::cout << "funct_multip_transform_tl::eval_tol_err() "
                        << "succeeded with 50-digit and 100-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf50_eval && cdf100_eval) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed last round: " << dtos(y_cdf100,0) << " "
                    << dtos(y_cdf50,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf100_eval && verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval is false): "
                    << dtos(y_cdf100,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform_tl::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval and "
                    << "cdf100_eval are both false)." << std::endl;
        }
      }
    
      /// Algorithm failed
      
      O2SCL_CONV2("Failed to compute with requested accuracy ",
                  "in funct_multip_transform_tl::eval_tol_err().",
                  o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Evaluate the function and uncertainty with a
        transformed upper limit with the default tolerance for the
        specified type
    */
    template<typename func_t, class fp_t>
    int eval_err_iu(func_t &&f, const fp_t &x,
                    fp_t &val, fp_t &err) const {
      return eval_tol_err('u',f,x,val,err);
    }
  
    /** \brief Evaluate the function and uncertainty with a
        transformed lower limit with the default tolerance for the
        specified type
    */
    template<typename func_t, class fp_t>
    int eval_err_il(func_t &&f, const fp_t &x,
                    fp_t &val, fp_t &err) const {
      return eval_tol_err('l',f,x,val,err);
    }
  
    /** \brief Evaluate the function and uncertainty with a
        transformed limits with the default tolerance for the
        specified type
    */
    template<typename func_t, class fp_t>
    int eval_err_i(func_t &&f, const fp_t &x,
                   fp_t &val, fp_t &err) const {
      return eval_tol_err('i',f,x,val,err);
    }
  
    /** \brief Evaluate the function with a transformed upper limit
        with the default tolerance for the specified type
    */
    template<typename func_t, class fp_t>
    fp_t eval_iu(func_t &&f, const fp_t &x) const {
      fp_t val;
      fp_t err;
      
      eval_err_iu(f,x,val,err);
    
      return val;
    }
      
    /** \brief Evaluate the function with a transformed lower limit
        with the default tolerance for the specified type
    */
    template<typename func_t, class fp_t>
    fp_t eval_il(func_t &&f, const fp_t &x) const {
      fp_t val;
      fp_t err;
      
      eval_err_il(f,x,val,err);
    
      return val;
    }
      
    /** \brief Evaluate the function with a transformed limits with
        the default tolerance for the specified type
    */
    template<typename func_t, class fp_t>
    fp_t eval_i(func_t &&f, const fp_t &x) const {
      fp_t val;
      fp_t err;
      
      eval_err_i(f,x,val,err);
    
      return val;
    }
      
  };

  /// Alias declarations using default types
  template <class fp_t> using funct_multip_transform=
    funct_multip_transform_tl<fp_t,o2fp_25,o2fp_35,o2fp_50,o2fp_100>;
  
  /// Alias declarations for \t cpp_dec_float types
  template <class fp_t> using funct_multip_transform_cdf=
    funct_multip_transform_tl<double,cpp_dec_float_25,cpp_dec_float_35,
                              cpp_dec_float_50,cpp_dec_float_100>;

#ifdef O2SCL_SET_MPFR  
  /// Alias declarations for \t mpfr types
  template <class fp_t> using funct_multip_transform_mpfr=
    funct_multip_transform_tl<double,mpfr_25,mpfr_35,mpfr_50,mpfr_100>;
#endif

  // end of #ifndef O2SCL_NO_BOOST_MULTIPRECISION
#endif
  
}

#endif
