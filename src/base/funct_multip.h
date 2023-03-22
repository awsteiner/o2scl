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
#ifndef O2SCL_FUNCT_MULTIP_H
#define O2SCL_FUNCT_MULTIP_H

/** \file funct.h
    \brief Multiprecisions extension to Function object classes
*/

#include <boost/multiprecision/cpp_dec_float.hpp>
#ifdef O2SCL_MPFR
#include <boost/multiprecision/mpfr.hpp>
#endif

namespace o2scl {

  /// One-dimensional function typedef in src/base/funct.h
  typedef std::function<long double(long double)> funct_ld;

  /** \brief One-dimensional function typedef in src/base/funct.h
   */
  typedef std::function<boost::multiprecision::number<
                          boost::multiprecision::cpp_dec_float<25> >
			(boost::multiprecision::number<
                         boost::multiprecision::cpp_dec_float<25> > )>
  funct_cdf25;
  
  /** \brief One-dimensional function typedef in src/base/funct.h
   */
  typedef std::function<boost::multiprecision::number<
                          boost::multiprecision::cpp_dec_float<35> >
			(boost::multiprecision::number<
                         boost::multiprecision::cpp_dec_float<35> > )>
  funct_cdf35;
  
  /** \brief One-dimensional function typedef in src/base/funct.h
   */
  typedef std::function<boost::multiprecision::cpp_dec_float_50
			(boost::multiprecision::cpp_dec_float_50)>
  funct_cdf50;

#ifdef O2SCL_MPFR
  
  /** \brief One-dimensional function typedef in src/base/funct.h
      
      This typedef is only present if O2SCL_MPFR is
      defined during compilation.
  */
  typedef std::function<boost::multiprecision::mpfr_float_50
			(boost::multiprecision::mpfr_float_50)>
  funct_mp50;
  
#endif
  
  /** \brief One-dimensional function typedef in src/base/funct.h
   */
  typedef std::function<boost::multiprecision::cpp_dec_float_100
			(boost::multiprecision::cpp_dec_float_100)>
  funct_cdf100;
  
  typedef std::function<int(long double,long double &)> funct_ret_ld;

  typedef std::function<int(boost::multiprecision::number<
                            boost::multiprecision::cpp_dec_float<25> >,
                            boost::multiprecision::number<
                            boost::multiprecision::cpp_dec_float<25> > &)>
  funct_ret_cdf25;
  
  typedef std::function<int(boost::multiprecision::number<
                            boost::multiprecision::cpp_dec_float<35> >,
                            boost::multiprecision::number<
                            boost::multiprecision::cpp_dec_float<35> > &)>
  funct_ret_cdf35;
  
  typedef std::function<int(boost::multiprecision::number<
                            boost::multiprecision::cpp_dec_float<50> >,
                            boost::multiprecision::number<
                            boost::multiprecision::cpp_dec_float<50> > &)>
  funct_ret_cdf50;
  
  typedef std::function<int(boost::multiprecision::number<
                            boost::multiprecision::cpp_dec_float<100> >,
                            boost::multiprecision::number<
                            boost::multiprecision::cpp_dec_float<100> > &)>
  funct_ret_cdf100;

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
  class funct_multip {

  protected:
    
    /// \name Typedefs for multiprecision types
    //@{
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> >
    cpp_dec_float_25;
  
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> >
    cpp_dec_float_35;
  
    typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
  
    typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;
    //@}

  public:

    funct_multip() {
      verbose=0;
      tol_rel=-1.0;
      err_nonconv=true;
    }

    ~funct_multip() {
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
          std::cout << "funct_multip::eval_tol_err(): "
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
          std::cout << "funct_multip::eval_tol_err(): "
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
          std::cout << "funct_multip::eval_tol_err(): "
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
            std::cout << "funct_multip::eval_tol_err() "
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
              std::cout << "funct_multip::eval_tol_err() "
                        << "succeeded with double and long double:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (d_eval && ld_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed first round: " << dtos(y_ld,0) << " "
                    << dtos(y_d,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (ld_eval && verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed first round (d_eval is false): "
                    << dtos(y_ld,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed first round "
                    << "(d_eval and ld_eval are both false)." << std::endl;
        }
      }
    
      /// Second pass, compare long double and 25-digit precision
      
      bool cdf25_eval=false;
      cpp_dec_float_25 x_cdf25=0, y_cdf25=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_25>::digits10)) {
        x_cdf25=static_cast<cpp_dec_float_25>(x);
        y_cdf25=f(x_cdf25);
        cdf25_eval=true;
      }

      if (ld_eval && cdf25_eval) {
        if (y_cdf25==0 && y_ld==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip::eval_tol_err() "
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
              std::cout << "funct_multip::eval_tol_err() "
                        << "succeeded with long double and 25-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (ld_eval && cdf25_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed second round: " << dtos(y_cdf25,0) << " "
                    << dtos(y_ld,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf25_eval && verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed second round (ld_eval is false): "
                    << dtos(y_cdf25,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed second round (ld_eval and "
                    << "cdf25_eval are both false)." << std::endl;
        }
      }
    
      /// Third pass, compare 25- and 35-digit precision

      bool cdf35_eval=false;
      cpp_dec_float_35 x_cdf35=0, y_cdf35=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_35>::digits10)) {
        x_cdf35=static_cast<cpp_dec_float_35>(x);
        y_cdf35=f(x_cdf35);
        cdf35_eval=true;
      }

      if (cdf25_eval && cdf35_eval) {
        if (y_cdf35==0 && y_cdf25==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip::eval_tol_err() "
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
              std::cout << "funct_multip::eval_tol_err() "
                        << "succeeded with 25-digit and 35-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf25_eval && cdf35_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed third round: " << dtos(y_cdf35,0) << " "
                    << dtos(y_cdf25,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf35_eval && verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval is false): "
                    << dtos(y_cdf35,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval and "
                    << "cdf35_eval are both false)." << std::endl;
        }
      }
    
      /// Fourth pass, compare 35- and 50-digit precision
      
      bool cdf50_eval=false;
      cpp_dec_float_50 x_cdf50=0, y_cdf50=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_50>::digits10)) {
        x_cdf50=static_cast<cpp_dec_float_50>(x);
        y_cdf50=f(x_cdf50);
        cdf50_eval=true;
      }

      if (cdf35_eval && cdf50_eval) {
        if (y_cdf50==0 && y_cdf35==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip::eval_tol_err() "
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
              std::cout << "funct_multip::eval_tol_err() "
                        << "succeeded with 35-digit and 50-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf35_eval && cdf50_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed fourth round: " << dtos(y_cdf50,0) << " "
                    << dtos(y_cdf35,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf50_eval && verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval is false): "
                    << dtos(y_cdf50,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval and "
                    << "cdf50_eval are both false)." << std::endl;
        }
      }
    
      /// Final pass, compare 50- and 100-digit precision
      
      bool cdf100_eval=false;
      cpp_dec_float_100 x_cdf100=0, y_cdf100=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_100>::digits10)) {
        x_cdf100=static_cast<cpp_dec_float_100>(x);
        y_cdf100=f(x_cdf100);
        cdf100_eval=true;
      }

      if (cdf100_eval && cdf50_eval) {
        if (y_cdf100==0 && y_cdf50==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip::eval_tol_err() "
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
              std::cout << "funct_multip::eval_tol_err() "
                        << "succeeded with 50-digit and 100-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf50_eval && cdf100_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed last round: " << dtos(y_cdf100,0) << " "
                    << dtos(y_cdf50,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf100_eval && verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval is false): "
                    << dtos(y_cdf100,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval and "
                    << "cdf100_eval are both false)." << std::endl;
        }
      }
    
      /// Algorithm failed
      
      O2SCL_CONV2("Failed to compute with requested accuracy ",
                  "in funct_multip::eval_tol_err().",
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

  template<class lim_fp_t> class funct_multip_transform {

  protected:
    
    /// \name Typedefs for multiprecision types
    //@{
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> >
    cpp_dec_float_25;
  
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> >
    cpp_dec_float_35;
  
    typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
  
    typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;
    //@}

  public:

    funct_multip_transform() {
      verbose=0;
      tol_rel=-1.0;
      err_nonconv=true;
    }

    ~funct_multip_transform() {
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

    /// Desc
    lim_fp_t upper_lim;
    
    /// Desc
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
          std::cout << "funct_multip_transform::eval_tol_err(): "
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
          std::cout << "funct_multip_transform::eval_tol_err(): "
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
          std::cout << "funct_multip_transform::eval_tol_err(): "
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
        long double x_ld;
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
            std::cout << "funct_multip_transform::eval_tol_err() "
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
              std::cout << "funct_multip_transform::eval_tol_err() "
                        << "succeeded with double and long double:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (d_eval && ld_eval) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed first round: " << dtos(y_ld,0) << " "
                    << dtos(y_d,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (ld_eval && verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed first round (d_eval is false): "
                    << dtos(y_ld,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed first round "
                    << "(d_eval and ld_eval are both false)." << std::endl;
        }
      }
    
      /// Second pass, compare long double and 25-digit precision
      
      bool cdf25_eval=false;
      cpp_dec_float_25 x_cdf25=0, y_cdf25=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_25>::digits10)) {
        cpp_dec_float_25 t_cdf25=static_cast<cpp_dec_float_25>(t);
        cpp_dec_float_25 x_cdf25;
        if (mode=='u') {
          x_cdf25=static_cast<cpp_dec_float_25>(lower_lim)+
            (1-t_cdf25)/t_cdf25;
          y_cdf25=f(x_cdf25)/t_cdf25/t_cdf25;
        } else if (mode=='l') {
          x_cdf25=static_cast<cpp_dec_float_25>(upper_lim)-
            (1-t_cdf25)/t_cdf25;
          y_cdf25=f(x_cdf25)/t_cdf25/t_cdf25;
        } else {
          x_cdf25=(1-t_cdf25)/t_cdf25;
          y_cdf25=(f(x_cdf25)+f(-x_cdf25))/t_cdf25/t_cdf25;
        }
        cdf25_eval=true;
      }

      if (ld_eval && cdf25_eval) {
        if (y_cdf25==0 && y_ld==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform::eval_tol_err() "
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
              std::cout << "funct_multip_transform::eval_tol_err() "
                        << "succeeded with long double and 25-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (ld_eval && cdf25_eval) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed second round: " << dtos(y_cdf25,0) << " "
                    << dtos(y_ld,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf25_eval && verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed second round (ld_eval is false): "
                    << dtos(y_cdf25,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed second round (ld_eval and "
                    << "cdf25_eval are both false)." << std::endl;
        }
      }
    
      /// Third pass, compare 25- and 35-digit precision

      bool cdf35_eval=false;
      cpp_dec_float_35 x_cdf35=0, y_cdf35=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_35>::digits10)) {
        cpp_dec_float_35 t_cdf35=static_cast<cpp_dec_float_35>(t);
        cpp_dec_float_35 x_cdf35;
        if (mode=='u') {
          x_cdf35=static_cast<cpp_dec_float_35>(lower_lim)+
            (1-t_cdf35)/t_cdf35;
          y_cdf35=f(x_cdf35)/t_cdf35/t_cdf35;
        } else if (mode=='l') {
          x_cdf35=static_cast<cpp_dec_float_35>(upper_lim)-
            (1-t_cdf35)/t_cdf35;
          y_cdf35=f(x_cdf35)/t_cdf35/t_cdf35;
        } else {
          x_cdf35=(1-t_cdf35)/t_cdf35;
          y_cdf35=(f(x_cdf35)+f(-x_cdf35))/t_cdf35/t_cdf35;
        }
        cdf35_eval=true;
      }

      if (cdf25_eval && cdf35_eval) {
        if (y_cdf35==0 && y_cdf25==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform::eval_tol_err() "
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
              std::cout << "funct_multip_transform::eval_tol_err() "
                        << "succeeded with 25-digit and 35-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf25_eval && cdf35_eval) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed third round: " << dtos(y_cdf35,0) << " "
                    << dtos(y_cdf25,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf35_eval && verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval is false): "
                    << dtos(y_cdf35,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval and "
                    << "cdf35_eval are both false)." << std::endl;
        }
      }
    
      /// Fourth pass, compare 35- and 50-digit precision
      
      bool cdf50_eval=false;
      cpp_dec_float_50 x_cdf50=0, y_cdf50=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_50>::digits10)) {
        cpp_dec_float_50 t_cdf50=static_cast<cpp_dec_float_50>(t);
        cpp_dec_float_50 x_cdf50;
        if (mode=='u') {
          x_cdf50=static_cast<cpp_dec_float_50>(lower_lim)+
            (1-t_cdf50)/t_cdf50;
          y_cdf50=f(x_cdf50)/t_cdf50/t_cdf50;
        } else if (mode=='l') {
          x_cdf50=static_cast<cpp_dec_float_50>(upper_lim)-
            (1-t_cdf50)/t_cdf50;
          y_cdf50=f(x_cdf50)/t_cdf50/t_cdf50;
        } else {
          x_cdf50=(1-t_cdf50)/t_cdf50;
          y_cdf50=(f(x_cdf50)+f(-x_cdf50))/t_cdf50/t_cdf50;
        }
        cdf50_eval=true;
      }

      if (cdf35_eval && cdf50_eval) {
        if (y_cdf50==0 && y_cdf35==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform::eval_tol_err() "
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
              std::cout << "funct_multip_transform::eval_tol_err() "
                        << "succeeded with 35-digit and 50-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf35_eval && cdf50_eval) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed fourth round: " << dtos(y_cdf50,0) << " "
                    << dtos(y_cdf35,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf50_eval && verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval is false): "
                    << dtos(y_cdf50,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval and "
                    << "cdf50_eval are both false)." << std::endl;
        }
      }
    
      /// Final pass, compare 50- and 100-digit precision
      
      bool cdf100_eval=false;
      cpp_dec_float_100 x_cdf100=0, y_cdf100=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_100>::digits10)) {
        cpp_dec_float_100 t_cdf100=static_cast<cpp_dec_float_100>(t);
        cpp_dec_float_100 x_cdf100;
        if (mode=='u') {
          x_cdf100=static_cast<cpp_dec_float_100>(lower_lim)+
            (1-t_cdf100)/t_cdf100;
          y_cdf100=f(x_cdf100)/t_cdf100/t_cdf100;
        } else if (mode=='l') {
          x_cdf100=static_cast<cpp_dec_float_100>(upper_lim)-
            (1-t_cdf100)/t_cdf100;
          y_cdf100=f(x_cdf100)/t_cdf100/t_cdf100;
        } else {
          x_cdf100=(1-t_cdf100)/t_cdf100;
          y_cdf100=(f(x_cdf100)+f(-x_cdf100))/t_cdf100/t_cdf100;
        }
        cdf100_eval=true;
      }

      if (cdf100_eval && cdf50_eval) {
        if (y_cdf100==0 && y_cdf50==0) {
          val=0;
          err=0;
          if (verbose>0) {
            std::cout << "funct_multip_transform::eval_tol_err() "
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
              std::cout << "funct_multip_transform::eval_tol_err() "
                        << "succeeded with 50-digit and 100-digit:\n  "
                        << "val,err: " << val << " " << err << std::endl;
            }
            return 0;
          }
        }
      }
      
      if (verbose>0) {
        if (cdf50_eval && cdf100_eval) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed last round: " << dtos(y_cdf100,0) << " "
                    << dtos(y_cdf50,0) << " "
                    << dtos(err,0) << " " << tol_loc << std::endl;
        } else if (cdf100_eval && verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval is false): "
                    << dtos(y_cdf100,0) << " "
                    << tol_loc << std::endl;
        } else if (verbose>1) {
          std::cout << "funct_multip_transform::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval and "
                    << "cdf100_eval are both false)." << std::endl;
        }
      }
    
      /// Algorithm failed
      
      O2SCL_CONV2("Failed to compute with requested accuracy ",
                  "in funct_multip_transform::eval_tol_err().",
                  o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Evaluate the function and return the error estimate
        with the default tolerance for the specified type
     */
    template<typename func_t, class fp_t>
    int eval_err_iu(func_t &&f, const fp_t &x,
                 fp_t &val, fp_t &err) const {
      return eval_tol_err('u',f,x,val,err);
    }
  
    template<typename func_t, class fp_t>
    int eval_err_il(func_t &&f, const fp_t &x,
                 fp_t &val, fp_t &err) const {
      return eval_tol_err('l',f,x,val,err);
    }
  
    template<typename func_t, class fp_t>
    int eval_err_i(func_t &&f, const fp_t &x,
                 fp_t &val, fp_t &err) const {
      return eval_tol_err('i',f,x,val,err);
    }
  
    /** \brief Evalulate the function without an error estimate
     */
    template<typename func_t, class fp_t>
    fp_t eval_iu(func_t &&f, const fp_t &x) const {
      fp_t val;
      fp_t err;
      
      eval_err_iu(f,x,val,err);
    
      return val;
    }
      
    template<typename func_t, class fp_t>
    fp_t eval_il(func_t &&f, const fp_t &x) const {
      fp_t val;
      fp_t err;
      
      eval_err_il(f,x,val,err);
    
      return val;
    }
      
    template<typename func_t, class fp_t>
    fp_t eval_i(func_t &&f, const fp_t &x) const {
      fp_t val;
      fp_t err;
      
      eval_err_i(f,x,val,err);
    
      return val;
    }
      
  };

  /** \brief Evaluate a one-dimensional function from a string
      at multiprecision

      \note Experimental.

      \warning This class only supports a limited number of data
      types, including double, long double, and cpp_dec_float types
      with 25, 35, 50, or 100 digits. It is designed to be used with
      the \ref funct_multip class.
   */
  class funct_multip_string {

  protected:
    
    /// \name Typedefs for multiprecision types
    //@{
    typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<25> > cpp_dec_float_25;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<35> > cpp_dec_float_35;
    typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
    typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;
    //@}

    /// \name The function evaluation objects
    //@{
    calc_utf8<double> c;
    calc_utf8<long double> c_ld;
    calc_utf8<cpp_dec_float_25> c_25;
    calc_utf8<cpp_dec_float_35> c_35;
    calc_utf8<cpp_dec_float_50> c_50;
    calc_utf8<cpp_dec_float_100> c_100;
    //@}

    /// \name The unit conversion objects
    //@{
    convert_units<double> cu;
    convert_units<long double> cu_ld;
    convert_units<cpp_dec_float_25> cu_25;
    convert_units<cpp_dec_float_35> cu_35;
    convert_units<cpp_dec_float_50> cu_50;
    convert_units<cpp_dec_float_100> cu_100;
    //@}

    /// \name The variable lists
    //@{
    std::map<std::string,double> vars;
    std::map<std::string,long double> vars_ld;
    std::map<std::string,cpp_dec_float_25> vars_25;
    std::map<std::string,cpp_dec_float_35> vars_35;
    std::map<std::string,cpp_dec_float_50> vars_50;
    std::map<std::string,cpp_dec_float_100> vars_100;
    //@}

    /** \brief If true, then the most recent function has been
        compiled in all of the function evaluation objects
    */
    bool compiled;
    
    /// The expression to be evaluated
    std::string st_form;
    
    /// The variable
    std::string st_var;

  public:
    
    funct_multip_string() {
      verbose=0;
      err_nonconv=true;
      compiled=false;
    }

    virtual ~funct_multip_string() {
    }
    
    /** \brief Set the function to compute
     */
    int set_function(std::string expr, std::string var) {
      st_form=expr;
      st_var=var;
      compiled=false;
      return 0;
    }

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief If true, call the error handler if the function
        evaluation fails
     */
    bool err_nonconv;

    /** \brief Compute the function at the value \c x 
     */
    template<class fp_t> fp_t operator()(fp_t x) {
    
      if (compiled==false) {

        if (verbose>3) {
          c.verbose=2;
          c_ld.verbose=2;
          c_25.verbose=2;
          c_35.verbose=2;
          c_50.verbose=2;
          c_100.verbose=2;
        } else if (verbose>2) {
          c.verbose=1;
          c_ld.verbose=1;
          c_25.verbose=1;
          c_35.verbose=1;
          c_50.verbose=1;
          c_100.verbose=1;
        }
        
        if (verbose>1) {
          std::cout << "funct_multip_string::operator() "
                    << "compiling with function "
                    << st_form << " and variable " << st_var
                    << std::endl;
        }
        
        c.compile(st_form.c_str());
        c_ld.compile(st_form.c_str());
        c_25.compile(st_form.c_str());
        c_35.compile(st_form.c_str());
        c_50.compile(st_form.c_str());
        c_100.compile(st_form.c_str());
        
        std::vector<std::u32string> vs=c.get_var_list();
        
        // If there are undefined variables, then attempt to get them
        // from the constant database
        if (vs.size()!=0) {
          
          for(size_t i=0;i<vs.size();i++) {
            
            std::string vsi2;
            char32_to_utf8(vs[i],vsi2);

            if (vsi2!=st_var) {
            
              if (verbose>1) {
                std::cout << "funct_multip_string::operator() "
                          << "trying to find constant " << vsi2
                          << std::endl;
              }
              
              std::vector<typename find_constants<
                double>::const_entry> matches;
              int fret=cu.find_nothrow(vsi2,"mks",matches);
              
              if (fret==find_constants<
                  double>::one_exact_match_unit_match ||
                  fret==find_constants<
                  double>::one_pattern_match_unit_match) {
                
                find_constants<double>::const_entry &fcl=matches[0];
                vars.insert(std::make_pair(vsi2,fcl.val));
                
                std::vector<typename
                            find_constants<long double>::const_entry>
                  matches_ld;
                cu_ld.find_nothrow(vsi2,"mks",matches_ld);
                vars_ld.insert(std::make_pair(vsi2,matches_ld[0].val));
                
                std::vector<typename
                            find_constants<cpp_dec_float_25>::const_entry>
                  matches_25;
                cu_25.find_nothrow(vsi2,"mks",matches_25);
                vars_25.insert(std::make_pair(vsi2,matches_25[0].val));
                
                std::vector<typename
                            find_constants<cpp_dec_float_35>::const_entry>
                  matches_35;
                cu_35.find_nothrow(vsi2,"mks",matches_35);
                vars_35.insert(std::make_pair(vsi2,matches_35[0].val));
                
                std::vector<typename
                            find_constants<cpp_dec_float_50>::const_entry>
                  matches_50;
                cu_50.find_nothrow(vsi2,"mks",matches_50);
                vars_50.insert(std::make_pair(vsi2,matches_50[0].val));
                
                std::vector<typename
                            find_constants<cpp_dec_float_100>::const_entry>
                  matches_100;
                cu_100.find_nothrow(vsi2,"mks",matches_100);
                vars_100.insert(std::make_pair(vsi2,matches_100[0].val));
                
              } else {
                std::cerr << "Cannot find constant " << vsi2
                          << "in funct_multip_string::operator()."
                          << std::endl;
                O2SCL_ERR2("Cannot find constant in ",
                           "funct_multip_string.",o2scl::exc_efailed);
              }
            }
          }
        }
        
        compiled=true;
      }

      // AWS, 7/1/22: This is a hack to determine the type so we can
      // get the right convert_units object.
      
      int d10=std::numeric_limits<fp_t>::digits10;
      if (verbose>1) {
        std::cout << "funct_multip_string::operator(): input is "
                  << x << " and d10 is " << d10 << std::endl;
      }
      if (d10==15) {
        vars[st_var]=static_cast<double>(x);
        fp_t ret=static_cast<fp_t>(c.eval(&vars));
        if (verbose>1) {
          std::cout << "funct_multip_string::operator(): double "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==18) {
        vars_ld[st_var]=static_cast<long double>(x);
        fp_t ret=static_cast<fp_t>(c_ld.eval(&vars_ld));
        if (verbose>1) {
          std::cout << "funct_multip_string::operator(): long double "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==25) {
        vars_25[st_var]=static_cast<cpp_dec_float_25>(x);
        fp_t ret=static_cast<fp_t>(c_25.eval(&vars_25));
        if (verbose>1) {
          std::cout << "funct_multip_string::operator(): 25-digit "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==35) {
        vars_35[st_var]=static_cast<cpp_dec_float_35>(x);
        fp_t ret=static_cast<fp_t>(c_35.eval(&vars_35));
        if (verbose>1) {
          std::cout << "funct_multip_string::operator(): 35-digit "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==50) {
        vars_50[st_var]=static_cast<cpp_dec_float_50>(x);
        fp_t ret=static_cast<fp_t>(c_50.eval(&vars_50));
        if (verbose>1) {
          std::cout << "funct_multip_string::operator(): 50-digit "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==100) {
        vars_100[st_var]=static_cast<cpp_dec_float_100>(x);
        fp_t ret=static_cast<fp_t>(c_100.eval(&vars_100));
        if (verbose>1) {
          std::cout << "funct_multip_string::operator(): 100-digit "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      }

      O2SCL_ERR("Unexpected type in funct_multip_string.",
                o2scl::exc_einval);
      return o2scl::exc_einval;
    }

  private:

    funct_multip_string(const funct_multip_string &);
    funct_multip_string& operator=(const funct_multip_string&);
    
  };
    
  
}

#endif
