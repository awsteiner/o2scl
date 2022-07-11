/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_FUNCT_H
#define O2SCL_FUNCT_H

/** \file funct.h
    \brief Function object classes for one-dimensional functions
*/

#include <string>
#include <functional>

#include <gsl/gsl_math.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#ifdef O2SCL_MPFR
#include <boost/multiprecision/mpfr.hpp>
#endif

#include <o2scl/err_hnd.h>
#include <o2scl/calc_utf8.h>
#include <o2scl/lib_settings.h>

#ifdef O2SCL_PYTHON
#include <Python.h>
#endif

namespace o2scl {

  // Experimental template typedef
  //template<class fp_t> using funct2 = std::function<fp_t(fp_t)>;
  
  /// One-dimensional function typedef in src/base/funct.h
  typedef std::function<double(double)> funct;

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
  
  /** \brief One-dimensional function from a string
      
      For example,
      \code
      funct_string f("pi*r^2","r");
      f.set_parm("pi",o2scl_const::pi);
      for(double r=1.0;r<=2.0;r+=0.1) {
      cout << f(x) << endl;
      }
      \endcode
      will print out the area of circles having radii between 1 and 2.
  */
  class funct_string {
    
  public:
    
    /** \brief Specify the string and the parameters
     */
    funct_string(std::string expr, std::string var) {
      calc.compile(expr.c_str(),&vars);
      st_form=expr;
      st_var=var;
    }

    virtual ~funct_string() {
    };

  
    /** \brief Specify the string and the parameters
     */
    int set_function(std::string expr, std::string var) {
      calc.compile(expr.c_str(),&vars);
      st_form=expr;
      st_var=var;
      return 0;
    }

    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms in the constructor
    */
    int set_parm(std::string name, double val) {
      if (name==st_var) {
	O2SCL_ERR2("A parameter cannot have the same name as ",
		   "the variable in funct_string::set_parm().",
		   o2scl::exc_einval);
      }
      vars[name]=val;
      return 0;
    }
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      vars[st_var]=x;
      return calc.eval(&vars);
    }

  protected:

    /// The object for evaluating strings
    mutable o2scl::calc_utf8<> calc;

    /// Parameter map
    mutable std::map<std::string,double> vars;
    
    /// The expr
    std::string st_form;
    /// The variable
    std::string st_var;

    funct_string() {};

  private:

    funct_string(const funct_string &);
    funct_string& operator=(const funct_string&);

  };

  /** \brief A wrapper to specify \ref o2scl::funct objects 
      to GSL
  */
  class funct_gsl : public gsl_function {

  protected:
    
    /// The function wrapper
    static double funct_wrap(double x, void *params) {
      funct *fp=(funct *)params;
      return (*fp)(x);
    }

  public:

    /// Create an object based on the specified function, \c f
    funct_gsl(funct &f) {
      function=&funct_wrap;
      params=&f;
    }
    
  };

  /** \brief Two-dimensional function from a string
   */
  class funct2_string {
    
  public:
    
    /** \brief Specify the string and the parameters
     */
    funct2_string(std::string expr, std::string var1, std::string var2) {
      calc.compile(expr.c_str(),&vars);
      st_form=expr;
      st_var1=var1;
      st_var2=var2;
    }

    virtual ~funct2_string() {
    };

  
    /** \brief Specify the string and the parameters
     */
    int set_function(std::string expr, std::string var1,
		     std::string var2) {
      calc.compile(expr.c_str(),&vars);
      st_form=expr;
      st_var1=var1;
      st_var2=var2;
      return 0;
    }

    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms in the constructor
    */
    int set_parm(std::string name, double val) {
      if (name==st_var1 || name==st_var2) {
	O2SCL_ERR2("A parameter cannot have the same name as ",
		   "a variable in funct_string::set_parm().",
		   o2scl::exc_einval);
      }
      vars[name]=val;
      return 0;
    }
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x, double y) const {
      vars[st_var1]=x;
      vars[st_var2]=y;
      return calc.eval(&vars);
    }

  protected:

    /// The object for evaluating strings
    mutable o2scl::calc_utf8<> calc;

    /// Parameter map
    mutable std::map<std::string,double> vars;
    
    /// The expr
    std::string st_form;
    /// The variable
    std::string st_var1;
    /// The variable
    std::string st_var2;

    funct2_string() {};

  private:

    funct2_string(const funct2_string &);
    funct2_string& operator=(const funct2_string&);

  };

#ifdef O2SCL_PYTHON
  
  /** \brief One-dimensional function from a python function
   */
  class funct_python {

  protected:

    /// Python unicode object containing function name
    PyObject *pName;
    
    /// Python module containing function
    PyObject *pModule;
    
    /// Function arguments
    PyObject *pArgs;

    /// Python function
    PyObject *pFunc;

    /// Verbosity parameter
    int verbose;
    
  public:
    
    /** \brief Specify the python and the parameters
     */
    funct_python(std::string module, std::string func, int v=0);
    
    virtual ~funct_python();
  
    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string func);
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const;

  protected:

    funct_python() {};

  private:

    funct_python(const funct_python &);
    funct_python& operator=(const funct_python&);

  };

  /** \brief One-dimensional function from a python function
   */
  class funct_python_method {

  protected:

    /// Python unicode object containing function name
    PyObject *pName;
    
    /// Python module containing function
    PyObject *pModule;
    
    /// Function arguments
    PyObject *pArgs;

    /// Python function
    PyObject *pFunc;

    /// Python class instance
    PyObject *pInstance;

    /// Python class
    PyObject *pClass;

    /// Verbosity parameter
    int verbose;
    
  public:
    
    /** \brief Specify the python and the parameters
     */
    funct_python_method(std::string module, std::string class_name,
                        std::string func, int v=0);
    
    virtual ~funct_python_method();
  
    /** \brief Specify the python and the parameters

        This function is called by the constructor and thus
        cannot be virtual.
    */
    int set_function(std::string module, std::string class_name,
                     std::string func);
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const;

  protected:

    funct_python_method() {};

  private:

    funct_python_method(const funct_python_method &);
    funct_python_method& operator=(const funct_python_method&);

  };

#endif

  /** \brief Use multiprecision to automatically evaluate a function to
      a specified level of precision

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

      \note The algorithm attempts not to be wasteful, but is not
      necessarily optimized for speed. One way to improve it would
      be to more intelligently choose the number of digits used
      in the boost multiprecision numbers based on the tolerance
      which was specified.
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
      if (tol_loc>pow(10.0,-std::numeric_limits<double>::max_digits10)) {
        x_d=static_cast<double>(x);
        y_d=f(x_d);
        d_eval=true;
      }
      
      bool ld_eval=false;
      long double x_ld=0, y_ld=0;
      if (tol_loc>pow(10.0,-std::numeric_limits<long double>::max_digits10)) {
        x_ld=static_cast<long double>(x);
        y_ld=f(x_ld);
        ld_eval=true;
      }

      if (d_eval && ld_eval) {
        if (y_ld==0 && y_d==0) {
          val=0;
          err=0;
          return 0;
        }
      
        if (y_ld!=0) {
          err=static_cast<fp_t>(abs(y_ld-y_d)/abs(y_ld));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_ld);
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
        } else if (ld_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed first round (d_eval is false): "
                    << dtos(y_ld,0) << " "
                    << tol_loc << std::endl;
        } else {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed first round "
                    << "(d_eval and ld_eval are both false)." << std::endl;
        }
      }
    
      /// Second pass, compare long double and 25-digit precision
      
      bool cdf25_eval=false;
      cpp_dec_float_25 x_cdf25=0, y_cdf25=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_25>::max_digits10)) {
        x_cdf25=static_cast<cpp_dec_float_25>(x);
        y_cdf25=f(x_cdf25);
        cdf25_eval=true;
      }

      if (ld_eval && cdf25_eval) {
        if (y_cdf25==0 && y_ld==0) {
          val=0;
          err=0;
          return 0;
        }
        if (y_cdf25!=0) {
          err=static_cast<fp_t>(abs(y_cdf25-y_ld)/abs(y_cdf25));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf25);
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
        } else if (cdf25_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed second round (ld_eval is false): "
                    << dtos(y_cdf25,0) << " "
                    << tol_loc << std::endl;
        } else {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed second round (ld_eval and "
                    << "cdf25_eval are both false)." << std::endl;
        }
      }
    
      /// Third pass, compare 25- and 35-digit precision

      bool cdf35_eval=false;
      cpp_dec_float_35 x_cdf35=0, y_cdf35=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_35>::max_digits10)) {
        x_cdf35=static_cast<cpp_dec_float_35>(x);
        y_cdf35=f(x_cdf35);
        cdf35_eval=true;
      }

      if (cdf25_eval && cdf35_eval) {
        if (y_cdf35==0 && y_cdf25==0) {
          val=0;
          err=0;
          return 0;
        }
        if (y_cdf35!=0) {
          err=static_cast<fp_t>(abs(y_cdf35-y_cdf25)/abs(y_cdf35));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf35);
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
        } else if (cdf35_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval is false): "
                    << dtos(y_cdf35,0) << " "
                    << tol_loc << std::endl;
        } else {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed third round (cdf25_eval and "
                    << "cdf35_eval are both false)." << std::endl;
        }
      }
    
      /// Fourth pass, compare 35- and 50-digit precision
      
      bool cdf50_eval=false;
      cpp_dec_float_50 x_cdf50=0, y_cdf50=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_50>::max_digits10)) {
        x_cdf50=static_cast<cpp_dec_float_50>(x);
        y_cdf50=f(x_cdf50);
        cdf50_eval=true;
      }

      if (cdf35_eval && cdf50_eval) {
        if (y_cdf50==0 && y_cdf35==0) {
          val=0;
          err=0;
          return 0;
        }
        if (y_cdf50!=0) {
          err=static_cast<fp_t>(abs(y_cdf50-y_cdf35)/abs(y_cdf50));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf50);
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
        } else if (cdf50_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval is false): "
                    << dtos(y_cdf50,0) << " "
                    << tol_loc << std::endl;
        } else {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed fourth round (cdf35_eval and "
                    << "cdf50_eval are both false)." << std::endl;
        }
      }
    
      /// Final pass, compare 50- and 100-digit precision
      
      bool cdf100_eval=false;
      cpp_dec_float_100 x_cdf100=0, y_cdf100=0;
      if (tol_loc>pow(10.0,-std::numeric_limits
                      <cpp_dec_float_100>::max_digits10)) {
        x_cdf100=static_cast<cpp_dec_float_100>(x);
        y_cdf100=f(x_cdf100);
        cdf100_eval=true;
      }

      if (cdf100_eval && cdf50_eval) {
        if (y_cdf100==0 && y_cdf50==0) {
          val=0;
          err=0;
          return 0;
        }
        if (y_cdf100!=0) {
          err=static_cast<fp_t>(abs(y_cdf100-y_cdf50)/abs(y_cdf100));
          if (err<tol_loc) {
            val=static_cast<fp_t>(y_cdf100);
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
        } else if (cdf100_eval) {
          std::cout << "funct_multip::eval_tol_err():\n  "
                    << "Failed last round (cdf50_eval is false): "
                    << dtos(y_cdf100,0) << " "
                    << tol_loc << std::endl;
        } else {
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

  /// One-dimensional function typedef in src/base/funct.h
  typedef std::function<int(double,double &)> funct_ret;

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

}

#endif
