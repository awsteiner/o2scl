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

  /** \brief A wrapper to provide a template operator() function
      from a set of type-specific function objects
   */
  class funct_multip_wrapper {

  protected:
    
    funct f;
    funct_ld f_ld;
    funct_cdf25 f_cdf25;
    funct_cdf35 f_cdf35;
    funct_cdf50 f_cdf50;
    funct_cdf100 f_cdf100;
  
    /** \brief The function for double
     */
    double eval(double x) {
      return f(x);
    }
  
    /** \brief The function for long double
     */
    long double eval(long double x) {
      return f_ld(x);
    }

    /** \brief The function for 25-digit numbers
     */
    boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<25>> eval
    (boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<25>> x) {
      return f_cdf25(x);
    }
  
    /** \brief The function for 35-digit numbers
     */
    boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<35>> eval
    (boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<35>> x) {
      return f_cdf35(x);
    }
  
    /** \brief The function for 50-digit numbers
     */
    boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<50>> eval
    (boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<50>> x) {
      return f_cdf50(x);
    }
  
    /** \brief The function for 100-digit numbers
     */
    boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<100>> eval
    (boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<100>> x) {
      return f_cdf100(x);
    }
  
  public:

    /** \brief Create a wrapper from the given function objects
    */
    funct_multip_wrapper(funct &f_x,
                         funct_ld &f_ld_x,
                         funct_cdf25 &f_cdf25_x,
                         funct_cdf35 &f_cdf35_x,
                         funct_cdf50 &f_cdf50_x,
                         funct_cdf100 &f_cdf100_x) {
      f=f_x;
      f_ld=f_ld_x;
      f_cdf25=f_cdf25_x;
      f_cdf35=f_cdf35_x;
      f_cdf50=f_cdf50_x;
      f_cdf100=f_cdf100_x;
      return;
    }

    /** \brief The template operator() function
     */
    template<class fp_t> fp_t operator()(fp_t x) {
      return eval(x);
    }

  };

  /** \brief Use multiprecision to automatically evaluate a function to
      a specified level of precision
  */
  template<class func_t=funct_multip_wrapper> class funct_multip {
  
  protected:

    /** \brief Desc
     */
    func_t &f;

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

    /** \brief Create a function evaluator based on the operator()
        given in the class \c fx of type \c func_t
     */
    funct_multip(func_t &fx) : f(fx) {
      verbose=0;
      tol_rel=-1.0;
    }

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief Evalulate the function with the specified tolerance
        and provide an error estimate
     */
    template<class fp_t> int eval_tol_err(fp_t x, fp_t &val,
                                          fp_t &err, double tol_loc=-1) {
    
      if (tol_loc<=0.0 && tol_rel<=0.0) {
        tol_loc=pow(10.0,-std::numeric_limits<fp_t>::digits10);
      } else {
        tol_loc=tol_rel;
      }
      if (verbose>0) {
        std::cout << "Set tol to: " << tol_loc << std::endl;
      }
    
      double x_d=static_cast<double>(x);
      long double x_ld=static_cast<long double>(x);
      double y_d=f(x_d);
      long double y_ld=f(x_ld);

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
      if (verbose>0) {
        std::cout << "Failed 1: " << dtos(y_ld,0) << " "
                  << dtos(y_d,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      cpp_dec_float_25 x_cdf25=static_cast<cpp_dec_float_25>(x);
      cpp_dec_float_25 y_cdf25=f(x_cdf25);

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
      if (verbose>0) {
        std::cout << "Failed 2: " << dtos(y_cdf25,0) << " "
                  << dtos(y_ld,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      cpp_dec_float_35 x_cdf35=static_cast<cpp_dec_float_35>(x);
      cpp_dec_float_35 y_cdf35=f(x_cdf35);
        
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
      if (verbose>0) {
        std::cout << "Failed 3: " << dtos(y_cdf35,0) << " "
                  << dtos(y_cdf25,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      cpp_dec_float_50 x_cdf50=static_cast<cpp_dec_float_50>(x);
      cpp_dec_float_50 y_cdf50=f(x_cdf50);
    
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
      if (verbose>0) {
        std::cout << "Failed 4: " << dtos(y_cdf50,0) << " "
                  << dtos(y_cdf35,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      cpp_dec_float_100 x_cdf100=static_cast<cpp_dec_float_100>(x);
      cpp_dec_float_100 y_cdf100=f(x_cdf100);
    
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
      if (verbose>0) {
        std::cout << "Failed 5: " << dtos(y_cdf100,0) << " "
                  << dtos(y_cdf50,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in funct_multip::eval_tol_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

    /** \brief Evaluate the function and return the error estimate
        with the default tolerance for the specified type
     */
    template<class fp_t> int eval_err(fp_t x, fp_t &val, fp_t &err) {
      return eval_tol_err(x,val,err);
    }
  
    /** \brief Evalulate the function without an error estimate
     */
    template<class fp_t> fp_t operator()(fp_t x) {
      fp_t val, err;
      
      eval_err(x,val,err);
    
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

  /** \brief A wrapper to provide a template operator() function
      from a set of type-specific function objects
   */
  class funct_ret_multip_wrapper {

  protected:
    
    funct_ret f;
    funct_ret_ld f_ld;
    funct_ret_cdf25 f_cdf25;
    funct_ret_cdf35 f_cdf35;
    funct_ret_cdf50 f_cdf50;
    funct_ret_cdf100 f_cdf100;
  
    /** \brief The function for double
     */
    int eval(double x, double &y) {
      return f(x,y);
    }
  
    /** \brief The function for long double
     */
    int eval(long double x, long double &y) {
      return f_ld(x,y);
    }

    /** \brief The function for 25-digit numbers
     */
    int eval
    (boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<25>> x,
     boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<25>> &y) {
      return f_cdf25(x,y);
    }
  
    /** \brief The function for 35-digit numbers
     */
    int eval
    (boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<35>> x,
     boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<35>> &y) {
      return f_cdf35(x,y);
    }
    
    /** \brief The function for 50-digit numbers
     */
    int eval
    (boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<50>> x,
     boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<50>> &y) {
      return f_cdf50(x,y);
    }
  
    /** \brief The function for 100-digit numbers
     */
    int eval
    (boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<100>> x,
     boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<100>> &y) {
      return f_cdf100(x,y);
    }
  
  public:

    /** \brief Create a wrapper from the given function objects
    */
    funct_ret_multip_wrapper(funct_ret &f_x,
                             funct_ret_ld &f_ld_x,
                             funct_ret_cdf25 &f_cdf25_x,
                             funct_ret_cdf35 &f_cdf35_x,
                             funct_ret_cdf50 &f_cdf50_x,
                             funct_ret_cdf100 &f_cdf100_x) {
      f=f_x;
      f_ld=f_ld_x;
      f_cdf25=f_cdf25_x;
      f_cdf35=f_cdf35_x;
      f_cdf50=f_cdf50_x;
      f_cdf100=f_cdf100_x;
      return;
    }

    /** \brief The template operator() function
     */
    template<class fp_t> int operator()(fp_t x, fp_t &y) {
      return eval(x,y);
    }

  };

  /** \brief Use multiprecision to automatically evaluate a function to
      a specified level of precision
  */
  template<class func_t=funct_ret_multip_wrapper> class funct_ret_multip {
  
  protected:

    /** \brief Desc
     */
    func_t &f;

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

    /** \brief Create a function evaluator based on the operator()
        given in the class \c fx of type \c func_t
     */
    funct_ret_multip(func_t &fx) : f(fx) {
      verbose=0;
      tol_rel=-1.0;
    }

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief Evalulate the function with the specified tolerance
        and provide an error estimate
     */
    template<class fp_t> int eval_tol_err(fp_t x, fp_t &val,
                                          fp_t &err, double tol_loc=-1.0) {
      
      if (tol_loc<=0.0 && tol_rel<=0.0) {
        tol_loc=pow(10.0,-std::numeric_limits<fp_t>::digits10);
      } else {
        tol_loc=tol_rel;
      }
      if (verbose>0) {
        std::cout << "Set tol to: " << tol_loc << std::endl;
      }
    
      double x_d=static_cast<double>(x), y_d;
      long double x_ld=static_cast<long double>(x), y_ld;
      int r_d, r_ld;
      r_d=f(x_d,y_d);
      r_ld=f(x_ld,y_ld);

      if (r_d==0 && r_ld==0 && y_ld==0 && y_d==0) {
        val=0;
        err=0;
        return 0;
      }
      if (r_d==0 && r_ld==0 && y_ld!=0) {
        err=static_cast<fp_t>(abs(y_ld-y_d)/abs(y_ld));
        if (err<tol_loc) {
          val=static_cast<fp_t>(y_ld);
          return 0;
        }
      }
      if (verbose>0) {
        std::cout << "Failed 1: " << r_ld << " " << r_d << " "
                  << dtos(y_ld,0) << " " << dtos(y_d,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      cpp_dec_float_25 x_cdf25=static_cast<cpp_dec_float_25>(x);
      cpp_dec_float_25 y_cdf25;
      int r_cdf25=f(x_cdf25,y_cdf25);

      if (r_cdf25==0 && r_ld==0 && y_cdf25==0 && y_ld==0) {
        val=0;
        err=0;
        return 0;
      }
      if (r_cdf25==0 && r_ld==0 && y_cdf25!=0) {
        err=static_cast<fp_t>(abs(y_cdf25-y_ld)/abs(y_cdf25));
        if (err<tol_loc) {
          val=static_cast<fp_t>(y_cdf25);
          return 0;
        }
      }
      if (verbose>0) {
        std::cout << "Failed 2: " << r_cdf25 << " " << r_ld << " "
                  << dtos(y_cdf25,0) << " " << dtos(y_ld,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      cpp_dec_float_35 x_cdf35=static_cast<cpp_dec_float_35>(x);
      cpp_dec_float_35 y_cdf35;
      int r_cdf35=f(x_cdf35,y_cdf35);
        
      if (r_cdf35==0 && r_cdf25==0 && y_cdf35==0 && y_cdf25==0) {
        val=0;
        err=0;
        return 0;
      }
      if (r_cdf35==0 && r_cdf25==0 && y_cdf35!=0) {
        err=static_cast<fp_t>(abs(y_cdf35-y_cdf25)/abs(y_cdf35));
        if (err<tol_loc) {
          val=static_cast<fp_t>(y_cdf35);
          return 0;
        }
      }
      if (verbose>0) {
        std::cout << "Failed 3: " << r_cdf35 << " " << r_cdf25 << " "
                  << dtos(y_cdf35,0) << " " << dtos(y_cdf25,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      cpp_dec_float_50 x_cdf50=static_cast<cpp_dec_float_50>(x);
      cpp_dec_float_50 y_cdf50;
      int r_cdf50=f(x_cdf50,y_cdf50);
    
      if (r_cdf50==0 && r_cdf35==0 && y_cdf50==0 && y_cdf35==0) {
        val=0;
        err=0;
        return 0;
      }
      if (r_cdf50==0 && r_cdf35==0 && y_cdf50!=0) {
        err=static_cast<fp_t>(abs(y_cdf50-y_cdf35)/abs(y_cdf50));
        if (err<tol_loc) {
          val=static_cast<fp_t>(y_cdf50);
          return 0;
        }
      }
      if (verbose>0) {
        std::cout << "Failed 4: " << r_cdf50 << " " << r_cdf35 << " "
                  << dtos(y_cdf50,0) << " " << dtos(y_cdf35,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      cpp_dec_float_100 x_cdf100=static_cast<cpp_dec_float_100>(x);
      cpp_dec_float_100 y_cdf100;
      int r_cdf100=f(x_cdf100,y_cdf100);
    
      if (r_cdf100==0 && r_cdf50==0 && y_cdf100==0 && y_cdf50==0) {
        val=0;
        err=0;
        return 0;
      }
      if (r_cdf100==0 && r_cdf50==0 && y_cdf100!=0) {
        err=static_cast<fp_t>(abs(y_cdf100-y_cdf50)/abs(y_cdf100));
        if (err<tol_loc) {
          val=static_cast<fp_t>(y_cdf100);
          return 0;
        }
      }
      if (verbose>0) {
        std::cout << "Failed 5: " << r_cdf100 << " " << r_cdf50 << " "
                  << dtos(y_cdf100,0) << " " << dtos(y_cdf50,0) << " "
                  << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in funct_ret_multip::eval_tol_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

    /** \brief Evaluate the function and return the error estimate
        with the default tolerance for the specified type
     */
    template<class fp_t> int eval_err(fp_t x, fp_t &val, fp_t &err) {
      return eval_tol_err(x,val,err);
    }
  
    /** \brief Evalulate the function without an error estimate
     */
    template<class fp_t> fp_t operator()(fp_t x) {
      fp_t val,err;
      
      eval_err(x,val,err);
    
      return val;
    }
    
  };

  
}

#endif
