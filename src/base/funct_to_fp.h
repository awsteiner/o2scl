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
#ifndef O2SCL_FUNCT_TO_FP_H
#define O2SCL_FUNCT_TO_FP_H
#include <iostream>
#include <string>

#include <o2scl/convert_units.h>
#include <o2scl/find_constants.h>
#include <o2scl/string_conv.h>
#include <o2scl/calc_utf8.h>
#include <o2scl/funct_multip.h>

/** \file funct_to_fp.h
    \brief Functions which convert strings to floating point numbers
*/

namespace o2scl {

  /** \brief Array of multi-dimensional functions in an array of strings
   */
  class mm_funct_strings {

  public:
      
    /** \brief Specify the strings
     */
    template<class vec_string_t=std::vector<std::string> >
      mm_funct_strings(int nv, vec_string_t &exprs,
			 vec_string_t &var_arr) {

      st_nv=nv;
      st_forms.resize(nv);
      st_vars.resize(nv);
      calc.resize(nv);
      for (int i=0;i<nv;i++) {
        calc[i]=new calc_utf8<>;
	calc[i]->compile(exprs[i].c_str(),&vars);
	st_vars[i]=var_arr[i];
	st_forms[i]=exprs[i];
      }
    }
      
    virtual ~mm_funct_strings() {
      for (size_t i=0;i<calc.size();i++) {
        delete calc[i];
      }
      calc.clear();
    };
      
    /** \brief Set the values of the auxilliary parameters that were
	specified in 'parms' in the constructor
    */
    int set_parm(std::string name, double val) {
      vars[name]=val;
      return 0;
    }
      
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    template<class vec_t=boost::numeric::ublas::vector<double> >
      int operator()(size_t nv, const vec_t &x, vec_t &y) {

      for(size_t i=0;i<nv;i++) {
	vars[st_vars[i]]=x[i];
      }
      for(size_t i=0;i<nv;i++) {
	y[i]=calc[i]->eval(&vars);
      }
      return 0;
    }
      
    /// Set the functions
    template<class vec_string_t=std::vector<std::string> >
      void set_function(int nv, vec_string_t &exprs,
			vec_string_t &var_arr) {

      st_nv=nv;
      st_forms.resize(nv);
      st_vars.resize(nv);
      calc.resize(nv);
      for (int i=0;i<nv;i++) {
        calc[i]=new calc_utf8<>;
	calc[i]->compile(exprs[i],&vars);
	st_vars[i]=var_arr[i];
	st_forms[i]=exprs[i];
      }

      return;
    }
      
#ifndef DOXYGEN_INTERNAL
      
  protected:
      
    /// The function parsers
    std::vector<calc_utf8<> *> calc;
      
    /// External variables to include in the function parsing
    std::map<std::string,double> vars;
      
    /// The expressions
    std::vector<std::string> st_forms;
      
    /// The variables
    std::vector<std::string> st_vars;
      
    /// The number of variables
    int st_nv;
      
    mm_funct_strings() {};
      
  private:
      
    mm_funct_strings(const mm_funct_strings &);
    mm_funct_strings& operator=(const mm_funct_strings&);
      
#endif
      
  };

  /** \brief A multi-dimensional function from a string
   */
  class multi_funct_strings {
    
  public:
  
    /** \brief Specify the string and the parameters
     */
    template<class vec_string_t=std::vector<std::string> >
    multi_funct_strings(std::string expr, int nv,
                        vec_string_t &var_arr) {
    
      st_nv=nv;
      st_funct=expr;
      st_vars.resize(nv);
      for (int i=0;i<nv;i++) {
	calc.compile(expr.c_str(),&vars);
	st_vars[i]=var_arr[i];
      }
    }
  
    /** \brief Specify the string and the parameters
     */
    template<class vec_string_t=std::vector<std::string> >
    void set_function(std::string expr, int nv, vec_string_t &var_arr) {

      st_nv=nv;
      st_funct=expr;
      st_vars.resize(nv);
      for (int i=0;i<nv;i++) {
	calc.compile(expr.c_str(),&vars);
	st_vars[i]=var_arr[i];
      }
      return;
    }

    virtual ~multi_funct_strings() {
    };
  
    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms  in the constructor
    */
    int set_parm(std::string name, double val) {
      vars[name]=val;
      return 0;
    }

    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
    template<class vec_t=boost::numeric::ublas::vector<double> >
    double operator()(size_t nv, const vec_t &x) {

      for(size_t i=0;i<nv;i++) {
	vars[st_vars[i]]=x[i];
      }

      return calc.eval(&vars);
    }

  protected:

    /// The function parser
    calc_utf8<> calc;

    /// External variables to include in the function parsing
    std::map<std::string,double> vars;

    /// The number of variables
    int st_nv;

    /// The function string
    std::string st_funct;
    
    /// The variable string
    std::vector<std::string> st_vars;
  
    multi_funct_strings() {}
  
  private:

    multi_funct_strings(const multi_funct_strings &);
    multi_funct_strings& operator=(const multi_funct_strings&);

  };

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

#ifndef O2SCL_NO_BOOST_MULTIPRECISION

  /** \brief Evaluate a one-dimensional function from a string
      at multiprecision

      \note Experimental.

      \warning This class only supports a limited number of data
      types, including double, long double, and cpp_dec_float types
      with 25, 35, 50, or 100 digits. It is designed to be used with
      the \ref funct_multip class.
   */
  template<class fp_25_t, class fp_35_t, class fp_50_t, class fp_100_t>
  class funct_multip_string_tl {

  protected:
    
    /// \name The function evaluation objects
    //@{
    calc_utf8<double> c;
    calc_utf8<long double> c_ld;
    calc_utf8<fp_25_t> c_25;
    calc_utf8<fp_35_t> c_35;
    calc_utf8<fp_50_t> c_50;
    calc_utf8<fp_100_t> c_100;
    //@}

    /// \name The unit conversion objects
    //@{
    convert_units<double> cu;
    convert_units<long double> cu_ld;
    convert_units<fp_25_t> cu_25;
    convert_units<fp_35_t> cu_35;
    convert_units<fp_50_t> cu_50;
    convert_units<fp_100_t> cu_100;
    //@}

    /// \name The variable lists
    //@{
    std::map<std::string,double> vars;
    std::map<std::string,long double> vars_ld;
    std::map<std::string,fp_25_t> vars_25;
    std::map<std::string,fp_35_t> vars_35;
    std::map<std::string,fp_50_t> vars_50;
    std::map<std::string,fp_100_t> vars_100;
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
    
    funct_multip_string_tl() {
      verbose=0;
      err_nonconv=true;
      compiled=false;
    }

    virtual ~funct_multip_string_tl() {
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
          std::cout << "funct_multip_string_tl::operator() "
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
                std::cout << "funct_multip_string_tl::operator() "
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
                            find_constants<fp_25_t>::const_entry>
                  matches_25;
                cu_25.find_nothrow(vsi2,"mks",matches_25);
                vars_25.insert(std::make_pair(vsi2,matches_25[0].val));
                
                std::vector<typename
                            find_constants<fp_35_t>::const_entry>
                  matches_35;
                cu_35.find_nothrow(vsi2,"mks",matches_35);
                vars_35.insert(std::make_pair(vsi2,matches_35[0].val));
                
                std::vector<typename
                            find_constants<fp_50_t>::const_entry>
                  matches_50;
                cu_50.find_nothrow(vsi2,"mks",matches_50);
                vars_50.insert(std::make_pair(vsi2,matches_50[0].val));
                
                std::vector<typename
                            find_constants<fp_100_t>::const_entry>
                  matches_100;
                cu_100.find_nothrow(vsi2,"mks",matches_100);
                vars_100.insert(std::make_pair(vsi2,matches_100[0].val));
                
              } else {
                std::cerr << "Cannot find constant " << vsi2
                          << "in funct_multip_string_tl::operator()."
                          << std::endl;
                O2SCL_ERR2("Cannot find constant in ",
                           "funct_multip_string_tl.",o2scl::exc_efailed);
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
        std::cout << "funct_multip_string_tl::operator(): input is "
                  << x << " and d10 is " << d10 << std::endl;
      }
      if (d10==15) {
        vars[st_var]=static_cast<double>(x);
        fp_t ret=static_cast<fp_t>(c.eval(&vars));
        if (verbose>1) {
          std::cout << "funct_multip_string_tl::operator(): double "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==18) {
        vars_ld[st_var]=static_cast<long double>(x);
        fp_t ret=static_cast<fp_t>(c_ld.eval(&vars_ld));
        if (verbose>1) {
          std::cout << "funct_multip_string_tl::operator(): long double "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==25) {
        vars_25[st_var]=static_cast<fp_25_t>(x);
        fp_t ret=static_cast<fp_t>(c_25.eval(&vars_25));
        if (verbose>1) {
          std::cout << "funct_multip_string_tl::operator(): 25-digit "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==35) {
        vars_35[st_var]=static_cast<fp_35_t>(x);
        fp_t ret=static_cast<fp_t>(c_35.eval(&vars_35));
        if (verbose>1) {
          std::cout << "funct_multip_string_tl::operator(): 35-digit "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==50) {
        vars_50[st_var]=static_cast<fp_50_t>(x);
        fp_t ret=static_cast<fp_t>(c_50.eval(&vars_50));
        if (verbose>1) {
          std::cout << "funct_multip_string_tl::operator(): 50-digit "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      } else if (d10==100) {
        vars_100[st_var]=static_cast<fp_100_t>(x);
        fp_t ret=static_cast<fp_t>(c_100.eval(&vars_100));
        if (verbose>1) {
          std::cout << "funct_multip_string_tl::operator(): 100-digit "
                    << "precision returning " << ret << std::endl;
        }
        return ret;
      }

      O2SCL_ERR("Unexpected type in funct_multip_string_tl.",
                o2scl::exc_einval);
      return o2scl::exc_einval;
    }

  private:

    funct_multip_string_tl(const funct_multip_string_tl &);
    funct_multip_string_tl& operator=(const funct_multip_string_tl&);
    
  };
    
  /** \brief Typedef for default floating point type
   */
  typedef funct_multip_string_tl<o2fp_25,o2fp_35,o2fp_50,o2fp_100>
  funct_multip_string;

  /** \brief Typedef for cpp_dec_float types
   */
  typedef funct_multip_string_tl<cpp_dec_float_25,cpp_dec_float_35,
                          cpp_dec_float_50,cpp_dec_float_100>
  funct_multip_string_cdf;

#ifdef O2SCL_SET_MPFR
  /** \brief Typedef for mpfr types
   */
  typedef funct_multip_string_tl<mpfr_25,mpfr_35,mpfr_50,mpfr_100>
  funct_multip_string_mpfr;
#endif

#endif
  
  /** \brief Convert a formula to a floating point number and 
      return an integer to indicate success or failure
      
      This is an alternate version of \ref function_to_fp()
      which does not call the error handler and returns a non-zero
      integer when it fails.
  */
  template<class fp_t=double>
  int function_to_fp_nothrow(std::string s, fp_t &result,
                             convert_units<fp_t> &cu,
                             int verbose=0, rng<> *r=0) {
    
    std::string s2;
    // Remove quotes and apostrophes
    for(size_t i=0;i<s.length();i++) {
      if (s[i]!='\"' && s[i]!='\'') {
        s2+=s[i];
      }
    }
    
    calc_utf8<fp_t> calc;
    if (r!=0) {
      calc.set_rng(*r);
    }

    if (verbose>1) {
      calc.verbose=verbose;
    }
    int ret=calc.compile_nothrow(s2.c_str(),0);
    if (ret!=0) return ret;

    if (verbose>=2) {
      std::cout << calc.RPN_to_string() << std::endl;
    }

    std::vector<std::u32string> vs=calc.get_var_list();

    // If there are undefined variables, then attempt to get them
    // from the constant database
    if (vs.size()!=0) {
    
      find_constants<fp_t> &fc=cu.fc;
    
      std::map<std::string,fp_t> vars;
      std::vector<typename find_constants<fp_t>::const_entry> matches;
      
      for(size_t i=0;i<vs.size();i++) {

        std::string vsi2;
        char32_to_utf8(vs[i],vsi2);

        if (verbose>2) {
          std::cout << "Function function_to_fp_nothrow(): "
                    << "trying to find constant " << vsi2 << std::endl;
        }

        cu.verbose=verbose;
        int fret=cu.find_nothrow(vsi2,"mks",matches);
      
        if (fret==find_constants<fp_t>::one_exact_match_unit_match ||
            fret==find_constants<fp_t>::one_pattern_match_unit_match) {

          typename find_constants<fp_t>::const_entry &fcl=matches[0];
          
          vars.insert(std::make_pair(vsi2,fcl.val));
          if (verbose>1) {
            std::cout << "Function function_to_fp_nothrow(): "
                      << "Found constant " << vsi2
                      << " value " << fcl.val << std::endl;
          }
        
        } else {
        
          if (verbose>=2) {
            std::cout << "Variable " << vsi2
                      << " not uniquely specified in constant list ("
                      << fret << ")." << std::endl;
          }
        
          return 1;
        }
      }

      // Evaluate the expression with the variables assigned above
      int ret2=calc.eval_nothrow(&vars,result);
      if (ret2!=0) return ret2;
    
    } else {

      // Evaluate the expression (no variables necessary)
      if (verbose>=3) {
        std::cout << "In funct_to_fp(), calling calc.eval_nothrow()."
                  << std::endl;
      }
      int ret2=calc.eval_nothrow(0,result);
      if (ret2!=0) return ret2;
      if (verbose>=3) {
        std::cout << "In funct_to_fp(), done with calc.eval_nothrow()."
                  << std::endl;
      }
    }
  
    return 0;
  }

  /** \brief Convert a formula to a double 
      
      This function removes all quotes and apostrophes from the string
      and then uses \ref o2scl::calculator to convert strings like
      "-1.0e-3", "pi/3.0" and "exp(cos(-1.0e-2))" to floating point
      numbers. This function uses the \o2 constant database from
      \ref lib_settings_class::get_find_constants() to interpret
      constant values.
  */
  template<class fp_t=double>
  fp_t function_to_fp(std::string s, int verbose=0) {
    fp_t res;
    convert_units<fp_t> cu;
    int ret=function_to_fp_nothrow<fp_t>(s,res,cu,verbose);
    if (ret!=0) {
      O2SCL_ERR("Function function_to_fp() failed.",ret);
    }
    return res;
  }

  /** \brief Convert a formula to a double 
      
      This function removes all quotes and apostrophes from the string
      and then uses \ref o2scl::calculator to convert strings like
      "-1.0e-3", "pi/3.0" and "exp(cos(-1.0e-2))" to floating point
      numbers. This function uses the \o2 constant database from
      \ref lib_settings_class::get_find_constants() to interpret
      constant values.
  */
  int function_to_double_nothrow(std::string s, double &result,
                                 int verbose=0, rng<> *r=0);
    
  /** \brief Convert a formula to a double 
      
      This function removes all quotes and apostrophes from the string
      and then uses \ref o2scl::calculator to convert strings like
      "-1.0e-3", "pi/3.0" and "exp(cos(-1.0e-2))" to floating point
      numbers. This function uses the \o2 constant database from
      \ref lib_settings_class::get_find_constants() to interpret
      constant values.
  */
  double function_to_double(std::string s, int verbose=0);

}

#endif
