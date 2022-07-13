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
#ifndef O2SCL_LIB_SETTINGS_H
#define O2SCL_LIB_SETTINGS_H
#include <iostream>
#include <string>

#include <o2scl/convert_units.h>
#include <o2scl/find_constants.h>

/** \file lib_settings.h
    \brief Library settings class and global settings object
*/

/** \brief The main \o2 namespace
    
    By default, all \o2 classes and functions which are not listed as
    being in one of \o2's smaller specialized namespaces are in this
    namespace. 

    \comment
    \htmlonly
    For a full list of all the 
    O<span style='position: relative; top: 0.3em; font-size: 0.8em'>2</span>scl 
    classes, see
    <a href="annotated.html">Data Structures</a>.
    \endhtmlonly
    \endcomment
    
    This namespace documentation is in the file 
    <tt>src/base/lib_settings.h</tt>
*/
/*
  \comment
  For a full list of all the
  O<span style='position: relative; top: 0.3em; font-size: 0.8em'>2</span>scl 
  global objects which are not
  classes, see <a href="globals.html">Files::Globals</a>.
  \endcomment
*/
namespace o2scl {
}

namespace o2scl {

  /** \brief A class to manage global library settings

      This class reports global \o2 settings such as the current
      version, whether or not sub-libraries were installed and what
      the current parent directory for \o2 data files is.

      A global object of this type is defined in
      <tt>lib_settings.h</tt> named \ref o2scl_settings .

      This class should not typically be instantiated by the end-user.
  */
  class lib_settings_class {

  public:

    lib_settings_class();

    ~lib_settings_class();

    /// \name Library settings
    //@{
    /** \brief Return the data directory */
    std::string get_data_dir() {
      return data_dir;
    }

    /** \brief Set the data directory */
    int set_data_dir(std::string dir) {
      data_dir=dir;
      return 0;
    }
    
    /** \brief Return the doc directory */
    std::string get_doc_dir() {
      return doc_dir;
    }

    /** \brief Set the doc directory */
    int set_doc_dir(std::string dir) {
      doc_dir=dir;
      return 0;
    }
    
    /// Return true if \o2 was installed with OpenMP support
    bool openmp_support();

    /// Return true if \o2 was installed with readline support
    bool readline_support();

    /// Return true if \o2 was installed with mpfr support
    bool mpfr_support();

    /// Return true if \o2 was installed with ncurses support
    bool ncurses_support();

    /// Return true if \o2 was installed with support for GSL V2.0+
    bool gsl2_support();

    /// Return true if \o2 was installed with Armadillo support
    bool armadillo_support();

    /// Return true if \o2 was installed with Eigen support
    bool eigen_support();

    /// Return true if \o2 was installed with FFTW support
    bool fftw_support();

    /// Return true if \o2 was installed with CUBATURE support
    bool cubature_support();

    /// Return true if \o2 was installed with POLYLOGARITHM support
    bool polylogarithm_support();

    // Return true if \o2 was installed with Python support
    bool python_support();

    /// Return true if \o2 was installed with HDF5 compression support
    bool hdf5_compression_support();

    /** \brief Return system type determined by autoconf

	Returns either "OSX", "Linux" or "unknown".
    */
    std::string system_type();
    
    /** \brief Return true if range checking was turned on during 
	installation (default true)
    */
    bool range_check();

    /// Return the time \o2 was compiled
    std::string time_compiled();

    /// Return the date \o2 was compiled
    std::string date_compiled();

    /// Return the library version
    std::string o2scl_version();

    /// Report some of the settings from config.h
    void config_h_report();

    /// Obtain HDF5 version from the library compiled with O2scl
    void hdf5_lib_version(unsigned &maj, unsigned &min, unsigned &rel);
    
    /// Obtain HDF5 version from the current headers
    void hdf5_header_version(unsigned &maj, unsigned &min, unsigned &rel);
    //@}
    
    /// \name Miscellaneous config.h string properties
    //@{
    std::string o2scl_name();
    std::string o2scl_package();
    std::string o2scl_bugreport();
    std::string o2scl_string();
    std::string o2scl_tarname();
    //@}

    /// \name Unit conversion objects
    //@{
    /// Default convert_units object
    convert_units<double> def_cu;
    
    /// Get the global convert_units object
    convert_units<double> &get_convert_units() {
      return *cup;
    }
    //@}

    // AWS: 2/22/21: I was originally thinking of using this to
    // control OpenMP threads, but I think for now the best is just to
    // use export OMP_NUM_THREADS to control this
    //size_t omp_num_threads;

    /// True if Python has been initialized (default false)
    bool py_initialized;

    /// Initialize the python interface
    int py_init_nothrow(int verbose=0);

    /// Initialize the python interface
    void py_init(int verbose=0);
    
    /// Finalize the python interface
    int py_final_nothrow(int verbose=0);

    /// Finalize the python interface
    void py_final(int verbose=0);

    /// Add path \c path to the python system search path
    void add_python_path(std::string path, int verbose=0);
    
  protected:

#ifndef DOXYGEN_INTERNAL

    /// \name Internal data set in the constructor [protected]
    //@{
    /// The present data directory
    std::string data_dir;

    /// The present documentation directory
    std::string doc_dir;

    /// Pointer to current \ref convert_units object
    convert_units<double> *cup;

    /// Pointer to current \ref find_constants object
    find_constants<> *fcp;
    //@}

#endif
  
  };

  /** \brief The global library settings object

      This global object is used by \o2 classes to store the global
      constant database, the global unit conversion database, and the
      directory for \o2p and \o2e data files. It may also be used by
      the end-user to probe details of the \o2 installation.
  */
  extern lib_settings_class o2scl_settings;
  
  /** \brief Convert a formula to a floating point number and 
      return an integer to indicate success or failure
      
      This is an alternate version of \ref function_to_double()
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
  
    int ret=calc.compile_nothrow(s2.c_str(),0);
    if (ret!=0) return ret;

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
      int ret2=calc.eval_nothrow(0,result);
      if (ret2!=0) return ret2;
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
      O2SCL_ERR("Function function_to_double() failed.",ret);
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

  /** \brief Find constant named \c name with unit \c unit and
      return the associated value
  */
  template<class fp_t=double>
  fp_t find_constant(std::string name, std::string unit) {
    o2scl::convert_units<fp_t> &cu=o2scl_settings.get_convert_units();
    return cu.find_unique(name,unit);
  }

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

      O2SCL_ERR("Unexpected type in funct_multip_strings.",
                o2scl::exc_einval);
      return o2scl::exc_einval;
    }
    
  };
    
}

extern "C" {
  void *o2scl_get_o2scl_settings();
}

#endif
