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
#ifndef O2SCL_FUNCT_TO_FP_H
#define O2SCL_FUNCT_TO_FP_H
#include <iostream>
#include <string>

#include <o2scl/convert_units.h>
#include <o2scl/find_constants.h>

/** \file funct_to_fp.h
    \brief Desc
*/

namespace o2scl {

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

}

#endif
