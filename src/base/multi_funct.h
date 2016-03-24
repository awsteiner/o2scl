/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#ifndef O2SCL_MULTI_FUNCT_H
#define O2SCL_MULTI_FUNCT_H

/** \file multi_funct.h
    \brief Function object classes for a multi-dimensional function
*/

#include <string>
#include <functional>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/err_hnd.h>
#include <o2scl/shunting_yard.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// Multi-dimensional function typedef
  typedef std::function<
    double(size_t,const boost::numeric::ublas::vector<double> &)>
    multi_funct11;
  
  /** \brief A multi-dimensional function from a string
   */
  class multi_funct11_strings {
    
  public:
  
    /** \brief Specify the string and the parameters
     */
    template<class vec_string_t=std::vector<std::string> >
      multi_funct11_strings(std::string expr, int nv,
			    vec_string_t &var_arr) {
    
      st_nv=nv;
      st_form=expr;
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
      st_form=expr;
      st_vars.resize(nv);
      for (int i=0;i<nv;i++) {
	calc.compile(expr.c_str(),&vars);
	st_vars[i]=var_arr[i];
      }
      return;
    }

    virtual ~multi_funct11_strings() {
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

      for(int i=0;i<nv;i++) {
	vars[st_vars[i]]=x[i];
      }

      return calc.eval(&vars);
    }

#ifndef DOXYGEN_INTERNAL

  protected:

    calculator calc;
    /// Desc
    std::map<std::string,double> vars;

    /// The number of variables
    int st_nv;

    /// The expr string
    std::string st_form;
    /// The variable string
    std::vector<std::string> st_vars;
  
    multi_funct11_strings() {}
  
#ifndef DOXYGEN_NO_O2NS
#endif

  private:

    multi_funct11_strings(const multi_funct11_strings &);
    multi_funct11_strings& operator=(const multi_funct11_strings&);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
