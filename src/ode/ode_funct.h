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
#ifndef O2SCL_ODE_FUNCT_H
#define O2SCL_ODE_FUNCT_H

/** \file ode_funct.h
    \brief Function object classes for ODE functions
*/

#include <string>
#include <map>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/calc_utf8.h>

namespace o2scl {

  /** \brief Ordinary differential equation function in src/ode/ode_funct.h
   */
  typedef std::function<int(double,size_t,
			    const boost::numeric::ublas::vector<double> &,
			    boost::numeric::ublas::vector<double> &)> 
    ode_funct;

  /** \brief Simple wrapper to translate Boost format to O2scl format
   */
  template<class vec_t=std::vector<double>,
	   class func_t=std::function<int(double,size_t,const vec_t &,
					  vec_t &)>,
           class fp_t=double> class ode_funct_boost {
    
  protected:

    /// Pointer to function object
    func_t *fp;

    /// Number of variables
    size_t nv2;
    
  public:

    /** \brief Create a wrapper
     */
    ode_funct_boost(func_t &f, size_t nv) {
      fp=&f;
      nv2=nv;
    }
    
    /** \brief Boost interface calling O2scl format
     */
    void operator()(const vec_t &y, vec_t &dydt, fp_t x) {
      (*fp)(x,nv2,y,dydt);
      return;
    }
    
  };
  
  /** \brief One-dimensional function from strings
      \nothing
  */
  class ode_funct_strings {
    
  public:
    
    /** \brief Specify the string and the parameters
     */
    template<class vec_string_t=std::vector<std::string> >
      ode_funct_strings(size_t nv, vec_string_t &exprs,
			  vec_string_t &funcs, std::string var) {
      
      calc.resize(nv);
      st_forms.resize(nv);
      st_funcs.resize(nv);
      for(size_t i=0;i<nv;i++) {
	st_forms[i]=exprs[i];
	calc[i].compile(exprs[i].c_str(),&vars);
	st_funcs[i]=funcs[i];
      }
      st_nv=nv;
      st_var=var;
    }
    
    virtual ~ode_funct_strings() {
    }
    
    /** \brief Specify the string and the parameters
     */
    template<class vec_string_t=std::vector<std::string> >
      void set_function(size_t nv, vec_string_t &exprs,
			vec_string_t &funcs, std::string var) {
      calc.resize(nv);
      st_forms.resize(nv);
      st_funcs.resize(nv);
      for(size_t i=0;i<nv;i++) {
	st_forms[i]=exprs[i];
	st_funcs[i]=funcs[i];
      }
      st_nv=nv;
      st_var=var;
      return;
    }

    /** \brief Set the values of the auxilliary parameters that were
	specified in 'parms' in the constructor
    */
    int set_parm(std::string name, double val) {
      vars[name]=val;
      return 0;
    }
  
    template <class vec_y_t=boost::numeric::ublas::vector<double>,
      class vec_dydx_t=vec_y_t>
      int operator()(double x, size_t nv, const vec_y_t &y, 
		     vec_dydx_t &dydx) {

      for(size_t i=0;i<st_nv;i++) {
	vars[st_funcs[i]]=y[i];
	vars[st_var]=x;
      }
      for(size_t i=0;i<st_nv;i++) {
	dydx[i]=calc[i].eval(&vars);
      }
      return 0;
    }

    protected:

    /// The function parser
    std::vector<o2scl::calc_utf8<> > calc;
      
    /// List of variables and values
    std::map<std::string,double> vars;

    /// The number of variables
    size_t st_nv;
    /// The expressions
    std::vector<std::string> st_forms;
    /// The variables
    std::string st_var;
    /// The function names
    std::vector<std::string> st_funcs;

    ode_funct_strings() {};

    private:

    ode_funct_strings(const ode_funct_strings &);
    ode_funct_strings& operator=(const ode_funct_strings&);

  };
  
}

#endif
