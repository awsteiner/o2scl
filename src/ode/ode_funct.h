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
#ifndef O2SCL_ODE_FUNCT_H
#define O2SCL_ODE_FUNCT_H

/** \file ode_funct.h
    \brief Function object classes for ODE functions
*/

#include <string>

#include <boost/numeric/ublas/vector.hpp>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Ordinary differential equation function
   */
  typedef std::function<int(double,size_t,
			    const boost::numeric::ublas::vector<double> &,
			    boost::numeric::ublas::vector<double> &)> 
    ode_funct11;

#ifdef O2SCL_NEVER_DEFINED
  
  /** \brief One-dimensional function from strings
      \nothing
  */
  template <class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t>
    class ode_funct11_strings {
    public:

    /** \brief Specify the string and the parameters
     */
    ode_funct11_strings(size_t nv, std::string *formulas, std::string funcs, 
			std::string var, size_t np=0, std::string parms="") {
      size_t i;
      fpw=new FunctionParser[nv];
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],funcs+","+var);
	}
	st_np=0;
	st_parms="";
      } else {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],funcs+","+var+","+parms);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_forms=formulas;
      st_var=var;
      st_funcs=funcs;
      st_nv=nv;
    }
  
    virtual ~ode_funct11_strings() {
      if (st_np>0) {
	delete[] arr;
      }
      delete[] fpw;
    }
  
    /** \brief Specify the string and the parameters
     */
    int set_function(size_t nv, std::string *formulas, std::string funcs, 
		     std::string var, size_t np=0, std::string parms="") {
      size_t i;
      if (nv!=st_nv) {
	delete[] fpw;
	fpw=new FunctionParser[nv];
      }
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],funcs+","+var);
	}
	st_np=0;
	st_parms="";
      } else {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],funcs+","+var+","+parms);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_forms=formulas;
      st_var=var;
      st_funcs=funcs;
      st_nv=nv;
      return 0;
    }

    /** \brief Set the values of the auxilliary parameters that were
	specified in 'parms' in the constructor
    */
    template<class vec_p_t> void set_parms(vec_p_t &p) {
      for(size_t i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return;
    }
  
    virtual int operator()(double x, size_t nv, const vec_y_t &y, 
			   vec_dydx_t &dydx) {
      size_t i;
      if(st_np<1) {
	double *all=new double[nv+1];
	for(i=0;i<nv;i++) all[i]=y[i];
	all[nv]=x;
	for(i=0;i<st_nv;i++) {
	  dydx[i]=fpw[i].Eval(all);
	}
	delete[] all;
      } else {
	double *all=new double[st_np+st_nv+1];
	for(i=0;i<st_nv;i++) all[i]=y[i];
	all[st_nv]=x;
	for(i=st_nv+1;i<st_np+st_nv+1;i++) all[i]=arr[i-st_nv-1];
	for(i=0;i<st_nv;i++) {
	  dydx[i]=fpw[i].Eval(all);
	}
	delete[] all;
      }
      return 0;
    }

#ifndef DOXYGEN_INTERNAL

    protected:

    /// The formula parser
    FunctionParser *fpw;
    /// The number of parameters
    size_t st_np;
    /// The number of variables
    size_t st_nv;
    /// The parameters
    double *arr;
    /// The formulas
    std::string *st_forms;
    /// The variables
    std::string st_var;
    /// The parameters
    std::string st_parms;
    /// The function names
    std::string st_funcs;

    ode_funct11_strings() {};

    private:

    ode_funct11_strings(const ode_funct11_strings &);
    ode_funct11_strings& operator=(const ode_funct11_strings&);

#endif

  };
  
#endif

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
