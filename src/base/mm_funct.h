/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_MM_FUNCT_H
#define O2SCL_MM_FUNCT_H

/** \file mm_funct.h
    \brief Function object classes for multi-dimensional functions
*/

#include <string>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/fparser.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// Array of multi-dimensional functions typedef
  typedef std::function<
    int(size_t,const boost::numeric::ublas::vector<double> &,
	boost::numeric::ublas::vector<double> &) > mm_funct11;

  /** \brief Array of multi-dimensional functions in an array of strings
   */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class mm_funct11_strings {
    public:
    
    /** \brief Specify the strings
     */
    mm_funct11_strings(int nv, std::string *formulas, std::string vars, 
		     int np=0, std::string parms="") {
      int i;
      fpw=new FunctionParser[nv];
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],vars);
	}
	st_np=0;
	st_parms="";
      } else {
	std::string all=vars+","+parms;
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],all);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_forms=formulas;
      st_vars=vars;
      st_nv=nv;
    }
    
    virtual ~mm_funct11_strings() {
      if (st_np>0) {
	delete[] arr;
      }
      delete[] fpw;
    };
    
    /** \brief Set the values of the auxilliary parameters that were
	specified in 'parms' in the constructor
    */
    int set_parms(const vec_t &p) {
      for(int i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return 0;
    }
    
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y) {
      int i;
      if(st_np<1) {
	for(i=0;i<st_nv;i++) {
	  y[i]=fpw[i].Eval(x);
	}
      } else {
	double *all=new double[st_np+st_nv];
	for(i=0;i<st_nv;i++) all[i]=x[i];
	for(i=st_nv;i<st_np+st_nv;i++) {
	  all[i]=arr[i-st_nv];
	}
	for(i=0;i<st_nv;i++) {
	  y[i]=fpw[i].Eval(all);
	}
	delete[] all;
      }
      return 0;
    }

    /// Set the functions
    int set_function(int nv, std::string *formulas, std::string vars, 
		     int np=0, std::string parms="") {
      int i;
      if (nv!=st_nv) {
	delete[] fpw;
	fpw=new FunctionParser[nv];
      }
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],vars);
	}
	st_np=0;
	st_parms="";
      } else {
	std::string all=vars+","+parms;
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],all);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np+nv];
      }
      st_forms=formulas;
      st_vars=vars;
      st_nv=nv;
      return 0;
    }
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The function parser
    FunctionParser *fpw;

    /// The number of parameters
    int st_np;

    /// The number of variables
    int st_nv;

    /// The arguments to the function parser
    double *arr;

    /// The formulas
    std::string *st_forms;

    /// The variables
    std::string st_vars;

    /// The parameters
    std::string st_parms;
    
    mm_funct11_strings() {};

    private:
    
    mm_funct11_strings(const mm_funct11_strings &);
    mm_funct11_strings& operator=(const mm_funct11_strings&);
    
#endif
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
