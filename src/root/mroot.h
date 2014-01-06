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
/** \file mroot.h
    \brief File for multi-dimensional solver base class
*/
#ifndef O2SCL_MROOT_H
#define O2SCL_MROOT_H

#include <o2scl/mm_funct.h>
#include <o2scl/jacobian.h>

namespace o2scl {

  /** \brief Multidimensional root-finding [abstract base]
      
      <b>The template parameters:</b>
      The template parameter \c func_t specifies the functions to 
      solve and should be a class containing a definition 
      \code
      func_t::operator()(size_t nv, const vec_t &x, vec_t &y);
      \endcode
      where \c y is the value of the functions at \c x with parameter \c pa
      and \c x and \c y are a array-like classes defining \c operator[] 
      of size \c nv. If the Jacobian matrix is to be specified by the
      user, then the parameter \c jfunc_t specifies the jacobian
      and should contain the definition
      \code
      jfunc_t::operator(size_t nv, vec_t &x, vec_t &y, mat_t &j);
      \endcode
      where \c x is the independent variables, \c y is the array of 
      function values, and \c j is the Jacobian matrix. This template
      parameter can be ignored if only the function msolve() will be used.

      \warning Many of the routines assume that the scale of the
      functions and their variables is of order unity. The solution
      routines may lose precision if this is not the case.

      There is an example for the usage of the multidimensional solver
      classes given in <tt>examples/ex_mroot.cpp</tt>, see \ref
      ex_mroot_sect .

      \future Change ntrial to size_t?
  */
#ifndef O2SCL_CPP11 
  template<class func_t, class vec_t=boost::numeric::ublas::vector<double>,
    class jfunc_t=jac_funct<vec_t,boost::numeric::ublas::matrix<double> > > 
    class mroot
#else
    template<class func_t=mm_funct11, 
    class vec_t=boost::numeric::ublas::vector<double>,
    class jfunc_t=jac_funct11 > class mroot
#endif
    {
      
    public:
    
    mroot() {
      ntrial=100;
      tol_rel=1.0e-8;
      verbose=0;
      tol_abs=1.0e-12;
      last_ntrial=0;
      err_nonconv=true;
    }

    virtual ~mroot() {}
  
    /// The maximum value of the functions for success (default 1.0e-8)
    double tol_rel;
    
    /// The minimum allowable stepsize (default 1.0e-12)
    double tol_abs;
    
    /// Output control (default 0)
    int verbose;
    
    /// Maximum number of iterations (default 100)
    int ntrial;

    /// The number of iterations for in the most recent minimization
    int last_ntrial;
    
    /** \brief If true, call the error handler if msolve() or
	msolve_de() does not converge (default true)
    */
    bool err_nonconv;
    
    /// Return the type, \c "mroot".
    virtual const char *type() { return "mroot"; }
      
    /// Solve \c func using \c x as an initial guess, returning \c x.
    virtual int msolve(size_t n, vec_t &x, func_t &func)=0;
    
    /** \brief Solve \c func with derivatives \c dfunc using \c x as 
	an initial guess, returning \c x.
	
	By default, this function just calls msolve() and ignores the 
	last argument.
    */
    virtual int msolve_de(size_t n, vec_t &x, func_t &func,
			  jfunc_t &dfunc) {
      return msolve(n,x,func);
    }

    /** \brief Print out iteration information.
	  
	Depending on the value of the variable verbose, this prints out
	the iteration information. If verbose=0, then no information is
	printed, while if verbose>1, then after each iteration, the
	present values of x and y are output to std::cout along with the
	iteration number. If verbose>=2 then each iteration waits for a
	character.  

	This is implemented as a template class using a new vector
	type because sometimes the internal vector class is distinct
	from the user-specified vector class (e.g. in \ref
	o2scl::mroot_hybrids).
    */
    template<class vec2_t, class vec3_t>
    int print_iter(size_t n, const vec2_t &x, const vec3_t &y, 
		   int iter, double value=0.0, double limit=0.0,
		   std::string comment="") 
    {
      if (verbose<=0) return o2scl::success;
      
      char ch;
	
      std::cout << comment << " Iteration: " << iter << std::endl;
      for(size_t i=0;i<n;i++) {
	std::cout.width(3);
	std::cout << i;
	if (x[i]>=0.0) {
	  std::cout << "  " << x[i];
	} else {
	  std::cout << " " << x[i];
	}
	if (y[i]>=0.0) {
	  std::cout << "  " << y[i] << std::endl;
	} else {
	  std::cout << " " << y[i] << std::endl;
	}
      }
      std::cout << "Val: " << value << " Lim: " << limit << std::endl;
      if (verbose>1) {
	std::cout << "Press a key and type enter to continue. ";
	std::cin >> ch;
      }
  
      return o2scl::success;
	
    }
      
  };
  
}

#endif
