/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#ifndef O2SCL_COMP_GEN_INTE_H
#define O2SCL_COMP_GEN_INTE_H

/** \file inte_gen_comp.h
    \brief File defining \ref o2scl::inte_gen_comp
*/

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/funct.h>
#include <o2scl/inte.h>
#include <o2scl/inte_gen.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Naive generalized multi-dimensional integration 

      Naively combine several one-dimensional integration objects from
      class inte in order to perform a multi-dimensional integration.
      The integration routines are specified in the function
      set_ptrs().

      The integration routines are called in order of the index
      specified in the function set_oned_inte(). For
      <tt>n-</tt>dimensional integration, <tt>n</tt> one-dimensional
      integration objects should be specified, with indexes <tt>0</tt>
      through <tt>n-1</tt>. The integration routines are called in
      order of their index, so that the outermost integration is done
      by the routine specified with index 0. The integral 
      performed is:
      \f[ 
      \int_{x_0=a_0}^{x_0=b_0} f(x_0) \int_{x_1=a_1(x_0)}^{x_1=b_1(x_0)} 
      f(x_0, x_1) ...
      \int_{x_{\mathrm{n}-1}=a_{\mathrm{n}-1}(x_0,x_1,..,x_{\mathrm{n}-2})}^
      {x_{\mathrm{n}-1}=b_{\mathrm{n}-1}(x_0,x_1,..,x_{\mathrm{n}-2})} 
      f(x_0,x_1,...,x_{\mathrm{n-1}})~d x_{\mathrm{n}-1}~...~d x_1~d x_0
      \f]

      This class is particularly useful if \f$ f_0 \f$ is time-consuming
      to evaluate, and separable from \f$ f_{n-1} \f$ .

      See the discussion about the functions \c func, \c lower and \c
      upper in the documentation for the class inte_gen.

      No error estimate is performed. Error estimation for multiple
      dimension integrals is provided by the Monte Carlo integration
      classes (see \ref mcarlo_section).

      \future Provide an example of usage for this class.
  */
  template<class func_t, class lfunc_t=func_t, 
    class ufunc_t=func_t, class vec_t=boost::numeric::ublas::vector<double> >
    class inte_gen_comp : 
    public inte_gen<func_t,lfunc_t,ufunc_t,vec_t> {

#ifndef DOXYGEN_INTERNAL

  protected:

  /// The independent variable vector
  vec_t *cx;
  /// The function specifying the lower limits
  lfunc_t *lowerp;
  /// The function specifying the upper limits
  ufunc_t *upperp;
  /// The function to be integrated
  func_t *mf;
  /// The number of dimensions
  size_t ndim;

#endif

  public:

  inte_gen_comp() {
      
    nint=0;
    max_dim=100;
  }
    
  virtual ~inte_gen_comp() {
  }
    
  /** \brief Set the one-dimensional integration object with 
      index \c i.
  */
  int set_oned_inte(inte<> &it, size_t i) {

    if (i>=max_dim) {
      O2SCL_ERR("Index >= max_dim in inte_multi_comp::set_oned_inte().",
		  exc_einval);
    }

    if (nint==0) {

      // Create new space
      nint=i+1;
      iptrs=new inte<> *[nint];
      tptrs=new bool[nint];

    } else if (i>nint-1) {
	
      // Create new space and copy old info over
      size_t nint_new=i+1;

      inte<> **iptrs_new=
	new inte<> *[nint_new];
      bool *tptrs_new=new bool[nint_new];

      for(size_t j=0;j<nint;j++) {
	iptrs_new[j]=iptrs[j];
	tptrs_new[j]=tptrs[j];
      }

      delete[] iptrs;
      delete[] tptrs;

      iptrs=iptrs_new;
      tptrs=tptrs_new;
      nint=nint_new;

    }

    iptrs[i]=&it;
    tptrs[i]=true;

    return success;
  }

  /** \brief Integrate function \c func from 
      \f$ {\ell}_i=f_i(x_i) \f$ to \f$ u_i=g_i(x_i) \f$ for 
      \f$ 0<i<\mathrm{n}-1 \f$.
  */
  virtual double ginteg(func_t &func, size_t n, func_t &lower, 
			func_t &upper) {
      
    // Test to make sure the 1-d integrators were set
    bool enough_oned=true;
    if (n>nint) enough_oned=false;
    for(size_t i=0;i<n;i++) {
      if (tptrs[i]==false) enough_oned=false;
    }
      
    if (enough_oned==false) {
      O2SCL_ERR("Too few objects specified with set_ptrs() in minteg().",
	      exc_einval);
      return 0.0;
    }
      
    // Perform integration
    vec_t c(n);

    cx=&c;
    lowerp=&lower;
    upperp=&upper;
    mf=&func;
    ndim=n;
    size_t ix=0;
      
    funct11 fmn=
    std::bind(std::mem_fn<double(double,size_t &)>
	      (&inte_gen_comp<func_t,lfunc_t,ufunc_t,vec_t>::odfunc),
	      this,std::placeholders::_1,ix);
    double res=iptrs[0]->integ(fmn,lower(0,c),upper(0,c));
      
    return res;
  }
  
  /// Return string denoting type ("inte_gen_comp")
  virtual const char *type() { return "inte_gen_comp"; }

  /// The maxiumum number of integration dimensions (default 100)
  size_t max_dim;

#ifndef DOXYGEN_INTERNAL
    
  protected:
    
  /// The one-dimensional integration function
  double odfunc(double x, size_t &ix) {

    double res;

    (*cx)[ix]=x;
    if (ix==ndim-1) {
      res=(*mf)(ndim,(*cx));
    } else {
      size_t ix_next=ix+1;
      res=(*mf)(ix_next,*cx);
      
      funct11 fmn=
	std::bind(std::mem_fn<double(double,size_t &)>
		  (&inte_gen_comp<func_t,lfunc_t,ufunc_t,vec_t>::odfunc),
		  this,std::placeholders::_1,ix_next);

      res*=iptrs[ix]->integ(fmn,(*lowerp)(ix_next,*cx),
			    (*upperp)(ix_next,*cx));
    }
    return res;
  }
  
  /// The size of the integration object arrays
  size_t nint;

  /// Pointers to the integration objects
  inte<> **iptrs;

  /// Flag indicating if integration object has been set
  bool *tptrs;

#endif


  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

