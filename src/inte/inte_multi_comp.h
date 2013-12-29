/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef O2SCL_COMPOSITE_INTE_H
#define O2SCL_COMPOSITE_INTE_H

#include <iostream>

#include <o2scl/inte_multi.h>
#include <o2scl/inte.h>
#include <o2scl/funct.h>
#include <o2scl/multi_funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Naive multi-dimensional integration over a hypercube

      Naively combine several one-dimensional integration objects from
      class inte in order to perform a multi-dimensional integration
      over a region defined by constant limits. For more general regions
      of integration, use children of the class \ref inte_gen.
      
      The 1-dimensional integration routines are specified in the function
      set_oned_inte().

      The integration routines are called in order of the index
      specified in the function set_oned_inte(). For
      <tt>n-</tt>dimensional integration, <tt>n</tt> one-dimensional
      integration objects should be specified, with indexes <tt>0</tt>
      through <tt>n-1</tt>. The integration routines are called in
      order of their index, so that the outermost integration is done
      by the routine specified with index 0.
      \f[ 
      \int_{x_0=a_0}^{x_0=b_0} \int_{x_1=a_1}^{x_1=b_1} ...
      \int_{x_{\mathrm{n}-1}=a_{\mathrm{n}-1}}^
      {x_{\mathrm{n}-1}=b_{\mathrm{n}-1}} 
      f(x_0,x_1,...,x_{\mathrm{n}})
      \f]
      
      No error estimate is performed. Error estimation for multiple
      dimension integrals is provided by the Monte Carlo integration
      classes (see \ref o2scl::mcarlo).

      \future Create a function to set an entire array of
      one-dimensional integration objects at once?
      \future Convert the inte<funct> ** to a std::vector<inte<funct> *>
  */
  template<class func_t=multi_funct<>, 
    class vec_t=boost::numeric::ublas::vector<double> 
    > class inte_multi_comp : 
    public inte_multi<func_t,vec_t> {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
  /// The user-specified upper limits
  const vec_t *ax;
  /// The user-specified lower limits
  const vec_t *bx;
  /// The independent variable vector
  vec_t *cx;
  /// The user-specified function
  func_t *mf;
  /// The user-specified number of dimensions
  size_t ndim;
  
#endif

  public:

  inte_multi_comp() {
    nint=0;
    max_dim=100;
  }
    
  virtual ~inte_multi_comp() {
    if (nint>0) {
      delete[] iptrs;
      delete[] tptrs;
    }
  }

  /** \brief Set the one-dimensional integration object with 
      index \c i.
  */
  int set_oned_inte(inte<funct> &it, size_t i) {

    if (i>=max_dim) {
      O2SCL_ERR_RET("Index >= max_dim in inte_multi_comp::set_oned_inte().",
		  exc_einval);
    }

    if (nint==0) {

      // Create new space
      nint=i+1;
      iptrs=new inte<funct> *[nint];
      tptrs=new bool[nint];

    } else if (i>nint-1) {
	
      // Create new space and copy old info over
      size_t nint_new=i+1;

      inte<funct> **iptrs_new=new inte<funct> *[nint_new];
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

  /** \brief Integrate function \c func over the hypercube from 
      \f$ x_i=a_i \f$ to \f$ x_i=b_i \f$ for 
      \f$ 0<i< \f$ ndim-1  
      
      The function set_oned_inte() must be used first to set the 
      one-dimensional routines to use.
  */
  virtual int minteg_err(func_t &func, size_t n, const vec_t &a, 
			 const vec_t &b, double &res, double &err) {
    err=0.0;

    // Test to make sure the 1-d integrators were set
    bool enough_oned=true;
    if (n>nint) enough_oned=false;
    for(size_t i=0;i<n;i++) {
      if (tptrs[i]==false) enough_oned=false;
    }
      
    if (enough_oned==false) {
      O2SCL_ERR_RET("Too few objects specified with set_ptrs() in minteg().",
		    exc_einval);
    }

    // Perform integration
    vec_t c;
    c.resize(n);

    ax=&a;
    bx=&b;
    cx=&c;
    mf=&func;
    ndim=n;
    size_t ix=0;
	
    funct_mfptr_param<inte_multi_comp<func_t,vec_t>,size_t>
    fmn(this,&inte_multi_comp<func_t,vec_t>::odfunc,ix);
    
    res=iptrs[0]->integ(fmn,a[0],b[0]);
    
    return success;
    
  }
    
  /// Return string denoting type ("inte_multi_comp")
  virtual const char *type() { return "inte_multi_comp"; }

  /// The maxiumum number of integration dimensions (default 100)
  size_t max_dim;

  protected:

#ifndef DOXYGEN_INTERNAL

  /// The one-dimensional integration function
  double odfunc(double x, size_t &ix) {

    double res;
	
    (*cx)[ix]=x;
    if (ix==ndim-1) {
      res=(*mf)(ndim,(*cx));
    } else {
      size_t ix_next=ix+1;

      /// This function to send to the integrators
      funct_mfptr_param<inte_multi_comp<func_t,vec_t>,size_t> 
	fmn(this,&inte_multi_comp<func_t,vec_t>::odfunc,ix_next);

      res=iptrs[ix]->integ(fmn,(*ax)[ix+1],(*bx)[ix+1]);
    }
    return res;
  }
    
  /// The size of the integration object arrays
  size_t nint;

  /// Pointers to the integration objects
  inte<funct> **iptrs;

  /// Flag indicating if integration object has been set
  bool *tptrs;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

