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
/** \file min_cern.h
    \brief File defining \ref o2scl::min_cern
*/
#ifndef O2SCL_CERN_MINIMIZE_H
#define O2SCL_CERN_MINIMIZE_H

#include <o2scl/funct.h>
#include <o2scl/min.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief One-dimensional minimization (CERNLIB)

      The golden section search is applied in the interval \f$ (a,b) \f$ 
      using a fixed number \f$ n \f$ of function evaluations where 
      \f[
      n=\left[ 2.08 \ln(|a-b|/\mathrm{tol\_abs})+\frac{1}{2}\right]+1 
      \f]
      
      The accuracy depends on the function. A choice of \f$
      \mathrm{tol\_abs}>10^{-8} \f$ usually results in a relative error of
      $x$ which is smaller than or of the order of \f$ \mathrm{tol\_abs}
      \f$ .
      
      This routine strictly searches the interval \f$ (a,b) \f$ . If
      the function is nowhere flat in this interval, then min_bkt()
      will return either \f$ a \f$ or \f$ b \f$ and \ref min_type
      is set to 1. Note that if there are more than 1 local minima
      in the specified interval, there is no guarantee that this
      method will find the global minimum. 
      
      \note The number of function evaluations can be quite large if
      mmin::tol_abs is sufficiently small. If mmin::tol_abs is 
      exactly zero, then the error handler will be called.

      \comment
      If mmin::tol_abs is exactly zero, then \f$ 10^{-8} \f$ 
      will be used instead.
      [5/23/11 - Changed this consistent with the policy of
      throwing on incorrect input.]
      \endcomment

      Based on the CERNLIB routines RMINFC and DMINFC, which was 
      based on \ref Fletcher87, and \ref Krabs83 and is documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d503/top.html
  */
#ifdef DOXYGEN_NO_O2NS
  template<class func_t=funct11> class min_cern : 
  public min_bkt_base
#else
    template<class func_t=funct11> class min_cern :
  public min_bkt_base<func_t> 
#endif
    {
      
#ifndef DOXGYEN_INTERNAL

      protected:

      /// C analog of Fortran's "Nearest integer" function 
      inline int nint(double x) {
	if (x<0.0) return ((int)(x-0.5));
	else return ((int)(x+0.5));
      }
    
#endif      

    public:
      
      min_cern() {
	delta_set=false;
	min_type=0;
      }
    
      /** \brief Calculate the minimum \c min of \c func between
	  \c a and \c b
	  
	  The initial value of \c x is ignored.
	  
	  If there is no minimum in the given interval, then on exit
	  \c x will be equal to either \c a or \c b and \ref min_type
	  will be set to 1 instead of zero. The error handler is not
	  called, as this need not be interpreted as an error.
      */
      virtual int min_bkt(double &x, double a, double b, double &y,
			  func_t &func) {
      
	// The square root of 5
	static const double w5=2.23606797749979;
	static const double hv=(3.0-w5)/2.0, hw=(w5-1.0)/2.0;
	double c, d, v=0.0, fv=0.0, w=0.0, fw=0.0, h;
      
	if (!delta_set) delta=10.0*this->tol_abs;
	
	int n=-1;
	if (this->tol_abs==0.0) {
	  O2SCL_ERR("Tolerance zero in min_cern::min_bkt().",
			exc_einval);
	} else {
	  if (a!=b) n=nint(2.0*log(fabs((a-b)/this->tol_abs)));
	}
	this->last_ntrial=n;
	c=a;
	d=b;
	if (a>b) {
	  c=b;
	  d=a;
	}
	bool llt=true, lge=true;
      
	// This is never an infinite loop so long as 'n' is finite. We
	// should maybe check somehow that 'n' is not too large
	while (true) {
	  h=d-c;
	  if (this->verbose>0) {
	    x=(c+d)/2.0;
	    y=func(x);
	    this->print_iter(x,y,this->last_ntrial-n,((double)n),this->tol_abs,
			     "min_cern");
	  }
	  if (n<0) {
	    x=(c+d)/2.0;
	    y=func(x);
	    if ((fabs(x-a)>delta && fabs(x-b)>delta)) min_type=0;
	    min_type=1;
	    return 0;
	  }
	  if (llt) {
	    v=c+hv*h;
	    fv=func(v);
	  }
	  if (lge) {
	    w=c+hw*h;
	    fw=func(w);
	  }
	  if (fv<fw) {
	    llt=true;
	    lge=false;
	    d=w;
	    w=v;
	    fw=fv;
	  } else {
	    llt=false;
	    lge=true;
	    c=v;
	    v=w;
	    fv=fw;
	  }
	  n--;
	}

	O2SCL_ERR("Failed sanity check in min_cern::min_bkt().",
		  exc_esanity);
	return exc_esanity;
      }
    
      /** \brief Set the value of \f$ \delta \f$ 
	  
	  If this is not called before min_bkt() is used, then the
	  suggested value \f$ \delta=10 \mathrm{tol_abs} \f$ is used.
      */
      int set_delta(double d) {
	if (d>0.0) {
	  delta=d;
	  delta_set=true;
	}
	return 0;
      }
    
      /// Return string denoting type ("min_cern")
      virtual const char *type() { return "min_cern"; }

      /// Type of minimum found
      int min_type;

    protected:

#ifndef DOXYGEN_INTERNAL

      /// The value of delta as specified by the user
      double delta;

      /// True if the value of delta has been set
      bool delta_set;

#endif

    };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
