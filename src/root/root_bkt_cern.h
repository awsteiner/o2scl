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
/** \file root_bkt_cern.h
    \brief File defining \ref o2scl::root_bkt_cern
*/
#ifndef O2SCL_ROOT_BKT_CERN_H
#define O2SCL_ROOT_BKT_CERN_H

#include <o2scl/root.h>

namespace o2scl {

  /** \brief One-dimensional root-finding routine (CERNLIB)
      
      This class attempts to find \f$ x_1 \f$ and \f$ x_2 \f$ in \f$
      [a,b] \f$ such that: 
      - \f$ f(x_1) f(x_2) \leq 0 \f$, 
      - \f$ |f(x_1)| \leq|f(x_2)| \f$, and 
      - \f$ | x_1-x_2| \leq 2~\mathrm{tol\_abs}~(1+|x_0|) \f$. 

      The function solve_bkt() requires inputs \c x1 and \c x2 such
      that the first condition, \f$ f(x_1) f(x_2) \leq 0 \f$, already
      holds.
      
      The variable root::tol_abs defaults to \f$ 10^{-8} \f$ and
      root::ntrial defaults to 200.
      \comment 
      I'm not sure why I chose these values for tol_abs and ntrial above,
      as it doesn't seem to me that these values are suggested in
      CERNLIB. However, they're reasonable so I'll leave them in for
      now. 
      \endcomment
	
      The function solve_bkt() will call the error handler if the root
      is not initially bracketed. If \ref root::err_nonconv is true
      (as it is by default), then the error handler will also be
      called if the number function evaluations is greater than
      root::ntrial.

      After a call to \ref solve_bkt(), \ref root::last_ntrial 
      contains the total number of iterations which were used
      
      \verbatim embed:rst
      See the :ref:`One-dimensional solvers` section of the User's
      guide for general information about O2scl solvers.

      This class is Based on the CERNLIB routines RZEROX and DZEROX,
      which was based on [Bus75]_ and was documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/c200/top.html
      \endverbatim
  */
  template<class func_t=funct> class root_bkt_cern : 
    public root_bkt<func_t> {
    
#ifndef DOXYGEN_INTERNAL

  protected:
      
  /** \brief Internal storage for the mode

      This internal variable is actually defined to be smaller by
      1 than the "mode" as it is defined in the CERNLIB
      documentation in order to avoid needless subtraction in
      solve_bkt().
  */
  int mode;
      
  /// FORTRAN-like function for sign
  inline double sign(double a, double b) {
    if (b>=0.0) return fabs(a);
    return -fabs(a);
  }
      
#endif
      
  public:
      
  root_bkt_cern() {
    mode=0;
    this->tol_abs=1.0e-8;
    this->ntrial=200;
  }
      
  /** \brief Set mode of solution (1 or 2)
         
      - \c 1 should be used for simple functions where the cost is
      inexpensive in comparison to one iteration of solve_bkt(),
      or functions which have a pole near the root (this is the
      default).
      - \c 2 should be used for more time-consuming functions.
 
      If an integer other than \c 1 or \c 2 is specified, 
      the error handler is called.
  */
  int set_mode(int m) {
    if (m!=1 && m!=2) {
      O2SCL_ERR("Invalid mode in root_bkt_cern::set_mode().",
		    o2scl::exc_einval);
    }
    mode=m-1;
    return 0;
  }

  /// Return the type, \c "root_bkt_cern".
  virtual const char *type() { return "root_bkt_cern"; }

  /** \brief Solve \c func in region \f$ x_1<x<x_2 \f$ returning  
      \f$ x_1 \f$.
  */
  virtual int solve_bkt(double &x1, double x2, func_t &func) {
	
    double im1[2]={2.0,3.0}, im2[2]={-1.0,3.0}, c=0.0, fa, fb;
    double atl, a, b, mft;
    double fc=0.0, d=0.0, fd=0.0, tol, h, hb, w, p, q, fdb, fda, f=0.0;
    bool lmt[2];
    int ie=0, loop, mf;
    char ch;
    
    fb=func(x1);
    fa=func(x2);
	
    if (fa*fb>0.0) {
      O2SCL_ERR2("Endpoints don't bracket function in ",
		     "root_bkt_cern::solve_bkt().",exc_einval);
    }
	
    atl=fabs(this->tol_abs);
    b=x1;
    a=x2;
    lmt[1]=true;
    mf=2;
    loop=1;
    do {
      if (loop==1) {
	c=a;
	fc=fa;
	ie=0;
      } else if (loop==2) {
	ie=0;
      }
      if (fabs(fc)<fabs(fb)) {
	if (c!=a) {
	  d=a;
	  fd=fa;
	}
	a=b;
	b=c;
	c=a;
	fa=fb;
	fb=fc;
	fc=fa;
      }
      tol=atl*(1.0+fabs(c));
      h=0.5*(c+b);
      hb=h-b;
      
      if (this->verbose>0) {
	this->print_iter(c,fc,mf-2,fabs(hb),tol,"root_bkt_cern");
      }

      if (fabs(hb)>tol) {
	if (ie>im1[mode]) {
	  w=hb;
	} else {
	  tol*=sign(1.0,hb);
	  p=(b-a)*fb;
	  lmt[0]=(ie<=1);
	  if (lmt[mode]) {
	    q=fa-fb;
	    lmt[1]=false;
	  } else {
	    fdb=(fd-fb)/(d-b);
	    fda=(fd-fa)/(d-a);
	    p*=fda;
	    q=fdb*fa-fda*fb;
	  }
	  if (p<0.0) {
	    p=-p;
	    q=-q;
	  }
	  if (ie==im2[mode]) p+=p;
	  if (p==0.0 || p<=q*tol) {
	    w=tol;
	  } else if (p<hb*q) {
	    w=p/q;
	  } else {
	    w=hb;
	  }
	}
	d=a;
	a=b;
	fd=fa;
	fa=fb;
	b+=w;
	mf++;
	if (mf>this->ntrial) {
	  this->last_ntrial=this->ntrial;
	  O2SCL_CONV2_RET("Too many function calls in ",
			  "root_bkt_cern::solve_bkt().",exc_emaxiter,
			  this->err_nonconv);
	}
	    
	fb=func(b);

	if (fb==0.0 || sign(1.0,fc)==sign(1.0,fb)) {
	  loop=1;
	} else if (w==hb) {
	  loop=2;
	} else {
	  ie++;
	  loop=3;
	}
      } else {
	loop=0;
      }
    } while (loop>0);
    x1=c;
    this->last_ntrial=mf;
    return o2scl::success;
  }
            
  };

}
  
#endif
