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
/* roots/brent.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, 
 * Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */

#ifndef O2SCL_ROOT_BRENT_GSL_H
#define O2SCL_ROOT_BRENT_GSL_H

#include <limits>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <o2scl/funct.h>
#include <o2scl/root.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief One-dimensional root-finding (GSL)

      This class finds the root of a user-specified function. If \ref
      test_form is 0, then solve_bkt() stops when the size of the
      bracket is smaller than \ref root::tol_abs. If \ref test_form is
      1, then the function stops when the residual is less than \ref
      root::tol_rel. If test_form is 2, then both tests are applied.

      See the \ref onedsolve_subsect section of the User's guide for
      general information about \o2 solvers. An example demonstrating
      the usage of this class is given in
      <tt>examples/ex_fptr.cpp</tt> and the \ref ex_fptr_sect .

      \future There is some duplication in the variables \c x_lower, 
      \c x_upper, \c a, and \c b, which could be removed. Some
      better variable names would also be helpful.

      \future Create a meaningful enum list for \ref
      o2scl::root_brent_gsl::test_form.
      
      \comment
      Note that I used \c instead of \ref to refer to variables above
      since the variables a and b are protected, and not available if
      compiling the documentation without the internal portion.

      Also, at first I got confused that GSL didn't require
      lower<upper, but it turns out that this is indeed a requirement
      in GSL, but I didn't see it because it was in roots/fsolver.c
      rather than in roots/brent.c . Thus, everything looks fine now.
      \endcomment
  */
  template<class func_t=funct > class root_brent_gsl : 
  public root_bkt<func_t> {
    
  public:
  
  root_brent_gsl() {
    test_form=0;
  }
    
  /// Return the type, \c "root_brent_gsl".
  virtual const char *type() { return "root_brent_gsl"; }

  /** \brief Perform an iteration

      This function currently always returns \ref success.
  */
  int iterate(func_t &f) {
      
    double tol, m;
	
    int ac_equal=0;
	
    if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {
      ac_equal=1;
      c=a;
      fc=fa;
      d=b-a;
      e=b-a;
    }
  
    if (fabs(fc) < fabs(fb)) {
      ac_equal=1;
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
  
#ifdef O2SCL_CPP11
    tol=0.5*fabs(b)*std::numeric_limits<double>::epsilon();
#else
    tol=0.5*GSL_DBL_EPSILON*fabs(b);
#endif
    m=0.5*(c-b);
  
    if (fb == 0) {
      root=b;
      x_lower=b;
      x_upper=b;
    
      return o2scl::success;
    }
    if (fabs(m) <= tol) {
      root=b;
    
      if (b < c) {
	x_lower=b;
	x_upper=c;
      } else {
	x_lower=c;
	x_upper=b;
      }
    
      return o2scl::success;
    }
  
    if (fabs(e) < tol || fabs(fa) <= fabs(fb)) {
      // use bisection 
      d=m;            
      e=m;
    } else {

      // use inverse cubic interpolation 
      double p, q, r;   
      double s=fb/fa;
    
      if (ac_equal) {
	p=2*m*s;
	q=1-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2*m*q*(q-r)-(b-a)*(r-1));
	q=(q-1)*(r-1)*(s-1);
      }
      
      if (p > 0) {
	q=-q;
      } else {
	p=-p;
      }
      double dtmp;
      if (3*m*q-fabs(tol*q)<fabs(e*q)) dtmp=3*m*q-fabs(tol*q);
      else dtmp=fabs(e*q);
      if (2*p<dtmp) {
	e=d;
	d=p/q;
      } else {
	// Interpolation failed, fall back to bisection.
	d=m;
	e=m;
      }
    }
  
    a=b;
    fa=fb;
  
    if (fabs(d) > tol) {
      b+=d;
    } else {
      b+=(m > 0 ? +tol : -tol);
    }
	
    fb=f(b);
  
    // Update the best estimate of the root and bounds on each
    // iteration
	
    root=b;
  
    if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {
      c=a;
    }
    if (b < c) {
      x_lower=b;
      x_upper=c;
    } else {
      x_lower=c;
      x_upper=b;
    }
  
    return o2scl::success;
  }
      
  /// Solve \c func in region \f$ x_1<x<x_2 \f$ returning \f$ x_1 \f$.
  virtual int solve_bkt(double &x1, double x2, func_t &f) {
	
    int status, iter=0;
	
    if (set(f,x1,x2)!=0) {
      O2SCL_ERR2_RET("Function set() failed in",
		     "root_brent_gsl::solve_bkt().",exc_einval);
    }
  
    if (test_form==0) {

      // Test the bracket size

      status=gsl_continue;
      while (status==gsl_continue && iter<this->ntrial) {
      
	iter++;
	iterate(f);
	status=gsl_root_test_interval(x_lower,x_upper,0.0,this->tol_abs);
      
	if (this->verbose>0) {
	  double y;

	  y=f(root);
	    
	  this->print_iter(root,y,iter,fabs(x_upper-x_lower),this->tol_abs,
			   "root_brent_gsl");
	}
      }
	
    } else if (test_form==1) {

      // Test the residual

      status=gsl_continue;
      while (status==gsl_continue && iter<this->ntrial) {
      
	iter++;
	iterate(f);

	double y=f(root);

	if (fabs(y)<this->tol_rel) status=o2scl::success;
      
	status=gsl_root_test_interval(x_lower,x_upper,0.0,this->tol_abs);
      
	if (this->verbose>0) {
	  this->print_iter(root,y,iter,fabs(x_upper-x_lower),this->tol_abs,
			   "root_brent_gsl");
	}
      }


    } else {

      // Test the bracket size and the residual

      status=gsl_continue;
      while (status==gsl_continue && iter<this->ntrial) {
      
	iter++;
	iterate(f);
	status=gsl_root_test_interval(x_lower,x_upper,0.0,this->tol_abs);
      
	if (status==o2scl::success) {
	  double y=f(root);
	  if (fabs(y)>=this->tol_rel) status=gsl_continue;
	}
	    
	if (this->verbose>0) {
	  double y=f(root);
	  this->print_iter(root,y,iter,fabs(x_upper-x_lower),this->tol_abs,
			   "root_brent_gsl");
	}
      }

    }

    x1=root;
  
    if (status!=o2scl::success) {
      int ret=o2scl::err_hnd->get_errno();
      return ret;
    }

    if (iter>=this->ntrial) {
      O2SCL_CONV2_RET("Function solve_bkt() exceeded maximum number ",
		      "of iterations.",o2scl::exc_emaxiter,
		      this->err_nonconv);
    }
  
    return o2scl::success;
  }

  /// The type of convergence test applied: 0, 1, or 2 (default 0)
  int test_form;
     
  /// Get the most recent value of the root
  double get_root() { return root; }
      
  /// Get the lower limit
  double get_lower() { return x_lower; }
      
  /// Get the upper limit
  double get_upper() { return x_upper; }
    
  /** \brief Set the information for the solver

      This function currently always returns \ref success.
  */
  int set(func_t &ff, double lower, double upper) {
      
    if (lower > upper) {
      double tmp=lower;
      lower=upper;
      upper=tmp;
    }
	
    x_lower=lower;
    x_upper=upper;
    root=0.5*(lower+upper);
  
    double f_lower, f_upper;
  
    f_lower=ff(x_lower);
    f_upper=ff(x_upper);
	
    a=x_lower;
    fa=f_lower;
	
    b=x_upper;
    fb=f_upper;
	
    c=x_upper;
    fc=f_upper;
	
    d=x_upper-x_lower;
    e=x_upper-x_lower;
	
    if ((f_lower<0.0 && f_upper<0.0) || 
	(f_lower>0.0 && f_upper>0.0)) {
      O2SCL_ERR2_RET("Endpoints don't straddle y=0 in ",
		     "root_brent_gsl::set().",o2scl::exc_einval);
    }
	
    return o2scl::success;
	
  }

#ifndef DOXYGEN_INTERNAL
      
  protected:
      
  /// The present solution estimate
  double root;
  /// The present lower limit
  double x_lower;
  /// The present upper limit
  double x_upper;

  /// \name Storage for solver state
  //@{
  double a, b, c, d, e;
  double fa, fb, fc;
  //@}
      
#endif
 
  };
   
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
