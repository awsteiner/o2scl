/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

/** \file root_brent_gsl.h
    \brief File defining \ref o2scl::root_brent_gsl 
*/

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
  template<class func_t=funct, class fp_t=double> class root_brent_gsl : 
    public root_bkt<func_t,func_t,fp_t> {

  protected:

    /** \brief Floating point-type agnostic version of
	\c gsl_root_test_interval() .
     */
    int test_interval(fp_t xx_lower, fp_t xx_upper, fp_t epsabs,
		      fp_t epsrel) {
      fp_t abs_lower, abs_upper;
	
      if (xx_lower<0.0) abs_lower=-xx_lower;
      else abs_lower=xx_lower;
      if (xx_upper<0.0) abs_upper=-xx_upper;
      else abs_upper=xx_upper;
      
      fp_t min_abs, tolerance;
      if (epsrel<0.0) {
	O2SCL_ERR2("Relative tolerance is negative in ",
		   "root_brent_gsl::test_interval().",o2scl::exc_ebadtol);
      }
      if (epsabs<0.0) {
	O2SCL_ERR2("Absolute tolerance is negative in ",
		   "root_brent_gsl::test_interval().",o2scl::exc_ebadtol);
      }
      if (xx_lower>xx_upper) {
	O2SCL_ERR2("Lower bound larger than upper bound in ",
		   "root_brent_gsl::test_interval().",o2scl::exc_einval);
      }

      if ((xx_lower>0.0 && xx_upper>0.0) ||
	  (xx_lower<0.0 && xx_upper<0.0)) {
	if (abs_lower<abs_upper) min_abs=abs_lower;
	else min_abs=abs_upper;
      } else {
	min_abs=0;
      }

      tolerance=epsabs+epsrel*min_abs;

      // AWS: I could combine this if statement but this form ensures
      // success in case floating point problems imply that
      // xx_lower<xx_upper and xx_lower>=xx_upper are both false.
      
      if (xx_lower<xx_upper) {
	if (xx_upper-xx_lower<tolerance) {
	  return o2scl::success;
	}
      } else {
	if (xx_lower-xx_upper<tolerance) {
	  return o2scl::success;
	}
      }
      
      return o2scl::gsl_continue;
    }
    
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
      
    fp_t tol, m, two=2;
	
    int ac_equal=0;
	
    if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {
      ac_equal=1;
      c=a;
      fc=fa;
      d=b-a;
      e=b-a;
    }
  
    if (o2scl::o2abs(fc) < o2scl::o2abs(fb)) {
      ac_equal=1;
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    
    tol=o2scl::o2abs(b)*std::numeric_limits<fp_t>::epsilon()/two;
    m=(c-b)/two;
  
    if (fb == 0) {
      root=b;
      x_lower=b;
      x_upper=b;
    
      return o2scl::success;
    }
    if (o2scl::o2abs(m) <= tol) {
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
  
    if (o2scl::o2abs(e) < tol || o2scl::o2abs(fa) <= o2scl::o2abs(fb)) {
      // [GSL] Use bisection 
      d=m;            
      e=m;
    } else {

      // [GSL] Use inverse cubic interpolation 
      fp_t p, q, r;   
      fp_t s=fb/fa;
    
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
      fp_t dtmp;
      if (3*m*q-o2scl::o2abs(tol*q)<o2scl::o2abs(e*q)) {
	dtmp=3*m*q-o2scl::o2abs(tol*q);
      } else {
	dtmp=o2scl::o2abs(e*q);
      }
      if (2*p<dtmp) {
	e=d;
	d=p/q;
      } else {
	// [GSL] Interpolation failed, fall back to bisection.
	d=m;
	e=m;
      }
    }
  
    a=b;
    fa=fb;
  
    if (o2scl::o2abs(d) > tol) {
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
  virtual int solve_bkt(fp_t &x1, fp_t x2, func_t &f) {
	
    int status, iter=0;

    int set_ret=set(f,x1,x2);
    if (set_ret!=0) return set_ret;
    
    if (test_form==0) {

      // Test the bracket size

      status=gsl_continue;
      while (status==gsl_continue && iter<this->ntrial) {
      
	iter++;
	iterate(f);
	status=test_interval(x_lower,x_upper,
			     this->tol_abs,this->tol_rel);
      
	if (this->verbose>0) {
	  fp_t y;

	  y=f(root);

	  // This additional temporary seems to be required
	  // for boost::multiprecision types
	  fp_t x_diff=x_upper-x_lower;
	  this->print_iter(root,y,iter,o2scl::o2abs(x_diff),
			   this->tol_abs,"root_brent_gsl (abs)");
	}
      }
	
    } else if (test_form==1) {

      // Test the residual

      status=gsl_continue;
      while (status==gsl_continue && iter<this->ntrial) {
      
	iter++;
	iterate(f);

	fp_t y=f(root);

	if (o2scl::o2abs(y)<this->tol_rel) status=o2scl::success;
      
	if (this->verbose>0) {
	  this->print_iter(root,y,iter,o2scl::o2abs(y),this->tol_rel,
			   "root_brent_gsl (rel)");
	}
      }


    } else {

      // Test the bracket size and the residual

      status=gsl_continue;
      while (status==gsl_continue && iter<this->ntrial) {
      
	iter++;
	iterate(f);
	status=test_interval(x_lower,x_upper,
				      this->tol_abs,this->tol_rel);
      
	if (status==o2scl::success) {
	  fp_t y=f(root);
	  if (o2scl::o2abs(y)>=this->tol_rel) status=gsl_continue;
	  if (this->verbose>0) {
	    this->print_iter(root,y,iter,o2scl::o2abs(y),this->tol_rel,
			     "root_brent_gsl (rel)");
	  }
	} else {
	  if (this->verbose>0) {
	    fp_t y=f(root);
	    // This additional temporary seems to be required
	    // for boost::multiprecision types
	    fp_t x_diff=x_upper-x_lower;
	    this->print_iter(root,y,iter,o2scl::o2abs(x_diff),
			     this->tol_abs,"root_brent_gsl (abs)");
	  }
	}
      }

    }

    x1=root;
  
    if (iter>=this->ntrial) {
      O2SCL_CONV2_RET("Function root_brent_gsl::solve_bkt() exceeded ",
		      "maximum number of iterations.",o2scl::exc_emaxiter,
		      this->err_nonconv);
    }
  
    if (status!=o2scl::success) {
      return status;
    }

    return o2scl::success;
  }

  /// The type of convergence test applied: 0, 1, or 2 (default 0)
  int test_form;
     
  /// Get the most recent value of the root
  fp_t get_root() { return root; }
      
  /// Get the lower limit
  fp_t get_lower() { return x_lower; }
      
  /// Get the upper limit
  fp_t get_upper() { return x_upper; }
    
  /** \brief Set the information for the solver

      This function currently always returns \ref success.
  */
  int set(func_t &ff, fp_t lower, fp_t upper) {
      
    if (lower > upper) {
      fp_t tmp=lower;
      lower=upper;
      upper=tmp;
    }
	
    x_lower=lower;
    x_upper=upper;
    fp_t two=2;
    root=(lower+upper)/two;
  
    fp_t f_lower, f_upper;
    
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
      O2SCL_CONV2_RET("Endpoints don't straddle y=0 in ",
		      "root_brent_gsl::set().",o2scl::exc_einval,
		      this->err_nonconv);
    }
	
    return o2scl::success;
	
  }

#ifndef DOXYGEN_INTERNAL
      
  protected:
      
  /// The present solution estimate
  fp_t root;
  /// The present lower limit
  fp_t x_lower;
  /// The present upper limit
  fp_t x_upper;

  /// \name Storage for solver state
  //@{
  fp_t a, b, c, d, e;
  fp_t fa, fb, fc;
  //@}
      
#endif
 
  };
   
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
