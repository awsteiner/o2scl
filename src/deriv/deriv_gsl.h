/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
/* deriv/deriv.c
 * 
 * Copyright (C) 2004, 2007 Brian Gough
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

#ifndef O2SCL_DERIV_GSL_H
#define O2SCL_DERIV_GSL_H

/** \file deriv_gsl.h
    \brief File defining \ref o2scl::deriv_gsl
*/

#include <iostream>
#include <cmath>
#include <limits>

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

#include <o2scl/deriv.h>
#include <o2scl/funct.h>
#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Numerical differentiation (GSL)

      This class computes the numerical derivative of a function. The
      stepsize \ref h should be specified before use. If similar
      functions are being differentiated in succession, the user may
      be able to increase the speed of later derivatives by setting
      the new stepsize equal to the optimized stepsize from the
      previous differentiation, by setting \ref h to \ref h_opt.

      The results will be incorrect for sufficiently difficult
      functions or if the step size is not properly chosen.

      Some successive derivative computations can be made more
      efficient by using the optimized stepsize in \ref
      deriv_gsl::h_opt , which is set by the most recent last
      derivative computation.

      If the function returns a non-finite value, or if \ref func_max
      is greater than zero and the absolute value of the function is
      larger than \ref func_max, then this class attempts to decrease
      the step size by a factor of 10 in order to compute the
      derivative. The class gives up after 20 reductions of the
      step size. 

      If \ref h is negative or zero, the initial step size is chosen
      to be \f$ 10^{-4} |x| \f$ or if \f$x=0\f$, then the initial step
      size is chosen to be \f$ 10^{-4} \f$ .

      Setting \ref deriv_base::verbose to a number greater than zero
      results in output for each call to \ref central_deriv() which 
      looks like:
      \verbatim
      deriv_gsl: 
      step: 1.000000e-04
      abscissas: 4.999500e-01 4.999000e-01 5.000500e-01 5.001000e-01
      ordinates: 4.793377e-01 4.793816e-01 4.794694e-01 4.795132e-01
      res: 8.775825e-01 trc: 1.462163e-09 rnd: 7.361543e-12
      \endverbatim
      where the last line contains the result (<tt>res</tt>), the
      truncation error (<tt>trc</tt>) and the rounding error
      (<tt>rnd</tt>). If \ref deriv_base::verbose is greater than 1, a
      keypress is required after each iteration.

      If the function always returns a finite value, then computing
      first derivatives requires either 1 or 2 calls to \ref
      central_deriv() and thus either 4 or 8 function calls.

      \note Second and third derivatives are computed by naive nested
      applications of the formula for the first derivative. No
      uncertainty for these derivatives is provided.
      
      An example demonstrating the usage of this class is given in
      <tt>examples/ex_deriv.cpp</tt> and the \ref ex_deriv_sect .

      \future Include the forward and backward GSL derivatives. 
      These would be useful for EOS classes which run in to 
      trouble for negative densities.
  */
  template<class func_t=funct, class fp_t=double> class deriv_gsl : 
    public deriv_base<func_t,fp_t> {
    
  public:
  
    deriv_gsl() {
      h=0.0;
      h_opt=0.0;
      func_max=-1.0;
    }

    virtual ~deriv_gsl() {}
    
    /** \brief Initial stepsize 
	
	This should be specified before a call to deriv() or
	deriv_err(). If it is less than or equal to zero, then \f$ x
	10^{-4} \f$ will used, or if \c x is zero, then \f$ 10^{-4} \f$
	will be used.
    */
    fp_t h;

    /** \brief Maximum absolute value of function, or 
	a negative value for no maximum (default -1)
    */
    fp_t func_max;

    /** \brief The last value of the optimized stepsize

	This is initialized to zero in the constructor and set by
	deriv_err() to the most recent value of the optimized stepsize.
    */
    fp_t h_opt;
  
    /** \brief Calculate the first derivative of \c func  w.r.t. x and 
	uncertainty
    */
    virtual int deriv_err(fp_t x, func_t &func, fp_t &dfdx, fp_t &err) {
      return deriv_tlate<func_t>(x,func,dfdx,err);
    }

    /// Return string denoting type ("deriv_gsl")
    virtual const char *type() { return "deriv_gsl"; }

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Internal template version of the derivative function
     */
    template<class func2_t> int deriv_tlate(fp_t x, func2_t &func, 
					    fp_t &dfdx, fp_t &err) {
      fp_t hh;
      if (h<=0.0) {
	if (x==0.0) hh=1.0e-4;
	else if (x<0.0) h=-1.0e-4*x;
	else hh=1.0e-4*x;
      } else {
	hh=h;
      }

      fp_t r_0, round, trunc, error;
      // Ensure all floating-point constants are initialized by
      // integers
      fp_t one=1, two=2, three=3, ten=10;
      
      size_t it_count=0;
      bool fail=true;
      while (fail && it_count<20) {
      
	fail=false;
      
	int cret=central_deriv(x,hh,r_0,round,trunc,func);
	if (cret!=0) fail=true;

	error=round+trunc;
      
	if (fail==false && round < trunc && (round > 0 && trunc > 0)) {
	  fp_t r_opt, round_opt, trunc_opt, error_opt;
	
	  /* Compute an optimised stepsize to minimize the total error,
	     using the scaling of the truncation error (O(h^2)) and
	     rounding error (O(1/h)). */
	
	  h_opt=hh*pow(round/(two*trunc),one/three);
	  cret=central_deriv(x,h_opt,r_opt,round_opt,trunc_opt,func);
	  if (cret!=0) fail=true;
	  error_opt=round_opt+trunc_opt;
	
	  /* Check that the new error is smaller, and that the new derivative
	     is consistent with the error bounds of the original estimate. */
	
	  if (fail==false && error_opt < error &&
	      abs(r_opt-r_0) < two*two*error) {
	    r_0=r_opt;
	    error=error_opt;
	  }
	}

	it_count++;
	if (fail==true) {
	  hh/=ten;
	  if (this->verbose>0) {
	    std::cout << "Function deriv_gsl::deriv_tlate out of range. "
		      << "Decreasing step." << std::endl;
	  }
	}
      }

      if (fail==true || it_count>=20) {
	if (this->err_nonconv) {
	  O2SCL_ERR2("Failed to find finite derivative in ",
		     "deriv_gsl::deriv_tlate<>.",o2scl::exc_efailed);
	}
	return o2scl::exc_efailed;
      }
      
      dfdx=r_0;
      err=error;
      
      return 0;
    }
  
    /** \brief Internal version of calc_err() for second
	and third derivatives
    */
    virtual int deriv_err_int
    (fp_t x, typename deriv_base<func_t,fp_t>::internal_func_t &func,
     fp_t &dfdx, fp_t &err) {
      return deriv_tlate<>(x,func,dfdx,err);
    }
    
    /** \brief Compute derivative using 5-point rule
	
	Compute the derivative using the 5-point rule (x-h, x-h/2, x,
	x+h/2, x+h) and the error using the difference between the
	5-point and the 3-point rule (x-h,x,x+h). Note that the
	central point is not used for either.

	This must be a class template because it is used by
	both deriv_err() and deriv_err_int().
    */
    template<class func2_t> 
    int central_deriv(fp_t x, fp_t hh, fp_t &result, 
		      fp_t &abserr_round, fp_t &abserr_trunc, 
		      func2_t &func) {
      
      fp_t fm1, fp1, fmh, fph;
    
      // Ensure all floating-point constants are initialized by
      // integers
      fp_t two=2, three=3, four=4, one=1;

      fp_t eps=std::numeric_limits<fp_t>::epsilon();
      
      fm1=func(x-hh);
      fp1=func(x+hh);

      fmh=func(x-hh/two);
      fph=func(x+hh/two);

      if (this->verbose>0) {
	std::cout << "deriv_gsl: " << std::endl;
	std::cout << "step: " << hh << std::endl;
	std::cout << "abscissas: " << x-hh/two << " " << x-hh << " " 
		  << x+hh/two << " " << x+hh << std::endl;
	std::cout << "ordinates: " << fm1 << " " << fmh << " " << fph << " " 
		  << fp1 << std::endl;
      }

      if (!isfinite(fm1) ||
	  !isfinite(fp1) ||
	  !isfinite(fmh) ||
	  !isfinite(fph) ||
	  (func_max>0.0 && (abs(fm1)>func_max ||
			    abs(fp1)>func_max ||
			    abs(fmh)>func_max ||
			    abs(fph)>func_max))) {
	return 1;
      }

      fp_t r3=(fp1-fm1)/two;
      fp_t r5=(four/three)*(fph-fmh)-(one/three)*r3;
      
      fp_t e3=(abs(fp1)+abs(fm1))*eps;
      fp_t e5=two*(abs(fph)+abs(fmh))*eps+e3;
      
      /* The next term is due to finite precision in x+h=O (eps*x) */
      
      fp_t dy=std::max(abs(r3/hh),abs(r5/hh))*abs(x/hh)*eps;
      
      /* The truncation error in the r5 approximation itself is O(h^4).
	 However, for safety, we estimate the error from r5-r3, which is
	 O(h^2).  By scaling h we will minimise this estimated error, not
	 the actual truncation error in r5. 
      */
      
      result=r5/hh;
      /* Estimated truncation error O(h^2) */
      abserr_trunc=abs((r5-r3)/hh); 
      /* Rounding error (cancellations) */
      abserr_round=abs(e5/hh)+dy;   
      
      if (this->verbose>0) {
	std::cout << "res: " << result << " trc: " << abserr_trunc 
		  << " rnd: " << abserr_round << std::endl;
	if (this->verbose>1) {
	  char ch;
	  std::cin >> ch;
	}
      }
      
      return 0;
    }
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



