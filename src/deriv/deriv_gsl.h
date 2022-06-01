/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

#include <o2scl/misc.h>
#include <o2scl/deriv.h>
#include <o2scl/funct.h>
#include <o2scl/err_hnd.h>

namespace o2scl {

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
      
      \verbatim embed:rst
      An example demonstrating the usage of this class is 
      given in the :ref:`Differentiation example`.
      \endverbatim

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

      // Set equal to 0 to avoid uninitialized variable warnings
      fp_t hh=0;
      if (h<=0.0) {
	if (x==0.0) hh=1.0e-4;
	else hh=o2scl::o2abs(x)*1.0e-4;
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

	  fp_t tdiff=r_opt-r_0;
	  if (fail==false && error_opt < error &&
	      o2scl::o2abs(tdiff) < two*two*error) {
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
  
    /** \brief Internal version of deriv_err() for second
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

      fp_t xph=x+hh;
      fp_t xmh=x-hh;
      fm1=func(xmh);
      fp1=func(xph);

      fp_t xmh2=x-hh/two;
      fp_t xph2=x+hh/two;
      fmh=func(xmh2);
      fph=func(xph2);

      if (this->verbose>0) {
	std::cout << "deriv_gsl: " << std::endl;
	std::cout << "step: " << hh << std::endl;
	std::cout << "abscissas: " << x-hh/two << " " << x-hh << " " 
		  << x+hh/two << " " << x+hh << std::endl;
	std::cout << "ordinates: " << fm1 << " " << fmh << " " << fph << " " 
		  << fp1 << std::endl;
      }
      
      if (!o2isfinite(fm1) || !o2isfinite(fp1) ||
	  !o2isfinite(fmh) || !o2isfinite(fph) ||
	  (func_max>0.0 && (o2scl::o2abs(fm1)>func_max ||
			    o2scl::o2abs(fp1)>func_max ||
			    o2scl::o2abs(fmh)>func_max ||
			    o2scl::o2abs(fph)>func_max))) {
	return 1;
      }

      fp_t r3=(fp1-fm1)/two;
      fp_t r5=(four/three)*(fph-fmh)-(one/three)*r3;
      
      fp_t e3=(o2scl::o2abs(fp1)+o2scl::o2abs(fm1))*eps;
      fp_t e5=two*(o2scl::o2abs(fph)+o2scl::o2abs(fmh))*eps+e3;
      
      /* The next term is due to finite precision in x+h=O (eps*x) */
      fp_t trat0=x/hh;
      fp_t trat1=r3/hh;
      fp_t trat2=r5/hh;
      fp_t dy=std::max(o2scl::o2abs(trat1),o2scl::o2abs(trat2))*
	o2scl::o2abs(trat0)*eps;
      
      /* The truncation error in the r5 approximation itself is O(h^4).
	 However, for safety, we estimate the error from r5-r3, which is
	 O(h^2).  By scaling h we will minimise this estimated error, not
	 the actual truncation error in r5. 
      */
      
      result=r5/hh;
      /* Estimated truncation error O(h^2) */
      fp_t tdiff2=r5-r3;
      fp_t trat3=tdiff2/hh;
      abserr_trunc=o2scl::o2abs(trat3);
      /* Rounding error (cancellations) */
      fp_t trat4=e5/hh;
      abserr_round=o2scl::o2abs(trat4)+dy;
      
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

  /** \brief Evalulate a derivative to within a requested tolerance
      using multiprecision if necessary

      Experimental.

      This class uses \ref funct_multip to ensure the function
      evaluations are sufficiently accurate and then ensures that the
      derivative evaluation is below the requested relative tolerance.
      If the relative tolerance is not specified, then \f$ 10^{-d} \f$
      is used where \f$ d \f$ is the value reported by
      <tt>numeric_limits::digits10</tt> for the input floating point
      type.

      \note Derivatives near zero can be particularly troublesome,
      even for simple functions, since this class only uses relative
      tolerances.

      \verbatim embed:rst
      .. todo:: 

         In class deriv_multip_gsl:

         - More carefully optimize pow_tol_func and the 
           stepsize guesses.
         - Allow the user to specify the value of 'h' 
           for the individual derivative objects.
         - Allow the user to determine how many function evaluations
           or precision required for the last derivative.

      \endverbatim
  */
  class deriv_multip2_gsl {
    
  protected:
    
    typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<50>> cpp_dec_float_50;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

    /** \brief Calculate the first derivative of \c func  w.r.t. x and 
	uncertainty
    */
    template <typename func_t, class fp_t>
    int deriv_err_int(func_t &&f, fp_t x, fp_t &dfdx, fp_t &err,
                      double tol) {
      // Set equal to 0 to avoid uninitialized variable warnings
      double hh=0;
      if (h<=0.0) {
	if (x==0.0) hh=1.0e-4;
	else hh=static_cast<double>(o2scl::o2abs(x))*1.0e-4;
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
      
	int cret=central_deriv(f,x,hh,r_0,round,trunc,tol);
	if (cret!=0) fail=true;

	error=round+trunc;
      
	if (fail==false && round < trunc && (round > 0 && trunc > 0)) {
	  fp_t r_opt, round_opt, trunc_opt, error_opt;
	
	  /* Compute an optimised stepsize to minimize the total error,
	     using the scaling of the truncation error (O(h^2)) and
	     rounding error (O(1/h)). */
	
	  h_opt=hh*static_cast<double>(pow(round/(two*trunc),one/three));
	  cret=central_deriv(f,x,h_opt,r_opt,round_opt,trunc_opt,tol);
	  if (cret!=0) fail=true;
	  error_opt=round_opt+trunc_opt;
	
	  /* Check that the new error is smaller, and that the new derivative
	     is consistent with the error bounds of the original estimate. */

	  fp_t tdiff=r_opt-r_0;
	  if (fail==false && error_opt < error &&
	      o2scl::o2abs(tdiff) < two*two*error) {
	    r_0=r_opt;
	    error=error_opt;
	  }
	}

	it_count++;
	if (fail==true) {
	  hh/=static_cast<double>(ten);
	  if (this->verbose>0) {
	    std::cout << "Function deriv_gsl::deriv_tlate out of range. "
		      << "Decreasing step." << std::endl;
	  }
	}
      }

      if (fail==true || it_count>=20) {
        /*
	if (this->err_nonconv) {
	  O2SCL_ERR2("Failed to find finite derivative in ",
		     "deriv_gsl::deriv_tlate<>.",o2scl::exc_efailed);
	}
        */
	return o2scl::exc_efailed;
      }
      
      dfdx=r_0;
      err=error;
      
      return 0;
    }

    /** \brief Compute derivative using 5-point rule
	
	Compute the derivative using the 5-point rule (x-h, x-h/2, x,
	x+h/2, x+h) and the error using the difference between the
	5-point and the 3-point rule (x-h,x,x+h). Note that the
	central point is not used for either.

	This must be a class template because it is used by
	both deriv_err() and deriv_err_int().
    */
    template <typename func_t, class fp_t>
    int central_deriv(func_t &&func, fp_t x, double hh, fp_t &result, 
		      fp_t &abserr_round, fp_t &abserr_trunc, double tol) {

      funct_multip2 fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=tol;
      
      fp_t fm1, fp1, fmh, fph;
    
      // Ensure all floating-point constants are initialized by
      // integers
      fp_t two=2, three=3, four=4, one=1;

      fp_t eps=std::numeric_limits<fp_t>::epsilon();

      fp_t err;
      
      fp_t xph=x+hh;
      fp_t xmh=x-hh;

      int fm2_ret1=fm2.eval_tol_err(func,xmh,fm1,err);
      if (fm2_ret1!=0) return 1;
      int fm2_ret2=fm2.eval_tol_err(func,xph,fp1,err);
      if (fm2_ret2!=0) return 2;
      
      fp_t xmh2=x-hh/two;
      fp_t xph2=x+hh/two;
      
      int fm2_ret3=fm2.eval_tol_err(func,xmh2,fmh,err);
      if (fm2_ret3!=0) return 3;
      int fm2_ret4=fm2.eval_tol_err(func,xph2,fph,err);
      if (fm2_ret4!=0) return 4;
      
      if (this->verbose>0) {
	std::cout << "deriv_gsl: " << std::endl;
	std::cout << "step: " << hh << std::endl;
	std::cout << "abscissas: " << x-hh/two << " " << x-hh << " " 
		  << x+hh/two << " " << x+hh << std::endl;
	std::cout << "ordinates: " << fm1 << " " << fmh << " " << fph << " " 
		  << fp1 << std::endl;
      }
      
      if (!o2isfinite(fm1) || !o2isfinite(fp1) ||
	  !o2isfinite(fmh) || !o2isfinite(fph)) {
        return 5;
      }

      fp_t r3=(fp1-fm1)/two;
      fp_t r5=(four/three)*(fph-fmh)-(one/three)*r3;
      
      fp_t e3=(o2scl::o2abs(fp1)+o2scl::o2abs(fm1))*eps;
      fp_t e5=two*(o2scl::o2abs(fph)+o2scl::o2abs(fmh))*eps+e3;
      
      /* The next term is due to finite precision in x+h=O (eps*x) */
      fp_t trat0=x/hh;
      fp_t trat1=r3/hh;
      fp_t trat2=r5/hh;
      fp_t dy=std::max(o2scl::o2abs(trat1),o2scl::o2abs(trat2))*
	o2scl::o2abs(trat0)*eps;
      
      /* The truncation error in the r5 approximation itself is O(h^4).
	 However, for safety, we estimate the error from r5-r3, which is
	 O(h^2).  By scaling h we will minimise this estimated error, not
	 the actual truncation error in r5. 
      */
      
      result=r5/hh;
      /* Estimated truncation error O(h^2) */
      fp_t tdiff2=r5-r3;
      fp_t trat3=tdiff2/hh;
      abserr_trunc=o2scl::o2abs(trat3);
      /* Rounding error (cancellations) */
      fp_t trat4=e5/hh;
      abserr_round=o2scl::o2abs(trat4)+dy;
      
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
    
  public:

    /** \brief Initial stepsize 
	
	This should be specified before a call to deriv() or
	deriv_err(). If it is less than or equal to zero, then \f$ x
	10^{-4} \f$ will used, or if \c x is zero, then \f$ 10^{-4} \f$
	will be used.
    */
    double h;

    /** \brief The last value of the optimized stepsize

	This is initialized to zero in the constructor and set by
	deriv_err() to the most recent value of the optimized stepsize.
    */
    double h_opt;
  
    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief Power for tolerance of function evaluations 
        (default 1.33)
     */
    double pow_tol_func;

    /** \brief Verbosity parameter
     */
    int verbose;

    deriv_multip2_gsl() {
      tol_rel=-1.0;
      verbose=0;
      pow_tol_func=1.33;
      h=0.0;
      h_opt=0.0;
    }

    /** \brief Calculate the first derivative of \c func  w.r.t. x and 
	uncertainty
    */
    template<typename func_t, class fp_t>
    int deriv_err(func_t &&f, fp_t x, fp_t &dfdx, fp_t &err,
                  double tol_loc=-1.0) {
      
      if (tol_loc<=0.0) {
        if (tol_rel<=0.0) {
          tol_loc=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          tol_loc=tol_rel;
        }
      } 
      
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err(): set "
                  << "tolerance to: " << tol_loc << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      double tol=pow(tol_loc,pow_tol_func);
      
      double dfdx_d, err_d;
      long double dfdx_ld, err_ld;
      
      deriv_err_int(f,static_cast<double>(x),dfdx_d,err_d,tol);
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "double h_opt is: " << h_opt << std::endl;
      }

      // If the tolerance is large, then just return the result
      // quickly
      if (err_d<tol_loc) {
        dfdx=static_cast<fp_t>(dfdx_d);
        err=static_cast<fp_t>(err_d);
        return 0;
      }

      deriv_err_int(f,static_cast<long double>(x),dfdx_ld,err_ld,tol);
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "long double h_opt is: " << h_opt << std::endl;
      }

      /*
        Sometimes the deriv_gsl class overestimates the uncertainty,
        so we don't always check for, e.g., err_d<tol_loc and
        err_ld<tol_loc.
      */
      if (dfdx_d==0 && dfdx_ld==0 && err_d<tol_loc && err_ld<tol_loc) {
        dfdx=0;
        err=static_cast<fp_t>(err_d);
        if (err<err_ld) err=static_cast<fp_t>(err_ld);
        return 0;
      }
      
      if (dfdx_ld!=0) {
        err=static_cast<fp_t>(abs(dfdx_ld-dfdx_d)/abs(dfdx_ld));
        if (err<tol_loc) {
          dfdx=static_cast<fp_t>(dfdx_ld);
          return 0;
        }
      }
      
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "failed after long double:\n  "
                  << dtos(dfdx_d,0) << " "
                  << dtos(dfdx_ld,0) << " "
                  << dtos(err_d,0) << " "
                  << dtos(err_ld,0) << " "
                  << dtos(err,0) << " " 
                  << tol_loc << std::endl;
      }
    
      // Determine the next stepsize by the previous optimal stepsize
      h=h_opt/10;
      
      cpp_dec_float_25 dfdx_cdf25, err_cdf25;
      deriv_err_int(f,static_cast<cpp_dec_float_25>(x),
                dfdx_cdf25,err_cdf25,tol);
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "cpp_dec_float_25 h_opt is: "
                  << h_opt << std::endl;
      }

      if (dfdx_ld==0 && dfdx_cdf25==0 && err_ld<tol_loc &&
          err_cdf25<tol_loc) {
        dfdx=0;
        err=static_cast<fp_t>(err_ld);
        if (err<err_cdf25) err=static_cast<fp_t>(err_cdf25);
        return 0;
      }
      
      if (dfdx_cdf25!=0) {
        err=static_cast<fp_t>(abs(dfdx_cdf25-dfdx_ld)/abs(dfdx_cdf25));
        if (err<tol_loc) {
          dfdx=static_cast<fp_t>(dfdx_cdf25);
          return 0;
        }
      }
      
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "failed after cpp_dec_float_25:\n  "
                  << dtos(dfdx_ld,0) << " "
                  << dtos(dfdx_cdf25,0) << " "
                  << dtos(err_ld,0) << " "
                  << dtos(err_cdf25,0) << " "
                  << dtos(err,0) << " " 
                  << tol_loc << std::endl;
      }
    
      // Determine the next stepsize by the previous optimal stepsize
      h=h_opt/100;

      cpp_dec_float_35 dfdx_cdf35, err_cdf35;
      deriv_err_int(f,static_cast<cpp_dec_float_35>(x),
                dfdx_cdf35,err_cdf35,tol);
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "cpp_dec_float_35 h_opt is: "
                  << h_opt << std::endl;
      }

      if (dfdx_cdf25==0 && dfdx_cdf35==0 && err_cdf25<tol_loc &&
          err_cdf35<tol_loc) {
        dfdx=0;
        err=static_cast<fp_t>(err_cdf25);
        if (err<err_cdf35) err=static_cast<fp_t>(err_cdf35);
        return 0;
      }
      
      if (dfdx_cdf35!=0) {
        err=static_cast<fp_t>(abs(dfdx_cdf35-dfdx_cdf25)/abs(dfdx_cdf35));
        if (err<tol_loc) {
          dfdx=static_cast<fp_t>(dfdx_cdf35);
          return 0;
        }
      }
      
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "failed after cpp_dec_float_35:\n  "
                  << dtos(dfdx_cdf25,0) << " "
                  << dtos(dfdx_cdf35,0) << " "
                  << dtos(err_cdf25,0) << " "
                  << dtos(err_cdf35,0) << " "
                  << dtos(err,0) << " " 
                  << tol_loc << std::endl;
      }
    
      // Determine the next stepsize by the previous optimal stepsize
      h=h_opt/1e4;
      
      cpp_dec_float_50 dfdx_cdf50, err_cdf50;
      deriv_err_int(f,static_cast<cpp_dec_float_50>(x),
                dfdx_cdf50,err_cdf50,tol);
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "cpp_dec_float_50 h_opt is: "
                  << h_opt << std::endl;
      }

      if (dfdx_cdf35==0 && dfdx_cdf50==0 && err_cdf35<tol_loc &&
          err_cdf50<tol_loc) {
        dfdx=0;
        err=static_cast<fp_t>(err_cdf35);
        if (err<err_cdf50) err=static_cast<fp_t>(err_cdf50);
        return 0;
      }
      
      if (dfdx_cdf50!=0) {
        err=static_cast<fp_t>(abs(dfdx_cdf50-dfdx_cdf35)/abs(dfdx_cdf50));
        if (err<tol_loc) {
          dfdx=static_cast<fp_t>(dfdx_cdf50);
          return 0;
        }
      }
      
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "failed after cpp_dec_float_50:\n  "
                  << dtos(dfdx_cdf35,0) << " "
                  << dtos(dfdx_cdf50,0) << " "
                  << dtos(err_cdf35,0) << " "
                  << dtos(err_cdf50,0) << " "
                  << dtos(err,0) << " " 
                  << tol_loc << std::endl;
      }
    
      // Determine the next stepsize by the previous optimal stepsize
      h=h_opt/1e8;
      
      cpp_dec_float_100 dfdx_cdf100, err_cdf100;
      deriv_err_int(f,static_cast<cpp_dec_float_100>(x),
                dfdx_cdf100,err_cdf100,tol);
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "cpp_dec_float_100 h_opt is: "
                  << h_opt << std::endl;
      }

      if (dfdx_cdf50==0 && dfdx_cdf100==0 && err_cdf50<tol_loc &&
          err_cdf100<tol_loc) {
        dfdx=0;
        err=static_cast<fp_t>(err_cdf50);
        if (err<err_cdf100) err=static_cast<fp_t>(err_cdf100);
        return 0;
      }
      
      if (dfdx_cdf100!=0) {
        err=static_cast<fp_t>(abs(dfdx_cdf100-dfdx_cdf50)/abs(dfdx_cdf100));
        if (err<tol_loc) {
          dfdx=static_cast<fp_t>(dfdx_cdf100);
          return 0;
        }
      }
      
      if (verbose>0) {
        std::cout << "Function deriv_multi_gsl::deriv_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << dtos(dfdx_cdf50,0) << " "
                  << dtos(dfdx_cdf100,0) << " "
                  << dtos(err_cdf50,0) << " "
                  << dtos(err_cdf100,0) << " "
                  << dtos(err,0) << " " 
                  << tol_loc << std::endl;
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in deriv_multip_gsl::deriv_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

  };

}

#endif
