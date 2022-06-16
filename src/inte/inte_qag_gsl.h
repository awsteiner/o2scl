/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Jerry Gagelman and Andrew W. Steiner
  
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
/* 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#ifndef O2SCL_GSL_INTE_QAG_H
#define O2SCL_GSL_INTE_QAG_H

/** \file inte_qag_gsl.h
    \brief File defining \ref o2scl::inte_qag_gsl
*/
#include <o2scl/inte.h>
#include <o2scl/inte_kronrod_gsl.h>
#include <o2scl/funct.h>
#include <o2scl/string_conv.h>
#include <o2scl/smooth_func.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Adaptive numerical integration of a function (without 
      singularities) on a bounded interval (GSL)
      
      Adaptive integration of a univariate function requires two main
      procedures: approximating the integral on a bounded interval,
      and estimating the approximation error. The algorithm
      recursively refines the interval, computing the integral and its
      error estimate on each subinterval, until the total error
      estimate over \b all subintervals falls within the
      user-specified tolerance. The value returned is the sum of the
      (approximated) integrals over all subintervals.
         
      \verbatim embed:rst
      See :ref:`GSL-based integration details` in the User's 
      guide for general information about the GSL integration classes.
      \endverbatim

      \future There are a few fine-tuned parameters which should
      be re-expressed as data members in the convergence tests.
      \future Should QUADPACK parameters in round-off tests be subject 
      to adjustments by the end user?
      \future Add functions to examine the contents of the workspace
      to detect regions where the integrand may be problematic;
      possibly call these functions automatically depending on
      verbosity settings.
  */
  template<class func_t=funct> class inte_qag_gsl : 
    public inte_kronrod_gsl<func_t> {
    
  public:

  /// Create an integrator with the specified rule
  inte_qag_gsl() {
  }
  
  virtual ~inte_qag_gsl() {}
    
  /** \brief Integrate function \c func from \c a to \c b and place
      the result in \c res and the error in \c err
  */
  virtual int integ_err(func_t &func, double a, double b, 
			double &res, double &err) {
    return qag(func,a,b,this->tol_abs,this->tol_rel,&res,&err);
  }

#ifndef DOXYGEN_INTERNAL

  protected:

  /** \brief Perform an adaptive integration given the coefficients,
      and returning \c result

      \future Just move this function to integ_err().
  */
  int qag(func_t &func, const double a, const double b, 
	  const double l_epsabs, const double l_epsrel, 
	  double *result, double *abserr) {
    
    double area, errsum;
    double result0, abserr0, resabs0, resasc0;
    double tolerance;
    size_t iteration = 0;
    int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;
      
    double round_off;
      
    /* Initialize results */
    
    this->w->initialise(a,b);
    
    *result = 0;
    *abserr = 0;
    
    double dbl_eps=std::numeric_limits<double>::epsilon();
    
    if (l_epsabs <= 0 && 
	(l_epsrel < 50 * dbl_eps || l_epsrel < 0.5e-28)) {
      this->last_iter=0;
      std::string estr="Tolerance cannot be achieved with given ";
      estr+="value of tol_abs, "+dtos(l_epsabs)+", and tol_rel, "+
	dtos(l_epsrel)+", in inte_qag_gsl::qag().";
      O2SCL_ERR(estr.c_str(),exc_ebadtol);
    }
	
    /* perform the first integration */
    
    this->gauss_kronrod(func,a,b,&result0,&abserr0,&resabs0,&resasc0);
    
    this->w->set_initial_result(result0,abserr0);
      
    /* Test on accuracy */
      
    tolerance = GSL_MAX_DBL(l_epsabs, l_epsrel * fabs (result0));
      
    /* need IEEE rounding here to match original quadpack behavior */
	  
    round_off=gsl_coerce_double(50 * dbl_eps * resabs0);
	  
    if (abserr0 <= round_off && abserr0 > tolerance) {

      *result = result0;
      *abserr = abserr0;

      // We start with 1 here, because an integration
      // was already performed above
      this->last_iter=1;
	
      std::string estr="Cannot reach tolerance because of roundoff ";
      estr+="error on first attempt in inte_qag_gsl::qag().";
      O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);

    } else if ((abserr0 <= tolerance && 
		abserr0 != resasc0) || abserr0 == 0.0) {
      *result = result0;
      *abserr = abserr0;
	  
      // We start with 1 here, because an integration
      // was already performed above
      this->last_iter=1;

      return success;

    } else if (this->w->limit == 1) {

      *result = result0;
      *abserr = abserr0;
	
      // We start with 1 here, because an integration
      // was already performed above
      this->last_iter=1;

      O2SCL_CONV2_RET("A maximum of 1 iteration was insufficient ",
		      "in inte_qag_gsl::qag().",
		      exc_emaxiter,this->err_nonconv);
    }
      
    area = result0;
    errsum = abserr0;
      
    iteration = 1;
    
    do {
      
      double a1, b1, a2, b2;
      double a_i, b_i, r_i, e_i;
      double area1 = 0, area2 = 0, area12 = 0;
      double error1 = 0, error2 = 0, error12 = 0;
      double resasc1, resasc2;
      double resabs1, resabs2;
	  
      /* Bisect the subinterval with the largest error estimate */
	  
      this->w->retrieve (&a_i, &b_i, &r_i, &e_i);
	  
      a1 = a_i;
      b1 = 0.5 * (a_i + b_i);
      a2 = b1;
      b2 = b_i;
      
      this->gauss_kronrod(func,a1,b1,&area1,&error1,&resabs1,&resasc1);
      this->gauss_kronrod(func,a2,b2,&area2,&error2,&resabs2,&resasc2);
      
      area12 = area1 + area2;
      error12 = error1 + error2;
	  
      errsum += (error12 - e_i);
      area += area12 - r_i;
	  
      if (resasc1 != error1 && resasc2 != error2) {
	double delta = r_i - area12;
	  
	if (fabs (delta) <= 1.0e-5 * fabs (area12) && 
	    error12 >= 0.99 * e_i) {
	  roundoff_type1++;
	}
	if (iteration >= 10 && error12 > e_i) {
	  roundoff_type2++;
	}
      }
	
      tolerance = GSL_MAX_DBL (l_epsabs, l_epsrel * fabs (area));
	  
      if (errsum > tolerance) {
	if (roundoff_type1 >= 6 || roundoff_type2 >= 20) {
	  // round off error
	  error_type = 2;
	}
	/* set error flag in the case of bad integrand behaviour at
	   a point of the integration range */
	  
	if (this->w->subinterval_too_small (a1, a2, b2)) {
	  error_type = 3;
	}
      }
	
      this->w->update (a1, b1, area1, error1, a2, b2, area2, error2);
	  
      this->w->retrieve (&a_i, &b_i, &r_i, &e_i);
	  
      if (this->verbose>0) {
	std::cout << "inte_qag_gsl Iter: " << iteration;
	std::cout.setf(std::ios::showpos);
	std::cout << " Res: " << area;
	std::cout.unsetf(std::ios::showpos);
	std::cout << " Err: " << errsum
		  << " Tol: " << tolerance << std::endl;
	if (this->verbose>1) {
	  char ch;
	  std::cout << "Press a key and type enter to continue. " ;
	  std::cin >> ch;
	}
      }
	  
      iteration++;
	
    } while (iteration < this->w->limit && !error_type && 
	     errsum > tolerance);
	  
    *result = this->w->sum_results();
    *abserr = errsum;
      
    this->last_iter=iteration;
      
    if (errsum <= tolerance) {
      return success;
    } else if (error_type == 2) {
      std::string estr="Roundoff error prevents tolerance ";
      estr+="from being achieved in inte_qag_gsl::qag().";
      O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);
    } else if (error_type == 3) {
      std::string estr="Bad integrand behavior ";
      estr+=" in inte_qag_gsl::qag().";
      O2SCL_CONV_RET(estr.c_str(),exc_esing,this->err_nonconv);
    } else if (iteration == this->w->limit) {
      std::string estr="Maximum number of subdivisions ("+itos(iteration);
      estr+=") reached in inte_qag_gsl::qag().";
      O2SCL_CONV_RET(estr.c_str(),exc_emaxiter,this->err_nonconv);
    } else {
      // AWS, 2/21/22: This error also sometimes occurs when the
      // function returns a finite (but still very large) result.
      // I'm changing this to a convergence error for now.
      std::string estr="Could not integrate function in inte_qag_gsl::";
      estr+="qag() (it may have returned a non-finite result).";
      O2SCL_CONV_RET(estr.c_str(),exc_efailed,this->err_nonconv);
    }
      
    // No return statement needed since the above if statement
    // always forces a return, but some compilers like having one
    // anyway.
    return o2scl::success;
  }
    
#endif

  public:
  
  /// Return string denoting type ("inte_qag_gsl")
  const char *type() { return "inte_qag_gsl"; }
  
  };

  /** \brief Experimental
   */
  template<class func_t=funct> class inte_qag_smooth : 
    public inte<func_t,double> {
    
  protected:
    
    gauss_filter<func_t> gf;
    
  public:

    inte_qag_gsl<func_t> qag;
    
    inte_qag_smooth() {
      qag.err_nonconv=false;
    }
    
    virtual int integ_err(func_t &func, double a, double b, 
			  double &res, double &err) {
      
      int iret=qag.integ_err(func,a,b,res,err);
      //if (iret==0) return 0;
      std::cout << "iqs: " << res << " " << err << " " << iret << std::endl;
      
      funct f1=std::bind(std::mem_fn<double(double)>
                         (&gauss_filter<>::operator()),&gf,
                         std::placeholders::_1);

      for(double factor=1.0e13;factor>1.0e3/1.0001;factor/=1.0e2) {
      
        gf.set_func(func);
        gf.h_rel=(b-a)/factor;
        gf.set_alpha(3.0);
        gf.set_K(10);

        iret=qag.integ_err(f1,a,b,res,err);
        std::cout << factor << " " << res << " " << err << " "
                  << res/factor << " " << iret << std::endl;
        if (iret==0) {
          err=sqrt(err*err+pow(res/factor,2.0));
          if (err<this->tol_abs || err<this->tol_rel*fabs(res)) {
            //return 0;
          }
        }
      }

      for(double factor=1.0e13;factor>1.0e3/1.0001;factor/=1.0e2) {
      
        gf.set_func(func);
        gf.h_rel=(b-a)/factor;
        gf.set_alpha(3.0);
        gf.set_K(100);
        
        iret=qag.integ_err(f1,a,b,res,err);
        
        std::cout << factor << " " << res << " " << err << " "
                  << res/factor << " " << iret << std::endl;
        if (iret==0) {
          err=sqrt(err*err+pow(res/factor,2.0));
          if (err<this->tol_abs || err<this->tol_rel*fabs(res)) {
            //return 0;
          }
        }
      }

      exit(-1);
      return 0;
    }
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
