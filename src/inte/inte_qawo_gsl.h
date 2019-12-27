 /*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Jerry Gagelman and Andrew W. Steiner
  
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
#ifndef O2SCL_GSL_INTE_QAWO_H
#define O2SCL_GSL_INTE_QAWO_H

/** \file inte_qawo_gsl.h
    \brief File defining \ref o2scl::inte_qawo_gsl_sin and
    \ref o2scl::inte_qawo_gsl_cos
*/

#include <o2scl/inte.h>
#include <o2scl/inte_qawc_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Adaptive integration for oscillatory integrals (GSL)

      The integral 
      \f[
      \int_a^b f(x) \sin (\omega x)~dx
      \f]
      is computed for some frequency parameter \f$ \omega \f$,
      stored in \ref inte_qawo_gsl_sin::omega .
      
      An adaptive algorithm together with an series-acceleration
      method like that of \ref inte_qags_gsl is used. Those
      subintervals with "large" widths \f$ d \equiv b-a \f$ where \f$
      d\omega > 4 \f$ are computed using a 25-point Clenshaw-Curtis
      integration rule to handle the oscillatory behavior. In order to
      work efficiently, the Chebyshev moments for the particular
      weight function \f$ W \f$ are computed in advance. 
      
      See \ref gslinte_subsect in the User's guide for general
      information about the GSL integration classes.
  */
  template<class func_t> class inte_qawo_gsl_sin : 
  public inte_cheb_gsl<func_t> {
    
  public:

    inte_qawo_gsl_sin() {
      n_levels=10;
      omega=1.0;
    }

    virtual ~inte_qawo_gsl_sin() {}
      
    /// The user-specified frequency (default 1.0)
    double omega;
    
    /// The number of bisection levels (default 10)
    size_t n_levels;
    
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, double a, double b, 
			  double &res, double &err) {

      otable=gsl_integration_qawo_table_alloc
	(omega,b-a,GSL_INTEG_SINE,n_levels);
      
      int status=qawo(func,a,this->tol_abs,this->tol_rel,
		      this->w,this->otable,&res,&err);
      
      gsl_integration_qawo_table_free(otable);
      
      return status;
    }

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The integration workspace
    gsl_integration_qawo_table *otable;

    /** \brief The full GSL integration routine called by integ_err()
	
        \future Remove goto statements.
    */
    int qawo(func_t &func, const double a, const double epsabs, 
	     const double epsrel, inte_workspace_gsl *loc_w, 
	     gsl_integration_qawo_table *wf, double *result, double *abserr) {
      
      double area, errsum;
      double res_ext, err_ext;
      double result0=0.0, abserr0=0.0, resabs0=0.0, resasc0=0.0;
      double tolerance;
      
      double ertest = 0;
      double error_over_large_intervals = 0;
      double reseps = 0, abseps = 0, correc = 0;
      size_t ktmin = 0;
      int roundoff_type1 = 0, roundoff_type2 = 0, roundoff_type3 = 0;
      int error_type = 0, error_type2 = 0;
      
      size_t iteration = 0;
      
      int positive_integrand = 0;
      int extrapolate = 0;
      int extall = 0;
      int disallow_extrapolation = 0;
      
      typename inte_singular_gsl<func_t>::extrap_table table;
      
      double b = a + wf->L ;
      double abs_omega = fabs (wf->omega) ;
      
      /* Initialize results */
      
      loc_w->initialise(a,b);
      
      *result = 0;
      *abserr = 0;

      size_t limit=this->w->limit;
      
      /* Test on accuracy */
      
      double dbl_eps=std::numeric_limits<double>::epsilon();

      if (epsabs <= 0 && (epsrel < 50 * dbl_eps || epsrel < 0.5e-28)) {
	this->last_iter=0;
	std::string estr="Tolerance cannot be achieved with given ";
	estr+="value of tol_abs, "+dtos(epsabs)+", and tol_rel, "+
	  dtos(epsrel)+", in inte_qawo_gsl_sin::qawo().";
	O2SCL_ERR(estr.c_str(),exc_ebadtol);
      }
      
      /* Perform the first integration */

      this->qc25f(func, a, b, wf, 0, 
		  &result0, &abserr0, &resabs0, &resasc0);
      
      loc_w->set_initial_result (result0, abserr0);

      tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (result0));

      if (this->verbose>0) {
	std::cout << "inte_qawo_gsl Iter: " << 1;
	std::cout.setf(std::ios::showpos);
	std::cout << " Res: " << result0;
	std::cout.unsetf(std::ios::showpos);
	std::cout << " Err: " << abserr0
		  << " Tol: " << tolerance << std::endl;
	if (this->verbose>1) {
	  char ch;
	  std::cout << "Press a key and type enter to continue. " ;
	  std::cin >> ch;
	}
      }

      if (abserr0 <= 100 * GSL_DBL_EPSILON * resabs0 && 
	  abserr0 > tolerance) {
	*result = result0;
	*abserr = abserr0;
	  
	this->last_iter=1;
	std::string estr="Cannot reach tolerance because of roundoff error ";
	estr+="on first attempt in inte_qawo_gsl_sin::qawo().";
	O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);
      } else if ((abserr0 <= tolerance && abserr0 != resasc0) || 
		 abserr0 == 0.0) {
	*result = result0;
	*abserr = abserr0;

	this->last_iter=1;
	return success;
      } else if (limit == 1) {
	*result = result0;
	*abserr = abserr0;

	this->last_iter=1;
	std::string estr="A maximum of 1 iteration was insufficient ";
	estr+="in inte_qawo_gsl_sin::qawo().";
	O2SCL_CONV_RET(estr.c_str(),exc_emaxiter,this->err_nonconv);
      }

      /* Initialization */

      this->initialise_table(&table);

      if (0.5 * abs_omega * fabs(b - a) <= 2) {
	this->append_table (&table, result0);
	extall = 1;
      }
      
      area = result0;
      errsum = abserr0;

      res_ext = result0;
      err_ext = GSL_DBL_MAX;
      
      positive_integrand = this->test_positivity (result0, resabs0);

      iteration = 1;

      do {

	size_t current_level;
	double a1, b1, a2, b2;
	double a_i, b_i, r_i, e_i;
	double area1 = 0, area2 = 0, area12 = 0;
	double error1 = 0, error2 = 0, error12 = 0;
	double resasc1=0.0, resasc2=0.0;
	double resabs1=0.0, resabs2=0.0;
	double last_e_i;

	/* Bisect the subinterval with the largest error estimate */

	loc_w->retrieve (&a_i, &b_i, &r_i, &e_i);

	current_level = loc_w->level[loc_w->i] + 1;

	if (current_level >= wf->n) {
	  error_type = -1 ; /* exceeded limit of table */
	  break ;
	}
	
	a1 = a_i;
	b1 = 0.5 * (a_i + b_i);
	a2 = b1;
	b2 = b_i;

	iteration++;

	qc25f(func, a1, b1, wf, current_level, 
	      &area1, &error1, &resabs1, &resasc1);
	qc25f(func, a2, b2, wf, current_level, 
	      &area2, &error2, &resabs2, &resasc2);

	area12 = area1 + area2;
	error12 = error1 + error2;
	last_e_i = e_i;

	/* Improve previous approximations to the integral and test for
	   accuracy.

	   We write these expressions in the same way as the original
	   QUADPACK code so that the rounding errors are the same, which
	   makes testing easier. */

	errsum = errsum + error12 - e_i;
	area = area + area12 - r_i;

	tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area));

	if (this->verbose>0) {
	  std::cout << "inte_qawo_gsl Iter: " << iteration;
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

	if (resasc1 != error1 && resasc2 != error2) {

	  double delta = r_i - area12;

	  if (fabs (delta) <= 1.0e-5 * fabs (area12) && 
	      error12 >= 0.99 * e_i) {

	    if (!extrapolate) {
	      roundoff_type1++;
	    } else {
	      roundoff_type2++;
	    }
	  }

	  if (iteration > 10 && error12 > e_i) {
	    roundoff_type3++;
	  }
	}
	
	/* Test for roundoff and eventually set error flag */
	
	if (roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20) {
	  error_type = 2;       /* round off error */
	}
	
	if (roundoff_type2 >= 5) {
	  error_type2 = 1;
	}
	
	/* set error flag in the case of bad integrand behaviour at
	   a point of the integration range */
	  
	if (loc_w->subinterval_too_small (a1, a2, b2)) {
	  error_type = 4;
	}

	/* append the newly-created intervals to the list */

	loc_w->update (a1, b1, area1, error1, a2, b2, area2, error2);

	if (errsum <= tolerance) {
	  goto compute_result;
	}
	
	if (error_type) {
	  break;
	}
	
	if (iteration >= limit - 1) {
	  error_type = 1;
	  break;
	}
	
	/* set up variables on first iteration */

	if (iteration == 2 && extall) {
	  error_over_large_intervals = errsum;
	  ertest = tolerance;
	  this->append_table (&table, area);
	  continue;
	}
	
	if (disallow_extrapolation) {
	  continue;
	}
	
	if (extall) {
	  error_over_large_intervals += -last_e_i;
          
	  if (current_level < loc_w->maximum_level)
	    {
	      error_over_large_intervals += error12;
	    }
	  
	  if (extrapolate) {
	    goto label70;
	  }
	}
	
	if (this->large_interval(loc_w)) {
	  continue;
	}
	
	if (extall) {
	  extrapolate = 1;
	  loc_w->nrmax = 1;
	} else {
	  /* test whether the interval to be bisected next is the
	     smallest interval. */
	  size_t i = loc_w->i;
	  double width = loc_w->blist[i] - loc_w->alist[i];
          
	  if (0.25 * fabs(width) * abs_omega > 2) {
	    continue;
	  }
          
	  extall = 1;
	  error_over_large_intervals = errsum;
	  ertest = tolerance;
	  continue;
	}

      label70:

	if (!error_type2 && error_over_large_intervals > ertest) {
	  if (this->increase_nrmax (loc_w)) {
	    continue;
	  }
	}
	
	/* Perform extrapolation */

	this->append_table (&table, area);

	if (table.n < 3) {
	  this->reset_nrmax(loc_w);
	  extrapolate = 0;
	  error_over_large_intervals = errsum;
	  continue;
	}

	this->qelg (&table, &reseps, &abseps);

	ktmin++;

	if (ktmin > 5 && err_ext < 0.001 * errsum) {
	  error_type = 5;
	}

	if (abseps < err_ext) {
	  ktmin = 0;
	  err_ext = abseps;
	  res_ext = reseps;
	  correc = error_over_large_intervals;
	  ertest = GSL_MAX_DBL (epsabs, epsrel * fabs (reseps));
	  if (err_ext <= ertest) {
	    break;
	  }
	}

	/* Prepare bisection of the smallest interval. */

	if (table.n == 1) {
	  disallow_extrapolation = 1;
	}

	if (error_type == 5) {
	  break;
	}

	/* work on interval with largest error */

	this->reset_nrmax (loc_w);
	extrapolate = 0;
	error_over_large_intervals = errsum;

      } while (iteration < limit);

      *result = res_ext;
      *abserr = err_ext;

      if (err_ext == GSL_DBL_MAX)
	goto compute_result;

      if (error_type || error_type2) {

	if (error_type2) {
	  err_ext += correc;
	}

	if (error_type == 0)
	  error_type = 3;

	if (result != 0 && area != 0) {
	  if (err_ext / fabs (res_ext) > errsum / fabs (area)) {
	    goto compute_result;
	  }
	} else if (err_ext > errsum) {
	  goto compute_result;
	} else if (area == 0.0) {
	  goto return_error;
	}
      }
      
      /*  Test on divergence. */

      {
	double max_area = GSL_MAX_DBL (fabs (res_ext), fabs (area));

	if (!positive_integrand && max_area < 0.01 * resabs0) {
	  goto return_error;
	}
      }

      {
	double ratio = res_ext / area;

	if (ratio < 0.01 || ratio > 100 || errsum > fabs (area)) {
	  error_type = 6;
	}
      }

      goto return_error;

    compute_result:

      *result = loc_w->sum_results();
      *abserr = errsum;

    return_error:

      if (error_type > 2) {
	error_type--;
      }
      
      this->last_iter=iteration;

      if (error_type == 0) {
	return success;
      } else if (error_type == 1) {
	std::string estr="Number of iterations was insufficient ";
	estr+=" in inte_qawo_gsl_sin::qawo().";
	O2SCL_CONV_RET(estr.c_str(),exc_emaxiter,this->err_nonconv);
      } else if (error_type == 2) {
	std::string estr="Roundoff error prevents tolerance ";
	estr+="from being achieved in inte_qawo_gsl_sin::qawo().";
	O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);
      } else if (error_type == 3) {
	std::string estr="Bad integrand behavior ";
	estr+=" in inte_qawo_gsl_sin::qawo().";
	O2SCL_CONV_RET(estr.c_str(),exc_esing,this->err_nonconv);
      } else if (error_type == 4) {
	std::string estr="Roundoff error detected in extrapolation table ";
	estr+="in inte_qawo_gsl_sin::qawo().";
	O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);
      } else if (error_type == 5) {
	std::string estr="Integral is divergent or slowly convergent ";
	estr+="in inte_qawo_gsl_sin::qawo().";
	O2SCL_CONV_RET(estr.c_str(),exc_ediverge,this->err_nonconv);
      } else if (error_type == -1) {
	std::string estr="Exceeded limit of trigonometric table ";
	estr+="inte_qawo_gsl_sin::qawo()";
	O2SCL_ERR(estr.c_str(),exc_etable);
      } else {
	std::string estr="Could not integrate function in inte_qawo_gsl";
	estr+="::qawo() (it may have returned a non-finite result).";
	O2SCL_ERR(estr.c_str(),exc_efailed);
      }
      
      // No return statement needed since the above if statement
      // always forces a return, but some compilers like having one
      // anyway.
      return o2scl::success;
    }
    
    /// 25-point quadrature for oscillating functions
    void qc25f(func_t &func, double a, double b, 
	       gsl_integration_qawo_table *wf, size_t level, 
	       double *result, double *abserr, double *resabs,
	       double *resasc) {

      const double center = 0.5 * (a + b);
      const double half_length = 0.5 * (b - a);
      
      const double par = omega * half_length;
      
      if (fabs (par) < 2) {
	
	this->gauss_kronrod(func,a,b,result,abserr,resabs,resasc);
      
	return;

      } else {
	
	double *moment;
	double cheb12[13], cheb24[25];
	double result_abs, res12_cos, res12_sin, res24_cos, res24_sin;
	double est_cos, est_sin;
	double c, s;
	size_t i;
	  
	this->inte_cheb_series(func, a, b, cheb12, cheb24);

	if (level >= wf->n) {
	  /* table overflow should not happen, check before calling */
	  O2SCL_ERR("Table overflow in inte_qawo_gsl::qc25f().",
		    exc_esanity);
	  return;
	}

	/* obtain moments from the table */

	moment = wf->chebmo + 25 * level;

	res12_cos = cheb12[12] * moment[12];
	res12_sin = 0 ;

	for (i = 0; i < 6; i++) {
	  size_t k = 10 - 2 * i;
	  res12_cos += cheb12[k] * moment[k];
	  res12_sin += cheb12[k + 1] * moment[k + 1];
	}

	res24_cos = cheb24[24] * moment[24];
	res24_sin = 0 ;

	result_abs = fabs(cheb24[24]) ;

	for (i = 0; i < 12; i++) {
	  size_t k = 22 - 2 * i;
	  res24_cos += cheb24[k] * moment[k];
	  res24_sin += cheb24[k + 1] * moment[k + 1];
	  result_abs += fabs(cheb24[k]) + fabs(cheb24[k+1]);
	}

	est_cos = fabs(res24_cos - res12_cos);
	est_sin = fabs(res24_sin - res12_sin);

	c = half_length * cos(center * omega);
	s = half_length * sin(center * omega);

	if (wf->sine == GSL_INTEG_SINE) {
	  *result = c * res24_sin + s * res24_cos;
	  *abserr = fabs(c * est_sin) + fabs(s * est_cos);
	} else {
	  *result = c * res24_cos - s * res24_sin;
	  *abserr = fabs(c * est_cos) + fabs(s * est_sin);
	}
      
	*resabs = result_abs * half_length;
	*resasc = GSL_DBL_MAX;

	return;
      }
    }
    
    /// Add the oscillating part to the integrand
  virtual double transform(double t, func_t &func) {
    return func(t)*sin(this->omega*t);
    }

#endif
  
    /// Return string denoting type ("inte_qawo_gsl_sin")
    const char *type() { return "inte_qawo_gsl_sin"; }
  
  };

  /** \brief Adaptive integration a function with finite limits of 
      integration (GSL)

      The integral 
      \f[
      \int_a^b f(x) \cos (\omega x)~dx
      \f]
      is computed for some frequency parameter \f$ \omega \f$ .

      This class is exactly analogous to \ref inte_qawo_gsl_sin .
      See that class documentation for more details.
  */
  template<class func_t> class inte_qawo_gsl_cos : 
  public inte_qawo_gsl_sin<func_t> {
    
  public:

    inte_qawo_gsl_cos() {
    }

    virtual ~inte_qawo_gsl_cos() {}
      
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, double a, double b, 
			  double &res, double &err) {
      
      this->otable=gsl_integration_qawo_table_alloc
	(this->omega,b-a,GSL_INTEG_COSINE,this->n_levels);

      int status=this->qawo(func,a,this->tol_abs,this->tol_rel,
			    this->w,this->otable,&res,&err);
      
      gsl_integration_qawo_table_free(this->otable);
      
      return status;
    }

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// Add the oscillating part to the integrand
    virtual double transform(double t, func_t &func) {
      return func(t)*cos(this->omega*t);
    }

#endif
  
    /// Return string denoting type ("inte_qawo_gsl_cos")
    const char *type() { return "inte_qawo_gsl_cos"; }
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
