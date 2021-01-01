/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Jerry Gagelman and Andrew W. Steiner
  
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
#ifndef O2SCL_GSL_INTE_SINGULAR_H
#define O2SCL_GSL_INTE_SINGULAR_H

/** \file inte_singular_gsl.h
    \brief File defining \ref o2scl::inte_singular_gsl
*/

#include <cmath>

#include <gsl/gsl_integration.h>

#include <o2scl/string_conv.h>
#include <o2scl/inte.h>
#include <o2scl/inte_kronrod_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Base class for integrating a function with a 
      singularity (GSL)

      This class contains the extrapolation table mechanics and the
      base integration function for singular integrals from GSL. 
      \verbatim embed:rst
      The
      casual end-user should use the classes described in the
      :ref:`One-dimensional Integration based on GSL, CERNLIB, and Boost`
      section of the User's guide.
      \endverbatim

      \future Some of the functions inside this class could 
      be moved out of header files?
  */
  template<class func_t> class inte_singular_gsl : 
  public inte_kronrod_gsl<func_t> {

  public:
      
    /** \brief A structure for extrapolation for \ref o2scl::inte_qags_gsl

	\future Move this to a new class, with qelg() as a method
    */
    typedef struct extrapolation_table {
      /// Index of new element in the first column
      size_t n;
      /// Lower diagonals of the triangular epsilon table
      double rlist2[52];
      /// Number of calls
      size_t nres;
      /// Three most recent results
      double res3la[3];
    } extrap_table;

  protected:

    /// Initialize the table
    void initialise_table(struct extrapolation_table *table) {
      table->n = 0;
      table->nres = 0;
      return;
    }
      
    /// Append a result to the table
    void append_table(struct extrapolation_table *table, double y) {
      size_t n;
      n = table->n;
      table->rlist2[n] = y;
      table->n++;
      return;
    }
    
    /** \brief Test if the integrand satisfies \f$ f = |f| \f$
     */
    inline int test_positivity(double result, double resabs) {
      double dbl_eps=std::numeric_limits<double>::epsilon();
      int status=(fabs(result) >= (1-50*dbl_eps)*resabs);
      return status;
    }
      
    /** \brief Determines the limit of a given sequence 
	of approximations

	For certain convergent series \f$ \sum_k a_k \f$ whose error
	term \f$ E_n = \sum_{k=n}^\infty a_k \f$ is well behaved, it
	is possible to find a transformation of the sequence that
	yields a faster converging series to the same limit. This
	method of extrapolation applies to some sequences of
	adaptive-approximation and error-estimation for numerical
	integration. 

	\verbatim embed:rst
	This function implements the 
	:math:`\varepsilon`-algorithm ([Wynn56]_, [Piessens83]_) for an
	extrapolation table stored in table.
	\endverbatim

	Quadpack documentation
	\verbatim
	c
	c   list of major variables
	c   -----------------------
	c   e0     - the 4 elements on which the computation of a new
	c   e1       element in the epsilon table is based
	c   e2
	c   e3                 e0
	c                e3    e1    new
	c                      e2
	c   newelm - number of elements to be computed in the new
	c            diagonal
	c   error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
	c   result - the element in the new diagonal with least value
	c            of error
	c
	c   machine dependent constants
	c   ---------------------------
	c
	c   epmach is the largest relative spacing.
	c   oflow is the largest positive magnitude.
	c   limexp is the maximum number of elements the epsilon
	c   table can contain. if this number is reached, the upper
	c   diagonal of the epsilon table is deleted.
	c
	\endverbatim
    */
    void qelg(struct extrapolation_table *table, double *result,
	      double *abserr) {
      
      double *epstab = table->rlist2;
      double *res3la = table->res3la;
      const size_t n = table->n - 1;
	
      const double current = epstab[n];
	
      double absolute = GSL_DBL_MAX;
      double relative = 5 * GSL_DBL_EPSILON * fabs (current);
	
      const size_t newelm = n / 2;
      const size_t n_orig = n;
      size_t n_final = n;
      size_t i;
	
      const size_t nres_orig = table->nres;
	
      *result = current;
      *abserr = GSL_DBL_MAX;
	
      if (n < 2) {
	*result = current;
	*abserr = GSL_MAX_DBL (absolute, relative);
	return;
      }
      
      epstab[n + 2] = epstab[n];
      epstab[n] = GSL_DBL_MAX;
	
      for (i = 0; i < newelm; i++) {
	double res = epstab[n - 2 * i + 2];
	double e0 = epstab[n - 2 * i - 2];
	double e1 = epstab[n - 2 * i - 1];
	double e2 = res;
	    
	double e1abs = fabs (e1);
	double delta2 = e2 - e1;
	double err2 = fabs (delta2);
	double tol2 = GSL_MAX_DBL (fabs (e2), e1abs) * GSL_DBL_EPSILON;
	double delta3 = e1 - e0;
	double err3 = fabs (delta3);
	double tol3 = GSL_MAX_DBL (e1abs, fabs (e0)) * GSL_DBL_EPSILON;
	    
	double e3, delta1, err1, tol1, ss;
	    
	if (err2 <= tol2 && err3 <= tol3) {
	  /* If e0, e1 and e2 are equal to within machine accuracy,
	     convergence is assumed.  */
	
	  *result = res;
	  absolute = err2 + err3;
	  relative = 5 * GSL_DBL_EPSILON * fabs (res);
	  *abserr = GSL_MAX_DBL (absolute, relative);
	  return;
	}
      
	e3 = epstab[n - 2 * i];
	epstab[n - 2 * i] = e1;
	delta1 = e1 - e3;
	err1 = fabs (delta1);
	tol1 = GSL_MAX_DBL (e1abs, fabs (e3)) * GSL_DBL_EPSILON;
	    
	/* If two elements are very close to each other, omit a part of
	   the table by adjusting the value of n */
	    
	if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
	  n_final = 2 * i;
	  break;
	}
	    
	ss = (1 / delta1 + 1 / delta2) - 1 / delta3;
	    
	/* Test to detect irregular behaviour in the table, and
	   eventually omit a part of the table by adjusting the value of
	   n. */
	if (fabs (ss * e1) <= 0.0001) {
	  n_final = 2 * i;
	  break;
	}
	/* Compute a new element and eventually adjust the value of
	   result. */
	    
	res = e1 + 1 / ss;
	epstab[n - 2 * i] = res;
	    
	{
	  const double error = err2 + fabs (res - e2) + err3;
	      
	  if (error <= *abserr) {
	    *abserr = error;
	    *result = res;
	  }
	}
      }
	
      /* Shift the table */
	
      {
	const size_t limexp = 50 - 1;
	  
	if (n_final == limexp) {
	  n_final = 2 * (limexp / 2);
	}
      }
    
      if (n_orig % 2 == 1) {
	for (i = 0; i <= newelm; i++) {
	  epstab[1 + i * 2] = epstab[i * 2 + 3];
	}
      } else {
	for (i = 0; i <= newelm; i++) {
	  epstab[i * 2] = epstab[i * 2 + 2];
	}
      }
      if (n_orig != n_final) {
	for (i = 0; i <= n_final; i++) {
	  epstab[i] = epstab[n_orig - n_final + i];
	}
      }
	
      table->n = n_final + 1;
	
      if (nres_orig < 3) {

	res3la[nres_orig] = *result;
	*abserr = GSL_DBL_MAX;

      } else {                           
	/* Compute error estimate */
	*abserr = (fabs (*result - res3la[2]) + fabs (*result - res3la[1])
		   + fabs (*result - res3la[0]));
	    
	res3la[0] = res3la[1];
	res3la[1] = res3la[2];
	res3la[2] = *result;
      }
	
      /* In QUADPACK the variable table->nres is incremented at the top of
	 qelg, so it increases on every call. This leads to the array
	 res3la being accessed when its elements are still undefined, so I
	 have moved the update to this point so that its value more
	 useful. */
	
      table->nres = nres_orig + 1;
	
      *abserr = GSL_MAX_DBL (*abserr, 5 * GSL_DBL_EPSILON * fabs (*result));
	
      return;
    }
      
    /// Determine if an interval is large
    int large_interval (inte_workspace_gsl * workspace) {
      size_t i = workspace->i ;
      const size_t * level = workspace->level;
	
      if (level[i] < workspace->maximum_level) {
	return 1;
      } else {
	return 0;
      }
    }
      
    /// Reset workspace to work on the interval with the largest error
    inline void reset_nrmax (inte_workspace_gsl * workspace) {
      workspace->nrmax = 0;
      workspace->i = workspace->order[0] ;
    }
      
    /// Increase workspace
    int increase_nrmax (inte_workspace_gsl * workspace) {
      int k;
      int id = workspace->nrmax;
      int jupbnd;
	
      const size_t * level = workspace->level;
      const size_t * order = workspace->order;
	
      size_t limit = workspace->limit ;
      size_t last = workspace->size - 1 ;
	
      if (last > (1 + limit / 2)) {
	jupbnd = limit + 1 - last;
      } else {
	jupbnd = last;
      }
      
      for (k = id; k <= jupbnd; k++) {
	size_t i_max = order[workspace->nrmax];
	
	workspace->i = i_max ;
	
	if (level[i_max] < workspace->maximum_level) {
	  return 1;
	}
	
	workspace->nrmax++;
	
      }
      return 0;
    }

    /** \brief Integration function

        \future Remove goto statements. Before this is done, it might
	be best to add some tests which fail in the various ways.
    */
    int qags(func_t &func, const double a, const double b,
	     const double l_epsabs, const double l_epsrel,
	     double *result, double *abserr) {
      
      double area, errsum;
      double res_ext, err_ext;
      double result0, abserr0, resabs0, resasc0;
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
      int disallow_extrapolation = 0;
	  
      struct extrapolation_table table;
	  
      /* Initialize results */
      
      this->w->initialise(a,b);
      
      *result = 0;
      *abserr = 0;
	  
      size_t limit=this->w->limit;
	  
      /* Test on accuracy */

      if (this->tol_abs <= 0 && (this->tol_rel < 50 * GSL_DBL_EPSILON || 
				 this->tol_rel < 0.5e-28)) {
	this->last_iter=0;
	
	std::string estr="Tolerance cannot be achieved with given ";
	estr+="value of tol_abs, "+dtos(l_epsabs)+", and tol_rel, "+
	  dtos(l_epsrel)+", in inte_singular_gsl::qags().";
	O2SCL_ERR(estr.c_str(),exc_ebadtol);
      }
      
      /* Perform the first integration */
      
      this->gauss_kronrod(func,a,b,&result0,&abserr0,&resabs0,&resasc0);
      
      this->w->set_initial_result (result0, abserr0);
	  
      tolerance = GSL_MAX_DBL (this->tol_abs, this->tol_rel * fabs (result0));
	  
      if (abserr0 <= 100 * GSL_DBL_EPSILON * resabs0 && 
	  abserr0 > tolerance) {

	*result = result0;
	*abserr = abserr0;
	    
	this->last_iter=1;

	std::string estr="Cannot reach tolerance because of roundoff error ";
	estr+="on first attempt in inte_singular_gsl::qags().";
	O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);

      } else if ((abserr0 <= tolerance && 
		  abserr0 != resasc0) || abserr0 == 0.0) {

	*result = result0;
	*abserr = abserr0;
	this->last_iter=1;
	return success;

      } else if (limit == 1) {

	*result = result0;
	*abserr = abserr0;
	    
	this->last_iter=1;
	O2SCL_CONV2_RET("A maximum of 1 iteration was insufficient ",
			"in inte_singular_gsl::qags().",
			exc_emaxiter,this->err_nonconv);
      }
	  
      /* Initialization */
	  
      initialise_table (&table);
      append_table (&table, result0);
	  
      area = result0;
      errsum = abserr0;
	  
      res_ext = result0;
      err_ext = GSL_DBL_MAX;
	  
      positive_integrand = this->test_positivity (result0, resabs0);
	  
      iteration = 1;
	  
      do {

	// Output iteration information
	if (this->verbose>0) {
	  std::cout << this->type();
	  std::cout << " Iter: " << iteration;
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

	size_t current_level;
	double a1, b1, a2, b2;
	double a_i, b_i, r_i, e_i;
	double area1 = 0, area2 = 0, area12 = 0;
	double error1 = 0, error2 = 0, error12 = 0;
	double resasc1, resasc2;
	double resabs1, resabs2;
	double last_e_i;
	    
	/* Bisect the subinterval with the largest error estimate */
	    
	this->w->retrieve (&a_i, &b_i, &r_i, &e_i);
	    
	current_level = this->w->level[this->w->i] + 1;
	    
	a1 = a_i;
	b1 = 0.5 * (a_i + b_i);
	a2 = b1;
	b2 = b_i;
	    
	iteration++;
	
	this->gauss_kronrod(func,a1,b1,&area1,&error1,&resabs1,&resasc1);
	this->gauss_kronrod(func,a2,b2,&area2,&error2,&resabs2,&resasc2);
	
	area12 = area1 + area2;
	error12 = error1 + error2;
	last_e_i = e_i;
	    
	/* Improve previous approximations to the integral and test for
	   accuracy.
	       
	   We write these expressions in the same way as the original
	   QUADPACK code so that the rounding errors are the same, which
	   makes testing easier. 
	*/
	
	errsum = errsum + error12 - e_i;
	area = area + area12 - r_i;
	    
	tolerance = GSL_MAX_DBL (this->tol_abs, this->tol_rel * fabs (area));
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
	    
	// Test for roundoff and eventually set error flag
	    
	if (roundoff_type1 + roundoff_type2 >= 10 || 
	    roundoff_type3 >= 20) {
	  /* round off error */
	  error_type = 2;       
	}
	    
	if (roundoff_type2 >= 5) {
	  error_type2 = 1;
	}
	    
	// Set error flag in the case of bad integrand behaviour at
	// a point of the integration range 
	    
	if (this->w->subinterval_too_small (a1, a2, b2)) {
	  error_type = 4;
	}
	    
	/* append the newly-created intervals to the list */
	    
	this->w->update(a1,b1,area1,error1,a2,b2,area2,error2);
	    
	if (errsum <= tolerance) {
	  
	  // Output final iteration information
	  if (this->verbose>0) {
	    std::cout << this->type();
	    std::cout << " Iter: " << iteration;
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
	  
	  goto compute_result;
	}
	    
	if (error_type) {
	  break;
	}
	    
	if (iteration >= limit - 1) {
	  error_type = 1;
	  break;
	}
	    
	if (iteration == 2) {
	  error_over_large_intervals = errsum;
	  ertest = tolerance;
	  append_table (&table,area);
	  continue;
	}
	    
	if (disallow_extrapolation) {
	  continue;
	}
	    
	error_over_large_intervals += -last_e_i;
	    
	if (current_level < this->w->maximum_level) {
	  error_over_large_intervals += error12;
	}
	    
	if (!extrapolate) {

	  /* test whether the interval to be bisected next is the
	     smallest interval. */
		
	  if (large_interval (this->w)) {
	    continue;
	  }
		
	  extrapolate = 1;
	  this->w->nrmax = 1;
	}

	if (!error_type2 && error_over_large_intervals > ertest) {
	  if (increase_nrmax (this->w)) {
	    continue;
	  }
	}
	    
	// Perform extrapolation
	    
	append_table(&table,area);
	    
	qelg(&table,&reseps,&abseps);
	    
	ktmin++;
	    
	if (ktmin > 5 && err_ext < 0.001 * errsum) {
	  error_type = 5;
	}
	    
	if (abseps < err_ext) {
	  ktmin = 0;
	  err_ext = abseps;
	  res_ext = reseps;
	  correc = error_over_large_intervals;
	  ertest = GSL_MAX_DBL (this->tol_abs,
				this->tol_rel * fabs (reseps));
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
	    
	reset_nrmax (this->w);
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
	      
	if (res_ext != 0.0 && area != 0.0) {
	  if (err_ext / fabs (res_ext) > errsum / fabs (area))
	    goto compute_result;
	} else if (err_ext > errsum) {
	  goto compute_result;
	} else if (area == 0.0) {
	  goto return_error;
	}
      }
	  
      /*  Test on divergence. */
	  
      {
	double max_area = GSL_MAX_DBL (fabs (res_ext),fabs (area));
	    
	if (!positive_integrand && max_area < 0.01 * resabs0)
	  goto return_error;
      }
	  
      {
	double ratio = res_ext / area;
	    
	if (ratio < 0.01 || ratio > 100.0 || errsum > fabs (area))
	  error_type = 6;
      }
	  
      goto return_error;
	  
    compute_result:
	  
      *result = this->w->sum_results();
      *abserr = errsum;
	  
    return_error:
	  
      if (error_type > 2) {
	error_type--;
      }
	  	  
      this->last_iter=iteration;

      if (error_type == 0) {
	return success;
      } else if (error_type == 1) {
	std::string estr="Maximum number of subdivisions ("+itos(iteration);
	estr+=") reached in inte_singular_gsl::qags().";
	O2SCL_CONV_RET(estr.c_str(),exc_emaxiter,this->err_nonconv);
      } else if (error_type == 2) {
	std::string estr="Roundoff error prevents tolerance ";
	estr+="from being achieved in inte_singular_gsl::qags().";
	O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);
      } else if (error_type == 3) {
	std::string estr="Bad integrand behavior ";
	estr+="in inte_singular_gsl::qags().";
	O2SCL_CONV_RET(estr.c_str(),exc_esing,this->err_nonconv);
      } else if (error_type == 4) {
	std::string estr="Roundoff error detected in extrapolation table ";
	estr+="in inte_singular_gsl::qags().";
	O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);
      } else if (error_type == 5) {
	std::string estr="Integral is divergent or slowly convergent ";
	estr+="in inte_singular_gsl::qags().";
	O2SCL_CONV_RET(estr.c_str(),exc_ediverge,this->err_nonconv);
      } 

      std::string estr="Could not integrate function in inte_kronrod_gsl";
      estr+="::qags() (it may have returned a non-finite result).";
      O2SCL_ERR(estr.c_str(),exc_efailed);

      return exc_efailed;
    }                                               

  };

  /** \brief Integrate a function with a singularity (GSL) 
      [abstract base]

      This class contains the GSL-based integration function for 
      applying transformations to the user-defined integrand. 
      \verbatim embed:rst
      The
      casual end-user should use the classes described in the
      :ref:`One-dimensional Integration based on GSL, CERNLIB, and Boost`
      section of the User's guide.
      \endverbatim

  */
  template<class func_t=funct> class inte_transform_gsl : 
  public inte_singular_gsl<func_t> {
    
  public:
  
  /// The transformation to apply to the user-supplied function
  virtual double transform(double t, func_t &func)=0;
  
  /** \brief Integration wrapper for internal transformed function
      type
  */
  virtual void gauss_kronrod
  (func_t &func, double a, double b, 
   double *result, double *abserr, double *resabs, double *resasc) {

    funct fmp=std::bind(std::mem_fn<double(double,func_t &)>
			(&inte_transform_gsl<func_t>::transform),
			this,std::placeholders::_1,std::ref(func));
    
    return this->gauss_kronrod_base
      (fmp,a,b,result,abserr,resabs,resasc);
  }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
