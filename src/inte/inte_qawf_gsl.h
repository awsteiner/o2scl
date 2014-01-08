/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Jerry Gagelman
  and Andrew W. Steiner
  
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
#ifndef O2SCL_GSL_INTE_QAWF_H
#define O2SCL_GSL_INTE_QAWF_H

/** \file inte_qawf_gsl.h
    \brief File defining \ref o2scl::inte_qawf_gsl_sin and
    \ref o2scl::inte_qawf_gsl_cos
*/

#include <o2scl/inte.h>
#include <o2scl/inte_qawo_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Adaptive integration for oscillatory integrals (GSL)
      
      The Fourier integral
      \f[
      \int_a^{\infty} f(x) \sin(\omega x)~dx
      \f]
      is computed for some frequency parameter \f$ \omega \f$,
      stored in \ref inte_qawo_gsl_sin::omega .

      The integral is computed using the same method as \ref
      inte_qawo_gsl_sin and \ref inte_qawo_gsl_cos over each of the
      subintervals,
      \f{eqnarray*}{
      C_1 &=& [a, a+c] \\
      C_2 &=& [a+c, a+2c] \\
      &\vdots & \\
      C_k &=& [a +(k-1)c,\, a+kc],
      \f}
      where \f$ c = (2\mathrm{floor}(|\omega|)+1)\pi/|\omega|\f$. This
      width is chosen to cover an odd number of periods so that the
      contributions from the intervals alternate in sign and are
      monotonically decreasing when \f$ f \f$ is positive and
      monotonically decreasing. The sum of this sequence of
      contributions is accelerated using the \f$ \varepsilon \f$
      algorithm.
      
      The algorithm uses zero for the relative tolerance \ref
      inte::tol_rel and attempts to compute the integral to an
      overall absolute tolerance set by \ref inte::tol_abs. The
      following strategy is used: on each interval \f$ C_k\f$, the
      algorithm tries to achieve the tolerance
      \f[ 
      \mathrm{TOL}_k = u_k\cdot \epsilon_{\mathrm{abs}}
      \f]
      where \f$ u_k = (1-p)p^{k-1} \f$ and \f$ p = 0.9\f$. The sum of
      the geometric series of contributions from each interval gives
      an overall tolerance of \f$ \epsilon_{\mathrm{abs}}\f$. If the
      integration of a subinterval leads to difficulties then the accu
      racy requirement for subsequent intervals is relaxed,
      \f[ 
      \mathrm{TOL}_k = 
      u_k\cdot \max\{\epsilon_{\mathrm{abs}}, E_1, \ldots, E_{k-1} \}
      \f]
      where \f$ E_k\f$ is the estimated error on the interval \f$ C_k\f$.
      
      See \ref gslinte_subsect in the User's guide for general
      information about the GSL integration classes.

      When verbose output is enabled, this class outputs information
      from both the subintegrations performed by \ref
      inte_qawo_gsl_sin and the overall integration progress in this
      class.

      \todo More documentation and examples for the
      qawf, qawo and qawc integrators.
  */
  template<class func_t> class inte_qawf_gsl_sin : 
  public inte_qawo_gsl_sin<func_t> {
    
  public:
    
    inte_qawf_gsl_sin() {
    }
    
    virtual ~inte_qawf_gsl_sin() {}
    
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, double a, double b, 
			  double &res, double &err) {
      
      this->otable=gsl_integration_qawo_table_alloc
	(this->omega,1.0,GSL_INTEG_SINE,this->n_levels);
      this->cyclew=new inte_workspace_gsl;
      this->cyclew->allocate(this->w->limit);
      
      int status=qawf(func,a,this->tol_abs,&res,&err);
      
      gsl_integration_qawo_table_free(this->otable);
      this->cyclew->free();
      delete this->cyclew;
      
      return status;
    }

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The integration workspace
    inte_workspace_gsl *cyclew;
    
    /** \brief The full GSL integration routine called by integ_err()
     */
    int qawf(func_t &func, const double a, 
	     const double epsabs, double *result, double *abserr) {

      double area, errsum;
      double res_ext, err_ext;
      double correc, total_error = 0.0, truncation_error;

      size_t ktmin = 0;
      size_t iteration = 0;

      typename inte_singular_gsl<func_t>::extrap_table table;

      double cycle;
      //double omega = this->otable->omega;

      const double p = 0.9;
      double factor = 1;
      double initial_eps, eps;
      int error_type = 0;

      /* Initialize results */

      this->w->initialise(a,a);

      *result = 0;
      *abserr = 0;

      size_t limit=this->w->limit;
      /*
	if (limit > this->w->limit) {
	std::string estr="Iteration limit exceeds workspace ";
	estr+="in inte_qawf_gsl::qawf().";
	O2SCL_ERR_RET(estr.c_str(),exc_einval);
	}
      */

      /* Test on accuracy */

      if (epsabs <= 0) {
	std::string estr="The absolute tolerance must be positive ";
	estr+="in inte_qawf_gsl::qawf().";
	O2SCL_ERR_RET(estr.c_str(),exc_ebadtol);
      }

      if (this->omega == 0.0) {
	if (this->otable->sine == GSL_INTEG_SINE) {
	  /* The function sin(w x) f(x) is always zero for w = 0 */

	  *result = 0;
	  *abserr = 0;

	  return success;
	} else {
	  /* The function cos(w x) f(x) is always f(x) for w = 0 */

	  inte_qagiu_gsl<func_t> iu;
	      
	  int status=iu.integ_err(func,a,0.0,*result,*abserr);

	  return status;
	}
      }

      if (epsabs > GSL_DBL_MIN / (1 - p)) {
	eps = epsabs * (1 - p);
      } else {
	eps = epsabs;
      }

      initial_eps = eps;

      area = 0;
      errsum = 0;

      res_ext = 0;
      err_ext = GSL_DBL_MAX;
      correc = 0;
      
      cycle = (2 * floor (fabs (this->omega)) + 1) * 
	M_PI / fabs (this->omega);

      gsl_integration_qawo_table_set_length (this->otable, cycle);

      this->initialise_table (&table);

      for (iteration = 0; iteration < limit; iteration++) {
	double area1, error1, reseps, erreps;

	double a1 = a + iteration * cycle;
	double b1 = a1 + cycle;

	double epsabs1 = eps * factor;

	int status=this->qawo(func,a1,epsabs1,0.0,cyclew,this->otable,
			      &area1,&error1);
	  
	this->w->append_interval (a1, b1, area1, error1);
	  
	factor *= p;

	area = area + area1;
	errsum = errsum + error1;

	/* estimate the truncation error as 50 times the final term */

	truncation_error = 50 * fabs (area1);

	total_error = errsum + truncation_error;

	if (total_error < epsabs && iteration > 4) {
	  goto compute_result;
	}

	if (error1 > correc) {
	  correc = error1;
	}

	if (status) {
	  eps = GSL_MAX_DBL (initial_eps, correc * (1.0 - p));
	}

	if (status && total_error < 10 * correc && iteration > 3) {
	  goto compute_result;
	}

	this->append_table (&table, area);

	if (table.n < 2) {
	  continue;
	}

	this->qelg (&table, &reseps, &erreps);

	ktmin++;

	if (ktmin >= 15 && err_ext < 0.001 * total_error) {
	  error_type = 4;
	}

	if (erreps < err_ext) {
	  ktmin = 0;
	  err_ext = erreps;
	  res_ext = reseps;

	  if (err_ext + 10 * correc <= epsabs)
	    break;
	  if (err_ext <= epsabs && 10 * correc >= epsabs)
	    break;
	}

	if (this->verbose>0) {
	  std::cout << "inte_qawf_gsl Iter: " << iteration;
	  std::cout.setf(std::ios::showpos);
	  std::cout << " Res: " << area;
	  std::cout.unsetf(std::ios::showpos);
	  std::cout << " Err: " << total_error 
		    << " Tol: " << epsabs << std::endl;
	  if (this->verbose>1) {
	    char ch;
	    std::cout << "Press a key and type enter to continue. " ;
	    std::cin >> ch;
	  }
	}

      }

      if (iteration == limit) error_type = 1;

      if (err_ext == GSL_DBL_MAX)
	goto compute_result;

      err_ext = err_ext + 10 * correc;

      *result = res_ext;
      *abserr = err_ext;

      if (error_type == 0) {
	return success ;
      }

      if (res_ext != 0.0 && area != 0.0) {
	if (err_ext / fabs (res_ext) > errsum / fabs (area))
	  goto compute_result;
      } else if (err_ext > errsum) {
	goto compute_result;
      } else if (area == 0.0) {
	goto return_error;
      }

      if (error_type == 4) {
	err_ext = err_ext + truncation_error;
      }

      goto return_error;

    compute_result:

      *result = area;
      *abserr = total_error;

    return_error:

      if (error_type > 2)
	error_type--;

      if (error_type == 0) {
	return success;
      } else if (error_type == 1) {
	std::string estr="Number of iterations was insufficient ";
	estr+=" in inte_qawf_gsl::qawf().";
	O2SCL_ERR_RET(estr.c_str(),exc_emaxiter);
      } else if (error_type == 2) {
	std::string estr="Roundoff error prevents tolerance ";
	estr+="from being achieved in inte_qawf_gsl::qawf().";
	O2SCL_ERR_RET(estr.c_str(),exc_eround);
      } else if (error_type == 3) {
	std::string estr="Bad integrand behavior ";
	estr+=" in inte_qawf_gsl::qawf().";
	O2SCL_ERR_RET(estr.c_str(),exc_esing);
      } else if (error_type == 4) {
	std::string estr="Roundoff error detected in extrapolation table ";
	estr+="in inte_qawf_gsl::qawf().";
	O2SCL_ERR_RET(estr.c_str(),exc_eround);
      } else if (error_type == 5) {
	std::string estr="Integral is divergent or slowly convergent ";
	estr+="in inte_qawf_gsl::qawf().";
	O2SCL_ERR_RET(estr.c_str(),exc_ediverge);
      } else {
	std::string estr="Could not integrate function in inte_qawf_gsl";
	estr+="::qawf() (it may have returned a non-finite result).";
	O2SCL_ERR_RET(estr.c_str(),exc_efailed);
      }
    }

    /// Add the oscillating part to the integrand
  virtual double transform(double t, func_t &func) {
      return func(t)*sin(this->omega*t);
    }

#endif
  
    /// Return string denoting type ("inte_qawf_gsl_sin")
    const char *type() { return "inte_qawf_gsl_sin"; }
  
  };

  /** \brief Adaptive integration a function with finite limits of 
      integration (GSL)

      The Fourier integral
      \f[
      \int_a^{\infty} f(x) \cos(\omega x)~dx
      \f]
      is computed for some frequency parameter \f$ \omega \f$ .

      This class is exactly analogous to \ref inte_qawf_gsl_sin .
      See that class documentation for more details.
  */
  template<class func_t> class inte_qawf_gsl_cos : 
  public inte_qawf_gsl_sin<func_t> {
    
  public:

    inte_qawf_gsl_cos() {
    }

    virtual ~inte_qawf_gsl_cos() {}
      
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, double a, double b, 
			  double &res, double &err) {

      this->otable=gsl_integration_qawo_table_alloc
	(this->omega,b-a,GSL_INTEG_COSINE,this->n_levels);
      this->cyclew=new inte_workspace_gsl;
      this->cyclew->allocate(this->w->limit);
      
      int status=this->qawf(func,a,this->tol_abs,&res,&err);
      
      gsl_integration_qawo_table_free(this->otable);
      this->cyclew->free();
      delete this->cyclew;
      
      return status;
      
    }

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// Add the oscillating part to the integrand
  virtual double transform(double t, func_t &func) {
      return func(t)*cos(this->omega*t);
    }

#endif
  
    /// Return string denoting type ("inte_qawf_gsl_cos")
    const char *type() { return "inte_qawf_gsl_cos"; }
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
