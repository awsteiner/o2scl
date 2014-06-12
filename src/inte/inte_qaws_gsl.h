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
#ifndef GSL_INTE_QAWS_H
#define GSL_INTE_QAWS_H

/** \file inte_qaws_gsl.h
    \brief File defining \ref o2scl::inte_qaws_gsl
*/

#include <o2scl/inte_qawc_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
	
  /** \brief Adaptive integration with with algebraic-logarithmic
      singularities at the end-points (GSL)

      This class computes the weighted integral
      \f[
      \int_a^b f(x)(x - a)^\alpha (b - x)^\beta \log^\mu(x - a) 
      \log^\nu(b - x)~dx
      \f]
      where the parameters of the weight function must satisfy
      \f[
      \alpha > -1, \quad \beta > -1, \quad
      \mu \in \{0, 1\}, \quad \nu \in \{0, 1\},
      \f]
      and which are set by \ref set_weight(). Note that setting \f$
      \mu=0 \f$ or \f$ \nu=0 \f$ removes the respective factor \f$
      \log^mu(\ldots) \f$ or \f$ \log^\nu(\ldots) \f$ from the weight.
	 
      The adaptive refinement algorithm described for \ref
      inte_qag_gsl is used. When a subinterval contains one of the
      endpoints, a special 25-point modified Clenshaw-Curtis rule is
      used to control the singularities. For subintervals which do not
      include the endpoints, a Gauss-Kronrod integration rule is used.

      See \ref gslinte_subsect in the User's guide for general
      information about the GSL integration classes.
  */
  template<class func_t=funct11> class inte_qaws_gsl : 
  public inte_cheb_gsl<func_t> {
    
#ifndef DOXYGEN_INTERNAL
		
  protected:

  /** \name Data from \c gsl_integration_qaws_table
   */
  //@{
  double alpha;
  double beta;
  int mu;
  int nu;
	
  double ri[25];
  double rj[25];
  double rg[25];
  double rh[25];
  //@}
		
  /** \brief Set the array values \c ri, \c rj, \c rg, \c rh from the 
      current values \c alpha and \c beta. 
		 
      This is the function from the GSL source code \c integration/qmomo.c 
      that initializes \c gsl_integration_qaws_table. 
  */
  void initialise_qaws_table() {

    const double alpha_p1 = this->alpha + 1.0;
    const double beta_p1 = this->beta + 1.0;
		
    const double alpha_p2 = this->alpha + 2.0;
    const double beta_p2 = this->beta + 2.0;
		
    const double r_alpha = pow (2.0, alpha_p1);
    const double r_beta = pow (2.0,beta_p1);
		
    size_t i;
		
    double an, anm1;
		
    this->ri[0] = r_alpha / alpha_p1;
    this->ri[1] = this->ri[0] * this->alpha / alpha_p2;
		
    an = 2.0;
    anm1 = 1.0;
		
    for (i = 2; i < 25; i++) {
      this->ri[i] = -(r_alpha + an * (an - alpha_p2) * this->ri[i - 1])
	/ (anm1 * (an + alpha_p1));
      anm1 = an;
      an = an + 1.0;
    }
		
    rj[0] = r_beta / beta_p1;
    rj[1] = rj[0] * this->beta / beta_p2;
		
    an = 2.0;
    anm1 = 1.0;
		
    for (i = 2; i < 25; i++) {
      rj[i] = (-(r_beta + an * (an - beta_p2) * rj[i - 1])
	       / (anm1 * (an + beta_p1)));
      anm1 = an;
      an = an + 1.0;
    }
		
    this->rg[0] = -this->ri[0] / alpha_p1;
    this->rg[1] = -this->rg[0] - 2.0 * r_alpha / (alpha_p2 * alpha_p2);
		
    an = 2.0;
    anm1 = 1.0;
		
    for (i = 2; i < 25; i++) {
      this->rg[i] = (-(an * (an - alpha_p2) * 
		       this->rg[i - 1] - an * this->ri[i - 1]
		       + anm1 * this->ri[i]) 
		     / (anm1 * (an + alpha_p1)));
      anm1 = an;
      an = an + 1.0;
    }
		
    this->rh[0] = -this->rj[0] / beta_p1;
    this->rh[1] = -this->rh[0] - 2.0 * r_beta / (beta_p2 * beta_p2);
		
    an = 2.0;
    anm1 = 1.0;
		
    for (i = 2; i < 25; i++) {
      this->rh[i] = (-(an * (an - beta_p2) * 
		       this->rh[i - 1] - an * this->rj[i - 1]
		       + anm1 * this->rj[i]) 
		     / (anm1 * (an + beta_p1)));
      anm1 = an;
      an = an + 1.0;
    }
		
    for (i = 1; i < 25; i += 2) {
      this->rj[i] *= -1;
      this->rh[i] *= -1;
    }
    
    return;
  }
		
  /** \brief True if algebraic-logarithmic singularity is present at the 
      right endpoint in the definition \c f_trans.
  */
  bool fn_qaws_R;
  
  /** \brief True if algebraic-logarithmic singularity is present at the 
      left endpoint in the definition \c f_trans.
  */
  bool fn_qaws_L;
		
  /// Left endpoint in definition of \c f_trans
  double left_endpoint;

  /// Right endpoint in definition of \c f_trans.
  double right_endpoint;
		
  /** \brief Weighted integrand. 
   */
  virtual double transform(double t, func_t &func) {

    double factor = 1.0,y;
			
    if (fn_qaws_L) {
      if (alpha != 0.0) {
	factor *= pow(t - left_endpoint,alpha);
      }
      if (mu == 1) {
	factor *= log(t - left_endpoint);
      }
    }
			
    if (fn_qaws_R) {
      if (beta != 0.0) {
	factor *= pow(right_endpoint - t,beta);
      }
      if (nu == 1) {
	factor *= log(right_endpoint - t);
      }
    }
    
    return func(t)*factor;
  }
		
  /** \brief Clenshaw-Curtis 25-point integration and error estimator
      for functions with an algebraic-logarithmic singularity at the
      endpoint(s).
  */
  void qc25s(func_t &func, double a, double b, double a1, double b1,
	     double &result, double &abserr, int &err_reliable) {
    
    // Transformed function object for inte_cheb_series()
    funct11 fmp=
    std::bind(std::mem_fn<double(double,func_t &)>
	      (&inte_transform_gsl<func_t>::transform),
	      this,std::placeholders::_1,func);

    this->left_endpoint = a;
    this->right_endpoint = b;
    
    if (a1 == a && (this->alpha != 0.0 || this->mu != 0)) {

      double cheb12[13], cheb24[25];
	
      double factor = pow(0.5 * (b1 - a1),this->alpha + 1.0);
	
      // weighted_function.function = &fn_qaws_R;
      this->fn_qaws_R = true;
      this->fn_qaws_L = false;
      
      this->inte_cheb_series(fmp,a1,b1,cheb12,cheb24);
			
      if (this->mu == 0) {
	double res12 = 0,res24 = 0;
	double u = factor;
				
	this->compute_result(this->ri,cheb12,cheb24,res12,res24);
				
	result = u * res24;
	abserr = fabs(u * (res24 - res12));

      } else {

	double res12a = 0,res24a = 0;
	double res12b = 0,res24b = 0;
	
	double u = factor * log(b1 - a1);
	double v = factor;
				
	this->compute_result(this->ri,cheb12,cheb24,res12a,res24a);
	this->compute_result(this->rg,cheb12,cheb24,res12b,res24b);
				
	result = u * res24a + v * res24b;
	abserr = fabs(u*(res24a - res12a)) + fabs(v*(res24b - res12b));
      }
			
      err_reliable = 0;
      return;

    } else if (b1 == b && (this->beta != 0.0 || this->nu != 0)) {

      double cheb12[13], cheb24[25];
      double factor = pow(0.5 * (b1 - a1), this->beta + 1.0);
			
      // weighted_function.function = &fn_qaws_L;
      this->fn_qaws_L = true;
      this->fn_qaws_R = false;
			
      this->inte_cheb_series(fmp,a1,b1,cheb12,cheb24);
			
      if (this->nu == 0) {

	double res12 = 0, res24 = 0;
	double u = factor;
				
	this->compute_result(this->rj, cheb12,cheb24,res12,res24);
				
	result = u * res24;
	abserr = fabs(u * (res24 - res12));

      } else {

	double res12a = 0, res24a = 0;
	double res12b = 0, res24b = 0;
				
	double u = factor * log(b1 - a1);
	double v = factor;
				
	this->compute_result(this->rj,cheb12,cheb24,res12a,res24a);
	this->compute_result(this->rh,cheb12,cheb24,res12b,res24b);
				
	result = u * res24a + v * res24b;
	abserr = fabs(u*(res24a - res12a)) + fabs(v*(res24b - res12b));
      }
			
      err_reliable = 0;
      return;

    } else {

      double resabs, resasc;
			
      // weighted_function.function = &fn_qaws;
      this->fn_qaws_R = true;
      this->fn_qaws_L = true;
      
      this->gauss_kronrod(func,a1,b1,&result,&abserr,&resabs,&resasc);
			
      if (abserr == resasc) {
	err_reliable = 0;
      } else {
	err_reliable = 1;
      }
	  
      return;
    }
  }
		
  /** \brief Compute the 13-point and 25-point approximations from
      the Chebyshev moments and coefficients. 
  */
  void compute_result(double *r, double *cheb12, double *cheb24,
		      double &result12, double &result24) {
    
    result12=0.0;
    result24=0.0;

    size_t i;
    for (i = 0; i < 13; i++) {
      result12 += r[i] * cheb12[i];
    }
			
    for (i = 0; i < 25; i++) {
      result24 += r[i] * cheb24[i];
    }
  }
		
#endif	
  
  public:
		
  /** \brief Initialize the adptive workspace as with the constructor
      \ref inte_qag_gsl::inte_qag_gsl. 
		 
      The default paramters \f$ \alpha, \beta, \mu, \nu \f$ of the weight
      function are all zero. 
  */
  inte_qaws_gsl() : inte_cheb_gsl<func_t>() {
    set_weight(0.0,0.0,0,0);
  }
  
  ~inte_qaws_gsl() {}
	
  /** \brief Sets the exponents of singularites of the weight function.
		 
      The parameters determine the exponents of the weight function
      \f[
      W(x) = (x-a)^\alpha (b-x)^\beta \log^\mu(x-a) \log^\nu(b-x),
      \f]
      and must satsify
      \f[
      \alpha > -1, \quad \beta > -1, \quad
      \mu \in \{0, 1\}, \quad \nu \in \{0, 1\}.
      \f]
      In order for the adaptive algorithm to run quickly, a table of
      Chebyshev weights for the particular parameters are computed in
      advance.
  */
  int set_weight(double u_alpha, double u_beta, int u_mu, int u_nu) {

    if (u_alpha < -1.0) {
      std::string estr=((std::string)"Variable alpha must be ")+
      "greater than -1.0 in inte_qaws_gsl().";
      O2SCL_ERR(estr.c_str(),exc_einval);
    }
    if (u_beta < -1.0) {
      std::string estr=((std::string)"Variable beta must be ")+
      "greater than -1.0 in inte_qaws_gsl().";
      O2SCL_ERR(estr.c_str(),exc_einval);
    }
    if (u_mu != 0 && u_mu != 1) {
      std::string estr=((std::string)"Variable mu must be 0 or 1 ")+
      "in inte_qaws_gsl().";
      O2SCL_ERR(estr.c_str(),exc_einval);
    }
    if (u_nu != 0 && u_nu != 1) {
      std::string estr=((std::string)"Variable nu must be 0 or 1 ")+
      "in inte_qaws_gsl().";
      O2SCL_ERR(estr.c_str(),exc_einval);
    }
    
    this->alpha = u_alpha;
    this->beta = u_beta;
    this->mu = u_mu;
    this->nu = u_nu;
    
    initialise_qaws_table();
    return success;
  }
		
  /** \brief Returns the current values (via reference) of the 
      weight-function's parameters.
  */
  void get_weight(double &u_alpha, double &u_beta, int &u_mu, int &u_nu) {
    u_alpha = this->alpha;
    u_beta = this->beta;
    u_mu = this->mu;
    u_nu = this->nu;
  }

  /** \brief Integrate the function \c func on the interval (\c a, \c b)
      returning the \c result and error estimate \c abserr.
  */
  virtual int integ_err(func_t &func, double a, double b, 
			double &result, double &abserr) {

    double area, errsum;
    double result0, abserr0;
    double tolerance;
    this->last_iter = 0;
    int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;
		
    /* Initialize results */
		
    this->w->initialise(a,b);
		
    result = 0;
    abserr = 0;
    
    size_t limit=this->w->limit;
		
    if (b <= a) {
      std::string estr="Integration limits, a="+dtos(a);
      estr+=" and b="+dtos(b)+", must satisfy a < b";
      estr+=" in inte_qaws_gsl::gsl_qaws().";
      O2SCL_ERR(estr.c_str(),exc_einval);			
    }
		
#ifndef O2SCL_NO_CPP11
      double dbl_eps=std::numeric_limits<double>::epsilon();
#else 
      double dbl_eps=GSL_DBL_EPSILON;
#endif

    if (this->tol_abs <= 0 && (this->tol_rel < 50 * dbl_eps || 
			       this->tol_rel < 0.5e-28)) {
      this->last_iter=0;
      std::string estr="Tolerance cannot be achieved with given ";
      estr+="value of tol_abs, "+dtos(this->tol_abs)+", and tol_rel, "+
      dtos(this->tol_rel)+", in inte_qaws_gsl::integ_err().";
      O2SCL_ERR(estr.c_str(),exc_ebadtol);
    }
		
    /* perform the first integration */
		
    {
      double area1, area2;
      double error1, error2;
      int err_reliable1, err_reliable2;
      double a1 = a;
      double b1 = 0.5 * (a + b);
      double a2 = b1;
      double b2 = b;
			
      this->qc25s(func, a, b, a1, b1, area1, error1, err_reliable1);
      this->qc25s(func, a, b, a2, b2, area2, error2, err_reliable2);
			
      this->last_iter = 2;
			
      if (error1 > error2) {
	this->w->append_interval(a1, b1, area1, error1);
	this->w->append_interval(a2, b2, area2, error2);
      } else {
	this->w->append_interval(a2, b2, area2, error2);
	this->w->append_interval(a1, b1, area1, error1);
      }
      
      result0 = area1 + area2;
      abserr0 = error1 + error2;
    }
		
    /* Test on accuracy */
		
    tolerance = GSL_MAX_DBL (this->tol_abs, this->tol_rel * fabs (result0));
		
    /* Test on accuracy, use 0.01 relative error as an extra safety
       margin on the first iteration (ignored for subsequent iterations) */
		
    if (abserr0 < tolerance && abserr0 < 0.01 * fabs(result0)) {
      result = result0;
      abserr = abserr0;
      return success;
    } else if (limit == 1) {
      result = result0;
      abserr = abserr0;
      
      std::string estr = "A maximum of 1 iteration was insufficient ";
      estr += "in inte_qaws_gsl::gsl_qaws().";
      O2SCL_CONV_RET(estr.c_str(), exc_emaxiter, this->err_nonconv);
    }
    
    area = result0;
    errsum = abserr0;
		
    do {
	double a1, b1, a2, b2;
	double a_i, b_i, r_i, e_i;
	double area1 = 0, area2 = 0, area12 = 0;
	double error1 = 0, error2 = 0, error12 = 0;
	int err_reliable1, err_reliable2;
			
	/* Bisect the subinterval with the largest error estimate */
	this->w->retrieve(&a_i,&b_i,&r_i,&e_i);
			
	a1 = a_i; 
	b1 = 0.5 * (a_i + b_i);
	a2 = b1;
	b2 = b_i;
			
	qc25s(func, a, b, a1, b1, area1, error1, err_reliable1);
	qc25s(func, a, b, a2, b2, area2, error2, err_reliable2);
			
	area12 = area1 + area2;
	error12 = error1 + error2;
			
	errsum += (error12 - e_i);
	area += area12 - r_i;
			
	if (err_reliable1 && err_reliable2) {

	  double delta = r_i - area12;
	  
	  if (fabs(delta) <= 1.0e-5 * fabs (area12) 
	      && error12 >= 0.99 * e_i) {
	    roundoff_type1++;
	  }
	  if (this->last_iter >= 10 && error12 > e_i) {
	    roundoff_type2++;
	  }
	}
			
	tolerance = GSL_MAX_DBL (this->tol_abs, this->tol_rel * fabs (area));
			
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
			
	this->w->update(a1, b1, area1, error1, a2, b2, area2, error2);
	this->w->retrieve(&a_i,&b_i,&r_i,&e_i);
	this->last_iter++;
			
      } while (this->last_iter < this->w->limit 
	       && !error_type && errsum > tolerance);
		
    result = this->w->sum_results();
    abserr = errsum;
    this->interror = abserr;
		
    if (errsum <= tolerance) {
      return success;
    } else if (error_type == 2) {
      std::string estr="Round-off error prevents tolerance ";
      estr+="from being achieved in inte_qaws_gsl::gsl_qaws().";
      O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);
    } else if (error_type == 3) {
      std::string estr="Bad integrand behavior ";
      estr+=" in inte_qaws_gsl::gsl_qaws().";
      O2SCL_CONV_RET(estr.c_str(),exc_esing,this->err_nonconv);
    } else if (this->last_iter == limit) {
      std::string estr="Maximum number of subdivisions ("+itos(limit);
      estr+=") reached in inte_qaws_gsl::gsl_qaws().";
      O2SCL_CONV_RET(estr.c_str(),exc_emaxiter,this->err_nonconv);
    } else {
      std::string estr="Could not integrate function in ";
      estr+="inte_qaws_gsl::gsl_qaws().";
      O2SCL_ERR(estr.c_str(),exc_efailed);
    }
		
    // No return statement needed since the above if statement
    // always forces a return
  }
		
  /// Return string denoting type ("inte_qaws_gsl")
  const char *type() { return "inte_qaws_gsl"; }
		
  };
	
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
