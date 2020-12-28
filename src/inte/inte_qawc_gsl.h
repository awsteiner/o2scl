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
#ifndef O2SCL_GSL_INTE_QAWC_H
#define O2SCL_GSL_INTE_QAWC_H

/** \file inte_qawc_gsl.h
    \brief File defining \ref o2scl::inte_qawc_gsl
*/

#include <o2scl/err_hnd.h>
#include <o2scl/inte.h>
#include <o2scl/inte_singular_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Chebyshev integration base class (GSL)

      This class provides the basic Chebyshev integration functions
      for use in the GSL-based integration classes which
      require them. See \ref gslinte_subsect in the User's 
      guide for general information about the GSL integration classes.
  */
  template<class func_t> class inte_cheb_gsl : 
    public inte_transform_gsl<func_t> {
    
  protected:
    
    /// Compute the Chebyshev moments
    void compute_moments(double cc, double *moment) {
      size_t k;
      
      double a0 = log (fabs ((1.0 - cc) / (1.0 + cc)));
      double a1 = 2 + a0 * cc;

      moment[0] = a0;
      moment[1] = a1;

      for (k = 2; k < 25; k++) {
	double a2;
	
	if ((k % 2) == 0) {
	  a2 = 2.0 * cc * a1 - a0;
	} else {
	  const double km1 = k - 1.0;
	  a2 = 2.0 * cc * a1 - a0 - 4.0 / (km1 * km1 - 1.0);
	}
	
	moment[k] = a2;
	
	a0 = a1;
	a1 = a2;
      }
    }

    /** \brief Compute Chebyshev series expansion using a FFT method

	The Chebyshev coefficients for the truncated expansions,
	\f[ 
	f(x) = 
	\frac{a_0}{2}T_0(x) + \frac{a_d}{2}T_d(x) + 
	\sum_{k=}^{d-1} a_k^{(d)}T_k(x),
	\f]
	are computed for \f$ d=12 \f$ and \f$ d=24 \f$ using an FFT
	algorithm.
	\verbatim embed:rst
	The FFT algorithm, from [Tolstov62]_, is adapted so that the both
	sets of coefficients are computed simultaneously.
	\endverbatim

	Given the function specified in \c f, this function computes
	the 13 Chebyshev coefficients, \f$ C^{12}_{k} \f$ of degree 12
	and 25 Chebyshev coefficients of degree 24, \f$ C^{24}_{k}
	\f$, for the interval \f$ [a,b] \f$ using a FFT method.
	
	These coefficients are constructed to approximate
	the original function with
	\f[
	f = \sum_{k=1}^{13} C^{12}_{k} T_{k-1}(x)
	\f]
	and 
	\f[
	f = \sum_{k=1}^{25} C^{24}_{k} T_{k-1}(x)
	\f]
	where \f$ T_{k-1}(x) \f$ is the Chebyshev polynomial of
	degree \f$ k-1 \f$ evaluated at the point \f$ x \f$.

	It is assumed that memory for \c cheb12 and \c cheb24 has
	been allocated beforehand.

	Originally written in QUADPACK by R. Piessens and E. de
	Doncker, translated into C for GSL by Brian Gough, and then
	rewritten for \o2.
    */
    template<class func2_t> 
      void inte_cheb_series(func2_t &f, double a, double b, 
			    double *cheb12, double *cheb24) {
      size_t i;
      double fval[25], v[12];
      
      /* These are the values of cos(pi*k/24) for k=1..11 needed for the
	 Chebyshev expansion of f(x) */
      
      const double x[11] = { 0.9914448613738104,     
			     0.9659258262890683,
			     0.9238795325112868,     
			     0.8660254037844386,
			     0.7933533402912352,     
			     0.7071067811865475,
			     0.6087614290087206,     
			     0.5000000000000000,
			     0.3826834323650898,     
			     0.2588190451025208,
			     0.1305261922200516 };
  
      const double center = 0.5 * (b + a);
      const double half_length =  0.5 * (b - a);

      double y1, y2, y3;
      y1=f(b);
      y2=f(center);
      y3=f(a);
      fval[0] = 0.5 * y1;
      fval[12] = y2;
      fval[24] = 0.5 * y3;

      for (i = 1; i < 12; i++) {
	const size_t j = 24 - i;
	const double u = half_length * x[i-1];
	double yp, ym;
	yp=f(center+u);
	ym=f(center-u);
	fval[i] = yp;
	fval[j] = ym;
      }

      for (i = 0; i < 12; i++) {
	const size_t j = 24 - i;
	v[i] = fval[i] - fval[j];
	fval[i] = fval[i] + fval[j];
      }

      {
	const double alam1 = v[0] - v[8];
	const double alam2 = x[5] * (v[2] - v[6] - v[10]);

	cheb12[3] = alam1 + alam2;
	cheb12[9] = alam1 - alam2;
      }

      {
	const double alam1 = v[1] - v[7] - v[9];
	const double alam2 = v[3] - v[5] - v[11];
	{
	  const double alam = x[2] * alam1 + x[8] * alam2;

	  cheb24[3] = cheb12[3] + alam;
	  cheb24[21] = cheb12[3] - alam;
	}

	{
	  const double alam = x[8] * alam1 - x[2] * alam2;
	  cheb24[9] = cheb12[9] + alam;
	  cheb24[15] = cheb12[9] - alam;
	}
      }

      {
	const double part1 = x[3] * v[4];
	const double part2 = x[7] * v[8];
	const double part3 = x[5] * v[6];
    
	{
	  const double alam1 = v[0] + part1 + part2;
	  const double alam2 = x[1] * v[2] + part3 + x[9] * v[10];
      
	  cheb12[1] = alam1 + alam2;
	  cheb12[11] = alam1 - alam2;
	}
    
	{
	  const double alam1 = v[0] - part1 + part2;
	  const double alam2 = x[9] * v[2] - part3 + x[1] * v[10];
	  cheb12[5] = alam1 + alam2;
	  cheb12[7] = alam1 - alam2;
	}
      }

      {
	const double alam = (x[0] * v[1] + x[2] * v[3] + x[4] * v[5]
			     + x[6] * v[7] + x[8] * v[9] + x[10] * v[11]);
	cheb24[1] = cheb12[1] + alam;
	cheb24[23] = cheb12[1] - alam;
      }

      {
	const double alam = (x[10] * v[1] - x[8] * v[3] + x[6] * v[5] 
			     - x[4] * v[7] + x[2] * v[9] - x[0] * v[11]);
	cheb24[11] = cheb12[11] + alam;
	cheb24[13] = cheb12[11] - alam;
      }

      {
	const double alam = (x[4] * v[1] - x[8] * v[3] - x[0] * v[5] 
			     - x[10] * v[7] + x[2] * v[9] + x[6] * v[11]);
	cheb24[5] = cheb12[5] + alam;
	cheb24[19] = cheb12[5] - alam;
      }

      {
	const double alam = (x[6] * v[1] - x[2] * v[3] - x[10] * v[5] 
			     + x[0] * v[7] - x[8] * v[9] - x[4] * v[11]);
	cheb24[7] = cheb12[7] + alam;
	cheb24[17] = cheb12[7] - alam;
      }

      for (i = 0; i < 6; i++) {
	const size_t j = 12 - i;
	v[i] = fval[i] - fval[j];
	fval[i] = fval[i] + fval[j];
      }

      {
	const double alam1 = v[0] + x[7] * v[4];
	const double alam2 = x[3] * v[2];

	cheb12[2] = alam1 + alam2;
	cheb12[10] = alam1 - alam2;
      }

      cheb12[6] = v[0] - v[4];

      {
	const double alam = x[1] * v[1] + x[5] * v[3] + x[9] * v[5];
	cheb24[2] = cheb12[2] + alam;
	cheb24[22] = cheb12[2] - alam;
      }

      {
	const double alam = x[5] * (v[1] - v[3] - v[5]);
	cheb24[6] = cheb12[6] + alam;
	cheb24[18] = cheb12[6] - alam;
      }

      {
	const double alam = x[9] * v[1] - x[5] * v[3] + x[1] * v[5];
	cheb24[10] = cheb12[10] + alam;
	cheb24[14] = cheb12[10] - alam;
      }

      for (i = 0; i < 3; i++) {
	const size_t j = 6 - i;
	v[i] = fval[i] - fval[j];
	fval[i] = fval[i] + fval[j];
      }

      cheb12[4] = v[0] + x[7] * v[2];
      cheb12[8] = fval[0] - x[7] * fval[2];

      {
	const double alam = x[3] * v[1];
	cheb24[4] = cheb12[4] + alam;
	cheb24[20] = cheb12[4] - alam;
      }

      {
	const double alam = x[7] * fval[1] - fval[3];
	cheb24[8] = cheb12[8] + alam;
	cheb24[16] = cheb12[8] - alam;
      }

      cheb12[0] = fval[0] + fval[2];

      {
	const double alam = fval[1] + fval[3];
	cheb24[0] = cheb12[0] + alam;
	cheb24[24] = cheb12[0] - alam;
      }

      cheb12[12] = v[0] - v[2];
      cheb24[12] = cheb12[12];

      for (i = 1; i < 12; i++) {
	cheb12[i] *= 1.0 / 6.0;
      }

      cheb12[0] *= 1.0 / 12.0;
      cheb12[12] *= 1.0 / 12.0;

      for (i = 1; i < 24; i++) {
	cheb24[i] *= 1.0 / 12.0;
      }

      cheb24[0] *= 1.0 / 24.0;
      cheb24[24] *= 1.0 / 24.0;
    }

  };

  /** \brief Adaptive Cauchy principal value integration (GSL)

      The Cauchy principal value of the integral of 
      \f[
      \int_a^b \frac{f(x)}{x-c}~dx =
      \lim_{\epsilon\to 0^+}
      \left\{ \int_a^{c-\epsilon} \frac{f(x)}{x-c}~dx +
      \int_{c+\epsilon}^b \frac{f(x)}{x-c}~dx \right\}.
      \f]
      over \f$ (a,b), \f$ with a singularity at \f$ c, \f$ is
      computed. The adaptive refinement algorithm described for
      inte_qag_gsl is used with modifications to ensure that
      subdivisions do not occur at the singular point \f$ x = c\f$ .
      When a subinterval contains the point \f$ x = c \f$ or is close
      to it, a special 25-point modified Clenshaw-Curtis rule is used
      to control the singularity. Further away from the singularity
      the algorithm uses a Gauss-Kronrod integration rule.
      
      The location of the singularity must be specified before-hand in
      inte_qawc_gsl::s, and the singularity must not be at one of the
      endpoints. Note that when integrating a function of the form \f$
      \frac{f(x)}{(x-s)} \f$, the denominator \f$ (x-s) \f$ must not
      be specified in the argument \c func to integ(). Note that this
      is different from how the \ref inte_cauchy_cern operates.

      See \ref gslinte_subsect in the User's guide for general
      information about the GSL integration classes.

      \future Make inte_cauchy_cern and this class consistent in the
      way which they require the user to provide the denominator
      in the integrand
  */
  template<class func_t> class inte_qawc_gsl : 
  public inte_cheb_gsl<func_t> {
    
  public:

    inte_qawc_gsl() {
    }

    virtual ~inte_qawc_gsl() {}
    
    /// The singularity
    double s;

    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, double a, double b, 
			  double &res, double &err) {
      
      return this->qawc(func,a,b,s,this->tol_abs,this->tol_rel,&res,&err);
    }

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief The full GSL integration routine called by integ_err()
     */
    int qawc(func_t &func, const double a, const double b, const double c,
	     const double epsabs, const double epsrel, 
	     double *result, double *abserr) {
      
      double area, errsum;
      double result0, abserr0;
      double tolerance;
      size_t iteration = 0;
      int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;
      int err_reliable;
      int sign = 1;
      double lower, higher;

      /* Initialize results */

      *result = 0;
      *abserr = 0;

      size_t limit=this->w->limit;

      if (b < a)  {
	lower = b; 
	higher = a;
	sign = -1;
      } else {
	lower = a;
	higher = b;
      }

      this->w->initialise(lower,higher);

      double dbl_eps=std::numeric_limits<double>::epsilon();
      
      if (epsabs <= 0 && (epsrel < 50 * dbl_eps || epsrel < 0.5e-28)) {
	this->last_iter=0;
	std::string estr="Tolerance cannot be achieved with given ";
	estr+="value of tol_abs, "+dtos(epsabs)+", and tol_rel, "+
	  dtos(epsrel)+", in inte_qawc_gsl::qawc().";
	O2SCL_ERR(estr.c_str(),exc_ebadtol);
      }

      if (c == a || c == b) {
	this->last_iter=0;
	std::string estr="Cannot integrate with singularity on endpoint ";
	estr+="in inte_qawc_gsl::qawc().";
	O2SCL_ERR(estr.c_str(),exc_einval);
      }      

      /* perform the first integration */
      
      this->qc25c(func,lower,higher,c,&result0,&abserr0, &err_reliable);

      this->w->set_initial_result (result0, abserr0);

      /* Test on accuracy, use 0.01 relative error as an extra safety
	 margin on the first iteration (ignored for subsequent iterations) 
      */
      tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (result0));

      if (abserr0 < tolerance && abserr0 < 0.01 * fabs(result0)) {

	this->last_iter=1;
	*result = sign * result0;
	*abserr = abserr0;
	return success;

      } else if (limit == 1) {

	*result = sign * result0;
	*abserr = abserr0;

	this->last_iter=1;
	std::string estr="A maximum of 1 iteration was insufficient ";
	estr+="in inte_qawc_gsl::qawc().";
	O2SCL_CONV_RET(estr.c_str(),exc_emaxiter,this->err_nonconv);
      }

      area = result0;
      errsum = abserr0;

      iteration = 1;

      do {

	double a1, b1, a2, b2;
	double a_i, b_i, r_i, e_i;
	double area1 = 0, area2 = 0, area12 = 0;
	double error1 = 0, error2 = 0, error12 = 0;
	int err_reliable1, err_reliable2;

	/* Bisect the subinterval with the largest error estimate */

	this->w->retrieve (&a_i, &b_i, &r_i, &e_i);

	a1 = a_i; 
	b1 = 0.5 * (a_i + b_i);
	a2 = b1;
	b2 = b_i;

	if (c > a1 && c <= b1) {
	  b1 = 0.5 * (c + b2) ;
	  a2 = b1;
	} else if (c > b1 && c < b2) {
	  b1 = 0.5 * (a1 + c) ;
	  a2 = b1;
	}

	qc25c (func, a1, b1, c, &area1, &error1, &err_reliable1);
	qc25c (func, a2, b2, c, &area2, &error2, &err_reliable2);

	area12 = area1 + area2;
	error12 = error1 + error2;

	errsum += (error12 - e_i);
	area += area12 - r_i;

	if (err_reliable1 && err_reliable2) {
	  double delta = r_i - area12;

	  if (fabs (delta) <= 1.0e-5 * fabs (area12) && 
	      error12 >= 0.99 * e_i) {
	    roundoff_type1++;
	  }
	  if (iteration >= 10 && error12 > e_i) {
	    roundoff_type2++;
	  }
	}

	tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area));

	if (errsum > tolerance) {
	  if (roundoff_type1 >= 6 || roundoff_type2 >= 20) {
	    error_type = 2;   /* round off error */
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
	  std::cout << "inte_qawc_gsl Iter: " << iteration;
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

      } while (iteration < limit && !error_type && errsum > tolerance);

      *result = sign * this->w->sum_results();
      *abserr = errsum;

      this->last_iter=iteration;
      if (errsum <= tolerance) {
	return success;
      } else if (error_type == 2) {
	std::string estr="Roundoff error prevents tolerance ";
	estr+="from being achieved in inte_qawc_gsl::qawc().";
	O2SCL_CONV_RET(estr.c_str(),exc_eround,this->err_nonconv);
      } else if (error_type == 3) {
	std::string estr="Bad integrand behavior ";
	estr+=" in inte_qawc_gsl::qawc().";
	O2SCL_CONV_RET(estr.c_str(),exc_esing,this->err_nonconv);
      } else if (iteration == limit) {
	std::string estr="Maximum number of subdivisions ("+itos(iteration);
	estr+=") reached in inte_qawc_gsl::qawc().";
	O2SCL_CONV_RET(estr.c_str(),exc_emaxiter,this->err_nonconv);
      } else {
	std::string estr="Could not integrate function in inte_qawc_gsl::";
	estr+="qawc() (it may have returned a non-finite result).";
	O2SCL_ERR(estr.c_str(),exc_efailed);
      }

      // No return statement needed since the above if statement
      // always forces a return, but some compilers like having one
      // anyway.
      return o2scl::success;
    }

    /// 25-point quadrature for Cauchy principal values
    void qc25c(func_t &func, double a, double b, double c, 
	       double *result, double *abserr, int *err_reliable) {

      double cc = (2 * c - b - a) / (b - a);
      
      if (fabs (cc) > 1.1) {
	double resabs, resasc;
	    
	this->gauss_kronrod(func,a,b,result,abserr,&resabs,&resasc);
      
	if (*abserr == resasc) {
	  *err_reliable = 0;
	} else {
	  *err_reliable = 1;
	}

	return;

      } else {

	double cheb12[13], cheb24[25], moment[25];
	double res12 = 0, res24 = 0;
	size_t i;
	this->inte_cheb_series(func, a, b, cheb12, cheb24);
	this->compute_moments (cc, moment);
	  
	for (i = 0; i < 13; i++) {
	  res12 += cheb12[i] * moment[i];
	}
	  
	for (i = 0; i < 25; i++) {
	  res24 += cheb24[i] * moment[i];
	}
	  
	*result = res24;
	*abserr = fabs(res24 - res12) ;
	*err_reliable = 0;

	return;
      }
    }

    /// Add the singularity to the function
    virtual double transform(double t, func_t &func) {
      double y;
      y=func(t);
      return y/(t-s);
    }

#endif
  
    /// Return string denoting type ("inte_qawc_gsl")
    const char *type() { return "inte_qawc_gsl"; }
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
