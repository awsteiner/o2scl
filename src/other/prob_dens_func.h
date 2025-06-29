/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2012-2025, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
/** \file prob_dens_func.h
    \brief File for probability density functions
*/
#ifndef O2SCL_PROB_DENS_FUNC_H
#define O2SCL_PROB_DENS_FUNC_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

// For solving quadratics in bivariate gaussians
#include <gsl/gsl_poly.h>

#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <o2scl/columnify.h>
#include <o2scl/rng.h>
#include <o2scl/search_vec.h>
#include <o2scl/cholesky.h>
#include <o2scl/lu.h>
#include <o2scl/vec_stats.h>
#include <o2scl/constants.h>

namespace o2scl {

  /** \brief A one-dimensional probability density function

      This class is experimental.

      \verbatim embed:rst
      .. todo:: 

         Future: Give functions for mean, median, mode, variance, etc?

      \endverbatim

      \comment
      For now, there aren't any pure virtual functions,
      since this causes problems in creating an
      std::vector<prob_dens_func> object below (especially
      with intel compilers)
      \endcomment
  */
  class prob_dens_func {
    
  public:

    virtual ~prob_dens_func() {
    }
    
    /// Sample from the specified density
    virtual double operator()() const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0.0;
    }

    /// Sample from the specified density
    virtual double sample() const {
      return (*this)();
    }
    
    /// The normalized density 
    virtual double pdf(double x) const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0.0;
    }
    
    /// The log of the normalized density 
    virtual double log_pdf(double x) const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0.0;
    }
    
    /// The cumulative distribution function (from the lower tail)
    virtual double cdf(double x) const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0.0;
    }
    
    /// The inverse cumulative distribution function
    virtual double invert_cdf(double cdf) const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0.0;
    }

    /// Entropy of the distribution (\f$ - \int f \ln f \f$ )
    virtual double entropy() const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0.0;
    }
    
  };
  
  /** \brief A one-dimensional Gaussian probability density

      The distribution
      \f[
      P(x)=\frac{1}{\sigma \sqrt{2 \pi}} 
      e^{-\frac{\left(x-x_0\right)^2}{2\sigma^2}}
      \f]

      This class is experimental.
  */
  class prob_dens_gaussian : public prob_dens_func {
    
  protected:
    
    /** \brief Central value
     */
    double cent_;

    /** \brief Width parameter 

	A value of -1 indicates it is yet unspecified.
    */
    double sigma_;

    /// Base random number generator
    rng<> r2;
    
    /// C++ base normal distribution
    mutable std::normal_distribution<double> nd;
    
  public:
    
    /** \brief Create a standard normal distribution
     */
    prob_dens_gaussian() : nd(0.0,1.0) {
      cent_=0.0;
      sigma_=1.0;
    }

    /** \brief Create a Gaussian distribution with width \c sigma

	The value of \c sigma must be larger than zero.
    */
    prob_dens_gaussian(double cent, double sigma) : nd(cent,sigma) {
      if (sigma<0.0) {
	O2SCL_ERR2("Tried to create a Gaussian dist. with sigma",
		   "<0 in prob_dens_gaussian::prob_dens_gaussian().",
		   exc_einval);
      }
      cent_=cent;
      sigma_=sigma;
    }
    
    virtual ~prob_dens_gaussian() {
    }

    /// Copy constructor
    prob_dens_gaussian(const prob_dens_gaussian &pdg) : prob_dens_func() {
      cent_=pdg.cent_;
      sigma_=pdg.sigma_;
      r2=pdg.r2;
      return;
    }

    /// Copy constructor with operator=
    prob_dens_gaussian &operator=(const prob_dens_gaussian &pdg) {
      // Check for self-assignment
      if (this!=&pdg) {
	cent_=pdg.cent_;
	sigma_=pdg.sigma_;
	r2=pdg.r2;
      }
      return *this;
    }

    /// Set the seed
    void set_seed(unsigned long int s) {
      r2.set_seed(s);
      return;
    }

    /// Set the center
    void set_center(double cent) {
      cent_=cent;
      nd=std::normal_distribution<double>(cent,sigma_);
      return;
    }

    /// Set the Gaussian width (must be positive)
    void set_sigma(double sigma) {
      if (sigma<0.0) {
	O2SCL_ERR2("Tried to set negative sigma ",
		   "in prob_dens_gaussian::prob_dens_gaussian().",
		   exc_einval);
      }
      sigma_=sigma;
      nd=std::normal_distribution<double>(cent_,sigma);
      return;
    }

    /// Get the center
    double mean() {
      return cent_;
    }

    /// Get the Gaussian width
    double stddev() {
      if (sigma_<0.0) {
	O2SCL_ERR2("Width not set in prob_dens_gaussian::",
		   "get_sigma().",exc_einval);
      }
      return sigma_;
    }

    /// Sample from the specified density
    virtual double operator()() const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Width not set in prob_dens_gaussian::",
		   "operator().",exc_einval);
      }
      return nd(r2.def_engine);
    }
    
    /// The normalized density 
    virtual double pdf(double x) const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Width not set in prob_dens_gaussian::",
		   "pdf().",exc_einval);
      }
      return gsl_ran_gaussian_pdf(x-cent_,sigma_);
    }

    /// The log of the normalized density 
    virtual double log_pdf(double x) const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Width not set in prob_dens_gaussian::",
		   "pdf().",exc_einval);
      }
      return log(gsl_ran_gaussian_pdf(x-cent_,sigma_));
    }

    /// The cumulative distribution function (from the lower tail)
    virtual double cdf(double x) const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Width not set in prob_dens_gaussian::",
		   "cdf().",exc_einval);
      }
      return gsl_cdf_gaussian_P(x-cent_,sigma_);
    }
    
    /// The inverse cumulative distribution function
    virtual double invert_cdf(double in_cdf) const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Width not set in prob_dens_gaussian::",
		   "invert_cdf().",exc_einval);
      }
      if (in_cdf<0.0 || in_cdf>1.0) {
	O2SCL_ERR((((std::string)"Requested cdf inverse of ")+
                   o2scl::dtos(in_cdf)+" which is outside of [0,1] "+
                   "in prob_dens_gaussian::invert_cdf().").c_str(),
                  exc_einval);
      }
      return gsl_cdf_gaussian_Pinv(in_cdf,sigma_)+cent_;
    }

    /// Entropy of the distribution (\f$ - \int f \ln f \f$ )
    virtual double entropy() const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Width not set in prob_dens_gaussian::",
		   "invert_cdf().",exc_einval);
      }
      return log(2.0*o2scl_const::pi*exp(1.0)*sigma_*sigma_);
    }
    
  };

  /** \brief A one-dimensional probability density over
      a finite range

      This class is experimental.
  */
  class prob_dens_frange : public prob_dens_func {

  public:

    virtual ~prob_dens_frange() {
    }
    
    /// Lower limit of the range
    virtual double lower_limit() const=0;

    /// Uower limit of the range
    virtual double upper_limit() const=0;

  };
  
  /** \brief A uniform one-dimensional probability density
      over a finite range

      A flat distribution given by \f$ P(x)=1/(b-a) \f$ for \f$ a<x<b
      \f$, where \f$ a \f$ is the lower limit and \f$ b \f$ is the
      upper limit.

      This class is experimental.
  */
  class prob_dens_uniform : public prob_dens_frange {
    
  protected:

    /// Lower limit
    double ll;

    /// Upper limit
    double ul;

    /// The GSL random number generator
    rng<> r;
    
  public:

    /** \brief Create a blank uniform distribution
     */
    prob_dens_uniform() {
      ll=1.0;
      ul=0.0;
    }

    /** \brief Create a uniform distribution from \f$ a<x<b \f$ 
     */
    prob_dens_uniform(double a, double b) {
      // Ensure a<b
      if (a>b) {
	double tmp=a;
	a=b;
	b=tmp;
      }
      ll=a;
      ul=b;
    }

    virtual ~prob_dens_uniform() {
    }
    
    /// Copy constructor
    prob_dens_uniform(const prob_dens_uniform &pdg) : prob_dens_frange() {
      ll=pdg.ll;
      ul=pdg.ul;
    }

    /// Copy constructor with operator=
    prob_dens_uniform &operator=(const prob_dens_uniform &pdg) {
      // Check for self-assignment
      if (this==&pdg) return *this;
      ll=pdg.ll;
      ul=pdg.ul;
      return *this;
    }

    /// Set the seed
    void set_seed(unsigned long int s) {
      r.set_seed(s);
      return;
    }

    /** \brief Set the limits of the uniform distribution
     */
    void set_limits(double a, double b) {
      // Ensure a<b
      if (a>b) {
	double tmp=a;
	a=b;
	b=tmp;
      }
      ll=a;
      ul=b;
      return;
    }
    
    /// Lower limit of the range
    virtual double lower_limit() const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "lower_limit().",exc_einval);
      }
      return ll;
    }

    /// Uower limit of the range
    virtual double upper_limit() const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "upper_limit().",exc_einval);
      }
      return ul;
    }

    /// Operator from the specified density
    virtual double operator()() const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "operator().",exc_einval);
      }
      double ret=r.random()*(ul-ll)+ll;
      return ret;
    }
    
    /// The normalized density 
    virtual double pdf(double x) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "pdf().",exc_einval);
      }
      if (x<ll || x>ul) return 0.0;
      return 1.0/(ul-ll);
    }
    
    /// The log of the normalized density 
    virtual double log_pdf(double x) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "pdf().",exc_einval);
      }
      if (x<ll || x>ul) return 0.0;
      return log(1.0/(ul-ll));
    }
    
    /// The cumulative distribution function (from the lower tail)
    virtual double cdf(double x) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "cdf().",exc_einval);
      }
      if (x<ll) return 0.0;
      if (x>ul) return 1.0;
      return (x-ll)/(ul-ll);
    }
    
    /// The inverse cumulative distribution function
    virtual double invert_cdf(double in_cdf) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "invert_cdf().",exc_einval);
      }
      if (in_cdf<0.0 || in_cdf>1.0) {
	O2SCL_ERR((((std::string)"Requested cdf inverse of ")+
                   o2scl::dtos(in_cdf)+" which is outside of [0,1] "+
                   "in prob_dens_gaussian::invert_cdf().").c_str(),
                  exc_einval);
      }
      if (in_cdf==1.0) return ul;
      if (in_cdf==0.0) return ll;
      return in_cdf*(ul-ll)+ll;
    }

    /// Entropy of the distribution (\f$ - \int f \ln f \f$ )
    virtual double entropy() const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "entropy().",exc_einval);
      }
      return log(ul-ll);
    }

  };

  /** \brief A one-dimensional probability density over the
      positive real numbers

      This class is experimental.
  */
  class prob_dens_positive : public prob_dens_func {

  public:
    
    virtual ~prob_dens_positive() {
    }
    
  };

  /** \brief Lognormal density function

      The distribution
      \f[
      P(x)=\frac{1}{x \sigma \sqrt{2 \pi}} 
      \exp \left[-\frac{\left(\ln x-\mu\right)^2}{2\sigma^2}\right]
      \f]

      This class is experimental.
  */
  class prob_dens_lognormal : public prob_dens_positive {

  protected:
    
    /** \brief Width parameter 

	A value of -1 indicates it is yet unspecified.
    */
    double sigma_;
    
    /** \brief Central value

	A value of -1 indicates it is yet unspecified.
    */
    double mu_;

    /// The GSL random number generator
    rng<> r;

    /// C++ base normal distribution
    std::lognormal_distribution<double> lnd;
    
  public:

    /** \brief Create a blank lognormal distribution
     */
    prob_dens_lognormal() : lnd() {
      sigma_=-1.0;
      mu_=0.0;
    }
    
    /** \brief Create lognormal distribution with mean parameter \c mu
	and width parameter \c sigma

	The value of \c sigma must be larger than zero.
    */
    prob_dens_lognormal(double mu, double sigma) : lnd(mu,sigma) {
      if (sigma<0.0) {
	O2SCL_ERR2("Tried to create log normal dist. with mu or sigma",
		   "<0 in prob_dens_lognormal::prob_dens_lognormal().",
		   exc_einval);
      }
      mu_=mu;
      sigma_=sigma;
    }

    virtual ~prob_dens_lognormal() {
    }

    /// Copy constructor
    prob_dens_lognormal(const prob_dens_lognormal &pdg) : prob_dens_positive() {
      mu_=pdg.mu_;
      sigma_=pdg.sigma_;
      lnd=std::lognormal_distribution<double>(mu_,sigma_);
    }

    /// Copy constructor with operator=
    prob_dens_lognormal &operator=(const prob_dens_lognormal &pdg) {
      // Check for self-assignment
      if (this==&pdg) return *this;
      mu_=pdg.mu_;
      sigma_=pdg.sigma_;
      lnd=pdg.lnd;
      return *this;
    }

    /** \brief Set the maximum and width of the lognormal distribution
     */
    void set_mu_sigma(double mu, double sigma) {
      if (sigma<0.0) {
	O2SCL_ERR2("Tried to set mu or sigma negative",
		   "in prob_dens_lognormal::prob_dens_lognormal().",
		   exc_einval);
      }
      mu_=mu;
      sigma_=sigma;
      lnd=std::lognormal_distribution<double>(mu_,sigma_);
      return;
    }

    /// Set the seed
    void set_seed(unsigned long int s) {
      r.set_seed(s);
      return;
    }

    /// Sample from the specified density
    virtual double operator()() const {
      double u, r2=0.0;
      do {
        u=-1.0+2.0*r.random();
        double v=-1.0+2.0*r.random();
        r2=u*u+v*v;
      } while (r2>1.0 || r2==0.0);
      double normal=u*sqrt(-2.0*log(r2)/r2);
      double z=exp(sigma_*normal+mu_);
      return z;
    }
    
    /// The normalized density 
    virtual double pdf(double x) const {
      if (x<0.0) {
	return 0.0;
      }
      double u=(log(x)-mu_)/sigma_;
      double p=1.0/(x*fabs(sigma_)*sqrt(2.0*o2scl_const::pi))*
        exp(-u*u/2.0);

      return p;
    }
    
    /// The log of the normalized density 
    virtual double log_pdf(double x) const {
      if (x<0.0) {
	return 0.0;
      }
      return log(this->pdf(x));
    }
    
    /// The cumulative distribution function (from the lower tail)
    virtual double cdf(double x) const {
      if (x<0.0) {
	return 0.0;
      }
      return gsl_cdf_lognormal_P(x,mu_,sigma_);
    }
    
    /// The inverse cumulative distribution function
    virtual double invert_cdf(double in_cdf) const {
      if (in_cdf<0.0 || in_cdf>1.0) {
	O2SCL_ERR((((std::string)"Requested cdf inverse of ")+
                   o2scl::dtos(in_cdf)+" which is outside of [0,1] "+
                   "in prob_dens_gaussian::invert_cdf().").c_str(),
                  exc_einval);
      }
      return gsl_cdf_lognormal_Pinv(in_cdf,mu_,sigma_);
    }

    /// Entropy of the distribution (\f$ - \int f \ln f \f$ )
    virtual double entropy() const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Parameters not set in prob_dens_lognormal::",
		   "entropy().",exc_einval);
      }
      return 0.5+0.5*log(2.0*o2scl_const::pi*sigma_*sigma_)+mu_;
    }
    
  };
  
  /** \brief A multi-dimensional probability density function
      
      \note This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class prob_dens_mdim {
    
  public:

    prob_dens_mdim() {
      verbose=0;
    }
    
    virtual ~prob_dens_mdim() {
    }
    
    /** \brief Verbosity parameter
     */
    int verbose;
    
    /// Return the dimensionality
    virtual size_t dim() const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0;
    }
  
    /// The normalized density 
    virtual double pdf(const vec_t &x) const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0.0;
    }
  
    /// The log of the normalized density 
    virtual double log_pdf(const vec_t &x) const {
      double val=pdf(x);
      if (!std::isfinite(val) || val<0.0) {
        O2SCL_ERR2("PDF not finite or negative in ",
                   "prob_dens_mdim::log_pdf().",o2scl::exc_efailed);
      }
      double val2=log(pdf(x));
      if (!std::isfinite(val2)) {
        std::cout << val << " " << val2 << std::endl;
        O2SCL_ERR2("Log of PDF not finite in ",
                   "prob_dens_mdim::log_pdf().",o2scl::exc_efailed);
      }
      return val2;
    }
  
    /// Sample the distribution
    virtual void operator()(vec_t &x) const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return;
    }

  };

  /** \brief A multidimensional distribution formed by the product
      of several one-dimensional distributions
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class prob_dens_mdim_factor : public prob_dens_mdim<vec_t> {
    
  protected:
    
    /// Vector of one-dimensional distributions
    std::vector<prob_dens_func> list;
    
  public:
    
    prob_dens_mdim_factor(std::vector<prob_dens_func> &p_list) {
      list=p_list;
    }
  
    /// Copy constructor
    prob_dens_mdim_factor(const prob_dens_mdim_factor &pdmf) {
      list=pdmf.list;
    }
  
    /// Copy constructor with operator=
    prob_dens_mdim_factor &operator=
    (const prob_dens_mdim_factor &pdmf) {
      // Check for self-assignment
      if (this!=&pdmf) {
        list=pdmf.list;
      }
      return *this;
    }
  
    /// Return the dimensionality
    virtual size_t dim() const {
      return list.size();
    }
  
    /// The normalized density 
    virtual double pdf(const vec_t &x) const {
      double ret=1.0;
      for(size_t i=0;i<list.size();i++) ret*=list[i].pdf(x[i]);
      return ret;
    }

    /// The log of the normalized density 
    virtual double log_pdf(const vec_t &x) const {
      double ret=0.0;
      for(size_t i=0;i<list.size();i++) ret+=list[i].log_pdf(x[i]);
      return ret;
    }
    
    /// Sample the distribution
    virtual void operator()(vec_t &x) const {
      for(size_t i=0;i<list.size();i++) x[i]=list[i]();
      return;
    }
  
  };

  /** \brief A bivariate gaussian probability distribution

      For a two-dimensional gaussian, given a mean 
      \f$ ( \mu_x, \mu_y ) \f$ and a covariance matrix
      \f[
      \Sigma = \left( 
      \begin{array}{cc}
      \sigma_x^2 & \rho \sigma_x \sigma_y \\
      \rho \sigma_x \sigma_y & \sigma_y^2 \\
      \end{array}
      \right)
      \f]
      the PDF is
      \f[
      pdf(x,y) = \left(2 \pi \sigma_x \sigma_y \sqrt{1-\rho^2}\right)^{-1}
      \exp \left\{ - \frac{1}{2 (1-\rho^2)} 
      \left[ \frac{(x-\mu_x)^2}{\sigma_x^2} + 
      \frac{(y-\mu_y)^2}{\sigma_y^2} - 
      \frac{2 \rho (x-\mu_x)(y-\mu_y)}{\sigma_x \sigma_y} \right]
      \right\}
      \f]
      (taken from the Wikipedia page on the "Multivariate normal
      distribution").
      
      The function \ref o2scl::prob_dens_mdim_biv_gaussian::contour()
      gives a point on the contour line for a fixed value of the
      PDF given an angle \f$ \theta \f$. In particular,
      it solves 
      \f[
      c = \mathrm{pdf}(r \cos \theta + \mu_x, r \sin \theta + \mu_y )
      \f]
      for the radius \f$ r \f$ and then stores the values 
      \f$ r \cos \theta + \mu_x \f$ and \f$ r \sin \theta + \mu_y \f$
      in the reference parameters named \c x and \c y . Thus
      this function can be used to map out the full contour
      by selecting values for \f$ \theta \in [0,2 \pi] \f$.
      
      The function \ref
      o2scl::prob_dens_mdim_biv_gaussian::level_fixed_integral() gives
      the value of the PDF for which the integral inside the
      corresponding contour is some fraction of the total integral
      (which is always 1). Given a fraction \f$ f \f$, the argument of
      the exponential is related to the inverse of the cumulative
      distribution function for the chi-squared probability
      distribution for two degrees of freedom for \f$ 1-f \f$. For a
      fraction \f$ f \f$, the value \f$ \chi^2 \f$ (i.e. the
      Mahalanobis distance) is \f$ \chi^2 = -2 \log (1-f) \f$ and then
      the value of the PDF for the corresponding contour is \f$
      pdf(x,y) = \left(2 \pi \sigma_x \sigma_y
      \sqrt{1-\rho^2}\right)^{-1} \exp (-\chi^2/2) \f$ .
      
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class prob_dens_mdim_biv_gaussian : public prob_dens_mdim<vec_t> {

  private:

    /// The x coordinate of the centroid
    double x0;

    /// The y coordinate of the centroid
    double y0;

    /// The x standard deviation
    double sig_x;

    /// The y standard deviation
    double sig_y;

    /// The covariance (always between -1 and 1)
    double rho;
  
  public:
  
    prob_dens_mdim_biv_gaussian() {
    }
  
    virtual ~prob_dens_mdim_biv_gaussian() {
    }
  
    /// Copy constructor
    prob_dens_mdim_biv_gaussian(const prob_dens_mdim_biv_gaussian &pdmbg) {
      x0=pdmbg.x0;
      y0=pdmbg.y0;
      sig_x=pdmbg.sig_x;
      sig_y=pdmbg.sig_y;
      rho=pdmbg.rho;
    }
  
    /// Copy constructor with operator=
    prob_dens_mdim_biv_gaussian &operator=
    (const prob_dens_mdim_biv_gaussian &pdmbg) {
      // Check for self-assignment
      if (this!=&pdmbg) {
        x0=pdmbg.x0;
        y0=pdmbg.y0;
        sig_x=pdmbg.sig_x;
        sig_y=pdmbg.sig_y;
        rho=pdmbg.rho;
      }
      return *this;
    }
  
    /** \brief Set the properties of the distribution

        \note If \f$ |\rho|\geq 1 \f$ this function will
        call the error handler.
    */
    void set(double x_cent, double y_cent, double x_std, double y_std,
             double covar) {
      if (fabs(covar)>=1.0) {
        O2SCL_ERR2("Covariance cannot have magnitude equal or larger than ",
                   "1 in prob_dens_mdim_biv_gaussian::set().",
                   o2scl::exc_einval);
      }
      x0=x_cent;
      y0=y_cent;
      sig_x=x_std;
      sig_y=y_std;
      rho=covar;
      return;
    }

    /** \brief Get the properties of the distribution
     */
    void get(double &x_cent, double &y_cent, double &x_std, double &y_std,
             double &covar) {
      
      x_cent=x0;
      y_cent=y0;
      x_std=sig_x;
      y_std=sig_y;
      covar=rho;
      
      return;
    }

    /** \brief Compute the normalized probability density
     */
    virtual double pdf(const vec_t &v) const {
      double x=v[0], y=v[1];
      double arg=-((x-x0)*(x-x0)/sig_x/sig_x+
                   (y-y0)*(y-y0)/sig_y/sig_y-
                   2.0*rho*(x-x0)*(y-y0)/sig_x/sig_y)/2.0/(1.0-rho*rho);
      double ret=exp(arg)/2.0/o2scl_const::pi/sig_x/sig_y/
        sqrt(1.0-rho*rho);
      return ret;
    }
  
    /** \brief Return the contour level corresponding to a fixed
        integral

        This function returns the value of the PDF (as returned by
        pdf()) for which the integral inside the contour line for that
        value is equal to the specified integral. 

        For example, in the standard two-dimensional normal case, if
        \f[
        \int_0^{a}~r dr~\exp\left( - \frac{r^2}{2} \right) 
        \f]
        is equal to the value given in \c integral, then this 
        function returns the value
        \f[
        \frac{1}{2 \pi} \exp \left( - \frac{a^2}{2} \right) \, .
        \f]
    */
    virtual double level_fixed_integral(double integral) {
      // This comes from inverting the cumulative distribution function
      // for the chi-squared distribution for two degrees of of freedom,
      // i.e. exp(-x/2)
      double arg=-2.0*log(1.0-integral);
      // Now compute the pdf for the fixed value of the
      // squared Mahalanobis distance
      return exp(-0.5*arg)/2.0/o2scl_const::pi/sig_x/
        sig_y/sqrt(1.0-rho*rho);
    }
  
    /** \brief Return the properties of an ellipse based on the 
        integral
    */
    void ellipse_frac_integral(double integral, double &x_cent,
                               double &y_cent, double &x_wid,
                               double &y_wid, double &angle) {
                               
      if (integral<=0.0 || integral>=1.0) {
        O2SCL_ERR("Invalid fraction in ellipse_frac_integral().",
                  o2scl::exc_einval);
      }
      x_cent=x0;
      y_cent=y0;
      // Compute the level for the specified integral
      double level=level_fixed_integral(integral);
      // Compute the associated argument to the exponential,
      // this is (x-x0)^2/sig_x^2+(y-y0)^2/sig_y^2 
      double arg=-2.0*log(level*sig_x*sig_y*sqrt(1.0-rho*rho));
      x_wid=sig_x*sqrt(arg);
      y_wid=sig_y*sqrt(arg);
      angle=atan(rho);
      return;
    }
    
    /** \brief Return a point on the contour for a specified level
        given an angle
    */
    virtual void contour(double level, double theta, vec_t &x) {
      if (level<0.0) {
        O2SCL_ERR2("Cannot produce contours for negative values in ",
                   "prob_dens_mdim_biv_gaussian::contour().",
                   o2scl::exc_einval);
      }
      double max=0.5/sig_x/sig_y/o2scl_const::pi/sqrt(1.0-rho*rho);
      if (level>max) {
        O2SCL_ERR2("Cannot produce contours larger than maximum in ",
                   "prob_dens_mdim_biv_gaussian::contour().",
                   o2scl::exc_einval);
      }
      double arg=-log(level*2.0*o2scl_const::pi*sig_x*sig_y*
                      sqrt(1.0-rho*rho))*2.0*(1.0-rho*rho);
      double r2=arg/(cos(theta)*cos(theta)/sig_x/sig_x+
                     sin(theta)*sin(theta)/sig_y/sig_y-
                     2.0*rho/sig_x/sig_y*cos(theta)*sin(theta));
      x[0]=sqrt(r2)*cos(theta)+x0;
      x[1]=sqrt(r2)*sin(theta)+y0;
      return;
    }
  
  };
  
  /** \brief A multi-dimensional Gaussian probability density function
      using a Cholesky decomposition

      Given a (square) covariance matrix, \f$ \Sigma \f$, and a mean
      vector \f$ \mu \f$ the PDF is
      \f[
      P(x) = \det \left( 2 \pi \Sigma \right)^{-1/2}
      \exp \left[ -\frac{1}{2} (x-\mu)^T \Sigma^{-1} (x-\mu) \right]
      \f]
      
      Given the Cholesky decomposition \f$ A A^{T} = \Sigma \f$,
      and a vector, \f$ z \f$ of samples from the standard Gaussian
      with 0 mean and unit variance, one can create a sample 
      \f$ x \f$ from \f$ x = \mu + A z \f$ .

      \note This class inverts the matrix, necessary for computing the
      pdf, but not for sampling the distribution, so for large
      matrices the inversion can be a waste of computation if the pdf
      is not needed.

      A separate class for the two-dimensional case is in \ref
      prob_dens_mdim_biv_gaussian .
      
      \note Note that, for example, a LU decomposition does not work
      for this class because a Cholesky decomposition (or a spectral
      decomposition) is required for sampling. For this reason, we
      cannot use a generic \ref o2scl_linalg::matrix_invert_det
      object. However, we still could use Cholesky decompositions from
      armadillo or Eigen.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_t=boost::numeric::ublas::matrix<double> >
  class prob_dens_mdim_gaussian : public prob_dens_mdim<vec_t> {
    
  protected:

    /// Cholesky decomposition
    mat_t chol;

    /// Inverse of the covariance matrix
    mat_t covar_inv;

    /// Location of the peak
    vec_t peak;

    /// Normalization factor, \f$ \det ( 2 \pi \Sigma)^{-1/2} \f$
    double norm;

    /// Number of dimensions
    size_t ndim;

  public:

    /// Verbosity parameter (default 0)
    int verbose;
    
    /** \brief Read the distribution from an input file
     */
    virtual int read_generic(std::istream &fin) {

      std::string stemp;
      
      // Read first line and into list
      fin >> ndim;
      if (verbose>1) {
        std::cout << "prob_dens_mdim_gaussian::read_generic():" << std::endl;
        std::cout << "  ndim: " << ndim << std::endl;
      }
      peak.resize(ndim);
      std::cout << "  peak: ";
      for(size_t i=0;i<ndim;i++) {
        fin >> peak[i];
        std::cout << peak[i] << " ";
      }
      std::cout << std::endl;
      fin >> stemp;
      
      mat_t covar(ndim,ndim);
      
      if (stemp[0]=='s' || stemp[0]=='S') {
        vec_t stds(ndim);
        std::cout << "  width: ";
        for(size_t i=0;i<ndim;i++) {
          fin >> stds[i];
          std::cout << stds[i] << " ";
        }
        std::cout << std::endl;
        for(size_t i=0;i<ndim;i++) {
          for(size_t j=0;j<=i;j++) {
            if (i==j) {
              covar(i,j)=stds[i]*stds[i];
            } else {
              fin >> covar(i,j);
              covar(i,j)*=stds[i]*stds[j];
              covar(j,i)=covar(i,j);
            }
          }
        }
      } else if (stemp[0]=='c' || stemp[0]=='C') {
        covar(0,0)=o2scl::stod(stemp);
        for(size_t i=0;i<ndim;i++) {
          for(size_t j=0;j<=i;j++) {
            if (i>0) {
              fin >> covar(i,j);
            }
            covar(j,i)=covar(i,j);
          }
        }
      } else {
        norm=o2scl::stod(stemp);
        for(size_t i=0;i<ndim;i++) {
          for(size_t j=0;j<=i;j++) {
            if (i>0) {
              fin >> chol(i,j);
            }
            chol(j,i)=chol(i,j);
          }
        }
        for(size_t i=0;i<ndim;i++) {
          for(size_t j=0;j<=i;j++) {
            if (i>0) {
              fin >> covar_inv(i,j);
            }
            covar_inv(j,i)=covar_inv(i,j);
          }
        }
      }

      set_covar(ndim,peak,covar);
      
      return 0;
    }

    /** \brief Write the Gaussian to a generic output file
     */
    virtual int write_generic(std::ostream &fout) {
      fout << ndim << std::endl;
      o2scl::vector_out(fout,peak,true);
      fout << norm << std::endl;
      o2scl::matrix_out(fout,ndim,ndim,chol);
      o2scl::matrix_out(fout,ndim,ndim,covar_inv);
      return 0;
    }
    
    /** \brief Standard normal
        \comment
        This has to be public so the user can set the random seed,
        or we have to create a new set_seed() function.
        \endcomment
    */
    o2scl::prob_dens_gaussian pdg;

    /** \brief Get the Cholesky decomposition 
     */
    const mat_t &get_chol() const {
      return chol;
    }

    /** \brief Get the covariance matrix

        The covariance matrix is not stored. This function computes it
        directly from the Cholesky decomposition.
     */
    template<class mat2_t> void get_covar(mat2_t &m) const {
      mat_t chol2(ndim,ndim);
      o2scl::matrix_transpose(ndim,ndim,chol,chol2);
      o2scl_cblas::dgemm(o2scl_cblas::o2cblas_RowMajor,
                         o2scl_cblas::o2cblas_NoTrans,
                         o2scl_cblas::o2cblas_NoTrans,
                         ndim,ndim,ndim,1.0,chol,chol2,0.0,m);
      return;
    }
    
    /** \brief Get the inverse of the covariance matrix
     */
    const mat_t &get_covar_inv() const {
      return covar_inv;
    }

    /** \brief Get the peak location
     */
    const vec_t &get_peak() const {
      return peak;
    }

    /** \brief Get the normalization
     */
    double get_norm() const {
      return norm;
    }

    /// The dimensionality
    virtual size_t dim() const {
      return ndim;
    }

    /// Create an empty distribution
    prob_dens_mdim_gaussian() {
      ndim=0;
      err_nonconv=true;
      verbose=0;
    }

    virtual ~prob_dens_mdim_gaussian() {
      ndim=0;
    }

    /// Copy constructor
    prob_dens_mdim_gaussian(const prob_dens_mdim_gaussian &pdmg_loc) {
      ndim=pdmg_loc.ndim;
      peak=pdmg_loc.peak;
      chol=pdmg_loc.chol;
      covar_inv=pdmg_loc.covar_inv;
      norm=pdmg_loc.norm;
      verbose=pdmg_loc.verbose;
    }

    /** \brief If true, call the error handler when convergence fails
     */
    bool err_nonconv;
    
    /** \brief Create a distribution from the covariance matrix
     */
    prob_dens_mdim_gaussian(size_t p_ndim, vec_t &p_peak, mat_t &covar) {
      set(p_ndim,p_peak,covar);
    }

    /// Copy constructor with operator=
    prob_dens_mdim_gaussian &operator=
    (const prob_dens_mdim_gaussian &pdmg_loc) {
      // Check for self-assignment
      if (this!=&pdmg_loc) {
        ndim=pdmg_loc.ndim;
        peak=pdmg_loc.peak;
        chol=pdmg_loc.chol;
        covar_inv=pdmg_loc.covar_inv;
        norm=pdmg_loc.norm;
        verbose=pdmg_loc.verbose;
      }
      return *this;
    }

    /// \name Set functions
    //@{
    /** \brief Create a distribution from a set of weighted samples from a 
        multidimensional Gaussian, returning the peak values and
        covariance matrix
      
        The matrix \c pts should have a size of \c n_pts in the first
        index and \c p_mdim in the second index
    */
    template<class mat2_t, class vec2_t,
             class mat2_col_t=const_matrix_column_gen<mat2_t> >
    int set_ret_wgts(size_t p_mdim, size_t n_pts, const mat2_t &pts,
                     const vec2_t &wgts, vec_t &peak_arg, mat_t &covar_arg) {
    
      // Set peak with average and diagonal elements in covariance
      // matrix with variance
      for(size_t i=0;i<p_mdim;i++) {
        const mat2_col_t col(pts,i);
        peak_arg[i]=o2scl::wvector_mean<mat2_col_t>(n_pts,col,wgts);
        // Square standard deviation
        covar_arg(i,i)=o2scl::wvector_stddev<mat2_col_t>(n_pts,col,wgts);
        covar_arg(i,i)*=covar_arg(i,i);
      }
      
      // Setup off-diagonal covariance matrix
      for(size_t i=0;i<p_mdim;i++) {
        mat2_col_t col_i(pts,i);
        for(size_t j=i+1;j<p_mdim;j++) {
          const mat2_col_t col_j(pts,j);
          double cov=o2scl::wvector_covariance(n_pts,col_i,col_j,wgts);
          covar_arg(i,j)=cov;
          covar_arg(j,i)=cov;
        }
      }
      set_covar(p_mdim,peak_arg,covar_arg);
      return 0;
    }

    /** \brief Create a distribution from a set of samples from a 
        multidimensional Gaussian, returning the peak values and
        covariance matrix
      
        The matrix \c pts should have a size of \c n_pts in the first
        index and \c p_mdim in the second index
    */
    template<class mat2_t, 
             class mat2_col_t=const_matrix_column_gen<mat2_t> >
    int set_ret(size_t p_mdim, size_t n_pts, const mat2_t &pts,
                vec_t &peak_arg, mat_t &covar_arg) {
      
      // Set peak with average and diagonal elements in covariance
      // matrix with variance
      for(size_t i=0;i<p_mdim;i++) {
        const mat2_col_t col(pts,i);
        peak_arg[i]=o2scl::vector_mean<mat2_col_t,double>(n_pts,col);
        // Square standard deviation
        covar_arg(i,i)=o2scl::vector_stddev<mat2_col_t>(n_pts,col);
        covar_arg(i,i)*=covar_arg(i,i);
      }
      
      // Setup off-diagonal covariance matrix
      for(size_t i=0;i<p_mdim;i++) {
        mat2_col_t col_i(pts,i);
        for(size_t j=i+1;j<p_mdim;j++) {
          const mat2_col_t col_j(pts,j);
          double cov=o2scl::vector_covariance(n_pts,col_i,col_j);
          covar_arg(i,j)=cov;
          covar_arg(j,i)=cov;
        }
      }

      if (this->verbose>0) {
        std::cout << "prob_dens_mdim_gaussian::set_ret():"
                  << std::endl;
        std::cout << "  peak: ";
        vector_out(std::cout,peak_arg,true);
        std::cout << "  covar: " << std::endl;
        matrix_out(std::cout,p_mdim,p_mdim,covar_arg,"  ");
      }
      
      set_covar(p_mdim,peak_arg,covar_arg);
      
      return 0;
    }
  
    /** \brief Create a distribution from a set of weighted samples from a 
        multidimensional Gaussian
      
        The matrix \c pts should have a size of \c n_pts in the first
        index and \c p_mdim in the second index
    */
    template<class mat2_t, class vec2_t,
             class mat2_col_t=const_matrix_column_gen<mat2_t> >
    int set_wgts(size_t p_mdim, size_t n_pts, const mat2_t &pts,
                 const vec2_t &wgts) {
      
      vec_t peak_arg(p_mdim);
      mat_t covar_arg(p_mdim,p_mdim);

      set_ret_wgts<mat2_t,vec2_t,mat2_col_t>(p_mdim,n_pts,pts,wgts,
                                             peak_arg,covar_arg);

      return 0;
    }
    
    /** \brief Create a distribution from a set of samples from a 
        multidimensional Gaussian
      
        The matrix \c pts should have a size of \c n_pts in the first
        index and \c p_mdim in the second index
    */
    template<class mat2_t, 
             class mat2_col_t=const_matrix_column_gen<mat2_t> >
    int set(size_t p_mdim, size_t n_pts, const mat2_t &pts) {
      
      vec_t peak_arg(p_mdim);
      mat_t covar_arg(p_mdim,p_mdim);

      set_ret<mat2_t,mat2_col_t>(p_mdim,n_pts,pts,peak_arg,
                                       covar_arg);
      
      return 0;
    }
    
    /** \brief Set the peak and covariance matrix for the distribution

        \note This function is called in constructors and thus 
        should not be virtual.
    */
    int set_covar(size_t p_ndim, vec_t &p_peak, mat_t &covar) {
      
      if (p_ndim==0) {
        O2SCL_ERR("Zero dimension in prob_dens_mdim_gaussian::set().",
                       o2scl::exc_einval);
      }
      ndim=p_ndim;
      norm=1.0;
      peak.resize(ndim);
      for(size_t i=0;i<ndim;i++) peak[i]=p_peak[i];

      if (this->verbose>0) {
        std::cout << "prob_dens_mdim_gaussian::set_covar():" << std::endl;
        std::cout << "  peak: ";
        vector_out(std::cout,ndim,peak,true);
        std::cout << "  covar: " << std::endl;
        matrix_out(std::cout,ndim,ndim,covar,"  ");
      }

      double sqrt_det;
      
      if (true) {
        
        // Perform the Cholesky decomposition of the covariance matrix
        chol=covar;
        o2scl_linalg::cholesky_decomp(ndim,chol);

        // The choleksy_decomp() function, by default, stores the
        // transpose in the upper triangular part, but the upper
        // triangular part needs to be zero for proper sampling, so we
        // set it to zero here.
        matrix_make_lower(ndim,ndim,chol);
        
        // Find the inverse
        covar_inv=chol;
        o2scl_linalg::cholesky_invert<mat_t>(ndim,covar_inv);

        // Compute the determinant of the Cholesky decomposision
        sqrt_det=1.0;
        for(size_t i=0;i<ndim;i++) {
          if (!std::isfinite(chol(i,i))) {
            O2SCL_CONV2_RET("An entry of the Cholesky decomposition was ",
                            "not finite in prob_dens_mdim_gaussian::set().",
                            o2scl::exc_einval,err_nonconv);
          }
          sqrt_det*=chol(i,i);
        }
        
      } else {

        o2scl_linalg::matrix_invert_det_cholesky<mat_t> mi;
        mi.invert_det(ndim,covar,covar_inv,sqrt_det);
        sqrt_det=sqrt(sqrt_det);
        
      }
      
      if (this->verbose>0) {
        std::cout << "  chol: " << std::endl;
        matrix_out(std::cout,ndim,ndim,chol,"  ");
        std::cout << "  covar_inv: " << std::endl;
        matrix_out(std::cout,ndim,ndim,covar_inv,"  ");
      }

      // Compute normalization
      norm=pow(2.0*o2scl_const::pi,-((double)ndim)/2.0)/sqrt_det;
      if (!std::isfinite(norm)) {
        O2SCL_CONV2_RET("Normalization not finite in ",
                        "prob_dens_mdim_gaussian::set().",
                        o2scl::exc_einval,err_nonconv);
      }
      
      return 0;
    }

    /** \brief Set the probability distribution from a 
        bivariate Gaussian
     */
    int set_from_biv(prob_dens_mdim_biv_gaussian<vec_t> &pdmbg) {

      double x_cent, y_cent, x_std, y_std, covar;
      pdmbg.get(x_cent,y_cent,x_std,y_std,covar);
      
      vec_t peak2(2);
      peak2[0]=x_cent;
      peak2[1]=y_cent;

      mat_t m_covar(2,2);
      m_covar(0,0)=x_std*x_std;
      m_covar(0,1)=x_std*y_std*m_covar;
      m_covar(1,0)=m_covar(0,1);
      m_covar(1,1)=y_std*y_std;

      set_covar(2,peak2,m_covar);
      
      return 0;
    }
    
    /** \brief Alternate set function for use when covariance matrix
        has already been decomposed and inverted
    */
    void set_alt(size_t p_ndim, vec_t &p_peak, mat_t &p_chol,
                 mat_t &p_covar_inv, double p_norm) {
      ndim=p_ndim;
      peak=p_peak;
      chol=p_chol;
      covar_inv=p_covar_inv;
      norm=p_norm;
      return;
    }

    /** \brief Given a data set and a covariance function, construct
        probability distribution based on a Gaussian process
        
        \note The type <tt>mat_col_t</tt> is a matrix column type for
        the internal object matrix type <tt>mat_t</tt>, and not
        associated with the data type <tt>vec_vec_t</tt>. Since the
        default matrix type is <tt>boost::numeric::ublas::matrix &lt;
        double &gt; </tt> a good matrix column type for this function
        is <tt>boost::numeric::ublas::matrix_column &lt;
        boost::numeric::ublas::matrix &lt; double &gt; &gt;</tt> .
        This matrix column type is needed for the LU decomposition and
        inversion.

        \future Clarify the relationship between this and
        interpm_krige.
    */
    template<class vec_vec_t, class mat_col_t, class func_t> 
    void set_gproc(size_t n_dim, size_t n_init, 
                   vec_vec_t &x, vec_t &y, func_t &fcovar) {
      
      // Construct the four covariance matrices
    
      mat_t KXsX(n_dim,n_init);
      for(size_t irow=n_init;irow<n_dim+n_init;irow++) {
        for(size_t icol=0;icol<n_init;icol++) {
          KXsX(irow-n_init,icol)=fcovar(x[irow],x[icol]);
        }
      }
    
      mat_t KXXs=boost::numeric::ublas::trans(KXsX);
    
      mat_t KXX(n_init,n_init);
      for(size_t irow=0;irow<n_init;irow++) {
        for(size_t icol=0;icol<n_init;icol++) {
          if (irow>icol) {
            KXX(irow,icol)=KXX(icol,irow);
          } else {
            KXX(irow,icol)=fcovar(x[irow],x[icol]);
          }
        }
      }
    
      mat_t KXsXs(n_dim,n_dim);
      for(size_t irow=n_init;irow<n_dim+n_init;irow++) {
        for(size_t icol=n_init;icol<n_dim+n_init;icol++) {
          if (irow>icol) {
            KXsXs(irow-n_init,icol-n_init)=KXsXs(icol-n_init,irow-n_init);
          } else {
            KXsXs(irow-n_init,icol-n_init)=fcovar(x[irow],x[icol]);
          }
        }
      }

      // Construct the inverse of KXX
      mat_t inv_KXX(n_init,n_init);
      vec_t mean(n_dim);
      mat_t covar(n_dim,n_dim);
      
      if (true) {
      
        o2scl::permutation p;
        int signum;
        o2scl_linalg::LU_decomp(n_init,KXX,p,signum);
        if (o2scl_linalg::diagonal_has_zero(n_dim,KXX)) {
          O2SCL_ERR2("KXX matrix is singular in ",
                     "prob_dens_mdim_gaussian::set_gproc().",
                     o2scl::exc_efailed);
        }
        o2scl_linalg::LU_invert<mat_t,mat_t,mat_col_t>(n_init,KXX,p,inv_KXX);
        
        // Compute the mean vector
        vec_t prod(n_init);
        boost::numeric::ublas::axpy_prod(inv_KXX,y,prod,true);
        boost::numeric::ublas::axpy_prod(KXsX,prod,mean,true);
        
        // Compute the covariance matrix
        mat_t prod2(n_init,n_dim), prod3(n_dim,n_dim);
        boost::numeric::ublas::axpy_prod(inv_KXX,KXXs,prod2,true);
        boost::numeric::ublas::axpy_prod(KXsX,prod2,prod3,true);
        covar=KXsXs-prod3;
        
      } else {
        
        o2scl_linalg::matrix_invert_det_cholesky<mat_t> mi;
        mi.invert(n_init,KXX,inv_KXX);

        // Compute the mean vector
        vec_t prod(n_init);
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           n_init,n_init,1.0,inv_KXX,y,0.0,prod);
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           n_init,n_init,1.0,KXsX,prod,0.0,mean);

        // Compute the covariance matrix
        mat_t prod2(n_init,n_dim), prod3(n_dim,n_dim);
        o2scl_cblas::dgemm(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           o2scl_cblas::o2cblas_NoTrans,
                           n_init,n_dim,n_dim,1.0,inv_KXX,KXXs,0.0,prod2);
        o2scl_cblas::dgemm(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           o2scl_cblas::o2cblas_NoTrans,
                           n_dim,n_init,n_dim,1.0,KXsX,prod2,0.0,prod3);
        covar=KXsXs-prod3;
        
      }
      
      // Now use set() in the parent class
      this->set(n_dim,mean,covar);
    
    }
    //@}

    /// \name Generic methods for multidimensional prob. dists.
    //@{
    /// The normalized density 
    virtual double pdf(const vec_t &x) const {
      if (ndim==0) {
        O2SCL_ERR2("Distribution not set in prob_dens_mdim_gaussian::",
                   "pdf().",o2scl::exc_einval);
      }
      double ret=norm;
      vec_t qq(ndim), vtmpx(ndim);
      for(size_t i=0;i<ndim;i++) {
        qq[i]=x[i]-peak[i];
      }

      o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                         o2scl_cblas::o2cblas_NoTrans,
                         ndim,ndim,1.0,covar_inv,qq,0.0,vtmpx);

      double ip=0.0;
      for(size_t i=0;i<ndim;i++) {
        ip+=qq[i]*vtmpx[i];
      }
      ret*=exp(-0.5*ip);
      
      return ret;
    }

    /// The log of the normalized density 
    virtual double log_pdf(const vec_t &x) const {
      if (ndim==0) {
        O2SCL_ERR2("Distribution not set in prob_dens_mdim_gaussian::",
                   "pdf().",o2scl::exc_einval);
      }
      double ret=log(norm);
      vec_t qq(ndim), vtmpx(ndim);
      for(size_t i=0;i<ndim;i++) qq[i]=x[i]-peak[i];
      
      o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                         o2scl_cblas::o2cblas_NoTrans,
                         ndim,ndim,1.0,covar_inv,qq,0.0,vtmpx);

      double ip=0.0;
      for(size_t i=0;i<ndim;i++) {
        ip+=qq[i]*vtmpx[i];
      }
      ret+=-0.5*ip;
      
      return ret;
    }

    /// Sample the distribution
    virtual void operator()(vec_t &x) const {
      if (ndim==0) {
        O2SCL_ERR2("Distribution not set in prob_dens_mdim_gaussian::",
                   "operator().",o2scl::exc_einval);
      }
      vec_t qq(ndim), vtmpx(ndim);
      for(size_t i=0;i<ndim;i++) qq[i]=pdg();

      o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                         o2scl_cblas::o2cblas_NoTrans,
                         ndim,ndim,1.0,chol,qq,0.0,vtmpx);

      for(size_t i=0;i<ndim;i++) x[i]=peak[i]+vtmpx[i];
      return;
    }
    //@}

    /** \brief Create a bivariate Gaussian probability distribution
     */
    prob_dens_mdim_biv_gaussian<vec_t> make_biv() const {

      if (ndim!=2) {
        O2SCL_ERR2("Distribution not two-dimensional in ",
                   "prob_dens_mdim_gaussian::make_biv().",
                   o2scl::exc_einval);
      }
      
      prob_dens_mdim_biv_gaussian<vec_t> pdmbg;
      
      double x_cent=peak[0], y_cent=peak[1];

      // The determinant of the inverse covariance matrix
      double det=covar_inv(0,0)*covar_inv(1,1)-
        covar_inv(0,1)*covar_inv(1,0);

      // Obtain the covariances by directly inverting the
      // inverse covariance matrix
      double x_std=covar_inv(1,1)/det, y_std=covar_inv(0,0)/det;
      double covar=-covar_inv(1,0)/det/x_std/y_std;

      // Set the bivariate Gaussian 
      pdmbg.set(x_cent,y_cent,x_std,y_std,covar);
      
      return pdmbg;
    }
    
  };

  /** \brief Gaussian distribution bounded by a hypercube

      \note This class naively resamples the Gaussian until
      a sample is within bounds. This is a temporary hack and
      can be very slow depending on the size of the volume
      excluded. 

      \warning The PDF is not yet properly normalized.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_t=boost::numeric::ublas::matrix<double> >
  class prob_dens_mdim_bound_gaussian :
    public prob_dens_mdim_gaussian<vec_t,mat_t> {
    
  protected:
  
    /** \brief Lower limits
     */
    vec_t low;
  
    /** \brief Upper limits
     */
    vec_t high;

  public:

    /** \brief Maximum number of samples (default \f$ 10^5 \f$)
     */
    size_t samp_max;
  
    /** \brief Create an empty distribution
     */
    prob_dens_mdim_bound_gaussian() {
      samp_max=100000;
    }
  
    /** \brief Create a distribution with the specified peak, covariance
        matrix, lower limits, and upper limits
    */
    prob_dens_mdim_bound_gaussian(size_t p_ndim, vec_t &p_peak, mat_t &covar,
                                  vec_t &p_low, vec_t &p_high) {
      set(p_ndim,p_peak,covar,p_low,p_high);
      samp_max=100000;
    }
  
    /** \brief Set the peak, covariance matrix, lower limits, and upper
        limits

        \note This function is called in constructors and thus 
        should not be virtual.
    */
    void set(size_t p_ndim, vec_t &p_peak, mat_t &covar,
             vec_t &p_low, vec_t &p_high) {
      prob_dens_mdim_gaussian<vec_t,mat_t>::set(p_ndim,p_peak,covar);
      low=p_low;
      high=p_high;
      return;
    }
  
    /** \brief Compute the probability density function (arbitrary
        normalization)
    */
    virtual double pdf(const vec_t &x) const {
      for(size_t i=0;i<this->ndim;i++) {
        if (x[i]<low[i]) {
          O2SCL_ERR("Parameter too small in pdf().",
                    o2scl::exc_einval);
        }
        if (x[i]>high[i]) {
          O2SCL_ERR("Parameter too large in pdf().",
                    o2scl::exc_einval);
        }
      }
      return prob_dens_mdim_gaussian<vec_t,mat_t>::pdf(x);
    }
  
    /** \brief Compute the natural log of the probability density function
        (arbitrary normalization)
    */
    virtual double log_pdf(const vec_t &x) const {
      for(size_t i=0;i<this->ndim;i++) {
        if (x[i]<low[i]) {
          O2SCL_ERR("Parameter too small in log_pdf().",
                    o2scl::exc_einval);
        }
        if (x[i]>high[i]) {
          O2SCL_ERR("Parameter too large in log_pdf().",
                    o2scl::exc_einval);
        }
      }
      return prob_dens_mdim_gaussian<vec_t,mat_t>::log_pdf(x);
    }
  
    /** \brief Sample the distribution
     */
    virtual void operator()(vec_t &x) const {
      bool done=false;
      size_t j=0;
      while (done==false) {
        done=true;
        prob_dens_mdim_gaussian<vec_t,mat_t>::operator()(x);
        j++;
        for(size_t i=0;i<this->ndim;i++) {
          if (x[i]<low[i]) {
            done=false;
            i=this->ndim;
          } else if (x[i]>high[i]) {
            done=false;
            i=this->ndim;
          }
        }
        if (j>samp_max) {
          O2SCL_ERR2("Sampling failed in ",
                     "prob_dens_mdim_bound_gaussian::operator().",
                     o2scl::exc_einval);
        }
      }
      return;
    }
  
  };
  
  /** \brief A multi-dimensional conditional probability density function

      Note that conditional probabilities are typically written \f$
      P(A|B) \f$, i.e. the probability of \f$ A \f$ given \f$ B \f$.
      \o2 arranges the function parameters for the functions \ref
      o2scl::prob_cond_mdim::pdf, \ref o2scl::prob_cond_mdim::log_pdf
      \ref o2scl::prob_cond_mdim::operator()(), so that \f$ B \f$ is
      given first, and \f$ A \f$ is second.

      \ref
      o2scl::prob_cond_mdim::log_metrop_hast is a vector from \f$ B
      \f$ as denoted above.
      
      This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class prob_cond_mdim {
    
  public:

    virtual ~prob_cond_mdim() {
    }
    
    /// The dimensionality
    virtual size_t dim() const {
      return 0;
    }
  
    /** \brief The conditional probability of x_A given x_B, 
        i.e. \f$ P(A|B) \f$
    */
    virtual double pdf(const vec_t &x_B, const vec_t &x_A) const=0;
  
    /** \brief The log of the conditional probability of x_A given x_B
        i.e. \f$ \log [P(A|B)] \f$
    */
    virtual double log_pdf(const vec_t &x_B, const vec_t &x_A) const=0;
  
    /// Sample the distribution
    virtual void operator()(const vec_t &x_B, vec_t &x_A) const=0;

    /** \brief Sample the distribution and return the 
        log of the Metropolis-Hastings ratio

        The Metropolis-Hastings ratio for a step beginning at \f$ x \f$
        and ending at \f$ x^{\prime} \f$ is 
        obeys
        \f[
        \frac{P(x^{\prime})g(x|x^{\prime})}{P(x)g(x^{\prime}|x)}
        \f]
        taking the log, this gives 
        \f[
        \log[P(x^{\prime})] - \log[P(x)] + 
        \log \left[ \frac{g(x|x^{\prime})}{g(x^{\prime}|x)} \right]
        \f]
        thus this function computes 
        \f[
        \log \left[ g(x|x^{\prime}) \right]
        - \log \left[ g(x^{\prime}|x) \right]
        \f]
        and thus, to keep a similar notation to 
        \ref prob_cond_mdim::pdf() where \f$ g(x^{\prime}|x) \f$
        is obtained from 
        \code
        pdf(x,x_prime)
        \endcode
        this function computes
        \code
        h(x,x_prime) = log_pdf(x_prime,x)-log_pdf(x,x_prime);
        \endcode

        To check this, in the limit that \f$ g(x|x^{\prime}) 
        \rightarrow P(x) \f$ this function returns
        \f[
        \log \left[ \frac{P(x)}{P(x^{\prime})} \right]
        \f]

    */
    virtual double log_metrop_hast(const vec_t &x, vec_t &x_prime) const {
      operator()(x,x_prime);
      double val1=log_pdf(x_prime,x);
      double val2=log_pdf(x,x_prime);
      if (!std::isfinite(val1)) {
        std::cout << "val1: " << val1 << std::endl;
        O2SCL_ERR("Log pdf not finite in log_metrop_hast 1.",
                  o2scl::exc_efailed);
      }
      if (!std::isfinite(val2)) {
        std::cout << "val2: " << val2 << std::endl;
        O2SCL_ERR("Log pdf not finite in log_metrop_hast 2.",
                  o2scl::exc_efailed);
      }
      return val1-val2;
    }
  
  };

  /** \brief Conditional probability for a random walk inside a
      hypercube

      \comment
      I had previously used std::uniform_real_distribution
      instead of rng, but this caused problems with
      intel compilers.
      \endcomment

      This conditional probability is most useful in providing a
      Metropolis-Hastings distribution with a fixed step size which
      properly handles a boundary. The Metropolis-Hastings step is
      accepted if \f$ r \in [0,1] \f$ obeys
      \f[
      r < \frac{P(x^{\prime})g(x|x^{\prime})}
      {P(x)g(x^{\prime}|x)}
      \f]
      The function \f$ g(x^{\prime}|x) = g(x^{\prime},x)/g(x) \f$, and
      because of the fixed step size these probabilities are just
      proportional to the inverse allowed volume, i.e. \f$ V(x)^{-1}
      V^{-1}(x^{\prime}) / V(x)^{-1} = V^{-1}(x^{\prime}) \f$ . If \f$
      x^{\prime} \f$ is near a boundary then \f$ V(x^{\prime}) \f$ is
      decreased and the conditional probability increases accordingly.
      If the distance between \f$ x \f$ and \f$ x^{\prime} \f$ is
      unreachable in a step, then the PDF is zero.

      \note This class is experimental.

      \note Const functions are not thread-safe because the
      class contains an internal mutable random number generator
      object.

      \comment
      If we do not include the g ratio, then the edges
      will be undersampled because we don't properly include
      the rejections 

      For \f$ g(x^{\prime}|x) \f$, if x is near the edge, then the
      cond. prob. is larger, thus the g ratio is smaller than 1,
      encouraging more rejections.
      \endcomment
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class prob_cond_mdim_fixed_step : public prob_cond_mdim<vec_t> {

  protected:

    /** \brief Step sizes
     */
    std::vector<double> u_step;

    /** \brief Lower limits
     */
    std::vector<double> u_low;

    /** \brief Upper limits
     */
    std::vector<double> u_high;

    /** \brief Internal random number generator
     */
    rng<> rg;
  
    /** \brief Internal set function

        \comment
        This can't be virtual because it needs to be called
        by the constructor
        \endcomment
    */
    int set_internal(size_t sz, vec_t &step, vec_t &low, vec_t &high) {
      u_step.resize(sz);
      u_low.resize(sz);
      u_high.resize(sz);
      for(size_t i=0;i<sz;i++) {
        u_step[i]=step[i];

        if (!std::isfinite(low[i]) || !std::isfinite(high[i])) {
          O2SCL_ERR2("Limit not finite in prob_cond_mdim_fixed_step::",
                     "set_internal().",o2scl::exc_einval);
        }
      
        // Force low and high to be properly ordered
        if (low[i]>high[i]) {
          u_high[i]=low[i];
          u_low[i]=high[i];
        } else {
          u_low[i]=low[i];
          u_high[i]=high[i];
        }
      }
      return 0;
    }

  public:
  
    prob_cond_mdim_fixed_step() {
    }

    /// Copy constructor
    prob_cond_mdim_fixed_step(const prob_cond_mdim_fixed_step &pcmfs) {
      u_step=pcmfs.u_step;
      u_low=pcmfs.u_low;
      u_high=pcmfs.u_high;
    }
  
    /// Copy constructor with operator=
    prob_cond_mdim_fixed_step &operator=
    (const prob_cond_mdim_fixed_step &pcmfs) {

      // Check for self-assignment
      if (this!=&pcmfs) {
        u_step=pcmfs.u_step;
        u_low=pcmfs.u_low;
        u_high=pcmfs.u_high;
      }
      return *this;
    }

    /// Virtual destructor
    virtual ~prob_cond_mdim_fixed_step() {
    }
  
    /** \brief Set the random number generator seed
     */
    void set_seed(unsigned long int s) {
      rg.set_seed(s);
      return;
    }

    /** \brief Create a conditional probability object
        with specified step sizes and limits
    */
    template<class=vec_t> prob_cond_mdim_fixed_step
    (vec_t &step, vec_t &low, vec_t &high) {
      if (step.size()!=low.size()) {
        O2SCL_ERR2("Vectors 'step' and 'low' mismatched in ",
                   "prob_cond_mdim_fixed_step constructor.",
                   o2scl::exc_einval);
      }
      if (step.size()!=high.size()) {
        O2SCL_ERR2("Vectors 'step' and 'high' mismatched in ",
                   "prob_cond_mdim_fixed_step constructor.",
                   o2scl::exc_einval);
      }
      set_internal(step.size(),step,low,high);
    }
  
    /** \brief Set step sizes and limits
     */
    virtual int set(vec_t &step, vec_t &low, vec_t &high) {
      if (step.size()!=low.size()) {
        O2SCL_ERR2("Vectors 'step' and 'low' mismatched in ",
                   "prob_cond_mdim_fixed_step::set().",
                   o2scl::exc_einval);
      }
      if (step.size()!=high.size()) {
        O2SCL_ERR2("Vectors 'step' and 'high' mismatched in ",
                   "prob_cond_mdim_fixed_step::set().",
                   o2scl::exc_einval);
      }
      set_internal(step.size(),step,low,high);
      return 0;
    }

    /// The dimensionality
    virtual size_t dim() const {
      return u_step.size();
    }
  
    /** \brief The conditional probability of x_A given x_B, 
        i.e. \f$ P(A|B) \f$
    */
    virtual double pdf(const vec_t &x_B, const vec_t &x_A) const {
      double vol1=1.0;
      for(size_t i=0;i<u_step.size();i++) {

        if (fabs(x_A[i]-x_B[i])>u_step[i]) return 0.0;
      
        if (x_B[i]-u_step[i]<u_low[i]) {
          if (x_B[i]+u_step[i]>u_high[i]) {
            // If x+step is too large and x-step is too small
            vol1*=u_high[i]-u_low[i];
          } else {
            // If x-step is too small
            vol1*=x_B[i]+u_step[i]-u_low[i];
          }
        } else {
          if (x_B[i]+u_step[i]>u_high[i]) {
            // If x+step is too large
            vol1*=u_high[i]-(x_B[i]-u_step[i]);
          } else {
            // The normal case, where the volumes are both inside
            // of the boundaries
            vol1*=2.0*u_step[i];
          }
        }
      }
      return 1.0/vol1;
    }
  
    /** \brief The log of the conditional probability of x_A given x_B
        i.e. \f$ \log [P(A|B)] \f$
    */
    virtual double log_pdf(const vec_t &x_B, const vec_t &x_A) const {
      return log(pdf(x_B,x_A));
    }
  
    /// Sample the distribution
    virtual void operator()(const vec_t &x_B, vec_t &x_A) const {
      size_t nv=u_step.size();
      for(size_t i=0;i<nv;i++) {
        if (x_B[i]<u_low[i] || x_B[i]>u_high[i]) {
          std::cout << "Input out of bounds in fixed_step::operator(): "
                    << i << " " << x_B[i] << " "
                    << u_low[i] << " " << u_high[i] << std::endl;
          O2SCL_ERR("Input out of bounds in fixed_step::operator().",
                    o2scl::exc_einval);
        }
        do {
          x_A[i]=x_B[i]+u_step[i]*(rg.random()*2.0-1.0);
        } while (x_A[i]<u_low[i] || x_A[i]>u_high[i]);
      }
      return;
    }
  
  };

  /** \brief A multi-dimensional conditional probability density function
      independent of the input

      The conditional probability, \f$ P(A|B) = P(A,B)/P(B) \f$. If
      the joint probability is factorizable because the events \f$ A
      \f$ and \f$ B \f$ are independent, i.e. \f$ P(A,B) = P(A) P(B)
      \f$, then \f$ P(A|B) = P(A) \f$ and is independent of \f$ B \f$.
      This class handles that particular case.

      \note This class stores a shared pointer of the underlying
      probability distribution, so copies created by the copy 
      constructor point to the same object. If this object is not
      thread-safe, then copies of this class are also not
      thread-safe.
      
      This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class prob_cond_mdim_indep : public prob_cond_mdim<vec_t> {

  protected:

    /// The underlying probability distribution
    std::shared_ptr<prob_dens_mdim<vec_t> > base;
    
  public:

    /** \brief Create a conditional probability distribution
        based on a default multidimensional Gaussian
    */
    prob_cond_mdim_indep() : 
      base(new prob_dens_mdim_gaussian<vec_t>) {
    }
    
    /** \brief Create a conditional probability distribution
        based on the specified probability distribution
    */
    prob_cond_mdim_indep(std::shared_ptr<prob_dens_mdim<vec_t>> out) : 
      base(out) {
    }
    
    /// Copy constructor
    prob_cond_mdim_indep(const prob_cond_mdim_indep &pcmi) {
      base=pcmi.base;
    }
  
    /// Copy constructor with operator=
    prob_cond_mdim_indep &operator=(const prob_cond_mdim_indep &pcmi) {
      // Check for self-assignment
      if (this!=&pcmi) {
        base=pcmi.base;
      }
      return *this;
    }

    /// Set the base probability distribution
    void set_base(std::shared_ptr<prob_dens_mdim<vec_t>> b) {
      base=b;
      return;
    }
  
    /// The dimensionality
    virtual size_t dim() const {
      return base->dim();
    }
  
    /** \brief The conditional probability of x_A given x_B, 
        i.e. \f$ P(A|B) \f$
    */
    virtual double pdf(const vec_t &x_B, const vec_t &x_A) const {
      return base->pdf(x_A);
    }
  
    /** \brief The log of the conditional probability of x_A given x_B
        i.e. \f$ \log [P(A|B)] \f$
    */
    virtual double log_pdf(const vec_t &x_B, const vec_t &x_A) const {
      double r=base->log_pdf(x_A);
      return r;
    }
  
    /// Sample the distribution
    virtual void operator()(const vec_t &x_B, vec_t &x_A) const {
      (*base)(x_A);
      return;
    }

    virtual const char *type() { return "prob_cons_mdim_indep"; }
  
  };
  
  /** \brief A multi-dimensional Gaussian conditional probability
      density function

      This class is experimental.

      \note Const functions are not thread-safe because
      mutable storage is used.

      \verbatim embed:rst
      
      .. todo::

         In class prob_cond_mdim_gaussian:

         - This should be a symmetric conditional probability, 
           i.e. \f$ P(x|y) = P(y|x) \f$. Test this.

      \endverbatim
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_t=boost::numeric::ublas::matrix<double> >
  class prob_cond_mdim_gaussian : public prob_cond_mdim<vec_t> {
    
  protected:

    /// Cholesky decomposition
    mat_t chol;

    /// Inverse of the covariance matrix
    mat_t covar_inv;

    /// Normalization factor
    double norm;

    /// Number of dimensions
    size_t ndim;

    /// Temporary storage 1
    mutable vec_t q;

    /// Temporary storage 2
    mutable vec_t vtmp;

    /// Standard normal
    mutable o2scl::prob_dens_gaussian pdg;
    
  public:

    /** \brief Create an empty distribution 
     */
    prob_cond_mdim_gaussian() {
      ndim=0;
    }
  
    /// Copy constructor
    prob_cond_mdim_gaussian(const prob_cond_mdim_gaussian &pcmg) {
      ndim=pcmg.ndim;
      chol=pcmg.chol;
      covar_inv=pcmg.covar_inv;
      norm=pcmg.norm;
      q.resize(ndim);
      vtmp.resize(ndim);
    }
  
    /// Copy constructor with operator=
    prob_cond_mdim_gaussian &operator=(const prob_cond_mdim_gaussian &pcmg) {
      // Check for self-assignment
      if (this!=&pcmg) {
        ndim=pcmg.ndim;
        chol=pcmg.chol;
        covar_inv=pcmg.covar_inv;
        norm=pcmg.norm;
        q.resize(ndim);
        vtmp.resize(ndim);
      }
      return *this;
    }
  
    /** \brief Create a distribution from the covariance matrix
     */
    prob_cond_mdim_gaussian(size_t p_ndim, mat_t &covar) {
      set(p_ndim,covar);
    }
  
    /// The dimensionality
    virtual size_t dim() const {
      return ndim;
    }

    /// Set the seed
    void set_seed(unsigned long int s) {
      pdg.set_seed(s);
      return;
    }

    /** \brief Clear the distribution
     */
    void clear() {
      covar_inv.clear();
      chol.clear();
      ndim=0;
      return;
    }
    
    /** \brief Set the covariance matrix for the distribution
     */
    void set(size_t p_ndim, mat_t &covar) {
      if (p_ndim==0) {
        O2SCL_ERR("Zero dimension in prob_cond_mdim_gaussian::set().",
                  o2scl::exc_einval);
      }
      ndim=p_ndim;
      norm=1.0;
      q.resize(ndim);
      vtmp.resize(ndim);

      // Perform the Cholesky decomposition
      chol=covar;
      o2scl_linalg::cholesky_decomp(ndim,chol);
      
      // The choleksy_decomp() function, by default, stores the
      // transpose in the upper triangular part, but the
      // upper triangular part needs to be zero for proper
      // sampling, so we set it to zero here.
      matrix_make_lower(ndim,ndim,chol);
        
      // Find the inverse
      covar_inv=chol;
      o2scl_linalg::cholesky_invert<mat_t>(ndim,covar_inv);
      
      // Compute the determinant of the Cholesky decomposision
      double sqrt_det=1.0;
      for(size_t i=0;i<ndim;i++) {
        if (!std::isfinite(chol(i,i))) {
          O2SCL_ERR2("An entry of the Cholesky decomposition was not finite ",
                     "in prob_cond_mdim_gaussian::set().",o2scl::exc_einval);
        }
        sqrt_det*=chol(i,i);
      }
    
      // Compute normalization
      norm=pow(2.0*o2scl_const::pi,-((double)ndim)/2.0)/sqrt_det;
      if (!std::isfinite(norm)) {
        O2SCL_ERR2("Normalization not finite in ",
                   "prob_dens_mdim_gaussian::set().",o2scl::exc_einval);
      }
    }

    /** \brief The conditional probability of x_A given x_B, 
        i.e. \f$ P(A|B) \f$
    */
    virtual double pdf(const vec_t &x_B, const vec_t &x_A) const {
      if (ndim==0) {
        O2SCL_ERR2("Distribution not set in prob_cond_mdim_gaussian::",
                   "pdf().",o2scl::exc_einval);
      }
      double ret=norm;
      for(size_t i=0;i<ndim;i++) q[i]=x_A[i]-x_B[i];
      vtmp=prod(covar_inv,q);
      ret*=exp(-0.5*inner_prod(q,vtmp));
      return ret;
    }

    /** \brief The log of the conditional probability of x_A given x_B
        i.e. \f$ \log [P(A|B)] \f$
    */
    virtual double log_pdf(const vec_t &x_B, const vec_t &x_A) const {
      if (ndim==0) {
        O2SCL_ERR2("Distribution not set in prob_cond_mdim_gaussian::",
                   "pdf().",o2scl::exc_einval);
      }
      double ret=log(norm);
      for(size_t i=0;i<ndim;i++) q[i]=x_A[i]-x_B[i];
      vtmp=prod(covar_inv,q);
      ret+=-0.5*inner_prod(q,vtmp);
      /*
        std::cout << "pdmg lp: ";
        for (size_t i=0;i<ndim;i++) std::cout << x_A[i] << " ";
        for (size_t i=0;i<ndim;i++) std::cout << q[i] << " ";
        for (size_t i=0;i<ndim;i++) std::cout << vtmp[i] << " ";
        std::cout << ret << std::endl;
      */
      return ret;
    }

    /// Sample the distribution
    virtual void operator()(const vec_t &x_B, vec_t &x_A) const {
      if (ndim==0) {
        O2SCL_ERR2("Distribution not set in prob_cond_mdim_gaussian::",
                   "operator().",o2scl::exc_einval);
      }
      for (size_t i=0;i<ndim;i++) q[i]=pdg();
      vtmp=prod(chol,q);
      for (size_t i=0;i<ndim;i++) x_A[i]=x_B[i]+vtmp[i];
      /*
        std::cout << "pdmg op: ";
        for (size_t i=0;i<ndim;i++) std::cout << x_A[i] << " ";
        std::cout << std::endl;
      */
      return;
    }

  };

  /** \brief A probability density distribution from a Gaussian 
      mixture model
   */
  template<class gauss_vec_t=boost::numeric::ublas::vector<double>,
           class gauss_mat_t=boost::numeric::ublas::matrix<double>>
  class prob_dens_mdim_gmm : public prob_dens_mdim<gauss_vec_t> {
    
  public:

    /** \brief Base random number generator

        This is used to automatically generate initial means for the
        Gaussians if a user-specified guess is not provided.
     */
    rng<> r2;
    
    /// The Gaussians
    std::vector<o2scl::prob_dens_mdim_gaussian
                <gauss_vec_t,gauss_mat_t>> pdmg;
    
    typedef boost::numeric::ublas::vector<double> internal_vec_t;
    
    /// The weights (must add to 1) 
    internal_vec_t weights;

    /// Verbosity parameter
    int verbose;
    
    prob_dens_mdim_gmm() {
      verbose=0;
    }
    
    /// Copy constructor
    prob_dens_mdim_gmm(const prob_dens_mdim_gmm &pdmg_loc) {
      pdmg=pdmg_loc.pdmg;
      weights=pdmg_loc.weights;
      verbose=pdmg_loc.verbose;
    }

    /// Copy constructor with operator=
    prob_dens_mdim_gmm &operator=
    (const prob_dens_mdim_gmm &pdmg_loc) {
      // Check for self-assignment
      if (this!=&pdmg_loc) {
        pdmg=pdmg_loc.pdmg;
        weights=pdmg_loc.weights;
        verbose=pdmg_loc.verbose;
      }
      return *this;
    }

    /** \brief Clear the distribution
     */
    void clear() {
      weights.clear();
      pdmg.clear();
      return;
    }

    /** \brief Return the number of Gaussians, or 0 if
        the object is empty
    */
    size_t get_n_components() {
      return pdmg.size();
    }
    
    /** \brief Read the Gaussian mixture from an input file
     */
    virtual int read_generic(std::istream &fin) {
      
      double data;
      std::string line;
      std::string cname;
      
      size_t nd_in;
      // Read first line and into list
      fin >> nd_in;
      weights.resize(nd_in);
      for(size_t i=0;i<nd_in;i++) {
        fin >> weights[i];
      }
      pdmg.resize(nd_in);
      for(size_t i=0;i<nd_in;i++) {
        pdmg[i].read_generic(fin);
      }
      return 0;
    }

    /** \brief Write the Gaussian mixture to an output file
     */
    virtual int write_generic(std::ostream &fout) {

      double data;
      std::string line;
      std::string cname;

      fout << weights.size() << std::endl;
      o2scl::vector_out(fout,weights,true);
      for(size_t i=0;i<weights.size();i++) {
        pdmg[i].write_generic(fout);
      }
      return 0;
    }
    
    /** \brief Sample the distribution
     */
    virtual void operator()(gauss_vec_t &x) const {
      if (weights.size()==0) {
        O2SCL_ERR2("No Gaussians specified in ",
                  "prob_dens_mdim_gmm::operator().",o2scl::exc_einval);
      }
      internal_vec_t partial_sums(weights.size());
      for(size_t i=0;i<weights.size();i++) {
        if (i==0) {
          partial_sums[0]=weights[0];
        } else {
          partial_sums[i]=partial_sums[i-1]+weights[i];
        }
      }

      double v=r2.random();
      for(size_t k=0;k<weights.size();k++) {
        if (v<partial_sums[k] || k==weights.size()-1) {
          pdmg[k](x);
          return;
        }
      }
      O2SCL_ERR2("Weight arithmetic problem in ",
                 "prob_dens_mdim_gmm::operator().",
                 o2scl::exc_esanity);
      return;
    }

    /// Return the dimensionality
    virtual size_t dim() const {
      if (weights.size()==0) return 0;
      return pdmg[0].dim();
    }

    /** \brief Compute the normalized probability density
     */
    virtual double pdf(const gauss_vec_t &x) const {
      double ret=0.0;
      for(size_t k=0;k<weights.size();k++) {
        ret+=weights[k]*pdmg[k].pdf(x);
      }
      return ret;
    }
    
    /** \brief The log of the normalized density
     */
    virtual double log_pdf(const gauss_vec_t &x) const {
      if (weights.size()==0) {
        O2SCL_ERR2("No Gaussians specified in ",
                  "prob_dens_mdim_gmm::operator().",o2scl::exc_einval);
      }
      return log(pdf(x));
    }
    
  };

  /** \brief A vector of multidimensional probability distributions
   */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class vec_prob_dens_mdim {

  protected:

    /// The internal vector of pointers
    std::vector<prob_dens_mdim<vec_t> *> list;
    
  public:

    virtual ~vec_prob_dens_mdim() {
      free();
    }

    /// Return the vector size
    size_t size() {
      return list.size();
    }
    
    /// Clear all of the memory
    void free() {
      for(size_t i=0;i<list.size();i++) {
        delete list[i];
      }
      list.clear();
    }
    
    /// Return a const reference
    virtual const prob_dens_mdim<vec_t> &operator()(size_t ix) const {
      return (*(list[ix]));
    }
    
    /// Return a non-const reference
    virtual prob_dens_mdim<vec_t> &operator()(size_t ix) {
      return (*(list[ix]));
    }
    
    /// Add a distribution of a template type and return a reference
    template<class prob_mdim_t> prob_mdim_t &add() {
      prob_mdim_t *ptr=new prob_mdim_t;
      list.push_back(ptr);
      return *ptr;
    }
    
  };
  
  /** \brief A vector of conditional probability distributions
   */
  template<class vec_t=boost::numeric::ublas::vector<double> >
  class vec_prob_cond_mdim {

  protected:

    /// The internal vector of pointers
    std::vector<prob_cond_mdim<vec_t> *> list;
    
  public:

    virtual ~vec_prob_cond_mdim() {
      free();
    }

    /// Return the vector size
    size_t size() {
      return list.size();
    }
    
    /// Clear all of the memory
    void free() {
      for(size_t i=0;i<list.size();i++) {
        delete list[i];
      }
      list.clear();
    }
    
    /// Return a const reference
    virtual const prob_cond_mdim<vec_t> &operator()(size_t ix) const {
      return (*(list[ix]));
    }
    
    /// Return a non-const reference
    virtual prob_cond_mdim<vec_t> &operator()(size_t ix) {
      return (*(list[ix]));
    }
    
    /// Add a distribution of a template type and return a reference
    template<class cond_mdim_t> cond_mdim_t &add() {
      cond_mdim_t *ptr=new cond_mdim_t;
      list.push_back(ptr);
      return *ptr;
    }
    
    /// Add a prob_cond_mdim_indep distribution
    prob_cond_mdim_indep<vec_t> &add_cond_mdim_indep
    (prob_dens_mdim<vec_t> &base) {
      prob_cond_mdim_indep<vec_t> *ptr=new
        prob_cond_mdim_indep<vec_t>(base);
      list.push_back(ptr);
    }
    
  };
  
#ifdef O2SCL_NEVER_DEFINED
  
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_t=const_matrix_view_table<> >
  class prob_dens_mdim_kde : public prob_dens_mdim<vec_t> {

  protected:

    size_t n_dim;
    size_t n_points;
    mat_t *dp;
    prob_dens_gaussian pdg;
    double bw;
    
  public:

    prob_dens_mdim_kde() {
      n_dim=0;
      n_points=0;
    }

    int set_data(mat_t &data) {
      n_dim=data.size1();
      n_points=data.size2();
      dp=&data;
      return 0;
    }

    /// Return the dimensionality
    virtual size_t dim() const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0;
    }
  
    /// The normalized density 
    virtual double pdf(const vec_t &x) const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return 0.0;
    }
  
    /// The log of the normalized density 
    virtual double log_pdf(const vec_t &x) const {
      double val=pdf(x);
      if (!std::isfinite(val) || val<0.0) {
        O2SCL_ERR2("PDF not finite or negative in ",
                   "prob_dens_mdim::log_pdf().",o2scl::exc_efailed);
      }
      double val2=log(pdf(x));
      if (!std::isfinite(val2)) {
        std::cout << val << " " << val2 << std::endl;
        O2SCL_ERR2("Log of PDF not finite in ",
                   "prob_dens_mdim::log_pdf().",o2scl::exc_efailed);
      }
      return val2;
    }
  
    /// Sample the distribution
    virtual void operator()(vec_t &x) const {
      O2SCL_ERR("Executing blank parent function.",o2scl::exc_eunimpl);
      return;
    }

    
  };
#endif
  
}

#endif
