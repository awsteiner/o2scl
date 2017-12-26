/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2018, Andrew W. Steiner
  
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
/** \file prob_dens_func.h
    \brief File for probability density functions
*/
#ifndef O2SCL_PROB_DENS_FUNC_H
#define O2SCL_PROB_DENS_FUNC_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <o2scl/hist.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/search_vec.h>
#include <o2scl/cholesky.h>
#include <o2scl/lu.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A one-dimensional probability density function

      This class is experimental.

      \future Give functions for mean, median, mode, variance, etc?

      \comment
      For now, there aren't any pure virtual functions,
      since this causes problems in creating an
      std::vector<prob_dens_func> object below (especially
      with intel compilers)
      \endcomment
  */
  class prob_dens_func {
    
  public:
    
    /// Sample from the specified density
    virtual double operator()() const {
      return 0.0;
    }
    
    /// The normalized density 
    virtual double pdf(double x) const {
      return 0.0;
    }
    
    /// The log of the normalized density 
    virtual double log_pdf(double x) const {
      return 0.0;
    }
    
    /// The cumulative distribution function (from the lower tail)
    virtual double cdf(double x) const {
      return 0.0;
    }
    
    /// The inverse cumulative distribution function
    virtual double invert_cdf(double cdf) const {
      return 0.0;
    }

    /// Entropy of the distribution (\f$ - \int f \ln f \f$ )
    virtual double entropy() const {
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
    
    /// Base GSL random number generator
    o2scl::rng_gsl r;
    
  public:
    
    /** \brief Create a standard normal distribution
     */
    prob_dens_gaussian() {
      cent_=0.0;
      sigma_=1.0;
    }

    /** \brief Create a Gaussian distribution with width \c sigma

	The value of \c sigma must be larger than zero.
    */
    prob_dens_gaussian(double cent, double sigma) {
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
      r=pdg.r;
    }

    /// Copy constructor with operator=
    prob_dens_gaussian &operator=(const prob_dens_gaussian &pdg) {
      // Check for self-assignment
      if (this!=&pdg) {
	cent_=pdg.cent_;
	sigma_=pdg.sigma_;
	r=pdg.r;
      }
      return *this;
    }

    /// Set the seed
    void set_seed(unsigned long int s) {
      r.set_seed(s);
      return;
    }

    /// Set the center
    void set_center(double cent) {
      cent_=cent;
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
      return cent_+gsl_ran_gaussian(&r,sigma_);
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
	O2SCL_ERR2("Requested cdf inverse outside of [0,1] in ",
		   "prob_dens_gaussian::invert_cdf().",exc_einval);
      }
      return gsl_cdf_gaussian_Pinv(in_cdf,sigma_)+cent_;
    }

    /// The inverse cumulative distribution function
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
    rng_gsl r;
    
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
      return gsl_ran_flat(&r,ll,ul);
    }
    
    /// The normalized density 
    virtual double pdf(double x) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "pdf().",exc_einval);
      }
      if (x<ll || x>ul) return 0.0;
      return gsl_ran_flat_pdf(x,ll,ul);
    }
    
    /// The log of the normalized density 
    virtual double log_pdf(double x) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "pdf().",exc_einval);
      }
      if (x<ll || x>ul) return 0.0;
      return log(gsl_ran_flat_pdf(x,ll,ul));
    }
    
    /// The cumulative distribution function (from the lower tail)
    virtual double cdf(double x) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "cdf().",exc_einval);
      }
      if (x<ll) return 0.0;
      if (x>ul) return 1.0;
      return gsl_cdf_flat_P(x,ll,ul);
    }
    
    /// The inverse cumulative distribution function
    virtual double invert_cdf(double in_cdf) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "invert_cdf().",exc_einval);
      }
      if (in_cdf<0.0 || in_cdf>1.0) {
	O2SCL_ERR2("Requested cdf inverse outside of [0,1] in ",
		   "prob_dens_uniform::invert_cdf().",exc_einval);
      }
      return gsl_cdf_flat_Pinv(in_cdf,ll,ul);
    }

    /// The inverse cumulative distribution function
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
    rng_gsl r;
    
  public:

    /** \brief Create a blank lognormal distribution
     */
    prob_dens_lognormal() {
      sigma_=-1.0;
      mu_=0.0;
    }
    
    /** \brief Create lognormal distribution with mean parameter \c mu
	and width parameter \c sigma

	The value of \c sigma must be larger than zero.
    */
    prob_dens_lognormal(double mu, double sigma) {
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
    }

    /// Copy constructor with operator=
    prob_dens_lognormal &operator=(const prob_dens_lognormal &pdg) {
      // Check for self-assignment
      if (this==&pdg) return *this;
      mu_=pdg.mu_;
      sigma_=pdg.sigma_;
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
    }

    /// Set the seed
    void set_seed(unsigned long int s) {
      r.set_seed(s);
    }

    /// Sample from the specified density
    virtual double operator()() const {
      return gsl_ran_lognormal(&r,mu_,sigma_);
    }
    
    /// The normalized density 
    virtual double pdf(double x) const {
      if (x<0.0) {
	return 0.0;
      }
      return gsl_ran_lognormal_pdf(x,mu_,sigma_);
    }
    
    /// The log of the normalized density 
    virtual double log_pdf(double x) const {
      if (x<0.0) {
	return 0.0;
      }
      return log(gsl_ran_lognormal_pdf(x,mu_,sigma_));
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
	O2SCL_ERR2("Requested cdf inverse outside of [0,1] in ",
		   "prob_dens_lognormal::invert_cdf().",exc_einval);
      }
      return gsl_cdf_lognormal_Pinv(in_cdf,mu_,sigma_);
    }

    /// The inverse cumulative distribution function
    virtual double entropy() const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Parameters not set in prob_dens_lognormal::",
		   "entropy().",exc_einval);
      }
      return 0.5+0.5*log(2.0*o2scl_const::pi*sigma_*sigma_)+mu_;
    }
    
  };
  
  /** \brief Probability density function based on a histogram

      This class is experimental.
  */
  class prob_dens_hist : public prob_dens_frange {
    
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

  protected:
    
    /// Search through the partial sums
    search_vec<ubvector> sv;
  
    /// Number of original histogram bins
    size_t n;
  
    /** \brief Normalized partial sum of histogram bins
	
	This vector has size \ref n plus one.
    */
    ubvector sum;
  
    /** \brief Vector specifying original histogram bins
	
	This vector has size \ref n plus one.
    */
    ubvector range;
  
    /// Random number generator
    mutable rng_gsl rng;
  
  public:
  
    prob_dens_hist();
  
    ~prob_dens_hist();

    /// Initialize with histogram \c h
    void init(hist &h);

    /// Get reference to partial sums
    const ubvector &partial_sums() { return sum; }
    
    /// Get reference to bin ranges
    const ubvector &bin_ranges() { return range; }
    
    /// Generate a sample
    virtual double operator()() const;

    /// Lower limit of the range
    virtual double lower_limit() const;
    
    /// Uower limit of the range
    virtual double upper_limit() const;

    /// The normalized density 
    virtual double pdf(double x) const;
    
    /// The log of the normalized density 
    virtual double log_pdf(double x) const;
    
    /// Cumulative distribution function (from the lower tail)
    virtual double cdf(double x) const;

    /// Inverse cumulative distribution function (from the lower tail)
    virtual double invert_cdf(double x) const;

    /// Inverse cumulative distribution function (from the lower tail)
    virtual double entropy() const {
      return 0.0;
    }
    
  };

  /** \brief A multi-dimensional probability density function
      
      This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class prob_dens_mdim {
    
  public:
  
  /// Return the dimensionality
  virtual size_t dim() const {
    return 0;
  }
  
  /// The normalized density 
  virtual double pdf(const vec_t &x) const {
    return 0.0;
  }
  
  /// The log of the normalized density 
  virtual double log_pdf(const vec_t &x) const {
    return 0.0;
  }
  
  /// Sample the distribution
  virtual void operator()(vec_t &x) const {
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
    double ret=1.0;
    for(size_t i=0;i<list.size();i++) ret*=list[i].pdf(x[i]);
    return log(ret);
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
      c = pdf(r \cos \theta + \mu_x, r \sin \theta + \mu_y )
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

  /// The covariance
  double rho;
  
  public:
  
  prob_dens_mdim_biv_gaussian() {
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

      This class is experimental.

      \future Create alternate versions based on other
      decompositions?
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

  /// Normalization factor
  double norm;

  /// Number of dimensions
  size_t ndim;

  /// Temporary storage 1
  mutable vec_t q;

  /// Temporary storage 2
  mutable vec_t vtmp;

  /// Standard normal
  o2scl::prob_dens_gaussian pdg;
    
  public:
  
  /// The dimensionality
  virtual size_t dim() const {
    return ndim;
  }

  /// Create an empty distribution
  prob_dens_mdim_gaussian() {
    ndim=0;
  }

  /** \brief Create a distribution from the covariance matrix
   */
  prob_dens_mdim_gaussian(size_t p_ndim, vec_t &p_peak, mat_t &covar) {
    set(p_ndim,p_peak,covar);
  }

  /** \brief Set the peak and covariance matrix for the distribution
   */
  void set(size_t p_ndim, vec_t &p_peak, mat_t &covar) {
    if (p_ndim==0) {
      O2SCL_ERR("Zero dimension in prob_dens_mdim_gaussian::set().",
		o2scl::exc_einval);
    }
    ndim=p_ndim;
    norm=1.0;
    peak.resize(ndim);
    for(size_t i=0;i<ndim;i++) peak[i]=p_peak[i];
    q.resize(ndim);
    vtmp.resize(ndim);

    // Perform the Cholesky decomposition of the covariance matrix
    chol=covar;
    o2scl_linalg::cholesky_decomp(ndim,chol);
      
    // Find the inverse
    covar_inv=chol;
    o2scl_linalg::cholesky_invert<mat_t>(ndim,covar_inv);

    // Force chol to be lower triangular and compute the determinant
    double det=1.0;
    for(size_t i=0;i<ndim;i++) {
      det*=chol(i,i);
      for(size_t j=0;j<ndim;j++) {
	if (i<j) chol(i,j)=0.0;
      }
    }
    det*=det;

    // Compute normalization
    norm=pow(2.0*o2scl_const::pi,-((double)ndim)/2.0)/sqrt(det);
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
      probability distribution based on a Gaussian process which
      includes noise
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
    vec_t prod(n_init), mean(n_dim);
    boost::numeric::ublas::axpy_prod(inv_KXX,y,prod,true);
    boost::numeric::ublas::axpy_prod(KXsX,prod,mean,true);
    
    // Compute the covariance matrix
    mat_t covar(n_dim,n_dim), prod2(n_init,n_dim), prod3(n_dim,n_dim);
    boost::numeric::ublas::axpy_prod(inv_KXX,KXXs,prod2,true);
    boost::numeric::ublas::axpy_prod(KXsX,prod2,prod3,true);
    covar=KXsXs-prod3;
    
    // Now use set() in the parent class
    this->set(n_dim,mean,covar);
    
  }

  /// The normalized density 
  virtual double pdf(const vec_t &x) const {
    if (ndim==0) {
      O2SCL_ERR2("Distribution not set in prob_dens_mdim_gaussian::",
		 "pdf().",o2scl::exc_einval);
    }
    double ret=norm;
    for(size_t i=0;i<ndim;i++) q[i]=x[i]-peak[i];
    vtmp=prod(covar_inv,q);
    ret*=exp(-0.5*inner_prod(q,vtmp));
    return ret;
  }

  /// The log of the normalized density 
  virtual double log_pdf(const vec_t &x) const {
    if (ndim==0) {
      O2SCL_ERR2("Distribution not set in prob_dens_mdim_gaussian::",
		 "pdf().",o2scl::exc_einval);
    }
    double ret=log(norm);
    for(size_t i=0;i<ndim;i++) q[i]=x[i]-peak[i];
    vtmp=prod(covar_inv,q);
    ret+=-0.5*inner_prod(q,vtmp);
    return ret;
  }

  /// Sample the distribution
  virtual void operator()(vec_t &x) const {
    if (ndim==0) {
      O2SCL_ERR2("Distribution not set in prob_dens_mdim_gaussian::",
		 "operator().",o2scl::exc_einval);
    }
    for(size_t i=0;i<ndim;i++) q[i]=pdg();
    vtmp=prod(chol,q);
    for(size_t i=0;i<ndim;i++) x[i]=peak[i]+vtmp[i];
    return;
  }

    };

  /** \brief A multi-dimensional conditional probability density function
      
      This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class prob_cond_mdim {
    
  public:
  
  /// The dimensionality
  virtual size_t dim() const {
    return 0;
  }
  
  /// The conditional probability
  virtual double pdf(const vec_t &x, const vec_t &x2) const=0;
  
  /// The log of the conditional probability
  virtual double log_pdf(const vec_t &x, const vec_t &x2) const=0;
  
  /// Sample the distribution
  virtual void operator()(const vec_t &x, vec_t &x2) const=0;

  /** \brief Sample the distribution and return the 
      log of the Metropolis-Hastings ratio
  */
  virtual double metrop_hast(const vec_t &x, vec_t &x2) const {
    operator()(x,x2);
    return log_pdf(x,x2)-log_pdf(x2,x);
  }
  
  };

  /** \brief A constrained random walk in the shape of 
      a hypercube

      \comment
      I had previously used std::uniform_real_distribution
      instead of rng_gsl, but this caused problems with
      intel compilers.

      The Metropolis step is
      P(x')g(xt|x')
      -------------
      P(xt)g(x'|xt)

      and if we do not include the g ratio, then the edges
      will be undersampled because we don't properly include
      the rejections 

      For g(x'|xt), if xt is near the edge, then the cond. prob.
      is larger, thus the g ratio is smaller than 1, encouraging
      more rejections.
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
  mutable rng_gsl rg;
  
  /** \brief Internal set function

      \comment
      This can't be virtual because it needs to be called
      by the constructor
      \endcomment
  */
  int set_internal(vec_t &step, vec_t &low, vec_t &high) {
    u_step.resize(step.size());
    u_low.resize(step.size());
    u_high.resize(step.size());
    for(size_t i=0;i<step.size();i++) {
      u_step[i]=step[i];
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
    set_internal(step,low,high);
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
    set_internal(step,low,high);
    return 0;
  }

  /// The dimensionality
  virtual size_t dim() const {
    return u_step.size();
  }
  
  /// The conditional probability
  virtual double pdf(const vec_t &x, const vec_t &x2) const {
    double vol1=1.0;
    for(size_t i=0;i<u_step.size();i++) {
      if (fabs(x2[i]-x[i])>u_step[i]) return 0.0;
      if (x[i]-u_step[i]<u_low[i]) {
	if (x[i]+u_step[i]>u_high[i]) {
	  vol1*=u_high[i]-u_low[i];
	} else {
	  vol1*=x[i]+u_step[i]-u_low[i];
	}
      } else {
	if (x[i]+u_step[i]>u_high[i]) {
	  vol1*=u_high[i]-(x[i]-u_step[i]);
	} else {
	  vol1*=2.0*u_step[i];
	}
      }
    }
    return 1.0/vol1;
  }
  
  /// The log of the conditional probability 
  virtual double log_pdf(const vec_t &x, const vec_t &x2) const {
    return log(pdf(x,x2));
  }
  
  /// Sample the distribution
  virtual void operator()(const vec_t &x, vec_t &x2) const {
    size_t nv=u_step.size();
    for(size_t i=0;i<nv;i++) {
      if (x[i]<u_low[i] || x[i]>u_high[i]) {
	std::cout << "Herex: " << i << " " << x[i] << " "
		  << u_low[i] << " " << u_high[i] << std::endl;
	O2SCL_ERR("Input out of bounds in fixed_step::operator().",
		  o2scl::exc_einval);
      }
      x2[i]=x[i]+u_step[i]*(rg.random()*2.0-1.0);
      while (x2[i]<u_low[i] || x2[i]>u_high[i]) {
	x2[i]=x[i]+u_step[i]*(rg.random()*2.0-1.0);
      }
    }
    return;
  }
  
  };

  /** \brief A multi-dimensional conditional probability density function
      independent of the input

      This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class prob_cond_mdim_invar {

  protected:
  
  prob_dens_mdim<vec_t> &base;
  
  public:

  prob_cond_mdim_invar(prob_dens_mdim<vec_t> &out) : base(out) {
  }
  
  /// The dimensionality
  virtual size_t dim() const {
    return base.dim();
  }
  
  /// The conditional probability 
  virtual double pdf(const vec_t &x, const vec_t &x2) const {
    return base.pdf(x2);
  }
  
  /// The log of the conditional probability 
  virtual double log_pdf(const vec_t &x, const vec_t &x2) const {
    return base.log_pdf(x2);
  }
  
  /// Sample the distribution
  virtual void operator()(const vec_t &x, vec_t &x2) const {
    return base(x2);
  }
  
  };
  
  /** \brief A multi-dimensional Gaussian conditional probability
      density function

      This class is experimental.
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
  o2scl::prob_dens_gaussian pdg;
    
  public:

  /** \brief Create an empty distribution 
   */
  prob_cond_mdim_gaussian() {
    ndim=0;
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
      
    // Find the inverse
    covar_inv=chol;
    o2scl_linalg::cholesky_invert<mat_t>(ndim,covar_inv);
      
    // Force chol to be lower triangular and compute the determinant
    double det=1.0;
    for(size_t i=0;i<ndim;i++) {
      det*=chol(i,i);
      for(size_t j=0;j<ndim;j++) {
	if (i<j) chol(i,j)=0.0;
      }
    }
    det*=det;

    // Compute normalization
    norm=pow(2.0*o2scl_const::pi,-((double)ndim)/2.0)/sqrt(det);
  }

  /// The conditional probability 
  virtual double pdf(const vec_t &x, const vec_t &x2) const {
    if (ndim==0) {
      O2SCL_ERR2("Distribution not set in prob_cond_mdim_gaussian::",
		 "pdf().",o2scl::exc_einval);
    }
    double ret=norm;
    for(size_t i=0;i<ndim;i++) q[i]=x2[i]-x[i];
    vtmp=prod(covar_inv,q);
    ret*=exp(-0.5*inner_prod(q,vtmp));
    return ret;
  }

  /// The log of the conditional probability 
  virtual double log_pdf(const vec_t &x, const vec_t &x2) const {
    if (ndim==0) {
      O2SCL_ERR2("Distribution not set in prob_cond_mdim_gaussian::",
		 "pdf().",o2scl::exc_einval);
    }
    double ret=log(norm);
    for(size_t i=0;i<ndim;i++) q[i]=x2[i]-x[i];
    vtmp=prod(covar_inv,q);
    ret+=-0.5*inner_prod(q,vtmp);
    return ret;
  }

  /// Sample the distribution
  virtual void operator()(const vec_t &x, vec_t &x2) const {
    if (ndim==0) {
      O2SCL_ERR2("Distribution not set in prob_cond_mdim_gaussian::",
		 "operator().",o2scl::exc_einval);
    }
    for(size_t i=0;i<ndim;i++) q[i]=pdg();
    vtmp=prod(chol,q);
    for(size_t i=0;i<ndim;i++) x2[i]=x[i]+vtmp[i];
    return;
  }

    };

#ifdef O2SCL_NEVER_DEFINED  
  /** \brief A multidimensional normal distribution from
      a Gaussian process
      
      \future The linear algebra only works with ublas and is
      not optimized.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double>,
    class mat_col_t=boost::numeric::ublas::matrix_column<mat_t> >
    class prob_dens_mdim_gproc :
    public o2scl::prob_dens_mdim_gaussian<vec_t> {
    
  public:
  
    };
#endif
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
