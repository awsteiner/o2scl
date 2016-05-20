/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2016, Andrew W. Steiner
  
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

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/hist.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/search_vec.h>
#include <o2scl/cholesky.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A one-dimensional probability density function

      This class is experimental.

      \future Give functions for mean, median, mode, variance, etc?
   */
  class prob_dens_func {
    //: public funct {
    
  public:
    
    /// Sample from the specified density
    virtual double operator()() const=0;
    
    /// The normalized density 
    virtual double function(double x) const=0;
    
    /// The cumulative distribution function (from the lower tail)
    virtual double cdf(double x) const=0;
    
    /// The inverse cumulative distribution function
    virtual double invert_cdf(double cdf) const=0;

    /// Entropy of the distribution (\f$ - \int f \ln f \f$ )
    virtual double entropy() const=0;
    
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
    gsl_rng *r;
    
  public:
    
    /** \brief Create a standard normal distribution
     */
    prob_dens_gaussian() {
      cent_=0.0;
      sigma_=1.0;
      r=gsl_rng_alloc(gsl_rng_mt19937);
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
      r=gsl_rng_alloc(gsl_rng_mt19937);
    }

    virtual ~prob_dens_gaussian() {
      gsl_rng_free(r);
    }

    /// Copy constructor
  prob_dens_gaussian(const prob_dens_gaussian &pdg) : prob_dens_func() {
      cent_=pdg.cent_;
      sigma_=pdg.sigma_;
      r=gsl_rng_alloc(gsl_rng_mt19937);
    }

    /// Copy constructor with operator=
    prob_dens_gaussian &operator=(const prob_dens_gaussian &pdg) {
      // Check for self-assignment
      if (this==&pdg) return *this;
      cent_=pdg.cent_;
      sigma_=pdg.sigma_;
      return *this;
    }

    /// Set the seed
    void set_seed(unsigned long int s) { 
      gsl_rng_set(r,s);
    }

    /// Set the center
    void set_center(double cent) {
      cent_=cent;
    }

    /// Set the Gaussian width
    void set_sigma(double sigma) {
      if (sigma<0.0) {
	O2SCL_ERR2("Tried to set sigma negative",
		   "in prob_dens_gaussian::prob_dens_gaussian().",
		   exc_einval);
      }
      sigma_=sigma;
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
      return cent_+gsl_ran_gaussian(r,sigma_);
    }
    
    /// The normalized density 
    virtual double function(double x) const {
      if (sigma_<0.0) {
	O2SCL_ERR2("Width not set in prob_dens_gaussian::",
		   "function().",exc_einval);
      }
      return gsl_ran_gaussian_pdf(x-cent_,sigma_);
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
    gsl_rng *r;
    
  public:

    /** \brief Create a blank uniform distribution
     */
    prob_dens_uniform() {
      ll=1.0;
      ul=0.0;
      r=gsl_rng_alloc(gsl_rng_mt19937);
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
      r=gsl_rng_alloc(gsl_rng_mt19937);
    }

    virtual ~prob_dens_uniform() {
      gsl_rng_free(r);
    }

    /// Copy constructor
  prob_dens_uniform(const prob_dens_uniform &pdg) : prob_dens_frange() {
      ll=pdg.ll;
      ul=pdg.ul;
      r=gsl_rng_alloc(gsl_rng_mt19937);
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
      gsl_rng_set(r,s);
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
      return gsl_ran_flat(r,ll,ul);
    }
    
    /// The normalized density 
    virtual double function(double x) const {
      if (ll>ul) {
	O2SCL_ERR2("Limits not set in prob_dens_uniform::",
		   "function().",exc_einval);
      }
      if (x<ll || x>ul) return 0.0;
      return gsl_ran_flat_pdf(x,ll,ul);
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
    gsl_rng *r;
    
  public:

    /** \brief Create a blank lognormal distribution
     */
    prob_dens_lognormal() {
      sigma_=-1.0;
      mu_=0.0;
      r=gsl_rng_alloc(gsl_rng_mt19937);
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
      r=gsl_rng_alloc(gsl_rng_mt19937);
    }

    virtual ~prob_dens_lognormal() {
      gsl_rng_free(r);
    }

    /// Copy constructor
  prob_dens_lognormal(const prob_dens_lognormal &pdg) : prob_dens_positive() {
      mu_=pdg.mu_;
      sigma_=pdg.sigma_;
      r=gsl_rng_alloc(gsl_rng_mt19937);
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
      gsl_rng_set(r,s);
    }

    /// Sample from the specified density
    virtual double operator()() const {
      return gsl_ran_lognormal(r,mu_,sigma_);
    }
    
    /// The normalized density 
    virtual double function(double x) const {
      if (x<0.0) {
	return 0.0;
      }
      return gsl_ran_lognormal_pdf(x,mu_,sigma_);
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

    /// Generate a sample
    virtual double operator()() const;

    /// Lower limit of the range
    virtual double lower_limit() const;
    
    /// Uower limit of the range
    virtual double upper_limit() const;

    /// The normalized density 
    virtual double function(double x) const;
    
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
    
    /// Return the probability density
    virtual double function(vec_t &x) const=0;
    
  /// Sample the distribution
  virtual void operator()(vec_t &x) const=0;

  };

  /** \brief A multidimensional distribution formed by the product
      of several one-dimensional distributions
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class prob_dens_mdim_factor : public prob_dens_mdim<vec_t> {
    
    /// Vector of one-dimensional distributions
    std::vector<prob_dens_func> list;
    
  public:

  prob_dens_mdim_factor(std::vector<prob_dens_func> &p_list) {
    list=p_list;
  }
  
  /// Return the probability density
    virtual double function(vec_t &x) const {
      double ret=1.0;
      for(size_t i=0;i<list.size();i++) ret*=list[i].function(x[i]);
      return ret;
    }
  
  /// Sample the distribution
    virtual void operator()(vec_t &x) const {
      for(size_t i=0;i<list.size();i++) x[i]=list[i]();
      return;
    }
  
  };
  
  /** \brief A multi-dimensional probability density function

      This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double> >
    class prob_dens_mdim_gauss : public prob_dens_mdim<vec_t> {
    
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
  
  /** \brief Create a distribution from the covariance matrix
   */
    prob_dens_mdim_gauss(size_t p_ndim, vec_t &p_peak, mat_t &covar) {
      ndim=p_ndim;
      norm=1.0;
      peak.resize(ndim);
      for(size_t i=0;i<ndim;i++) peak[i]=p_peak[i];
      q.resize(ndim);
      vtmp.resize(ndim);

      // Perform the Cholesky decomposition
      chol=covar;
      o2scl_linalg::cholesky_decomp(ndim,chol);
      
      // Find the inverse
      covar_inv=chol;
      o2scl_linalg::cholesky_invert<mat_t>(ndim,covar_inv);
      
      // Force chol to be lower triangular
      for(size_t i=0;i<ndim;i++) {
	for(size_t j=0;j<ndim;j++) {
	  if (i<j) chol(i,j)=0.0;
	}
      }
    }

    /// Return the probability density
    virtual double function(vec_t &x) const {
      double ret=norm;
      for(size_t i=0;i<ndim;i++) q[i]=x[i]-peak[i];
      vtmp=prod(covar_inv,q);
      ret*=exp(-0.5*inner_prod(q,vtmp));
      return ret;
    }

    /// Sample the distribution
    virtual void operator()(vec_t &x) const {
      for(size_t i=0;i<ndim;i++) q[i]=pdg();
      vtmp=prod(chol,q);
      for(size_t i=0;i<ndim;i++) x[i]=peak[i]+vtmp[i];
      return;
    }

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
