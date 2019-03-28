/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#ifndef O2SCL_VEC_STATS_H
#define O2SCL_VEC_STATS_H

/** \file vec_stats.h
    \brief Statistical functions for vector types

    This file contains several function templates for computing
    statistics of vectors of double-precision data. It includes mean,
    median, variance, standard deviation, covariance, correlation, and
    other functions.

    No additional range checking is done on the vectors.

    \future Consider generalizing to other data types.
*/

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// \name Vector mean, standard deviation, and variance
  //@{
  /** \brief Compute the mean of the first \c n elements of a vector

      This function produces the same results
      as <tt>gsl_stats_mean()</tt>.

      If \c n is zero, this function will return zero.
  */
  template<class vec_t> double vector_mean(size_t n, const vec_t &data) {
    long double mean=0.0;
    for(size_t i=0;i<n;i++) {
      mean+=(data[i]-mean)/(i+1);
    }
    return mean;
  }

  /** \brief Compute the mean of all of the vector elements

      This function uses <tt>size()</tt> to determine the vector size
      and produces the same results as <tt>gsl_stats_mean()</tt>.

      If the vector size is zero, this function will return zero.
  */
  template<class vec_t> double vector_mean(const vec_t &data) {
    return vector_mean(data.size(),data);
  }

  /** \brief Compute variance with specified mean known in advance
      
      This function computes
      \f[
      \frac{1}{N} \sum_{i} \left( x_i - \mu \right)^2
      \f]
      where the value of \f$ \mu \f$ is given in \c mean. 

      This function produces the same results as
      <tt>gsl_stats_variance_with_fixed_mean()</tt>. IF \c n is zero,
      this function returns zero.
  */
  template<class vec_t>
    double vector_variance_fmean(size_t n, const vec_t &data, double mean) {
    long double var=0.0;
    for(size_t i=0;i<n;i++) {
      long double delta=(data[i]-mean);
      var+=(delta*delta-var)/(i+1);
    }
    return var;
  }

  /** \brief Compute variance with specified mean known in advance
      
      This function computes
      \f[
      \frac{1}{N} \sum_{i} \left( x_i - \mu \right)^2
      \f]
      where the value of \f$ \mu \f$ is given in \c mean. 

      This function produces the same results as
      <tt>gsl_stats_variance_with_fixed_mean()</tt>. If the vector
      size is zero, this function will return zero.
  */
  template<class vec_t>
    double vector_variance_fmean(const vec_t &data, double mean) {
    return vector_variance_fmean(data.size(),data,mean);
  }

  /** \brief Compute the variance with specified mean

      This function computes
      \f[
      \frac{1}{N-1} \sum_{i} \left( x_i - \mu \right)^2
      \f]
      where the value of \f$ \mu \f$ is given in \c mean.
      
      This function produces the same results
      as <tt>gsl_stats_variance_m</tt>.

      If \c n is 0 or 1, this function will call the error
      handler.
  */
  template<class vec_t>
    double vector_variance(size_t n, const vec_t &data, double mean) {

    if (n<2) {
      O2SCL_ERR2("Cannot compute variance with less than 2 elements",
		     " in vector_variance().",exc_einval);
    }

    double var=vector_variance_fmean<vec_t>(n,data,mean);
    return var*n/(n-1);
  }

  /** \brief Compute the variance with specified mean

      This function computes
      \f[
      \frac{1}{N-1} \sum_{i} \left( x_i - \mu \right)^2
      \f]
      where the value of \f$ \mu \f$ is given in \c mean.
      
      This function produces the same results
      as <tt>gsl_stats_variance_m</tt>.

      If \c n is 0 or 1, this function will call the error
      handler.
  */
  template<class vec_t>
    double vector_variance(const vec_t &data, double mean) {
    return vector_variance(data.size(),data,mean);
  }

  /** \brief Compute the variance 

      This function computes
      \f[
      \frac{1}{N-1} \sum_{i} \left( x_i - \mu \right)^2
      \f]
      where \f$ \mu \f$ is the mean computed with \ref vector_mean().
      
      This function produces the same results
      as <tt>gsl_stats_variance</tt>.

      If \c n is 0 or 1, this function will call the error handler.
  */
  template<class vec_t> double vector_variance(size_t n, const vec_t &data) {

    if (n<2) {
      O2SCL_ERR2("Cannot compute variance with less than 2 elements",
		     " in vector_variance().",exc_einval);
    }
    
    double mean=vector_mean<vec_t>(n,data);
    double var=vector_variance_fmean<vec_t>(n,data,mean);
    return var*n/(n-1);
  }

  /** \brief Compute the variance 

      This function computes
      \f[
      \frac{1}{N-1} \sum_{i} \left( x_i - \mu \right)^2
      \f]
      where \f$ \mu \f$ is the mean computed with \ref vector_mean().
      
      This function produces the same results
      as <tt>gsl_stats_variance</tt>.

      If \c n is 0 or 1, this function will call the error handler.
  */
  template<class vec_t> double vector_variance(const vec_t &data) {
    return vector_variance(data.size(),data);
  }

  /** \brief Standard deviation with specified mean known in advance

      This function computes
      \f[
      \sqrt{\frac{1}{N} \sum_{i} \left( x_i - \mu \right)^2}
      \f]
      where the value of \f$ \mu \f$ is given in \c mean. 

      This function produces the same results
      as <tt>gsl_stats_sd_with_fixed_mean()</tt>.

      If \c n is zero, this function will return zero without calling
      the error handler.
  */
  template<class vec_t>
    double vector_stddev_fmean(size_t n, const vec_t &data, 
			       double mean) {
    double sd=vector_variance_fmean<vec_t>(n,data,mean);
    return std::sqrt(sd);
  }

  /** \brief Standard deviation with specified mean known in advance

      This function computes
      \f[
      \sqrt{\frac{1}{N} \sum_{i} \left( x_i - \mu \right)^2}
      \f]
      where the value of \f$ \mu \f$ is given in \c mean. 

      This function produces the same results
      as <tt>gsl_stats_sd_with_fixed_mean()</tt>.

      If \c n is zero, this function will return zero without calling
      the error handler.
  */
  template<class vec_t>
    double vector_stddev_fmean(const vec_t &data, double mean) {
    return vector_stddev_fmean(data.size(),data,mean);
  }

  /** \brief Standard deviation with specified mean

      This function computes
      \f[
      \sqrt{\frac{1}{N-1} \sum_{i} \left( x_i - \mu \right)^2}
      \f]
      where \f$ \mu \f$ is the mean computed with \ref vector_mean().

      This function produces the same results
      as <tt>gsl_stats_sd()</tt>.

      If \c n is 0 or 1, this function will call the error handler.
  */
  template<class vec_t> double vector_stddev(size_t n, const vec_t &data) {
    
    if (n<2) {
      O2SCL_ERR2("Cannot compute std. dev. with less than 2 elements",
		     " in vector_stddev().",exc_einval);
    }
    
    double mean=vector_mean<vec_t>(n,data);
    double var=vector_variance_fmean<vec_t>(n,data,mean);
    return std::sqrt(var*n/(n-1));
  }

  /** \brief Standard deviation with specified mean

      This function computes
      \f[
      \sqrt{\frac{1}{N-1} \sum_{i} \left( x_i - \mu \right)^2}
      \f]
      where \f$ \mu \f$ is the mean computed with \ref vector_mean().

      This function produces the same results
      as <tt>gsl_stats_sd()</tt>.

      If \c n is 0 or 1, this function will call the error handler.
  */
  template<class vec_t> double vector_stddev(const vec_t &data) {
    return vector_stddev(data.size(),data);
  }

  /** \brief Standard deviation with specified mean

      This function computes
      \f[
      \sqrt{\frac{1}{N-1} \sum_{i} \left( x_i - \mu \right)^2}
      \f]
      where the value of \f$ \mu \f$ is given in \c mean. 

      This function produces the same results
      as <tt>gsl_stats_sd_m()</tt>.

      If \c n is 0 or 1, this function will call the error
      handler.
  */
  template<class vec_t> double vector_stddev(size_t n, const vec_t &data, 
					     double mean) {

    if (n<2) {
      O2SCL_ERR2("Cannot compute std. dev. with less than 2 elements",
		     " in vector_stddev().",exc_einval);
    }
    
    double sd=vector_variance_fmean<vec_t>(n,data,mean);
    return std::sqrt(sd*n/(n-1));
  }

  /** \brief Standard deviation with specified mean

      This function computes
      \f[
      \sqrt{\frac{1}{N-1} \sum_{i} \left( x_i - \mu \right)^2}
      \f]
      where the value of \f$ \mu \f$ is given in \c mean. 

      This function produces the same results
      as <tt>gsl_stats_sd_m()</tt>.

      If \c n is 0 or 1, this function will call the error
      handler.
  */
  template<class vec_t> double vector_stddev(const vec_t &data, double mean) {
    return vector_stddev(data.size(),data,mean);
  }
  
  /** \brief Absolute deviation from the specified mean

      This function computes
      \f[
      \sum_i | x_i - \mu |
      \f]
      where the value of \f$ \mu \f$ is given in \c mean. 

      This function produces the same results
      as <tt>gsl_stats_absdev_m()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t> double vector_absdev(size_t n, const vec_t &data, 
					     double mean) {
    
    if (n==0) return 0.0;

    long double sum=0.0;
    for(size_t i=0;i<n;i++) {
      sum+=fabs(data[i]-mean);
    }
    return sum/n;
  }
  //@}

  /// \name Vector absolute deviation, skewness, and kurtosis
  //@{
  /** \brief Absolute deviation from the specified mean

      This function computes
      \f[
      \sum_i | x_i - \mu |
      \f]
      where the value of \f$ \mu \f$ is given in \c mean. 

      This function produces the same results
      as <tt>gsl_stats_absdev_m()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t> double vector_absdev(const vec_t &data, 
					     double mean) {
    return vector_absdev(data.size(),data,mean);
  }
  
  /** \brief Absolute deviation from the computed mean

      This function computes
      \f[
      \sum_i | x_i - \mu |
      \f]
      where the value of \f$ \mu \f$ is mean as computed
      from \ref vector_mean().

      This function produces the same results
      as <tt>gsl_stats_absdev()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t>
    double vector_absdev(size_t n, const vec_t &data) {
    double mean=vector_mean<vec_t>(n,data);
    return vector_absdev(n,data,mean);
  }

  /** \brief Absolute deviation from the computed mean

      This function computes
      \f[
      \sum_i | x_i - \mu |
      \f]
      where the value of \f$ \mu \f$ is mean as computed
      from \ref vector_mean().

      This function produces the same results
      as <tt>gsl_stats_absdev()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t>
    double vector_absdev(const vec_t &data) {
    return vector_absdev(data.size(),data);
  }

  /** \brief Skewness with specified mean and standard deviation

      This function computes 
      \f[
      \frac{1}{N} \sum_i \left[ 
      \frac{ \left(x_i - \mu \right)}{ \sigma }\right]^3
      \f]
      where the values of \f$ \mu \f$ and \f$ \sigma \f$ 
      are given in \c mean and \c stddev.

      This function produces the same results
      as <tt>gsl_stats_skew_m_sd()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t> double vector_skew(size_t n, const vec_t &data, 
					   double mean, double stddev) {
    long double skew=0.0;
    for(size_t i=0;i<n;i++) {
      long double x=(data[i]-mean)/stddev;
      skew+=(x*x*x-skew)/(i+1);
    }
    return skew;
  }

  /** \brief Skewness with specified mean and standard deviation

      This function computes 
      \f[
      \frac{1}{N} \sum_i \left[ 
      \frac{ \left(x_i - \mu \right)}{ \sigma }\right]^3
      \f]
      where the values of \f$ \mu \f$ and \f$ \sigma \f$ 
      are given in \c mean and \c stddev.

      This function produces the same results
      as <tt>gsl_stats_skew_m_sd()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t> double vector_skew(const vec_t &data, 
					   double mean, double stddev) {
    return vector_skew(data.size(),data,mean,stddev);
  }

  /** \brief Skewness with computed mean and standard deviation

      This function computes 
      \f[
      \frac{1}{N} \sum_i \left[ 
      \frac{ \left(x_i - \mu \right)}{ \sigma }\right]^3
      \f]
      where the values of \f$ \mu \f$ and \f$ \sigma \f$ 
      are computed using \ref vector_mean() and \ref vector_stddev().

      This function produces the same results
      as <tt>gsl_stats_skew()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t> double vector_skew(size_t n, const vec_t &data) {
    double mean=vector_mean<vec_t>(n,data);
    double sd=vector_stddev<vec_t>(n,data,mean);
    return vector_skew(n,data,mean,sd);
  }

  /** \brief Skewness with computed mean and standard deviation

      This function computes 
      \f[
      \frac{1}{N} \sum_i \left[ 
      \frac{ \left(x_i - \mu \right)}{ \sigma }\right]^3
      \f]
      where the values of \f$ \mu \f$ and \f$ \sigma \f$ 
      are computed using \ref vector_mean() and \ref vector_stddev().

      This function produces the same results
      as <tt>gsl_stats_skew()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t> double vector_skew(const vec_t &data) {
    return vector_skew(data.size(),data);
  }

  /** \brief Kurtosis with specified mean and standard deviation

      This function computes 
      \f[
      -3 + \frac{1}{N} \sum_i \left[ 
      \frac{ \left(x_i - \mu \right)}{ \sigma }\right]^4
      \f]
      where the values of \f$ \mu \f$ and \f$ \sigma \f$ 
      are given in \c mean and \c stddev.

      This function produces the same results
      as <tt>gsl_stats_kurtosis_m_sd()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t>
    double vector_kurtosis(size_t n, const vec_t &data, double mean,
			   double stddev) {
    long double avg=0.0;
    for(size_t i=0;i<n;i++) {
      long double x=(data[i]-mean)/stddev;
      avg+=(x*x*x*x-avg)/(i+1);
    }
    return avg-3.0;
  }

  /** \brief Kurtosis with specified mean and standard deviation

      This function computes 
      \f[
      -3 + \frac{1}{N} \sum_i \left[ 
      \frac{ \left(x_i - \mu \right)}{ \sigma }\right]^4
      \f]
      where the values of \f$ \mu \f$ and \f$ \sigma \f$ 
      are given in \c mean and \c stddev.

      This function produces the same results
      as <tt>gsl_stats_kurtosis_m_sd()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t>
    double vector_kurtosis(const vec_t &data, double mean,
			   double stddev) {
    return vector_kurtosis(data.size(),data,mean,stddev);
  }

  /** \brief Kurtosis with computed mean and standard deviation

      This function computes 
      \f[
      -3 + \frac{1}{N} \sum_i \left[ 
      \frac{ \left(x_i - \mu \right)}{ \sigma }\right]^4
      \f]
      where the values of \f$ \mu \f$ and \f$ \sigma \f$ 
      are computed using \ref vector_mean() and \ref vector_stddev().

      This function produces the same results
      as <tt>gsl_stats_kurtosis()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t> double vector_kurtosis(size_t n, const vec_t &data) {
    double mean=vector_mean<vec_t>(n,data);
    double sd=vector_stddev<vec_t>(n,data,mean);
    return vector_kurtosis(n,data,mean,sd);
  }

  /** \brief Kurtosis with computed mean and standard deviation

      This function computes 
      \f[
      -3 + \frac{1}{N} \sum_i \left[ 
      \frac{ \left(x_i - \mu \right)}{ \sigma }\right]^4
      \f]
      where the values of \f$ \mu \f$ and \f$ \sigma \f$ 
      are computed using \ref vector_mean() and \ref vector_stddev().

      This function produces the same results
      as <tt>gsl_stats_kurtosis()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t> double vector_kurtosis(const vec_t &data) {
    return vector_kurtosis(data.size(),data);
  }
  //@}

  /// \name Other vector functions
  //@{
  /** \brief Compute the covariance of two vectors
      
      This function computes
      \f[
      \frac{1}{n-1} \sum_i \left(x_i - {\bar{x}}\right)
      \left(y_i - {\bar{y}}\right)
      \f]
      where \f$ {\bar{x}} \f$ and \f$ {\bar{y}} \f$ are specified
      in \c mean1 and \c mean2, respectively.

      This function produces the same results
      as <tt>gsl_stats_covariance_m()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t, class vec2_t>
    double vector_covariance(size_t n, const vec_t &data1, const vec2_t &data2,
			     double mean1, double mean2) {
    double covar=0.0;
    for(size_t i=0;i<n;i++) {
      double delta1=(data1[i]-mean1);
      double delta2=(data2[i]-mean2);
      covar+=(delta1*delta2-covar)/(i+1);
    }
    return covar*n/(n-1);
  }

  /** \brief Compute the covariance of two vectors

      This function computes
      \f[
      \frac{1}{n-1} \sum_i \left(x_i - {\bar{x}}\right)
      \left(y_i - {\bar{y}}\right)
      \f]
      where \f$ {\bar{x}} \f$ and \f$ {\bar{y}} \f$ are specified
      in \c mean1 and \c mean2, respectively.

      This function produces the same results
      as <tt>gsl_stats_covariance_m()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t, class vec2_t>
    double vector_covariance(const vec_t &data1, const vec2_t &data2,
			     double mean1, double mean2) {
    return vector_covariance(data1.size(),data1,data2,mean1,mean2);
  }

  /** \brief Compute the covariance of two vectors

      This function computes
      \f[
      \frac{1}{n-1} \sum_i \left(x_i - {\bar{x}}\right)
      \left(y_i - {\bar{y}}\right)
      \f]
      where \f$ {\bar{x}} \f$ and \f$ {\bar{y}} \f$ are 
      the averages of \c data1 and \c data2 and are computed
      automatically using \ref vector_mean().

      This function produces the same
      results as <tt>gsl_stats_covariance()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t, class vec2_t>
    double vector_covariance(size_t n, const vec_t &data1, 
			     const vec2_t &data2) {
    double covar=0.0;
    double mean1=vector_mean<vec_t>(n,data1);
    double mean2=vector_mean<vec_t>(n,data2);
    for(size_t i=0;i<n;i++) {
      long double delta1=(data1[i]-mean1);
      long double delta2=(data2[i]-mean2);
      covar+=(delta1*delta2-covar)/(i+1);
    }
    return covar*n/(n-1);
  }
  
  /** \brief Compute the covariance of two vectors

      This function computes
      \f[
      \frac{1}{n-1} \sum_i \left(x_i - {\bar{x}}\right)
      \left(y_i - {\bar{y}}\right)
      \f]
      where \f$ {\bar{x}} \f$ and \f$ {\bar{y}} \f$ are 
      the averages of \c data1 and \c data2 and are computed
      automatically using \ref vector_mean().

      This function produces the same
      results as <tt>gsl_stats_covariance()</tt>.

      If \c n is zero, this function will return zero
      without calling the error handler.
  */
  template<class vec_t, class vec2_t>
    double vector_covariance(const vec_t &data1, 
			     const vec2_t &data2) {
    return vector_covariance(data1.size(),data1,data2);
  }
  
  /** \brief Pearson's correlation

      This function computes the Pearson correlation coefficient 
      between \c data1 and \c data2 .
      
      This function produces the same
      results as <tt>gsl_stats_correlation()</tt>.

      \comment
      r = cov(x, y) / (\Hat\sigma_x \Hat\sigma_y)
      = {1/(n-1) \sum (x_i - \Hat x) (y_i - \Hat y)
      \over
      \sqrt{1/(n-1) \sum (x_i - \Hat x)^2} \sqrt{1/(n-1) 
      \sum (y_i - \Hat y)^2}
      }
      \endcomment

      If \c n is zero, this function will call the error handler.
  */
  template<class vec_t, class vec2_t>
    double vector_correlation(size_t n, const vec_t &data1, 
			      const vec2_t &data2) {
    size_t i;

    if (n<1) {
      O2SCL_ERR2("Cannot compute correlation with no elements",
		     " in vector_correlation().",exc_einval);
    }

    double sum_xsq=0.0;
    double sum_ysq=0.0;
    double sum_cross=0.0;
    double ratio;
    double delta_x, delta_y;
    double mean_x, mean_y;
    double r;

    /*
     * Compute:
     * sum_xsq = Sum [ (x_i - mu_x)^2 ],
     * sum_ysq = Sum [ (y_i - mu_y)^2 ] and
     * sum_cross = Sum [ (x_i - mu_x) * (y_i - mu_y) ]
     * using the above relation from Welford's paper
     */

    mean_x=data1[0];
    mean_y=data2[0];

    for (i=1; i < n; ++i) {
      ratio=i / (i + 1.0);
      delta_x=data1[i] - mean_x;
      delta_y=data2[i] - mean_y;
      sum_xsq += delta_x * delta_x * ratio;
      sum_ysq += delta_y * delta_y * ratio;
      sum_cross += delta_x * delta_y * ratio;
      mean_x += delta_x / (i + 1.0);
      mean_y += delta_y / (i + 1.0);
    }
    
    r=sum_cross / (std::sqrt(sum_xsq) * std::sqrt(sum_ysq));
    
    return r;
  }

  /** \brief Pearson's correlation

      This function computes the Pearson correlation coefficient 
      between \c data1 and \c data2 .
      
      This function produces the same
      results as <tt>gsl_stats_correlation()</tt>.

      \comment
      r = cov(x, y) / (\Hat\sigma_x \Hat\sigma_y)
      = {1/(n-1) \sum (x_i - \Hat x) (y_i - \Hat y)
      \over
      \sqrt{1/(n-1) \sum (x_i - \Hat x)^2} \sqrt{1/(n-1) 
      \sum (y_i - \Hat y)^2}
      }
      \endcomment

      If \c n is zero, this function will call the error handler.
  */
  template<class vec_t, class vec2_t>
    double vector_correlation(const vec_t &data1, 
			      const vec2_t &data2) {
    return vector_correlation(data1.size(),data1,data2);
  }

  /** \brief The pooled variance of two vectors

      This function computes
      \f[
      s_{p}^2 = \frac{(n_1-1)s_1^2+(n_2-1)s_2^2}{n_1+n_2-2}
      \f]
      where \f$ n_i \f$ is the number of elements in vector \f$ i \f$
      and \f$ s_i^2 \f$ is the variance of vector \f$ i \f$. 

      From http://en.wikipedia.org/wiki/Pooled_variance, "Under the
      assumption of equal population variances, the pooled sample
      variance provides a higher precision estimate of variance than
      the individual sample variances."

      This function produces the same
      results as <tt>gsl_stats_pvariance()</tt>.
  */
  template<class vec_t, class vec2_t>
    double vector_pvariance(size_t n1, const vec_t &data1, 
			    size_t n2, const vec2_t &data2) {
    double var1=vector_variance<vec_t>(n1,data1);
    double var2=vector_variance<vec2_t>(n2,data2);
    return (((n1-1)*var1)+((n2-1)*var2))/(n1+n2-2);
  }

  /** \brief The pooled variance of two vectors

      This function computes
      \f[
      s_{p}^2 = \frac{(n_1-1)s_1^2+(n_2-1)s_2^2}{n_1+n_2-2}
      \f]
      where \f$ n_i \f$ is the number of elements in vector \f$ i \f$
      and \f$ s_i^2 \f$ is the variance of vector \f$ i \f$. 

      From http://en.wikipedia.org/wiki/Pooled_variance, "Under the
      assumption of equal population variances, the pooled sample
      variance provides a higher precision estimate of variance than
      the individual sample variances."

      This function produces the same
      results as <tt>gsl_stats_pvariance()</tt>.
  */
  template<class vec_t, class vec2_t>
    double vector_pvariance(const vec_t &data1, 
			    const vec2_t &data2) {
    return vector_pvariance(data1.size(),data1,data2.size(),data2);
  }

  /** \brief Quantile from sorted data (ascending only)

      This function returns the quantile \c f of data which
      has already been sorted in ascending order. The quantile,
      \f$ q \f$ , is
      found by interpolation using 
      \f[
      q = \left(1-\delta\right) x_i \delta x_{i+1}
      \f]
      where \f$ i = \mathrm{floor}[ (n-1)f ] \f$ and 
      \f$ \delta = (n-1)f -i \f$ .

      This function produces the same
      results as <tt>gsl_stats_quantile_from_sorted_data()</tt>.

      No checks are made to ensure the data is sorted, or to ensure
      that \f$ 0 \leq 0 \leq 1 \f$. If \c n is zero, this function
      will return zero without calling the error handler.
  */
  template<class vec_t>
    double vector_quantile_sorted(size_t n, const vec_t &data, 
				  const double f) {

    double index=f*(n-1);
    size_t lhs=((size_t)index);
    double delta=index-lhs;
    if (n==0) return 0.0;
    if (lhs==n-1) return data[lhs];
    return (1-delta)*data[lhs]+delta*data[lhs+1];
  }
  
  /** \brief Quantile from sorted data (ascending only)

      This function returns the quantile \c f of data which
      has already been sorted in ascending order. The quantile,
      \f$ q \f$ , is
      found by interpolation using 
      \f[
      q = \left(1-\delta\right) x_i \delta x_{i+1}
      \f]
      where \f$ i = \mathrm{floor}[ (n-1)f ] \f$ and 
      \f$ \delta = (n-1)f -i \f$ .

      This function produces the same
      results as <tt>gsl_stats_quantile_from_sorted_data()</tt>.

      No checks are made to ensure the data is sorted, or to ensure
      that \f$ 0 \leq 0 \leq 1 \f$. If \c n is zero, this function
      will return zero without calling the error handler.
  */
  template<class vec_t>
    double vector_quantile_sorted(const vec_t &data, const double f) {
    return vector_quantile_sorted<vec_t>(data.size(),data,f);
  }
  
  /** \brief Return the median of sorted (ascending or descending) data

      This function returns the median of sorted data (either
      ascending or descending), assuming the data has already been
      sorted. When the data set has an odd number of elements, the
      median is the value of the element at index \f$ (n-1)/2 \f$,
      otherwise, the median is taken to be the average of the elements
      at indices \f$ (n-1)/2 \f$ and \f$ n/2 \f$ .

      This function produces the same
      results as <tt>gsl_stats_median_from_sorted_data()</tt>.

      No checks are made to ensure the data is sorted. If \c n is
      zero, this function will return zero without calling the error
      handler.
  */
  template<class vec_t>
    double vector_median_sorted(size_t n, const vec_t &data) {
    
    if (n==0) return 0.0;
    
    size_t lhs=(n-1)/2;
    size_t rhs=n/2;
    
    if (lhs==rhs) return data[lhs];

    return (data[lhs]+data[rhs])/2.0;
  }

  /** \brief Return the median of sorted (ascending or descending) data

      This function returns the median of sorted data (either
      ascending or descending), assuming the data has already been
      sorted. When the data set has an odd number of elements, the
      median is the value of the element at index \f$ (n-1)/2 \f$,
      otherwise, the median is taken to be the average of the elements
      at indices \f$ (n-1)/2 \f$ and \f$ n/2 \f$ .

      This function produces the same
      results as <tt>gsl_stats_median_from_sorted_data()</tt>.

      No checks are made to ensure the data is sorted. If \c n is
      zero, this function will return zero without calling the error
      handler.
  */
  template<class vec_t>
    double vector_median_sorted(const vec_t &data) {
    return vector_median_sorted<vec_t>(data.size(),data);
  }

  /** \brief Compute the chi-squared statistic

      This function computes
      \f[
      \sum_i \left( \frac{\mathrm{obs}_i - \mathrm{exp}_i}
      {\mathrm{err}_i}\right)^2
      \f]
      where \f$ \mathrm{obs} \f$ are the observed values,
      \f$ \mathrm{exp} \f$ are the expected values, and 
      \f$ \mathrm{err} \f$ are the errors.
   */
  template<class vec_t, class vec2_t, class vec3_t>
    double vector_chi_squared(size_t n, const vec_t &obs, const vec2_t &exp,
			      const vec3_t &err) {
    double chi2=0.0;
    for(size_t i=0;i<n;i++) {
      chi2+=pow((obs[i]-exp[i])/err[i],2.0);
    }
    return chi2;
  }

  /** \brief Compute the chi-squared statistic

      This function computes
      \f[
      \sum_i \left( \frac{\mathrm{obs}_i - \mathrm{exp}_i}
      {\mathrm{err}_i}\right)^2
      \f]
      where \f$ \mathrm{obs} \f$ are the observed values,
      \f$ \mathrm{exp} \f$ are the expected values, and 
      \f$ \mathrm{err} \f$ are the errors.
   */
  template<class vec_t, class vec2_t, class vec3_t>
    double vector_chi_squared(const vec_t &obs, const vec2_t &exp,
			      const vec3_t &err) {
    return vector_chi_squared<vec_t,vec2_t,vec3_t>(obs.size(),obs,exp,err);
  }

  /** \brief Optimal bin size using Scott's method for the
      first <tt>n</tt> elements
  */
  template<class vec_t> double vector_bin_size_scott
    (size_t n, const vec_t &v) {
    if (n<=1) return 0.0;
    double ret=3.5*vector_stddev(n,v)/cbrt(((double)n));
    return ret;
  }

  /** \brief Optimal bin size using Scott's method
      
      This function computes the optimal bin size \f$ \Delta_b \f$ of 
      a histogram using the expression
      \f[
      \Delta_b = \frac{3.5 \sigma}{n^{1/3}}
      \f]

      From \ref Scott79 .

      \note If <tt>n</tt> is less than or equal to 1, this
      function returns 0.0 without calling the error handler.
  */
  template<class vec_t> double vector_bin_size_scott
    (const vec_t &v) {
    return vector_bin_size_scott(v.size(),v);
  }
  
  /** \brief Obtain a quantile from a sorted vector

      This is a generic version of 
      <tt>gsl_stats_quantile_from_sorted_data()</tt>.

      If <tt>f</tt> is less than 0 or greater than 1, the error
      handler is called. If <tt>n</tt> is zero, this function returns
      zero without calling the error handler.
   */
  template<class vec_t> double vector_sorted_quantile
    (size_t n, const vec_t &v, double f) {

    if (f<0.0 || f>1.0) {
      O2SCL_ERR("Invalid fraction for vector_sorted_quantile",
		o2scl::exc_einval);
    }
    
    double index=f*(n-1);
    size_t lhs=(int)index;
    double delta=index-lhs;
    double result;
    
    if (n == 0) {
      return 0.0;
    }
    
    if (lhs == n-1) {
      result=v[lhs];
    } else {
      result=(1-delta)*v[lhs]+delta*v[(lhs+1)];
    }
    
    return result;
  }
  
  /** \brief Optimal bin size using the Freedman-Diaconis rule
      for the first <tt>n</tt> elements      
  */
  template<class vec_t> double vector_bin_size_freedman
    (size_t n, vec_t &v) {
    vector_sort<vec_t,double>(n,v);
    double ret=2.0*(vector_sorted_quantile(n,v,0.75)-
		    vector_sorted_quantile(n,v,0.25))/cbrt(((double)n));
    return ret;
  }
  
  /** \brief Optimal bin size using the Freedman-Diaconis rule
      
      This function computes the optimal bin size \f$ \Delta_b \f$ of 
      a histogram using the expression
      \f[
      \Delta_b = \frac{2\left(q_{0.75}-q_{0.25}\right)}{n^{1/3}}
      \f]
      where \f$ q_{i} \f$ is the \f$ i \f$ quantile 
      of the data (note this is quantile not quartile).
      This function sorts the vector in order to obtain
      the result.

      From \ref Freedman81 .

      \note If <tt>n</tt> is less than or equal to 1, this
      function returns 0.0 without calling the error handler.
  */
  template<class vec_t> double vector_bin_size_freedman
    (vec_t &v) {
    return vector_bin_size_freedman(v.size(),v);
  }
  //@}

  /// \name Weighted vector mean, standard deviation, and variance
  //@{
  /** \brief Compute the mean of weighted data

      This function computes 
      \f[
      \left( \sum_i w_i x_i \right) \left( \sum_i w_i \right)^{-1}
      \f]

      This function produces the same results
      as <tt>gsl_stats_wmean()</tt>.

      \comment
      M(n) = M(n-1) + (data[n] - M(n-1)) (w(n)/(W(n-1) + w(n)))
      W(n) = W(n-1) + w(n)
      \endcomment
  */
  template<class vec_t, class vec2_t>
    double wvector_mean(size_t n, const vec_t &data, const vec2_t &weights) {

    long double wmean=0.0;
    long double W=0.0;
    for(size_t i=0;i<n;i++) {
      double wi=weights[i];
      if (wi>0.0) {
	W+=wi;
	wmean+=(data[i]-wmean)*(wi/W);
      }
    }
    
    return wmean;
  }

  /** \brief Compute the mean of weighted data

      This function computes 
      \f[
      \left( \sum_i w_i x_i \right) \left( \sum_i w_i \right)^{-1}
      \f]

      This function produces the same results
      as <tt>gsl_stats_wmean()</tt>.

      \comment
      M(n) = M(n-1) + (data[n] - M(n-1)) (w(n)/(W(n-1) + w(n)))
      W(n) = W(n-1) + w(n)
      \endcomment
  */
  template<class vec_t, class vec2_t>
    double wvector_mean(const vec_t &data, const vec2_t &weights) {
    return wvector_mean<vec_t,vec2_t>(data.size(),data,weights);
  }

  /** \brief Compute a normalization factor for weighted data

      This function is used internally in \ref wvector_variance(size_t
      n, vec_t &data, const vec2_t &weights, double wmean) and \ref
      wvector_stddev(size_t n, vec_t &data, const vec2_t &weights, double
      wmean) .
  */
  template<class vec_t> double wvector_factor(size_t n, const vec_t &weights) {
    
    long double a=0.0;
    long double b=0.0;
    long double factor;
    for(size_t i=0;i<n;i++) {
      double wi=weights[i];
      if (wi>0.0) {
	a+=wi;
	b+=wi*wi;
      }
    }
    factor=a*a/(a*a-b);
    return factor;
  }

  /** \brief Compute a normalization factor for weighted data

      This function is used internally in \ref wvector_variance(size_t
      n, vec_t &data, const vec2_t &weights, double wmean) and \ref
      wvector_stddev(size_t n, vec_t &data, const vec2_t &weights, double
      wmean) .
  */
  template<class vec_t> double wvector_factor(const vec_t &weights) {
    return wvector_factor<vec_t>(weights.size(),weights);
  }
  
  /** \brief Compute the variance of a weighted vector with a mean
      known in advance

      This function computes
      \f[
      \left[ \sum_i w_i \left(x_i-\mu\right)^2 \right] 
      \left[ \sum_i w_i \right]^{-1}
      \f]

      This function produces the same results
      as <tt>gsl_stats_wvariance_with_fixed_mean()</tt>.

  */
  template<class vec_t, class vec2_t>
    double wvector_variance_fmean(size_t n, const vec_t &data,
				  const vec2_t &weights, double wmean) {
    long double wvariance=0.0;
    long double W=0.0;
    for(size_t i=0;i<n;i++) {
      double wi=weights[i];
      if (wi>0.0) {
	const long double delta=data[i]-wmean;
	W+=wi;
	wvariance+=(delta*delta-wvariance)*(wi/W);
      }
    }

    return wvariance;
  }

  /** \brief Compute the variance of a weighted vector with a mean
      known in advance

      This function computes
      \f[
      \left[ \sum_i w_i \left(x_i-\mu\right)^2 \right] 
      \left[ \sum_i w_i \right]^{-1}
      \f]

      This function produces the same results
      as <tt>gsl_stats_wvariance_with_fixed_mean()</tt>.

  */
  template<class vec_t, class vec2_t>
    double wvector_variance_fmean(const vec_t &data,
				  const vec2_t &weights, double wmean) {
    return wvector_variance_fmean(data.size(),data,weights,wmean);
  }

  /** \brief Compute the variance of a weighted vector with
      specified mean

      This function produces the same results
      as <tt>gsl_stats_wvariance_m()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_variance(size_t n, const vec_t &data,
			    const vec2_t &weights, double wmean) {

    const double variance=wvector_variance_fmean
      (n,data,weights,wmean);
    const double scale=wvector_factor(n,weights);
    const double wvar=scale*variance;
    return wvar;
  }

  /** \brief Compute the variance of a weighted vector with
      specified mean

      This function produces the same results
      as <tt>gsl_stats_wvariance_m()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_variance(const vec_t &data,
			    const vec2_t &weights, double wmean) {
    return wvector_variance<vec_t,vec2_t>(data.size(),data,weights,wmean);
  }

  /** \brief Compute the variance of a weighted vector where mean
      is computed automatically

      This function produces the same results
      as <tt>gsl_stats_wvariance()</tt>.
   */
  template<class vec_t, class vec2_t>
    double wvector_variance(size_t n, const vec_t &data,
			    const vec2_t &weights) {

    double wmean=wvector_mean(n,data,weights);
    return wvector_variance<vec_t,vec2_t>(n,data,weights,wmean);
  }

  /** \brief Compute the variance of a weighted vector where mean
      is computed automatically

      This function produces the same results
      as <tt>gsl_stats_wvariance()</tt>.
   */
  template<class vec_t, class vec2_t>
    double wvector_variance(const vec_t &data, const vec2_t &weights) {
    return wvector_variance(data.size(),data,weights);
  }

  /** \brief Compute the standard deviation of a weighted vector 
      with a mean known in advance

      This function produces the same results
      as <tt>gsl_stats_wsd_with_fixed_mean()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_stddev_fmean(size_t n, const vec_t &data,
				const vec2_t &weights, double wmean) {
    return sqrt(wvector_variance_fmean(n,data,weights,wmean));
  }

  /** \brief Compute the standard deviation of a weighted vector 
      with a mean known in advance

      This function produces the same results
      as <tt>gsl_stats_wsd_with_fixed_mean()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_stddev_fmean(const vec_t &data,
				const vec2_t &weights, double wmean) {
    return wvector_stddev_fmean<vec_t,vec2_t>
      (data.size(),data,weights,wmean);
  }

  /** \brief Compute the standard deviation of a weighted vector where mean
      is computed automatically

      This function produces the same results
      as <tt>gsl_stats_wsd()</tt>.
   */
  template<class vec_t, class vec2_t>
    double wvector_stddev(size_t n, const vec_t &data,
			  const vec2_t &weights) {
    double wmean=wvector_mean(n,data,weights);
    return sqrt(wvector_variance(n,data,weights,wmean));
  }

  /** \brief Compute the standard deviation of a weighted vector where mean
      is computed automatically

      This function produces the same results
      as <tt>gsl_stats_wsd()</tt>.
   */
  template<class vec_t, class vec2_t>
    double wvector_stddev(const vec_t &data, const vec2_t &weights) {
    return wvector_stddev(data.size(),data,weights);
  }

  /** \brief Compute the standard deviation of a weighted vector with
      specified mean

      This function produces the same results
      as <tt>gsl_stats_wsd_m()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_stddev(size_t n, const vec_t &data,
			  const vec2_t &weights, double wmean) {
    const double variance=wvector_variance_fmean
      (n,data,weights,wmean);
    const double scale=wvector_factor(n,weights);
    double wvar=scale*variance;
    return sqrt(wvar);
  }

  /** \brief Compute the standard deviation of a weighted vector with
      specified mean

      This function produces the same results
      as <tt>gsl_stats_wsd_m()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_stddev(const vec_t &data,
			  const vec2_t &weights, double wmean) {
    return wvector_stddev<vec_t,vec2_t>(data.size(),data,weights,wmean);
  }
  //@}

  /// \name Other weighted vector functions
  //@{
  /** \brief The weighted covariance of two vectors

      \note Experimental
  */
  template<class vec_t, class vec2_t, class vec3_t>
    double wvector_covariance(size_t n, const vec_t &data1,
			      const vec2_t &data2,
			      const vec3_t &weights) {
    double mean1=wvector_mean(n,data1,weights);
    double mean2=wvector_mean(n,data2,weights);
    double covar=0.0;
    double W=0.0;
    for(size_t i=0;i<n;i++) {
      double wi=weights[i];
      if (wi>0.0) {
	W+=wi;
	double delta1=(data1[i]-mean1);
	double delta2=(data2[i]-mean2);
	covar+=(wi/W)*(delta1*delta2-covar);
      }
    }
    double scale=wvector_factor(n,weights);
    return covar*scale;
  }

  /** \brief The weighted covariance of two vectors

      \note Experimental
  */
  template<class vec_t, class vec2_t, class vec3_t>
    double wvector_covariance(const vec_t &data1, const vec2_t &data2,
			      const vec3_t &weights) {
    return wvector_covariance<vec_t,vec2_t,vec3_t>
      (data1.size(),data1,data2,weights);
  }

  /** \brief Compute the weighted sum of squares of data about the 
      specified weighted mean

      This function produces the same results
      as <tt>gsl_stats_wtss_m()</tt>.
   */
  template<class vec_t, class vec2_t>
    double wvector_sumsq(size_t n, const vec_t &data,
			 const vec2_t &weights, double wmean) {
    long double wtss=0.0;
    for(size_t i=0;i<n;i++) {
      double wi=weights[i];
      if (wi>0.0) {
	const long double delta=data[i]-wmean;
	wtss+=wi*delta*delta;
      }
    }
    
    return wtss;
  }

  /** \brief Compute the weighted sum of squares of data about the 
      specified weighted mean

      This function produces the same results
      as <tt>gsl_stats_wtss_m()</tt>.
   */
  template<class vec_t, class vec2_t>
    double wvector_sumsq(const vec_t &data,
			 const vec2_t &weights, double wmean) {
    return wvector_sumsq<vec_t,vec2_t>(data.size(),data,weights,wmean);
  }

  /** \brief Compute the weighted sum of squares of data about the 
      weighted mean

      This function produces the same results
      as <tt>gsl_stats_wtss()</tt>.
   */
  template<class vec_t, class vec2_t>
    double wvector_sumsq(size_t n, const vec_t &data,
			 const vec2_t &weights) {
    
    double wmean=wvector_mean(n,data,weights);
    return wvector_sumsq(n,data,weights,wmean);
  }

  /** \brief Compute the weighted sum of squares of data about the 
      weighted mean

      This function produces the same results
      as <tt>gsl_stats_wtss()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_sumsq(const vec_t &data, const vec2_t &weights) {
    return wvector_sumsq<vec_t,vec2_t>(data.size(),data,weights);
  }

  /** \brief Compute the absolute deviation of data about a specified mean

      This function produces the same results
      as <tt>gsl_stats_wabsdev_m()</tt>.
   */
  template<class vec_t, class vec2_t> 
    double wvector_absdev(size_t n, const vec_t &data, const vec2_t &weights, 
			  double wmean) {
    long double wabsdev=0.0;
    long double W=0.0;
    for(size_t i=0;i<n;i++) {
      double wi=weights[i];
      if (wi>0.0) {
	const long double delta=fabs(data[i]-wmean);
	W+=wi;
	wabsdev+=(delta-wabsdev)*(wi/W);
      }
    }
    return wabsdev;
  }

  /** \brief Compute the absolute deviation of data about a specified mean

      This function produces the same results
      as <tt>gsl_stats_wabsdev_m()</tt>.
   */
  template<class vec_t, class vec2_t> 
    double wvector_absdev(const vec_t &data, const vec2_t &weights, 
			  double wmean) {
    return wvector_absdev<vec_t,vec2_t>(data.size(),data,weights,wmean);
  }

  /** \brief Compute the absolute deviation of data about a specified mean
      
      This function produces the same results
      as <tt>gsl_stats_wabsdev()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_absdev(size_t n, const vec_t &data,
			  const vec2_t &weights) {
    
    double wmean=wvector_mean(n,data,weights);
    return wvector_absdev(n,data,weights,wmean);
  }

  /** \brief Compute the absolute deviation of data about a specified mean
      
      This function produces the same results
      as <tt>gsl_stats_wabsdev()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_absdev(const vec_t &data, const vec2_t &weights) {
    return wvector_absdev<vec_t,vec2_t>(data.size(),data,weights);
  }

  /** \brief Compute the skewness of data with specified mean
      and standard deviation

      This function produces the same results
      as <tt>gsl_stats_wskew_m_sd()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_skew(size_t n, const vec_t &data, const vec2_t &weights,
			double wmean, double wsd) {
    long double wskew=0.0;
    long double W=0.0;
    for(size_t i=0;i<n;i++) {
      double wi=weights[i];
      if (wi>0.0) {
	const long double x=(data[i]-wmean)/wsd;
	W+=wi;
	wskew+=(x*x*x-wskew)*(wi/W);
      }
    }
    return wskew;
  }
  
  /** \brief Compute the skewness of data with specified mean
      and standard deviation

      This function produces the same results
      as <tt>gsl_stats_wskew_m_sd()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_skew(const vec_t &data, const vec2_t &weights,
			double wmean, double wsd) {
    return wvector_skew<vec_t,vec2_t>(data.size(),data,weights,wmean,wsd);
  }
  
  /** \brief Compute the skewness of data with specified mean
      and standard deviation
      
      This function produces the same results
      as <tt>gsl_stats_wskew()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_skew(size_t n, const vec_t &data, const vec2_t &weights) {
    double wmean=wvector_mean(n,data,weights);
    double wsd=wvector_stddev(n,data,weights,wmean);
    return wvector_skew(n,data,weights,wmean,wsd);
  }

  /** \brief Compute the skewness of data with specified mean
      and standard deviation
      
      This function produces the same results
      as <tt>gsl_stats_wskew()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_skew(const vec_t &data, const vec2_t &weights) {
    return wvector_skew<vec_t,vec2_t>(data.size(),data,weights);
  }

  /** \brief Compute the kurtosis of data with specified mean
      and standard deviation

      This function produces the same results
      as <tt>gsl_stats_wkurtosis_m_sd()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_kurtosis(size_t n, const vec_t &data, const vec2_t &weights,
			    double wmean, double wsd) {
    long double wavg=0.0;
    long double W=0.0;
    for(size_t i=0;i<n;i++) {
      double wi=weights[i];
      if (wi>0.0) {
	const long double x=(data[i]-wmean)/wsd;
	W+=wi;
	wavg+=(x*x*x*x-wavg)*(wi/W);
      }
    }
    return wavg-3.0;
  }

  /** \brief Compute the kurtosis of data with specified mean
      and standard deviation

      This function produces the same results
      as <tt>gsl_stats_wkurtosis_m_sd()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_kurtosis(const vec_t &data, const vec2_t &weights,
			    double wmean, double wsd) {
    return wvector_kurtosis<vec_t,vec2_t>
      (data.size(),data,weights,wmean,wsd);
  }

  /** \brief Compute the kurtosis of data with specified mean
      and standard deviation
      
      This function produces the same results
      as <tt>gsl_stats_wkurtosis()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_kurtosis(size_t n, const vec_t &data, 
			    const vec2_t &weights) {
    double wmean=wvector_mean(n,data,weights);
    double wsd=wvector_stddev(n,data,weights,wmean);
    return wvector_kurtosis(n,data,weights,wmean,wsd);
  }

  /** \brief Compute the kurtosis of data with specified mean
      and standard deviation
      
      This function produces the same results
      as <tt>gsl_stats_wkurtosis()</tt>.
  */
  template<class vec_t, class vec2_t>
    double wvector_kurtosis(const vec_t &data, const vec2_t &weights) {
    return wvector_kurtosis<vec_t,vec2_t>(data,weights);
  }
  //@}

  // This section has to appear after wvector_mean()
  /// \name Vector autocorrelation
  //@{
  /** \brief Lag-1 autocorrelation

      This function computes
      \f[
      \left[
      \sum_i \left(x_i - \mu\right) \left(x_{i-1} - \mu \right)
      \right] \left[ 
      \sum_i \left(x_i - \mu\right)^2 
      \right]^{-1}
      \f]

      This function produces the same results
      as <tt>gsl_stats_lag1_autocorrelation_m()</tt>.

      If \c n is less than 2, this function will call the error handler.
  */
  template<class vec_t>
    double vector_lag1_autocorr(size_t n, const vec_t &data, double mean) {
    
    if (n<2) {
      O2SCL_ERR2("Cannot compute lag1 with less than 2 elements",
		     " in vector_lag1_autocorr().",exc_einval);
    }

    long double q=0.0;
    long double v=(data[0]-mean)*(data[0]-mean);
    for(size_t i=1;i<n;i++) {
      long double delta0=data[i-1]-mean;
      long double delta1=data[i]-mean;
      q+=(delta0*delta1-q)/(i+1);
      v+=(delta1*delta1-v)/(i+1);
    }

    return q/v;
  }

  /** \brief Lag-1 autocorrelation

      This function computes
      \f[
      \left[
      \sum_i \left(x_i - \mu\right) \left(x_{i-1} - \mu \right)
      \right] \left[ 
      \sum_i \left(x_i - \mu\right)^2 
      \right]^{-1}
      \f]

      This function produces the same results
      as <tt>gsl_stats_lag1_autocorrelation_m()</tt>.

      If \c n is less than 2, this function will call the error handler.
  */
  template<class vec_t>
    double vector_lag1_autocorr(const vec_t &data, double mean) {
    return vector_lag1_autocorr(data.size(),data,mean);
  }

  /** \brief Lag-1 autocorrelation

      This function computes
      \f[
      \left[
      \sum_i \left(x_i - \mu\right) \left(x_{i-1} - \mu \right)
      \right] \left[ 
      \sum_i \left(x_i - \mu\right)^2 
      \right]^{-1}
      \f]
      
      This function produces the same results
      as <tt>gsl_stats_lag1_autocorrelation()</tt>.

      If \c n is less than 2, this function will call the error handler.
  */
  template<class vec_t> double vector_lag1_autocorr
    (size_t n, const vec_t &data) {
    double mean=vector_mean<vec_t>(n,data);
    return vector_lag1_autocorr(n,data,mean);
  }

  /** \brief Lag-1 autocorrelation

      This function computes
      \f[
      \left[
      \sum_i \left(x_i - \mu\right) \left(x_{i-1} - \mu \right)
      \right] \left[ 
      \sum_i \left(x_i - \mu\right)^2 
      \right]^{-1}
      \f]
      
      This function produces the same results
      as <tt>gsl_stats_lag1_autocorrelation()</tt>.

      If \c n is less than 2, this function will call the error handler.
  */
  template<class vec_t> double vector_lag1_autocorr(const vec_t &data) {
    return vector_lag1_autocorr(data.size(),data);
  }

  /** \brief Lag-k autocorrelation

      This function computes
      \f[
      \left[
      \sum_i \left(x_i - \mu\right) \left(x_{i-k} - \mu \right)
      \right] \left[ 
      \sum_i \left(x_i - \mu\right)^2 
      \right]^{-1}
      \f]

      If <tt>n<=k</tt>, this function will call the error handler.
  */
  template<class vec_t>
    double vector_lagk_autocorr(size_t n, const vec_t &data, size_t k,
				double mean) {
    
    if (n<=k) {
      O2SCL_ERR2("Not enough elements ",
		 "in vector_lagk_autocorr().",exc_einval);
    }

    long double q=0.0, v=0.0;
    for(size_t i=0;i<k;i++) {
      v+=(data[i]-mean)*(data[i]-mean)/(i+1);
    }
    for(size_t i=k;i<n;i++) {
      long double delta0=data[i-k]-mean;
      long double delta1=data[i]-mean;
      q+=(delta0*delta1-q)/(i+1);
      v+=(delta1*delta1-v)/(i+1);
    }
    return q/v;
  }

  /** \brief Lag-k autocorrelation

      This function computes
      \f[
      \left[
      \sum_i \left(x_i - \mu\right) \left(x_{i-k} - \mu \right)
      \right] \left[ 
      \sum_i \left(x_i - \mu\right)^2 
      \right]^{-1}
      \f]

      If <tt>n<=k</tt>, this function will call the error handler.
  */
  template<class vec_t>
    double vector_lagk_autocorr(const vec_t &data, size_t k,
				double mean) {
    return vector_lagk_autocorr(data.size(),k,mean);
  }

  /** \brief Lag-k autocorrelation

      This function computes
      \f[
      \left[
      \sum_i \left(x_i - \mu\right) \left(x_{i-k} - \mu \right)
      \right] \left[ 
      \sum_i \left(x_i - \mu\right)^2 
      \right]^{-1}
      \f]

      If <tt>n<=k</tt>, this function will call the error handler.
  */
  template<class vec_t> double vector_lagk_autocorr
    (size_t n, const vec_t &data, size_t k) {
    double mean=vector_mean<vec_t>(n,data);
    return vector_lagk_autocorr(n,data,k,mean);
  }

  /** \brief Lag-k autocorrelation

      This function computes
      \f[
      \left[
      \sum_i \left(x_i - \mu\right) \left(x_{i-k} - \mu \right)
      \right] \left[ 
      \sum_i \left(x_i - \mu\right)^2 
      \right]^{-1}
      \f]

      If <tt>n<=k</tt>, this function will call the error handler.
  */
  template<class vec_t> double vector_lagk_autocorr
    (const vec_t &data, size_t k) {
    return vector_lagk_autocorr(data.size(),data,k);
  }

  /** \brief Construct an autocorrelation vector

      This constructs a vector \c ac_vec for which the kth entry
      stores the lag-k autocorrelation. This function chooses \f$
      k_{\mathrm{max}} =n/2 \f$ where \f$ n \f$ is the length of the
      \c data vector. The vector \c ac_vec is resized to accomodate
      exactly \f$ k_{\mathrm{max}} \f$ values, from 0 to 
      \f$ k_{\mathrm{max}}-1 \f$.
  */
  template<class vec_t, class resize_vec_t> void vector_autocorr_vector
    (const vec_t &data, resize_vec_t &ac_vec) {
    size_t kmax=data.size()/2;
    double mean=vector_mean(data);
    ac_vec.resize(kmax);
    ac_vec[0]=1.0;
    for(size_t k=1;k<kmax;k++) {
      ac_vec[k]=vector_lagk_autocorr(data.size(),data,k,mean);
    }
    return;
  }

  /** \brief Use the Goodman method to compute the
      autocorrelation length

      Representing the lag-k correlation coefficient by
      \f$ \hat{C}(k) \f$, Goodman defines
      \f[
      \hat{\tau}(M) = 1 + 2 \sum_{s=1}^{M} \frac{\hat{C}(k)}{\hat{C}(0)}
      \f]
      Then the autocorrelation length is the value of 
      \f$ \hat{\tau}(M) \f$ for which 
      \f[
      5 \hat{\tau}(M)/M \leq 1
      \f]

      This function computes the value of \f$ 5 \hat{\tau}(M)/M \f$
      and stores it in the <tt>five_tau_over_M</tt> vector and then
      returns the first value of \f$ M \f$ for which the vector is
      less than or equal to 1.0. If this function returns 0, then \f$
      5 \hat{\tau}(M)/M \f$ is greater than 1.0 for all \f$ M \f$, and
      this can be a sign that the autocorrelation length is too long
      to accurately resolve.

      On completion, the vector \c five_tau_over_m will have
      one less element than the vector \c ac_vec .
  */
  template<class vec_t, class resize_vec_t> size_t vector_autocorr_tau
    (const vec_t &ac_vec, resize_vec_t &five_tau_over_M) {
    five_tau_over_M.resize(0);
    size_t len=0;
    bool len_set=false;
    for (size_t M=1;M<ac_vec.size();M++) {
      double sum=0.0;
      for(size_t s=1;s<=M;s++) {
	sum+=ac_vec[s];
      }
      double val=(1.0+2.0*sum)/((double)M)*5.0;
      if (len_set==false && val<=1.0) {
	len=M;
	len_set=true;
      }
      five_tau_over_M.push_back(val);
    }
    return len;
  }
  
  /** \brief Lag-k autocorrelation for the first
      \c n elements with a vector multiplier given 
      the mean
   */
  template<class vec_t, class vec2_t>
    double vector_lagk_autocorr_mult(size_t n, const vec_t &data,
				     const vec2_t &mult, size_t k,
				     double mean) {
    
    size_t n2=0;
    for(size_t i=0;i<n;i++) {
      size_t m=((size_t)(mult[i]*(1.0+1.0e-10)));
      if (m==0) {
	O2SCL_ERR2("Mult vector is zero ",
		   "in vector_lagk_autocorr_mult().",exc_einval);
      }
      n2+=m;
    }
    
    if (n2<=k) {
      O2SCL_ERR2("Not enough elements ",
		 "in vector_lagk_autocorr_mult().",exc_einval);
    }

    long double q=0.0, v=0.0;
    size_t im=0, ix=0, im2=0, ix2=0;
    for(size_t i=0;i<k;i++) {
      v+=(data[ix]-mean)*(data[ix]-mean)/(i+1);
      im++;
      if (im>=((size_t)(mult[ix]*(1.0+1.0e-10)))) {
	im=0;
	ix++;
      }
    }
    for(size_t i=k;i<n2;i++) {
      long double delta0=data[ix2]-mean;
      long double delta1=data[ix]-mean;
      q+=(delta0*delta1-q)/(i+1);
      v+=(delta1*delta1-v)/(i+1);
      im++;
      if (im>=((size_t)(mult[ix]*(1.0+1.0e-10)))) {
	im=0;
	ix++;
      }
      im2++;
      if (im2>=((size_t)(mult[ix2]*(1.0+1.0e-10)))) {
	im2=0;
	ix2++;
      }
    }
    return q/v;
  }

  /** \brief Lag-k autocorrelation for the first
      \c n elements with a vector multiplier
   */
  template<class vec_t, class vec2_t>
    double vector_lagk_autocorr_mult(size_t n, const vec_t &data,
				     const vec2_t &mult, size_t k) {
    double mean=wvector_mean(n,mult,data);
    return vector_lagk_autocorr_mult(n,data,mult,k,mean);
  }

  /** \brief Lag-k autocorrelation with a vector multiplier
      given the mean
   */
  template<class vec_t, class vec2_t>
    double vector_lagk_autocorr_mult(const vec_t &data,
				     const vec2_t &mult, size_t k,
				     double mean) {
    return vector_lagk_autocorr_mult(data.size(),data,mult,k,mean);
  }
  
  /** \brief Lag-k autocorrelation with a vector multiplier
   */
  template<class vec_t, class vec2_t>
    double vector_lagk_autocorr_mult(const vec_t &data,
				     const vec2_t &mult, size_t k) {
    return vector_lagk_autocorr_mult(data.size(),data,mult,k);
  }
  
  /** \brief Construct an autocorrelation vector using a multiplier
      using the first \c n2 elements of vectors \c data and \c mult
   */
  template<class vec_t, class vec2_t, class resize_vec_t>
    void vector_autocorr_vector_mult
    (size_t n2, const vec_t &data, const vec2_t &mult, resize_vec_t &ac_vec) {

    size_t n=0;
    for(size_t i=0;i<n2;i++) {
      n+=((size_t)(mult[i]*(1.0+1.0e-10)));
    }
    
    size_t kmax=n/2;
    double mean=wvector_mean(data,mult);
    ac_vec.resize(kmax);
    ac_vec[0]=1.0;
    for(size_t k=1;k<kmax;k++) {
      ac_vec[k]=vector_lagk_autocorr_mult(n2,data,mult,k,mean);
    }
    return;
  }

  /** \brief Construct an autocorrelation vector using a multiplier
   */
  template<class vec_t, class vec2_t, class resize_vec_t>
    void vector_autocorr_vector_mult
    (const vec_t &data, const vec2_t &mult, resize_vec_t &ac_vec) {

    vector_autocorr_vector_mult(data.size(),data,mult,ac_vec);
    
    return;
  }
  //@}

  /** \brief Take a vector of data and convert it to a vector
      of bin edges automatically adjusting for increasing or
      decreasing and linear or logarithmic spacing
   */
  template<class vec_t, class vec2_t>
    void vector_to_bins(const vec_t &v_grid, vec2_t &v_bins,
			int verbose=1) {

    size_t n=v_grid.size();
    v_bins.resize(n+1);

    // Vector of differences
    std::vector<double> diffs;

    // Compute quality factor for linear bins (smaller number is
    // better)
    vector_diffs(v_grid,diffs);
    double mean=vector_mean(diffs);
    double std=vector_stddev(diffs);
    double qual_lin=std/mean;

    // Compute quality factor for logarithmic bins (smaller number is
    // better)
    std::vector<double> logs;
    for(size_t i=0;i<n;i++) logs[i]=log(v_grid[i]);
    vector_diffs(logs,diffs);
    mean=vector_mean(diffs);
    std=vector_stddev(diffs);
    double qual_log=std/mean;

    if (qual_log<qual_lin) {
      if (verbose>0) {
	std::cout << "Auto-detected log mode in vector_to_bins()."
		  << std::endl;
      }
      if (logs[1]>logs[0]) {
	// Increasing, log
	v_bins[0]=exp(logs[0]-(logs[1]-logs[0])/2.0);
	v_bins[n]=exp(logs[n-1]+(logs[n-1]-logs[n-2])/2.0);
      } else {
	// Decreasing, log
	v_bins[0]=exp(logs[0]+(logs[0]-logs[1])/2.0);
	v_bins[n]=exp(logs[n-1]-(logs[n-2]-logs[n-1])/2.0);
      }
      for(size_t i=1;i<n-1;i++) {
	v_bins[i]=exp((logs[i-1]+logs[i])/2.0);
      }
    } else {
      if (v_grid[1]>v_grid[0]) {
	// Increasing, linear
	v_bins[0]=v_grid[0]-(v_grid[1]-v_grid[0])/2.0;
	v_bins[n]=v_grid[n-1]+(v_grid[n-1]-v_grid[n-2])/2.0;
      } else {
	// Decreasing, linear
	v_bins[0]=v_grid[0]+(v_grid[0]-v_grid[1])/2.0;
	v_bins[n]=v_grid[n-1]-(v_grid[n-2]-v_grid[n-1])/2.0;
      }
      for(size_t i=1;i<n-1;i++) {
	v_bins[i]=(v_grid[i-1]+v_grid[i])/2.0;
      }
    }
    
    return;
  }

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
