/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/vec_stats.h>
#include <gsl/gsl_statistics.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  cout.setf(ios::scientific);

  const size_t N=10;
  double mean1=3.7;

  double x[N]={1,2,3,4,5,6,4,3,4,5};
  double x2[N]={3,1,4,1,5,9,2,6,5,3};
  double x3[12]={3,1,4,1,5,9,2,6,5,3,5,8};
  double w[N]={0.1,0.2,0.1,0.2,0.1,
	       0.1,0.2,0.1,0.2,0.1};
  
  t.test_rel(vector_max_value<double [N],double>(N,x),
	     gsl_stats_max(x,1,N),1.0e-8,"max");
  t.test_rel(vector_min_value<double [N],double>(N,x),
	     gsl_stats_min(x,1,N),1.0e-8,"min");
  t.test_rel(vector_sum<double [N],double>(N,x),
	     37.0,1.0e-8,"sum");
  t.test_rel(vector_mean(N,x),
	     gsl_stats_mean(x,1,N),1.0e-8,"mean");

  t.test_rel(vector_variance<double [N]>(N,x),
	     gsl_stats_variance(x,1,N),1.0e-8,"variance 1");
  t.test_rel(vector_variance<double [N]>(N,x,mean1),
	     gsl_stats_variance_m(x,1,N,mean1),1.0e-8,"variance 2");
  t.test_rel(vector_variance_fmean<double [N]>(N,x,mean1),
	     gsl_stats_variance_with_fixed_mean(x,1,N,mean1),1.0e-8,
	     "variance 3");

  double sdtmp=vector_stddev<double [N]>(N,x);
  t.test_rel(vector_stddev<double [N]>(N,x),
	     gsl_stats_sd(x,1,N),1.0e-8,"sd 1");
  t.test_rel(vector_stddev<double [N]>(N,x,mean1),
	     gsl_stats_sd_m(x,1,N,mean1),1.0e-8,"sd 2");
  t.test_rel(vector_stddev_fmean<double [N]>(N,x,mean1),
	     gsl_stats_sd_with_fixed_mean(x,1,N,mean1),1.0e-8,"sd 3");

  t.test_rel(vector_absdev<double [N]>(N,x),
	     gsl_stats_absdev(x,1,N),1.0e-8,"abs 1");
  t.test_rel(vector_absdev<double [N]>(N,x,mean1),
	     gsl_stats_absdev_m(x,1,N,mean1),1.0e-8,"abs 2");
    
  t.test_rel(vector_skew<double [N]>(N,x),
	     gsl_stats_skew(x,1,N),1.0e-8,"skew 1");
  t.test_rel(vector_skew<double [N]>(N,x,mean1,sdtmp),
	     gsl_stats_skew_m_sd(x,1,N,mean1,sdtmp),1.0e-8,"skew 2");
  
  t.test_rel(vector_kurtosis<double [N]>(N,x),
	     gsl_stats_kurtosis(x,1,N),1.0e-8,"kurt 1");
  t.test_rel(vector_kurtosis<double [N]>(N,x,mean1,sdtmp),
	     gsl_stats_kurtosis_m_sd(x,1,N,mean1,sdtmp),1.0e-8,"kurt 2");
  
  t.test_rel(vector_lag1_autocorr<double [N]>(N,x),
	     gsl_stats_lag1_autocorrelation(x,1,N),1.0e-8,"lag1 1");
  t.test_rel(vector_lag1_autocorr<double [N]>(N,x,mean1),
	     gsl_stats_lag1_autocorrelation_m(x,1,N,mean1),1.0e-8,"lag1 2");
  
  t.test_rel(vector_covariance<double [N]>(N,x,x2),
	     gsl_stats_covariance(x,1,x2,1,N),1.0e-8,"covariance 1");
  t.test_rel(vector_covariance<double [N]>
	     (N,x,x2,vector_mean<double [N]>(N,x),
	      vector_mean<double [N]>(N,x2)),
	     gsl_stats_covariance_m
	     (x,1,x2,1,N,vector_mean<double [N]>(N,x),
	      vector_mean<double [N]>(N,x2)),
	     1.0e-8,"covariance 2");

  t.test_rel(vector_correlation<double [N]>(N,x,x2),
	     gsl_stats_correlation(x,1,x2,1,N),1.0e-8,"correlation 1");

  t.test_rel(vector_pvariance<double [N],double [12]>(N,x,12,x3),
	     gsl_stats_pvariance(x,1,N,x3,1,12),1.0e-8,"pvariance 1");
  
  vector_sort<double [N],double>(N,x);
  t.test_rel(vector_median_sorted<double [N]>(N,x),
	     gsl_stats_median_from_sorted_data(x,1,N),1.0e-8,
	     "median from sorted");

  t.test_rel(vector_quantile_sorted<double [N]>(N,x,0.25),
	     gsl_stats_quantile_from_sorted_data(x,1,N,0.25),1.0e-8,
	     "quantile from sorted");
  t.test_rel(vector_quantile_sorted<double [N]>(N,x,0.9),
	     gsl_stats_quantile_from_sorted_data(x,1,N,0.9),1.0e-8,
	     "quantile from sorted");

  double wmean=wvector_mean(N,x,w);

  t.test_rel(wmean,gsl_stats_wmean(w,1,x,1,N),1.0e-8,"wgtd. mean");

  t.test_rel(wvector_variance<double [N]>(N,x,w),
	     gsl_stats_wvariance(w,1,x,1,N),1.0e-8,"wvariance 1");
  t.test_rel(wvector_variance<double [N]>(N,x,w,wmean),
	     gsl_stats_wvariance_m(w,1,x,1,N,wmean),1.0e-8,"wvariance 2");
  t.test_rel(wvector_variance_fmean<double [N],double [N]>(N,x,w,wmean),
	     gsl_stats_wvariance_with_fixed_mean(w,1,x,1,N,wmean),1.0e-8,
	     "wvariance 3");

  t.test_rel(wvector_stddev<double [N]>(N,x,w),
	     gsl_stats_wsd(w,1,x,1,N),1.0e-8,"wsd 1");
  t.test_rel(wvector_stddev<double [N]>(N,x,w,wmean),
	     gsl_stats_wsd_m(w,1,x,1,N,wmean),1.0e-8,"wsd 2");
  t.test_rel(wvector_stddev_fmean<double [N]>(N,x,w,wmean),
	     gsl_stats_wsd_with_fixed_mean(w,1,x,1,N,wmean),1.0e-8,"wsd 3");

  double wsdtmp=wvector_stddev(N,x,w);

  t.test_rel(wvector_sumsq<double [N]>(N,x,w),
	     gsl_stats_wtss(w,1,x,1,N),1.0e-8,"wtss 1");
  t.test_rel(wvector_sumsq<double [N]>(N,x,w,wmean),
	     gsl_stats_wtss_m(w,1,x,1,N,wmean),1.0e-8,"wtss 2");

  t.test_rel(wvector_absdev<double [N]>(N,x,w),
	     gsl_stats_wabsdev(w,1,x,1,N),1.0e-8,"wabs 1");
  t.test_rel(wvector_absdev<double [N]>(N,x,w,wmean),
	     gsl_stats_wabsdev_m(w,1,x,1,N,wmean),1.0e-8,"wabs 2");

  t.test_rel(wvector_skew<double [N]>(N,x,w),
	     gsl_stats_wskew(w,1,x,1,N),1.0e-8,"wskew 1");
  t.test_rel(wvector_skew<double [N]>(N,x,w,wmean,wsdtmp),
	     gsl_stats_wskew_m_sd(w,1,x,1,N,wmean,wsdtmp),1.0e-8,"wskew 2");
  
  t.test_rel(wvector_kurtosis<double [N]>(N,x,w),
	     gsl_stats_wkurtosis(w,1,x,1,N),1.0e-8,"wkurtosis 1");
  t.test_rel(wvector_kurtosis<double [N]>(N,x,w,wmean,wsdtmp),
	     gsl_stats_wkurtosis_m_sd(w,1,x,1,N,wmean,wsdtmp),1.0e-8,
	     "wkurtosis 2");

  t.report();
  
  return 0;
}

