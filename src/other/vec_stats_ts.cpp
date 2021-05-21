/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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

#ifdef O2SCL_FFTW
#include <fftw3.h>
#endif

#include <gsl/gsl_statistics.h>
#include <o2scl/test_mgr.h>
#include <o2scl/vec_stats.h>
#include <o2scl/prob_dens_func.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);
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

  if (true) {
    cout << endl;
    cout << "------------------------------------------------------------"
	 << endl;
    cout << "Testing vector_autocorr_vector_mult(): " << endl;
    cout << endl;
    
    /* 
       This is the test case in Goodman's original acor code
       which should report a correlation length of about 19.
     */
    vector<double> act;
    double aca=0.9;
    double acx=0.0;
    rng_gsl r;
    r.clock_seed();
    for(size_t i=0;i<4000000;i++) {
      acx=aca*acx+r.random();
      act.push_back(acx);
    }
    double mean, sig, tau;
    vector_acor(act.size(),act,mean,sig,tau,1);
    t.test_abs(tau,19.0,2.0,"acor");
    cout << "Results from acor: " << mean << " " << sig << " "
         << tau << endl;
    std::vector<double> ac, ac2, ftom, ftom2;
    o2scl::vector_autocorr_vector(act,ac,1000);
    size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
    cout << "ac_len: " << ac_len << endl;
    t.test_abs(((double)ac_len),19.0,2.0,"vector_autocorr_tau");
    cout << endl;
  }
  
  if (true) {
    /* 
       Test a simple data set with a known covariance length of 25
     */
    vector<double> act;
    rng_gsl r;
    r.clock_seed();
    double x0=r.random();
    static const size_t NN=100000;
    for(size_t i=0;i<NN;i++) {
      if (i%25==24) x0=r.random();
      act.push_back(x0);
    }
    double mean, sig, tau;
    vector_acor(act.size(),act,mean,sig,tau,0);
    cout << "Results from acor: " << mean << " " << sig << " "
         << tau << endl;
    t.test_abs(tau,25.0,5.0,"acor 2");
    std::vector<double> ac, ac2, ftom, ftom2;
    o2scl::vector_autocorr_vector(act,ac);
    size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
    t.test_abs(((double)ac_len),25.0,5.0,"vector_autocorr_tau 2");
    cout << "ac_len: " << ac_len << endl;
    cout << endl;

#ifdef O2SCL_FFTW

    /*
      From https://github.com/kaityo256/fftw_sample

      MIT License
      
      Copyright (c) 2012 H. Watanabe 
      
      Permission is hereby granted, free of charge, to any person
      obtaining a copy of this software and associated documentation
      files (the "Software"), to deal in the Software without
      restriction, including without limitation the rights to use,
      copy, modify, merge, publish, distribute, sublicense, and/or
      sell copies of the Software, and to permit persons to whom the
      Software is furnished to do so, subject to the following
      conditions:
      
      The above copyright notice and this permission notice shall be
      included in all copies or substantial portions of the Software.
      
      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
      OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
      NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
      HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
      WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
      FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
      OTHER DEALINGS IN THE SOFTWARE.
     */
    
    fftw_complex *in;
    fftw_complex *out;
    in=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*NN);
    out=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*NN);
    
    mean=vector_mean(act);
    sig=vector_stddev(act);
    
    for(size_t i=0;i<NN;i++){
      in[i][0]=(act[i]-mean)/sig;
      in[i][1]=0.0;
    }
    
    fftw_plan plan=fftw_plan_dft_1d(NN,in,out,
                                      FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    
    for(size_t i=0;i<NN;i++) {
      double re=out[i][0];
      double im=out[i][1];
      in[i][0]=(re*re+im*im)/((double)NN);
      in[i][1]=0.0;
    }
    fftw_plan plan2=fftw_plan_dft_1d(NN,in,out,
                                     FFTW_BACKWARD,FFTW_ESTIMATE);
    
    fftw_execute(plan2);
    
    vector<double> fft_out;
    for(size_t i=0;i<NN;i++){
      fft_out.push_back(out[i][0]/((double)NN));
    }

    for(size_t i=0;i<100;i++) {
      cout << i << " " << fft_out[i] << " " << ac[i] << endl;
    }
    
    fftw_free(in);
    fftw_free(out);
    
#endif

  }

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  
  ubmatrix covar_prior(1,1);
  covar_prior(0,0)=2.0;
  ubmatrix covar_post(1,1);
  covar_post(0,0)=3.0;
  ubvector mean_prior(1);
  mean_prior(0)=5.0;
  ubvector mean_post(1);
  mean_post(0)=7.0;
  double kl1=kl_div_gaussian(1,mean_prior,mean_post,
                             covar_prior,covar_post);
  double kl2=kl_div_gaussian(mean_prior(0),mean_post(0),
                             covar_prior(0,0),covar_post(0,0));
  t.test_rel(kl1,kl2,1.0e-12,"KL div");
  
  t.report();
  
  return 0;
}

