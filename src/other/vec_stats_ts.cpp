/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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

#ifdef O2SCL_SET_FFTW
#include <fftw3.h>
#endif

#include <gsl/gsl_statistics.h>
#include <o2scl/test_mgr.h>
#include <o2scl/vec_stats.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/invert.h>

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
  t.test_rel(vector_mean<double [N],double>(N,x),
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
	     (N,x,x2,vector_mean<double [N],double>(N,x),
	      vector_mean<double [N],double>(N,x2)),
	     gsl_stats_covariance_m
	     (x,1,x2,1,N,vector_mean<double [N],double>(N,x),
	      vector_mean<double [N],double>(N,x2)),
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
    rng<> r;
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
    o2scl::vector_autocorr_vector(act.size(),act,ac,1000);
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
    rng<> r;
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
    t.test_abs(tau,25.0,50.0,"acor 2");
    std::vector<double> ac, ac2, ftom, ftom2;
    o2scl::vector_autocorr_vector(act.size(),act,ac);
    size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
    t.test_abs(((double)ac_len),25.0,6.0,"vector_autocorr_tau 2");
    cout << "ac_len: " << ac_len << endl;
    cout << endl;

#ifdef O2SCL_SET_FFTW

    mean=vector_mean(act);
    sig=vector_stddev(act);
    vector<double> fft_out;
    vector_autocorr_vector_fftw(act,fft_out,mean,sig);
    // At this point, the vector ac has half the size of the fft_out
    // vector
    t.test_abs_vec(ac.size(),fft_out,ac,4.0e-2,"autocorr comparison");
    
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
  o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> mid;
  double kl1=kl_div_gaussian(1,mean_prior,mean_post,
                             covar_prior,covar_post,mid);
  double kl2=kl_div_gaussian(mean_prior(0),mean_post(0),
                             covar_prior(0,0),covar_post(0,0));
  t.test_rel(kl1,kl2,1.0e-12,"KL div");

#ifdef O2SCL_SET_FFTW

  {
    vector<double> sines, sines3;
    vector<complex<double>> sines2, fft, fft2, sines4;

    cout << "Data: " << endl;
    // We choose 26 here because FFTW is better at arrays with
    // certain sizes
    cout.setf(ios::showpos);
    for(double xq=0.0;xq<o2scl_const::pi*4.001;
        xq+=o2scl_const::pi*4.0/25.0) {
      sines.push_back(sin(2.0*xq)+sin(4.0*xq));
      complex<double> xc=sin(2.0*xq)+sin(4.0*xq);
      cout << xc.real() << endl;
      sines2.push_back(xc);
    }
    cout.unsetf(ios::showpos);
    cout << endl;

    // Test real forward and backward transformations
    vector_forward_fft(sines,fft);
    cout.setf(ios::showpos);
    for(size_t j=0;j<fft.size();j++) {
      cout.width(2);
      cout << j << " " << sines[j] << " ";
      cout << fft[j].real() << " " << fft[j].imag() << endl;
    }
    cout << endl;
    vector_backward_fft(fft,sines3);
    
    for(size_t j=0;j<sines.size();j++) {
      if (j<fft.size()) {
        cout.width(2);
        cout << j << " " << sines[j] << " " << fft[j].real() << " "
             << fft[j].imag() << " "
             << sines3[j] << " " << sines3[j]/sines[j] << endl;
      } else {
        cout.width(2);
        cout << j << " " << sines[j] << " " << "             " << " "
             << "             " << " " 
             << sines3[j] << " " << sines3[j]/sines[j] << endl;
      }
      if (j>0 && j<25) {
        t.test_rel(sines3[j]/sines[j],26.0,1.0e-6,"forward and backward");
      } else if (j==25) {
        t.test_rel(sines3[j]/sines[j],26.0,1.0e-2,"forward and backward");
      }
    }

    cout.unsetf(ios::showpos);
    cout << endl;
    
    // Compare real and complex FFTs
    vector_forward_complex_fft(sines2,fft2);
    cout << "Complex: " << sines2.size() << " " << fft2.size() << endl;
    cout.setf(ios::showpos);
    for(size_t j=0;j<fft2.size();j++) {
      cout.width(2);
      cout << j << " ";
      cout << sines2[j].real() << " " << sines2[j].imag() << " "
           << fft2[j].real() << " " << fft2[j].imag() << endl;
      if (j<fft.size()) {
        if (abs(fft2[j].real())>1.0e-4) {
          t.test_rel(fft[j].real(),fft2[j].real(),1.0e-10,"fft real part");
        }
        if (abs(fft2[j].imag())>1.0e-4) {
          t.test_rel(fft[j].imag(),fft2[j].imag(),1.0e-10,"fft imag part");
        }
      }
    }
    cout.unsetf(ios::showpos);
    cout << endl;

    // Compare complex forward and backward FFTs
    vector_backward_complex_fft(fft2,sines4);
    cout.setf(ios::showpos);
    for(size_t j=0;j<sines2.size();j++) {
      cout.width(2);
      cout << j << " ";
      cout << sines2[j].real() << " " << sines2[j].imag() << " "
           << sines4[j].real() << " " << sines4[j].imag() << " "
           << sines4[j].real()/sines2[j].real() << endl;
      if (j!=0 && j!=sines2.size()-1) {
        t.test_rel(sines4[j].real()/sines2[j].real(),26.0,
                   1.0e-6,"complex forward backward");
      }
    }
    cout.unsetf(ios::showpos);
    cout << endl;
  }
  
  if (true) {
    
    // Prepare data for 2D FFTs with a rectangular matrix
    size_t n_rows=11;
    size_t n_cols=14;
    vector<double> sines(n_rows*n_cols), sines3;
    vector<complex<double>> fft, sines2, fft2, sines4;
    vector_view_matrix<vector<double>,double>
      vvm_sines(sines,n_rows*n_cols,n_rows);
  
    cout << "Matrix data: " << endl;
    cout.setf(ios::showpos);
    for(size_t ii=0;ii<n_rows*n_cols;ii++) {
      size_t i, j;
      vvm_sines.matrix_indices(ii,i,j);

      double xq=sin(i*o2scl_const::pi/2)+sin(j*o2scl_const::pi/2)+
        sin(i*o2scl_const::pi/4)+sin((i+j)*o2scl_const::pi/8);
      vvm_sines(i,j)=xq;

      complex<double> xc=xq;
      sines2.push_back(xc);
      cout.width(2);
      cout << i << " ";
      cout.width(2);
      cout << j << " " << xq << endl;
    }
    cout.unsetf(ios::showpos);
    cout << endl;

    // Compare the forward and backward FFTs
    matrix_forward_fft(n_rows,n_cols,sines,fft);
    cout << "Matrix real: " << fft.size() << endl;
    cout.setf(ios::showpos);
    for(size_t j=0;j<fft.size();j++) {
      cout.width(2);
      cout << j << " ";
      cout << fft[j].real() << " " << fft[j].imag() << endl;
    }
    cout.unsetf(ios::showpos);
    cout << endl;
    cout << "1." << endl;
    matrix_backward_fft_copy(n_rows,(n_cols)/2+1,fft,sines3);
    cout << "2." << endl;
    cout.setf(ios::showpos);
    for(size_t j=0;j<fft.size();j++) {
      cout.width(2);
      cout << j << " ";
      cout << sines[j] << " " << sines3[j] << " "
           << sines3[j]/sines[j] << endl;
      if (fabs(sines[j])>1.0e-10) {
        t.test_rel(sines3[j]/sines[j],(double)(11*14),1.0e-6,"fw and bw");
      }
    }
    cout.unsetf(ios::showpos);
    cout << endl;

    // Compare real and complex FFTs
    cout << "Matrix complex: " << endl;
    matrix_forward_complex_fft(n_rows,n_cols,sines2,fft2);
    cout.setf(ios::showpos);
    for(size_t j=0;j<fft2.size();j++) {
      cout.width(2);
      cout << j << " ";
      cout << fft2[j].real() << " " << fft2[j].imag() << endl;
    }
    cout << endl;
    vector_view_matrix<vector<complex<double>>,complex<double>>
      vvm_fft(fft,fft.size(),n_rows);
    vector_view_matrix<vector<complex<double>>,complex<double>>
      vvm_fft2(fft2,fft2.size(),n_rows);
    for(size_t i=0;i<n_rows;i++) {
      for(size_t j=0;j<n_cols;j++) {
        cout.width(2);
        cout << i << " ";
        cout.width(2);
        cout << j << " ";
        if (j<vvm_fft.size2()) {
          cout << fft[j].real() << " " << fft[j].imag() << " ";
        } else {
          cout << "+x.xxxxxxe+xx +x.xxxxxxe+xx ";
        }
        cout << fft2[j].real() << " " << fft2[j].imag() << endl;
        if (j<vvm_fft.size2()) {
          if (abs(fft2[j].real())>1.0e-4) {
            t.test_rel(fft[j].real(),fft2[j].real(),1.0e-10,
                       "fft real part");
          }
          if (abs(fft2[j].imag())>1.0e-4) {
            t.test_rel(fft[j].imag(),fft2[j].imag(),1.0e-10,
                       "fft imag part");
          }
        }
      }
    }
    cout.unsetf(ios::showpos);
    cout << endl;

    // Compare complex forward and backward FFTs
    cout << "Matrix complex: " << endl;
    matrix_backward_complex_fft(n_rows,n_cols,fft2,sines4);
    cout.setf(ios::showpos);
    for(size_t j=0;j<n_rows*n_cols;j++) {
      cout.width(2);
      cout << j << " ";
      cout << sines2[j].real() << " " << sines2[j].imag() << " ";
      cout << sines4[j].real() << " " << sines4[j].imag() << " "
           << sines4[j].real()/sines2[j].real() << " "
           << sines4[j].real()/sines2[j].real() << endl;
    }

  }
    
#endif    

  if (true) {
    vector<double> x1, x4;
    rng<> r;
    r.clock_seed();
    for(size_t i=0;i<100;i++) {
      x1.push_back(2.0+r.random());
      x4.push_back(1.0/((double)(i+1))+2.0+r.random());
    }
    vector<vector<double>> xall;
    xall.push_back(x1);
    xall.push_back(x4);
    double gr=mult_vector_gelman_rubin<vector<double>,double>(xall);
    t.test_rel(gr,1.0,1.0,"gelman rubin");
  }
  
  t.report();
  
  return 0;
}

