/*
  -------------------------------------------------------------------
  
  Copyright (C) 2021-2023, Andrew W. Steiner and Jesse Farr
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/vec_stats.h>

using namespace std;
using namespace o2scl;

double o2scl::kl_div_gaussian(double mean_prior, double mean_post,
                              double covar_prior, double covar_post) {
  
  double covar_prior_inv=1.0/covar_prior;
  double prod1=covar_prior_inv*covar_post;
  double diff=mean_prior-mean_post;
  double prod2=diff*covar_prior_inv*diff;
  double div=0.5*(prod1+prod2-1.0+log(covar_prior/covar_post));
  
  return div;
}

void o2scl::vector_forward_fft(const std::vector<double> &data,
                               std::vector<std::complex<double>> &fft) {
  
#ifdef O2SCL_FFTW
  
  fft.resize(data.size()/2+1);
  
  // First, note that FFTW_ESTIMATE means that FFTW is estimating an
  // efficient plan, not that FFTW is estimating the FFT. The result
  // is exact. What it does mean is that the input data is not
  // modified, so we allow the input vector to be const even though
  // we have to cast that const-ness away.
  double *non_const=(double *)(&(data[0]));
  fftw_complex *fft2=reinterpret_cast<fftw_complex *>(&(fft[0]));
  fftw_plan plan=fftw_plan_dft_r2c_1d(data.size(),non_const,fft2,
                                      FFTW_ESTIMATE);
  fftw_execute(plan);
  
#else
  
  O2SCL_ERR("FFTW support not included in this O2scl installation.",
            o2scl::exc_eunsup);
  
#endif
  
  return;
}

void o2scl::vector_forward_complex_fft
(const std::vector<std::complex<double>> &data,
 std::vector<std::complex<double>> &fft) {
  
#ifdef O2SCL_FFTW
  
  fft.resize(data.size());
  
  // First, note that FFTW_ESTIMATE means that FFTW is estimating an
  // efficient plan, not that FFTW is estimating the FFT. The result
  // is exact. What it does mean is that the input data is not
  // modified, so we allow the input vector to be const even though
  // we have to cast that const-ness away.
  std::vector<std::complex<double>> *non_const=
    (std::vector<std::complex<double>> *)(&(data[0]));
  fftw_complex *data2=reinterpret_cast<fftw_complex *>(non_const);
  fftw_complex *fft2=reinterpret_cast<fftw_complex *>(&(fft[0]));
  fftw_plan plan=fftw_plan_dft_1d(data.size(),data2,fft2,
                                      FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  
#else
  
  O2SCL_ERR("FFTW support not included in this O2scl installation.",
            o2scl::exc_eunsup);
  
#endif
  
  return;
}

void o2scl::matrix_forward_fft
(size_t m, size_t n, const std::vector<double> &data,
 std::vector<std::complex<double>> &fft) {
  
#ifdef O2SCL_FFTW

  fft.resize(m*(n/2+1));
    
  // First, note that FFTW_ESTIMATE means that FFTW is estimating an
  // efficient plan, not that FFTW is estimating the FFT. The result
  // is exact. What it does mean is that the input data is not
  // modified, so we allow the input vector to be const even though
  // we have to cast that const-ness away.
  double *non_const=(double *)(&(data[0]));
  fftw_complex *fft2=reinterpret_cast<fftw_complex *>(&(fft[0]));
    
  // The forward FFT
  fftw_plan plan=fftw_plan_dft_r2c_2d(m,n,non_const,fft2,
                                      FFTW_ESTIMATE);
  fftw_execute(plan);

#else
    
  O2SCL_ERR("FFTW support not included in this O2scl installation.",
            o2scl::exc_eunsup);
    
#endif

  return;
}

void o2scl::matrix_forward_complex_fft
(size_t m, size_t n, const std::vector<std::complex<double>> &data,
 std::vector<std::complex<double>> &fft) {
  
#ifdef O2SCL_FFTW
  
  fft.resize(m*n);
    
  // First, note that FFTW_ESTIMATE means that FFTW is estimating an
  // efficient plan, not that FFTW is estimating the FFT. The result
  // is exact. What it does mean is that the input data is not
  // modified, so we allow the input vector to be const even though
  // we have to cast that const-ness away.
  std::vector<std::complex<double>> *non_const=
    (std::vector<std::complex<double>> *)(&(data[0]));
  fftw_complex *data2=reinterpret_cast<fftw_complex *>(non_const);
  fftw_complex *fft2=reinterpret_cast<fftw_complex *>(&(fft[0]));
    
  // The forward FFT
  fftw_plan plan=fftw_plan_dft_2d(m,n,data2,fft2,
                                  FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(plan);

#else
    
  O2SCL_ERR("FFTW support not included in this O2scl installation.",
            o2scl::exc_eunsup);
    
#endif

  return;
}

