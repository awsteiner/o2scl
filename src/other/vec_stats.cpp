/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2021-2025, Andrew W. Steiner and Jesse Farr
  
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

void o2scl::vector_forward_fft_cuda(const std::vector<double> &data,
                                    std::vector<std::complex<double>> &fft) {

#ifdef O2SCL_SET_CUDA
  
  int ret=vector_forward_fft_cuda_nothrow(data.size(),data,fft);
  if (ret!=0) {
    std::string errs=((std::string)"Error number ")+o2scl::itos(ret)+
      " in vector_forward_fft_cuda().";
    O2SCL_ERR(errs.c_str(),o2scl::exc_einval);
  }
  
#else
  
  O2SCL_ERR("CUDA support not included in this O2scl installation.",
            o2scl::exc_eunsup);
  
#endif
  
  return;
}

void o2scl::vector_backward_fft_cuda
(const std::vector<std::complex<double>> &fft,
 std::vector<double> &data) {

#ifdef O2SCL_SET_CUDA
  
  int ret=vector_backward_fft_cuda_nothrow(data.size(),data,fft);
  if (ret!=0) {
    std::string errs=((std::string)"Error number ")+o2scl::itos(ret)+
      " in vector_backward_fft_cuda().";
    O2SCL_ERR(errs.c_str(),o2scl::exc_einval);
  }

#else
  
  O2SCL_ERR("CUDA support not included in this O2scl installation.",
            o2scl::exc_eunsup);
  
#endif
  
  return;
}

void o2scl::vector_forward_fft(const std::vector<double> &data,
                               std::vector<std::complex<double>> &fft) {
  
#ifdef O2SCL_SET_FFTW
  
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

void o2scl::vector_backward_fft(const std::vector<std::complex<double>> &data,
                               std::vector<double> &fft) {  
  
#ifdef O2SCL_SET_FFTW
  
  fft.resize((data.size()-1)*2);
  
  // First, note that FFTW_ESTIMATE means that FFTW is estimating an
  // efficient plan, not that FFTW is estimating the FFT. The result
  // is exact. What it does mean is that the input data is not
  // modified, so we allow the input vector to be const even though
  // we have to cast that const-ness away.
  fftw_complex *data2=(fftw_complex *)(reinterpret_cast<const
                                       fftw_complex *>(&(data[0])));
  double *fft2=(double *)(&(fft[0]));
  // Note that the c2r transforms don't preserve input by default
  fftw_plan plan=fftw_plan_dft_c2r_1d(fft.size(),data2,fft2,
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
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
  
#ifdef O2SCL_SET_FFTW
  
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

void o2scl::vector_backward_complex_fft
(const std::vector<std::complex<double>> &data,
 std::vector<std::complex<double>> &fft) {
  
#ifdef O2SCL_SET_FFTW
  
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
                                      FFTW_BACKWARD,FFTW_ESTIMATE);
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
  
#ifdef O2SCL_SET_FFTW

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

void o2scl::matrix_backward_fft_copy
(size_t m, size_t n, const std::vector<std::complex<double>> &data,
  std::vector<double> &fft) {
  
#ifdef O2SCL_SET_FFTW
  
  fft.resize(m*(n-1)*2);

  // AWS, 5/11/23: This transform crashes unless I make a copy of the
  // input data, even if I use FFTW_PRESERVE_INPUT. Maybe some extra
  // array padding is required in that case? Until I resolve that
  // problem, this code works for now
  
  fftw_complex *data2=(fftw_complex *)fftw_malloc(m*n*sizeof(fftw_complex));
  for(size_t i=0;i<m*n;i++) {
    data2[i][0]=data[i].real();
    data2[i][1]=data[i].imag();
  }

  //fftw_complex *data2=(fftw_complex *)(reinterpret_cast<const
  //fftw_complex *>(&(data[0])));
  
  double *fft2=(double *)(&(fft[0]));
  
  // Note that FFTW_ESTIMATE means that FFTW is estimating an
  // efficient plan, not that FFTW is estimating the FFT.
  
  fftw_plan plan=fftw_plan_dft_c2r_2d(m,(n-1)*2,data2,fft2,
                                      FFTW_ESTIMATE);
  //| FFTW_PRESERVE_INPUT);
  
  fftw_execute(plan);

  fftw_free(data2);
    
#else

  O2SCL_ERR("FFTW support not included in this O2scl installation.",
            o2scl::exc_eunsup);
    
#endif

  return;
}

void o2scl::matrix_forward_complex_fft
(size_t m, size_t n, const std::vector<std::complex<double>> &data,
 std::vector<std::complex<double>> &fft) {
  
#ifdef O2SCL_SET_FFTW
  
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

void o2scl::matrix_backward_complex_fft
(size_t m, size_t n, 
 const std::vector<std::complex<double>> &data,
 std::vector<std::complex<double>> &fft) {  
#ifdef O2SCL_SET_FFTW
  
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
                                  FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(plan);

#else
    
  O2SCL_ERR("FFTW support not included in this O2scl installation.",
            o2scl::exc_eunsup);
    
#endif

  return;
}

