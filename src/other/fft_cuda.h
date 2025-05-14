/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2025, Andrew W. Steiner
  
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
/** \file fft_cuda.h
    \brief File for CUDA FFTs
*/
#ifndef O2SCL_FFT_CUDA_H
#define O2SCL_FFT_CUDA_H

#include <cufft.h>
#include <vector>
#include <complex>

namespace o2scl {

  /** \brief Desc
   */
  int vector_forward_fft_cuda_nothrow
  (size_t N, const std::vector<double> &data,
   std::vector<std::complex<double>> &fft);

  /** \brief Desc
   */
  int vector_backward_fft_cuda_nothrow
  (size_t N, const std::vector<std::complex<double>> &fft,
   std::vector<double> &data);
  
}

#endif

