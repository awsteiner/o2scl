/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef O2SCL_GSL_FFT_H
#define O2SCL_GSL_FFT_H

#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Real mixed-radix fast Fourier transform

      This is a simple wrapper for the GSL FFT functions which
      automatically allocates the necessary memory.

      \future Generalize to generic vector types. 
  */
  class gsl_fft {
  public:

    gsl_fft();

    virtual ~gsl_fft();

    /** \brief Perform the FFT transform
     */
    int transform(int n, double *x);

    /** \brief Perform the inverse FFT transform
     */
    int inverse_transform(int n, double *x);
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /// The current memory size
    int mem_size;

    /// Reallocate memory 
    int mem_resize(int new_size);

    /// The GSL workspace
    gsl_fft_real_workspace *work;

    /// The table for the forward transform
    gsl_fft_real_wavetable *real;

    /// The table for the inverse transform
    gsl_fft_halfcomplex_wavetable *hc;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
