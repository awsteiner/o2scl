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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/gsl_fft.h>

using namespace std;
using namespace o2scl;

gsl_fft::gsl_fft() {
  mem_size=0;
}

gsl_fft::~gsl_fft() {
  gsl_fft_real_wavetable_free(real);
  gsl_fft_halfcomplex_wavetable_free(hc);
  gsl_fft_real_workspace_free(work);
}

int gsl_fft::mem_resize(int new_size) {
  if (mem_size!=0) {
    gsl_fft_real_wavetable_free(real);
    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);
  }
  mem_size=new_size;
  work=gsl_fft_real_workspace_alloc(mem_size);
  real=gsl_fft_real_wavetable_alloc(mem_size);
  hc=gsl_fft_halfcomplex_wavetable_alloc(mem_size);
  return 0;
}

int gsl_fft::transform(int n, double *x) {
  if (n!=mem_size) mem_resize(n);
  return gsl_fft_real_transform(x,1,n,real,work);
}

int gsl_fft::inverse_transform(int n, double *x) {
  if (n!=mem_size) mem_resize(n);
  return gsl_fft_halfcomplex_inverse(x,1,n,hc,work);
}

