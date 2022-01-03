/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
/** \file lu.h
    \brief Header wrapper for \ref lu_base.h
*/
#ifndef O2SCL_LU_H
#define O2SCL_LU_H

#include <o2scl/err_hnd.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>

namespace o2scl_linalg {
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M(i,j)
#include <o2scl/lu_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2
  
}

namespace o2scl_linalg_bracket {

  /** \brief Specialized version of LU_decomp for C-style 2D arrays

      \note Note that \c N and \c n must be equal, by no checking
      is done to ensure that this is the case
   */
  template<size_t N>
    int LU_decomp_array_2d(const size_t n, double A[][N], 
			   o2scl::permutation &p, int &signum) {
    
    size_t i, j, k;

    signum=1;
    p.init();
  
    for (j = 0; j < N - 1; j++) {
    
      /* Find maximum in the j-th column */
      double ajj, max = fabs(A[j][j]);
      size_t i_pivot = j;
      
      for (i = j + 1; i < N; i++) {
	double aij = fabs (A[i][j]);
      
	if (aij > max) {
	  max = aij;
	  i_pivot = i;
	}
      }

      if (i_pivot != j) {

	// Swap rows j and i_pivot
	double temp;
	double *r1=&(A[j][0]);
	double *r2=&(A[i_pivot][0]);
	for (k=0;k<N;k++) {
	  temp=r1[k];
	  r1[k]=r2[k];
	  r2[k]=temp;
	}
	p.swap(j,i_pivot);
	signum=-signum;
      }
      
      ajj = A[j][j];
      
      if (ajj != 0.0) {
	for (i = j + 1; i < N; i++) {
	  double aij = A[i][j] / ajj;
	  A[i][j]=aij;
	  for (k = j + 1; k < N; k++) {
	    double aik = A[i][k];
	    double ajk = A[j][k];
	    A[i][k]=aik - aij * ajk;
	  }
	}
      }
    }
  
    return o2scl::success;
  }
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M[i][j]
#include <o2scl/lu_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2

}

#endif
