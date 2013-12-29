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
/* linalg/hh.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
/** \file hh_base.h
    \brief File for householder solver
*/
#include <o2scl/err_hnd.h>
#include <o2scl/cblas.h>
#include <o2scl/permutation.h>

#ifdef DOXYGENP
namespace o2scl_linalg {
#endif

  /** \brief Solve a linear system after Householder decomposition in place

      \future Handle memory allocation like the tridiagonal 
      functions.
  */
  template<class mat_t, class vec_t>
    int HH_svx(size_t N, size_t M, mat_t &A, vec_t &x) {
    size_t i, j, k;

    boost::numeric::ublas::vector<double> d(N);

    /* Perform Householder transformation. */

    for (i = 0; i < N; i++) {
      const double aii = O2SCL_IX2(A,i,i);
      double alpha;
      double f;
      double ak;
      double max_norm = 0.0;
      double r = 0.0;

      for (k = i; k < M; k++) {
	double aki = O2SCL_IX2(A,k,i);
	r += aki * aki;
      }

      if (r == 0.0) {
	/* Rank of matrix is less than size1. */
	O2SCL_ERR_RET("Matrix is rank deficient in HH_svx().",
		      o2scl::exc_esing);
      }

      alpha = sqrt (r) * GSL_SIGN (aii);

      ak = 1.0 / (r + alpha * aii);
      O2SCL_IX2(A,i,i)=aii+alpha;
      
      d[i] = -alpha;

      for (k = i + 1; k < N; k++) {
	double norm = 0.0;
	f = 0.0;
	for (j = i; j < M; j++) {
	  double ajk = O2SCL_IX2(A,j,k);
	  double aji = O2SCL_IX2(A,j,i);
	  norm += ajk * ajk;
	  f += ajk * aji;
	}
	max_norm = GSL_MAX (max_norm, norm);

	f *= ak;

	for (j = i; j < M; j++) {
	  double ajk = O2SCL_IX2(A,j,k);
	  double aji = O2SCL_IX2(A,j,i);
	  O2SCL_IX2(A,j,k)=ajk-f*aji;
	}
      }

#ifdef O2SCL_CPP11
      double dbl_eps=std::numeric_limits<double>::epsilon();
#else 
      double dbl_eps=GSL_DBL_EPSILON;
#endif
      
      if (fabs (alpha) < 2.0*dbl_eps*sqrt(max_norm)) {
	/* Apparent singularity. */
	O2SCL_ERR_RET("Apparent singularity in HH_svx().",o2scl::exc_esing);
      }

      /* Perform update of RHS. */

      f = 0.0;
      for (j = i; j < M; j++) {
	f += O2SCL_IX(x,j)*O2SCL_IX2(A,j,i);
      }
      f *= ak;
      for (j = i; j < M; j++) {
	double xj = O2SCL_IX(x,j);
	double aji = O2SCL_IX2(A,j,i);
	O2SCL_IX(x,j)=xj - f * aji;
      }
    }

    /* Perform back-substitution. */
    
    for (i = N; i-- > 0; ) {
      double xi = O2SCL_IX(x,i);
      double sum = 0.0;
      for (k = i + 1; k < N; k++) {
	sum += O2SCL_IX2(A,i,k)*O2SCL_IX(x,k);
      }
    
      O2SCL_IX(x,i)=(xi - sum) / d[i];
    }

    return o2scl::success;
  }
  
  /** \brief Solve linear system after Householder decomposition
  */
  template<class mat_t, class vec_t, class vec2_t>
    int HH_solve(size_t n, mat_t &A, const vec_t &b, vec2_t &x) {
    for(size_t i=0;i<n;i++) O2SCL_IX(x,i)=O2SCL_IX(b,i);
    return HH_svx(n,n,A,x);
  }

#ifdef DOXYGENP
}
#endif
