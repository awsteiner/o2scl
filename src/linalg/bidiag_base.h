/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2024, Andrew W. Steiner

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
/* linalg/bidiag.c
 * 
 * Copyright (C) 2001, 2007 Brian Gough
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
/** \file bidiag_base.h
    \brief File defining bidiagonalization functions
*/
#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Factor a matrix into bidiagonal form

      Factor matrix A of size <tt>(M,N)</tt> with \f$ M\geq N \f$ into
      \f$ A = U B V^T \f$ where U and V are orthogonal and B is upper
      bidiagonal.

      After the function call, the matrix \f$ B \f$ is stored the
      diagonal and first superdiagonal of \c A. The matrices \f$ U \f$
      and \f$ V \f$ are stored as packed sets of Householder
      transformations in the lower and upper triangular parts of \c A,
      respectively.

      \verbatim embed:rst
      Adapted from the GSL version which was based on algorithm 5.4.2
      in [Golub96]_.
      \endverbatim
  */
  template<class mat_t, class vec_t, class vec2_t> 
    int bidiag_decomp(size_t M, size_t N, mat_t &A, 
		      vec_t &tau_U, vec2_t &tau_V) {
    
    if (M<N) {
      O2SCL_ERR2("Bidiagonal decomposition requires M>=N in ",
		     "bidiag_comp().",o2scl::exc_ebadlen);
    }
    
    for (size_t i=0;i<N;i++) {

      // [GSL] Apply Householder transformation to current column
      double tau_i=householder_transform_subcol(A,i,i,M);

      // [GSL] Apply the transformation to the remaining columns
      if (i+1<N) {
	householder_hm_subcol(A,i,i+1,M,N,A,i,i,tau_i);
      }
      O2SCL_IX(tau_U,i)=tau_i;

      // [GSL] Apply Householder transformation to current row
      
      if (i+1<N) {
	double tau2_i=householder_transform_subrow(A,i,i+1,N);
	
	// [GSL] Apply the transformation to the remaining rows
	if (i+1<M) {
	  householder_mh_subrow(A,i+1,i+1,M,N,A,i,i+1,tau2_i);
	}
	O2SCL_IX(tau_V,i)=tau2_i;
      }

    }
    
    return o2scl::success;
  }

  /** \brief Unpack a matrix \c A with the bidiagonal decomposition
      and create matrices \c U, \c V, diagonal \c diag and
      superdiagonal \c superdiag

      Given a matrix \c A of size <tt>(M,N)</tt> with \f$ M \geq N \f$
      created by \ref bidiag_decomp(), this function creates the
      matrix \c U of size <tt>(M,N)</tt>, the matrix \c V of size
      <tt>(N,N)</tt>, the diagonal \c diag of size \c N and the
      super-diagonal \c superdiag of size \c N-1.
  */
  template<class mat_t, class vec_t, class mat2_t, class vec2_t, 
    class mat3_t, class vec3_t, class vec4_t>
    int bidiag_unpack(size_t M, size_t N, const mat_t &A, 
		      const vec_t &tau_U, mat2_t &U, const vec2_t &tau_V, 
		      mat3_t &V, vec3_t &diag, vec4_t &superdiag) {
    
    if (M<N) {
      O2SCL_ERR2("Matrix A must have M >= N in ",
		     "bidiag_unpack().",o2scl::exc_ebadlen);
    }
      
    const size_t K=N;
    
    // [GSL] Copy diagonal into diag
    
    for (size_t i=0;i<N;i++) {
      O2SCL_IX(diag,i)=O2SCL_IX2(A,i,i);
    }
    
    // [GSL] Copy superdiagonal into superdiag
    
    for (size_t i=0;i<N-1;i++) {
      O2SCL_IX(superdiag,i)=O2SCL_IX2(A,i,i+1);
    }
    
    // [GSL] Initialize V to the identity
    for(size_t i=0;i<N;i++) {
      for(size_t j=0;j<N;j++) {
	if (i==j) O2SCL_IX2(V,i,j)=1.0;
	else O2SCL_IX2(V,i,j)=0.0;
      }
    }
    
    for (size_t i=N-1;i-- > 0;) {

      // [GSL] Householder row transformation to accumulate V
      householder_hm_subrow(V,i+1,i+1,N,N,A,i,i+1,O2SCL_IX(tau_V,i));
    }
    
    // [GSL] Initialize U to the identity
    
    for(size_t i=0;i<M;i++) {
      for(size_t j=0;j<N;j++) {
	if (i==j) O2SCL_IX2(U,i,j)=1.0;
	else O2SCL_IX2(U,i,j)=0.0;
      }
    }
    
    for (size_t j=N;j-- > 0;) {
      householder_hm_subcol(U,j,j,M,N,A,j,j,O2SCL_IX(tau_U,j));
    }
    
    return o2scl::success;
  }

  /** \brief Unpack a matrix \c A with the bidiagonal decomposition
      and create matrix \c V
  */
  template<class mat_t, class vec_t, class vec2_t, class mat2_t> 
    int bidiag_unpack2(size_t M, size_t N, mat_t &A, vec_t &tau_U, 
		       vec2_t &tau_V, mat2_t &V) {

    if (M<N) {
      O2SCL_ERR2("Matrix A must have M >= N in ",
		     "bidiag_unpack2().",o2scl::exc_ebadlen);
    }

    const size_t K=M;

    // [GSL] Initialize V to the identity

    for(size_t i=0;i<N;i++) {
      for(size_t j=0;j<N;j++) {
	if (i==j) O2SCL_IX2(V,i,j)=1.0;
	else O2SCL_IX2(V,i,j)=0.0;
      }
    }

    for (size_t i=N-1;i-- > 0;) {

      // [GSL] Householder row transformation to accumulate V
      householder_hm_subrow(V,i+1,i+1,N,N,A,i,i+1,O2SCL_IX(tau_V,i));

    }
      
    // [GSL] Copy superdiagonal into tau_v
      
    for (size_t i=0;i<N-1;i++) {
      O2SCL_IX(tau_V,i)=O2SCL_IX2(A,i,i+1);
    }
      
    // [GSL] Allow U to be unpacked into the same memory as A, copy
    // diagonal into tau_U
      
    for (size_t j=N; j-- > 0;) {
      // [GSL] Householder column transformation to accumulate U
      double tj=O2SCL_IX(tau_U,j);
      O2SCL_IX(tau_U,j)=O2SCL_IX2(A,j,j);
      householder_hm1_sub(M,N,tj,A,j,j);
    }
      
    return o2scl::success;
  }

  /** \brief Unpack the diagonal and superdiagonal of the bidiagonal
      decomposition of \c A into \c diag and \c superdiag
  */
  template<class mat_t, class vec_t, class vec2_t>
    int bidiag_unpack_B(size_t M, size_t N, const mat_t &A, 
			vec_t &diag, vec2_t &superdiag) {

    size_t K=N;
    if (M<=N) K=M;

    // [GSL] Copy diagonal into diag 
    for (size_t i=0;i<K;i++) {
      O2SCL_IX(diag,i)=O2SCL_IX2(A,i,i);
    }
    
    // [GSL] Copy superdiagonal into superdiag
    for (size_t i=0;i<K-1;i++) {
      O2SCL_IX(superdiag,i)=O2SCL_IX2(A,i,i+1);
    }

    return o2scl::success;
  }

#ifdef DOXYGEN
}
#endif

