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
/* linalg/lu.c
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
/** \file lu_base.h 
    \brief Functions related to LU decomposition
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Return 1 if at least one of the elements in the
      diagonal is zero
   */
  template<class mat_t>
    int diagonal_has_zero(const size_t N, mat_t &A) {
    
    for(size_t i=0;i<N;i++) {
      if (O2SCL_IX2(A,i,i)==0.0) return 1;
    }
    return 0;
  }

  /** \brief Compute the LU decomposition of the matrix \c A

      On output the diagonal and upper triangular part of the input
      matrix A contain the matrix U. The lower triangular part of the
      input matrix (excluding the diagonal) contains L. The diagonal
      elements of L are unity, and are not stored.

      The permutation matrix P is encoded in the permutation p. The
      j-th column of the matrix P is given by the k-th column of the
      identity matrix, where k = p_j the j-th element of the
      permutation vector. The sign of the permutation is given by
      signum. It has the value (-1)^n, where n is the number of
      interchanges in the permutation.
      
      The algorithm used in the decomposition is Gaussian Elimination
      with partial pivoting (Golub & Van Loan, Matrix Computations,
      Algorithm 3.4.1).

      \future The "swap rows j and i_pivot" section could probably
      be made more efficient using a "matrix_row"-like object
      as done in GSL. (7/16/09 - I've tried this, and it doesn't
      seem to improve the speed significantly.)
  */
  template<class mat_t>
    int LU_decomp(const size_t N, mat_t &A, o2scl::permutation &p, 
		  int &signum) {
    
    size_t i, j, k;
  
    signum=1;
    p.init();
  
    for (j = 0; j < N - 1; j++) {
    
      /* Find maximum in the j-th column */
      if (j>=N) O2SCL_ERR("Fail 1.",o2scl::exc_einval);
      double ajj, max = fabs(O2SCL_IX2(A,j,j));
      size_t i_pivot = j;
      
      for (i = j + 1; i < N; i++) {
	if (i>=N || j>=N) O2SCL_ERR("Fail 2.",o2scl::exc_einval);
	double aij = fabs (O2SCL_IX2(A,i,j));
      
	if (aij > max) {
	  max = aij;
	  i_pivot = i;
	}
      }

      if (i_pivot != j) {

	// Swap rows j and i_pivot
	double temp;
	for (k=0;k<N;k++) {
	  if (k>=N || j>=N) O2SCL_ERR("Fail 3.",o2scl::exc_einval);
	  temp=O2SCL_IX2(A,j,k);
	  if (k>=N || i_pivot>=N) O2SCL_ERR("Fail 4.",o2scl::exc_einval);
	  O2SCL_IX2(A,j,k)=O2SCL_IX2(A,i_pivot,k);
	  O2SCL_IX2(A,i_pivot,k)=temp;
	}
	p.swap(j,i_pivot);
	signum=-signum;
      }
    
      ajj = O2SCL_IX2(A,j,j);
      
      if (ajj != 0.0) {
	for (i = j + 1; i < N; i++) {
	  if (j>=N || i>=N) O2SCL_ERR("Fail 5.",o2scl::exc_einval);
	  double aij = O2SCL_IX2(A,i,j) / ajj;
	  O2SCL_IX2(A,i,j)=aij;
	  for (k = j + 1; k < N; k++) {
	    if (j>=N || i>=N || k>=N) O2SCL_ERR("Fail 6.",o2scl::exc_einval);
	    double aik = O2SCL_IX2(A,i,k);
	    double ajk = O2SCL_IX2(A,j,k);
	    O2SCL_IX2(A,i,k)=aik - aij * ajk;
	  }
	}
      }
    }
  
    return o2scl::success;
  }

  /** \brief Solve a linear system after LU decomposition in place
      
      These functions solve the square system A x = b in-place using
      the LU decomposition of A into (LU,p). On input x should contain
      the right-hand side b, which is replaced by the solution on
      output.
  */
  template<class mat_t, class vec_t>
    int LU_svx(const size_t N, const mat_t &LU, 
	       const o2scl::permutation &p, vec_t &x) {
  
    if (diagonal_has_zero(N,LU)) {
      O2SCL_ERR("LU matrix is singular in LU_svx().",
		    o2scl::exc_edom);
    }

    /* Apply permutation to RHS */
    p.apply(x);
  
    /* Solve for c using forward-substitution, L c = P b */
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Lower,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_Unit,
		       N,N,LU,x);
  
    /* Perform back-substitution, U x = c */
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,
		       N,N,LU,x);
  
    return o2scl::success;
  }

  /** \brief An alternate form of LU decomposition with 
      matrix row objects

      \comment
      I was testing this out, but it doesn't seem much faster than
      the original LU_decomp().
      \comment
  */
  template<class mat_t, class mat_row_t>
    int LU_decomp_alt(const size_t N, mat_t &A, o2scl::permutation &p, 
		      int &signum) {
    
    size_t i, j, k;
  
    signum=1;
    p.init();
  
    for (j = 0; j < N - 1; j++) {
    
      /* Find maximum in the j-th column */
      double ajj, max = fabs(O2SCL_IX2(A,j,j));
      size_t i_pivot = j;
      
      for (i = j + 1; i < N; i++) {
	double aij = fabs (O2SCL_IX2(A,i,j));
      
	if (aij > max) {
	  max = aij;
	  i_pivot = i;
	}
      }

      if (i_pivot != j) {

	// Swap rows j and i_pivot
	double temp;
	mat_row_t r1(A,j);
	mat_row_t r2(A,i_pivot);
	for (k=0;k<N;k++) {
	  temp=O2SCL_IX(r1,k);
	  O2SCL_IX(r1,k)=O2SCL_IX(r2,k);
	  O2SCL_IX(r2,k)=temp;
	}
	p.swap(j,i_pivot);
	signum=-signum;
      }
    
      ajj = O2SCL_IX2(A,j,j);
      
      if (ajj != 0.0) {
	for (i = j + 1; i < N; i++) {
	  double aij = O2SCL_IX2(A,i,j) / ajj;
	  O2SCL_IX2(A,i,j)=aij;
	  for (k = j + 1; k < N; k++) {
	    double aik = O2SCL_IX2(A,i,k);
	    double ajk = O2SCL_IX2(A,j,k);
	    O2SCL_IX2(A,i,k)=aik - aij * ajk;
	  }
	}
      }
    }
  
    return o2scl::success;
  }

  /** \brief Solve a linear system after LU decomposition

      This function solve the square system A x = b using the LU
      decomposition of A into (LU, p) given by gsl_linalg_LU_decomp or
      gsl_linalg_complex_LU_decomp.
  */
  template<class mat_t, class vec_t, class vec2_t>
    int LU_solve(const size_t N, const mat_t &LU, const o2scl::permutation &p, 
		 const vec_t &b, vec2_t &x) {
    
    if (diagonal_has_zero(N,LU)) {
      O2SCL_ERR("LU matrix is singular in LU_solve().",
		    o2scl::exc_edom);
    }

    /* Copy x <- b */
    o2scl::vector_copy(N,b,x);
  
    /* Solve for x */
    return LU_svx(N,LU,p,x);
  }

  /** \brief Refine the solution of a linear system
      
      These functions apply an iterative improvement to x, the solution
      of A x = b, using the LU decomposition of A into (LU,p). The
      initial residual r = A x - b is also computed and stored in
      residual.
  */
  template<class mat_t, class mat2_t, class vec_t, class vec2_t, class vec3_t> 
    int LU_refine(const size_t N, const mat_t &A, const mat2_t &LU,
		  const o2scl::permutation &p, const vec_t &b, vec2_t &x,
		  vec3_t &residual) {
  
    if (diagonal_has_zero(N,LU)) {
      O2SCL_ERR("LU matrix is singular in LU_refine().",
		    o2scl::exc_edom);
    }

    /* Compute residual, residual = (A * x  - b) */
    o2scl::vector_copy(N,b,residual);
    o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_NoTrans,
		       N,N,1.0,A,x,-1.0,residual);
  
    /* Find correction, delta = - (A^-1) * residual, and apply it */
  
    int status=LU_svx(N,LU,p,residual);
    o2scl_cblas::daxpy(-1.0,5,residual,x);
  
    return status;
  }

  /** \brief Compute the inverse of a matrix from its LU decomposition

      These functions compute the inverse of a matrix A from its LU
      decomposition (LU,p), storing the result in the matrix inverse.
      The inverse is computed by solving the system A x = b for each
      column of the identity matrix. It is preferable to avoid direct
      use of the inverse whenever possible, as the linear solver
      functions can obtain the same result more efficiently and
      reliably.

      \future Could rewrite to avoid mat_col_t, (9/16/09 - However,
      the function may be faster if mat_col_t is left in, so it's
      unclear what's best.)
  */
  template<class mat_t, class mat2_t, class mat_col_t> 
    int LU_invert(const size_t N, const mat_t &LU, 
		  const o2scl::permutation &p, mat2_t &inverse) {

    if (diagonal_has_zero(N,LU)) {
      O2SCL_ERR("LU matrix is singular in LU_invert().",
		    o2scl::exc_edom);
    }

    size_t i;
  
    int status=o2scl::success;
  
    // Set matrix 'inverse' to the identity
    for(i=0;i<N;i++) {
      for(size_t j=0;j<N;j++) {
	if (j>=N || i>=N) O2SCL_ERR("Fail 7.",o2scl::exc_einval);
	if (i==j) O2SCL_IX2(inverse,i,j)=1.0;
	else O2SCL_IX2(inverse,i,j)=0.0;
      }
    }
  
    for (i = 0; i < N; i++) {
      mat_col_t c(inverse,i);
      int status_i=LU_svx(N,LU,p,c);
    
      if (status_i) {
	status = status_i;
      }
    }
  
    return status;
  }

  /** \brief Compute the determinant of a matrix from its LU decomposition
      
      Compute the determinant of a matrix A from its LU decomposition,
      LU. The determinant is computed as the product of the diagonal
      elements of U and the sign of the row permutation signum.
  */
  template<class mat_t> 
    double LU_det(const size_t N, const mat_t &LU, int signum) {

    size_t i;

    double det=(double)signum;
  
    for (i=0;i<N;i++) {
      det*=O2SCL_IX2(LU,i,i);
    }
  
    return det;
  }

  /** \brief Compute the logarithm of the absolute value of the
      determinant of a matrix from its LU decomposition
      
      Compute the logarithm of the absolute value of the determinant of
      a matrix A, \f$ \ln|\det(A)| \f$, from its LU decomposition, LU.
      This function may be useful if the direct computation of the
      determinant would overflow or underflow.
  */
  template<class mat_t> 
    double LU_lndet(const size_t N, const mat_t &LU) {
  
    size_t i;
    double lndet = 0.0;
  
    for (i = 0; i < N; i++) {
      lndet+=log(fabs(O2SCL_IX2(LU,i,i)));
    }
  
    return lndet;
  }

  /** \brief Compute the sign of the 
      determinant of a matrix from its LU decomposition
      
      Compute the sign or phase factor of the determinant of a matrix
      A, \f$ \det(A)/|\det(A)| \f$, from its LU decomposition, LU.
  */
  template<class mat_t> 
    int LU_sgndet(const size_t N, const mat_t &LU, int signum) {

    size_t i;
    int s=signum;

    for (i=0;i<N;i++) {
      double u=O2SCL_IX2(LU,i,i);
      if (u<0) {
	s*=-1;
      } else if (u == 0) {
	s=0;
	break;
      }
    }

    return s;
  }

#ifdef DOXYGEN
}
#endif
