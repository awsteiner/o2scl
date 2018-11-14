/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2018, Andrew W. Steiner

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
/* Cholesky Decomposition
 *
 * Copyright (C) 2000 Thomas Walter
 *
 * 03 May 2000: Modified for GSL by Brian Gough
 * 29 Jul 2005: Additions by Gerard Jungman
 * Copyright (C) 2000,2001, 2002, 2003, 2005, 2007 Brian Gough, 
 * Gerard Jungman
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 */
/** \file cholesky_base.h
    \brief File defining Cholesky decomposition
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Compute the in-place Cholesky decomposition of a symmetric
      positive-definite square matrix
    
      On input, the upper triangular part of A is ignored (only the 
      lower triangular part and diagonal are used). On output,
      the diagonal and lower triangular part contain the matrix L
      and the upper triangular part contains L^T. 
    
      If the matrix is not positive-definite, the error handler 
      will be called, unless \c err_on_fail is false, in which
      case a non-zero value will be returned.
  */
  template<class mat_t> int cholesky_decomp(const size_t M, mat_t &A,
					    bool err_on_fail=true) {
  
    size_t i,j,k;

    /* [GSL] Do the first 2 rows explicitly. It is simple, and faster.
       And one can return if the matrix has only 1 or 2 rows.
    */

    double A_00=O2SCL_IX2(A,0,0);
  
    // AWS: The GSL version stores GSL_NAN in L_00 and then throws
    // an error if A_00 <= 0. We throw the error first and then
    // the square root should always be safe?

    if (A_00<=0.0) {
      if (err_on_fail) {
	O2SCL_ERR2("Matrix not positive definite (A[0][0]<=0) in ",
		   "cholesky_decomp().",o2scl::exc_einval);
      } else {
	return 1;
      }
    }
  
    double L_00=sqrt(A_00);
    O2SCL_IX2(A,0,0)=L_00;
  
    if (M>1) {
      double A_10=O2SCL_IX2(A,1,0);
      double A_11=O2SCL_IX2(A,1,1);
          
      double L_10=A_10/L_00;
      double diag=A_11-L_10*L_10;
    
      if (diag<=0.0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Matrix not positive definite (diag<=0 for 2x2) in ",
		     "cholesky_decomp().",o2scl::exc_einval);
	} else {
	  return 2;
	}
      }
      double L_11=sqrt(diag);

      O2SCL_IX2(A,1,0)=L_10;
      O2SCL_IX2(A,1,1)=L_11;
    }
      
    for (k=2;k<M;k++) {
      double A_kk=O2SCL_IX2(A,k,k);
          
      for (i=0;i<k;i++) {
	double sum=0.0;

	double A_ki=O2SCL_IX2(A,k,i);
	double A_ii=O2SCL_IX2(A,i,i);

	// AWS: Should change to use a form of ddot() here
	if (i>0) {
	  sum=0.0;
	  for(j=0;j<i;j++) {
	    sum+=O2SCL_IX2(A,i,j)*O2SCL_IX2(A,k,j);
	  }
	}

	A_ki=(A_ki-sum)/A_ii;
	O2SCL_IX2(A,k,i)=A_ki;
      } 

      {
	double sum=dnrm2_subrow(A,k,0,k);
	double diag=A_kk-sum*sum;

	if (diag<=0.0) {
	  if (err_on_fail) {
	    O2SCL_ERR2("Matrix not positive definite (diag<=0) in ",
		       "cholesky_decomp().",o2scl::exc_einval);
	  } else {
	    return 3;
	  }
	}

	double L_kk=sqrt(diag);
      
	O2SCL_IX2(A,k,k)=L_kk;
      }
    }

    /* [GSL] Now copy the transposed lower triangle to the upper
       triangle, the diagonal is common.
    */
      
    for (i=1;i<M;i++) {
      for (j=0;j<i;j++) {
	double A_ij=O2SCL_IX2(A,i,j);
	O2SCL_IX2(A,j,i)=A_ij;
      }
    } 
  
    return 0;
  }

  /** \brief Compute the determinant of a matrix from its Cholesky decomposition
      
  */
  template<class mat_t> 
    double cholesky_det(const size_t M, const mat_t &A) {

    double det=1.0;
  
    for (size_t i=0;i<M;i++) {
      det*=O2SCL_IX2(A,i,i);
    }
  
    return det;
  }

  /** \brief Compute the logarithm of the absolute value of the
      determinant of a matrix from its Cholesky decomposition
  */
  template<class mat_t> 
    double cholesky_lndet(const size_t M, const mat_t &A) {
  
    double lndet = 0.0;
  
    for (size_t i = 0; i < M; i++) {
      lndet+=log(fabs(O2SCL_IX2(A,i,i)));
    }
  
    return lndet;
  }

  /** \brief Solve a symmetric positive-definite linear system after a 
      Cholesky decomposition

      Given the Cholesky decomposition of a matrix A in \c LLT, 
      this function solves the system <tt>A*x=b</tt>. 
  */
  template<class mat_t, class vec_t, class vec2_t>
    void cholesky_solve(const size_t N, const mat_t &LLT, 
			const vec_t &b, vec2_t &x) {
  
    // [GSL] Copy x <- b 
    o2scl::vector_copy(N,b,x);
  
    // [GSL] Solve for c using forward-substitution, L c=b 
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Lower, 
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,N,N,LLT,x);
  
    // [GSL] Perform back-substitution,U x=c 
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,N,N,LLT,x);
  
    return;
  }

  /** \brief Solve a linear system in place using a Cholesky decomposition
   */
  template<class mat_t, class vec_t> 
    void cholesky_svx(const size_t N, const mat_t &LLT, vec_t &x) {
  
    // [GSL] Solve for c using forward-substitution, L c=b 
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Lower,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,N,N,LLT,x);
  
    // [GSL] Perform back-substitution, U x=c 
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,N,N,LLT,x);
  
    return;
  }

  /** \brief Compute the inverse of a symmetric positive definite matrix
      given the Cholesky decomposition

      Given a Cholesky decomposition produced by \ref cholesky_decomp(),
      this function returns the inverse of that matrix in \c LLT.
  */
  template<class mat_t> void cholesky_invert(const size_t N, mat_t &LLT) {
  
    size_t i, j;
    double sum;

    // [GSL] invert the lower triangle of LLT
    for (i=0;i<N;++i) {
    
      j=N-i-1;
    
      O2SCL_IX2(LLT,j,j)=1.0/O2SCL_IX2(LLT,j,j);
      double ajj=-O2SCL_IX2(LLT,j,j);
    
      if (j<N-1) {

	// This section is just the equivalent of dtrmv() for a part of
	// the matrix LLT.
	{
	
	  size_t ix=N-j-2;
	  for (size_t ii=N-j-1;ii>0 && ii--;) {
	    double temp=0.0;
	    const size_t j_min=0;
	    const size_t j_max=ii;
	    size_t jx=j_min;
	    for (size_t jj=j_min;jj<j_max;jj++) {
	      temp+=O2SCL_IX2(LLT,jx+j+1,j)*O2SCL_IX2(LLT,ii+j+1,jj+j+1);
	      jx++;
	    }
	    O2SCL_IX2(LLT,ix+j+1,j)=temp+O2SCL_IX2(LLT,ix+j+1,j)*
	      O2SCL_IX2(LLT,ii+j+1,ii+j+1);
	    ix--;
	  }

	}

	o2scl_cblas::dscal_subcol(LLT,j+1,j,N,ajj);

      }
    }
  
    /*
      [GSL] The lower triangle of LLT now contains L^{-1}. Now compute
      A^{-1}=L^{-t} L^{-1}
    
      The (ij) element of A^{-1} is column i of L^{-1} dotted into
      column j of L^{-1}
    */
  
    for (i=0;i<N;++i) {

      for (j=i+1;j<N;++j) {

	// [GSL] Compute Ainv_{ij}=sum_k Linv_{ki} Linv_{kj}.

	// AWS: Should change to use a form of ddot() here
	sum=0.0;
	for(size_t k=j;k<N;k++) {
	  sum+=O2SCL_IX2(LLT,k,i)*O2SCL_IX2(LLT,k,j);
	}
      
	// [GSL] Store in upper triangle
	O2SCL_IX2(LLT,i,j)=sum;
      }
    
      // [GSL] now compute the diagonal element
    
      // AWS: Should change to use a form of ddot() here
      sum=0.0;
      for(size_t k=i;k<N;k++) {
	sum+=O2SCL_IX2(LLT,k,i)*O2SCL_IX2(LLT,k,i);
      }

      O2SCL_IX2(LLT,i,i)=sum;
    }
  
    // [GSL] Copy the transposed upper triangle to the lower triangle 

    for (j=1;j<N;j++) {
      for (i=0;i<j;i++) {
	O2SCL_IX2(LLT,j,i)=O2SCL_IX2(LLT,i,j);
      }
    } 
  
    return;
  }

  /** \brief Cholesky decomposition with unit-diagonal triangular parts.
      
      A = L D L^T, where diag(L) = (1,1,...,1).
      Upon exit, A contains L and L^T as for Cholesky, and
      the diagonal of A is (1,1,...,1). The vector Dis set
      to the diagonal elements of the diagonal matrix D.
  */
  template<class mat_t, class vec_t>
    int cholesky_decomp_unit(const size_t N, mat_t &A, vec_t &D) {
  
    size_t i, j;
  
    // [GSL] Initial Cholesky decomposition
    int stat_chol=cholesky_decomp(N,A);
  
    // [GSL] Calculate D from diagonal part of initial Cholesky
    for(i=0;i<N;++i) {
      const double C_ii=O2SCL_IX2(A,i,i);
      O2SCL_IX(D,i)=C_ii*C_ii;
    }
  
    // [GSL] Multiply initial Cholesky by 1/sqrt(D) on the right 
    for(i=0;i<N;++i) {
      for(j=0;j<N;++j) {
	O2SCL_IX2(A,i,j)=O2SCL_IX2(A,i,j)/sqrt(O2SCL_IX(D,j));
      }
    }
  
    /* [GSL] Because the initial Cholesky contained both L and
       transpose(L), the result of the multiplication is not symmetric
       anymore; but the lower triangle _is_ correct. Therefore we
       reflect it to the upper triangle and declare victory.
    */
    for(i=0;i<N;++i) {
      for(j=i+1;j<N;++j) {
	O2SCL_IX2(A,i,j)=O2SCL_IX2(A,j,i);
      }
    }
  
    return stat_chol;
  }

#ifdef DOXYGEN
}
#endif
