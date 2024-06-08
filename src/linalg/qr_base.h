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
/* linalg/qr.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, 
 * Brian Gough
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
/** \file qr_base.h
    \brief File for QR decomposition and associated solver
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Compute the QR decomposition of matrix \c A
  */
  template<class mat_t, class vec_t>
    void QR_decomp(size_t M, size_t N, mat_t &A, vec_t &tau) {
    
    size_t imax;
    if (M<N) imax=M;
    else imax=N;
    for(size_t i=0;i<imax;i++) {
      O2SCL_IX(tau,i)=householder_transform_subcol(A,i,i,M);
      if (i+1<N) {
	householder_hm_subcol(A,i,i+1,M,N,A,i,i,O2SCL_IX(tau,i));
      }
    }
    return;
  }

  /** \brief Form the product Q^T v from a QR factorized matrix
  */
  template<class mat_t, class vec_t, class vec2_t>
    void QR_QTvec(const size_t M, const size_t N, 
		 const mat_t &QR, const vec_t &tau, vec2_t &v) {
    
    // compute Q^T v 
    size_t imax;
    if (M<N) imax=M;
    else imax=N;
    for(size_t i=0;i<imax;i++) {
      householder_hv_subcol(QR,v,O2SCL_IX(tau,i),i,M);
    }
    return;
  }

  /** \brief Unpack the QR matrix to the individual Q and R components
  */
  template<class mat1_t, class mat2_t, class mat3_t, class vec_t>
    void QR_unpack(const size_t M, const size_t N,
		  const mat1_t &QR, const vec_t &tau, mat2_t &Q, 
		  mat3_t &R) {
    
    size_t i, j;
    
    // [GSL] Initialize Q to the identity 
    
    for(i=0;i<M;i++) {
      for(j=0;j<M;j++) {
	if (i==j) O2SCL_IX2(Q,i,j)=1.0;
	else O2SCL_IX2(Q,i,j)=0.0;
      }
    }
    
    size_t istart;
    if (M<N) istart=M;
    else istart=N;

    for (i=istart; i-- > 0;) {
      householder_hm_subcol(Q,i,i,M,M,QR,i,i,O2SCL_IX(tau,i));
    }
    
    // [GSL] Form the right triangular matrix R from a packed QR matrix
    
    for (i=0; i < M; i++) {
      for (j=0; j < i && j < N; j++) {
	O2SCL_IX2(R,i,j)=0.0;
      }
      for (j=i; j < N; j++) {
	O2SCL_IX2(R,i,j)=O2SCL_IX2(QR,i,j);
      }
    }
    
    return;
  }

  /** \brief Solve the system A x = b in place using the QR factorization
  */
  template<class mat_t, class vec_t, class vec2_t>
    void QR_svx(size_t M, size_t N, const mat_t &QR, const vec_t &tau, 
	       vec2_t &x) {
    QR_QTvec(M,N,QR,tau,x);
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,
		       M,N,QR,x);
    return;
  }

  /** \brief Solve the system A x = b using the QR factorization
  */
  template<class mat_t, class vec_t, class vec2_t, class vec3_t>
    void QR_solve(size_t N, const mat_t &QR, const vec_t &tau, 
		 const vec2_t &b, vec3_t &x) {
    o2scl::vector_copy(N,b,x);
    QR_svx(N,N,QR,tau,x);
    return;
  }
  
  /** \brief Update a QR factorisation for A= Q R,  A' = A + u v^T,

      The parameters \c M and \c N are the number of rows and columns 
      of the matrix \c R.

      \verbatim
      * Q' R' = QR + u v^T
      *       = Q (R + Q^T u v^T)
      *       = Q (R + w v^T)
      *
      * where w = Q^T u.
      *
      * Algorithm from Golub and Van Loan, "Matrix Computations", Section
      * 12.5 (Updating Matrix Factorizations, Rank-One Changes)
      \endverbatim
  */
  template<class mat1_t, class mat2_t, class vec1_t, class vec2_t> 
    void QR_update(size_t M, size_t N, mat1_t &Q, mat2_t &R, 
		  vec1_t &w, vec2_t &v) {
    
    size_t j;
    // Integer to avoid problems with decreasing loop below
    int k;
    double w0;
	 
    /* [GSL] Apply Given's rotations to reduce w to (|w|, 0, 0, ... , 0)
       
       J_1^T .... J_(n-1)^T w = +/- |w| e_1
       
       simultaneously applied to R,  H = J_1^T ... J^T_(n-1) R
       so that H is upper Hessenberg.  (12.5.2) 
    */
	 
    // Loop from k = M-1 to 1
    for (k=(int)(M-1); k > 0; k--) {

      double c, s;
      double wk=O2SCL_IX(w,k);
      double wkm1=O2SCL_IX(w,k-1);
	     
      o2scl_linalg::create_givens(wkm1,wk,c,s);
      apply_givens_vec(w,k-1,k,c,s);
      apply_givens_qr(M,N,Q,R,k-1,k,c,s);
    }
    
    w0=O2SCL_IX(w,0);
	 
    // [GSL] Add in w v^T  (Equation 12.5.3) 
	 
    for (j=0; j < N; j++) {
      double r0j=O2SCL_IX2(R,0,j);
      double vj=O2SCL_IX(v,j);
      O2SCL_IX2(R,0,j)=r0j+w0*vj;
    }
	 
    // [GSL] Apply Givens transformations R' = G_(n-1)^T ... G_1^T H
    // (Equation 12.5.4)
    int kmax;
    if (M<N+1) kmax=M;
    else kmax=N+1;
    for (k=1;k<kmax;k++) {

      double c, s;
      double diag=O2SCL_IX2(R,k-1,k-1);
      double offdiag=O2SCL_IX2(R,k,k-1);
      
      o2scl_linalg::create_givens(diag,offdiag,c,s);
      apply_givens_qr(M,N,Q,R,k-1,k,c,s);
    
      O2SCL_IX2(R,k,k-1)=0.0;
    }
	 
    return;
  }
  
  /** \brief Compute the unpacked QR decomposition of matrix \c A

      If \o2 is compiled with Armadillo support, this is specialized
      for <tt>arma::mat</tt> to use <tt>arma::qr_econ</tt>. If \o2 is
      compiled with Eigen support, this is specialized for
      <tt>Eigen::MatrixXd</tt>.
   */
  template<class mat_t, class mat2_t, class mat3_t>
    void QR_decomp_unpack(const size_t M, const size_t N,
			  mat_t &A, mat2_t &Q, mat3_t &R) {
    
    boost::numeric::ublas::vector<double> tau(M);
    QR_decomp(M,N,A,tau);
    QR_unpack(M,N,A,tau,Q,R);
    
    return;
  }

#ifdef DOXYGEN
}
#endif
