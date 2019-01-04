/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
/* linalg/householder.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 
 * 2007 Gerard Jungman, Brian Gough
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
/** \file householder_base.h
    \brief File for Householder transformations

    \todo Better documentation for the Householder functions.
*/

#include <gsl/gsl_machine.h>

#include <o2scl/err_hnd.h>
#include <o2scl/cblas.h>
#include <o2scl/permutation.h>

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif
  
  /** \brief Replace the vector \c v with a Householder vector and a
      coefficient tau that annihilates <tt>v[1]</tt> through
      <tt>v[n-1]</tt> (inclusive)
      
      On exit, this function returns the value of \f$ \tau = 2/ (v^{T}
      v) \f$. If \c n is less than or equal to 1 then this function
      returns zero without calling the error handler.
  */
  template<class vec_t> 
    double householder_transform(const size_t n, vec_t &v) {

    if (n <= 1) {
      
      // tau=0 
      return 0.0;
      
    } else { 
      
      double alpha, beta, tau;
      double xnorm=O2SCL_CBLAS_NAMESPACE::dnrm2_subvec(n,v,1);
    
      if (xnorm == 0) {
	// tau=0 
	return 0.0;
      }
    
      alpha=O2SCL_IX(v,0);
      beta=-(alpha >= 0.0 ? +1.0 : -1.0)*hypot(alpha,xnorm);
      tau=(beta-alpha)/beta;
    
      double dbl_eps=std::numeric_limits<double>::epsilon();
      double dbl_min=std::numeric_limits<double>::min();

      double s=(alpha-beta);
      if (fabs(s)>dbl_min) {
	O2SCL_CBLAS_NAMESPACE::dscal_subvec(1.0/s,n,v,1);
	O2SCL_IX(v,0)=beta;
      } else {
	O2SCL_CBLAS_NAMESPACE::dscal_subvec(dbl_eps/s,n,v,1);
	O2SCL_CBLAS_NAMESPACE::dscal_subvec(1.0/dbl_eps,n,v,1);
	O2SCL_IX(v,0)=beta;
      }
      
      return tau;
    }
  }

  /** \brief Compute the Householder transform of a vector
      formed with \c n rows of a column of a matrix
      
      This performs a Householder transform of a vector defined by a
      column of a matrix \c A which starts at element
      <tt>A(ir,ic)</tt> and ends at element <tt>A(M-1,ic)</tt>.
      If <tt>M-1</tt> is equal to <tt>ir+1</tt>, this function quietly
      does nothing.

      Used in \ref QR_decomp() and \ref SV_decomp_mod().
  */
  template<class mat_t> 
    double householder_transform_subcol(mat_t &A, const size_t ir,
					const size_t ic, const size_t M) {
    
    if (M == ir+1) {
      // [GSL] tau=0 
      return 0.0;
    } else { 
      double alpha, beta, tau;
      
      double xnorm=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,ir+1,ic,M);
      
      if (xnorm == 0) {
	// [GSL] tau=0 
	return 0.0;
      }
    
      alpha=O2SCL_IX2(A,ir,ic);
      beta=-(alpha >= 0.0 ? +1.0 : -1.0)*hypot(alpha,xnorm);
      tau=(beta-alpha)/beta;
      
      double dbl_eps=std::numeric_limits<double>::epsilon();
      double dbl_min=std::numeric_limits<double>::min();

      double s=(alpha-beta);
      if (fabs(s)>dbl_min) {
	O2SCL_CBLAS_NAMESPACE::dscal_subcol(A,ir+1,ic,M,1.0/s);
	O2SCL_IX2(A,ir,ic)=beta;
      } else {
	O2SCL_CBLAS_NAMESPACE::dscal_subcol(A,ir+1,ic,M,dbl_eps/s);
	O2SCL_CBLAS_NAMESPACE::dscal_subcol(A,ir+1,ic,M,1.0/dbl_eps);
	O2SCL_IX2(A,ir,ic)=beta;
      }
    
      return tau;
    }
  }

  /** \brief Compute the Householder transform of a vector
      formed with the last \c n columns of a row of a matrix

      This performs a Householder transform of a vector 
      defined by a row of a matrix \c A which starts at element
      <tt>A(ir,ic)</tt> and ends at element <tt>A(ir,N-1)</tt>
      If <tt>N-1</tt> is equal to <tt>ic</tt>, this function quietly
      does nothing.
  */
  template<class mat_t> 
    double householder_transform_subrow(mat_t &A, const size_t ir,
					const size_t ic, const size_t N) {
    
    if (N == ic+1) {
      // [GSL] tau=0 
      return 0.0;
    } else { 
      double alpha, beta, tau;
      
      double xnorm=O2SCL_CBLAS_NAMESPACE::dnrm2_subrow(A,ir,ic+1,N);
      
      if (xnorm == 0) {
	// [GSL] tau=0 
	return 0.0;
      }
    
      alpha=O2SCL_IX2(A,ir,ic);
      beta=-(alpha >= 0.0 ? +1.0 : -1.0)*hypot(alpha,xnorm);
      tau=(beta-alpha)/beta;
      
      double dbl_eps=std::numeric_limits<double>::epsilon();
      double dbl_min=std::numeric_limits<double>::min();

      double s=(alpha-beta);
      if (fabs(s)>dbl_min) {
	O2SCL_CBLAS_NAMESPACE::dscal_subrow(A,ir,ic+1,N,1.0/s);
	O2SCL_IX2(A,ir,ic)=beta;
      } else {
	O2SCL_CBLAS_NAMESPACE::dscal_subrow(A,ir,ic+1,N,dbl_eps/s);
	O2SCL_CBLAS_NAMESPACE::dscal_subrow(A,ir,ic+1,N,1.0/dbl_eps);
	O2SCL_IX2(A,ir,ic)=beta;
      }
    
      return tau;
    }
  }

  /** \brief Apply a Householder transformation \f$ (v,\tau) \f$ to 
      matrix \f$ A \f$ of size <tt>(M,N)</tt>
      
      The vector \c v must have at least \c N entries, with the
      exception that the vector element <tt>v[0]</tt> is never
      referenced by this function.
  */
  template<class vec_t, class mat_t>
    void householder_hm(const size_t M, const size_t N, 
			double tau, const vec_t &v, mat_t &A) {
    
    if (tau == 0.0) {
      return;
    }
    
    for (size_t j=0;j<N;j++) {
      
      // [GSL] Compute wj=Akj vk 
      double wj=O2SCL_IX2(A,0,j);
      
      for (size_t i=1;i<M;i++) {
	wj+=O2SCL_IX2(A,i,j)*O2SCL_IX(v,i);
      }
      
      // [GSL] Aij=Aij-tau vi wj 
      
      // [GSL] i=0 
      O2SCL_IX2(A,0,j)-=tau*wj;
      
      // [GSL] i=1 .. M-1 
      for (size_t i=1;i<M;i++) {
	O2SCL_IX2(A,i,j)-=tau*wj*O2SCL_IX(v,i);
      }
    }
    return;
  }

  /** \brief Apply a Householder transformation to the lower-right
      part of \c M when the transformation is stored in a column of 
      \c M2
      
      This applies a householder transformation <tt>(v,tau)</tt> to a
      lower-right submatrix of \c M. The submatrix has <tt>nr-ir</tt>
      rows and <tt>nc-ic</tt> columns and starts at row \c ir of
      column \c ic of the original matrix \c M. The vector containing
      the transformation is taken from a column of \c M2 starting at
      row \c ir2 and column \c ic2. The matrix \c M2 must have at
      least <tt>ic2+1</tt> columns and at least <tt>nr-ir+ir2</tt>
      rows.

      This function is used in \ref QR_decomp() and \ref QR_unpack() .
  */
  template<class mat_t>
    void householder_hm_subcol(mat_t &M, const size_t ir,
			       const size_t ic, const size_t nr,
			       const size_t nc, const mat_t &M2,
			       const size_t ir2, const size_t ic2,
			       double tau) {

    if (tau == 0.0) {
      return;
    }
    
    for (size_t j=ic;j<nc;j++) {
      
      // [GSL] Compute wj=Akj vk 
      double wj=O2SCL_IX2(M,ir,j);
      
      for (size_t i=ir+1;i<nr;i++) {
	wj+=O2SCL_IX2(M,i,j)*O2SCL_IX2(M2,i-ir+ir2,ic2);
      }
      
      // [GSL] Aij=Aij-tau vi wj 
      
      // [GSL] i=0 
      O2SCL_IX2(M,ir,j)-=tau*wj;
      
      // [GSL] i=1 .. M-1 
      for (size_t i=ir+1;i<nr;i++) {
	O2SCL_IX2(M,i,j)-=tau*wj*O2SCL_IX2(M2,i-ir+ir2,ic2);
      }
    }
    return;
  }

  /** \brief Apply a Householder transformation to the lower-right
      part of \c M when the transformation is stored in a row of \c M2
      
      This applies a householder transformation <tt>(v,tau)</tt> to a
      lower-right submatrix of \c M. The submatrix has <tt>nr-ir</tt>
      rows and <tt>nc-ic</tt> columns and starts at row \c ir of
      column \c ic of the original matrix \c M. The vector containing
      the transformation is taken from a row of \c M2 starting at row
      \c ir2 and column \c ic2. The matrix \c M2 must have
      at least <tt>ir2+1</tt> rows and <tt>nr-ir+ic2</tt> columns.

      Used in \ref bidiag_unpack().
  */
  template<class mat_t>
    void householder_hm_subrow(mat_t &M, const size_t ir,
			       const size_t ic, const size_t nr,
			       const size_t nc, const mat_t &M2,
			       const size_t ir2, const size_t ic2,
			       double tau) {
    
    if (tau == 0.0) {
      return;
    }
    
    for (size_t j=ic;j<nc;j++) {
      
      // [GSL] Compute wj=Akj vk 
      double wj=O2SCL_IX2(M,ir,j);
      
      for (size_t i=ir+1;i<nr;i++) {
        wj+=O2SCL_IX2(M,i,j)*O2SCL_IX2(M2,ir2,i-ir+ic2);
      }
      
      // [GSL] Aij=Aij-tau vi wj 
      
      // [GSL] i=0 
      O2SCL_IX2(M,ir,j)-=tau*wj;
      
      // [GSL] i=1 .. M-1 
      for (size_t i=ir+1;i<nr;i++) {
        O2SCL_IX2(M,i,j)-=tau*wj*O2SCL_IX2(M2,ir2,i-ic+ic2);
      }
    }
    return;
  }

  /** \brief Apply a Householder transformation \c v to vector \c w 
   */
  template<class vec_t, class vec2_t>
    void householder_hv(const size_t N, double tau, const vec_t &v, 
		       vec2_t &w) {
    
    if (tau==0) return;
    
    // compute d=v'w
    
    double d0=O2SCL_IX(w,0);
    double d1, d;
      
    d1=O2SCL_CBLAS_NAMESPACE::ddot_subvec(N,v,w,1);
      
    d=d0+d1;
      
    // compute w=w-tau(v)(v'w)
    O2SCL_IX(w,0)-=tau*d;
    
    O2SCL_CBLAS_NAMESPACE::daxpy_subvec(-tau*d,N,v,w,1);

    return;
  }

  /** \brief Apply a Householder transformation \c v to vector \c w 
      where \c v is stored as a column in a matrix \c A
      
      Used in \ref QR_QTvec().
  */
  template<class mat_t, class vec_t>
    void householder_hv_subcol(const mat_t &A, vec_t &w, double tau,
			       const size_t ie, const size_t N) {
    
    if (tau==0) return;
    
    // compute d=v'w
    
    double d0=O2SCL_IX(w,ie);
    double d1, d;
      
    d1=O2SCL_CBLAS_NAMESPACE::ddot_subcol(N,A,ie+1,ie,w);
      
    d=d0+d1;
      
    // compute w=w-tau(v)(v'w)
    O2SCL_IX(w,ie)-=tau*d;
    
    O2SCL_CBLAS_NAMESPACE::daxpy_subcol(-tau*d,N,A,ie+1,ie,w);

    return;
  }

  /** \brief Apply a Householder transformation \f$ (v,\tau) \f$ to a
      matrix being build up from the identity matrix, using the first
      column of A as a Householder vector 
  */
  template<class mat_t> 
    void householder_hm1(const size_t M, const size_t N, 
			 double tau, mat_t &A) {
    
    if (tau == 0) {
      O2SCL_IX2(A,0,0)=1.0;
      for (size_t j=1;j<N;j++) {
	O2SCL_IX2(A,0,j)=0.0;
      }
      for (size_t i=1;i<M;i++) {
	O2SCL_IX2(A,i,0)=0.0;
      }
      return;
    }

    // [GSL] w=A' v 

    for (size_t j=1;j<N;j++) {
      double wj=0.0;
      for (size_t i=1;i<M;i++) {
	wj+=O2SCL_IX2(A,i,j)*O2SCL_IX2(A,i,0);
      }
      
      // [GSL] A=A-tau v w' 
      O2SCL_IX2(A,0,j)=-tau*wj;
      
      for (size_t i=1;i<M;i++) {
	O2SCL_IX2(A,i,j)-=tau*wj*O2SCL_IX2(A,i,0);
      }
    }
    
    for (size_t i=1;i<M;i++) {
      O2SCL_IX2(A,i,0)=-tau*O2SCL_IX2(A,i,0);
    }
    O2SCL_IX2(A,0,0)=1.0-tau;
    
    return;
  }
  
  /** \brief Apply a Householder transformation \f$ (v,\tau) \f$ to a
      matrix being build up from the identity matrix, using the first
      column of A as a Householder vector 

      Used in \ref SV_decomp_mod() and \ref bidiag_unpack2().
  */
  template<class mat_t> 
    void householder_hm1_sub(const size_t M, const size_t N,
			     double tau, mat_t &A, size_t ir, size_t ic) {
    
    size_t irp1=ir+1;
    size_t icp1=ic+1;

    if (tau == 0) {
      O2SCL_IX2(A,ir,ic)=1.0;
      for (size_t j=icp1;j<N;j++) {
	O2SCL_IX2(A,ir,j)=0.0;
      }
      for (size_t i=irp1;i<M;i++) {
	O2SCL_IX2(A,i,ic)=0.0;
      }
      return;
    }

    // [GSL] w=A' v 

    for (size_t j=icp1;j<N;j++) {
      double wj=0.0;
      for (size_t i=irp1;i<M;i++) {
	wj+=O2SCL_IX2(A,i,j)*O2SCL_IX2(A,i,ic);
      }
      
      // [GSL] A=A-tau v w' 
      O2SCL_IX2(A,ir,j)=-tau*wj;
      
      for (size_t i=irp1;i<M;i++) {
	O2SCL_IX2(A,i,j)-=tau*wj*O2SCL_IX2(A,i,ic);
      }
    }
    
    for (size_t i=irp1;i<M;i++) {
      O2SCL_IX2(A,i,ic)=-tau*O2SCL_IX2(A,i,ic);
    }
    O2SCL_IX2(A,ir,ic)=1.0-tau;
    
    return;
  }
  
  /** \brief Apply the Householder transformation <tt>(v,tau)</tt> to
      the right-hand side of the matrix \c A.
  */
  template<class vec_t, class mat_t>
    void householder_mh(const size_t M, const size_t N,
			double tau, const vec_t &v, mat_t &A) {
    
    if (tau==0.0) {
      return;
    }

    // [GSL] A=A-tau w  v' 
    size_t i, j;
    
    for (i=0;i<M;i++) {
      
      double wi=O2SCL_IX2(A,i,0);

      // [GSL] note, computed for v(0)=1 above 
      for (j=1;j<N;j++)  {
	wi+=O2SCL_IX2(A,i,j)*O2SCL_IX(v,j);
      }
      
      // [GSL] j=0 
      O2SCL_IX2(A,i,0)-=tau*wi;
      
      // [GSL] j=1 .. N-1 
      for (j=1;j<N;j++) {
	O2SCL_IX2(A,i,j)-=tau*wi*O2SCL_IX(v,j);
      }
    }
    
    return;
  }
  
  /** \brief Apply the Householder transformation <tt>(v,tau)</tt> to
      the right-hand side of the matrix \c A.

      Used in \ref bidiag_decomp().
  */
  template<class mat_t, class mat2_t>
    void householder_mh_subrow
    (mat_t &M, const size_t ir, const size_t ic,
     const size_t nr, const size_t nc, const mat2_t &M2,
     const size_t ir2, const size_t ic2, double tau) {
    
    if (tau==0.0) {
      return;
    }

    // [GSL] A=A-tau w  v' 
    size_t i, j, last=nc;
    
    for (i=ir;i<nr;i++) {
      
      double wi=O2SCL_IX2(M,i,ic);

      // [GSL] note, computed for v(0)=1 above 
      for (j=ic+1;j<last;j++)  {
	wi+=O2SCL_IX2(M,i,j)*O2SCL_IX2(M2,ir2,j-ic+ic2);
      }
      
      // [GSL] j=0 
      O2SCL_IX2(M,i,ic)-=tau*wi;
      
      // [GSL] j=1 .. N-1 
      for (j=ic+1;j<last;j++) {
	O2SCL_IX2(M,i,j)-=tau*wi*O2SCL_IX2(M2,ir2,j-ic+ic2);
      }
    }

    return;
  }
  
#ifdef DOXYGEN
}
#endif
