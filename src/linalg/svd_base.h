/*
  -------------------------------------------------------------------

  Copyright (C) 2010-2015, Andrew W. Steiner

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
/* linalg/svd.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard
 * Jungman, Brian Gough
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
/** \file svd_base.h
    \brief File for SVD composition
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Factorise a general matrix into its SV
      decomposition using the Golub-Reinsch algorithm

      This factors matrix \c A of size <tt>(M,N)</tt> into
      \f[
      A = U~D~V^T
      \f]
      where \c U is a column-orthogonal matrix of size <tt>(M,N)</tt>
      (stored in \c A), \c D is a diagonal matrix of size
      <tt>(N,N)</tt> (stored in the vector \c S of size \c N), and \c
      V is a orthogonal matrix of size <tt>(N,N)</tt>. The vector \c
      work is a workspace vector of size \c N. The matrices \c U and
      \c V are constructed so that
      \f[
      U^T~U = I \qquad \mathrm{and} \qquad V^T~V = V~V^T = I
      \f]

      This algorithm requres \f$ M \geq N \f$. 

      \todo Test N=1 case, N=2 case, and non-square matrices.
  */
  template<class mat_t, class mat2_t, class vec_t, class vec2_t>
    void SV_decomp(size_t M, size_t N, 
		   mat_t &A, mat2_t &V, vec_t &S, vec2_t &work) {
    
    if (M<N) {
      O2SCL_ERR2("SVD decomposition of MxN matrix with M<N, ",
		 "is not implemented in SV_decomp().",o2scl::exc_eunimpl);
    }
    
    // K should be the smaller of M and N, which is always N
    const size_t K=N;
    
    // Handle the case of N = 1 (SVD of a column vector) 
    
    if (N==1) {
      double norm=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,0,0,M);
      O2SCL_IX(S,0)=norm;
      O2SCL_IX2(V,0,0)=1.0;
      
      if (norm != 0.0) {
	O2SCL_CBLAS_NAMESPACE::dscal_subcol(A,0,0,M,1.0/norm);
      }
	
      return;
    }
    
    // [GSL] bidiagonalize matrix A, unpack A into U S V */
    bidiag_decomp(M,N,A,S,work);
    
    //std::cout << "A: " << A(0,0) << " " << A(M-1,N-1) << std::endl;
    //std::cout << "S: " << S[0] << " " << S[S.size()-1] << std::endl;

    bidiag_unpack2(M,N,A,S,work,V);

    //std::cout << "S2: " << S[0] << " " << S[S.size()-1] << std::endl;

    // [GSL] apply reduction steps to B=(S,Sd)

    chop_small_elements(N,S,work);

    //std::cout << "S3: " << S[0] << " " << S[S.size()-1] << std::endl;

    // [GSL] Progressively reduce the matrix until it is diagonal
    
    size_t b=N-1;
    size_t iter=0;
    
    while (b>0) {
      
      double fbm1=O2SCL_IX(work,b-1);
      if (fbm1==0.0 || gsl_isnan(fbm1)) {
	b--;
	continue;
      }
      //std::cout << "b,fbm1: " << b << " " << fbm1 << std::endl;
      
      // [GSL] Find the largest unreduced block (a,b) starting from b
      // and working backwards
      
      size_t a=b-1;
      
      while (a>0) {
	
	double fam1=O2SCL_IX(work,a-1);
	if (fam1==0.0 || gsl_isnan(fam1)) {
	  break;
	}
	a--;

	//std::cout << "a,fam1: " << a << " " << fam1 << std::endl;
      }
      
      iter++;
      if (iter>100*N) {
	O2SCL_ERR("SV decomposition failed to converge in SV_decomp().",
		  o2scl::exc_emaxiter);
      }

      int rescale=0;
      double scale=1.0;
      double norm=0.0;
      size_t n_block=b-a+1;

      // [GSL] Find the maximum absolute values of the diagonal
      // and subdiagonal

      for(size_t i=0;i<n_block;i++) {
	double aa=fabs(O2SCL_IX(S,a+i));
	if (aa>norm) norm=aa;
	//std::cout << "aa: " << aa << std::endl;
      }

      for(size_t i=0;i<n_block-1;i++) {
	double aa=fabs(O2SCL_IX(work,a+i));
	if (aa>norm) norm=aa;
	//std::cout << "aa2: " << aa << std::endl;
      }

      // [GSL] Temporarily scale the submatrix if necessary

      if (norm>GSL_SQRT_DBL_MAX) {
	scale=norm/GSL_SQRT_DBL_MAX;
	rescale=1;
      } else if (norm<GSL_SQRT_DBL_MIN && norm>0.0) {
	scale=norm/GSL_SQRT_DBL_MIN;
	rescale=1;
      }

      //std::cout << "rescale: " << rescale << std::endl;

      if (rescale) {
	O2SCL_CBLAS_NAMESPACE::dscal_subvec(1.0/scale,N,S,a);
	O2SCL_CBLAS_NAMESPACE::dscal_subvec(1.0/scale,N,work,a);
      }

      // [GSL] Perform the implicit QR step

      /*
	for(size_t ii=0;ii<M;ii++) {
	for(size_t jj=0;jj<N;jj++) {
	std::cout << ii << "." << jj << "." << A(ii,jj) << std::endl;
	}
	}
	for(size_t ii=0;ii<N;ii++) {
	for(size_t jj=0;jj<N;jj++) {
	std::cout << "V: " << ii << "." << jj << "." << V(ii,jj) << std::endl;
	}
	}
      */
      
      // At this point, we want to perform qrstep on the arguments 
      // assuming S is a vector of size n_block,
      // f is a vector of size n_block, U is matrix of size (M,n_block),
      // and V is a matrix of size (N,n_block). 
      qrstep_sub(M,N,n_block+a,S,work,A,V,a);

      /*
	for(size_t ii=0;ii<M;ii++) {
	for(size_t jj=0;jj<N;jj++) {
	std::cout << ii << " " << jj << " " << A(ii,jj) << std::endl;
	}
	}
	for(size_t ii=0;ii<N;ii++) {
	for(size_t jj=0;jj<N;jj++) {
	std::cout << "V: " << ii << " " << jj << " " << V(ii,jj) << std::endl;
	}
	}
      */
      
      // [GSL] Remove any small off-diagonal elements
      
      chop_small_elements(N,S,work);

      // [GSL] Undo the scaling if needed
      if (rescale) {
	O2SCL_CBLAS_NAMESPACE::dscal_subvec(scale,N,S,a);
	O2SCL_CBLAS_NAMESPACE::dscal_subvec(scale,N,work,a);
      }
    }

    // [GSL] Make singular values positive by reflections if necessary 
  
    for (size_t j=0;j<K;j++) {
      double Sj=O2SCL_IX(S,j);
      if (Sj<0.0) {
	for (size_t i=0;i<N;i++) {
	  O2SCL_IX2(V,i,j)=-O2SCL_IX2(V,i,j);
	}
	O2SCL_IX(S,j)=-Sj;
      }
    }
    //std::cout << "Here8" << std::endl;
  
    // [GSL] Sort singular values into decreasing order 
    
    for (size_t i=0;i<K;i++) {
      
      double S_max=O2SCL_IX(S,i);
      size_t i_max=i;
      
      for (size_t j=i+1;j<K;j++) {
	double Sj=O2SCL_IX(S,j);
	if (Sj>S_max) {
	  S_max=Sj;
	  i_max=j;
	}
      }
      
      if (i_max != i) {
	// [GSL] swap eigenvalues 
	o2scl::vector_swap<vec_t,double>(S,i,i_max);

	// [GSL] swap eigenvectors 
	o2scl::matrix_swap_cols<mat_t,double>(M,A,i,i_max);
	o2scl::matrix_swap_cols<mat2_t,double>(M,V,i,i_max);
      }
    }
    //std::cout << "Here9" << std::endl;
  
    return;
  }
  
  /** \brief SV decomposition by the modified Golub-Reinsch 
      algorithm which is better for \f$ M \gg N \f$

      This factors matrix \c A of size <tt>(M,N)</tt> into
      \f[
      A = U~D~V^T
      \f]
      where \c U is a column-orthogonal matrix of size <tt>(M,N)</tt>
      (stored in \c A), \c D is a diagonal matrix of size
      <tt>(N,N)</tt> (stored in the vector \c S of size \c N), and \c
      V is a orthogonal matrix of size <tt>(N,N)</tt>. The vector \c
      work is a workspace vector of size \c N and the matrix
      \c X is a workspace of size <tt>(N,N)</tt>. 
   */
  template<class mat_t, class mat2_t, class mat3_t, class vec_t, 
    class vec2_t> void SV_decomp_mod
    (size_t M, size_t N, mat_t &A, mat2_t &X, mat3_t &V, vec_t &S, 
     vec2_t &work) {
      
      if (N == 1) {
	double norm=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,0,0,M);
	O2SCL_IX(S,0)=norm;
	O2SCL_IX2(V,0,0)=1.0;
	
	if (norm != 0.0) {
	  O2SCL_CBLAS_NAMESPACE::dscal_subcol(A,0,0,M,1.0/norm);
	}
	
	return;
      }

      // AWS: This next loop is just a QR decomposition. Replace
      // with QR_decomp()?
      
      // [GSL] Convert A into an upper triangular matrix R
      for (size_t i=0;i<N;i++) {
	double tau_i=householder_transform_subcol(A,i,i,M);
	// [GSL] Apply the transformation to the remaining columns
	if (i+1<N) {
	  householder_hm_subcol(A,i,i+1,M,N,A,i,i,tau_i);
	}
	O2SCL_IX(S,i)=tau_i;
      }
      
      // [GSL] Copy the upper triangular part of A into X 
      
      for (size_t i=0;i<N;i++) {
	for (size_t j=0;j<i;j++) {
	  O2SCL_IX2(X,i,j)=0.0;
	}

	O2SCL_IX2(X,i,i)=O2SCL_IX2(A,i,i);
	for (size_t j=i+1;j<N;j++) {
	  O2SCL_IX2(X,i,j)=O2SCL_IX2(A,i,j);
	}
      }
      
      // [GSL] Convert A into an orthogonal matrix L 
      
      for (size_t j=N;j-- > 0;) {
	// [GSL] Householder column transformation to accumulate L 
	double tj=O2SCL_IX(S,j);
	householder_hm1_sub(M,N,tj,A,j,j);
      }
      
      // [GSL] unpack R into X V S 
      
      SV_decomp(N,N,X,V,S,work);
      
      // [GSL] Multiply L by X, to obtain U = L X, stored in U 
      
      for (size_t i=0;i<M;i++) {

	for(size_t j=0;j<N;j++) {
	  O2SCL_IX(work,j)=0.0;
	}

	for (size_t j=0;j<N;j++) {
	  double Lij=O2SCL_IX2(A,i,j);
	  O2SCL_CBLAS_NAMESPACE::daxpy_subrow(Lij,N,X,j,0,work);
	}

	for (size_t j=0;j<N;j++) {
	  O2SCL_IX2(A,i,j)=O2SCL_IX(work,j);
	}

      }
      
      return;
    }
  
  /** \brief Solve the system A x = b using the SV decomposition

      Solves a linear system using the output of \ref SV_decomp().
      Only non-zero singular values are used in computing solution. In
      the over-determined case, \f$ M>N \f$, the system is solved in
      the least-squares sense.
   */
  template<class mat_t, class mat2_t, class vec_t, class vec2_t, class vec3_t>
    void SV_solve(size_t M, size_t N, 
		  mat_t &U, mat2_t &V, vec_t &S, vec2_t &b, vec3_t &x) {

    //if (U->size1 != b->size)
    //else if (U->size2 != S->size)
    //else if (V->size1 != V->size2)
    //else if (S->size != V->size1)
    //else if (V->size2 != x->size)

    double *w=new double[N];
    
    O2SCL_CBLAS_NAMESPACE::dgemv
      (O2SCL_CBLAS_NAMESPACE::o2cblas_RowMajor,
       O2SCL_CBLAS_NAMESPACE::o2cblas_Trans,M,N,1.0,U,b,0.0,w);
    
    for (size_t i=0;i<N;i++) {
      double alpha=O2SCL_IX(S,i);
      if (alpha != 0) alpha=1.0/alpha;
      w[i]*=alpha;
    }
    
    O2SCL_CBLAS_NAMESPACE::dgemv
      (O2SCL_CBLAS_NAMESPACE::o2cblas_RowMajor,
       O2SCL_CBLAS_NAMESPACE::o2cblas_NoTrans,N,N,1.0,V,w,0.0,x);
    
    delete[] w;
    
    return;
  }

  /** \brief SV decomposition using one-sided Jacobi orthogonalization
      
      This factors matrix \c A of size <tt>(M,N)</tt> into
      \f[
      A = U~D~V^T
      \f]
      where \c U is a column-orthogonal matrix of size <tt>(M,N)</tt>
      (stored in \c A), \c D is a diagonal matrix of size
      <tt>(N,N)</tt> (stored in the vector \c S of size \c N), and \c
      V is a orthogonal matrix of size <tt>(N,N)</tt>. 

      This function computes singular values to higher relative
      accuracy than \ref SV_decomp() and \ref SV_decomp_mod().

      \comment
      Algorithm due to J.C. Nash, Compact Numerical Methods for
      Computers (New York: Wiley and Sons, 1979), chapter 3.
      See also Algorithm 4.1 in
      James Demmel, Kresimir Veselic, "Jacobi's Method is more
      accurate than QR", Lapack Working Note 15 (LAWN15), October 1989.
      Available from netlib.
      
      Based on code by Arthur Kosowsky, Rutgers University
      kosowsky@physics.rutgers.edu  
      
      Another relevant paper is, P.P.M. De Rijk, "A One-Sided Jacobi
      Algorithm for computing the singular value decomposition on a
      vector computer", SIAM Journal of Scientific and Statistical
      Computing, Vol 10, No 2, pp 359-371, March 1989.
      \endcomment

      \future There were originally a few GSL_COERCE_DBL calls which
      have been temporarily removed and could be restored.
  */
  template<class mat_t, class mat2_t, class vec_t>
    void SV_decomp_jacobi(size_t M, size_t N, mat_t &A, mat2_t &Q, vec_t &S) {

    //if (A->size1<A->size2)
    //else if (Q->size1 != A->size2)
    //else if (Q->size1 != Q->size2)
    //else if (S->size != A->size2)
    //const size_t M=A->size1;
    //const size_t N=A->size2;

    // [GSL] Initialize the rotation counter and the sweep counter. 
    int count=1;
    int sweep=0;
    int sweepmax=5*N;

      double dbl_eps=std::numeric_limits<double>::epsilon();

    double tolerance=10*M*dbl_eps;

    // [GSL] Always do at least 12 sweeps. 
    if (sweepmax<12) sweepmax=12;

    // [GSL] Set Q to the identity matrix. 
    for(size_t i=0;i<N;i++) {
      for(size_t j=0;j<N;j++) {
	if (i==j) O2SCL_IX2(Q,i,j)=1.0;
	else O2SCL_IX2(Q,i,j)=0.0;
      }
    }

    // [GSL] Store the column error estimates in S, pfor use during the
    // orthogonalization 

    for (size_t j=0;j<N;j++) {
      double sj=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,0,j,M);
      O2SCL_IX(S,j)=dbl_eps*sj;
    }
    
    // [GSL] Orthogonalize A by plane rotations. 

    while (count > 0 && sweep <= sweepmax) {

      // [GSL] Initialize rotation counter. 
      count=N*(N-1)/2;

      for (size_t j=0;j<N-1;j++) {
	for (size_t k=j+1;k<N;k++) {

	  double a=0.0;
	  double b=0.0;
	  double p=0.0;
	  double q=0.0;
	  double cosine, sine;
	  double v;
	  double abserr_a, abserr_b;
	  int sorted, orthog, noisya, noisyb;

	  //gsl_blas_ddot(&cj.vector,&ck.vector,&p);
	  for(size_t ii=0;ii<M;ii++) {
	    p+=O2SCL_IX2(A,ii,j)*O2SCL_IX2(A,ii,k);
	  }
	  // [GSL] equation 9a: p = 2 x.y 
	  p*=2.0; 

	  a=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,0,j,M);
	  b=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,0,k,M);

	  q=a*a-b*b;
	  v=hypot(p,q);

	  // [GSL] test for columns j,k orthogonal, or dominant errors 
	  
	  abserr_a=O2SCL_IX(S,j);
	  abserr_b=O2SCL_IX(S,k);

	  //sorted=(GSL_COERCE_DBL(a) >= GSL_COERCE_DBL(b));
	  //orthog=(fabs (p) <= tolerance*GSL_COERCE_DBL(a*b));
	  sorted=(a>=b);
	  orthog=(fabs (p) <= tolerance*a*b);
	  noisya=(a<abserr_a);
	  noisyb=(b<abserr_b);

	  if (sorted && (orthog || noisya || noisyb)) {
	    count--;
	    continue;
	  }

	  // [GSL] calculate rotation angles 
	  if (v == 0 || !sorted) {
	    cosine=0.0;
	    sine=1.0;
	  } else {
	    cosine=sqrt((v+q)/(2.0*v));
	    sine=p/(2.0*v*cosine);
	  }

	  // [GSL] apply rotation to A 
	  for (size_t i=0;i<M;i++) {
	    const double Aij=O2SCL_IX2(A,i,j);
	    const double Aik=O2SCL_IX2(A,i,k);
	    O2SCL_IX2(A,i,j)=Aij*cosine+Aik*sine;
	    O2SCL_IX2(A,i,k)=-Aij*sine+Aik*cosine;
	  }
	  
	  O2SCL_IX(S,j)=fabs(cosine)*abserr_a+fabs(sine)*abserr_b;
	  O2SCL_IX(S,k)=fabs(sine)*abserr_a+fabs(cosine)*abserr_b;

	  // [GSL] apply rotation to Q 
	  for (size_t i=0;i<N;i++) {
	    const double Qij=O2SCL_IX2(Q,i,j);
	    const double Qik=O2SCL_IX2(Q,i,k);
	    O2SCL_IX2(Q,i,j)=Qij*cosine+Qik*sine;
	    O2SCL_IX2(Q,i,k)=-Qij*sine+Qik*cosine;
	  }
	}
      }

      // [GSL] Sweep completed. 
      sweep++;
    }

    // [GSL] Orthogonalization complete. Compute singular values.

    {
      double prev_norm=-1.0;

      for (size_t j=0;j<N;j++) {

	double norm=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,0,j,M);

	// [GSL] Determine if singular value is zero, according to the
	// criteria used in the main loop above (i.e. comparison
	// with norm of previous column). 

	if (norm == 0.0 || prev_norm == 0.0 
	    || (j > 0 && norm <= tolerance*prev_norm)) {

	  // [GSL] singular 
	  O2SCL_IX(S,j)=0.0;
	  // [GSL] annihilate column 
	  for(size_t ii=0;ii<M;ii++) {
	    O2SCL_IX2(A,ii,j)=0.0;
	  }

	  prev_norm=0.0;

	} else {

	  // [GSL] non-singular 
	  O2SCL_IX(S,j)=norm;
	  // [GSL] normalize column 
	  O2SCL_CBLAS_NAMESPACE::dscal_subcol(A,0,j,M,1.0/norm);

	  prev_norm=norm;

	}
      }
    }

    if (count > 0) {
      // [GSL] reached sweep limit 
      O2SCL_ERR("Jacobi iterations did not reach desired tolerance",
		o2scl::exc_etol);
    }

    return;
  }

  /** \brief  Balance a general matrix A by scaling the columns
      by the diagonal matrix D
  */
  template<class mat_t, class vec_t>
    void balance_columns(size_t M, size_t N, mat_t &A, vec_t &D) {

    for(size_t j=0;j<N;j++) O2SCL_IX(D,j)=1.0;

    for(size_t j=0;j<N;j++) {
      double s=O2SCL_CBLAS_NAMESPACE::dasum_subcol(A,0,j,M);
      double f=1.0;
      if (s==0.0 || !o2scl::is_finite(s)) {
	O2SCL_IX(D,j)=f;
	continue;
      }
      while (s>1.0) {
	s/=2.0;
	f*=2.0;
      }
      while (s<0.5) {
	s*=2.0;
	f/=2.0;
      }
      O2SCL_IX(D,j)=f;
      if (f!=1.0) {
	O2SCL_CBLAS_NAMESPACE::dscal_subcol(A,0,j,M,1.0/f);
      }
    }

    return;
  }

#ifdef DOXYGEN
}
#endif
