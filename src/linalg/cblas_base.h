/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
/*
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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
/** \file cblas_base.h
    \brief O2scl basic linear algebra function templates

    See \ref o2scl_cblas for more documentation on
    these functions.
    
    \future Add float and complex versions?
    \future There are some Level-1 BLAS functions which are already
    present in vector.h. Should we move all of them to one place or
    the other? The ones in vector.h are generic in the sense that they
    can use doubles or floats, but the ones here can use either () or
    [].

*/

#ifdef DOXYGEN
namespace o2scl_cblas {
#endif
  
  /// Matrix order, either column-major or row-major
  enum o2cblas_order {o2cblas_RowMajor=101, o2cblas_ColMajor=102};

  /// Transpose operations
  enum o2cblas_transpose {o2cblas_NoTrans=111, o2cblas_Trans=112, 
			  o2cblas_ConjTrans=113};

  /// Upper- or lower-triangular
  enum o2cblas_uplo {o2cblas_Upper=121, o2cblas_Lower=122};

  /// Unit or generic diagonal
  enum o2cblas_diag {o2cblas_NonUnit=131, o2cblas_Unit=132};

  /// Left or right sided operation
  enum o2cblas_side {o2cblas_Left=141, o2cblas_Right=142};

  /// \name Test if matrix elements are finite
  //@{
  /** \brief Test if the first \c n elements of a matrix are finite

      If \c n is zero, this will return true without throwing
      an exception.
  */
  template<class mat_t>
  bool matrix_is_finite(size_t m, size_t n, mat_t &data) {
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (!std::isfinite(O2SCL_IX2(data,i,j))) return false;
      }
    }
    return true;
  }

  /** \brief Test if a matrix is finite

      If \c n is zero, this will return true without throwing
      an exception.
  */
  template<class mat_t> bool matrix_is_finite(mat_t &data) {
    return matrix_is_finite(data.size1(),data.size2(),data);
  }
  //@}

  /// \name Level-1 BLAS functions
  //@{
  /** \brief Compute the absolute sum of vector elements

      If \c alpha is zero, this function returns and performs
      no computations. 
  */
  template<class vec_t> double dasum(const size_t N, const vec_t &X) {
    double r=0.0;
    for(size_t i=0;i<N;i++) {
      r+=fabs(X[i]);
    }
    return r;
  }

  /** \brief Compute \f$ y= \alpha x+y \f$

      If \c alpha is zero, this function returns and performs
      no computations. 
  */
  template<class vec_t, class vec2_t>
  void daxpy(const double alpha, const size_t N, const vec_t &X, 
	     vec2_t &Y) {
    
    size_t i;
    
    if (alpha == 0.0) {
      return;
    }
    
    const size_t m=N % 4;
    
    for (i=0;i<m;i++) {
      O2SCL_IX(Y,i)+=alpha*O2SCL_IX(X,i);
    }
    
    for (i=m;i+3<N;i+=4) {
      O2SCL_IX(Y,i)+=alpha*O2SCL_IX(X,i);
      O2SCL_IX(Y,i+1)+=alpha*O2SCL_IX(X,i+1);
      O2SCL_IX(Y,i+2)+=alpha*O2SCL_IX(X,i+2);
      O2SCL_IX(Y,i+3)+=alpha*O2SCL_IX(X,i+3);
    }
  }
  
  /// Compute \f$ r=x \cdot y \f$
  template<class vec_t, class vec2_t> 
  double ddot(const size_t N, const vec_t &X, const vec2_t &Y) {

    double r=0.0;
    size_t i;

    const size_t m=N % 4;
      
    for (i=0;i<m;i++) {
      r+=O2SCL_IX(X,i)*O2SCL_IX(Y,i);
    }
      
    for (i=m;i+3<N;i+=4) {
      r+=O2SCL_IX(X,i)*O2SCL_IX(Y,i);
      r+=O2SCL_IX(X,i+1)*O2SCL_IX(Y,i+1);
      r+=O2SCL_IX(X,i+2)*O2SCL_IX(Y,i+2);
      r+=O2SCL_IX(X,i+3)*O2SCL_IX(Y,i+3);
    }

    return r;
  }
  
  /** \brief Compute the norm of the vector \c X
      
      \note The suffix "2" on the function name indicates that this
      computes the "2-norm", not that the norm is squared.
      
      If \c N is less than or equal to zero, this function returns
      zero without calling the error handler.

      This function works only with vectors which hold \c double. For
      the norm of a general floating point vector, see \ref
      vector_norm().
  */
  template<class vec_t> double dnrm2(const size_t N, const vec_t &X) {
    
    double scale=0.0;
    double ssq=1.0;
    size_t i;
    
    if (N == 0) {
      return 0;
    } else if (N == 1) {
      return fabs(O2SCL_IX(X,0));
    }
    
    for (i=0;i<N;i++) {
      const double x=O2SCL_IX(X,i);
      
      if (x != 0.0) {
	const double ax=fabs(x);
	
	if (scale<ax) {
	  ssq=1.0+ssq*(scale/ax)*(scale/ax);
	  scale=ax;
	} else {
	  ssq+=(ax/scale)*(ax/scale);
	}
      }
      
    }
    
    return scale*sqrt(ssq);
  }

  /** \brief Compute \f$ x=\alpha x \f$
   */
  template<class vec_t> 
  void dscal(const double alpha, const size_t N, vec_t &X) {

    size_t i;
    const size_t m=N % 4;

    for (i=0;i<m;i++) {
      O2SCL_IX(X,i)*=alpha;
    }
    
    for (i=m;i+3<N;i+=4) {
      O2SCL_IX(X,i)*=alpha;
      O2SCL_IX(X,i+1)*=alpha;
      O2SCL_IX(X,i+2)*=alpha;
      O2SCL_IX(X,i+3)*=alpha;
    }
  }
  //@}

  /// \name Level-2 BLAS functions
  //@{
  /** \brief Compute \f$ y=\alpha \left[\mathrm{op}(A)\right] x+
      \beta y \f$.

      If \c M or \c N is zero, or if \c alpha is zero and \c beta is
      one, this function performs no calculations and returns without
      calling the error handler.
  */
  template<class mat_t, class vec_t, class vec2_t>
  void dgemv(const enum o2cblas_order order, 
	     const enum o2cblas_transpose TransA, const size_t M, 
	     const size_t N, const double alpha, const mat_t &A,
	     const vec_t &X, const double beta, vec2_t &Y) {
    
    size_t i, j;
    size_t lenX, lenY;
    
    // If conjugate transpose is requested, just assume plain transpose
    const int Trans=(TransA != o2cblas_ConjTrans) ? TransA : o2cblas_Trans;
    
    if (M == 0 || N == 0) {
      return;
    }

    if (alpha == 0.0 && beta == 1.0) {
      return;
    }

    if (Trans == o2cblas_NoTrans) {
      lenX=N;
      lenY=M;
    } else {
      lenX=M;
      lenY=N;
    }

    /* form  y := beta*y */
    if (beta == 0.0) {
      size_t iy=0;
      for (i=0;i<lenY;i++) {
	O2SCL_IX(Y,iy)=0.0;
	iy++;
      }
    } else if (beta != 1.0) {
      size_t iy=0;
      for (i=0;i<lenY;i++) {
	O2SCL_IX(Y,iy) *= beta;
	iy++;
      }
    }

    if (alpha == 0.0) {
      return;
    }

    if ((order == o2cblas_RowMajor && Trans == o2cblas_NoTrans) ||
	(order == o2cblas_ColMajor && Trans == o2cblas_Trans)) {

      /* form  y := alpha*A*x+y */
      size_t iy=0;
      for (i=0;i<lenY;i++) {
	double temp=0.0;
	size_t ix=0;
	for (j=0;j<lenX;j++) {
	  temp+=O2SCL_IX(X,ix)*O2SCL_IX2(A,i,j);
	  ix++;
	}
	O2SCL_IX(Y,iy)+=alpha*temp;
	iy++;
      }

    } else if ((order == o2cblas_RowMajor && Trans == o2cblas_Trans) ||
	       (order == o2cblas_ColMajor && Trans == o2cblas_NoTrans)) {

      /* form  y := alpha*A'*x+y */
      size_t ix=0;
      for (j=0;j<lenX;j++) {
	const double temp=alpha*O2SCL_IX(X,ix);
	if (temp != 0.0) {
	  size_t iy=0;
	  for (i=0;i<lenY;i++) {
	    O2SCL_IX(Y,iy)+=temp*O2SCL_IX2(A,j,i);
	    iy++;
	  }
	}
	ix++;
      }

    } else {
      O2SCL_ERR("Unrecognized operation in dgemv().",o2scl::exc_einval);
    }
    return;
  }
  
  /** \brief Compute \f$ x=\mathrm{op} (A)^{-1} x \f$
      
      If \c N is zero, this function does nothing and returns zero.
  */
  template<class mat_t, class vec_t> 
  void dtrsv(const enum o2cblas_order order, 
	     const enum o2cblas_uplo Uplo,
	     const enum o2cblas_transpose TransA, 
	     const enum o2cblas_diag Diag,
	     const size_t M, const size_t N, const mat_t &A, vec_t &X) {

    const int nonunit=(Diag == o2cblas_NonUnit);
    int ix, jx;
    int i, j;
    const int Trans=(TransA != o2cblas_ConjTrans) ? TransA : o2cblas_Trans;

    if (N == 0) {
      return;
    }

    /* form  x := inv( A )*x */
    
    if ((order == o2cblas_RowMajor && Trans == o2cblas_NoTrans && 
	 Uplo ==  o2cblas_Upper) ||
	(order == o2cblas_ColMajor && Trans == o2cblas_Trans && 
	 Uplo == o2cblas_Lower)) {
      
      /* backsubstitution */

      // O2scl: Note that subtraction of 1 from size_t N here is ok 
      // since we already handled the N=0 case above
      ix=(int)(N-1);
      if (nonunit) {
	O2SCL_IX(X,ix)=O2SCL_IX(X,ix)/O2SCL_IX2(A,N-1,N-1);
      }
      ix--;

      for (i=(int)(N-1);i>0 && i--;) {
	double tmp=O2SCL_IX(X,ix);
	jx=ix +1;
	for (j=i+1;j<((int)N);j++) {
	  const double Aij=O2SCL_IX2(A,i,j);
	  tmp-=Aij*O2SCL_IX(X,jx);
	  jx++;
	}
	if (nonunit) {
	  O2SCL_IX(X,ix)=tmp/O2SCL_IX2(A,i,i);
	} else {
	  O2SCL_IX(X,ix)=tmp;
	}
	ix--;
      }

    } else if ((order == o2cblas_RowMajor && Trans == o2cblas_NoTrans && 
		Uplo == o2cblas_Lower) || 
	       (order == o2cblas_ColMajor && Trans == o2cblas_Trans && 
		Uplo == o2cblas_Upper)) {
	
      /* forward substitution */
      ix=0;
      if (nonunit) {
	O2SCL_IX(X,ix)=O2SCL_IX(X,ix)/O2SCL_IX2(A,0,0);
      }
      ix++;
      for (i=1;i<((int)N);i++) {
	double tmp=O2SCL_IX(X,ix);
	jx=0;
	for (j=0;j<i;j++) {
	  const double Aij=O2SCL_IX2(A,i,j);
	  tmp-=Aij*O2SCL_IX(X,jx);
	  jx++;
	}
	if (nonunit) {
	  O2SCL_IX(X,ix)=tmp/O2SCL_IX2(A,i,i);
	} else {
	  O2SCL_IX(X,ix)=tmp;
	}
	ix++;
      }

    } else if ((order == o2cblas_RowMajor && Trans == o2cblas_Trans && 
		Uplo == o2cblas_Upper) ||
	       (order == o2cblas_ColMajor && Trans == o2cblas_NoTrans && 
		Uplo == o2cblas_Lower)) {
	
      /* form  x := inv( A' )*x */
	
      /* forward substitution */
      ix=0;
      if (nonunit) {
	O2SCL_IX(X,ix)=O2SCL_IX(X,ix)/O2SCL_IX2(A,0,0);
      }
      ix++;
      for (i=1;i<((int)N);i++) {
	double tmp=O2SCL_IX(X,ix);
	jx=0;
	for (j=0;j<i;j++) {
	  const double Aji=O2SCL_IX2(A,j,i);
	  tmp-=Aji*O2SCL_IX(X,jx);
	  jx++;
	}
	if (nonunit) {
	  O2SCL_IX(X,ix)=tmp/O2SCL_IX2(A,i,i);
	} else {
	  O2SCL_IX(X,ix)=tmp;
	}
	ix++;
      }

    } else if ((order == o2cblas_RowMajor && Trans == o2cblas_Trans && 
		Uplo == o2cblas_Lower) ||
	       (order == o2cblas_ColMajor && Trans == o2cblas_NoTrans && 
		Uplo == o2cblas_Upper)) {
	
      /* backsubstitution */
      // O2scl: Note that subtraction of 1 from size_t N here is ok 
      // since we already handled the N=0 case above
      ix=(int)(N-1);
      if (nonunit) {
	O2SCL_IX(X,ix)=O2SCL_IX(X,ix)/O2SCL_IX2(A,N-1,N-1);
      }
      ix--;
      for (i=(int)(N-1);i>0 && i--;) {
	double tmp=O2SCL_IX(X,ix);
	jx=ix+1;
	for (j=i+1;j<((int)N);j++) {
	  const double Aji=O2SCL_IX2(A,j,i);
	  tmp-=Aji*O2SCL_IX(X,jx);
	  jx++;
	}
	if (nonunit) {
	  O2SCL_IX(X,ix)=tmp/O2SCL_IX2(A,i,i);
	} else {
	  O2SCL_IX(X,ix)=tmp;
	}
	ix--;
      }

    } else {
      O2SCL_ERR("Unrecognized operation in dtrsv().",o2scl::exc_einval);
    }
    return;
  }
  
  /** \brief Compute \f$ x=op(A) x \f$ for the triangular matrix \c A
   */
  template<class mat_t, class vec_t>
  void dtrmv(const enum o2cblas_order Order, 
	     const enum o2cblas_uplo Uplo,
	     const enum o2cblas_transpose TransA,
	     const enum o2cblas_diag Diag, const size_t N,
	     const mat_t &A, vec_t &x) {

    int i, j;

    const int nonunit=(Diag == o2cblas_NonUnit);
    const int Trans=(TransA != o2cblas_ConjTrans) ? TransA : o2cblas_Trans;
    
    if ((Order == o2cblas_RowMajor && Trans == o2cblas_NoTrans && 
	 Uplo == o2cblas_Upper) || 
	(Order == o2cblas_ColMajor && Trans == o2cblas_Trans && 
	 Uplo == o2cblas_Lower)) {

      /* form  x := A*x */

      for (i=0;i<N;i++) {
	double temp=0.0;
	const size_t j_min=i+1;
	const size_t j_max=N;
	size_t jx=j_min;
	for (j=j_min;j<j_max;j++) {
	  temp+=O2SCL_IX(x,jx)*O2SCL_IX2(A,i,j);
	  jx++;
	}
	if (nonunit) {
	  O2SCL_IX(x,i)=temp+O2SCL_IX(x,i)*O2SCL_IX2(A,i,i);
	} else {
	  O2SCL_IX(x,i)+=temp;
	}
      }

    } else if ((Order == o2cblas_RowMajor && Trans == o2cblas_NoTrans && 
		Uplo == o2cblas_Lower) || 
	       (Order == o2cblas_ColMajor && Trans == o2cblas_Trans && 
		Uplo == o2cblas_Upper)) {
	       
      // O2scl: Note that subtraction of 1 from size_t N here is ok 
      // since we already handled the N=0 case above
      size_t ix=N-1;
      for (i=N;i>0 && i--;) {
	double temp=0.0;
	const size_t j_min=0;
	const size_t j_max=i;
	size_t jx=j_min;
	for (j=j_min;j<j_max;j++) {
	  temp+=O2SCL_IX(x,jx)*O2SCL_IX2(A,i,j);
	  jx++;
	}
	if (nonunit) {
	  O2SCL_IX(x,ix)=temp+O2SCL_IX(x,ix)*O2SCL_IX2(A,i,i);
	} else {
	  O2SCL_IX(x,ix)+=temp;
	}
	ix--;
      }

    } else if ((Order == o2cblas_RowMajor && Trans == o2cblas_Trans && 
		Uplo == o2cblas_Upper) || 
	       (Order == o2cblas_ColMajor && Trans == o2cblas_NoTrans && 
		Uplo == o2cblas_Lower)) {
	       
      /* form  x := A'*x */
      size_t ix=N-1;
      for (i=N;i>0 && i--;) {
	double temp=0.0;
	const size_t j_min=0;
	const size_t j_max=i;
	size_t jx=j_min;
	for (j=j_min;j<j_max;j++) {
	  temp+=O2SCL_IX(x,jx)*O2SCL_IX2(A,j,i);
	  jx++;
	}
	if (nonunit) {
	  O2SCL_IX(x,ix)=temp+O2SCL_IX(x,ix)*O2SCL_IX2(A,i,i);
	} else {
	  O2SCL_IX(x,ix)+=temp;
	}
	ix--;
      }

    } else if ((Order == o2cblas_RowMajor && Trans == o2cblas_Trans && 
		Uplo == o2cblas_Lower) ||
	       (Order == o2cblas_ColMajor && Trans == o2cblas_NoTrans && 
		Uplo == o2cblas_Upper)) { 

      for (i=0;i<N;i++) {
	double temp=0.0;
	const size_t j_min=i+1;
	const size_t j_max=N;
	size_t jx=i+1;
	for (j=j_min;j<j_max;j++) {
	  temp+=O2SCL_IX(x,jx)*O2SCL_IX2(A,j,i);
	  jx++;
	}
	if (nonunit) {
	  O2SCL_IX(x,i)=temp+O2SCL_IX(x,i)*O2SCL_IX2(A,i,i);
	} else {
	  O2SCL_IX(x,i)+=temp;
	}
      }

    } else {
      O2SCL_ERR("Unrecognized operation in dtrmv().",
		o2scl::exc_einval);
    }

    return;
  }
  //@}

  /// \name Level-3 BLAS functions
  //@{
  /** \brief Compute \f$ y=\alpha \mathrm{op}(A) \mathrm{op}(B) +
      \beta C \f$
      
      When both \c TransA and \c TransB are \c NoTrans, this function
      operates on the first M rows and K columns of matrix A, and the
      first K rows and N columns of matrix B to produce a matrix C
      with M rows and N columns. 

      This function works for all values of \c Order, \c TransA, and
      \c TransB.
  */
  template<class mat_t>
  void dgemm(const enum o2cblas_order Order, 
	     const enum o2cblas_transpose TransA,
	     const enum o2cblas_transpose TransB, const size_t M, 
	     const size_t N, const size_t K, const double alpha, 
	     const mat_t &A, const mat_t &B, const double beta, mat_t &C) {
    
    size_t i, j, k;
    size_t n1, n2;
    int TransF, TransG;

    if (alpha == 0.0 && beta == 1.0) {
      return;
    }

    /*
      This is a little more complicated than the original in GSL,
      which assigned the matrices A and B to variables *F and *G which
      then allowed some code duplication. We can't really do that
      here, since we don't have that kind of type info on A and B, so
      we just handle the two cases separately.
    */
    
    if (Order == o2cblas_RowMajor) {

      n1=M;
      n2=N;

      /* form  y := beta*y */
      if (beta == 0.0) {
	for (i=0;i<n1;i++) {
	  for (j=0;j<n2;j++) {
	    O2SCL_IX2(C,i,j)=0.0;
	  }
	}
      } else if (beta != 1.0) {
	for (i=0;i<n1;i++) {
	  for (j=0;j<n2;j++) {
	    O2SCL_IX2(C,i,j)*=beta;
	  }
	}
      }

      if (alpha == 0.0) {
	return;
      }

      TransF=(TransA == o2cblas_ConjTrans) ? o2cblas_Trans : TransA;
      TransG=(TransB == o2cblas_ConjTrans) ? o2cblas_Trans : TransB;
      
      if (TransF == o2cblas_NoTrans && TransG == o2cblas_NoTrans) {

	/* form  C := alpha*A*B+C */

	for (k=0;k<K;k++) {
	  for (i=0;i<n1;i++) {
	    const double temp=alpha*O2SCL_IX2(A,i,k);
	    if (temp != 0.0) {
	      for (j=0;j<n2;j++) {
		O2SCL_IX2(C,i,j)+=temp*O2SCL_IX2(B,k,j);
	      }
	    }
	  }
	}

      } else if (TransF == o2cblas_NoTrans && TransG == o2cblas_Trans) {

	/* form  C := alpha*A*B'+C */

	for (i=0;i<n1;i++) {
	  for (j=0;j<n2;j++) {
	    double temp=0.0;
	    for (k=0;k<K;k++) {
	      temp+=O2SCL_IX2(A,i,k)*O2SCL_IX2(B,j,k);
	    }
	    O2SCL_IX2(C,i,j)+=alpha*temp;
	  }
	}

      } else if (TransF == o2cblas_Trans && TransG == o2cblas_NoTrans) {

	for (k=0;k<K;k++) {
	  for (i=0;i<n1;i++) {
	    const double temp=alpha*O2SCL_IX2(A,k,i);
	    if (temp != 0.0) {
	      for (j=0;j<n2;j++) {
		O2SCL_IX2(C,i,j)+=temp*O2SCL_IX2(B,k,j);
	      }
	    }
	  }
	}

      } else if (TransF == o2cblas_Trans && TransG == o2cblas_Trans) {
	
	for (i=0;i<n1;i++) {
	  for (j=0;j<n2;j++) {
	    double temp=0.0;
	    for (k=0;k<K;k++) {
	      temp+=O2SCL_IX2(A,k,i)*O2SCL_IX2(B,j,k);
	    }
	    O2SCL_IX2(C,i,j)+=alpha*temp;
	  }
	}

      } else {
	O2SCL_ERR("Unrecognized operation in dgemm().",o2scl::exc_einval);
      }

    } else {

      // Column-major case

      n1=N;
      n2=M;

      /* form  y := beta*y */
      if (beta == 0.0) {
	for (i=0;i<n1;i++) {
	  for (j=0;j<n2;j++) {
	    O2SCL_IX2(C,i,j)=0.0;
	  }
	}
      } else if (beta != 1.0) {
	for (i=0;i<n1;i++) {
	  for (j=0;j<n2;j++) {
	    O2SCL_IX2(C,i,j)*=beta;
	  }
	}
      }

      if (alpha == 0.0) {
	return;
      }

      TransF=(TransB == o2cblas_ConjTrans) ? o2cblas_Trans : TransB;
      TransG=(TransA == o2cblas_ConjTrans) ? o2cblas_Trans : TransA;

      if (TransF == o2cblas_NoTrans && TransG == o2cblas_NoTrans) {

	/* form  C := alpha*A*B+C */

	for (k=0;k<K;k++) {
	  for (i=0;i<n1;i++) {
	    const double temp=alpha*O2SCL_IX2(B,i,k);
	    if (temp != 0.0) {
	      for (j=0;j<n2;j++) {
		O2SCL_IX2(C,i,j)+=temp*O2SCL_IX2(A,k,j);
	      }
	    }
	  }
	}

      } else if (TransF == o2cblas_NoTrans && TransG == o2cblas_Trans) {

	/* form  C := alpha*A*B'+C */

	for (i=0;i<n1;i++) {
	  for (j=0;j<n2;j++) {
	    double temp=0.0;
	    for (k=0;k<K;k++) {
	      temp+=O2SCL_IX2(B,i,k)*O2SCL_IX2(A,j,k);
	    }
	    O2SCL_IX2(C,i,j)+=alpha*temp;
	  }
	}

      } else if (TransF == o2cblas_Trans && TransG == o2cblas_NoTrans) {

	for (k=0;k<K;k++) {
	  for (i=0;i<n1;i++) {
	    const double temp=alpha*O2SCL_IX2(B,k,i);
	    if (temp != 0.0) {
	      for (j=0;j<n2;j++) {
		O2SCL_IX2(C,i,j)+=temp*O2SCL_IX2(A,k,j);
	      }
	    }
	  }
	}

      } else if (TransF == o2cblas_Trans && TransG == o2cblas_Trans) {

	for (i=0;i<n1;i++) {
	  for (j=0;j<n2;j++) {
	    double temp=0.0;
	    for (k=0;k<K;k++) {
	      temp+=O2SCL_IX2(B,k,i)*O2SCL_IX2(A,j,k);
	    }
	    O2SCL_IX2(C,i,j)+=alpha*temp;
	  }
	}

      } else {
	O2SCL_ERR("Unrecognized operation in dgemm().",o2scl::exc_einval);
      }
    }

    return;
  }

  /** \brief Compute \f$ B=\alpha \mathrm{op}[\mathrm{inv}(A)] B
      \f$ where $A$ is triangular
      
      This function works for all values of \c Order, \c Side, \c Uplo,
      \c TransA, and \c Diag . The variable \c Side is \c Left
      when A is on the left 

      This function operates on the first M rows and N columns of the
      matrix B. If Side is Left, then this function operates on the
      first M rows and M columns of A. If Side is Right, then this
      function operates on the first N rows and N columns of A.
  */
  template<class mat_t>
  void dtrsm(const enum o2cblas_order Order,
	     const enum o2cblas_side Side, 
	     const enum o2cblas_uplo Uplo, 
	     const enum o2cblas_transpose TransA,
	     const enum o2cblas_diag Diag, 
	     const size_t M, const size_t N,
	     const double alpha, const mat_t &A, mat_t &B) {
    
    size_t i, j, k;
    size_t n1, n2;
    
    const int nonunit = (Diag == o2cblas_NonUnit);
    int side, uplo, trans;
    
    if (Order == o2cblas_RowMajor) {
      n1 = M;
      n2 = N;
      side = Side;
      uplo = Uplo;
      trans = (TransA == o2cblas_ConjTrans) ? o2cblas_Trans : TransA;
    } else {
      n1 = N;
      n2 = M;
      side = (Side == o2cblas_Left) ? o2cblas_Right : o2cblas_Left;
      uplo = (Uplo == o2cblas_Upper) ? o2cblas_Lower : o2cblas_Upper;
      trans = (TransA == o2cblas_ConjTrans) ? o2cblas_Trans : TransA;
    }
    
    if (side == o2cblas_Left && uplo == o2cblas_Upper &&
	trans == o2cblas_NoTrans) {
      
      /* form  B := alpha * inv(TriU(A)) *B */
      
      if (alpha != 1.0) {
	for (i = 0; i < n1; i++) {
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j)*=alpha;
	  }
	}
      }
      
      for (i = n1; i > 0 && i--;) {
	if (nonunit) {
	  double Aii = O2SCL_IX2(A,i,i);
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j)/=Aii;
	  }
	}
	
	for (k = 0; k < i; k++) {
	  const double Aki = O2SCL_IX2(A,k,i);
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,k,j)-=Aki*O2SCL_IX2(B,i,j);
	  }
	}
      }

    } else if (side == o2cblas_Left && uplo == o2cblas_Upper &&
	       trans == o2cblas_Trans) {

      /* form  B := alpha * inv(TriU(A))' *B */

      if (alpha != 1.0) {
	for (i = 0; i < n1; i++) {
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = 0; i < n1; i++) {
	if (nonunit) {
	  double Aii = O2SCL_IX2(A,i,i);
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) /= Aii;
	  }
	}

	for (k = i + 1; k < n1; k++) {
	  const double Aik = O2SCL_IX2(A,i,k);
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,k,j) -= Aik * O2SCL_IX2(B,i,j);
	  }
	}
      }

    } else if (side == o2cblas_Left && uplo == o2cblas_Lower &&
	       trans == o2cblas_NoTrans) {

      /* form  B := alpha * inv(TriL(A))*B */

      if (alpha != 1.0) {
	for (i = 0; i < n1; i++) {
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = 0; i < n1; i++) {
	if (nonunit) {
	  double Aii = O2SCL_IX2(A,i,i);
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) /= Aii;
	  }
	}

	for (k = i + 1; k < n1; k++) {
	  const double Aki = O2SCL_IX2(A,k,i);
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,k,j) -= Aki * O2SCL_IX2(B,i,j);
	  }
	}
      }


    } else if (side == o2cblas_Left && uplo == o2cblas_Lower &&
	       trans == o2cblas_Trans) {

      /* form  B := alpha * TriL(A)' *B */

      if (alpha != 1.0) {
	for (i = 0; i < n1; i++) {
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = n1; i > 0 && i--;) {
	if (nonunit) {
	  double Aii = O2SCL_IX2(A,i,i);
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) /= Aii;
	  }
	}

	for (k = 0; k < i; k++) {
	  const double Aik = O2SCL_IX2(A,i,k);
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,k,j) -= Aik * O2SCL_IX2(B,i,j);
	  }
	}
      }

    } else if (side == o2cblas_Right && uplo == o2cblas_Upper &&
	       trans == o2cblas_NoTrans) {

      /* form  B := alpha * B * inv(TriU(A)) */

      if (alpha != 1.0) {
	for (i = 0; i < n1; i++) {
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = 0; i < n1; i++) {
	for (j = 0; j < n2; j++) {
	  if (nonunit) {
	    double Ajj = O2SCL_IX2(A,j,j);
	    O2SCL_IX2(B,i,j) /= Ajj;
	  }

	  {
	    double Bij = O2SCL_IX2(B,i,j);
	    for (k = j + 1; k < n2; k++) {
	      O2SCL_IX2(B,i,k) -= O2SCL_IX2(A,j,k) * Bij;
	    }
	  }
	}
      }

    } else if (side == o2cblas_Right && uplo == o2cblas_Upper &&
	       trans == o2cblas_Trans) {

      /* form  B := alpha * B * inv(TriU(A))' */

      if (alpha != 1.0) {
	for (i = 0; i < n1; i++) {
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = 0; i < n1; i++) {
	for (j = n2; j > 0 && j--;) {

	  if (nonunit) {
	    double Ajj = O2SCL_IX2(A,j,j);
	    O2SCL_IX2(B,i,j) /= Ajj;
	  }

	  {
	    double Bij = O2SCL_IX2(B,i,j);
	    for (k = 0; k < j; k++) {
	      O2SCL_IX2(B,i,k) -= O2SCL_IX2(A,k,j) * Bij;
	    }
	  }
	}
      }

    } else if (side == o2cblas_Right && uplo == o2cblas_Lower &&
	       trans == o2cblas_NoTrans) {

      /* form  B := alpha * B * inv(TriL(A)) */

      if (alpha != 1.0) {
	for (i = 0; i < n1; i++) {
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = 0; i < n1; i++) {
	for (j = n2; j > 0 && j--;) {

	  if (nonunit) {
	    double Ajj = O2SCL_IX2(A,j,j);
	    O2SCL_IX2(B,i,j) /= Ajj;
	  }

	  {
	    double Bij = O2SCL_IX2(B,i,j);
	    for (k = 0; k < j; k++) {
	      O2SCL_IX2(B,i,k) -= O2SCL_IX2(A,j,k) * Bij;
	    }
	  }
	}
      }
      
    } else if (side == o2cblas_Right && uplo == o2cblas_Lower &&
	       trans == o2cblas_Trans) {

      /* form  B := alpha * B * inv(TriL(A))' */

      
      if (alpha != 1.0) {
	for (i = 0; i < n1; i++) {
	  for (j = 0; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = 0; i < n1; i++) {
	for (j = 0; j < n2; j++) {
	  if (nonunit) {
	    double Ajj = O2SCL_IX2(A,j,j);
	    O2SCL_IX2(B,i,j) /= Ajj;
	  }

	  {
	    double Bij = O2SCL_IX2(B,i,j);
	    for (k = j + 1; k < n2; k++) {
	      O2SCL_IX2(B,i,k) -= O2SCL_IX2(A,k,j) * Bij;
	    }
	  }
	}
      }

    } else {
      O2SCL_ERR("Bad operation in dtrsm().",o2scl::exc_einval);
    }
    
    return;
  }
  //@}
  
  /// \name Helper Level-1 BLAS functions - Subvectors
  //@{
  /** \brief Compute \f$ y=\alpha x+y \f$ beginning with index \c ie 
      and ending with index \c N-1
      
      This function is used in \ref householder_hv().

      If \c alpha is identical with zero or <tt>N==ie</tt>, this
      function will perform no calculations and return without calling
      the error handler.
      
      If <tt>ie</tt> is greater than <tt>N-1</tt> then the error 
      handler will be called if \c O2SCL_NO_RANGE_CHECK is not
      defined.
  */
  template<class vec_t, class vec2_t> 
  void daxpy_subvec(const double alpha, const size_t N, const vec_t &X,
		    vec2_t &Y, const size_t ie) {
    
    size_t i;

    if (alpha == 0.0) return;
#if O2SCL_NO_RANGE_CHECK
#else
    if (ie+1>N) {
      O2SCL_ERR("Invalid index in daxpy_subvec().",o2scl::exc_einval);
    }
#endif

    const size_t m=(N-ie) % 4;
    
    for (i=ie;i<ie+m;i++) {
      O2SCL_IX(Y,i)+=alpha*O2SCL_IX(X,i);
    }

    for (;i+3<N;i+=4) {
      O2SCL_IX(Y,i)+=alpha*O2SCL_IX(X,i);
      O2SCL_IX(Y,i+1)+=alpha*O2SCL_IX(X,i+1);
      O2SCL_IX(Y,i+2)+=alpha*O2SCL_IX(X,i+2);
      O2SCL_IX(Y,i+3)+=alpha*O2SCL_IX(X,i+3);
    }
  }

  /** \brief Compute \f$ r=x \cdot y \f$ beginning with index \c ie and
      ending with index \c N-1
      
      This function is used in \ref householder_hv().

      If <tt>ie</tt> is greater than <tt>N-1</tt> then the error 
      handler will be called if \c O2SCL_NO_RANGE_CHECK is not
      defined.
  */
  template<class vec_t, class vec2_t> 
  double ddot_subvec(const size_t N, const vec_t &X, const vec2_t &Y,
		     const size_t ie) {
    double r=0.0;
    size_t i;

#if O2SCL_NO_RANGE_CHECK
#else
    if (ie+1>N) {
      O2SCL_ERR("Invalid index in ddot_subvec().",o2scl::exc_einval);
    }
#endif
    
    const size_t m=(N-ie) % 4;

    for (i=ie;i<ie+m;i++) {
      r+=O2SCL_IX(X,i)*O2SCL_IX(Y,i);
    }

    for (;i+3<N;i+=4) {
      r+=O2SCL_IX(X,i)*O2SCL_IX(Y,i);
      r+=O2SCL_IX(X,i+1)*O2SCL_IX(Y,i+1);
      r+=O2SCL_IX(X,i+2)*O2SCL_IX(Y,i+2);
      r+=O2SCL_IX(X,i+3)*O2SCL_IX(Y,i+3);
    }

    return r;
  }

  /** \brief Compute the norm of the vector \c X beginning with 
      index \c ie and ending with index \c N-1
      
      Used in \ref householder_transform().
      
      \note The suffix "2" on the function name indicates that this
      computes the "2-norm", not that the norm is squared.
      
      If <tt>ie</tt> is greater than <tt>N-1</tt> then the error 
      handler will be called if \c O2SCL_NO_RANGE_CHECK is not
      defined. 
  */
  template<class vec_t> 
  double dnrm2_subvec(const size_t N, const vec_t &X, const size_t ie) {
    
    double scale=0.0;
    double ssq=1.0;
    
#if O2SCL_NO_RANGE_CHECK
#else
    if (ie+1>N) {
      O2SCL_ERR("Invalid index in dnrm2_subvec().",o2scl::exc_einval);
    }
#endif

    if (ie+1==N) {
      return fabs(O2SCL_IX(X,ie));
    }
    
    for (size_t i=ie;i<N;i++) {
      const double x=O2SCL_IX(X,i);
      
      if (x != 0.0) {
	const double ax=fabs(x);
	
	if (scale<ax) {
	  ssq=1.0+ssq*(scale/ax)*(scale/ax);
	  scale=ax;
	} else {
	  ssq+=(ax/scale)*(ax/scale);
	}
      }
      
    }
    
    return scale*sqrt(ssq);
  }

  /** \brief Compute \f$ x=\alpha x \f$ beginning with index \c ie and
      ending with index \c N-1
      
      This function is used in \ref householder_transform().

      If <tt>ie</tt> is greater than <tt>N-1</tt> then the error 
      handler will be called if \c O2SCL_NO_RANGE_CHECK is not
      defined. 
  */
  template<class vec_t> 
  void dscal_subvec(const double alpha, const size_t N, vec_t &X,
		    const size_t ie) {

#if O2SCL_NO_RANGE_CHECK
#else
    if (ie+1>N) {
      O2SCL_ERR("Invalid index in dscal_subvec().",o2scl::exc_einval);
    }
#endif

    size_t i;
    const size_t m=(N-ie) % 4;

    for (i=ie;i<ie+m;i++) {
      O2SCL_IX(X,i)*=alpha;
    }
    
    for (;i+3<N;i+=4) {
      O2SCL_IX(X,i)*=alpha;
      O2SCL_IX(X,i+1)*=alpha;
      O2SCL_IX(X,i+2)*=alpha;
      O2SCL_IX(X,i+3)*=alpha;
    }
  }
  //@}

  /// \name Helper Level-1 BLAS functions - Subcolums of a matrix
  //@{
  /** \brief Compute \f$ y=\alpha x+y \f$ for a subcolumn of a matrix
      
      Given the matrix \c X, define the vector \c x as the column with
      index \c ic. This function computes \f$ y=\alpha x+y \f$ 
      for elements in the vectors \c x and \c y from row \c ir to
      row \c <tt>M-1</tt> (inclusive). All other elements in \c
      x and \c y are not referenced.
      
      Used in householder_hv_sub().
  */
  template<class mat_t, class vec_t> 
  void daxpy_subcol(const double alpha, const size_t M, const mat_t &X,
		    const size_t ir, const size_t ic, vec_t &y) {
    
#if O2SCL_NO_RANGE_CHECK
#else
    if (ir+1>M) {
      O2SCL_ERR("Invalid index in daxpy_subcol().",o2scl::exc_einval);
    }
#endif

    if (alpha == 0.0) {
      return;
    }

    size_t i;
    const size_t m=(M-ir) % 4;
    
    for (i=ir;i<m+ir;i++) {
      O2SCL_IX(y,i)+=alpha*O2SCL_IX2(X,i,ic);
    }
    
    for (;i+3<M;i+=4) {
      O2SCL_IX(y,i)+=alpha*O2SCL_IX2(X,i,ic);
      O2SCL_IX(y,i+1)+=alpha*O2SCL_IX2(X,i+1,ic);
      O2SCL_IX(y,i+2)+=alpha*O2SCL_IX2(X,i+2,ic);
      O2SCL_IX(y,i+3)+=alpha*O2SCL_IX2(X,i+3,ic);
    }

    return;
  }

  /** \brief Compute \f$ r=x \cdot y \f$ for a subcolumn of a matrix
      
      Given the matrix \c X, define the vector \c x as the column with
      index \c ic. This function computes \f$ r=x \cdot y \f$ 
      for elements in the vectors \c x and \c y from row \c ir to
      row \c <tt>M-1</tt> (inclusive). All other elements in \c
      x and \c y are not referenced.
      
      Used in householder_hv_sub().
  */
  template<class mat_t, class vec_t> 
  double ddot_subcol(const size_t M, const mat_t &X, const size_t ir, 
		     const size_t ic, const vec_t &y) {
#if O2SCL_NO_RANGE_CHECK
#else
    if (ir+1>M) {
      O2SCL_ERR("Invalid index in ddot_subcol().",o2scl::exc_einval);
    }
#endif

    double r=0.0;
    size_t i;
    const size_t m=(M-ir) % 4;
    
    for (i=ir;i<m+ir;i++) {
      r+=O2SCL_IX2(X,i,ic)*O2SCL_IX(y,i);
    }
    
    for (;i+3<M;i+=4) {
      r+=O2SCL_IX2(X,i,ic)*O2SCL_IX(y,i);
      r+=O2SCL_IX2(X,i+1,ic)*O2SCL_IX(y,i+1);
      r+=O2SCL_IX2(X,i+2,ic)*O2SCL_IX(y,i+2);
      r+=O2SCL_IX2(X,i+3,ic)*O2SCL_IX(y,i+3);
    }

    return r;
  }

  /** \brief Compute the norm of a subcolumn of a matrix
      
      Given the matrix \c A, define the vector \c x as the column with
      index \c ic. This function computes the norm of the part of \c x
      from row \c ir to row \c <tt>M-1</tt> (inclusive). All other
      elements in \c x are not referenced.
      
      if \c M is zero, then this function silently returns zero
      without calling the error handler.
      
      This function is used in householder_transform_subcol().
      
      \note The suffix "2" on the function name indicates that
      this computes the "2-norm", not that the norm is squared.
  */
  template<class mat_t> 
  double dnrm2_subcol(const mat_t &A, const size_t ir, const size_t ic,
		      const size_t M) {
    
    double scale=0.0;
    double ssq=1.0;
    size_t i;
    
#if O2SCL_NO_RANGE_CHECK
#else
    if (ir+1>M) {
      O2SCL_ERR("Invalid index in dnrm2_subcol().",o2scl::exc_einval);
    }
#endif

    // Handle the one-element vector case separately
    if (ir+1 == M) {
      return fabs(O2SCL_IX2(A,ir,ic));
    }
    
    for (i=ir;i<M;i++) {
      const double x=O2SCL_IX2(A,i,ic);
      
      if (x != 0.0) {
	const double ax=fabs(x);
	
	if (scale<ax) {
	  ssq=1.0+ssq*(scale/ax)*(scale/ax);
	  scale=ax;
	} else {
	  ssq+=(ax/scale)*(ax/scale);
	}
      }
      
    }
    
    return scale*sqrt(ssq);
  }

  /** \brief Compute \f$ x=\alpha x \f$ for a subcolumn of a matrix

      Given the matrix \c A, define the vector \c x as the column with
      index \c ic. This function computes \f$ x= \alpha x \f$ for
      elements in the vectors \c x from row \c ir to row \c
      <tt>M-1</tt> (inclusive). All other elements in \c x are not
      referenced.
      
      Used in householder_transform_subcol().
  */
  template<class mat_t> 
  void dscal_subcol(mat_t &A, const size_t ir, const size_t ic,
		    const size_t M, const double alpha) {

#if O2SCL_NO_RANGE_CHECK
#else
    if (ir+1>M) {
      O2SCL_ERR("Invalid index in dscal_subcol().",o2scl::exc_einval);
    }
#endif

    size_t i;
    const size_t m=(M-ir) % 4;
    
    for (i=ir;i<m+ir;i++) {
      O2SCL_IX2(A,i,ic)*=alpha;
    }
    
    for (;i+3<M;i+=4) {
      O2SCL_IX2(A,i,ic)*=alpha;
      O2SCL_IX2(A,i+1,ic)*=alpha;
      O2SCL_IX2(A,i+2,ic)*=alpha;
      O2SCL_IX2(A,i+3,ic)*=alpha;
    }

    return;
  }

  /** \brief Compute \f$ x=\alpha x \f$ for a subcolumn of a matrix

      Given the matrix \c A, define the vector \c x as the column with
      index \c ic. This function computes \f$ x= \alpha x \f$ for
      elements in the vectors \c x from row \c ir to row \c
      <tt>M-1</tt> (inclusive). All other elements in \c x are not
      referenced.
      
      Used in householder_transform_subcol().
  */
  template<class mat_t> 
  double dasum_subcol(mat_t &A, const size_t ir, const size_t ic,
		      const size_t M) {
    
#if O2SCL_NO_RANGE_CHECK
#else
    if (ir+1>M) {
      O2SCL_ERR("Invalid index in dscal_subcol().",o2scl::exc_einval);
    }
#endif

    size_t i;
    double r=0.0;
    const size_t m=(M-ir) % 4;
    
    for (i=ir;i<m+ir;i++) {
      r+=fabs(O2SCL_IX2(A,i,ic));
    }
    
    for (;i+3<M;i+=4) {
      r+=fabs(O2SCL_IX2(A,i,ic));
      r+=fabs(O2SCL_IX2(A,i+1,ic));
      r+=fabs(O2SCL_IX2(A,i+2,ic));
      r+=fabs(O2SCL_IX2(A,i+3,ic));
    }

    return r;
  }
  //@}

  /// \name Helper Level-1 BLAS functions - Subrows of a matrix
  //@{
  /** \brief Compute \f$ y=\alpha x+y \f$ for a subrow of a matrix

      Given the matrix \c X, define the vector \c x as the row with
      index \c ir. This function computes \f$ y=\alpha x+y \f$ for
      elements in the vectors \c x from column \c ic to column \c
      <tt>N-1</tt> (inclusive). All other elements in \c x and 
      \c y are not referenced.

      If <tt>ic</tt> is greater than <tt>N-1</tt> then the error 
      handler will be called if \c O2SCL_NO_RANGE_CHECK is not
      defined. 

      Used in householder_hv_sub().
  */
  template<class mat_t, class vec_t> 
  void daxpy_subrow(const double alpha, const size_t N, const mat_t &X,
		    const size_t ir, const size_t ic, vec_t &Y) {
    
#if O2SCL_NO_RANGE_CHECK
#else
    if (ic+1>N) {
      O2SCL_ERR("Invalid index in daxpy_subrow().",o2scl::exc_einval);
    }
#endif

    if (alpha == 0.0) {
      return;
    }

    size_t i;
    const size_t m=(N-ic) % 4;
    
    for (i=ic;i<m+ic;i++) {
      O2SCL_IX(Y,i)+=alpha*O2SCL_IX2(X,ir,i);
    }
    
    for (;i+3<N;i+=4) {
      O2SCL_IX(Y,i)+=alpha*O2SCL_IX2(X,ir,i);
      O2SCL_IX(Y,i+1)+=alpha*O2SCL_IX2(X,ir,i+1);
      O2SCL_IX(Y,i+2)+=alpha*O2SCL_IX2(X,ir,i+2);
      O2SCL_IX(Y,i+3)+=alpha*O2SCL_IX2(X,ir,i+3);
    }

    return;
  }

  /** \brief Compute \f$ r=x \cdot y \f$ for a subrow of a matrix
      
      Given the matrix \c X, define the vector \c x as the row with
      index \c ir. This function computes \f$ r=x \cdot y \f$ for
      elements in the vectors \c x from column \c ic to column \c
      <tt>N-1</tt> (inclusive). All other elements in \c x and 
      \c y are not referenced.

      If <tt>ic</tt> is greater than <tt>N-1</tt> then the error 
      handler will be called if \c O2SCL_NO_RANGE_CHECK is not
      defined. 

      Used in householder_hv_sub().
  */
  template<class mat_t, class vec_t> 
  double ddot_subrow(const size_t N, const mat_t &X, const size_t ir, 
		     const size_t ic, const vec_t &Y) {

#if O2SCL_NO_RANGE_CHECK
#else
    if (ic+1>N) {
      O2SCL_ERR("Invalid index in ddot_subrow().",o2scl::exc_einval);
    }
#endif

    double r=0.0;
    size_t i;
    const size_t m=(N-ic) % 4;
    
    for (i=ic;i<m+ic;i++) {
      r+=O2SCL_IX2(X,ir,i)*O2SCL_IX(Y,i);
    }
    
    for (;i+3<N;i+=4) {
      r+=O2SCL_IX2(X,ir,i)*O2SCL_IX(Y,i);
      r+=O2SCL_IX2(X,ir,i+1)*O2SCL_IX(Y,i+1);
      r+=O2SCL_IX2(X,ir,i+2)*O2SCL_IX(Y,i+2);
      r+=O2SCL_IX2(X,ir,i+3)*O2SCL_IX(Y,i+3);
    }

    return r;
  }

  /** \brief Compute the norm of a subrow of a matrix
      
      Given the matrix \c X, define the vector \c x as the row with
      index \c ir. This function computes the norm of the part of \c x
      from column \c ic to column \c <tt>N-1</tt> (inclusive). All
      other elements in \c x are not referenced.

      \note The suffix "2" on the function name indicates that this
      computes the "2-norm", not that the norm is squared.
  */
  template<class mat_t> 
  double dnrm2_subrow(const mat_t &M, const size_t ir, const size_t ic,
		      const size_t N) {
    
    double scale=0.0;
    double ssq=1.0;
    size_t i;
    
    if (ic+1==N) {
      return fabs(O2SCL_IX2(M,ir,ic));
    }
    
    for (i=ic;i<N;i++) {
      const double x=O2SCL_IX2(M,ir,i);
      
      if (x != 0.0) {
	const double ax=fabs(x);
	
	if (scale<ax) {
	  ssq=1.0+ssq*(scale/ax)*(scale/ax);
	  scale=ax;
	} else {
	  ssq+=(ax/scale)*(ax/scale);
	}
      }
      
    }
    
    return scale*sqrt(ssq);
  }

  /** \brief Compute \f$ x=\alpha x \f$ for a subrow of a matrix

      Given the matrix \c A, define the vector \c x as the row with
      index \c ir. This function computes \f$ x = \alpha x \f$ for
      elements in the vectors \c x from column \c ic to column \c
      <tt>N-1</tt> (inclusive). All other elements in \c x and 
      \c y are not referenced.

      If <tt>ic</tt> is greater than <tt>N-1</tt> then the error 
      handler will be called if \c O2SCL_NO_RANGE_CHECK is not
      defined. 
  */
  template<class mat_t> 
  void dscal_subrow(mat_t &A, const size_t ir, const size_t ic,
		    const size_t N, const double alpha) {

#if O2SCL_NO_RANGE_CHECK
#else
    if (ic+1>N) {
      O2SCL_ERR("Invalid index in dscal_subrow().",o2scl::exc_einval);
    }
#endif

    size_t i;
    const size_t m=(N-ic) % 4;
    
    for (i=ic;i<m+ic;i++) {
      O2SCL_IX2(A,ir,i)*=alpha;
    }
    
    for (;i+3<N;i+=4) {
      O2SCL_IX2(A,ir,i)*=alpha;
      O2SCL_IX2(A,ir,i+1)*=alpha;
      O2SCL_IX2(A,ir,i+2)*=alpha;
      O2SCL_IX2(A,ir,i+3)*=alpha;
    }
    
    return;
  }
  //@}

#ifdef O2SCL_NEVER_DEFINED
  /// \name Helper Level-3 BLAS functions
  //@{
  /** \brief Compute \f$ y=\alpha \mathrm{op}(A) \mathrm{op}(B) +
      \beta C \f$ using only part of the matrices A, B, and C

      When both \c TransA and \c TransB are \c NoTrans, this function
      operates on the rows in \f$ [\mathrm{rstarta},M-1] \f$ and
      columns \f$ [\mathrm{cstarta,K-1] \f$ in matrix A and rows in
      \f$ [\mathrm{rstartb},K-1] \f$ and columns \f$
      [\mathrm{cstartb,N-1] \f$ in matrix B. The results are placed in
      rows \f$ [\mathrm{rstartc},M-1] \f$ and columns \f$
      [\mathrm{cstartc},N-1] \f$ .
      
      This function works for all values of \c Order, \c TransA, and
      \c TransB.
  */
  template<class mat_t>
  void dgemm_submat(const enum o2cblas_order Order, 
		    const enum o2cblas_transpose TransA,
		    const enum o2cblas_transpose TransB, const size_t M, 
		    const size_t N, const size_t K, const double alpha, 
		    const mat_t &A, const mat_t &B,
		    const double beta, mat_t &C, size_t rstarta,
		    size_t cstarta, size_t rstartb, size_t cstartb,
		    size_t rstartc, size_t cstartc) {
    
    size_t i, j, k;
    size_t n1, n2;
    int TransF, TransG;

    if (alpha == 0.0 && beta == 1.0) {
      return;
    }

    /*
      This is a little more complicated than the original in GSL,
      which assigned the matrices A and B to variables *F and *G which
      then allowed some code duplication. We can't really do that
      here, since we don't have that kind of type info on A and B, so
      we just handle the two cases separately.
    */
    
    if (Order == o2cblas_RowMajor) {

      n1=M;
      n2=N;

      /* form  y := beta*y */
      if (beta == 0.0) {
	for (i=rstartc;i<n1;i++) {
	  for (j=cstartc;j<n2;j++) {
	    O2SCL_IX2(C,i,j)=0.0;
	  }
	}
      } else if (beta != 1.0) {
	for (i=rstartc;i<n1;i++) {
	  for (j=cstartc;j<n2;j++) {
	    O2SCL_IX2(C,i,j)*=beta;
	  }
	}
      }

      if (alpha == 0.0) {
	return;
      }

      TransF=(TransA == o2cblas_ConjTrans) ? o2cblas_Trans : TransA;
      TransG=(TransB == o2cblas_ConjTrans) ? o2cblas_Trans : TransB;
      
      if (TransF == o2cblas_NoTrans && TransG == o2cblas_NoTrans) {

	/* form  C := alpha*A*B+C */

	for (k=cstarta;k<K;k++) {
	  for (i=rstarta;i<n1;i++) {
	    const double temp=alpha*O2SCL_IX2(A,i,k);
	    if (temp != 0.0) {
	      for (j=rstartc;j<n2;j++) {
		O2SCL_IX2(C,i,j)+=temp*O2SCL_IX2(B,k,j);
	      }
	    }
	  }
	}

      } else if (TransF == o2cblas_NoTrans && TransG == o2cblas_Trans) {

	/* form  C := alpha*A*B'+C */

	for (i=mstart;i<n1;i++) {
	  for (j=nstart;j<n2;j++) {
	    double temp=0.0;
	    for (k=kstart;k<K;k++) {
	      temp+=O2SCL_IX2(A,i,k)*O2SCL_IX2(B,j,k);
	    }
	    O2SCL_IX2(C,i,j)+=alpha*temp;
	  }
	}

      } else if (TransF == o2cblas_Trans && TransG == o2cblas_NoTrans) {

	for (k=kstart;k<K;k++) {
	  for (i=mstart;i<n1;i++) {
	    const double temp=alpha*O2SCL_IX2(A,k,i);
	    if (temp != 0.0) {
	      for (j=nstart;j<n2;j++) {
		O2SCL_IX2(C,i,j)+=temp*O2SCL_IX2(B,k,j);
	      }
	    }
	  }
	}

      } else if (TransF == o2cblas_Trans && TransG == o2cblas_Trans) {
	
	for (i=mstart;i<n1;i++) {
	  for (j=nstart;j<n2;j++) {
	    double temp=0.0;
	    for (k=kstart;k<K;k++) {
	      temp+=O2SCL_IX2(A,k,i)*O2SCL_IX2(B,j,k);
	    }
	    O2SCL_IX2(C,i,j)+=alpha*temp;
	  }
	}

      } else {
	O2SCL_ERR("Unrecognized operation in dgemm_submat().",
		  o2scl::exc_einval);
      }

    } else {

      // Column-major case

      n1=N;
      n2=M;

      /* form  y := beta*y */
      if (beta == 0.0) {
	for (i=nstart;i<n1;i++) {
	  for (j=mstart;j<n2;j++) {
	    O2SCL_IX2(C,i,j)=0.0;
	  }
	}
      } else if (beta != 1.0) {
	for (i=nstart;i<n1;i++) {
	  for (j=mstart;j<n2;j++) {
	    O2SCL_IX2(C,i,j)*=beta;
	  }
	}
      }

      if (alpha == 0.0) {
	return;
      }

      TransF=(TransB == o2cblas_ConjTrans) ? o2cblas_Trans : TransB;
      TransG=(TransA == o2cblas_ConjTrans) ? o2cblas_Trans : TransA;

      if (TransF == o2cblas_NoTrans && TransG == o2cblas_NoTrans) {

	/* form  C := alpha*A*B+C */

	for (k=kstart;k<K;k++) {
	  for (i=nstart;i<n1;i++) {
	    const double temp=alpha*O2SCL_IX2(B,i,k);
	    if (temp != 0.0) {
	      for (j=mstart;j<n2;j++) {
		O2SCL_IX2(C,i,j)+=temp*O2SCL_IX2(A,k,j);
	      }
	    }
	  }
	}

      } else if (TransF == o2cblas_NoTrans && TransG == o2cblas_Trans) {

	/* form  C := alpha*A*B'+C */

	for (i=nstart;i<n1;i++) {
	  for (j=mstart;j<n2;j++) {
	    double temp=0.0;
	    for (k=kstart;k<K;k++) {
	      temp+=O2SCL_IX2(B,i,k)*O2SCL_IX2(A,j,k);
	    }
	    O2SCL_IX2(C,i,j)+=alpha*temp;
	  }
	}

      } else if (TransF == o2cblas_Trans && TransG == o2cblas_NoTrans) {

	for (k=kstart;k<K;k++) {
	  for (i=nstart;i<n1;i++) {
	    const double temp=alpha*O2SCL_IX2(B,k,i);
	    if (temp != 0.0) {
	      for (j=mstart;j<n2;j++) {
		O2SCL_IX2(C,i,j)+=temp*O2SCL_IX2(A,k,j);
	      }
	    }
	  }
	}

      } else if (TransF == o2cblas_Trans && TransG == o2cblas_Trans) {

	for (i=nstart;i<n1;i++) {
	  for (j=mstart;j<n2;j++) {
	    double temp=0.0;
	    for (k=kstart;k<K;k++) {
	      temp+=O2SCL_IX2(B,k,i)*O2SCL_IX2(A,j,k);
	    }
	    O2SCL_IX2(C,i,j)+=alpha*temp;
	  }
	}

      } else {
	O2SCL_ERR("Unrecognized operation in dgemm_submat().",
		  o2scl::exc_einval);
      }
    }

    return;
  }

  /** \brief Compute \f$ B=\alpha \mathrm{op}[\mathrm{inv}(A)] B
      \f$ where $A$ is triangular using only part of the matrices
      A and B
      
      This function works for all values of \c Order, \c Side, \c Uplo,
      \c TransA, and \c Diag .
  */
  template<class mat_t>
  void dtrsm_submat(const enum o2cblas_order Order,
		    const enum o2cblas_side Side, 
		    const enum o2cblas_uplo Uplo, 
		    const enum o2cblas_transpose TransA,
		    const enum o2cblas_diag Diag, 
		    const size_t M, const size_t N, const double alpha, 
		    const mat_t &A, mat_t &B, size_t mstart, size_t nstart) {
    
    size_t i, j, k;
    size_t n1, n2;
    size_t istart, jstart;
    
    const int nonunit = (Diag == o2cblas_NonUnit);
    int side, uplo, trans;
    
    if (Order == o2cblas_RowMajor) {
      n1 = M;
      n2 = N;
      side = Side;
      uplo = Uplo;
      trans = (TransA == o2cblas_ConjTrans) ? o2cblas_Trans : TransA;
      istart=mstart;
      jstart=nstart;
    } else {
      n1 = N;
      n2 = M;
      side = (Side == o2cblas_Left) ? o2cblas_Right : o2cblas_Left;
      uplo = (Uplo == o2cblas_Upper) ? o2cblas_Lower : o2cblas_Upper;
      trans = (TransA == o2cblas_ConjTrans) ? o2cblas_Trans : TransA;
      istart=nstart;
      jstart=mstart;
    }
    
    if (side == o2cblas_Left && uplo == o2cblas_Upper &&
	trans == o2cblas_NoTrans) {
      
      /* form  B := alpha * inv(TriU(A)) *B */
      
      if (alpha != 1.0) {
	for (i = istart; i < n1; i++) {
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j)*=alpha;
	  }
	}
      }
      
      for (i = n1; i > istart && i--;) {
	if (nonunit) {
	  double Aii = O2SCL_IX2(A,i,i);
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j)/=Aii;
	  }
	}
	
	for (k = 0; k < i; k++) {
	  const double Aki = O2SCL_IX2(A,k,i);
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,k,j)-=Aki*O2SCL_IX2(B,i,j);
	  }
	}
      }

    } else if (side == o2cblas_Left && uplo == o2cblas_Upper &&
	       trans == o2cblas_Trans) {

      /* form  B := alpha * inv(TriU(A))' *B */

      if (alpha != 1.0) {
	for (i = istart; i < n1; i++) {
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = istart; i < n1; i++) {
	if (nonunit) {
	  double Aii = O2SCL_IX2(A,i,i);
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) /= Aii;
	  }
	}

	for (k = i + 1; k < n1; k++) {
	  const double Aik = O2SCL_IX2(A,i,k);
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,k,j) -= Aik * O2SCL_IX2(B,i,j);
	  }
	}
      }

    } else if (side == o2cblas_Left && uplo == o2cblas_Lower &&
	       trans == o2cblas_NoTrans) {

      /* form  B := alpha * inv(TriL(A))*B */

      if (alpha != 1.0) {
	for (i = istart; i < n1; i++) {
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = istart; i < n1; i++) {
	if (nonunit) {
	  double Aii = O2SCL_IX2(A,i,i);
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) /= Aii;
	  }
	}

	for (k = i + 1; k < n1; k++) {
	  const double Aki = O2SCL_IX2(A,k,i);
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,k,j) -= Aki * O2SCL_IX2(B,i,j);
	  }
	}
      }


    } else if (side == o2cblas_Left && uplo == o2cblas_Lower &&
	       trans == o2cblas_Trans) {

      /* form  B := alpha * TriL(A)' *B */

      if (alpha != 1.0) {
	for (i = istart; i < n1; i++) {
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = n1; i > istart && i--;) {
	if (nonunit) {
	  double Aii = O2SCL_IX2(A,i,i);
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) /= Aii;
	  }
	}

	for (k = 0; k < i; k++) {
	  const double Aik = O2SCL_IX2(A,i,k);
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,k,j) -= Aik * O2SCL_IX2(B,i,j);
	  }
	}
      }

    } else if (side == o2cblas_Right && uplo == o2cblas_Upper &&
	       trans == o2cblas_NoTrans) {

      /* form  B := alpha * B * inv(TriU(A)) */

      if (alpha != 1.0) {
	for (i = istart; i < n1; i++) {
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = istart; i < n1; i++) {
	for (j = jstart; j < n2; j++) {
	  if (nonunit) {
	    double Ajj = O2SCL_IX2(A,j,j);
	    O2SCL_IX2(B,i,j) /= Ajj;
	  }

	  {
	    double Bij = O2SCL_IX2(B,i,j);
	    for (k = j + 1; k < n2; k++) {
	      O2SCL_IX2(B,i,k) -= O2SCL_IX2(A,j,k) * Bij;
	    }
	  }
	}
      }

    } else if (side == o2cblas_Right && uplo == o2cblas_Upper &&
	       trans == o2cblas_Trans) {

      /* form  B := alpha * B * inv(TriU(A))' */

      if (alpha != 1.0) {
	for (i = istart; i < n1; i++) {
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = istart; i < n1; i++) {
	for (j = n2; j > jstart && j--;) {

	  if (nonunit) {
	    double Ajj = O2SCL_IX2(A,j,j);
	    O2SCL_IX2(B,i,j) /= Ajj;
	  }

	  {
	    double Bij = O2SCL_IX2(B,i,j);
	    for (k = 0; k < j; k++) {
	      O2SCL_IX2(B,i,k) -= O2SCL_IX2(A,k,j) * Bij;
	    }
	  }
	}
      }

    } else if (side == o2cblas_Right && uplo == o2cblas_Lower &&
	       trans == o2cblas_NoTrans) {

      /* form  B := alpha * B * inv(TriL(A)) */

      if (alpha != 1.0) {
	for (i = istart; i < n1; i++) {
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = istart; i < n1; i++) {
	for (j = n2; j > jstart && j--;) {

	  if (nonunit) {
	    double Ajj = O2SCL_IX2(A,j,j);
	    O2SCL_IX2(B,i,j) /= Ajj;
	  }

	  {
	    double Bij = O2SCL_IX2(B,i,j);
	    for (k = 0; k < j; k++) {
	      O2SCL_IX2(B,i,k) -= O2SCL_IX2(A,j,k) * Bij;
	    }
	  }
	}
      }
      
    } else if (side == o2cblas_Right && uplo == o2cblas_Lower &&
	       trans == o2cblas_Trans) {

      /* form  B := alpha * B * inv(TriL(A))' */

      
      if (alpha != 1.0) {
	for (i = istart; i < n1; i++) {
	  for (j = jstart; j < n2; j++) {
	    O2SCL_IX2(B,i,j) *= alpha;
	  }
	}
      }

      for (i = istart; i < n1; i++) {
	for (j = jstart; j < n2; j++) {
	  if (nonunit) {
	    double Ajj = O2SCL_IX2(A,j,j);
	    O2SCL_IX2(B,i,j) /= Ajj;
	  }

	  {
	    double Bij = O2SCL_IX2(B,i,j);
	    for (k = j + 1; k < n2; k++) {
	      O2SCL_IX2(B,i,k) -= O2SCL_IX2(A,k,j) * Bij;
	    }
	  }
	}
      }

    } else {
      O2SCL_ERR("Bad operation in dtrsm().",o2scl::exc_einval);
    }
    
    return;
  }
  //@}

#endif
  
#ifdef DOXYGEN
}
#endif

