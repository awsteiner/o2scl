/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
/** \file lanczos_base.h
    \brief File defining \ref o2scl_linalg::lanczos
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Lanczos diagonalization
      
      This class approximates the largest eigenvalues of a symmetric
      matrix.

      The vector and matrix types can be any type which provides access
      via \c operator[] or \c operator(), given a suitable vector
      allocation type. 
      
      The tridiagonalization routine was rewritten from the EISPACK
      routines \c imtql1.f (but uses \c std::hypot() instead of \c
      pythag.f). 

      \future The function eigen_tdiag() automatically sorts the
      eigenvalues, which may not be necessary.

      \future Do something better than the simple matrix-vector product.
      For example, use dgemm() and allow user to specify column or
      row-major.

      \future Rework memory allocation to perform as needed.

      \comment 
      We need the o2scl:: prefix in the default template parameters
      below because this class isn't in the o2scl namespace.
      \endcomment
  */
  template<class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class lanczos {
  
  public:
  
  lanczos() {
    td_iter=30;
    td_lasteval=0;
  }
  
  /** \brief Number of iterations for finding the eigenvalues of the
      tridiagonal matrix (default 30)
  */
  size_t td_iter;

  /** \brief The index for the last eigenvalue not determined if
      tridiagonalization fails
  */
  size_t td_lasteval;

  /** \brief Approximate the largest eigenvalues of a symmetric
      matrix \c mat using the Lanczos method

      Given a square matrix \c mat with size \c size, this function
      applies \c n_iter iterations of the Lanczos algorithm to
      produce \c n_iter approximate eigenvalues stored in \c
      eigen. As a by-product, this function also partially
      tridiagonalizes the matrix placing the result in \c diag and
      \c off_diag. Before calling this function, space must have
      already been allocated for \c eigen, \c diag, and \c
      off_diag. All three of these arrays must have at least enough
      space for \c n_iter elements.
	
      Choosing /c n_iter = \c size will produce all of the exact
      eigenvalues and the corresponding tridiagonal matrix, but this
      may be slower than diagonalizing the matrix directly.
  */
  int eigenvalues(size_t size, mat_t &mat, size_t n_iter, 
		  vec_t &eigen, vec_t &diag, vec_t &off_diag) {
    double t;
    bool cont=true;
    size_t i, j, k;

    vec_t v;
    vec_t w;
    vec_t b3;
    vec_t prod;

    v.resize(size);
    w.resize(size);
    b3.resize(size);
    prod.resize(size);

    // Pick a unit vector
    O2SCL_IX(w,0)=1.0;
    for(i=1;i<size;i++) O2SCL_IX(w,i)=0.0;

    for(i=0;i<size;i++) O2SCL_IX(v,i)=0;
    j=0;
    while (cont) {
      if (j!=0) {
	for(i=0;i<size;i++) {
	  t=O2SCL_IX(w,i);
	  O2SCL_IX(w,i)=O2SCL_IX(v,i)/O2SCL_IX(off_diag,j-1);
	  O2SCL_IX(v,i)=-O2SCL_IX(off_diag,j-1)*t;
	}
      }
      product(size,mat,w,prod);
      for(k=0;k<size;k++) O2SCL_IX(v,k)+=O2SCL_IX(prod,k);
      O2SCL_IX(diag,j)=0.0;
      O2SCL_IX(off_diag,j)=0.0;

      for(k=0;k<size;k++) O2SCL_IX(diag,j)+=O2SCL_IX(w,k)*O2SCL_IX(v,k);
      for(k=0;k<size;k++) O2SCL_IX(v,k)-=O2SCL_IX(diag,j)*O2SCL_IX(w,k);
      for(k=0;k<size;k++) O2SCL_IX(off_diag,j)+=O2SCL_IX(v,k)*O2SCL_IX(v,k);
      O2SCL_IX(off_diag,j)=sqrt(O2SCL_IX(off_diag,j));
      j++;
    
      if (j>=n_iter || O2SCL_IX(off_diag,j-1)==0.0) cont=false;
    
      if (j>0) {
	for(k=0;k<size-1;k++) {
	  O2SCL_IX(b3,k+1)=O2SCL_IX(off_diag,k);
	  O2SCL_IX(eigen,k)=O2SCL_IX(diag,k);
	}
	O2SCL_IX(eigen,size-1)=O2SCL_IX(diag,size-1);
	if (eigen_tdiag(j,eigen,b3)!=0) {
	  O2SCL_ERR2("Call to eigen_tdiag() in lanczos::",
			 "eigenvalues() failed.",o2scl::exc_efailed);
	}
      }
    }
  
    return 0;
  }
    
  /** \brief In-place diagonalization of a tridiagonal matrix

      On input, the vectors \c diag and \c off_diag should both be
      vectors of size \c n. The diagonal entries stored in \c diag,
      and the \f$ n-1 \f$ off-diagonal entries should be stored in
      \c off_diag, starting with \c off_diag[1].  The value in \c
      off_diag[0] is unused. The vector \c off_diag is destroyed by
      the computation.

      This uses an implict QL method from the EISPACK routine \c
      imtql1. The value of \c ierr from the original Fortran routine
      is stored in \ref td_lasteval.

  */
  int eigen_tdiag(size_t n, vec_t &diag, vec_t &off_diag) {

    // The variable 'i' is set to zero here because Cygwin complained
    // about uninitialized variables. This is probably ok, but it
    // would be nice to double check that there isn't a problem with
    // setting i=0 here.

    int i=0,j,l,m,mml;
    double b,c,f,g,p,r,s,tst1,tst2;

    if (n==1) return 0;
      
    for(size_t ij=1;ij<n;ij++) {
      O2SCL_IX(off_diag,ij-1)=O2SCL_IX(off_diag,ij);
    }
    O2SCL_IX(off_diag,n-1)=0.0;
  
    bool done=false;
  
    l=1;
    j=0;
  
    while (done==false && l<=((int)n)) {
    
      // Look for small sub-diagonal element
      bool idone=false;
      for(m=l;m<((int)n) && idone==false;m++) {
	tst1=fabs(O2SCL_IX(diag,m-1))+fabs(O2SCL_IX(diag,m));
	tst2=tst1+fabs(O2SCL_IX(off_diag,m-1));
	if (tst2==tst1) {
	  m--;
	  idone=true;
	}
      }
    
      p=O2SCL_IX(diag,l-1);

      if (m!=l && j==((int)td_iter)) {
	  
	// Set error. No convergence after td_iter iterations
	td_lasteval=l;
	O2SCL_ERR("No convergence in lanczos::eigen_tdiag()",
		      o2scl::exc_efailed);
      }

      if (m!=l) {

	j++;

	// Form shift
	g=(O2SCL_IX(diag,l)-p)/(2.0*O2SCL_IX(off_diag,l-1));
	r=std::hypot(g,1.0);
      
	g=O2SCL_IX(diag,m-1)-p+O2SCL_IX(off_diag,l-1)/
	  (g+(g>=0.0 ? fabs(r) : -fabs(r)));
	s=1.0;
	c=1.0;
	p=0.0;
	mml=m-l;
      
	for(int ii=1;ii<=mml;ii++) {

	  i=m-ii;
	  f=s*O2SCL_IX(off_diag,i-1);
	  b=c*O2SCL_IX(off_diag,i-1);
	  r=std::hypot(f,g);
	  O2SCL_IX(off_diag,i)=r;
	
	  if (r==0.0) {

	    // Recover from underflow
	    O2SCL_IX(diag,i)-=p;
	    O2SCL_IX(off_diag,m-1)=0.0;
	    ii=mml+1;

	  } else {

	    s=f/r;
	    c=g/r;
	    g=O2SCL_IX(diag,i)-p;
	    r=(O2SCL_IX(diag,i-1)-g)*s+2.0*c*b;
	    p=s*r;
	    O2SCL_IX(diag,i)=g+p;
	    g=c*r-b;

	  }

	}
      
	O2SCL_IX(diag,l-1)-=p;
	O2SCL_IX(off_diag,l-1)=g;
	O2SCL_IX(off_diag,m-1)=0.0;
      

      } else {

	// Order eigenvalues

	if (l==1) {

	  i=1;
	  O2SCL_IX(diag,i-1)=p;
	
	} else {

	  bool skip=false;
	  for(int ii=2;ii<=l;ii++) {
	    i=l+2-ii;
	    if (p>=O2SCL_IX(diag,i-2)) {
	      ii=l+1;
	      skip=true;
	    } else {
	      O2SCL_IX(diag,i-1)=O2SCL_IX(diag,i-2);
	    }
	  }
	
	  if (skip==false) i=1;
	  O2SCL_IX(diag,i-1)=p;
	}

	j=0;
	l++;
      }
    
    }
  
    return 0;
  }

#ifndef DOXYGEN_INTERNAL

  protected:
      
  /** \brief Simple matrix-vector product
	
      It is assumed that memory is already allocated for \c prod.
  */
  void product(size_t n, mat_t &a, vec_t &w, vec_t &prod) {
    size_t i, j;
    for(i=0;i<n;i++) {
      O2SCL_IX(prod,i)=0.0;
      for(j=0;j<n;j++) {
	O2SCL_IX(prod,i)+=O2SCL_IX2(a,i,j)*O2SCL_IX(w,j);
      }
    }
    return;
  }
    
#endif
  
  };

#ifdef DOXYGEN
}
#endif
