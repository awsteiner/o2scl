/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2010-2024, Andrew W. Steiner

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
#include <iostream>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/svd.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_cblas;
using namespace o2scl_linalg;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

namespace svd_ts {

  void create_givens (const double a, const double b, double *c, double *s) {
    if (b == 0)
      {
	*c=1;
	*s=0;
      }
    else if (fabs (b) > fabs (a))
      {
	double t=-a / b;
	double s1=1.0 / sqrt (1 + t * t);
	*s=s1;
	*c=s1 * t;
      }
    else
      {
	double t=-b / a;
	double c1=1.0 / sqrt (1 + t * t);
	*c=c1;
	*s=c1 * t;
      }
  }
  
  void chop_small_elements (gsl_vector * d, gsl_vector * f) {
    const size_t N=d->size;
    double d_i=gsl_vector_get (d, 0);

    size_t i;

    for (i=0; i < N - 1; i++) {
      double f_i=gsl_vector_get (f, i);
      double d_ip1=gsl_vector_get (d, i + 1);
      
      if (fabs(f_i) < GSL_DBL_EPSILON*(fabs(d_i)+fabs(d_ip1))) {
	gsl_vector_set(f, i, 0.0);
      }
      
      d_i=d_ip1;
    }
    return;
  }

  double trailing_eigenvalue (const gsl_vector * d, const gsl_vector * f) {
    const size_t n=d->size;

    double da=gsl_vector_get (d, n - 2);
    double db=gsl_vector_get (d, n - 1);
    double fa=(n > 2) ? gsl_vector_get (f, n - 3) : 0.0;
    double fb=gsl_vector_get (f, n - 2);

    double mu;

    /* We can compute mu more accurately than using the formula above
       since we know the roots cannot be negative.  This also avoids
       the possibility of NaNs in the formula above.

       The matrix is [ da^2 + fa^2,  da fb      ;
       da fb      , db^2 + fb^2 ]
       and mu is the eigenvalue closest to the bottom right element.
    */
    
    double ta=da * da + fa * fa;
    double tb=db * db + fb * fb;
    double tab=da * fb;
    
    double dt=(ta - tb) / 2.0;
    
    double S=ta + tb;
    double da2=da * da, db2=db * db;
    double fa2=fa * fa, fb2=fb * fb;
    double P=(da2 * db2) + (fa2 * db2) + (fa2 * fb2);
    double D=hypot(dt, tab);
    double r1=S/2 + D;
    
    if (dt >= 0)
      {
        /* tb < ta, choose smaller root */
        mu=(r1 > 0) ?  P / r1 : 0.0;
      }
    else 
      {
        /* tb > ta, choose larger root */
        mu=r1;
      }

    return mu;
  }

  void create_schur (double d0, double f0, double d1, double * c, 
		     double * s) {
    double apq=2.0 * d0 * f0;

    if (d0 == 0 || f0 == 0)
      {
	*c=1.0;
	*s=0.0;
	return;
      }

    /* Check if we need to rescale to avoid underflow/overflow */
    if (fabs(d0) < GSL_SQRT_DBL_MIN || fabs(d0) > GSL_SQRT_DBL_MAX
	|| fabs(f0) < GSL_SQRT_DBL_MIN || fabs(f0) > GSL_SQRT_DBL_MAX
	|| fabs(d1) < GSL_SQRT_DBL_MIN || fabs(d1) > GSL_SQRT_DBL_MAX)
      {
	double scale;
	int d0_exp, f0_exp;
	frexp(d0, &d0_exp);
	frexp(f0, &f0_exp);
	/* Bring |d0*f0| into the range GSL_DBL_MIN to GSL_DBL_MAX */
	scale=ldexp(1.0, -(d0_exp + f0_exp)/4);
	d0 *= scale;
	f0 *= scale;
	d1 *= scale;
	apq=2.0 * d0 * f0;
      }

    if (apq != 0.0)
      {
	double t;
	double tau=(f0*f0 + (d1 + d0)*(d1 - d0)) / apq;
      
	if (tau >= 0.0)
	  {
	    t=1.0/(tau + hypot(1.0, tau));
	  }
	else
	  {
	    t=-1.0/(-tau + hypot(1.0, tau));
	  }

	*c=1.0 / hypot(1.0, t);
	*s=t * (*c);
      }
    else
      {
	*c=1.0;
	*s=0.0;
      }
  }

  void svd2 (gsl_vector * d, gsl_vector * f, gsl_matrix * U, 
	     gsl_matrix * V) {
    size_t i;
    double c, s, a11, a12, a21, a22;

    const size_t M=U->size1;
    const size_t N=V->size1;

    double d0=gsl_vector_get (d, 0);
    double f0=gsl_vector_get (f, 0);
  
    double d1=gsl_vector_get (d, 1);

    if (d0 == 0.0)
      {
	/* Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] */

	create_givens (f0, d1, &c, &s);

	/* compute B <= G^T B X,  where X=[0,1;1,0] */

	gsl_vector_set (d, 0, c * f0 - s * d1);
	gsl_vector_set (f, 0, s * f0 + c * d1);
	gsl_vector_set (d, 1, 0.0);

	/* Compute U <= U G */

	for (i=0; i < M; i++)
	  {
	    double Uip=gsl_matrix_get (U, i, 0);
	    double Uiq=gsl_matrix_get (U, i, 1);
	    gsl_matrix_set (U, i, 0, c * Uip - s * Uiq);
	    gsl_matrix_set (U, i, 1, s * Uip + c * Uiq);
	  }

	/* Compute V <= V X */

	gsl_matrix_swap_columns (V, 0, 1);

	return;
      }
    else if (d1 == 0.0)
      {
	/* Eliminate off-diagonal element in [d0,f0;0,0] */

	create_givens (d0, f0, &c, &s);

	/* compute B <= B G */

	gsl_vector_set (d, 0, d0 * c - f0 * s);
	gsl_vector_set (f, 0, 0.0);

	/* Compute V <= V G */

	for (i=0; i < N; i++)
	  {
	    double Vip=gsl_matrix_get (V, i, 0);
	    double Viq=gsl_matrix_get (V, i, 1);
	    gsl_matrix_set (V, i, 0, c * Vip - s * Viq);
	    gsl_matrix_set (V, i, 1, s * Vip + c * Viq);
	  }

	return;
      }
    else
      {
	/* Make columns orthogonal, A=[d0, f0; 0, d1] * G */

	create_schur (d0, f0, d1, &c, &s);

	/* compute B <= B G */
      
	a11=c * d0 - s * f0;
	a21=- s * d1;
      
	a12=s * d0 + c * f0;
	a22=c * d1;
      
	/* Compute V <= V G */
      
	for (i=0; i < N; i++)
	  {
	    double Vip=gsl_matrix_get (V, i, 0);
	    double Viq=gsl_matrix_get (V, i, 1);
	    gsl_matrix_set (V, i, 0, c * Vip - s * Viq);
	    gsl_matrix_set (V, i, 1, s * Vip + c * Viq);
	  }
      
	/* Eliminate off-diagonal elements, bring column with largest
	   norm to first column */
      
	if (hypot(a11, a21) < hypot(a12,a22))
	  {
	    double t1, t2;

	    /* B <= B X */

	    t1=a11; a11=a12; a12=t1;
	    t2=a21; a21=a22; a22=t2;

	    /* V <= V X */

	    gsl_matrix_swap_columns(V, 0, 1);
	  } 

	create_givens (a11, a21, &c, &s);
      
	/* compute B <= G^T B */
      
	gsl_vector_set (d, 0, c * a11 - s * a21);
	gsl_vector_set (f, 0, c * a12 - s * a22);
	gsl_vector_set (d, 1, s * a12 + c * a22);
      
	/* Compute U <= U G */
      
	for (i=0; i < M; i++)
	  {
	    double Uip=gsl_matrix_get (U, i, 0);
	    double Uiq=gsl_matrix_get (U, i, 1);
	    gsl_matrix_set (U, i, 0, c * Uip - s * Uiq);
	    gsl_matrix_set (U, i, 1, s * Uip + c * Uiq);
	  }

	return;
      }
  }


  void chase_out_intermediate_zero (gsl_vector * d, gsl_vector * f, 
				    gsl_matrix * U, size_t k0) {

    const size_t M=U->size1;
    const size_t n=d->size;
    double c, s;
    double x, y;
    size_t k;

    x=gsl_vector_get (f, k0);
    y=gsl_vector_get (d, k0+1);

    for (k=k0; k < n - 1; k++)
      {
	create_givens (y, -x, &c, &s);
      
	/* Compute U <= U G */

	{
	  size_t i;

	  for (i=0; i < M; i++)
	    {
	      double Uip=gsl_matrix_get (U, i, k0);
	      double Uiq=gsl_matrix_get (U, i, k + 1);
	      //std::cout << "Uip,Uiq: " << Uip << " " << Uiq << std::endl;
	      gsl_matrix_set (U, i, k0, c * Uip - s * Uiq);
	      gsl_matrix_set (U, i, k + 1, s * Uip + c * Uiq);
	    }
	}
      
	/* compute B <= G^T B */
      
	gsl_vector_set (d, k + 1, s * x + c * y);

	if (k == k0)
	  gsl_vector_set (f, k, c * x - s * y );

	if (k < n - 2) 
	  {
	    double z=gsl_vector_get (f, k + 1);
	    gsl_vector_set (f, k + 1, c * z); 

	    x=-s * z ;
	    y=gsl_vector_get (d, k + 2); 
	  }
      }
  }

  void chase_out_trailing_zero (gsl_vector * d, gsl_vector * f, 
				gsl_matrix * V) {

    const size_t N=V->size1;
    const size_t n=d->size;
    double c, s;
    double x, y;
    size_t k;

    x=gsl_vector_get (d, n - 2);
    y=gsl_vector_get (f, n - 2);

    for (k=n - 1; k-- > 0;)
      {
	create_givens (x, y, &c, &s);

	/* Compute V <= V G where G=[c, s ; -s, c] */

	{
	  size_t i;
   
	  for (i=0; i < N; i++)
	    {
	      double Vip=gsl_matrix_get (V, i, k);
	      double Viq=gsl_matrix_get (V, i, n - 1);
	      gsl_matrix_set (V, i, k, c * Vip - s * Viq);
	      gsl_matrix_set (V, i, n - 1, s * Vip + c * Viq);
	    }
	}

	/* compute B <= B G */
      
	gsl_vector_set (d, k, c * x - s * y);

	if (k == n - 2)
	  gsl_vector_set (f, k, s * x + c * y );

	if (k > 0) 
	  {
	    double z=gsl_vector_get (f, k - 1);
	    gsl_vector_set (f, k - 1, c * z); 

	    x=gsl_vector_get (d, k - 1); 
	    y=s * z ;
	  }
      }
  }

  void qrstep (gsl_vector * d, gsl_vector * f, gsl_matrix * U, 
	       gsl_matrix * V) {

    const size_t M=U->size1;
    const size_t N=V->size1;
    const size_t n=d->size;
    double y, z;
    double ak, bk, zk, ap, bp, aq, bq;
    size_t i, k;

    //std::cout << "M,N,n: " << M << " " << N << " " << n << std::endl;

    if (n == 1)
      return;  /* shouldn't happen */

    /* Compute 2x2 svd directly */

    if (n == 2)
      {
	svd2 (d, f, U, V);
	return;
      }

    /* Chase out any zeroes on the diagonal */

    for (i=0; i < n - 1; i++)
      {
	double d_i=gsl_vector_get (d, i);
	//std::cout << "d_i: " << i << " " << n << " "
	//<< d_i << std::endl;
      
	if (d_i == 0.0)
	  {
	    chase_out_intermediate_zero (d, f, U, i);
	    return;
	  }
      }

    /* Chase out any zero at the end of the diagonal */

    {
      double d_nm1=gsl_vector_get (d, n - 1);
      //std::cout << "d_nm1: " << d_nm1 << std::endl;

      if (d_nm1 == 0.0) 
	{
	  chase_out_trailing_zero (d, f, V);
	  return;
	}
    }


    /* Apply QR reduction steps to the diagonal and offdiagonal */

    {
      double d0=gsl_vector_get (d, 0);
      double f0=gsl_vector_get (f, 0);
    
      double d1=gsl_vector_get (d, 1);
      double f1=gsl_vector_get (f, 1);
      //std::cout << "d0,f0,d1,f1: " << d0 << " " << f0 << " " << d1 << " "
      //<< f1 << std::endl;
    
      {
	double mu=trailing_eigenvalue (d, f);
    
	y=d0 * d0 - mu;
	z=d0 * f0;
      }
    
      /* Set up the recurrence for Givens rotations on a bidiagonal matrix */
    
      ak=0;
      bk=0;
    
      ap=d0;
      bp=f0;
    
      aq=d1;
      bq=f1;
    }

    for (k=0; k < n - 1; k++)
      {
	double c, s;
	create_givens (y, z, &c, &s);

	/* Compute V <= V G */

	for (i=0; i < N; i++)
	  {
	    double Vip=gsl_matrix_get (V, i, k);
	    double Viq=gsl_matrix_get (V, i, k + 1);
	    //std::cout << "Vip,Viq: " << Vip << " " << Viq << std::endl;
	    gsl_matrix_set (V, i, k, c * Vip - s * Viq);
	    gsl_matrix_set (V, i, k + 1, s * Vip + c * Viq);
	  }

	/* compute B <= B G */

	{
	  double bk1=c * bk - s * z;

	  double ap1=c * ap - s * bp;
	  double bp1=s * ap + c * bp;
	  double zp1=-s * aq;

	  double aq1=c * aq;

	  if (k > 0)
	    {
	      gsl_vector_set (f, k - 1, bk1);
	    }

	  ak=ap1;
	  bk=bp1;
	  zk=zp1;

	  ap=aq1;

	  if (k < n - 2)
	    {
	      bp=gsl_vector_get (f, k + 1);
	    }
	  else
	    {
	      bp=0.0;
	    }

	  y=ak;
	  z=zk;
	}

	create_givens (y, z, &c, &s);

	/* Compute U <= U G */

	for (i=0; i < M; i++)
	  {
	    double Uip=gsl_matrix_get (U, i, k);
	    double Uiq=gsl_matrix_get (U, i, k + 1);
	    //std::cout << "Uip2,Uiq2: " << Uip << " " << Uiq << std::endl;
	    gsl_matrix_set (U, i, k, c * Uip - s * Uiq);
	    gsl_matrix_set (U, i, k + 1, s * Uip + c * Uiq);
	  }

	/* compute B <= G^T B */
	
	//std::cout << "k,bk,ap2: " << k << " " << bk << " " << ap << std::endl;
	//std::cout << "ak,zk,bp: " << ak << " " << zk << " " 
	// << bp << std::endl;

	{
	  //std::cout << "prod1: " << c*ak << " " << s*zk << std::endl;
	  //std::cout << "prod2: " << c*bk << " " << s*ap << std::endl;
	  //std::cout << "prod3: " << s*bk << " " << c*ap << std::endl;
	  double ak1=c * ak - s * zk;
	  double bk1=c * bk - s * ap;
	  double zk1=-s * bp;

	  double ap1=s * bk + c * ap;
	  double bp1=c * bp;

	  gsl_vector_set (d, k, ak1);

	  ak=ak1;
	  bk=bk1;
	  zk=zk1;

	  ap=ap1;
	  bp=bp1;
	  //std::cout << "c,s: " << c << " " << s << std::endl;
	  //std::cout << "k,bk,ap: " << k << " " << bk << " " << ap << std::endl;

	  if (k < n - 2)
	    {
	      aq=gsl_vector_get (d, k + 2);
	    }
	  else
	    {
	      aq=0.0;
	    }

	  y=bk;
	  z=zk;
	}
      }

    gsl_vector_set (f, n - 2, bk);
    gsl_vector_set (d, n - 1, ap);
    //std::cout << "bk,ap: " << bk << " " << ap << std::endl;
  }

  int
  gsl_linalg_SV_decomp (gsl_matrix * A, gsl_matrix * V, gsl_vector * S, 
			gsl_vector * work)
  {
    size_t a, b, i, j, iter;

    const size_t M=A->size1;
    const size_t N=A->size2;
    size_t K;
    if (M<N) K=M;
    else K=N;

    if (M < N)
      {
	GSL_ERROR ("svd of MxN matrix, M<N, is not implemented", GSL_EUNIMPL);
      }
    else if (V->size1 != N)
      {
	GSL_ERROR ("square matrix V must match second dimension of matrix A",
		   GSL_EBADLEN);
      }
    else if (V->size1 != V->size2)
      {
	GSL_ERROR ("matrix V must be square", GSL_ENOTSQR);
      }
    else if (S->size != N)
      {
	GSL_ERROR ("length of vector S must match second dimension of matrix A",
		   GSL_EBADLEN);
      }
    else if (work->size != N)
      {
	GSL_ERROR ("length of workspace must match second dimension of matrix A",
		   GSL_EBADLEN);
      }

    /* Handle the case of N=1 (SVD of a column vector) */

    if (N == 1)
      {
	gsl_vector_view column=gsl_matrix_column (A, 0);
	double norm=gsl_blas_dnrm2 (&column.vector);

	gsl_vector_set (S, 0, norm); 
	gsl_matrix_set (V, 0, 0, 1.0);
      
	if (norm != 0.0)
	  {
	    gsl_blas_dscal (1.0/norm, &column.vector);
	  }

	return GSL_SUCCESS;
      }
  
    {
      gsl_vector_view f=gsl_vector_subvector (work, 0, K - 1);
    
      /* bidiagonalize matrix A, unpack A into U S V */
    
      gsl_linalg_bidiag_decomp (A, S, &f.vector);

      //std::cout << "A: " << gsl_matrix_get(A,0,0) << " "
      //<< gsl_matrix_get(A,M-1,N-1) << std::endl;
      //std::cout << "S: " << S->data[0] << " " 
      //<< S->data[S->size-1] 
      //<< std::endl;
    
      gsl_linalg_bidiag_unpack2 (A, S, &f.vector, V);

      //std::cout << "S2: " << S->data[0] << " " 
      //<< S->data[S->size-1] 
      //<< std::endl;
    
      /* apply reduction steps to B=(S,Sd) */
    
      chop_small_elements (S, &f.vector);
    
      //std::cout << "S3: " << S->data[0] << " " 
      //<< S->data[S->size-1] 
      //<< std::endl;
    
      /* Progressively reduce the matrix until it is diagonal */
    
      b=N - 1;
      iter=0;

      while (b > 0)
	{
	  double fbm1=gsl_vector_get (&f.vector, b - 1);

	  if (fbm1 == 0.0 || gsl_isnan (fbm1))
	    {
	      b--;
	      continue;
	    }

	  //std::cout << "b,fbm1: " << b << " " << fbm1 << std::endl;
        
	  /* Find the largest unreduced block (a,b) starting from b
	     and working backwards */

	  a=b - 1;

	  while (a > 0)
	    {
	      double fam1=gsl_vector_get (&f.vector, a - 1);

	      if (fam1 == 0.0 || gsl_isnan (fam1))
		{
		  break;
		}
            
	      a--;

	      //std::cout << "a,fam1: " << a << " " << fam1 << std::endl;
	    }

	  iter++;
        
	  if (iter > 100 * N) 
	    {
	      GSL_ERROR("SVD decomposition failed to converge", GSL_EMAXITER);
	    }

        
	  {
	    const size_t n_block=b - a + 1;
	    gsl_vector_view S_block=gsl_vector_subvector (S, a, n_block);
	    gsl_vector_view f_block=gsl_vector_subvector 
	      (&f.vector, a, n_block - 1);
          
	    gsl_matrix_view U_block =
	      gsl_matrix_submatrix (A, 0, a, A->size1, n_block);
	    gsl_matrix_view V_block =
	      gsl_matrix_submatrix (V, 0, a, V->size1, n_block);
          
	    int rescale=0;
	    double scale=1; 
	    double norm=0;

	    /* Find the maximum absolute values of the diagonal and 
               subdiagonal */

	    for (i=0; i < n_block; i++) 
	      {
		double s_i=gsl_vector_get (&S_block.vector, i);
		double ax=fabs(s_i);
		if (ax > norm) norm=ax;
		//std::cout << "aa: " << ax << std::endl;
	      }

	    for (i=0; i < n_block - 1; i++) 
	      {
		double f_i=gsl_vector_get (&f_block.vector, i);
		double ax=fabs(f_i);
		if (ax > norm) norm=ax;
		//std::cout << "aa2: " << ax << std::endl;
	      }

	    /* Temporarily scale the submatrix if necessary */

	    if (norm > GSL_SQRT_DBL_MAX)
	      {
		scale=(norm / GSL_SQRT_DBL_MAX);
		rescale=1;
	      }
	    else if (norm < GSL_SQRT_DBL_MIN && norm > 0)
	      {
		scale=(norm / GSL_SQRT_DBL_MIN);
		rescale=1;
	      }

	    //std::cout << "rescale: " << rescale << std::endl;

	    if (rescale) 
	      {
		gsl_blas_dscal(1.0 / scale, &S_block.vector);
		gsl_blas_dscal(1.0 / scale, &f_block.vector);
	      }

	    /* Perform the implicit QR step */

	    /*
	      for(size_t ii=0;ii<M;ii++) {
	      for(size_t jj=0;jj<N;jj++) {
	      std::cout << ii << "." << jj << "." 
	      << gsl_matrix_get(A,ii,jj) << std::endl;
	      }
	      }
	      for(size_t ii=0;ii<N;ii++) {
	      for(size_t jj=0;jj<N;jj++) {
	      std::cout << "V: " << ii << "." << jj << "." 
	      << gsl_matrix_get(V,ii,jj) << std::endl;
	      }
	      }
	    */

	    qrstep (&S_block.vector, &f_block.vector, &U_block.matrix, 
		    &V_block.matrix);

	    /*
	      for(size_t ii=0;ii<M;ii++) {
	      for(size_t jj=0;jj<N;jj++) {
	      std::cout << ii << " " << jj << " " 
	      << gsl_matrix_get(A,ii,jj) << std::endl;
	      }
	      }
	      for(size_t ii=0;ii<N;ii++) {
	      for(size_t jj=0;jj<N;jj++) {
	      std::cout << "V: " << ii << " " << jj << " " 
	      << gsl_matrix_get(V,ii,jj) << std::endl;
	      }
	      }
	    */

	    /* remove any small off-diagonal elements */
          
	    chop_small_elements (&S_block.vector, &f_block.vector);
          
	    /* Undo the scaling if needed */

	    if (rescale)
	      {
		gsl_blas_dscal(scale, &S_block.vector);
		gsl_blas_dscal(scale, &f_block.vector);
	      }
	  }
        
	}
    }

    /* Make singular values positive by reflections if necessary */
  
    for (j=0; j < K; j++)
      {
	double Sj=gsl_vector_get (S, j);
      
	if (Sj < 0.0)
	  {
	    for (i=0; i < N; i++)
	      {
		double Vij=gsl_matrix_get (V, i, j);
		gsl_matrix_set (V, i, j, -Vij);
	      }
          
	    gsl_vector_set (S, j, -Sj);
	  }
      }
  
    /* Sort singular values into decreasing order */
  
    for (i=0; i < K; i++)
      {
	double S_max=gsl_vector_get (S, i);
	size_t i_max=i;
      
	for (j=i + 1; j < K; j++)
	  {
	    double Sj=gsl_vector_get (S, j);
          
	    if (Sj > S_max)
	      {
		S_max=Sj;
		i_max=j;
	      }
	  }
      
	if (i_max != i)
	  {
	    /* swap eigenvalues */
	    gsl_vector_swap_elements (S, i, i_max);
          
	    /* swap eigenvectors */
	    gsl_matrix_swap_columns (A, i, i_max);
	    gsl_matrix_swap_columns (V, i, i_max);
	  }
      }
  
    return GSL_SUCCESS;
  }

};

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);
  cout.precision(6);

  test_mgr t;
  t.set_output_level(2);

  if (false) {
  
    static const size_t arr_size=12;

    gsl_matrix *gm1=gsl_matrix_alloc(arr_size,arr_size);
    gsl_matrix *gm2=gsl_matrix_alloc(arr_size,arr_size);
    gsl_matrix *gm3=gsl_matrix_alloc(arr_size,arr_size);
    gsl_matrix *gm4=gsl_matrix_alloc(arr_size,arr_size);
    gsl_matrix *mat_base=gsl_matrix_alloc(arr_size,arr_size);
    gsl_vector *gv1=gsl_vector_alloc(arr_size);
    gsl_vector *gv2=gsl_vector_alloc(arr_size);
    gsl_vector *gv3=gsl_vector_alloc(arr_size);
    gsl_vector *gv4=gsl_vector_alloc(arr_size);
    ubmatrix om1(arr_size,arr_size);
    ubmatrix om2(arr_size,arr_size);
    ubmatrix om3(arr_size,arr_size);
    ubvector ov1(arr_size);
    ubvector ov2(arr_size);
    ubvector ov3(arr_size);
    ubvector ov4(arr_size);

    // ----------------------------------------------------------
    // Tests of basic SV decomposition

    {

      // Setup arbitrary matrix
      for(size_t i=0;i<arr_size;i++) {
	for(size_t j=0;j<arr_size;j++) {
	  gsl_matrix_set(gm1,i,j,sqrt(((double)i))+sin(((double)j)));
	}
      }

      // Compute initial SVD
      gsl_linalg_SV_decomp(gm1,gm2,gv1,gv2);

      // Create one with better singular values. Note that because we
      // are ensuring that the all the singular values are distinct (the
      // non-degenerate case) then the left- and right-singular vectors
      // are unique up to a choice in sign.
      for(size_t j=0;j<arr_size;j++) {
	gsl_matrix_set(gm3,j,j,((double)(j+1)));
      }

      // Reconstruct a new matrix, stored in mat_base, for future
      // SV decompositions
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,gm1,gm3,0.0,gm4);
      gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,gm4,gm2,0.0,mat_base);

      // Copy things back over to gm1 and om1
      for(size_t i=0;i<arr_size;i++) {
	for(size_t j=0;j<arr_size;j++) {
	  gsl_matrix_set(gm1,i,j,gsl_matrix_get(mat_base,i,j));
	  om1(i,j)=gsl_matrix_get(mat_base,i,j);
	}
      }

      // Test decomposition and solve
      gsl_linalg_SV_decomp(gm1,gm2,gv1,gv2);
      SV_decomp(arr_size,arr_size,om1,om2,ov1,ov2);
      gsl_linalg_SV_solve(gm1,gm2,gv1,gv3,gv4);
      SV_solve(arr_size,arr_size,om1,om2,ov1,ov3,ov4);

      t.test_abs_mat(arr_size,arr_size,om1,gsl_matrix_wrap(gm1),5.0e-14,"m1");
      t.test_abs_mat(arr_size,arr_size,om2,gsl_matrix_wrap(gm2),5.0e-14,"m2");
      t.test_rel_vec(arr_size,ov1,gsl_vector_wrap(gv1),1.0e-12,"v1");
      t.test_rel_vec(arr_size-1,ov2,gsl_vector_wrap(gv2),1.0e-12,"v2");
      t.test_rel_vec(arr_size,ov4,gsl_vector_wrap(gv4),1.0e-12,"v4 (solve)");

#ifdef O2SCL_ARMA

      arma::mat am1(arr_size,arr_size);
      arma::mat am2(arr_size,arr_size);
      arma::mat am3(arr_size,arr_size);
      arma::colvec av1(arr_size); 
      for(size_t i=0;i<arr_size;i++) {
	for(size_t j=0;j<arr_size;j++) {
	  am1(i,j)=gsl_matrix_get(mat_base,i,j);
	}
      }
      svd(am2,av1,am3,am1);

      // The columns are defined only up to a sign. Find which
      // columns are given with the opposite sign, and flip
      // their signs back
      for(size_t j=0;j<arr_size;j++) {
	if (fabs((om1(0,j)-am2(0,j))/om1(0,j))>1.0) {
	  for(size_t i=0;i<arr_size;i++) {
	    am2(i,j)*=-1.0;
	    am3(i,j)*=-1.0;
	  }
	}
      }

      t.test_rel_mat(arr_size,arr_size,om1,am2,1.0e-11,"arma m1");
      t.test_rel_mat(arr_size,arr_size,om2,am3,1.0e-11,"arma m2");
      t.test_rel_vec(arr_size,ov1,av1,1.0e-14,"arma v2");

#endif

#ifdef O2SCL_SET_EIGEN

      Eigen::MatrixXd em1(arr_size,arr_size);
      Eigen::MatrixXd em2(arr_size,arr_size);
      Eigen::MatrixXd em3(arr_size,arr_size);
      Eigen::VectorXd ev1(arr_size); 
      for(size_t i=0;i<arr_size;i++) {
	for(size_t j=0;j<arr_size;j++) {
	  em1(i,j)=gsl_matrix_get(mat_base,i,j);
	}
      }
    
      Eigen::JacobiSVD<Eigen::MatrixXd> 
	svd(em1,Eigen::ComputeThinU | Eigen::ComputeThinV);
      ev1=svd.singularValues();
      em2=svd.matrixU();
      em3=svd.matrixV();

      // The columns are defined only up to a sign. Find which
      // columns are given with the opposite sign, and flip
      // their signs back
      for(size_t j=0;j<arr_size;j++) {
	if (fabs((om1(0,j)-em2(0,j))/om1(0,j))>1.0) {
	  for(size_t i=0;i<arr_size;i++) {
	    em2(i,j)*=-1.0;
	    em3(i,j)*=-1.0;
	  }
	}
      }

      t.test_rel_mat(arr_size,arr_size,om1,em2,1.0e-11,"eigen m1");
      t.test_rel_mat(arr_size,arr_size,om2,em3,1.0e-11,"eigen m2");
      t.test_rel_vec(arr_size,ov1,ev1,1.0e-14,"eigen v2");

#endif

    }

    // ----------------------------------------------------------
    // Testing SV_decomp_mod()

    {

      // Setup original matrix
      for(size_t i=0;i<arr_size;i++) {
	for(size_t j=0;j<arr_size;j++) {
	  gsl_matrix_set(gm1,i,j,gsl_matrix_get(mat_base,i,j));
	  om1(i,j)=gsl_matrix_get(mat_base,i,j);
	}
      }
    
      // Test decomposition
      gsl_linalg_SV_decomp_mod(gm1,gm2,gm3,gv1,gv2);
      SV_decomp_mod(arr_size,arr_size,om1,om2,om3,ov1,ov2);

      t.test_abs_mat(arr_size,arr_size,om1,gsl_matrix_wrap(gm1),5.0e-14,
		     "mod m1");
      t.test_abs_mat(arr_size,arr_size,om2,gsl_matrix_wrap(gm2),5.0e-14,
		     "mod m2");
      t.test_abs_mat(arr_size,arr_size,om3,gsl_matrix_wrap(gm3),5.0e-14,
		     "mod m3");
      t.test_rel_vec(arr_size,ov1,gsl_vector_wrap(gv1),1.0e-12,
		     "mod v1");
      t.test_rel_vec(arr_size-1,ov2,gsl_vector_wrap(gv2),1.0e-12,
		     "mod v2");
    }
  
    // ----------------------------------------------------------
    // Testing SV_decomp_jacobi()

    {
      // Setup original matrix
      for(size_t i=0;i<arr_size;i++) {
	for(size_t j=0;j<arr_size;j++) {
	  gsl_matrix_set(gm1,i,j,gsl_matrix_get(mat_base,i,j));
	  om1(i,j)=gsl_matrix_get(mat_base,i,j);
	}
      }
    
      // Test decomposition
      gsl_linalg_SV_decomp_jacobi(gm1,gm2,gv1);
      SV_decomp_jacobi(arr_size,arr_size,om1,om2,ov1);
    
      t.test_abs_mat(arr_size,arr_size,om1,gsl_matrix_wrap(gm1),
		     5.0e-14,"jacobi m1");
      t.test_abs_mat(arr_size,arr_size,om2,gsl_matrix_wrap(gm2),
		     5.0e-14,"jacobi m2");
      t.test_rel_vec(arr_size,ov1,gsl_vector_wrap(gv1),
		     1.0e-12,"jacobi v1");

    }

  }

  t.report();
  return 0;
}

