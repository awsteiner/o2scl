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
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_linalg.h>

#include <o2scl/test_mgr.h>
#include <o2scl/svdstep.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

// The GSL svd2() function is static, so we use a simple copy
// for testing here
namespace gsl_svdstep {

  void create_givens (const double a, const double b, double *c, double *s) {
    if (b == 0)
      {
	*c = 1;
	*s = 0;
      }
    else if (fabs (b) > fabs (a))
      {
	double t = -a / b;
	double s1 = 1.0 / sqrt (1 + t * t);
	*s = s1;
	*c = s1 * t;
      }
    else
      {
	double t = -b / a;
	double c1 = 1.0 / sqrt (1 + t * t);
	*c = c1;
	*s = c1 * t;
      }
  }
  
  void chop_small_elements (gsl_vector * d, gsl_vector * f) {
    const size_t N = d->size;
    double d_i = gsl_vector_get (d, 0);

    size_t i;

    for (i = 0; i < N - 1; i++)
      {
	double f_i = gsl_vector_get (f, i);
	double d_ip1 = gsl_vector_get (d, i + 1);

	if (fabs (f_i) < GSL_DBL_EPSILON * (fabs (d_i) + fabs (d_ip1)))
	  {
	    gsl_vector_set (f, i, 0.0);
	  }

	d_i = d_ip1;
      }

  }

  double trailing_eigenvalue (const gsl_vector * d, const gsl_vector * f) {
    const size_t n = d->size;

    double da = gsl_vector_get (d, n - 2);
    double db = gsl_vector_get (d, n - 1);
    double fa = (n > 2) ? gsl_vector_get (f, n - 3) : 0.0;
    double fb = gsl_vector_get (f, n - 2);

    double mu;

    /* We can compute mu more accurately than using the formula above
       since we know the roots cannot be negative.  This also avoids
       the possibility of NaNs in the formula above.

       The matrix is [ da^2 + fa^2,  da fb      ;
       da fb      , db^2 + fb^2 ]
       and mu is the eigenvalue closest to the bottom right element.
    */
    
    double ta = da * da + fa * fa;
    double tb = db * db + fb * fb;
    double tab = da * fb;
    
    double dt = (ta - tb) / 2.0;
    
    double S = ta + tb;
    double da2 = da * da, db2 = db * db;
    double fa2 = fa * fa, fb2 = fb * fb;
    double P = (da2 * db2) + (fa2 * db2) + (fa2 * fb2);
    double D = hypot(dt, tab);
    double r1 = S/2 + D;
    
    if (dt >= 0)
      {
        /* tb < ta, choose smaller root */
        mu = (r1 > 0) ?  P / r1 : 0.0;
      }
    else 
      {
        /* tb > ta, choose larger root */
        mu = r1;
      }

    return mu;
  }

  void create_schur (double d0, double f0, double d1, double * c, 
		     double * s) {
    double apq = 2.0 * d0 * f0;

    if (d0 == 0 || f0 == 0)
      {
	*c = 1.0;
	*s = 0.0;
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
	scale = ldexp(1.0, -(d0_exp + f0_exp)/4);
	d0 *= scale;
	f0 *= scale;
	d1 *= scale;
	apq = 2.0 * d0 * f0;
      }

    if (apq != 0.0)
      {
	double t;
	double tau = (f0*f0 + (d1 + d0)*(d1 - d0)) / apq;
      
	if (tau >= 0.0)
	  {
	    t = 1.0/(tau + hypot(1.0, tau));
	  }
	else
	  {
	    t = -1.0/(-tau + hypot(1.0, tau));
	  }

	*c = 1.0 / hypot(1.0, t);
	*s = t * (*c);
      }
    else
      {
	*c = 1.0;
	*s = 0.0;
      }
  }

  void svd2 (gsl_vector * d, gsl_vector * f, gsl_matrix * U, 
	     gsl_matrix * V) {
    size_t i;
    double c, s, a11, a12, a21, a22;

    const size_t M = U->size1;
    const size_t N = V->size1;

    double d0 = gsl_vector_get (d, 0);
    double f0 = gsl_vector_get (f, 0);
  
    double d1 = gsl_vector_get (d, 1);

    if (d0 == 0.0)
      {
	/* Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] */

	create_givens (f0, d1, &c, &s);

	/* compute B <= G^T B X,  where X = [0,1;1,0] */

	gsl_vector_set (d, 0, c * f0 - s * d1);
	gsl_vector_set (f, 0, s * f0 + c * d1);
	gsl_vector_set (d, 1, 0.0);

	/* Compute U <= U G */

	for (i = 0; i < M; i++)
	  {
	    double Uip = gsl_matrix_get (U, i, 0);
	    double Uiq = gsl_matrix_get (U, i, 1);
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

	for (i = 0; i < N; i++)
	  {
	    double Vip = gsl_matrix_get (V, i, 0);
	    double Viq = gsl_matrix_get (V, i, 1);
	    gsl_matrix_set (V, i, 0, c * Vip - s * Viq);
	    gsl_matrix_set (V, i, 1, s * Vip + c * Viq);
	  }

	return;
      }
    else
      {
	/* Make columns orthogonal, A = [d0, f0; 0, d1] * G */

	create_schur (d0, f0, d1, &c, &s);

	/* compute B <= B G */
      
	a11 = c * d0 - s * f0;
	a21 = - s * d1;
      
	a12 = s * d0 + c * f0;
	a22 = c * d1;
      
	/* Compute V <= V G */
      
	for (i = 0; i < N; i++)
	  {
	    double Vip = gsl_matrix_get (V, i, 0);
	    double Viq = gsl_matrix_get (V, i, 1);
	    gsl_matrix_set (V, i, 0, c * Vip - s * Viq);
	    gsl_matrix_set (V, i, 1, s * Vip + c * Viq);
	  }
      
	/* Eliminate off-diagonal elements, bring column with largest
	   norm to first column */
      
	if (hypot(a11, a21) < hypot(a12,a22))
	  {
	    double t1, t2;

	    /* B <= B X */

	    t1 = a11; a11 = a12; a12 = t1;
	    t2 = a21; a21 = a22; a22 = t2;

	    /* V <= V X */

	    gsl_matrix_swap_columns(V, 0, 1);
	  } 

	create_givens (a11, a21, &c, &s);
      
	/* compute B <= G^T B */
      
	gsl_vector_set (d, 0, c * a11 - s * a21);
	gsl_vector_set (f, 0, c * a12 - s * a22);
	gsl_vector_set (d, 1, s * a12 + c * a22);
      
	/* Compute U <= U G */
      
	for (i = 0; i < M; i++)
	  {
	    double Uip = gsl_matrix_get (U, i, 0);
	    double Uiq = gsl_matrix_get (U, i, 1);
	    gsl_matrix_set (U, i, 0, c * Uip - s * Uiq);
	    gsl_matrix_set (U, i, 1, s * Uip + c * Uiq);
	  }

	return;
      }
  }


  void chase_out_intermediate_zero (gsl_vector * d, gsl_vector * f, 
				    gsl_matrix * U, size_t k0) {

    const size_t M = U->size1;
    const size_t n = d->size;
    double c, s;
    double x, y;
    size_t k;

    x = gsl_vector_get (f, k0);
    y = gsl_vector_get (d, k0+1);

    for (k = k0; k < n - 1; k++)
      {
	create_givens (y, -x, &c, &s);
      
	/* Compute U <= U G */

	{
	  size_t i;

	  for (i = 0; i < M; i++)
	    {
	      double Uip = gsl_matrix_get (U, i, k0);
	      double Uiq = gsl_matrix_get (U, i, k + 1);
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
	    double z = gsl_vector_get (f, k + 1);
	    gsl_vector_set (f, k + 1, c * z); 

	    x = -s * z ;
	    y = gsl_vector_get (d, k + 2); 
	  }
      }
  }

  void chase_out_trailing_zero (gsl_vector * d, gsl_vector * f, 
				gsl_matrix * V) {

    const size_t N = V->size1;
    const size_t n = d->size;
    double c, s;
    double x, y;
    size_t k;

    x = gsl_vector_get (d, n - 2);
    y = gsl_vector_get (f, n - 2);

    for (k = n - 1; k-- > 0;)
      {
	create_givens (x, y, &c, &s);

	/* Compute V <= V G where G = [c, s ; -s, c] */

	{
	  size_t i;
   
	  for (i = 0; i < N; i++)
	    {
	      double Vip = gsl_matrix_get (V, i, k);
	      double Viq = gsl_matrix_get (V, i, n - 1);
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
	    double z = gsl_vector_get (f, k - 1);
	    gsl_vector_set (f, k - 1, c * z); 

	    x = gsl_vector_get (d, k - 1); 
	    y = s * z ;
	  }
      }
  }

  void qrstep (gsl_vector * d, gsl_vector * f, gsl_matrix * U, 
	       gsl_matrix * V) {

    const size_t M = U->size1;
    const size_t N = V->size1;
    const size_t n = d->size;
    double y, z;
    double ak, bk, zk, ap, bp, aq, bq;
    size_t i, k;

    if (n == 1)
      return;  /* shouldn't happen */

    /* Compute 2x2 svd directly */

    if (n == 2)
      {
	svd2 (d, f, U, V);
	return;
      }

    /* Chase out any zeroes on the diagonal */

    for (i = 0; i < n - 1; i++)
      {
	double d_i = gsl_vector_get (d, i);
      
	if (d_i == 0.0)
	  {
	    chase_out_intermediate_zero (d, f, U, i);
	    return;
	  }
      }

    /* Chase out any zero at the end of the diagonal */

    {
      double d_nm1 = gsl_vector_get (d, n - 1);

      if (d_nm1 == 0.0) 
	{
	  chase_out_trailing_zero (d, f, V);
	  return;
	}
    }


    /* Apply QR reduction steps to the diagonal and offdiagonal */

    {
      double d0 = gsl_vector_get (d, 0);
      double f0 = gsl_vector_get (f, 0);
    
      double d1 = gsl_vector_get (d, 1);
      double f1 = gsl_vector_get (f, 1);
    
      {
	double mu = trailing_eigenvalue (d, f);
    
	y = d0 * d0 - mu;
	z = d0 * f0;
      }
    
      /* Set up the recurrence for Givens rotations on a bidiagonal matrix */
    
      ak = 0;
      bk = 0;
    
      ap = d0;
      bp = f0;
    
      aq = d1;
      bq = f1;
    }

    for (k = 0; k < n - 1; k++)
      {
	double c, s;
	create_givens (y, z, &c, &s);

	/* Compute V <= V G */

	for (i = 0; i < N; i++)
	  {
	    double Vip = gsl_matrix_get (V, i, k);
	    double Viq = gsl_matrix_get (V, i, k + 1);
	    gsl_matrix_set (V, i, k, c * Vip - s * Viq);
	    gsl_matrix_set (V, i, k + 1, s * Vip + c * Viq);
	  }

	/* compute B <= B G */

	{
	  double bk1 = c * bk - s * z;

	  double ap1 = c * ap - s * bp;
	  double bp1 = s * ap + c * bp;
	  double zp1 = -s * aq;

	  double aq1 = c * aq;

	  if (k > 0)
	    {
	      gsl_vector_set (f, k - 1, bk1);
	    }

	  ak = ap1;
	  bk = bp1;
	  zk = zp1;

	  ap = aq1;

	  if (k < n - 2)
	    {
	      bp = gsl_vector_get (f, k + 1);
	    }
	  else
	    {
	      bp = 0.0;
	    }

	  y = ak;
	  z = zk;
	}

	create_givens (y, z, &c, &s);

	/* Compute U <= U G */

	for (i = 0; i < M; i++)
	  {
	    double Uip = gsl_matrix_get (U, i, k);
	    double Uiq = gsl_matrix_get (U, i, k + 1);
	    gsl_matrix_set (U, i, k, c * Uip - s * Uiq);
	    gsl_matrix_set (U, i, k + 1, s * Uip + c * Uiq);
	  }

	/* compute B <= G^T B */

	{
	  double ak1 = c * ak - s * zk;
	  double bk1 = c * bk - s * ap;
	  double zk1 = -s * bp;

	  double ap1 = s * bk + c * ap;
	  double bp1 = c * bp;

	  gsl_vector_set (d, k, ak1);

	  ak = ak1;
	  bk = bk1;
	  zk = zk1;

	  ap = ap1;
	  bp = bp1;

	  if (k < n - 2)
	    {
	      aq = gsl_vector_get (d, k + 2);
	    }
	  else
	    {
	      aq = 0.0;
	    }

	  y = bk;
	  z = zk;
	}
      }

    gsl_vector_set (f, n - 2, bk);
    gsl_vector_set (d, n - 1, ap);
  }

};

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  void *vpx=0;
  size_t tmp;

  cout.setf(ios::scientific);
  static const size_t arr_size=5;

  using namespace o2scl_cblas;
  using namespace o2scl_linalg;

  // GSL

  gsl_vector *gv1=gsl_vector_alloc(arr_size);
  gsl_vector *gv2=gsl_vector_alloc(arr_size);
  gsl_matrix *gm1=gsl_matrix_alloc(arr_size,arr_size);
  gsl_matrix *gm2=gsl_matrix_alloc(arr_size,arr_size);

  for(size_t i=0;i<arr_size;i++) {
    gsl_vector_set(gv1,i,i/100.0);
    gsl_vector_set(gv2,i,(i+1)/100.0);
    for(size_t j=0;j<arr_size;j++) {
      gsl_matrix_set(gm1,i,j,exp(-((double)(i-j)*(i-j))/2.0));
      gsl_matrix_set(gm2,i,j,exp(-((double)(i-j)*(i-j))/2.0));
    }
  }

  gsl_svdstep::svd2(gv1,gv2,gm1,gm2);

  // O2scl
  ubmatrix om1(arr_size,arr_size), om2(arr_size,arr_size);
  ubvector ov1(arr_size), ov2(arr_size);

  for(size_t i=0;i<arr_size;i++) {
    ov1[i]=i/100.0;
    ov2[i]=(i+1)/100.0;
    for(size_t j=0;j<arr_size;j++) {
      om1(i,j)=exp(-((double)(i-j)*(i-j))/2.0);
      om2(i,j)=exp(-((double)(i-j)*(i-j))/2.0);
    }
  }

  o2scl_linalg::svd2(arr_size,arr_size,ov1,ov2,om1,om2);
  
  t.test_rel_vec(arr_size,ov1,gsl_vector_wrap(gv1),1.0e-10,"svd2");
  t.test_rel_vec(arr_size,ov2,gsl_vector_wrap(gv2),1.0e-10,"svd2");
  t.test_rel_mat(arr_size,arr_size,om1,gsl_matrix_wrap(gm1),1.0e-10,"svd2");
  t.test_rel_mat(arr_size,arr_size,om2,gsl_matrix_wrap(gm2),1.0e-10,"svd2");

  for(size_t i=0;i<arr_size;i++) {
    gsl_vector_set(gv1,i,i/100.0);
    gsl_vector_set(gv2,i,(i+1)/100.0);
    for(size_t j=0;j<arr_size;j++) {
      gsl_matrix_set(gm1,i,j,exp(-((double)(i-j)*(i-j))/2.0));
      gsl_matrix_set(gm2,i,j,exp(-((double)(i-j)*(i-j))/2.0));
    }
  }

  gsl_svdstep::qrstep(gv1,gv2,gm1,gm2);

  for(size_t i=0;i<arr_size;i++) {
    ov1[i]=i/100.0;
    ov2[i]=(i+1)/100.0;
    for(size_t j=0;j<arr_size;j++) {
      om1(i,j)=exp(-((double)(i-j)*(i-j))/2.0);
      om2(i,j)=exp(-((double)(i-j)*(i-j))/2.0);
    }
  }

  o2scl_linalg::qrstep(arr_size,arr_size,arr_size,ov1,ov2,om1,om2);
  
  t.test_rel_vec(arr_size,ov1,gsl_vector_wrap(gv1),1.0e-10,"qrstep");
  t.test_rel_vec(arr_size,ov2,gsl_vector_wrap(gv2),1.0e-10,"qrstep");
  t.test_rel_mat(arr_size,arr_size,om1,gsl_matrix_wrap(gm1),1.0e-10,"qrstep");
  t.test_rel_mat(arr_size,arr_size,om2,gsl_matrix_wrap(gm2),1.0e-10,"qrstep");

  t.report();
  return 0;
}

