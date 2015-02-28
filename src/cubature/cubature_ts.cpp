/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015, Andrew W. Steiner
  
  This file is part of O2scl. It has been adapted from cubature 
  written by Steven G. Johnson. 
  
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
 * Copyright (c) 2005-2013 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <iostream>
#include <o2scl/cubature.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int count = 0;
int which_integrand;
const double radius = 0.50124145262344534123412; /* random */

/* Simple constant function */
double fconst (double x[], size_t dim, void *params) {
  return 1;
}

/*** f0, f1, f2, and f3 are test functions from the Monte-Carlo
     integration routines in GSL 1.6 (monte/test.c).  Copyright (c)
     1996-2000 Michael Booth, GNU GPL. ****/

/* Simple product function */
double f0 (unsigned dim, const double *x, void *params) {
  double prod = 1.0;
  unsigned int i;
  for (i = 0; i < dim; ++i)
    prod *= 2.0 * x[i];
  return prod;
}

#define K_2_SQRTPI 1.12837916709551257390

/* Gaussian centered at 1/2. */
double f1 (unsigned dim, const double *x, void *params) {
  double a = *(double *)params;
  double sum = 0.;
  unsigned int i;
  for (i = 0; i < dim; i++) {
    double dx = x[i] - 0.5;
    sum += dx * dx;
  }
  return (pow (K_2_SQRTPI / (2. * a), (double) dim) *
	  exp (-sum / (a * a)));
}

/* double gaussian */
double f2 (unsigned dim, const double *x, void *params) {
  double a = *(double *)params;
  double sum1 = 0.;
  double sum2 = 0.;
  unsigned int i;
  for (i = 0; i < dim; i++) {
    double dx1 = x[i] - 1. / 3.;
    double dx2 = x[i] - 2. / 3.;
    sum1 += dx1 * dx1;
    sum2 += dx2 * dx2;
  }
  return 0.5 * pow (K_2_SQRTPI / (2. * a), dim) 
    * (exp (-sum1 / (a * a)) + exp (-sum2 / (a * a)));
}

/* Tsuda's example */
double f3 (unsigned dim, const double *x, void *params) {
  double c = *(double *)params;
  double prod = 1.;
  unsigned int i;
  for (i = 0; i < dim; i++)
    prod *= c / (c + 1) * pow((c + 1) / (c + x[i]), 2.0);
  return prod;
}

/* test integrand from W. J. Morokoff and R. E. Caflisch, "Quasi=
   Monte Carlo integration," J. Comput. Phys 122, 218-230 (1995).
   Designed for integration on [0,1]^dim, integral = 1. */
static double morokoff(unsigned dim, const double *x, void *params) {
  double p = 1.0 / dim;
  double prod = pow(1 + p, dim);
  unsigned int i;
  for (i = 0; i < dim; i++)
    prod *= pow(x[i], p);
  return prod;
}

/*** end of GSL test functions ***/

int f_test(unsigned dim, const double *x, void *data_,
	   unsigned fdim, double *retval) {
  
  double val;
  unsigned i, j;
  ++count;
  (void) data_; /* not used */
  for (j = 0; j < 1; ++j) {
    double fdata = which_integrand == 6 ? (1.0+sqrt (10.0))/9.0 : 0.1;
    switch (which_integrand) {
    case 0: /* simple smooth (separable) objective: prod. cos(x[i]). */
      val = 1;
      for (i = 0; i < dim; ++i)
	val *= cos(x[i]);
      break;
    case 1: { /* integral of exp(-x^2), rescaled to (0,infinity) limits */
      double scale = 1.0;
      val = 0;
      for (i = 0; i < dim; ++i) {
	if (x[i] > 0) {
	  double z = (1 - x[i]) / x[i];
	  val += z * z;
	  scale *= K_2_SQRTPI / (x[i] * x[i]);
	}
	else {
	  scale = 0;
	  break;
	}
      }
      val = exp(-val) * scale;
      break;
    }
    case 2: /* discontinuous objective: volume of hypersphere */
      val = 0;
      for (i = 0; i < dim; ++i)
	val += x[i] * x[i];
      val = val < radius * radius;
      break;
    case 3:
      val = f0(dim, x, &fdata);
      break;
    case 4:
      val = f1(dim, x, &fdata);
      break;
    case 5:
      val = f2(dim, x, &fdata);
      break;
    case 6:
      val = f3(dim, x, &fdata);
      break;
    case 7:
      val = morokoff(dim, x, &fdata);
      break;
    default:
      cout << "Unknown integrand." << endl;
      exit(-1);
    }
    retval[j] = val;
  }
  return 0;
}

#define K_PI 3.14159265358979323846

/* surface area of n-dimensional unit hypersphere */
static double S(unsigned n) {
  double val;
  int fact = 1;
  if (n % 2 == 0) { /* n even */
    val = 2 * pow(K_PI, n * 0.5);
    n = n / 2;
    while (n > 1) fact *= (n -= 1);
    val /= fact;
  }
  else { /* n odd */
    val = (1 << (n/2 + 1)) * pow(K_PI, n/2);
    while (n > 2) fact *= (n -= 2);
    val /= fact;
  }
  return val;
}

static double exact_integral(int which, unsigned dim, const double *xmax) {
  unsigned i;
  double val;
  switch(which) {
  case 0:
    val = 1;
    for (i = 0; i < dim; ++i)
      val *= sin(xmax[i]);
    break;
  case 2:
    val = dim == 0 ? 1 : S(dim) * pow(radius * 0.5, dim) / dim;
    break;
  default:
    val = 1.0;
  }
  return val;
}

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr tmgr;
  tmgr.set_output_level(2);

  size_t dim=3;
  double xmin[3], xmax[3];
  for (size_t i=0;i<dim;++i) {
    xmin[i]=0.0;
    xmax[i]=1.0;
  }
  
  //typedef int (*cfunc)(unsigned, const double *, void *, unsigned, double *);
  //typedef int (*cfunc_v)(unsigned, size_t, const double *, void *,
  //unsigned, double *);
  //hcub_class<cfunc,cfunc_v> hc;
  
  for(size_t test_iand=0;test_iand<8;test_iand++) {

    double tol, val, err;
    unsigned maxEval;

    tol=1.0e-2;
    maxEval=0;

    //hcub_class<cfunc,cfunc_v>::error_norm en=
    //hcub_class<cfunc,cfunc_v>::ERROR_INDIVIDUAL;
    
    which_integrand = test_iand; 
    
    if (test_iand!=2) {

      hcubature(1,f_test,0,dim,xmin,xmax, 
		maxEval,0,tol,ERROR_INDIVIDUAL,&val,&err);
      //hc.integ(1,&f_test,0,dim,xmin,xmax, 
      //maxEval,0,tol,en,&val,&err);
      
      cout << "# " << which_integrand << " " 
	   << "integral " << val << " " << "est. error " << err << " " 
	   << "true error " 
	   << fabs(val-exact_integral(which_integrand,dim,xmax)) << endl;
      cout << "evals " << count << endl;
      tmgr.test_gen(fabs(val-exact_integral(which_integrand,dim,xmax))<
		    err*2.0,"hcub");
    }

    if (test_iand!=3) {

      pcubature(1,f_test,0,dim,xmin,xmax, 
		maxEval,0,tol,ERROR_INDIVIDUAL,&val,&err);
      
      cout << "# " << which_integrand << " " 
	   << "integral " << val << " " << "est. error " << err << " " 
	   << "true error " 
	   << fabs(val-exact_integral(which_integrand,dim,xmax)) << endl;
      cout << "evals " << count << endl;
      tmgr.test_gen(fabs(val-exact_integral(which_integrand,dim,xmax))<
		    err*2.0,"pcub");
    }
    
  }

  tmgr.report();
  return 0;
}
