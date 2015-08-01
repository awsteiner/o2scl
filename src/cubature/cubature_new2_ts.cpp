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
#include <o2scl/cubature_new.h>
#include <o2scl/test_mgr.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::vector_range<ubvector> ubvector_range;

using namespace std;
using namespace o2scl;

static const bool debug=true;

int cub_count = 0;
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

/* Simple product function */
double f02(unsigned dim, const ubvector &x, void *params) {
  double prod = 1.0;
  unsigned int i;
  for (i = 0; i < dim; ++i)
    prod *= 2.0 * x[i];
  return prod;
}

#define K_2_SQRTPI 1.12837916709551257390

/* Gaussian centered at 1/2. */
double f12(unsigned dim, const ubvector &x, void *params) {
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
double f22(unsigned dim, const ubvector &x, void *params) {
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
double f32(unsigned dim, const ubvector &x, void *params) {
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
static double morokoff2(unsigned dim, const ubvector &x, void *params) {
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
  ++cub_count;
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

int f_test2(unsigned dim, const ubvector_range &x, 
	    unsigned fdim, ubvector_range &retval) {
  
  cout << retval[0] << endl;
  cout << x[0] << endl;
  
  double val;
  unsigned i, j;
  ++cub_count;
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
      val = f02(dim, x, &fdata);
      break;
    case 4:
      val = f12(dim, x, &fdata);
      break;
    case 5:
      val = f22(dim, x, &fdata);
      break;
    case 6:
      val = f32(dim, x, &fdata);
      break;
    case 7:
      val = morokoff2(dim, x, &fdata);
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

int fv(unsigned ndim, size_t npt, const double *x, unsigned fdim,
       double *fval) {
  for (unsigned i = 0; i < npt; ++i) {
    if (f_test(ndim, x + i*ndim, 0, fdim, fval + i*fdim)) {
      return o2scl::gsl_failure;
    }
  }
  return o2scl::success;
}

int fv2(unsigned ndim, size_t npt, const ubvector &x, unsigned fdim,
	ubvector &fval) {
  for (unsigned i = 0; i < npt; ++i) {
    ubvector_range x_shift=vector_range(x,i*ndim,(i+1)*ndim);
    ubvector_range fval_shift=vector_range(fval,i*fdim,(i+1)*fdim);
    f_test2(ndim, x_shift, fdim, fval_shift);
  }
  if (debug) {
    cout << ndim << " " << fdim << " " << npt << endl;
    for(unsigned i=0;i<npt;i++) {
      ubvector x_shift=vector_range(x,i*ndim,(i+1)*ndim);
      ubvector fval_shift=vector_range(fval,i*fdim,(i+1)*fdim);
      for(size_t k=0;k<ndim;k++) {
	cout << x_shift[k] << " ";
      }
      for(size_t k=0;k<fdim;k++) {
	cout << fval_shift[k] << " ";
      }
      cout << endl;
    }
    char ch;
    cin >> ch;
  }
  return o2scl::success;
}

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr tmgr;
  tmgr.set_output_level(1);

  size_t dim=3;
  double xmin[3], xmax[3];
  ubvector xmin2(3), xmax2(3);
  for (size_t i=0;i<dim;++i) {
    xmin[i]=0.0;
    xmax[i]=1.0;
    xmin2[i]=0.0;
    xmax2[i]=1.0;
  }

  typedef std::function<
    int(unsigned,size_t,const double *,unsigned,double *)> cub_funct_arr;
  typedef std::function<
    int(unsigned,size_t,const ubvector &,unsigned,ubvector &)>
    cub_funct_arr2;
  inte_hcubature_new<cub_funct_arr> hc;
  inte_pcubature_new<cub_funct_arr2,ubvector > pc;

  cub_funct_arr cfa=fv;
  cub_funct_arr2 cfa2=fv2;

  /*std::function<int(unsigned,const double *,unsigned,double *)> cfa=
    std::bind(f_test,std::placeholders::_1,std::placeholders::_2,0,
    std::placeholders::_3,std::placeholders::_4);
  */
    
  int test_n[14]={33,125,693,4913,70785,33,3861,35937,3465,35937,297,
		  729,33,729};
  
  double test_vals[14][3]={{5.958229e-01,3.519922e-06,3.523658e-07},
			   {5.958236e-01,2.130785e-04,3.832854e-07},
			   {1.002290e+00,9.980917e-03,2.290472e-03},
			    {9.999119e-01,1.113448e-03,8.812269e-05},
			    {6.514615e-02,6.405123e-04,7.924271e-04},
			    {1.000000e+00,2.220446e-16,2.220446e-16},
			    {1.000753e+00,9.612568e-03,7.526466e-04},
			    {1.000000e+00,2.155111e-04,1.324296e-08},
			    {9.852783e-01,9.774575e-03,1.472168e-02},
			   {9.999963e-01,7.175992e-05,3.650226e-06},
			   {9.998328e-01,7.738486e-03,1.671812e-04},
			   {9.999948e-01,1.425689e-03,5.187945e-06},
			   {1.001055e+00,4.808302e-03,1.055387e-03},
			   {9.967782e-01,6.471054e-03,3.221771e-03}};
  
  int tcnt=0;
  for(size_t test_iand=0;test_iand<8;test_iand++) {

    double tol, val, err;
    unsigned maxEval;

    tol=1.0e-2;
    maxEval=0;
    
    inte_hcubature_new<cub_funct_arr>::error_norm enh=
      inte_hcubature_new<cub_funct_arr>::ERROR_INDIVIDUAL;
    inte_pcubature_new<cub_funct_arr2,ubvector >::error_norm enp=
      inte_pcubature_new<cub_funct_arr2,ubvector >::ERROR_INDIVIDUAL;
    
    which_integrand = test_iand; 
    
    if (test_iand!=2) {

      cub_count=0;
      hc.integ(1,cfa,dim,xmin,xmax,maxEval,0,tol,enh,&val,&err);
	       
      cout << "# " << which_integrand << " " 
	   << "integral " << val << " " << "est. error " << err << " " 
	   << "true error " 
	   << fabs(val-exact_integral(which_integrand,dim,xmax)) << endl;
      cout << "evals " << cub_count << endl;

      tmgr.test_gen(fabs(val-exact_integral(which_integrand,dim,xmax))<
		    err*2.0,"hcub 2");
      tmgr.test_gen(test_n[tcnt]==cub_count,"cub_count");
      tmgr.test_rel(val,test_vals[tcnt][0],5.0e-6,"val");
      tmgr.test_rel(err,test_vals[tcnt][1],5.0e-6,"err");
      tmgr.test_rel(fabs(val-exact_integral(which_integrand,dim,xmax)),
		    test_vals[tcnt][2],5.0e-6,"diff w/ exact");
      tcnt++;

      cout << endl;
    }

    if (test_iand!=3) {

      ubvector vval(1), verr(1);

      cub_count=0;
      pc.integ(1,cfa2,dim,xmin2,xmax2,maxEval,0,tol,enp,vval,verr);
	       
      cout << "# " << which_integrand << " " 
	   << "integral " << vval[0] << " " << "est. error " << verr[0] << " " 
	   << "true error " 
	   << fabs(vval[0]-exact_integral(which_integrand,dim,xmax)) << endl;
      cout << "evals " << cub_count << endl;

      tmgr.test_gen(fabs(vval[0]-exact_integral(which_integrand,dim,xmax))<
		    verr[0]*2.0,"pcub 2");
      tmgr.test_gen(test_n[tcnt]==cub_count,"cub_count");
      tmgr.test_rel(vval[0],test_vals[tcnt][0],5.0e-6,"val");
      tmgr.test_rel(verr[0],test_vals[tcnt][1],5.0e-6,"err");
      tmgr.test_rel(fabs(vval[0]-exact_integral(which_integrand,dim,xmax)),
		    test_vals[tcnt][2],5.0e-6,"diff w/ exact");
      tcnt++;
      
      cout << endl;
    }
    
  }

  tmgr.report();
  return 0;
}
