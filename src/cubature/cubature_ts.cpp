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
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/vector.h>
#include <o2scl/cubature.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::vector_range<ubvector> ubvector_range;
typedef boost::numeric::ublas::vector_range<ubvector_range>
ubvector_range_range;

typedef boost::numeric::ublas::vector_range<const ubvector> ubvector_crange;
typedef boost::numeric::ublas::vector_range<const ubvector_range>
ubvector_crange_range;
typedef boost::numeric::ublas::vector_range<const ubvector_crange>
ubvector_crange_crange;

using namespace std;
using namespace o2scl;

int cub_count = 0;
int which_integrand;
static const double radius = 0.50124145262344534123412;
static const double k_2_sqrtpi=2.0/sqrt(o2scl_const::pi);

/*** f0, f1, f2, and f3 are test functions from the Monte-Carlo
     integration routines in GSL 1.6 (monte/test.c).  Copyright (c)
     1996-2000 Michael Booth, GNU GPL. ****/

/* Simple product function */
double f0 (size_t dim, const double *x) {
  double prod = 1.0;
  size_t i;
  for (i = 0; i < dim; ++i) {
    prod *= 2.0 * x[i];
  }
  return prod;
}

/* Gaussian centered at 1/2. */
double f1 (size_t dim, const double *x, double a) {
  double sum = 0.;
  size_t i;
  for (i = 0; i < dim; i++) {
    double dx = x[i] - 0.5;
    sum += dx * dx;
  }
  return (pow (k_2_sqrtpi / (2. * a), (double) dim) *
	  exp (-sum / (a * a)));
}

/* double gaussian */
double f2 (size_t dim, const double *x, double a) {
  double sum1 = 0.;
  double sum2 = 0.;
  size_t i;
  for (i = 0; i < dim; i++) {
    double dx1 = x[i] - 1. / 3.;
    double dx2 = x[i] - 2. / 3.;
    sum1 += dx1 * dx1;
    sum2 += dx2 * dx2;
  }
  return 0.5 * pow (k_2_sqrtpi / (2. * a), dim) 
    * (exp (-sum1 / (a * a)) + exp (-sum2 / (a * a)));
}

/* Tsuda's example */
double f3 (size_t dim, const double *x, double c) {
  double prod = 1.;
  size_t i;
  for (i = 0; i < dim; i++) {
    prod *= c / (c + 1) * pow((c + 1) / (c + x[i]), 2.0);
  }
  return prod;
}

/* test integrand from W. J. Morokoff and R. E. Caflisch, "Quasi=
   Monte Carlo integration," J. Comput. Phys 122, 218-230 (1995).
   Designed for integration on [0,1]^dim, integral = 1. */
double morokoff(size_t dim, const double *x) {
  double p = 1.0 / dim;
  double prod = pow(1 + p, dim);
  size_t i;
  for (i = 0; i < dim; i++) {
    prod *= pow(x[i], p);
  }
  return prod;
}

int f_test(size_t dim, const double *x, size_t fdim, double *retval) {
  
  double val;
  ++cub_count;

  double fdata;
  if (which_integrand==6) {
    fdata=(1.0+sqrt(10.0))/9.0;
  } else {
    fdata=0.1;
  }
  switch (which_integrand) {
  case 0:
    /* simple smooth (separable) objective: prod. cos(x[i]). */
    val = 1;
    for (size_t i = 0; i < dim; ++i) {
      val *= cos(x[i]);
    }
    break;
  case 1:
    {
      /* integral of exp(-x^2), rescaled to (0,infinity) limits */
      double scale = 1.0;
      val = 0;
      for (size_t i = 0; i < dim; ++i) {
	if (x[i] > 0) {
	  double z = (1 - x[i]) / x[i];
	  val += z * z;
	  scale *= k_2_sqrtpi / (x[i] * x[i]);
	} else {
	  scale = 0;
	  break;
	}
      }
      val = exp(-val) * scale;
      break;
    }
  case 2: /* discontinuous objective: volume of hypersphere */
    val = 0;
    for (size_t i = 0; i < dim; ++i) {
      val += x[i] * x[i];
    }
    val = val < radius * radius;
    break;
  case 3:
    val = f0(dim, x);
    break;
  case 4:
    val = f1(dim, x, fdata);
    break;
  case 5:
    val = f2(dim, x, fdata);
    break;
  case 6:
    val = f3(dim, x, fdata);
    break;
  case 7:
    val = morokoff(dim, x);
    break;
  default:
    cout << "Unknown integrand." << endl;
    exit(-1);
  }

  retval[0] = val;
  return 0;
}

/* surface area of n-dimensional unit hypersphere */
double S(size_t n) {
  double val;
  int fact = 1;
  if (n % 2 == 0) {
    /* n even */
    val = 2 * pow(o2scl_const::pi, n * 0.5);
    n = n / 2;
    while (n > 1) fact *= (n -= 1);
    val /= fact;
  } else {
    /* n odd */
    val = (1 << (n/2 + 1)) * pow(o2scl_const::pi, n/2);
    while (n > 2) fact *= (n -= 2);
    val /= fact;
  }
  return val;
}

template<class vec_t>
double exact_integral(int which, size_t dim,
		      const vec_t &xmax) {
  double val;
  switch(which) {
  case 0:
    val = 1;
    for (size_t i = 0; i < dim; ++i) {
      val *= sin(xmax[i]);
    }
    break;
  case 2:
    val = dim == 0 ? 1 : S(dim) * pow(radius * 0.5, dim) / dim;
    break;
  default:
    val = 1.0;
  }
  return val;
}

int fv(size_t ndim, size_t npt, const double *x, size_t fdim,
       double *fval) {
  for (size_t i = 0; i < npt; ++i) {
    if (f_test(ndim, x + i*ndim, fdim, fval + i*fdim)) {
      return o2scl::gsl_failure;
    }
  }
  return o2scl::success;
}

int fv_new(size_t ndim, size_t npt, const ubvector_crange &x, size_t fdim,
	   ubvector_range &fval) {
  double *x2=new double[ndim];
  double *f2=new double[fdim];
  for (size_t i = 0; i < npt; ++i) {
    const double *xp=&(x[i*ndim]);
    vector_copy(ndim,xp,x2);
    if (f_test(ndim,x2,fdim,f2)) {
      return o2scl::gsl_failure;
    }
    double *fp=&(fval[i*fdim]);
    vector_copy(fdim,f2,fp);
  }
  delete[] x2;
  delete[] f2;
  return 0;
}

/** Test integrating a few functions at once
 */
int fv2(size_t ndim, size_t npt, const double *x, size_t fdim,
	double *fval) {
  for (size_t i=0;i<npt;i++) {
    const double *x2=x+i*ndim;
    double *f2=fval+i*fdim;
    f2[0]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)));
    f2[1]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)))*x2[0]*x2[0];
    f2[2]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)))*x2[0]*x2[0]*x2[1]*x2[1];
  }
  return 0;
}

int fv2_new(size_t ndim, size_t npt, const ubvector_crange &x, size_t fdim,
	    ubvector_range &fval) {
  
  for (size_t i=0;i<npt;i++) {
    const ubvector_crange_crange x2=const_vector_range(x,i*ndim,(i+1)*ndim);
    ubvector_range_range f2=vector_range(fval,i*fdim,(i+1)*fdim);
    f2[0]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)));
    f2[1]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)))*x2[0]*x2[0];
    f2[2]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)))*x2[0]*x2[0]*x2[1]*x2[1];
  }
  return 0;
}

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr tmgr;
  tmgr.set_output_level(1);

  size_t dim=3;
  ubvector xmin(3), xmax(3);
  for (size_t i=0;i<dim;++i) {
    xmin[i]=0.0;
    xmax[i]=1.0;
  }

  typedef std::function<
    int(size_t,size_t,const double *,size_t,double *)> cub_funct_arr;
  typedef std::function<
    int(size_t,size_t,const ubvector_crange &,size_t,
	ubvector_range &)> cub_funct_ub;
  inte_hcubature<cub_funct_arr> hc;
  inte_pcubature<cub_funct_arr,ubvector,ubvector_crange,
		 ubvector_range> pc;

  inte_cubature_base::error_norm en=inte_cubature_base::ERROR_INDIVIDUAL;
  
  cub_funct_arr cfa=fv;
  cub_funct_ub cfa_new=fv_new;

  // Test both hcubature and pcubature with several integrands
  // and compare with original cubature testing results
  
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

    double tol;
    size_t maxEval;
    ubvector vval(1), verr(1);

    tol=1.0e-2;
    maxEval=0;
    
    which_integrand = test_iand; 
    
    if (test_iand!=2) {

      cub_count=0;
      hc.integ(1,cfa,dim,xmin,xmax,maxEval,0,tol,en,vval,verr);
	       
      cout << "# " << which_integrand << " " 
	   << "integral " << vval[0] << " "
	   << "est. error " << verr[0] << " " 
	   << "true error " 
	   << fabs(vval[0]-exact_integral(which_integrand,dim,xmax)) << endl;
      cout << "evals " << cub_count << endl;

      tmgr.test_gen(fabs(vval[0]-exact_integral(which_integrand,dim,xmax))<
		    verr[0]*2.0,"hcub 2");
      tmgr.test_gen(test_n[tcnt]==cub_count,"cub_count");
      tmgr.test_rel(vval[0],test_vals[tcnt][0],5.0e-6,"val");
      tmgr.test_rel(verr[0],test_vals[tcnt][1],5.0e-6,"err");
      tmgr.test_rel(fabs(vval[0]-exact_integral(which_integrand,dim,xmax)),
		    test_vals[tcnt][2],5.0e-6,"diff w/ exact");
      tcnt++;
    }
    
    if (test_iand!=3) {

      cub_count=0;
      pc.integ(1,cfa,dim,xmin,xmax,maxEval,0,tol,en,vval,verr);
	       
      cout << "# " << which_integrand << " " << "integral " << vval[0]
	   << " " << "est. error " << verr[0] << " " << "true error " 
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
    }
    
  }

  // With parallelism for hcubature
  hc.use_parallel=1;

  if (false) {

    tcnt=0;
    for(size_t test_iand=0;test_iand<8;test_iand++) {

      double tol;
      size_t maxEval;
      ubvector vval(1), verr(1);

      tol=1.0e-2;
      maxEval=0;
    
      which_integrand = test_iand; 
    
      if (test_iand!=2) {

	cub_count=0;
	hc.integ(1,cfa,dim,xmin,xmax,maxEval,0,tol,en,vval,verr);
	       
	cout << "# " << which_integrand << " " 
	     << "integral " << vval[0] << " "
	     << "est. error " << verr[0] << " " 
	     << "true error " 
	     << fabs(vval[0]-exact_integral(which_integrand,dim,xmax)) << endl;
	cout << "evals " << cub_count << endl;

	tmgr.test_gen(fabs(vval[0]-exact_integral(which_integrand,dim,xmax))<
		      verr[0]*2.0,"hcub 2");
	tmgr.test_gen(test_n[tcnt]==cub_count,"cub_count");
	tmgr.test_rel(vval[0],test_vals[tcnt][0],5.0e-6,"val");
	tmgr.test_rel(verr[0],test_vals[tcnt][1],5.0e-6,"err");
	tmgr.test_rel(fabs(vval[0]-exact_integral(which_integrand,dim,xmax)),
		      test_vals[tcnt][2],5.0e-6,"diff w/ exact");
	tcnt++;
      }
    
      if (test_iand!=3) {
	tcnt++;
      }
    
    }

  }

  // Now run with dim=1 (without parallelism)

  hc.use_parallel=0;

  int test_n2[14]={15,5,75,17,257,15,45,33,105,33,15,9,15,3};
  
  double test_vals2[14][3]={{8.414710e-01,9.342205e-15,0.000000e+00},
			    {8.414712e-01,3.009270e-04,1.804358e-07},
			    {1.000000e+00,2.435248e-06,6.890044e-13},
			    {9.999706e-01,1.113513e-03,2.937509e-05},
			    {5.030680e-01,3.067965e-03,1.826509e-03},
			    {1.000000e+00,1.110223e-14,0.000000e+00},
			    {1.000000e+00,1.225522e-04,4.218679e-10},
			    {1.000000e+00,2.155111e-04,4.414323e-09},
			    {9.999988e-01,2.820725e-07,1.214234e-06},
			    {9.999988e-01,7.176009e-05,1.216743e-06},
			    {1.000000e+00,8.919302e-07,8.393286e-14},
			    {9.999983e-01,1.425694e-03,1.729318e-06},
			    {1.000000e+00,1.110223e-14,0.000000e+00},
			    {1.000000e+00,0.000000e+00,0.000000e+00}};
  
  dim=1;
  tcnt=0;
  
  for(size_t test_iand=0;test_iand<8;test_iand++) {

    double tol;
    size_t maxEval;
    ubvector vval(1), verr(1);

    tol=1.0e-2;
    maxEval=0;
    
    which_integrand = test_iand; 
    
    if (test_iand!=2) {

      cub_count=0;
      hc.integ(1,cfa,dim,xmin,xmax,maxEval,0,tol,en,vval,verr);
	       
      cout << "# " << which_integrand << " " 
	   << "integral " << vval[0] << " "
	   << "est. error " << verr[0] << " " 
	   << "true error " 
	   << fabs(vval[0]-exact_integral(which_integrand,dim,xmax)) << endl;
      cout << "evals " << cub_count << endl;

      if (test_iand!=5) {
	tmgr.test_gen(fabs(vval[0]-exact_integral(which_integrand,dim,xmax))<
		      verr[0]*2.0,"hcub 2");
      }
      tmgr.test_gen(test_n2[tcnt]==cub_count,"cub_count");
      tmgr.test_rel(vval[0],test_vals2[tcnt][0],5.0e-6,"val");
      tmgr.test_rel(verr[0],test_vals2[tcnt][1],5.0e-6,"err");
      tmgr.test_rel(fabs(vval[0]-exact_integral(which_integrand,dim,xmax)),
		    test_vals2[tcnt][2],5.0e-6,"diff w/ exact");
      tcnt++;
    }

    if (test_iand!=3) {

      cub_count=0;
      pc.integ(1,cfa,dim,xmin,xmax,maxEval,0,tol,en,vval,verr);
	       
      cout << "# " << which_integrand << " " 
	   << "integral " << vval[0] << " " << "est. error "
	   << verr[0] << " " << "true error " 
	   << fabs(vval[0]-exact_integral(which_integrand,dim,xmax)) << endl;
      cout << "evals " << cub_count << endl;

      if (test_iand!=7) {
	tmgr.test_gen(fabs(vval[0]-exact_integral(which_integrand,dim,xmax))<
		      verr[0]*2.0,"pcub 2");
      }
      tmgr.test_gen(test_n2[tcnt]==cub_count,"cub_count");
      tmgr.test_rel(vval[0],test_vals2[tcnt][0],5.0e-6,"val");
      tmgr.test_rel(verr[0],test_vals2[tcnt][1],5.0e-6,"err");
      tmgr.test_rel(fabs(vval[0]-exact_integral(which_integrand,dim,xmax)),
		    test_vals2[tcnt][2],5.0e-6,"diff w/ exact");
      tcnt++;
    }
    
  }

  // Run another test which tests integrating more than one
  // function at a time
  
  {
    vector<double> vlow(2), vhigh(2);
    ubvector vlow2(2), vhigh2(2);
    vlow[0]=-2.0;
    vlow[1]=-2.0;
    vhigh[0]=2.0;
    vhigh[1]=2.0;
    vlow2[0]=-2.0;
    vlow2[1]=-2.0;
    vhigh2[0]=2.0;
    vhigh2[1]=2.0;
    vector<double> dres(3), derr(3);
    ubvector dres2(3), derr2(3);
    cub_funct_arr cfa2=fv2;
    cub_funct_ub cfa2_new=fv2_new;
    int ret=hc.integ(3,cfa2,2,vlow,vhigh,10000,0.0,1.0e-4,en,dres,derr);
    tmgr.test_gen(ret==0,"hc mdim ret");
    tmgr.test_rel(3.067993,dres[0],1.0e-6,"hc mdim val 0");
    tmgr.test_rel(1.569270,dres[1],1.0e-6,"hc mdim val 1");
    tmgr.test_rel(1.056968,dres[2],1.0e-6,"hc mdim val 2");
    ret=pc.integ(3,cfa2,2,vlow2,vhigh2,10000,0.0,1.0e-4,en,dres2,derr2);
    tmgr.test_gen(ret==0,"pc mdim ret");
    tmgr.test_rel(3.067993,dres2[0],1.0e-6,"pc mdim val 0");
    tmgr.test_rel(1.569270,dres2[1],1.0e-6,"pc mdim val 1");
    tmgr.test_rel(1.056968,dres2[2],1.0e-6,"pc mdim val 2");
  }
    
  tmgr.report();
  return 0;
}
