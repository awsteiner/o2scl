/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// For gsl_blas_dnrm2
#include <gsl/gsl_blas.h>

#include <o2scl/test_mgr.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/fit_nonlin.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

struct data {
  size_t n;
  double *y;
  double *sigma;
};

int expb_f(const gsl_vector *x, void *params, gsl_vector *f) {
  
  size_t n=((struct data *)params)->n;
  double *y=((struct data *)params)->y;
  double *sigma=((struct data *) params)->sigma;
  
  double A=gsl_vector_get(x,0);
  double lambda=gsl_vector_get(x,1);
  double b=gsl_vector_get(x,2);

  /* Model Yi=A*exp(-lambda*i)+b */
  for (size_t i=0;i<n;i++) {
    double Yi=A*exp(-lambda*((double)i))+b;
    gsl_vector_set(f,i,(Yi-y[i])/sigma[i]);
  }
  
  return success;
}

int expb_df(const gsl_vector *x, void *params, gsl_matrix *J) {

  size_t n=((struct data *)params)->n;
  double *y=((struct data *)params)->y;
  double *sigma=((struct data *) params)->sigma;

  size_t i;
  double A=gsl_vector_get(x,0);
  double lambda=gsl_vector_get(x,1);

  for (i=0;i<n;i++) {
    double t=i;
    double s=sigma[i];
    double e=exp(-lambda*t);
    gsl_matrix_set(J,i,0,e/s);
    gsl_matrix_set(J,i,1,-t*A*e/s);
    gsl_matrix_set(J,i,2,1/s);
  }

  return success;
}

int expb_df2(const gsl_vector *x, void *params, gsl_matrix *J) {

  size_t n=((struct data *)params)->n;
  double *y=((struct data *)params)->y;
  double *sigma=((struct data *) params)->sigma;

  gsl_vector *f=gsl_vector_alloc(n);
  
  double eps=1.0e-4;
  //double eps=GSL_SQRT_DBL_EPSILON;

  for(size_t j=0;j<3;j++) {

    double xtemp=gsl_vector_get(x,j);
    double step=fabs(xtemp)*eps;
    if (step<1.0e-15) step=eps;

    for(size_t i=0;i<n;i++) {

      gsl_vector_set((gsl_vector *)x,j,xtemp+step);
      expb_f(x,params,f);
      double yhi=gsl_vector_get(f,i);

      gsl_vector_set((gsl_vector *)x,j,xtemp);
      expb_f(x,params,f);
      double ylo=gsl_vector_get(f,i);

      gsl_matrix_set(J,i,j,(yhi-ylo)/step);
    }
  }

  gsl_vector_free(f);
    
  return success;
}

int expb_fdf(const gsl_vector *x, void *params, gsl_vector *f, 
	     gsl_matrix *J) {
  expb_f(x,params,f);
  expb_df(x,params,J);
  return success;
}

int expb_fdf2(const gsl_vector *x, void *params, gsl_vector *f, 
	      gsl_matrix *J) {
  expb_f(x,params,f);
  expb_df2(x,params,J);
  return success;
}

double func(size_t np, const ubvector &p, double x) {
  return p[0]*exp(-p[1]*x)+p[2];
}

int main(void) {
  test_mgr tm;
  tm.set_output_level(1);

  cout.setf(ios::scientific);
  cout.precision(8);
  
  std::vector<double> x0_s, x1_s, x2_s;
  std::vector<double> x0_u, x1_u, x2_u;
  double chi2red_s, chi2red_u;

  /* This is the data to be fitted */

  prob_dens_gaussian gauss(0.0,0.1);
  const size_t n=40;
  double xdat[40], y[40], sigma[40];
  for (size_t i=0;i<n;i++) {
    xdat[i]=((double)i);
    y[i]=1.0+5*exp(-0.1*((double)i))+gauss.sample();
    sigma[i]=0.1;
  }
  struct data d={n,y,sigma};

  // A copy in ubvectors 
  ubvector axdat(40);
  ubvector ay(40);
  ubvector asigma(40);
  vector_copy(40,xdat,axdat);
  vector_copy(40,y,ay);
  vector_copy(40,sigma,asigma);
  
  fit_funct11 fff=func;
  chi_fit_funct<> cff(40,axdat,ay,asigma,fff);
  //cff.auto_jac.epsrel=GSL_SQRT_DBL_EPSILON;
  cff.auto_jac.epsrel=1.0e-4;

  // Copy of GSL covariance matrix for testing
  gsl_matrix *covar1=gsl_matrix_alloc(3,3);

  //----------------------------------------------------------------
  // GSL, scaled version

  if (true) {

    cout << "GSL scaled" << endl;

    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;

    int status;
    size_t i, iter=0;

    const size_t p=3;

    gsl_multifit_function_fdf f;

    double x_init[3]={ 1.0, 0.0, 0.0 };

    gsl_vector_view x=gsl_vector_view_array (x_init, p);

    f.f=&expb_f;
    f.df=&expb_df2;
    f.fdf=&expb_fdf2;
    f.n=n;
    f.p=p;
    f.params=&d;

    T=gsl_multifit_fdfsolver_lmsder;
    s=gsl_multifit_fdfsolver_alloc(T, n, p);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);

    do {
      iter++;
      status=gsl_multifit_fdfsolver_iterate(s);

      cout << gsl_vector_get(s->x,0) << " ";
      cout << gsl_vector_get(s->x,1) << " ";
      cout << gsl_vector_get(s->x,2) << endl;

      x0_s.push_back(gsl_vector_get(s->x,0));
      x1_s.push_back(gsl_vector_get(s->x,1));
      x2_s.push_back(gsl_vector_get(s->x,2));
      
      if (status)
	break;

      status=gsl_multifit_test_delta(s->dx,s->x,1e-4,1e-4);
				     
    } while (status == gsl_continue && iter < 500);
  
    x0_s.push_back(gsl_vector_get(s->x,0));
    x1_s.push_back(gsl_vector_get(s->x,1));
    x2_s.push_back(gsl_vector_get(s->x,2));

    gsl_multifit_covar(s->J,0.0,covar1);

    double chi=gsl_blas_dnrm2(s->f);
    chi2red_s=chi*chi/(n-p);

    gsl_multifit_fdfsolver_free(s);

  } 
  cout << endl;

  //----------------------------------------------------------------
  // GSL, unscaled version

  if (true) {

    cout << "GSL unscaled" << endl;

    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;

    int status;
    size_t i, iter=0;

    const size_t p=3;

    gsl_matrix *covar=gsl_matrix_alloc (p, p);

    gsl_multifit_function_fdf f;

    double x_init[3] = { 1.0, 0.0, 0.0 };

    gsl_vector_view x = gsl_vector_view_array (x_init, p);

    f.f=&expb_f;
    f.df=&expb_df2;
    f.fdf=&expb_fdf2;
    f.n=n;
    f.p=p;
    f.params=&d;

    T=gsl_multifit_fdfsolver_lmder;
    s=gsl_multifit_fdfsolver_alloc(T, n, p);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);

    do {
      iter++;
      status=gsl_multifit_fdfsolver_iterate(s);

      cout << gsl_vector_get(s->x,0) << " ";
      cout << gsl_vector_get(s->x,1) << " ";
      cout << gsl_vector_get(s->x,2) << endl;
      
      x0_u.push_back(gsl_vector_get(s->x,0));
      x1_u.push_back(gsl_vector_get(s->x,1));
      x2_u.push_back(gsl_vector_get(s->x,2));

      if (status)
	break;
      
      status=gsl_multifit_test_delta(s->dx,s->x,1e-4,1e-4);
      
    } while (status == gsl_continue && iter < 500);
  
    x0_u.push_back(gsl_vector_get(s->x,0));
    x1_u.push_back(gsl_vector_get(s->x,1));
    x2_u.push_back(gsl_vector_get(s->x,2));

    gsl_multifit_covar (s->J, 0.0, covar);

    double chi=gsl_blas_dnrm2(s->f);
    chi2red_u=chi*chi/(n-p);

    gsl_multifit_fdfsolver_free(s);

  } 
  cout << endl;

  //----------------------------------------------------------------
  // O2scl, scaled version

  if (true) {

    cout << "O2scl scaled fit()." << endl;
    
    fit_nonlin<> gf;
    
    fit_funct11 ff=func;

    double chi2;
    ubmatrix mycovar(3,3);

    double x[3]={1.0,0.0,0.0};
    ubvector ax(3);
    vector_copy(3,x,ax);
    
    gf.fit(3,ax,mycovar,chi2,cff);
    
    tm.test_rel(ax[0],x0_s[x0_s.size()-1],1.0e-8,"scaled x0");
    tm.test_rel(ax[1],x1_s[x1_s.size()-1],1.0e-8,"scaled x1");
    tm.test_rel(ax[2],x2_s[x2_s.size()-1],1.0e-8,"scaled x2");
    tm.test_rel(chi2red_s,chi2/(n-3),1.0e-10,"scaled chi2");

    for(size_t i=0;i<3;i++) {
      for(size_t j=0;j<3;j++) {
	tm.test_rel(mycovar(i,j),gsl_matrix_get(covar1,i,j),
		    1.0e-8,"covariance mat.");
      }
    }

    cout << endl;
  }

  //----------------------------------------------------------------
  // O2scl, scaled version with set() and iterate()

  if (true) {

    cout << "O2scl scaled set/iterate()." << endl;

    fit_nonlin<> gf;
  
    fit_funct11 ff=func;

    double chi2;
    ubmatrix mycovar(3,3);

    double x[3]={1.0,0.0,0.0};
    ubvector ax(3);
    vector_copy(3,x,ax);

    gf.set(3,ax,cff);

    for(size_t i=0;i<10;i++) {
      gf.iterate();
      //cout << ax << endl;
      tm.test_rel(ax[0],x0_s[i],1.0e-3,"scaled x0 set/iter");
      tm.test_rel(ax[1],x1_s[i],1.0e-3,"scaled x1 set/iter");
      tm.test_rel(ax[2],x2_s[i],1.0e-3,"scaled x2 set/iter");
      if (gf.test_delta_f(3,gf.dx_,ax,gf.tol_abs,
			  gf.tol_rel)!=gsl_continue) {
	i=10;
      } 
    }
    chi2=o2scl_cblas::dnrm2(n,gf.f_);
    chi2*=chi2;

    tm.test_rel(ax[0],x0_s[x0_s.size()-1],1.0e-8,"post scaled x0 set/iter");
    tm.test_rel(ax[1],x1_s[x1_s.size()-1],1.0e-8,"post scaled x1 set/iter");
    tm.test_rel(ax[2],x2_s[x2_s.size()-1],1.0e-8,"post scaled x2 set/iter");
    tm.test_rel(chi2red_s,chi2/(n-3),1.0e-10,"post scaled chi2");

    cout << endl;
  }

  //----------------------------------------------------------------
  // O2scl, unscaled version

  if (true) {

    cout << "O2scl unscaled fit()." << endl;

    fit_nonlin<> gf;
    gf.use_scaled=false;
    
    fit_funct11 ff=func;

    double chi2;
    ubmatrix mycovar(3,3);

    double x[3]={1.0,0.0,0.0};
    ubvector ax(3);
    vector_copy(3,x,ax);

    gf.fit(3,ax,mycovar,chi2,cff);

    //cout << ax << endl;

    tm.test_rel(ax[0],x0_u[x0_u.size()-1],1.0e-11,"unscaled x0");
    tm.test_rel(ax[1],x1_u[x1_u.size()-1],1.0e-11,"unscaled x1");
    tm.test_rel(ax[2],x2_u[x2_u.size()-1],1.0e-11,"unscaled x2");
    tm.test_rel(chi2red_u,chi2/(n-3),1.0e-10,"unscaled chi2");

    cout << endl;
  }

  //----------------------------------------------------------------
  // O2scl, unscaled version with set() and iterate()

  {

    cout << "O2scl unscaled set/iterate()." << endl;

    fit_nonlin<> gf;
    gf.use_scaled=false;
  
    fit_funct11 ff=func;

    double chi2;
    ubmatrix mycovar(3,3);

    double x[3]={1.0,0.0,0.0};
    ubvector ax(3);
    vector_copy(3,x,ax);

    gf.set(3,ax,cff);

    for(size_t i=0;i<10;i++) {
      gf.iterate();

      //cout << ax << endl;

      tm.test_rel(ax[0],x0_u[i],1.0e-4,"unscaled x0 set/iter");
      tm.test_rel(ax[1],x1_u[i],1.0e-4,"unscaled x1 set/iter");
      tm.test_rel(ax[2],x2_u[i],1.0e-4,"unscaled x2 set/iter");
      if (gf.test_delta_f(3,gf.dx_,ax,gf.tol_abs,
			  gf.tol_rel)!=gsl_continue) {
	i=10;
      } 
    }
    chi2=o2scl_cblas::dnrm2(n,gf.f_);
    chi2*=chi2;
    
    tm.test_rel(ax[0],x0_u[x0_u.size()-1],1.0e-12,"post unscaled x0 set/iter");
    tm.test_rel(ax[1],x1_u[x1_u.size()-1],1.0e-11,"post unscaled x1 set/iter");
    tm.test_rel(ax[2],x2_u[x2_u.size()-1],1.0e-11,"post unscaled x2 set/iter");
    tm.test_rel(chi2red_u,chi2/(n-3),1.0e-10,"post unscaled chi2 set/iter");

    cout << endl;
  }

  tm.report();
  return 0;
}


