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
#include <o2scl/multi_funct.h>
#include <o2scl/mmin_bfgs2.h>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

// Test functions for ubvectors
double minfun(size_t n, const ubvector &x) {
  double ret;
  ret=x[0]*x[0]+(x[1]-2.0)*(x[1]-2.0)+3.0;
  return(ret);
}

int minfund(size_t n, ubvector &x, ubvector &g) {
  double ret;
  g[0]=2.0*x[0];
  g[1]=2.0*(x[1]-2.0);
  return 0;
}

// Test functions in GSL form
double my_f(const gsl_vector *v, void *pa) {
  double x = gsl_vector_get(v,0);
  double y = gsl_vector_get(v,1);
  return x*x+(y-2.0)*(y-2.0)+3.0;
}

void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
  double x = gsl_vector_get(v,0);
  double y = gsl_vector_get(v,1);
  
  gsl_vector_set(df,0,2.0*x);
  gsl_vector_set(df,1,2.0*(y-2.0));

}

void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df) {
  *f = my_f(x,params);
  my_df(x,params,df);
}

// Updated spring function
double spring_two(size_t nv, const ubvector &x) {
  double theta=atan2(x[1],x[0]);
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double z=x[2];
  double tmz=theta-z;
  double fact=8.0-pow(sin(tmz+o2scl_const::pi/2.0)+1.0,3.0);
  double rm1=r-1.0;
  double ret=fact+exp(rm1*rm1)+z*z/(30.0);
  return ret;
}

// Gradient of the spring function
int sgrad(size_t nv, ubvector &x, ubvector &g) {

  double theta=atan2(x[1],x[0]);
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double z=x[2];
  double tmz=theta-z;
  double rm1=r-1.0;
  double fact=8.0-pow(sin(tmz+o2scl_const::pi/2.0)+1.0,3.0);

  double dtdx=-x[1]/r/r;
  double dtdy=x[0]/r/r;
  double drdx=x[0]/r;
  double drdy=x[1]/r;
  double dfdt=-3.0*pow(sin(tmz+o2scl_const::pi/2.0)+1.0,2.0)*
    cos(tmz+o2scl_const::pi/2.0);
  double dfdz=2.0*z/(30.0)+3.0*pow(sin(tmz+o2scl_const::pi/2.0)+1.0,2.0)*
    cos(tmz+o2scl_const::pi/2.0);
  double dfdr=2.0*rm1*exp(rm1*rm1);

  g[0]=dfdr*drdx+dfdt*dtdx;
  g[1]=dfdr*drdy+dfdt*dtdy;
  g[2]=dfdz;

  return 0;
}

// Test functions in GSL form
double my_f2(const gsl_vector *v, void *pa) {
  double theta=atan2(gsl_vector_get(v,1),gsl_vector_get(v,0));
  double r=sqrt(gsl_vector_get(v,0)*gsl_vector_get(v,0)+
		gsl_vector_get(v,1)*gsl_vector_get(v,1));
  double z=gsl_vector_get(v,2);
  double tmz=theta-z;
  double fact=8.0-pow(sin(tmz+o2scl_const::pi/2.0)+1.0,3.0);
  double rm1=r-1.0;
  double ret=fact+exp(rm1*rm1)+z*z/(30.0);
  return ret;
}

void my_df2(const gsl_vector *v, void *params, gsl_vector *df) {
  double theta=atan2(gsl_vector_get(v,1),gsl_vector_get(v,0));
  double r=sqrt(gsl_vector_get(v,0)*gsl_vector_get(v,0)+
		gsl_vector_get(v,1)*gsl_vector_get(v,1));
  double z=gsl_vector_get(v,2);
  double tmz=theta-z;
  double rm1=r-1.0;
  double fact=8.0-pow(sin(tmz+o2scl_const::pi/2.0)+1.0,3.0);

  double dtdx=-gsl_vector_get(v,1)/r/r;
  double dtdy=gsl_vector_get(v,0)/r/r;
  double drdx=gsl_vector_get(v,0)/r;
  double drdy=gsl_vector_get(v,1)/r;
  double dfdt=-3.0*pow(sin(tmz+o2scl_const::pi/2.0)+1.0,2.0)*
    cos(tmz+o2scl_const::pi/2.0);
  double dfdz=2.0*z/(30.0)+3.0*pow(sin(tmz+o2scl_const::pi/2.0)+1.0,2.0)*
    cos(tmz+o2scl_const::pi/2.0);
  double dfdr=2.0*rm1*exp(rm1*rm1);

  gsl_vector_set(df,0,dfdr*drdx+dfdt*dtdx);
  gsl_vector_set(df,1,dfdr*drdy+dfdt*dtdy);
  gsl_vector_set(df,2,dfdz);

  return;
}

void my_fdf2(const gsl_vector *x, void *params, double *f, gsl_vector *df) {
  *f = my_f2(x,params);
  my_df2(x,params,df);
}

// This demonstrates that the user can specify their own automatic
// gradient class
template<class func_t, class vec_t> class gradient_gsl_new :
  public gradient_gsl<func_t,vec_t> {
public:
  virtual ~gradient_gsl_new() {}
  virtual int operator()(size_t nv, vec_t &x, vec_t &g) {
    std::cout << '.' << std::flush;
    return gradient_gsl<func_t,vec_t>::operator()(nv,x,g);
  }
};

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  double min;
  ubvector x(2), x3(3);
  double xa[2], xa3[3];
  int vp=0;
  int ret;
  mmin_bfgs2<> g;

  // BFGS has trouble converging (especially to zero, since the
  // minimimum of x[0] is exactly at zero) if the derivative is not
  // very accurate.
  g.def_grad.epsrel=1.0e-8;

  cout.setf(ios::scientific);

  multi_funct mf=minfun;
  grad_funct mfd=minfund;
  multi_funct mfs=spring_two;
  grad_funct mfds=sgrad;

  // Normal version 
  
  x[0]=1.0;
  x[1]=1.0;
  ret=g.mmin(2,x,min,mf);
  cout << g.last_ntrial << endl;
  cout << ret << " " << g.last_ntrial << " " << min << " " 
       << x[0] << " " << x[1] << endl;
  t.test_gen(ret==0,"mmin_bfgs2 0");
  t.test_abs(x[0],0.0,1.0e-4,"mmin_bfgs2 1");
  t.test_rel(x[1],2.0,2.0e-4,"mmin_bfgs2 2");
  cout << endl;
  
  // Show that we can use a user-specified automatic gradient object
  mmin_bfgs2<multi_funct,ubvector,grad_funct,
    gradient<multi_funct,ubvector>,
	   gradient_gsl_new<multi_funct,ubvector> > g3;
 
// Emacs has trouble with tabification
#ifdef O2SCL_NEVER_DEFINED
} {
#endif
  
  g3.def_grad.epsrel=1.0e-8;

  // Normal version with user auto gradient

  x[0]=1.0;
  x[1]=1.0;
  ret=g3.mmin(2,x,min,mf);
  cout << g3.last_ntrial << endl;
  cout << ret << " " << g3.last_ntrial << " " << min << " " 
       << x[0] << " " << x[1] << endl;
  t.test_gen(ret==0,"mmin_bfgs2 grad 0");
  t.test_abs(x[0],0.0,1.0e-4,"mmin_bfgs2 grad 1");
  t.test_rel(x[1],2.0,2.0e-4,"mmin_bfgs2 grad 2");
  cout << endl;

  // Normal version with gradient

  x[0]=1.0;
  x[1]=1.0;
  ret=g.mmin_de(2,x,min,mf,mfd);
  cout << g.last_ntrial << endl;
  cout << ret << " " << g.last_ntrial << " " << min << " " 
       << x[0] << " " << x[1] << endl;
  t.test_gen(ret==0,"mmin_bfgs2 grad 0");
  t.test_abs(x[0],0.0,1.0e-4,"mmin_bfgs2 grad 1");
  t.test_rel(x[1],2.0,2.0e-4,"mmin_bfgs2 grad 2");
  cout << endl;

  // GSL version with gradient
  {
    int iter = 0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    gsl_vector *gx;
    gsl_multimin_function_fdf my_func;
    
    my_func.f = &my_f;
    my_func.df = &my_df;
    my_func.fdf = &my_fdf;
    my_func.n = 2;
    my_func.params = 0;
    
    gx = gsl_vector_alloc (2);
    gsl_vector_set(gx,0,1.0);
    gsl_vector_set(gx,1,1.0);
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, 2);
    
    gsl_multimin_fdfminimizer_set(s,&my_func,gx,0.01,1e-4);
    
    do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);

      if (status)
	break;
      
      status = gsl_multimin_test_gradient(s->gradient,1e-4);

    } while (status == GSL_CONTINUE && iter < 100);
    
    cout << iter << endl;
    cout << status << " " << s->f << " "
	 << gsl_vector_get(s->x,0) << " "
	 << gsl_vector_get(s->x,1) << endl;
    t.test_rel(s->f,min,1.0e-12,"bfgs2 compare 5");
    t.test_abs(gsl_vector_get(s->x,0),x[0],1.0e-15,"bfgs2 compare 6");
    t.test_rel(gsl_vector_get(s->x,1),x[1],1.0e-12,"bfgs2 compare 7");
    t.test_gen(iter==g.last_ntrial,"bfgs2 compare 8");
    
    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free(gx);
  }
  cout << endl;

  // Normal version of spring function with numerical gradient
  
  g.def_grad.epsrel=1.0e-8;

  x3[0]=1.0;
  x3[1]=0.0;
  x3[2]=7.0*o2scl_const::pi;
  ret=g.mmin(3,x3,min,mfs);
  cout << g.last_ntrial << endl;
  cout << ret << " " << min << " " 
       << x3[0] << " " << x3[1] << " "<< x3[2] << endl;
  cout << endl;

  // Normal version of spring function with gradient
  // with a slightly different initial guess. 
  // Both this, and the GSL version below, fail in the same way
  
  g.err_nonconv=false;
  x3[0]=1.0;
  x3[1]=1.0;
  x3[2]=7.0*o2scl_const::pi;
  ret=g.mmin_de(3,x3,min,mfs,mfds);
  cout << g.last_ntrial << endl;
  cout << ret << " " << min << " " 
       << x3[0] << " " << x3[1] << " "<< x3[2] << endl;
  cout << endl;

  // GSL version of spring function with gradient
  {
    int iter = 0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    gsl_vector *gx;
    gsl_multimin_function_fdf my_func;
    
    my_func.f = &my_f2;
    my_func.df = &my_df2;
    my_func.fdf = &my_fdf2;
    my_func.n = 3;
    my_func.params = 0;
    
    gx = gsl_vector_alloc (3);
    gsl_vector_set(gx,0,1.0);
    gsl_vector_set(gx,1,1.0);
    gsl_vector_set(gx,2,7.0*o2scl_const::pi);
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, 3);
    
    gsl_multimin_fdfminimizer_set(s,&my_func,gx,0.01,1e-4);
    
    do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);

      if (status)
	break;
      
      status = gsl_multimin_test_gradient(s->gradient,1e-4);
      
    } while (status == GSL_CONTINUE && iter < 100);
    
    cout << iter << endl;
    cout << status << " " << s->f << " "
	 << gsl_vector_get(s->x,0) << " "
	 << gsl_vector_get(s->x,1) << " "
	 << gsl_vector_get(s->x,2) << endl;
    
    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free(gx);
  }

  t.report();
  return 0;
}
 
