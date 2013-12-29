/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/multi_funct.h>
#include <o2scl/mmin_conp.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

// Test functions for ubvectors
double minfun(size_t n, const ubvector &x) {
  double ret;
  ret=x[0]*x[0]+(x[1]-2.0)*(x[1]-2.0)+3.0;
  return ret;
}

int minfund(size_t n, ubvector &x, ubvector &g) {
  double ret;
  g[0]=2.0*x[0];
  g[1]=2.0*(x[1]-2.0);
  return 0;
}

// Test functions in GSL form
double my_f(const gsl_vector *v, void *pa) {
  double x=gsl_vector_get(v,0);
  double y=gsl_vector_get(v,1);
  return x*x+(y-2.0)*(y-2.0)+3.0;
}

void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
  double x=gsl_vector_get(v,0);
  double y=gsl_vector_get(v,1);
  
  gsl_vector_set(df,0,2.0*x);
  gsl_vector_set(df,1,2.0*(y-2.0));

}

void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df) {
  *f=my_f(x,params);
  my_df(x,params,df);
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
  ubvector x(2);
  double xa[2];
  int vp=0;
  int ret;

  mmin_conp<> g;

  cout.setf(ios::scientific);

#ifdef O2SCL_CPP11
  multi_funct11 mf=minfun;
  grad_funct11 mfd=minfund;
#else
  multi_funct_fptr<> mf(minfun);
  grad_funct_fptr<> mfd(minfund);
#endif

  // Normal version without gradient
  
  x[0]=1.0;
  x[1]=1.0;
  ret=g.mmin(2,x,min,mf);
  cout << g.last_ntrial << endl;
  cout << ret << " " << min << " " << x[0] << " " << x[1] << endl;
  t.test_gen(ret==0,"mmin_conp 0");
  t.test_abs(x[0],0.0,1.0e-4,"mmin_conp 1");
  t.test_rel(x[1],2.0,2.0e-4,"mmin_conp 2");
  cout << endl;

  // Show that we can use a user-specified automatic gradient object
#ifndef O2SCL_CPP11
  mmin_conp<multi_funct<>,ubvector,grad_funct<ubvector>,
	    gradient<multi_funct<>,ubvector>,
	    gradient_gsl_new<multi_funct<>,ubvector> > g3;
#else
mmin_conp<multi_funct11,ubvector,grad_funct<ubvector>,
	  gradient<multi_funct11,ubvector>,
	  gradient_gsl_new<multi_funct11,ubvector> > g3;
#endif

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
  cout << ret << " " << min << " " << x[0] << " " << x[1] << endl;
  t.test_gen(ret==0,"mmin_conp grad 0");
  t.test_abs(x[0],0.0,1.0e-4,"mmin_conp grad 1");
  t.test_rel(x[1],2.0,2.0e-4,"mmin_conp grad 2");
  cout << endl;

  // Normal version with gradient

  x[0]=1.0;
  x[1]=1.0;
  ret=g.mmin_de(2,x,min,mf,mfd);
  cout << g.last_ntrial << endl;
  cout << ret << " " << min << " " << x[0] << " " << x[1] << endl;
  t.test_gen(ret==0,"mmin_conp grad 0");
  t.test_abs(x[0],0.0,1.0e-4,"mmin_conp grad 1");
  t.test_rel(x[1],2.0,2.0e-4,"mmin_conp grad 2");
  cout << endl;

  // GSL version with gradient
  {
    int iter=0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    gsl_vector *gx;
    gsl_multimin_function_fdf my_func;
    
    my_func.f=&my_f;
    my_func.df=&my_df;
    my_func.fdf=&my_fdf;
    my_func.n=2;
    my_func.params=0;
    
    gx=gsl_vector_alloc (2);
    gsl_vector_set(gx,0,1.0);
    gsl_vector_set(gx,1,1.0);
    T=gsl_multimin_fdfminimizer_conjugate_pr;
    s=gsl_multimin_fdfminimizer_alloc (T, 2);
    
    gsl_multimin_fdfminimizer_set(s,&my_func,gx,0.01,1e-4);
    
    do {
      iter++;
      status=gsl_multimin_fdfminimizer_iterate(s);
      /*
	cout << "Iter: " << iter << " " << gsl_vector_get(s->x,0) << " "
	<< gsl_vector_get(s->x,1) << " " 
	<< gsl_vector_get(s->gradient,0) << " "
	<< gsl_vector_get(s->gradient,1) << " "
	<< endl;
      */
	
      if (status)
	break;
      
      status=gsl_multimin_test_gradient(s->gradient,1e-4);
      
    } while (status == GSL_CONTINUE && iter < 100);
    
    cout << iter << endl;
    cout << status << " " << s->f << " "
	 << gsl_vector_get(s->x,0) << " "
	 << gsl_vector_get(s->x,1) << endl;
    t.test_rel(s->f,min,1.0e-12,"conp compare 5");
    t.test_abs(gsl_vector_get(s->x,0),x[0],1.0e-15,"conp compare 6");
    t.test_rel(gsl_vector_get(s->x,1),x[1],1.0e-12,"conp compare 7");
    t.test_gen(iter==g.last_ntrial,"conp compare 8");
    
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(gx);
  }

  t.report();
  return 0;
}
 
