/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2016, Andrew W. Steiner

  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*
  This compares the speed of integrating ODE's with the Cash-Karp
  stepper using various approaches. Using arrays is typically faster
  than using ovector objects, and the class gsl_rkck_fast is faster
  than gsl_rkck.  Note the relative performance will change somewhat
  in real world applications because the differential equations will
  be more complex.
*/

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include <o2scl/ovector_tlate.h>
#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/ode_funct.h>
#include <o2scl/gsl_rkck.h>

using namespace std;
using namespace o2scl;

int derivs(double x, size_t nv, const  ovector_base &y, 
	   ovector_base &dydx, int &pa) {
  dydx[0]=y[0];
  return 0;
}

int derivsu(double x, size_t nv, const  uvector_base &y, 
	    uvector_base &dydx, int &pa) {
  dydx[0]=y[0];
  return 0;
}

int derivs_arr(double x, size_t nv, const double y[], 
	       double dydx[], int &pa) {
  dydx[0]=y[0];
  return 0;
}

int derivs_gsl(double x, const double y[], double dydx[], void *pa) {
  dydx[0]=y[0];
  return 0;
}

int main(void) {
  const size_t N=100000;
  int i;
  double x, dx=1.0e-1;
  ovector y(2), dydx(2), yerr(2);
  uvector yu(2), dydxu(2), yerru(2);
  double ya[1], dydxa[1], yerra[1];
  int vp=0;

  gsl_rkck<int,ode_funct_fptr<int >,ovector_base,ovector,
    ovector_alloc> rk;
  gsl_rkck<int,ode_funct_fptr<int,uvector_base>,uvector_base,uvector,
    uvector_alloc> rku;
  gsl_rkck_fast<1,int,ode_funct_fptr<int >,ovector_base,ovector,
    ovector_alloc> rkf;
  gsl_rkck<int,ode_vfunct_fptr<int,1>,
    double [1], double [1],array_alloc<double [1]> > rka;
  gsl_rkck_fast<1,int,ode_vfunct_fptr<int,1>,
    double [1], double [1],array_alloc<double [1]> > rkaf;

  test_mgr t;
  t.set_output_level(1);
  size_t tmp, t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0;

  cout.setf(ios::scientific);

  ode_funct_fptr<int,ovector_base> od(derivs);
  ode_funct_fptr<int,uvector_base> odu(derivsu);
  ode_vfunct_fptr<int,1> od_arr(derivs_arr);
  
  cout << endl;
  {
    // Show that GSL range checking is off
    gsl_vector *gv=gsl_vector_alloc(2);
    double xx=gsl_vector_get(gv,3);
    cout << "Gsl   range checking: 0" << endl;
  }
  cout << "O2scl range checking: " << lib_settings.range_check() << endl;
  cout << endl;
  
  for(size_t ii=0;ii<5;ii++) {

    // With ovector

    tmp=clock();
    for(size_t kk=0;kk<N;kk++) {
      x=1.0;
      y[0]=1.0;
      derivs(x,1,y,dydx,vp);
      for(i=1;i<=40;i++) {
	rk.step(x,dx,1,y,dydx,y,yerr,dydx,vp,od);
	x+=dx;
      }
    }
    t1+=(clock()-tmp)/10000;
    t.test_rel(y[0],exp(x)/exp(1.0),1.0e-6,"y_calculated-y_exact");
    t.test_abs(yerr[0],0.0,1.0e-6,"y_err");

    // Ovector with gsl_rkck_fast

    tmp=clock();
    for(size_t kk=0;kk<N;kk++) {
      x=1.0;
      y[0]=1.0;
      derivs(x,1,y,dydx,vp);
      for(i=1;i<=40;i++) {
	rkf.step(x,dx,1,y,dydx,y,yerr,dydx,vp,od);
	x+=dx;
      }
    }
    t2+=(clock()-tmp)/10000;
    t.test_rel(y[0],exp(x)/exp(1.0),1.0e-6,"y_calculated-y_exact");
    t.test_abs(yerr[0],0.0,1.0e-6,"y_err");
  
    // Arrays (GSL)
  
    gsl_odeiv_step *s=gsl_odeiv_step_alloc(gsl_odeiv_step_rkck,1);
    gsl_odeiv_system sys={derivs_gsl,0,1,0};

    tmp=clock();
    for(size_t kk=0;kk<N;kk++) {
      x=1.0;
      ya[0]=1.0;
      derivs_arr(x,1,ya,dydxa,vp);
      for(i=1;i<=40;i++) {
	gsl_odeiv_step_apply(s,x,dx,ya,yerra,dydxa,dydxa,&sys);
	x+=dx;
      }
    }
    t3+=(clock()-tmp)/10000;
    t.test_rel(y[0],exp(x)/exp(1.0),1.0e-6,"y_calculated-y_exact");
    t.test_abs(yerra[0],0.0,1.0e-6,"y_err");
  
    gsl_odeiv_step_free(s);

    // Arrays

    tmp=clock();
    for(size_t kk=0;kk<N;kk++) {
      x=1.0;
      ya[0]=1.0;
      derivs_arr(x,1,ya,dydxa,vp);
      for(i=1;i<=40;i++) {
	rka.step(x,dx,1,ya,dydxa,ya,yerra,dydxa,vp,od_arr);
	x+=dx;
      }
    }
    t4+=(clock()-tmp)/10000;
    t.test_rel(ya[0],exp(x)/exp(1.0),1.0e-6,"y_calculated-y_exact");
    t.test_abs(yerra[0],0.0,1.0e-6,"y_err");

    // Arrays with gsl_rkck_fast

    tmp=clock();
    for(size_t kk=0;kk<N;kk++) {
      x=1.0;
      ya[0]=1.0;
      derivs_arr(x,1,ya,dydxa,vp);
      for(i=1;i<=40;i++) {
	rkaf.step(x,dx,1,ya,dydxa,ya,yerra,dydxa,vp,od_arr);
	x+=dx;
      }
    }
    t5+=(clock()-tmp)/10000;
    t.test_rel(ya[0],exp(x)/exp(1.0),1.0e-6,"y_calculated-y_exact");
    t.test_abs(yerra[0],0.0,1.0e-6,"y_err");
  
    // Arrays with gsl_rkck_fast and direct function call
  
    typedef int (*der_fn_t)(double x, size_t nv, const double y[], 
			    double dydx[], int &pa);
    der_fn_t fnp=&derivs_arr;
    gsl_rkck_fast<1,int,der_fn_t,
      double [1], double [1],array_alloc<double [1]> > rkaf2;

    tmp=clock();
    for(size_t kk=0;kk<N;kk++) {
      x=1.0;
      ya[0]=1.0;
      derivs_arr(x,1,ya,dydxa,vp);
      for(i=1;i<=40;i++) {
	rkaf2.step(x,dx,1,ya,dydxa,ya,yerra,dydxa,vp,fnp);
	x+=dx;
      }
    }
    t6+=(clock()-tmp)/10000;
    t.test_rel(ya[0],exp(x)/exp(1.0),1.0e-6,"y_calculated-y_exact");
    t.test_abs(yerra[0],0.0,1.0e-6,"y_err");

    // With uvector

    tmp=clock();
    for(size_t kk=0;kk<N;kk++) {
      x=1.0;
      yu[0]=1.0;
      derivsu(x,1,yu,dydxu,vp);
      for(i=1;i<=40;i++) {
	rku.step(x,dx,1,yu,dydxu,yu,yerru,dydxu,vp,odu);
	x+=dx;
      }
    }
    t7+=(clock()-tmp)/10000;
    t.test_rel(yu[0],exp(x)/exp(1.0),1.0e-6,"y_calculated-y_exact");
    t.test_abs(yerru[0],0.0,1.0e-6,"y_err");

  }
  
  cout.width(20);
  cout << "Ovectors: " << t1 << endl;
  cout.width(20);
  cout << "Ovectors (fast): " << t2 << endl;
  cout.width(20);
  cout << "GSL: " << t3 << endl;
  cout.width(20);
  cout << "Arrays: " << t4 << endl;
  cout.width(20);
  cout << "Arrays (fast): " << t5 << endl;
  cout.width(20);
  cout << "Arrays (direct fn): " << t6 << endl;
  cout.width(20);
  cout << "Uvectors: " << t7 << endl;
  cout << endl;

  t.report();

  return 0;
}
