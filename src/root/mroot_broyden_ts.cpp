/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#include <cmath>

#include <gsl/gsl_multiroots.h>

#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>

#include <o2scl/mroot_broyden.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int gsl_fn(const gsl_vector *x, void *pa, gsl_vector *f) {
  gsl_vector_set(f,0,sin(gsl_vector_get(x,1)-0.2));
  gsl_vector_set(f,1,sin(gsl_vector_get(x,0)-0.25));
  return 0;
}

int gfn(size_t nv, const ubvector &x, ubvector &y) {
  y[0]=sin(x[1]-0.2);
  y[1]=sin(x[0]-0.25);
  return 0;
}

int main(void) {
  ubvector x(2), y(2), dx(2);
  int i;
  int vp=0;
  size_t tmp;
  int r1, r2, r3;
  bool done;
  test_mgr t;
  std::vector<double> resid_test, resid_test2;

  cout.precision(10);
  cout.setf(ios::scientific);

  t.set_output_level(2);

  mroot_broyden<mm_funct11> gmb1;
  mm_funct11 fmf=gfn;

  // 1 - Normal execution using a member function
  x[0]=0.5;
  x[1]=0.5;
  gmb1.msolve(2,x,fmf);
  t.test_rel(x[0],0.25,1.0e-6,"1a");
  t.test_rel(x[1],0.2,1.0e-6,"1b");

  // 2 - Using the set(), iterate() interface
  
  x[0]=0.5;
  x[1]=0.5;
  dx[0]=0.1;
  dx[1]=0.0;
  gmb1.allocate(2);
  gmb1.set(fmf,2,x,y,dx);
  done=false;
  do {
    gmb1.iterate();
    double resid=fabs(y[0])+fabs(y[1]);
    cout << resid << endl;
    resid_test.push_back(resid);
    if (resid<gmb1.tol_rel) done=true;
  } while (done==false);
  t.test_rel(x[0],0.25,1.0e-6,"2a");
  t.test_rel(x[1],0.2,1.0e-6,"2b");

  // 7 - GSL version
  {
    gsl_multiroot_fsolver *s;
    
    int status;
    size_t iter=0;
    
    const size_t n=2;
    gsl_multiroot_function f={&gsl_fn, n, 0};
    double x_init[2]={0.5,0.5};
    gsl_vector *gx=gsl_vector_alloc (n);
    gsl_vector_set(gx,0,x_init[0]);
    gsl_vector_set(gx,1,x_init[1]);
    s=gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_broyden, 2);
    gsl_multiroot_fsolver_set(s, &f, gx);

    do {

      iter++;
      status=gsl_multiroot_fsolver_iterate(s);
      if (status) break;
      status=gsl_multiroot_test_residual (s->f, 1e-8);
      resid_test2.push_back(fabs(gsl_vector_get(s->f,0))+
			    fabs(gsl_vector_get(s->f,1)));
      cout << "resid: " << fabs(gsl_vector_get(s->f,0))+
	fabs(gsl_vector_get(s->f,1)) << endl;
      
    } while (status == GSL_CONTINUE && iter < 1000);
    
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(gx);
  }

  t.test_rel_vec(resid_test.size(),resid_test,resid_test2,1.0e-10,
		 "GSL vs. O2scl");

#ifdef NEVER_DEFINED

  // 1 - Normal execution using a member function
  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const Eigen::VectorXd &,Eigen::VectorXd &)>
     (&cl::mfn),&acl,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  mroot_broyden<> cr1;
  
  x[0]=0.5;
  x[1]=0.5;
  cr1.msolve(2,x,fmf);
  t.test_rel(x[0],0.25,1.0e-6,"1a");
  t.test_rel(x[1],0.2,1.0e-6,"1b");

  // 3 - Having specified the Jacobian
  jac_funct_mfptr<cl,ubvector,ubmatrix> fmfd(&acl,&cl::mfnd);
    
  x[0]=0.5;
  x[1]=0.5;
  cr1.msolve_de(2,x,fmf,fmfd);
  t.test_rel(x[0],0.25,1.0e-6,"3a");
  t.test_rel(x[1],0.2,1.0e-6,"3b");

  // 4 - Using the set_de(), iterate() interface
  mroot_broyden<> cr1c;

  x[0]=0.5;
  x[1]=0.5;
  cr1c.allocate(2);
  cr1c.set_de(2,x,fmf,fmfd);
  done=false;
  do {
    r3=cr1c.iterate();
    double resid=fabs(cr1c.f[0])+fabs(cr1c.f[1]);
    if (resid<cr1c.tol_rel || r3>0) done=true;
    resid_test.push_back(resid);
    cout << "resid: " << resid << endl;
  } while (done==false);
  t.test_rel(cr1c.x[0],0.25,1.0e-6,"4a");
  t.test_rel(cr1c.x[1],0.2,1.0e-6,"4b");
  cr1c.free();

  // 5 - Using arrays instead of ubvectors
  std::function<int(size_t,const arr_t &,arr_t &)> fmf2=std::bind
    (std::mem_fn<int(size_t,const arr_t &,arr_t &)>
     (&cl::mfna),&acl,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  //mm_funct_mfptr<cl,arr_t> fmf2(&acl,&cl::mfna);
  mroot_broyden<std::function<int(size_t,const arr_t &,arr_t &)>,arr_t,
    arr_t,array_alloc<arr_t> > cr2;

// Emacs had trouble with tabifying the following lines
// so I add some braces
#ifdef O2SCL_NEVER_DEFINED
}
{
#endif
  
  xa[0]=0.5;
  xa[1]=0.5;
  cr2.msolve(2,xa,fmf2);
  t.test_rel(xa[0],0.25,1.0e-6,"5a");
  t.test_rel(xa[1],0.2,1.0e-6,"5b");
  
  // using arrays with jacobian
  jac_funct_mfptr<cl,arr_t,mat_t> fmfd2(&acl,&cl::mfnda);
  mroot_broyden<std::function<int(size_t,const arr_t &,arr_t &)>,arr_t,
    arr_t,array_alloc<arr_t>,mat_t,mat_t,
    array_2d_alloc<mat_t>,jac_funct<arr_t,mat_t> > crx;

// Emacs had trouble with tabifying the following lines
// so I add some braces
#ifdef O2SCL_NEVER_DEFINED
}
{
#endif
  
  xa[0]=0.5;
  xa[1]=0.5;
  crx.msolve_de(2,xa,fmf2,fmfd2);
  t.test_rel(xa[0],0.25,1.0e-6,"8a");
  t.test_rel(xa[1],0.2,1.0e-6,"8b");
  
  // 6 - Using a global function pointer directly
  typedef int (*gfnt)(size_t, const ubvector &, ubvector &);
  mroot_broyden<gfnt,ubvector> cr4;
  gfnt gfnv=&gfn;

  x[0]=0.5;
  x[1]=0.5;
  cr4.msolve(2,x,gfnv);
  t.test_rel(x[0],0.25,1.0e-6,"6a");
  t.test_rel(x[1],0.2,1.0e-6,"6b");
  
  // 7 - GSL version
  {
    gsl_multiroot_fsolver *s;
    
    int status;
    size_t iter=0;
    
    const size_t n=2;
    gsl_multiroot_function f={&gsl_fn, n, 0};
    
    double x_init[2]={0.5,0.5};
    gsl_vector *gx=gsl_vector_alloc (n);
    
    gsl_vector_set(gx,0,x_init[0]);
    gsl_vector_set(gx,1,x_init[1]);
    
    s=gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_broyden, 2);
    
    gsl_multiroot_fsolver_set(s, &f, gx);

    do {

      iter++;
      status=gsl_multiroot_fsolver_iterate(s);
      
      if (status) break;
      
      status=gsl_multiroot_test_residual (s->f, 1e-8);
      resid_test2.push_back(fabs(gsl_vector_get(s->f,0))+
			    fabs(gsl_vector_get(s->f,1)));
      cout << "resid: " << fabs(gsl_vector_get(s->f,0))+
	fabs(gsl_vector_get(s->f,1)) << endl;
      
    } while (status == GSL_CONTINUE && iter < 1000);
    
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(gx);
  }
  
  t.test_rel_vec(resid_test.size(),resid_test,resid_test2,1.0e-2,
		 "GSL vs. O2scl");

#endif
  
  t.report();
  return 0;
}
