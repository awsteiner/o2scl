/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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

#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_hybrids.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int gsl_fn(const gsl_vector *x, void *pa, gsl_vector *f) {
  gsl_vector_set(f,0,sin(gsl_vector_get(x,1)-0.2));
  gsl_vector_set(f,1,sin(gsl_vector_get(x,0)-0.25));
  return 0;
}

int gfn(size_t nv, const ubvector &x, 
	ubvector &y) {
  y[0]=sin(x[1]-0.2);
  y[1]=sin(x[0]-0.25);
  return 0;
}

class cl {

public:

  int mfn(size_t nv, const ubvector &x, ubvector &y) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }

  int mfnd(size_t nx, ubvector &x, size_t ny, ubvector &y, ubmatrix &j) {
    j(0,0)=0.0;
    j(0,1)=cos(x[1]-0.2);
    j(1,0)=cos(x[0]-0.25);
    j(1,1)=0.0;
    return 0;
  }

#ifdef O2SCL_EIGEN

  int mfn_Eigen(size_t nv, const Eigen::VectorXd &x, Eigen::VectorXd &y) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }

  int mfnd_Eigen(size_t nx, Eigen::VectorXd &x, size_t ny, 
		 Eigen::VectorXd &y, Eigen::MatrixXd &j) {
    j(0,0)=0.0;
    j(0,1)=cos(x[1]-0.2);
    j(1,0)=cos(x[0]-0.25);
    j(1,1)=0.0;
    return 0;
  }

#endif

#ifdef O2SCL_ARMA

  int mfn_arma(size_t nv, const arma::rowvec &x, arma::rowvec &y) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }

  int mfnd_arma(size_t nx, arma::rowvec &x, size_t ny, 
		 arma::rowvec &y, arma::mat &j) {
    j(0,0)=0.0;
    j(0,1)=cos(x[1]-0.2);
    j(1,0)=cos(x[0]-0.25);
    j(1,1)=0.0;
    return 0;
  }

#endif
  
};

int main(void) {
  cout.setf(ios::scientific);
  cout.precision(10);

  test_mgr t;
  t.set_output_level(2);

  cl acl;
  ubvector x(2);
  std::vector<double> resid_test, resid_test2;

  // 1 - Normal execution using a member function
  mm_funct11 fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&cl::mfn),&acl,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  mroot_hybrids<mm_funct11,ubvector,ubmatrix,jac_funct11> cr1;
  
  x[0]=0.5;
  x[1]=0.5;
  cr1.msolve(2,x,fmf);
  t.test_rel(x[0],0.25,1.0e-6,"normal a");
  t.test_rel(x[1],0.2,1.0e-6,"normal b");

  // 2 - Using the set(), iterate() interface
  mroot_hybrids<mm_funct11,ubvector,ubmatrix,jac_funct11> cr2;
  
  x[0]=0.5;
  x[1]=0.5;
  cr2.allocate(2);
  cr2.set(2,x,fmf);
  bool done=false;
  do {
    int r3=cr2.iterate();
    double resid=fabs(cr2.f[0])+fabs(cr2.f[1]);
    if (resid<cr2.tol_rel || r3>0) done=true;
  } while (done==false);
  t.test_rel(cr2.x[0],0.25,1.0e-6,"set/iterate a");
  t.test_rel(cr2.x[1],0.2,1.0e-6,"set/iterate b");

  // 3 - Having specified the Jacobian
  jac_funct11 fmfd=
    std::bind(std::mem_fn<int(size_t,ubvector &,size_t,
			      ubvector &,ubmatrix &)>(&cl::mfnd),
    &acl,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,std::placeholders::_4,
    std::placeholders::_5);
    
  x[0]=0.5;
  x[1]=0.5;
  cr1.msolve_de(2,x,fmf,fmfd);
  t.test_rel(x[0],0.25,1.0e-6,"jac a");
  t.test_rel(x[1],0.2,1.0e-6,"jac b");

  // 4 - Using the set_de(), iterate() interface
  mroot_hybrids<mm_funct11,ubvector,ubmatrix,jac_funct11> cr4;

  x[0]=0.5;
  x[1]=0.5;
  cr4.allocate(2);
  cr4.set_de(2,x,fmf,fmfd);
  done=false;
  do {
    int r3=cr4.iterate();
    double resid=fabs(cr4.f[0])+fabs(cr4.f[1]);
    if (resid<cr4.tol_rel || r3>0) done=true;
    resid_test.push_back(resid);
    cout << "resid: " << resid << endl;
  } while (done==false);
  t.test_rel(cr4.x[0],0.25,1.0e-6,"set_de/iterate a");
  t.test_rel(cr4.x[1],0.2,1.0e-6,"set_de/iterate b");

// Emacs had trouble with tabifying the following lines
// so I add some braces
#ifdef O2SCL_NEVER_DEFINED
}
{
#endif

  // 6 - Using a global function pointer directly
  typedef int (*gfnt)(size_t, const ubvector &, ubvector &);
  mroot_hybrids<gfnt,ubvector,ubmatrix,jac_funct11> cr6;
  gfnt gfnv=&gfn;

  x[0]=0.5;
  x[1]=0.5;
  cr6.msolve(2,x,gfnv);
  t.test_rel(x[0],0.25,1.0e-6,"global fptr a");
  t.test_rel(x[1],0.2,1.0e-6,"global fptr b");
  
  // 7 - GSL version
  {
    gsl_multiroot_fsolver *s;
    
    int status;
    size_t iter = 0;
    
    const size_t n = 2;
    gsl_multiroot_function f = {&gsl_fn, n, 0};
    
    double x_init[2] = {0.5,0.5};
    gsl_vector *gx = gsl_vector_alloc (n);
    
    gsl_vector_set(gx,0,x_init[0]);
    gsl_vector_set(gx,1,x_init[1]);
    
    s = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 2);
    
    gsl_multiroot_fsolver_set(s, &f, gx);

    do {

      iter++;
      status = gsl_multiroot_fsolver_iterate(s);
      
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

#ifdef O2SCL_EIGEN

  // 8 - Using Eigen
  typedef std::function<int(size_t,const Eigen::VectorXd &,
			    Eigen::VectorXd &) > mm_funct_Eigen;
  typedef std::function<int(size_t,Eigen::VectorXd &,
			    size_t,Eigen::VectorXd &,
			    Eigen::MatrixXd &) > jac_funct_Eigen;

  mm_funct_Eigen fmf_Eigen=std::bind
    (std::mem_fn<int(size_t,const Eigen::VectorXd &,Eigen::VectorXd &)>
     (&cl::mfn_Eigen),&acl,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  
  jac_funct_Eigen fmfd_Eigen=
    std::bind(std::mem_fn<int(size_t,Eigen::VectorXd &,size_t,
			      Eigen::VectorXd &,Eigen::MatrixXd &)>
    (&cl::mfnd_Eigen),
    &acl,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,std::placeholders::_4,
    std::placeholders::_5);

  mroot_hybrids<mm_funct_Eigen,Eigen::VectorXd,
		Eigen::MatrixXd,jac_funct_Eigen> cr8;
  
  Eigen::VectorXd xE(2);
  xE[0]=0.5;
  xE[1]=0.5;
  cr8.msolve_de(2,xE,fmf_Eigen,fmfd_Eigen);
  t.test_rel(xE[0],0.25,1.0e-6,"eigen a");
  t.test_rel(xE[1],0.2,1.0e-6,"eigen b");

#endif
  
#ifdef O2SCL_ARMA
  
  // 9 - Using Armadillo

  typedef std::function<int(size_t,const arma::rowvec &,
			    arma::rowvec &) > mm_funct_arma;
  typedef std::function<int(size_t,arma::rowvec &,
			    size_t,arma::rowvec &,
			    arma::mat &) > jac_funct_arma;

  mm_funct_arma fmf_arma=std::bind
    (std::mem_fn<int(size_t,const arma::rowvec &,arma::rowvec &)>
     (&cl::mfn_arma),&acl,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  
  jac_funct_arma fmfd_arma=
    std::bind(std::mem_fn<int(size_t,arma::rowvec &,size_t,
			      arma::rowvec &,arma::mat &)>
    (&cl::mfnd_arma),
    &acl,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,std::placeholders::_4,
    std::placeholders::_5);

  mroot_hybrids<mm_funct_arma,arma::rowvec,
		arma::mat,jac_funct_arma> cr9;

  arma::rowvec xA(2);
  xA[0]=0.5;
  xA[1]=0.5;
  cr9.msolve_de(2,xA,fmf_arma,fmfd_arma);
  t.test_rel(xA[0],0.25,1.0e-6,"arma a");
  t.test_rel(xA[1],0.2,1.0e-6,"arma b");

#endif
  
  // 1a - Member function with new C++11 extensions
  mm_funct11 f_new=
    std::bind(std::mem_fn<int(size_t,const ubvector &,ubvector &)>(&cl::mfn),
	      acl,std::placeholders::_1,std::placeholders::_2,
	      std::placeholders::_3);

#ifdef O2SCL_NEVER_DEFINED
} {
#endif
  mroot_hybrids<> cr1a;
  
  x[0]=0.5;
  x[1]=0.5;
  cr1a.msolve(2,x,f_new);
  t.test_rel(x[0],0.25,1.0e-6,"normal c++11 a");
  t.test_rel(x[1],0.2,1.0e-6,"normal c++11 b");

  // 2 - Using the set(), iterate() interface
  mroot_hybrids<> cr2a;
  
  x[0]=0.5;
  x[1]=0.5;
  cr2a.allocate(2);
  cr2a.set(2,x,f_new);
  done=false;
  do {
    int r3=cr2a.iterate();
    double resid=fabs(cr2a.f[0])+fabs(cr2a.f[1]);
    if (resid<cr2a.tol_rel || r3>0) done=true;
  } while (done==false);
  t.test_rel(cr2a.x[0],0.25,1.0e-6,"set/iterate c++11 a");
  t.test_rel(cr2a.x[1],0.2,1.0e-6,"set/iterate c++11 b");

  // 3 - Having specified the Jacobian
  jac_funct11 df_new=
    std::bind(std::mem_fn<int(size_t,ubvector &,size_t,
			      ubvector &,ubmatrix &)>(&cl::mfnd),
    acl,std::placeholders::_1,std::placeholders::_2,
    std::placeholders::_3,std::placeholders::_4,
    std::placeholders::_5);
    
#ifdef O2SCL_NEVER_DEFINED
} {
#endif

  x[0]=0.5;
  x[1]=0.5;
  cr1a.msolve_de(2,x,f_new,df_new);
  t.test_rel(x[0],0.25,1.0e-6,"jac c++11 a");
  t.test_rel(x[1],0.2,1.0e-6,"jac c++11 b");

  // 4 - Using the set_de(), iterate() interface
  mroot_hybrids<> cr4a;

  x[0]=0.5;
  x[1]=0.5;
  cr4a.allocate(2);
  cr4a.set_de(2,x,f_new,df_new);
  done=false;
  do {
    int r3=cr4a.iterate();
    double resid=fabs(cr4a.f[0])+fabs(cr4a.f[1]);
    if (resid<cr4a.tol_rel || r3>0) done=true;
    resid_test.push_back(resid);
    cout << "resid: " << resid << endl;
  } while (done==false);
  t.test_rel(cr4a.x[0],0.25,1.0e-6,"set_de/iterate c++11 a");
  t.test_rel(cr4a.x[1],0.2,1.0e-6,"set_de/iterate c++11 b");

  // 6 - Using a global function pointer directly
  typedef int (*gfnt)(size_t, const ubvector &, ubvector &);
  mroot_hybrids<gfnt> cr6a;
  gfnt gfnva=&gfn;

  x[0]=0.5;
  x[1]=0.5;
  cr6a.msolve(2,x,gfnva);
  t.test_rel(x[0],0.25,1.0e-6,"global fptr c++11 a");
  t.test_rel(x[1],0.2,1.0e-6,"global fptr c++11 b");
  
  // 7 - GSL version
  {
    gsl_multiroot_fsolver *s;
    
    int status;
    size_t iter = 0;
    
    const size_t n = 2;
    gsl_multiroot_function f = {&gsl_fn, n, 0};
    
    double x_init[2] = {0.5,0.5};
    gsl_vector *gx = gsl_vector_alloc (n);
    
    gsl_vector_set(gx,0,x_init[0]);
    gsl_vector_set(gx,1,x_init[1]);
    
    s = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, 2);
    
    gsl_multiroot_fsolver_set(s, &f, gx);

    do {

      iter++;
      status = gsl_multiroot_fsolver_iterate(s);
      
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
  
  t.report();
  return 0;
}
