/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2012, Andrew W. Steiner
  
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

/* Example: ex_mroot.cpp
   -------------------------------------------------------------------
   Several ways to use an O2scl solver to solve a simple function
*/

#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/mroot_cern.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int gfn(size_t nv, const ubvector &x, ubvector &y) {
  y[0]=sin(x[1]-0.2);
  y[1]=sin(x[0]-0.25);
  return 0;
}

class cl {

public:

  // Store the number of function and derivative evaluations
  int nf, nd;

  int mfn(size_t nv, const ubvector &x, ubvector &y) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    nf++;
    return 0;
  }

  int operator()(size_t nv, const ubvector &x, ubvector &y) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    nf++;
    return 0;
  }
  
  int mfnd(size_t nx, ubvector &x, size_t ny, 
	   ubvector &y, ubmatrix &j) {
    j(0,0)=0.0;
    j(0,1)=cos(x[1]-0.2);
    j(1,0)=cos(x[0]-0.25);
    j(1,1)=0.0;
    nd++;
    return 0;
  }

  template<class vec_t>
  int mfn_tlate(size_t nv, const vec_t &x, vec_t &y) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    nf++;
    return 0;
  }

  template<class vec_t, class mat_t>
  int mfnd_tlate(size_t nx, vec_t &x, size_t ny, vec_t &y, mat_t &j) {
    j(0,0)=0.0;
    j(0,1)=cos(x[1]-0.2);
    j(1,0)=cos(x[0]-0.25);
    j(1,1)=0.0;
    nd++;
    return 0;
  }

};

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  cl acl;
  ubvector x(2);

  /*
    Using a member function with ublas vectors
  */
  mm_funct_mfptr<cl> f1(&acl,&cl::mfn);
  mroot_hybrids<mm_funct<>,ubvector,ubmatrix,jac_funct<> > cr1;
    
  x[0]=0.5;
  x[1]=0.5;
  acl.nf=0;
  int ret1=cr1.msolve(2,x,f1);
  cout << "GSL solver (numerical Jacobian): " << endl;
  cout << "Return value: " << ret1 << endl;
  cout << "Number of iterations: " << cr1.last_ntrial << endl;
  cout << "Number of function evaluations: " << acl.nf << endl;
  cout << endl;
  t.test_rel(x[0],0.25,1.0e-6,"1a");
  t.test_rel(x[1],0.2,1.0e-6,"1b");

  /*
    Using the CERNLIB solver
  */
  mroot_cern<mm_funct<> > cr2;
    
  x[0]=0.5;
  x[1]=0.5;
  acl.nf=0;
  int ret2=cr2.msolve(2,x,f1);
  cout << "CERNLIB solver (numerical Jacobian): " << endl;
  cout << "Return value: " << ret2 << endl;
  cout << "INFO parameter: " << cr2.get_info() << endl;
  cout << "Number of function evaluations: " << acl.nf << endl;
  cout << endl;
  t.test_rel(x[0],0.25,1.0e-6,"2a");
  t.test_rel(x[1],0.2,1.0e-6,"2b");

  /*
    Using a member function with \ref ovector objects, but
    using the GSL-like interface with set() and iterate().
  */
  mroot_hybrids<mm_funct<> > cr3;

  x[0]=0.5;
  x[1]=0.5;
  cr3.allocate(2);
  cr3.set(2,x,f1);
  bool done=false;
  do {
    double r3=cr3.iterate();
    double resid=fabs(cr3.f[0])+fabs(cr3.f[1]);
    if (resid<cr3.tol_rel || r3>0) done=true;
  } while (done==false);
  t.test_rel(cr3.x[0],0.25,1.0e-6,"3a");
  t.test_rel(cr3.x[1],0.2,1.0e-6,"3b");

  /*
    Now instead of using the automatic Jacobian, using 
    a user-specified Jacobian.
  */
  jac_funct_mfptr<cl,ubvector,ubmatrix> j4(&acl,&cl::mfnd);
    
  x[0]=0.5;
  x[1]=0.5;
  acl.nf=0;
  acl.nd=0;
  int ret4=cr1.msolve_de(2,x,f1,j4);
  cout << "GSL solver (analytic Jacobian): " << endl;
  cout << "Return value: " << ret4 << endl;
  cout << "Number of iterations: " << cr1.last_ntrial << endl;
  cout << "Number of function evaluations: " << acl.nf << endl;
  cout << "Number of Jacobian evaluations: " << acl.nd << endl;
  cout << endl;
  t.test_rel(x[0],0.25,1.0e-6,"4a");
  t.test_rel(x[1],0.2,1.0e-6,"4b");

  /*
    Using a user-specified Jacobian and the GSL-like interface
  */
  mroot_hybrids<mm_funct<>,ubvector,ubmatrix,jac_funct<> > cr5;

  x[0]=0.5;
  x[1]=0.5;
  cr5.allocate(2);
  cr5.set_de(2,x,f1,j4);
  done=false;
  do {
    double r3=cr5.iterate();
    double resid=fabs(cr5.f[0])+fabs(cr5.f[1]);
    if (resid<cr5.tol_rel || r3>0) done=true;
  } while (done==false);
  t.test_rel(cr5.x[0],0.25,1.0e-6,"5a");
  t.test_rel(cr5.x[1],0.2,1.0e-6,"5b");

  /*
    Using a class with an operator(). Note that there can be only one
    operator() function in each class.
  */
  mroot_hybrids<cl> cr9;

  x[0]=0.5;
  x[1]=0.5;
  cr9.msolve(2,x,acl);
  t.test_rel(x[0],0.25,1.0e-6,"9a");
  t.test_rel(x[1],0.2,1.0e-6,"9b");

  /*
    Using a function pointer to a global function.
  */
  typedef int (*gfnt)(size_t, const ubvector &, ubvector &);
  mroot_hybrids<gfnt> cr10;
  gfnt f10=&gfn;

  x[0]=0.5;
  x[1]=0.5;
  cr10.msolve(2,x,f10);
  t.test_rel(x[0],0.25,1.0e-6,"10a");
  t.test_rel(x[1],0.2,1.0e-6,"10b");

  /* 
     Using different vector types: std::vector<double>, 
     ublas::unbounded_array<double>, and
     ublas::bounded_array<double>
  */
  std::vector<double> svx(2);
  svx[0]=0.5;
  svx[1]=0.5;
  mm_funct_mfptr<cl,std::vector<double> > 
    f11(&acl,&cl::mfn_tlate<std::vector<double> >);
  mroot_hybrids<mm_funct<std::vector<double> >,std::vector<double>,
		    ubmatrix,jac_funct<std::vector<double>,
				       ubmatrix> > cr11;
  cr11.msolve(2,svx,f11);
  t.test_rel(x[0],0.25,1.0e-6,"11a");
  t.test_rel(x[1],0.2,1.0e-6,"11b");

  typedef boost::numeric::ublas::unbounded_array<double> uarr;
  uarr uad(2);
  uad[0]=0.5;
  uad[1]=0.5;
  mm_funct_mfptr<cl,uarr> f12(&acl,&cl::mfn_tlate<uarr>);
  mroot_hybrids<mm_funct<uarr >,uarr,ubmatrix,
		    jac_funct<uarr,ubmatrix> > cr12;
  cr12.msolve(2,uad,f12);
  t.test_rel(x[0],0.25,1.0e-6,"12a");
  t.test_rel(x[1],0.2,1.0e-6,"12b");

  typedef boost::numeric::ublas::bounded_array<double,2> barr;
  barr bad(2);
  bad[0]=0.5;
  bad[1]=0.5;
  mm_funct_mfptr<cl,barr> f13(&acl,&cl::mfn_tlate<barr>);
  mroot_hybrids<mm_funct<barr >,barr,ubmatrix,
		    jac_funct<barr,ubmatrix> > cr13;
  cr13.msolve(2,bad,f13);
  t.test_rel(x[0],0.25,1.0e-6,"13a");
  t.test_rel(x[1],0.2,1.0e-6,"13b");

#ifdef O2SCL_EIGEN

  // Using Eigen with an analytic Jacobian
  mm_funct_mfptr<cl,Eigen::VectorXd> f14
    (&acl,&cl::mfn_tlate<Eigen::VectorXd>);
  jac_funct_mfptr<cl,Eigen::VectorXd,Eigen::MatrixXd> 
    fd14(&acl,&cl::mfnd_tlate<Eigen::VectorXd>);
  mroot_hybrids<mm_funct<Eigen::VectorXd>,Eigen::VectorXd,
		    Eigen::MatrixXd,jac_funct<Eigen::VectorXd,
					      Eigen::MatrixXd> > cr14;

  Eigen::VectorXd xE(2);
  xE[0]=0.5;
  xE[1]=0.5;
  cr14.msolve_de(2,xE,f14,fd14);
  t.test_rel(xE[0],0.25,1.0e-6,"eigen a");
  t.test_rel(xE[1],0.2,1.0e-6,"eigen b");

#endif

  t.report();
  return 0;
}
// End of example
