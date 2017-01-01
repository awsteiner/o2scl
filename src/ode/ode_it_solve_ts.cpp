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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_ieee_utils.h>

#ifdef O2SCL_ARMA
#include <armadillo>
#endif
#ifdef O2SCL_EIGEN
#include <eigen3/Eigen/Dense>
#endif

#include <o2scl/ode_it_solve.h>
#include <o2scl/linear_solver.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_linalg;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
typedef boost::numeric::ublas::compressed_matrix<double> ubcomp_matrix;

// Simple two-dimensional system
class fc1 {
  
public:

  double left(size_t ieq, double x, ubmatrix_row &yleft) {
    return yleft[1]-1.0;
  }
  
  double right(size_t ieq, double x, ubmatrix_row &yright) {
    return yright[0]-2.0;
  }
  
  double derivs(size_t ieq, double x, ubmatrix_row &y) {
    if (ieq==0) return y[0]+y[1];
    return y[0];
  }

#ifdef O2SCL_ARMA

  double lefta(size_t ieq, double x, arma::rowvec &yleft) {
    return yleft[1]-1.0;
  }
  
  double righta(size_t ieq, double x, arma::rowvec &yright) {
    return yright[0]-2.0;
  }
  
  double derivsa(size_t ieq, double x, arma::rowvec &y) {
    if (ieq==0) return y[0]+y[1];
    return y[0];
  }

#endif

#ifdef O2SCL_EIGEN

  double lefte(size_t ieq, double x, Eigen::MatrixXd::RowXpr &yleft) {
    return yleft[1]-1.0;
  }
  
  double righte(size_t ieq, double x, Eigen::MatrixXd::RowXpr &yright) {
    return yright[0]-2.0;
  }
  
  double derivse(size_t ieq, double x, Eigen::MatrixXd::RowXpr &y) {
    if (ieq==0) return y[0]+y[1];
    return y[0];
  }

#endif
  
};

// Simple two-dimensional system with different boundary conditions
class fc2 {
  
public:
  
  double left(size_t ieq, double x, ubmatrix_row &yleft) {
    return yleft[0]-1.0;
  }
  
  double right(size_t ieq, double x, ubmatrix_row &yright) {
    return yright[1]-2.0;
  }
  
  double derivs(size_t ieq, double x, ubmatrix_row &y) {
    if (ieq==1) return y[0]+y[1];
    return y[1];
  }
  
};

// Three-dimensional system
class fc3 {
  
public:
  
  double left(size_t ieq, double x, ubmatrix_row &yleft) {
    if (ieq==0) return yleft[0]-1.0;
    return yleft[1]*yleft[1]+yleft[2]*yleft[2]-2.0;
  }
  
  double right(size_t ieq, double x, ubmatrix_row &yright) {
    return yright[1]-3.0;
  }
  
  double derivs(size_t ieq, double x, ubmatrix_row &y) {
    if (ieq==1) return y[0]+y[1];
    else if (ieq==2) return y[0]+y[2];
    return y[1];
  }
  
};

// Three-dimensional system with different boundary conditions
class fc4 {
  
public:
  
  double left(size_t ieq, double x, ubmatrix_row &yleft) {
    return yleft[0];
  }
  
  double right(size_t ieq, double x, ubmatrix_row &yright) {
    if (ieq==0) return yright[0]-1.0;
    return yright[1]*yright[1]-yright[2]*yright[2]-2.0;
  }
  
  double derivs(size_t ieq, double x, ubmatrix_row &y) {
    if (ieq==1) return y[0]+y[1];
    else if (ieq==2) return y[0]+y[2];
    return y[1];
  }
  
};

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  double s5=sqrt(5.0);

  // Solve system 1 on a small grid

  {
    ubvector x(5);
    ubmatrix y(5,2);
    for(int i=0;i<5;i++) {
      x[i]=((double)i)/4.0;
      y(i,0)=2.0*x[i];
      y(i,1)=1.0+x[i]/2;
    }
    
    ubmatrix A(10,10);
    ubvector rhs(10), dy(10);
    fc1 f1;

    ode_it_funct11 f_derivs=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::derivs),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_left=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::left),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_right=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::right),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       

    ode_it_solve<> oit;
    oit.solve(5,2,1,x,y,f_derivs,f_left,f_right,A,rhs,dy);
  }

  // Solve system 1 on a larger grid

  {
    ubvector x(11);
    ubmatrix y(11,2);
    for(int i=0;i<11;i++) {
      x[i]=((double)i)/10.0;
      y(i,0)=2.0*x[i];
      y(i,1)=1.0+x[i]/2;
    }
  
    ubmatrix A(22,22);
    ubvector rhs(22), dy(22);
    fc1 f1;

    ode_it_funct11 f_derivs=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::derivs),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_left=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::left),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_right=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::right),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       

    ode_it_solve<> oit;
    oit.solve(11,2,1,x,y,f_derivs,f_left,f_right,A,rhs,dy);

    // Compare to exact solution

    double sol1, sol2;
    for(int kk=0;kk<11;kk++) {
      double z=x[kk];
      sol1=2.0*exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-(-5.0+s5)*exp(s5/2.0)-s5*exp(0.5+s5)+
	 (5.0+s5)*exp(0.5*s5*(1.0+2.0*z))+
	 s5*exp(0.5+s5*z));
      sol2=exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-4.0*s5*exp(s5/2.0)+(5.0+s5)*exp(0.5+s5)+
	 4.0*s5*exp(0.5*s5*(1.0+2.0*z))-
	 (-5.0+s5)*exp(0.5+s5*z));
      if (kk==0) {
	t.test_rel(y(kk,0),sol1,7.0e-1,"sys1 o2scl 1");
	t.test_rel(y(kk,1),sol2,7.0e-1,"sys1 o2scl 2");
      } else {
	t.test_rel(y(kk,0),sol1,5.0e-2,"sys1 o2scl 3");
	t.test_rel(y(kk,1),sol2,5.0e-2,"sys1 o2scl 4");
      }
    }
  
  }

  // Try reversing the order of the equations,
  // solve system 2 on a large grid

  {
    ubvector x(11);
    ubmatrix y(11,2);
    for(int i=0;i<11;i++) {
      x[i]=((double)i)/10.0;
      y(i,0)=1.0+x[i]/2;
      y(i,1)=2.0*x[i];
    }
  
    ubmatrix A(22,22);
    ubvector rhs(22), dy(22);

    fc2 f2;

    ode_it_funct11 f_derivs=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc2::derivs),&f2,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_left=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc2::left),&f2,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_right=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc2::right),&f2,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       

    ode_it_solve<> oit;
    oit.solve(11,2,1,x,y,f_derivs,f_left,f_right,A,rhs,dy);

    // Compare to exact result
    for(int kk=0;kk<11;kk++) {
      double z=x[kk];
      double sol1, sol2;
      sol1=exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-4.0*s5*exp(s5/2.0)+(5.0+s5)*exp(0.5+s5)+
	 4.0*s5*exp(0.5*s5*(1.0+2.0*z))-
	 (-5.0+s5)*exp(0.5+s5*z));
      sol2=2.0*exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-(-5.0+s5)*exp(s5/2.0)-s5*exp(0.5+s5)+
	 (5.0+s5)*exp(0.5*s5*(1.0+2.0*z))+
	 s5*exp(0.5+s5*z));
      if (kk==0) {
	t.test_rel(y(kk,0),sol1,7.0e-1,"sys2 o2scl 1");
	t.test_rel(y(kk,1),sol2,7.0e-1,"sys2 o2scl 2");
      } else {
	t.test_rel(y(kk,0),sol1,5.0e-2,"sys2 o2scl 3");
	t.test_rel(y(kk,1),sol2,5.0e-2,"sys2 o2scl 4");
      }
    }

  }
  
  // Solve system 3

  {
    ubvector x4(5);
    ubmatrix y(5,3);
    for(int i=0;i<5;i++) {
      x4[i]=((double)i)/4.0;
      y(i,0)=1.0+x4[i]+1.0;
      y(i,1)=3.0*x4[i];
      y(i,2)=-0.1*x4[i]-1.4;
    }
  
    ubmatrix A4(15,15);
    ubvector rhs4(15), dy(15);
    fc3 f3;

    ode_it_funct11 ofm3d=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc3::derivs),&f3,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 ofm3l=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc3::left),&f3,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 ofm3r=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc3::right),&f3,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       

    ode_it_solve<> oit;
    oit.solve(5,3,2,x4,y,ofm3d,ofm3l,ofm3r,A4,rhs4,dy);

    t.test_rel(y(0,0),1.0,1.0e-3,"sys3 o2scl 1");
    t.test_rel(y(4,0),2.30666,4.0e-3,"sys3 o2scl 2");
    t.test_rel(y(0,1),0.259509,1.0e-1,"sys3 o2scl 3");
    t.test_rel(y(4,1),3.0,1.0e-3,"sys3 o2scl 4");
    t.test_rel(y(0,2),-1.3902,4.0e-3,"sys3 o2scl 5");
    t.test_rel(y(4,2),-1.48437,3.0e-2,"sys3 o2scl 6");
  }
  
  // Solve system 4
  
  {
    ubvector x(5);
    ubmatrix y(5,3);
    for(int i=0;i<5;i++) {
      x[i]=((double)i)/4.0;
      y(i,0)=x[i];
      y(i,1)=1.5*x[i]+0.5;
      y(i,2)=-x[i]-0.6;
    }
  
    ubmatrix A(15,15);
    ubvector rhs(15), dy(15);

    fc4 f4;

    ode_it_funct11 ofm4d=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc4::derivs),&f4,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 ofm4l=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc4::left),&f4,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 ofm4r=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc4::right),&f4,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       

    ode_it_solve<> oit;
    oit.solve(5,3,1,x,y,ofm4d,ofm4l,ofm4r,A,rhs,dy);
  
    t.test_rel(y(0,0),0.0,1.0e-3,"sys4 o2scl 1");
    t.test_rel(y(4,0),1.0,1.0e-3,"sys4 o2scl 2");
    t.test_rel(y(0,1),0.496445,4.0e-2,"sys4 o2scl 3");
    t.test_rel(y(4,1),1.88562,1.0e-2,"sys4 o2scl 4");
    t.test_rel(y(0,2),-0.656063,1.0e-3,"sys4 o2scl 5");
    t.test_rel(y(4,2),-1.24722,1.0e-2,"sys4 o2scl 6");
  }

  // System 1 with sparse matrix format
  {
    ubvector x(11);
    ubmatrix y(11,2);
    for(int i=0;i<11;i++) {
      x[i]=((double)i)/10.0;
      y(i,0)=2.0*x[i];
      y(i,1)=1.0+x[i]/2;
    }
  
    ubcomp_matrix A(22,22);
    ubvector rhs(22), dy(22);
    fc1 f1;

    ode_it_funct11 f_derivs=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::derivs),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_left=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::left),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_right=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::right),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       

    ode_it_solve<ode_it_funct11,ubvector,ubmatrix,
		 ubmatrix_row,ubvector,ubcomp_matrix> oit;
    oit.solve(11,2,1,x,y,f_derivs,f_left,f_right,A,rhs,dy);

    // Compare to exact solution
    double sol1, sol2;
    for(int kk=0;kk<11;kk++) {
      double z=x[kk];
      sol1=2.0*exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-(-5.0+s5)*exp(s5/2.0)-s5*exp(0.5+s5)+
	 (5.0+s5)*exp(0.5*s5*(1.0+2.0*z))+
	 s5*exp(0.5+s5*z));
      sol2=exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-4.0*s5*exp(s5/2.0)+(5.0+s5)*exp(0.5+s5)+
	 4.0*s5*exp(0.5*s5*(1.0+2.0*z))-
	 (-5.0+s5)*exp(0.5+s5*z));
      if (kk==0) {
	t.test_rel(y(kk,0),sol1,7.0e-1,"sys1 sparse 1");
	t.test_rel(y(kk,1),sol2,7.0e-1,"sys1 sparse 2");
      } else {
	t.test_rel(y(kk,0),sol1,5.0e-2,"sys1 sparse 3");
	t.test_rel(y(kk,1),sol2,5.0e-2,"sys1 sparse 4");
      }
    }
  
  }

#ifdef O2SCL_NEVER_DEFINED

#ifdef O2SCL_ARMA

  // Armadillo with dense matrices
  {
    fc1 f1;

    ode_it_funct11 f_derivs=std::bind
      (std::mem_fn<double(size_t,double,arma::rowvec &)>
       (&fc1::derivsa),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_left=std::bind
      (std::mem_fn<double(size_t,double,arma::rowvec &)>
       (&fc1::lefta),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_right=std::bind
      (std::mem_fn<double(size_t,double,arma::rowvec &)>
       (&fc1::righta),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       

    ode_it_solve<ode_it_funct11,arma::colvec,arma::mat,
		 arma::rowvec,arma::colvec,arma::mat> oit;
    
    arma::colvec x(11);
    arma::mat y(11,2);
    for(int i=0;i<11;i++) {
      x[i]=((double)i)/10.0;
      y(i,0)=2.0*x[i];
      y(i,1)=1.0+x[i]/2;
    }
    
    arma::mat A(22,22);
    arma::colvec rhsa(22), dy(22);
    linear_solver_arma<arma::colvec,arma::mat> sol;
    oit.set_solver(sol);
    oit.solve(11,2,1,x,y,f_derivs,f_left,f_right,A,rhsa,dy);

    // Compare to exact solution
    double sol1, sol2;
    for(int kk=0;kk<11;kk++) {
      double z=x[kk];
      sol1=2.0*exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-(-5.0+s5)*exp(s5/2.0)-s5*exp(0.5+s5)+
	 (5.0+s5)*exp(0.5*s5*(1.0+2.0*z))+
	 s5*exp(0.5+s5*z));
      sol2=exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-4.0*s5*exp(s5/2.0)+(5.0+s5)*exp(0.5+s5)+
	 4.0*s5*exp(0.5*s5*(1.0+2.0*z))-
	 (-5.0+s5)*exp(0.5+s5*z));
      if (kk==0) {
	t.test_rel(y(kk,0),sol1,7.0e-1,"arma dense 1");
	t.test_rel(y(kk,1),sol2,7.0e-1,"arma dense 2");
      } else {
	t.test_rel(y(kk,0),sol1,5.0e-2,"arma dense 3");
	t.test_rel(y(kk,1),sol2,5.0e-2,"arma dense 4");
      }
    }

  }
  
  // Armadillo with sparse matrices and the o2scl solver
  {
    fc1 f1;

    ode_it_funct f_derivs=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::derivsa),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_left=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::lefta),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct11 f_right=std::bind
      (std::mem_fn<double(size_t,double,ubmatrix_row &)>
       (&fc1::righta),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       

    ode_it_solve<ode_it_funct11,arma::colvec,arma::mat,
		 arma::rowvec,arma::colvec,arma::sp_mat> oit;
    
    arma::colvec x(11);
    arma::mat y(11,2);
    for(int i=0;i<11;i++) {
      x[i]=((double)i)/10.0;
      y(i,0)=2.0*x[i];
      y(i,1)=1.0+x[i]/2;
    }
    
    arma::sp_mat A(22,22);
    arma::colvec rhs(22), dy(22);
    oit.solve(11,2,1,x,y,f_derivs,f_left,f_right,A,rhs,dy);

    // Compare to exact solution
    double sol1, sol2;
    for(int kk=0;kk<11;kk++) {
      double z=x[kk];
      sol1=2.0*exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-(-5.0+s5)*exp(s5/2.0)-s5*exp(0.5+s5)+
	 (5.0+s5)*exp(0.5*s5*(1.0+2.0*z))+
	 s5*exp(0.5+s5*z));
      sol2=exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-4.0*s5*exp(s5/2.0)+(5.0+s5)*exp(0.5+s5)+
	 4.0*s5*exp(0.5*s5*(1.0+2.0*z))-
	 (-5.0+s5)*exp(0.5+s5*z));
      if (kk==0) {
	t.test_rel(y(kk,0),sol1,7.0e-1,"arma sparse 1");
	t.test_rel(y(kk,1),sol2,7.0e-1,"arma sparse 2");
      } else {
	t.test_rel(y(kk,0),sol1,5.0e-2,"arma sparse 3");
	t.test_rel(y(kk,1),sol2,5.0e-2,"arma sparse 4");
      }
    }

  }
  
#endif

#endif

#ifdef O2SCL_EIGEN

  // Eigen with dense matrices and QR solver with column pivoting
  {
    fc1 f1;

    typedef std::function
      <double(size_t,double,Eigen::MatrixXd::RowXpr &)> ode_it_funct_eigen;
    
    ode_it_funct_eigen f_derivs=std::bind
      (std::mem_fn<double(size_t,double,Eigen::MatrixXd::RowXpr &)>
       (&fc1::derivse),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct_eigen f_left=std::bind
      (std::mem_fn<double(size_t,double,Eigen::MatrixXd::RowXpr &)>
       (&fc1::lefte),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    ode_it_funct_eigen f_right=std::bind
      (std::mem_fn<double(size_t,double,Eigen::MatrixXd::RowXpr &)>
       (&fc1::righte),&f1,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);       
    
    ode_it_solve<ode_it_funct_eigen,Eigen::VectorXd,Eigen::MatrixXd,
		 Eigen::MatrixXd::RowXpr,Eigen::VectorXd,Eigen::MatrixXd> oit;
    
    Eigen::VectorXd x(11);
    Eigen::MatrixXd y(11,2);
    for(int i=0;i<11;i++) {
      x[i]=((double)i)/10.0;
      y(i,0)=2.0*x[i];
      y(i,1)=1.0+x[i]/2;
    }
    
    Eigen::MatrixXd A(22,22);
    Eigen::VectorXd rhs(22), dy(22);
    linear_solver_eigen_colQR<Eigen::VectorXd,Eigen::MatrixXd> sol;
    oit.set_solver(sol);
    oit.solve(11,2,1,x,y,f_derivs,f_left,f_right,A,rhs,dy);

    // Compare to exact solution
    double sol1, sol2;
    for(int kk=0;kk<11;kk++) {
      double z=x[kk];
      sol1=2.0*exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-(-5.0+s5)*exp(s5/2.0)-s5*exp(0.5+s5)+
	 (5.0+s5)*exp(0.5*s5*(1.0+2.0*z))+
	 s5*exp(0.5+s5*z));
      sol2=exp(0.5*(z-s5*z))/
	(5.0*sqrt(exp(1.0))+(5.0+s5)*exp(0.5+s5)-
	 sqrt(5.0*exp(1.0)))*
	(-4.0*s5*exp(s5/2.0)+(5.0+s5)*exp(0.5+s5)+
	 4.0*s5*exp(0.5*s5*(1.0+2.0*z))-
	 (-5.0+s5)*exp(0.5+s5*z));
      if (kk==0) {
	t.test_rel(y(kk,0),sol1,7.0e-1,"eigen dense 1");
	t.test_rel(y(kk,1),sol2,7.0e-1,"eigen dense 2");
      } else {
	t.test_rel(y(kk,0),sol1,5.0e-2,"eigen dense 3");
	t.test_rel(y(kk,1),sol2,5.0e-2,"eigen dense 4");
      }
    }

  }
  
#endif

  t.report();

  return 0;
}
