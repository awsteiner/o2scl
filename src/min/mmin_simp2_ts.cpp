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
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/multi_funct.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/test_mgr.h>
#ifdef O2SCL_EIGEN
#include <eigen3/Eigen/Dense>
#endif

using namespace std;
using namespace o2scl;
#ifdef O2SCL_EIGEN
using namespace Eigen;
#endif

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

double minfun(size_t n, const ubvector &x) {
  double ret;
  ret=x[0]*x[0]+(x[1]-2.0)*(x[1]-2.0)+3.0;
  return ret;
}

double minfun2(size_t n, const ubvector &x) {
  double ret;
  ret=x[0]*x[0]+(x[1]-2.0e-9)*(x[1]-2.0e-9)+3.0;
  return ret;
}

#ifdef O2SCL_EIGEN
double minfun_Eigen(size_t n, const VectorXd &x) {
  double ret;
  ret=x[0]*x[0]+(x[1]-2.0)*(x[1]-2.0)+3.0;
  return ret;
}
#endif

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  double min=0.0, min2;
  ubvector x(2), x2(2);
  mmin_simp2<multi_funct> g;
  
  multi_funct mf=minfun;
  multi_funct mf2=minfun2;
  
  // Standard function

  x[0]=1.1;
  x[1]=0.9;
  g.mmin(2,x,min,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(x[0],0.0,1.0e-4,"mmin_simp2 1");
  t.test_rel(x[1],2.0,2.0e-4,"mmin_simp2 2");

  // Function demonstrating bad scaling

  x[0]=1.1e-9;
  x[1]=9.0e-10;
  g.mmin(2,x,min,mf2);
  cout << min << " " << x[0] << " " << x[1] << endl;
  //t.test_abs(x[0],0.0,1.0e-4,"mmin_simp2 1");
  //t.test_rel(x[1],2.0e-9,2.0e-4,"mmin_simp2 2");

  // GSL-like interface

  x[0]=1.0;
  x[1]=1.0;
  ubvector step_size(2);
  step_size[0]=0.1;
  step_size[1]=0.1;
  g.set(mf,2,x,step_size);
  for(size_t it=0;it<100 && g.size>g.tol_abs;it++) {
    g.iterate();
  }
  x[0]=g.x[0];
  x[1]=g.x[1];
  min=g.fval;
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(x[0],0.0,1.0e-4,"mmin_simp2 set/iterate 1");
  t.test_rel(x[1],2.0,2.0e-4,"mmin_simp2 set/iterate 2");

  // Specify full simplex

  ubmatrix simp(3,2);
  simp(0,0)=1.0;
  simp(0,1)=1.0;
  simp(1,0)=1.1;
  simp(1,1)=1.1;
  simp(2,0)=2.0;
  simp(2,1)=1.0;
  g.mmin_simplex(2,simp,min,mf);
  cout << min << " " << x[0] << " " << x[1] << endl;
  t.test_abs(simp(0,0),0.0,1.0e-4,"mmin_simp2 full simplex 1");
  t.test_rel(simp(0,1),2.0,2.0e-4,"mmin_simp2 full simplex 2");

#ifdef O2SCL_EIGEN

  // Test with Eigen

  typedef std::function<double(size_t,const Eigen::VectorXd &) > 
    multi_funct_Eigen;

  mmin_simp2<multi_funct_Eigen,VectorXd> g_Eigen;
  multi_funct_Eigen mf_Eigen=minfun_Eigen;

  VectorXd x_Eigen(2);
  x_Eigen[0]=1.1;
  x_Eigen[1]=0.9;
  g_Eigen.mmin(2,x_Eigen,min,mf_Eigen);
  cout << min << " " << x_Eigen[0] << " " << x_Eigen[1] << endl;
  t.test_abs(x_Eigen[0],0.0,1.0e-4,"mmin_simp2 Eigen 1");
  t.test_rel(x_Eigen[1],2.0,2.0e-4,"mmin_simp2 Eigen 2");

#endif

  t.report();
  return 0;
}
 
