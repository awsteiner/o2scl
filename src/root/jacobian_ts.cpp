/*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2013, Andrew W. Steiner
  
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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/jacobian.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int tfun(size_t nv, const ubvector &x, ubvector &y) {
  y[0]=x[0]*x[0]+sin(x[1]);
  y[1]=pow(x[1],3.0)+tan(x[0]);
  return 0;
}

int main(void) {

  jacobian_exact<mm_funct_fptr<> > ej;
  jacobian_gsl<mm_funct_fptr<> > sj;

  test_mgr t;
  t.set_output_level(2);

  mm_funct_fptr<> mff(tfun);
  
  ej.set_function(mff);
  sj.set_function(mff);

  ubvector x(2), y(2);
  ubmatrix jac(2,2), jex(2,2);

  x[0]=2.0;
  x[1]=3.0;

  tfun(2,x,y);

  cout.setf(ios::scientific);
  cout.precision(6);

  jex(0,0)=2.0*x[0];
  jex(0,1)=cos(x[1]);
  jex(1,0)=1.0/cos(x[0])/cos(x[0]);
  jex(1,1)=3.0*x[1]*x[1];

  cout << jex(0,0) << " " << jex(0,1) << endl;
  cout << jex(1,0) << " " << jex(1,1) << endl;
  cout << endl;

  sj(2,x,2,y,jac);
  cout << jac(0,0) << " " << jac(0,1) << endl;
  cout << jac(1,0) << " " << jac(1,1) << endl;
  cout << endl;

  t.test_rel(jac(0,0),jex(0,0),1.0e-3,"simple");
  t.test_rel(jac(0,1),jex(0,1),1.0e-3,"simple");
  t.test_rel(jac(1,0),jex(1,0),1.0e-3,"simple");
  t.test_rel(jac(1,1),jex(1,1),1.0e-3,"simple");
  cout << endl;
  
  ej(2,x,2,y,jac);
  cout << jac(0,0) << " " << jac(0,1) << endl;
  cout << jac(1,0) << " " << jac(1,1) << endl;
  cout << endl;
  
  t.test_rel(jac(0,0),jex(0,0),1.0e-10,"exact");
  t.test_rel(jac(0,1),jex(0,1),1.0e-10,"exact");
  t.test_rel(jac(1,0),jex(1,0),1.0e-10,"exact");
  t.test_rel(jac(1,1),jex(1,1),1.0e-10,"exact");
  cout << endl;

  cout << 2.0*x[0] << " " << cos(x[1]) << endl;
  cout << 1.0/cos(x[0])/cos(x[0]) << " "
       << 3.0*x[1]*x[1] << endl;
  cout << endl;

  t.report();
  return 0;
}

