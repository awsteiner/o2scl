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

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/funct.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double tx, double &pa);

double testfun(double tx, double &a) {
  return -cos(1.0/(tx+a))/(a+tx)/(a+tx);
}

int main(void) {
  double a, calc, exact, diff, ei;
  inte_adapt_cern<funct> cg;
  test_mgr t;
  t.set_output_level(2);
  
  a=0.01;
  funct tf=std::bind(testfun,std::placeholders::_1,a);

  cout.setf(ios::scientific);
  cout.precision(10);
  
  cg.integ_err(tf,0.0,1.0,calc,ei);
  exact=sin(1.0/(1.0+a))-sin(1.0/a);
  t.test_rel(calc,exact,1.0e-8,"inte_adapt_cern");
  diff=fabs(calc-exact);
  cout << calc << " " << exact << " " << diff << " " << ei << endl;

  // This is a nasty function and takes many subdivisions (68)
  cout << cg.get_nsubdivisions() << endl;
  size_t n=cg.get_nsubdivisions();
  typedef boost::numeric::ublas::vector<double> ubvector;
  ubvector xlo(n), xhi(n), val(n), err(n);
  cg.get_subdivisions(xlo,xhi,val,err);
  for(size_t i=0;i<n;i+=10) {
    cout << xlo[i] << " " << xhi[i] << " ";
    cout.setf(ios::showpos);
    cout << val[i] << " ";
    cout.unsetf(ios::showpos);
    cout << err[i] << endl;
  }

  a=0.01;
  cg.verbose=1;
  cout.precision(6);
  cg.integ_err(tf,0.0,1.0,calc,ei);
  cout.precision(10);
  exact=sin(1.0/(1.0+a))-sin(1.0/a);
  t.test_rel(calc,exact,1.0e-8,"inte_adapt_cern");
  diff=fabs(calc-exact);
  cout << calc << " " << exact << " " << diff << " " << ei << endl;

  t.report();
  return 0;
}

