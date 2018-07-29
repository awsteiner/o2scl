/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <o2scl/funct.h>
#include <o2scl/inte_gauss_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double tx, double &a);

double testfun(double tx, double &a) {
  return (a*sin(tx)/(tx+0.01));
}

double testfun2(double tx) {
  return 4.0*std::sqrt(1.0-tx*tx);
}

long double testfun2_ld(long double tx) {
  return 4.0*std::sqrt(1.0-tx*tx);
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  {
    inte_gauss_cern<funct> cg;
    double a=3.0, calc, exact, diff;
  
    funct tf=std::bind(testfun,std::placeholders::_1,a);

    cout.setf(ios::scientific);
    cout.precision(10);
  
    calc=cg.integ(tf,0.0,1.0);
    exact=a*(0.900729064796877177);
    t.test_rel(calc,exact,1.0e-8,"inte_gauss_cern 1");
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << endl;

    funct tf2=testfun2;

    calc=cg.integ(tf2,0.0,1.0);
    exact=boost::math::constants::pi<double>();
    t.test_rel(calc,exact,1.0e-8,"inte_gauss_cern 2");
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << endl;
  }

  {
    inte_gauss_cern<funct_ld,long double,
		    o2scl::inte_gauss_cern_x_long_double,
		    o2scl::inte_gauss_cern_w_long_double> cg_ld;
    cg_ld.tol_rel=1.0e-20;
    cg_ld.tol_abs=1.0e-20;
    long double a=3.0, calc, exact, diff;

    funct_ld tf2=testfun2_ld;
    
    calc=cg_ld.integ(tf2,0.0,1.0);
    exact=boost::math::constants::pi<long double>();
    t.test_rel(calc,exact,1.0e-16L,"inte_gauss_cern ld");
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << endl;
  }
  
  t.report();
  return 0;
}

