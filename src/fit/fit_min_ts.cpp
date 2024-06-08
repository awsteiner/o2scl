/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/fit_min.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

double func2(size_t np, const ubvector &p, double x) {
  return p[0]+x*p[1];
}

int main(void) {
  cout.setf(ios::scientific);
  test_mgr t;
  t.set_output_level(1);

  //----------------------------------------------------------------
  // Test fit_funct_fptr
  //----------------------------------------------------------------

  fit_funct ff2(func2);
  ubvector x_init(2);
  ubvector xdat(4), y(4), sigma(4);
  fit_min<> mf;
  double chi2;
  ubmatrix mycovar(2,2);

  x_init[0]=0.0;
  x_init[1]=1.0;
  xdat[0]=0.0;
  xdat[1]=1.0;
  xdat[2]=2.0;
  xdat[3]=3.0;
  y[0]=1.1;
  y[1]=2.9;
  y[2]=5.1;
  y[3]=8.0;
  sigma[0]=0.01;
  sigma[1]=0.1;
  sigma[2]=0.2;
  sigma[3]=1.0;

  chi_fit_funct<> cff(4,xdat,y,sigma,ff2);

  mf.verbose=2;
  mf.fit(2,x_init,mycovar,chi2,cff);

  cout << x_init[0] << " " << x_init[1] << endl;
  cout << sqrt(mycovar(0,0)) << " " << sqrt(mycovar(1,1)) << endl;
  cout << chi2 << endl;
  cout << endl;
  t.test_rel(x_init[0],1.0,0.1,"x0");
  t.test_rel(x_init[1],2.0,0.1,"x1");

  //----------------------------------------------------------------
  // Clean up
  //----------------------------------------------------------------

  t.report();

  return 0;
}

