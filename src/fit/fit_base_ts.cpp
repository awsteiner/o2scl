/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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

#include <o2scl/fit_base.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double func(size_t np, const ubvector &p, double x);

double func(size_t np, const ubvector &p, double x) {
  return p[0]+x*p[1];
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  
  fit_funct11 f1=func;

  ubvector par(2);
  par[0]=0.5;
  par[1]=1.5;
  double x=3.2, y;

  y=f1(2,par,x);
  t.test_rel(y,5.3,1.0e-6,"fptr");

  t.report();
  return 0;
}

