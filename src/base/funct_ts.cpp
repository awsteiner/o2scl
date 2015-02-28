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
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  funct11_string f("pi*r^2","r",1,"pi");
  ubvector par(1);
  par[0]=o2scl_const::pi;
  f.set_parms(par);
  for(double r=1.0;r<=2.0;r+=0.1) {
    cout << f(r) << endl;
  }

  t.report();
  return 0;
}

