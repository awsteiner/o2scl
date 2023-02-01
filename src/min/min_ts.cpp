/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/min.h>
#include <o2scl/funct.h>
#include <o2scl/min_cern.h>

using namespace std;
using namespace o2scl;

double minfun(double x) {
  return -exp(-(x-0.5)*(x-0.5));
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  // test constraint functions
  for(double x=50.0;x<=150.01;x+=10.0) {
    cout << x << " ";
    cout << constraint(x,100.0,20.0,1.0) << " "
	 << cont_constraint(x,100.0,20.0,1.0) << " ";
    cout << lower_bound(x,100.0,20.0,1.0) << " "
	 << cont_lower_bound(x,100.0,20.0,1.0) << endl;
  }
  
  min_cern<funct> mi;
  funct mf=minfun;
  int vp=0;
  double x1, x2, x3=0.0, f1=0.0, f2=0.0, f3=0.0;
  int ret;
  x1=0.0;
  x2=1.0;
  ret=mi.bracket(x1,x2,x3,f1,f2,f3,mf);
  cout << ret << endl;
  t.test_gen(ret==0,"Bracket 1");
  cout << x1 << " " << x3 << " " << x2 << endl;
  cout << f1 << " " << f3 << " " << f2 << endl;

  x1=0.0;
  x2=2.0;
  ret=mi.bracket(x1,x2,x3,f1,f2,f3,mf);
  t.test_gen(ret==0,"Bracket 2");
  cout << ret << endl;
  cout << x1 << " " << x3 << " " << x2 << endl;
  cout << f1 << " " << f3 << " " << f2 << endl;
  
  x1=10.0;
  x2=11.0;
  ret=mi.bracket(x1,x2,x3,f1,f2,f3,mf);
  t.test_gen(ret==0,"Bracket 3");
  cout << ret << endl;
  cout << err_hnd->get_str() << endl;
  cout << x1 << " " << x3 << " " << x2 << endl;
  cout << f1 << " " << f3 << " " << f2 << endl;

  t.report();
  return 0;
}
 
