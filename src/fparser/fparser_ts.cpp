/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2014, Andrew W. Steiner
  
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
#include <o2scl/fparser.h>
#include <o2scl/test_mgr.h>

#include <iostream>
#include <string>

using namespace std;
using namespace o2scl;

int main() {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  FunctionParser fparser;
  fparser.AddConstant("pi",acos(-1));

  int res=fparser.Parse("sin(pi*x)","x");
  t.test_gen(res==-1,"parse");
  
  cout.setf(ios::showpos);
  
  double arr[1];
  for(arr[0]=1.0;arr[0]<=4.001;arr[0]+=0.1) {
    t.test_rel(fparser.Eval(arr),sin(acos(-1.0)*arr[0]),1.0e-12,"sin");
  }
  
  t.report();
  return 0;
}
