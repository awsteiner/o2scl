/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/umatrix_tlate.h>

using namespace std;
using namespace o2scl;

#include <o2scl/test_mgr.h>

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  umatrix a(3,3);
  a[0][0]=0.0;
  a[0][1]=1.0;
  a[0][2]=4.0;
  a[1][0]=1.0;
  a[1][1]=5.0;
  a[1][2]=9.0;
  a[2][0]=2.0;
  a[2][1]=6.0;
  a[2][2]=5.0;
  a.set(0,0,3.0);

  umatrix_row v1(a,1);
  t.test_rel(v1[0],1.0,1.0e-12,"row 1");
  t.test_rel(v1[1],5.0,1.0e-12,"row 2");
  t.test_rel(v1[2],9.0,1.0e-12,"row 3");
  
  t.report();
  return 0;
}
