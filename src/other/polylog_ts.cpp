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

#include <o2scl/test_mgr.h>
#include <o2scl/polylog.h>

using namespace std;
using namespace o2scl;

int main(void) {
  polylog po;
  test_mgr t;
  t.set_output_level(2);
  t.test_rel(po.li1(-0.99)-po.li1(-1.0),po.li1(-1.0)-po.li1(-1.01),
	     1.e-2,"Li1");
  t.test_rel(po.li2(-0.99)-po.li2(-1.0),po.li2(-1.0)-po.li2(-1.01),
	     1.e-2,"Li2");
  t.test_rel(po.li3(-0.99)-po.li3(-1.0),po.li3(-1.0)-po.li3(-1.01),
	     1.e-2,"Li3");
  t.test_rel(po.li4(-0.99)-po.li4(-1.0),po.li4(-1.0)-po.li4(-1.01),
	     1.e-2,"Li4");
  t.test_rel(po.li5(-0.99)-po.li5(-1.0),po.li5(-1.0)-po.li5(-1.01),
	     1.e-2,"Li5");
  t.test_rel(po.li6(-0.99)-po.li6(-1.0),po.li6(-1.0)-po.li6(-1.01),
	     1.e-2,"Li6");
  t.report();
  return 0;
}
