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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <complex>
#include <o2scl/umatrix_cx_tlate.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  umatrix_cx oc(2,2);
  oc.set(0,0,3.0,-1.0);
  oc.set(0,1,-4.0,1.0);
  oc.set(1,0,5.0,-9.0);
  oc.set(1,1,-2.0,6.0);
  gsl_complex g;
  g=oc(1,1);
  t.test_rel(GSL_REAL(g),-2.0,1.0e-6,"r1");
  t.test_rel(GSL_IMAG(g),6.0,1.0e-6,"c1");

  t.report();
  return 0;
}
