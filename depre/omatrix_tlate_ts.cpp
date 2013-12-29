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
#include <o2scl/omatrix_tlate.h>
#include <o2scl/columnify.h>
#include <o2scl/test_mgr.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  omatrix a(3,3);
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

  gsl_matrix *gm=(gsl_matrix *)(&a);
  t.test_rel(gsl_matrix_get(gm,0,0),3.0,1.0e-12,"get 1");
  t.test_rel(gsl_matrix_get(gm,1,1),5.0,1.0e-12,"get 2");
  t.test_rel(gsl_matrix_get(gm,2,2),5.0,1.0e-12,"get 3");

  gsl_matrix *gm2=gsl_matrix_alloc(3,3);
  gsl_matrix_set(gm2,1,2,9.0);
  omatrix *amp=(omatrix *)gm2;
  t.test_rel((*amp)[1][2],9.0,1.0e-12,"get 4");
  
  omatrix_row v1(a,1);
  t.test_rel(v1[0],1.0,1.0e-12,"row 1");
  t.test_rel(v1[1],5.0,1.0e-12,"row 2");
  t.test_rel(v1[2],9.0,1.0e-12,"row 3");

  omatrix_col v2(a,1);
  t.test_rel(v2[0],1.0,1.0e-12,"col 1");
  t.test_rel(v2[1],5.0,1.0e-12,"col 2");
  t.test_rel(v2[2],6.0,1.0e-12,"col 3");

  umatrix ua(3,3);
  ua[0][0]=0.0;
  ua[0][1]=1.0;
  ua[0][2]=4.0;
  ua[1][0]=1.0;
  ua[1][1]=5.0;
  ua[1][2]=9.0;
  ua[2][0]=2.0;
  ua[2][1]=6.0;
  ua[2][2]=5.0;
  umatrix_col uv2(ua,1);
  t.test_rel(uv2[0],1.0,1.0e-12,"col 1");
  t.test_rel(uv2[1],5.0,1.0e-12,"col 2");
  t.test_rel(uv2[2],6.0,1.0e-12,"col 3");

  const double *x=0;
  gsl_vector_const_view_array(x,3);

  xmatrix xm(5,5);
  xm[2][2]=2.0;
  
  t.report();
  return 0;
}
