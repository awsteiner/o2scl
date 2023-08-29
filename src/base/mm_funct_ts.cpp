/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/mm_funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/funct_to_fp.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  vector<string> list(2), funcs(2);
  list[0]="x1";
  list[1]="x2";
  funcs[0]="x1*x2-x1";
  funcs[1]="x1*x2-x2";
  mm_funct_strings f(2,funcs,list);
  ubvector x(2), y(2);
  x[0]=1.5;
  x[1]=2.1;
  f(2,x,y);
  t.test_rel(x[0]*x[1]-x[0],y[0],1.0e-12,"mm_funct_strings");
  t.test_rel(x[0]*x[1]-x[1],y[1],1.0e-12,"mm_funct_strings");

#ifdef O2SCL_PYTHON

  o2scl_settings.py_init();
  o2scl_settings.add_python_path("../../data/o2scl/python/");
  {
    // We use the brackets to force the mm_funct_python
    // destructor to run before py_final()
    
    ubvector xxa(2), yya(2), y2(2);
    y2[0]=2.0*o2scl_const::pi;
    y2[1]=3.0*o2scl_const::pi;

    cout << "Calling fun() in mm_funct_test.py:" << endl;
    mm_funct_python<ubvector> fp("mm_funct_test","fun");
    xxa[0]=2.0;
    xxa[1]=3.0;
    int mfp_ret=fp(2,xxa,yya);
    t.test_rel_vec(2,yya,y2,1.0e-12,"mm_funct_python");
    cout << endl;

    for(size_t j=0;j<10;j++) {
      xxa[0]=sin(j);
      xxa[1]=cos(j);
      fp(2,xxa,yya);
      t.test_rel(yya[0],xxa[0]*o2scl_const::pi,1.0e-10,"mm_funct_python 1");
      t.test_rel(yya[1],xxa[1]*o2scl_const::pi,1.0e-10,"mm_funct_python 2");
    }
    cout << endl;

    mm_funct_python_ndarray<std::vector<double> >
      fp4("mm_funct_test","fun_numpy","mft",0);
    cout << "Calling mft.fun_numpy() in mm_funct_test.py:" << endl;
    std::vector<double> sx(2), sy(2), sy2(2);
    sx[0]=2.0;
    sx[1]=3.0;
    sy2[0]=2.0*o2scl_const::pi;
    sy2[1]=3.0*o2scl_const::pi;
    int mfcp_ret=fp4(2,sx,sy);
    t.test_rel_vec(2,sy,sy2,1.0e-12,"mm_funct_python_ndarray");
    cout << endl;

    for(size_t j=0;j<10;j++) {
      sx[0]=sin(j);
      sx[1]=cos(j);
      fp4(2,sx,sy);
      t.test_rel(sy[0],sx[0]*o2scl_const::pi,1.0e-10,
                 "mm_funct_python_ndarray 1");
      t.test_rel(sy[1],sx[1]*o2scl_const::pi,1.0e-10,
                 "mm_funct_python_ndarray 2");
    }
    cout << endl;

    
  }
  o2scl_settings.py_final();
  
#endif
  
  t.report();
  return 0;
}

