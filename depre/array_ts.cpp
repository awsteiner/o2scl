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
#include <o2scl/array.h>
#include <o2scl/columnify.h>
#include <o2scl/test_mgr.h>
#include <climits>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  /// Data for arrays, pointers, const arrays, and const pointers
  double arr1[3]={3.0,1.0,4.0};
  double *arr2=new double[3];
  arr2[0]=3.0;
  arr2[1]=1.0;
  arr2[2]=4.0;
  const double arr3[3]={3.0,1.0,4.0};
  const double *arr4=arr2;
  
  /// Try reversing
  array_reverse<3> ar1(arr1);
  array_reverse<3> ar2(arr2);
  array_const_reverse<3> ar3(arr3);
  array_const_reverse<3> ar4(arr4);

  t.test_rel(ar1[0],4.0,1.0e-6,"rv 1");
  t.test_rel(ar1[1],1.0,1.0e-6,"rv 2");
  t.test_rel(ar1[2],3.0,1.0e-6,"rv 3");
  for(int i=0;i<3;i++) {
    cout << ar1[i] << " " << ar2[i] << " " << ar3[i] << " " << ar4[i] << endl;
    t.test_rel(ar1[i],ar2[i],1.0e-6,"1=2");
    t.test_rel(ar1[i],ar3[i],1.0e-6,"1=3");
    t.test_rel(ar1[i],ar4[i],1.0e-6,"1=4");
  }

  /// Try subvectors
  array_subvector as2(arr2,1,2);
  array_subvector as1(arr1,1,2);
  array_const_subvector as3(arr3,1,2);
  array_const_subvector as4(arr4,1,2);

  for(int i=0;i<2;i++) {
    cout << as1[i] << " " << as2[i] << " " << as3[i] << " " << as4[i] << endl;
    t.test_rel(as1[i],as2[i],1.0e-6,"1=2");
    t.test_rel(as1[i],as3[i],1.0e-6,"1=3");
    t.test_rel(as1[i],as4[i],1.0e-6,"1=4");
  }

  /// Try reversed subvectors
  array_subvector_reverse asr1(arr1,1,2);
  array_subvector_reverse asr2(arr2,1,2);
  array_const_subvector_reverse asr3(arr3,1,2);
  array_const_subvector_reverse asr4(arr4,1,2);

  for(int i=0;i<2;i++) {
    cout << asr1[i] << " " << asr2[i] << " " << asr3[i] << " " 
	 << asr4[i] << endl;
    t.test_rel(asr1[i],asr2[i],1.0e-6,"1=2");
    t.test_rel(asr1[i],asr3[i],1.0e-6,"1=3");
    t.test_rel(asr1[i],asr4[i],1.0e-6,"1=4");
  }
  
  delete arr2;
  
  cout.setf(ios::scientific);
  double a2d[3][3]={{3,-1,4},{-1,5,-9},{2,-6,5}};
  matrix_out(cout,a2d,3,3);
  cout << endl;
  
  // We need to test the sort in addition to printing out the results
  
  cout << endl;

  double mrc[5][5];
  for(size_t i=0;i<5;i++) {
    for(size_t j=0;j<5;j++) {
      mrc[i][j]=fabs(sin(((double)(i+1)))*cos(((double)(j+1))));
      cout << mrc[i][j] << " ";
    }
    cout << endl;
  }

  cout << "col: " << endl;
  array_2d_col<5,5> mrcc(mrc,2);
  for(size_t i=0;i<5;i++) {
    t.test_rel(mrcc[i],mrc[i][2],1.0e-12,"col");
    cout << mrcc[i] << " ";
  }
  cout << endl;
  
  cout << "row: " << endl;
  array_2d_row<double[5][5]> mrcr(mrc,2);
  for(size_t i=0;i<5;i++) {
    t.test_rel(mrcr[i],mrc[2][i],1.0e-12,"row");
    cout << mrcr[i] << " ";
  }
  cout << endl;

  {
    // This tests pointer_alloc garbage collection
    pointer_alloc<double> pa;
    double *p1, *p2;
    pa.allocate(p1,8);
    pa.allocate(p2,8);
    pa.free(p1);
  }

  {
    // This tests pointer_alloc garbage collection
    pointer_2d_alloc<double> pa;
    double **p1, **p2;
    pa.allocate(p1,8,6);
    pa.allocate(p2,8,4);
    pa.free(p1,8);
  }
  
  t.report();
  return 0;
}

