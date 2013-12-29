/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2012, Andrew W. Steiner

  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------
*/
/*
  Testing the speed of various vector types
  
  Typical output: (I think this is with range-checking on)
  array: 0
  array:    3508 -0.506366 -0.506366
  pointer:  3539 -0.506366 -0.506366
  restrict: 3506 -0.506366 -0.506366
  uvector:  3766 -0.506366 -0.506366
  ufvector: 3720 -0.506366 -0.506366
  ovector:  3719 -0.506366 -0.506366

*/

#include <iostream>
#include <o2scl/ovector_tlate.h>
#include <o2scl/uvector_tlate.h>

using namespace std;
using namespace o2scl;

int main(void) {

  // Timing
  const int N=20000;
  const int M=1000;
  
  // array
  double dv[M];
  double dv2[M];
  // pointer
  double *pv=new double[M];
  // pointer (with restrict)
  double * __restrict__ rv=new double[M];
  // uvector 1
  ufvector<M> tv;
  // uvector 2
  uvector uv(M);
  // ovector
  ovector ov(M);

  for(int k=0;k<1;k++) {

    // Note that this version and the one directly below it are
    // different. The compiler is apparently able to optimize this
    // away to zero!
    size_t s1x=clock();
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {
        dv2[(i+j)%M]=sin(((double)i+j));
      }
    }
    size_t s2x=clock();
    cout << "array: " << (s2x-s1x)/10000 << endl;

    // This one doesn't work that way, probably because we have to
    // print out a result later?
    size_t s1=clock();
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {
	dv[(i+j)%M]=sin(((double)((i+j)%M)));
      }
    }
    size_t s2=clock();
    cout.width(10);
    cout << "array: " << (s2-s1)/10000 << " " << dv[100] << " "
	 << sin(100.0) << endl;

    size_t s2b=clock();
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {
	pv[(i+j)%M]=sin(((double)((i+j)%M)));
      }
    }
    size_t s2c=clock();
    cout.width(10);
    cout << "pointer: " << (s2c-s2b)/10000 << " " << dv[100] << " "
	 << sin(100.0) << endl;

    size_t sxb=clock();
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {
	rv[(i+j)%M]=sin(((double)((i+j)%M)));
      }
    }
    size_t sxc=clock();
    cout.width(10);
    cout << "restrict: " << (sxc-sxb)/10000 << " " << dv[100] << " "
	 << sin(100.0) << endl;

    size_t s3=clock();
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {
	uv[(i+j)%M]=sin(((double)((i+j)%M)));
      }
    }
    size_t s4=clock();
    cout.width(10);
    cout << "uvector: " << (s4-s3)/10000 << " " << dv[100] << " "
	 << sin(100.0) << endl;

    size_t s5=clock();
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {
	tv[(i+j)%M]=sin(((double)((i+j)%M)));
      }
    }
    size_t s6=clock();
    cout.width(10);
    cout << "ufvector: " << (s6-s5)/10000 << " " << dv[100] << " "
	 << sin(100.0) << endl;

    size_t s7=clock();
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {
	ov[(i+j)%M]=sin(((double)((i+j)%M)));
      }
    }
    size_t s8=clock();
    cout.width(10);
    cout << "ovector: " << (s8-s7)/10000 << " " << dv[100] << " "
	 << sin(100.0) << endl;
  }

  delete[] pv;
  delete[] rv;
  return 0;
}

