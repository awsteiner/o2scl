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
#include <iostream>

using namespace std;

const int N=1000;

template<class vec_t> int func1(vec_t &x) {
  int i, j;
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      x[(i+j)%5]=((double)i+j);
    }
  }
  return 0;
}

int func2(double x[]) {
  int i, j;
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      x[(i+j)%5]=((double)i+j);
    }
  }
  return 0;
}

int func3(double *x) {
  int i, j;
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      x[(i+j)%5]=((double)i+j);
    }
  }
  return 0;
}

int main(void) {
  size_t s;
  double x[5];
  double *xh=new double[5];
  
  s=clock();
  for(int k=0;k<300;k++) {
    func1<double [5]>(x);
  }
  cout << clock()-s << endl;

  s=clock();
  for(int k=0;k<300;k++) {
    func1<double *>(xh);
  }
  cout << clock()-s << endl;

  s=clock();
  for(int k=0;k<300;k++) {
    func2(x);
  }
  cout << clock()-s << endl;
  
  s=clock();
  for(int k=0;k<300;k++) {
    func3(xh);
  }
  cout << clock()-s << endl;
  
  s=clock();
  for(int k=0;k<300;k++) {
    func3(xh);
  }
  cout << clock()-s << endl;
  
  delete[] xh;

  return 0;
}
