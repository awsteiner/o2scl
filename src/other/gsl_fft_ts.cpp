/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <o2scl/gsl_fft.h>

using namespace std;
using namespace o2scl;

int main (void) {
  test_mgr t;
  t.set_output_level(1);
  
  int i;
  const int n=100;
  double data[n];
  
  gsl_fft rf;

  for (i=0;i<n;i++) data[i]=0.0;
  for (i=n/3;i<2*n/3;i++) data[i]=1.0;
  
  rf.transform(100,data);

  for (i=11;i<n;i++) data[i]=0;

  rf.inverse_transform(100,data);

  for (i=0;i<n;i++) {
    cout << i << " " << data[i] << endl;
    if (i<15 || i>85) t.test_rel(data[i],0.0,1.0e-1,"zero");
    if (i<55 && i>45) t.test_rel(data[i],1.0,1.0e-1,"one");
  }

  t.report();
  return 0;
}
