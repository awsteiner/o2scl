/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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

/// For the definition of M_PI
#include <gsl/gsl_math.h>

#include <o2scl/series_acc.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  const int N=20;
  double arr[N], sum, err;
  series_acc gs(N);

  cout.setf(ios::scientific);
  cout.precision(14);

  const double zeta_2=M_PI*M_PI/6.0;
  for (int n=0;n<N;n++) {
    double np1=n+1.0;
    arr[n]=1.0/(np1*np1);
  }
  
  sum=gs.series_accel(N,arr,err);
  t.test_rel(sum,zeta_2,err,"lewin, N=20");
  cout << sum << " " << " " << err << " " << zeta_2 << endl;

  cout << "Test changing the size:" << endl;
  gs.set_size(N-10);
  sum=gs.series_accel(N-10,arr,err);
  t.test_rel(sum,zeta_2,err,"lewin, N=10");
  cout << sum << " " << " " << err << " " << zeta_2 << endl;

  t.report();
  return 0;
}
