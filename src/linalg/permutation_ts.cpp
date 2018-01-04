/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2018, Andrew W. Steiner

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

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/permutation.h>

using namespace std;
using namespace o2scl;

int main(void) {

  test_mgr t;
  t.set_output_level(1);
  
  cout.setf(ios::scientific);

  typedef boost::numeric::ublas::vector<double> ubvector;

  // Init

  ubvector a(5), c(5);
  permutation p(5);

  p.set(0,3);
  p.set(1,2);
  p.set(2,4);
  p.set(3,0);
  p.set(4,1);
  for(size_t i=0;i<5;i++) {
    cout << p[i] << " " << p(i) << endl;
  }

  // Apply permutation

  a[0]=3.0;
  a[1]=1.0;
  a[2]=2.0;
  a[3]=8.0;
  a[4]=5.0;

  p.apply(a);

  // Apply inverse permutation

  a[0]=3.0;
  a[1]=1.0;
  a[2]=2.0;
  a[3]=8.0;
  a[4]=5.0;

  p.apply_inverse(a);
  permutation ip=p.inverse();

  c[0]=3.0;
  c[1]=1.0;
  c[2]=2.0;
  c[3]=8.0;
  c[4]=5.0;

  ip.apply(c);
  t.test_rel_vec(5,a,c,1.0e-10,"inverse permutation");

  t.report();
  return 0;
}

