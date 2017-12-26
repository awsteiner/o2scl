/*
  -------------------------------------------------------------------
  
  Copyright (C) 2010-2018, Andrew W. Steiner
  
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
#include <o2scl/hist.h>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/convert_units.h>
#include <o2scl/rng_gsl.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;

  t.set_output_level(1);
  
  cout.setf(ios::scientific);

  hist h;
  h.set_bin_edges(uniform_grid_width<>(1.0,0.1,10));
  
  for(size_t i=0;i<10;i++) {
    cout << i << " " << h.get_bin_low_i(i) << " " 
	 << h.get_bin_high_i(i) << endl;
  }
  cout << endl;

  rng_gsl gr;
  for(size_t i=0;i<10000;i++) {
    h.update(gr.random()*gr.random()+1.0);
  }

  for(size_t i=0;i<10;i++) {
    cout << i << " " << h.get_rep_i(i) << " " << h.get_wgt_i(i) << endl;
  }
  cout << endl;
  
  for(double x=1.0;x<=2.001;x+=0.15) {
    cout << x << " " << h.get_rep(x) << " " 
	 << h(x) << " " << h.interp(x) << " "
	 << h.deriv(x) << endl;
    t.test_rel(h(x),h.interp(x),1.0e-10,"func and interp");
  }
  cout << endl;

  vector<double> x;
  for(size_t i=0;i<1000;i++) {
    x.push_back(sin(gr.random()));
  }
  hist h2(x.size(),x,10);
  for(size_t i=0;i<10;i++) {
    cout << i << " " << h2.get_rep_i(i) << " " << h2[i] << endl;
  }
  cout << h2.sum_wgts() << endl;
  
  t.report();
  return 0;
}

