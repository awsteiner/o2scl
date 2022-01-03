/*
  -------------------------------------------------------------------
  
  Copyright (C) 2010-2022, Andrew W. Steiner
  
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
#include <o2scl/hist_2d.h>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/convert_units.h>
#include <o2scl/rng.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  
  cout.setf(ios::scientific);

  hist_2d h;
  h.set_bin_edges(uniform_grid_width<>(1.0,0.1,10),
		  uniform_grid_width<>(0.0,1.0,10));
  
  for(size_t i=0;i<10;i++) {
    cout << i << " " << h.get_x_low_i(i) << " " << h.get_x_high_i(i) << " "
	 << h.get_y_low_i(i) << " " << h.get_y_high_i(i) << endl;
  }
  cout << endl;

  rng<> gr;
  for(size_t i=0;i<10000;i++) {
    h.update(gr.random()*gr.random()+1.0,gr.random()*gr.random()*9.0);
  }
  
  t.report();
  return 0;
}

