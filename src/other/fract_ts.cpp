/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020-2021, Andrew W. Steiner
  
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
#include <o2scl/fract.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  uniform_grid<double> grid=uniform_grid_end<double>(-2,2,199);
  
  fract f;
  table3d t3d;
  std::vector<double> roots_x, roots_y, min, max;

  roots_x.push_back(-1.0);
  roots_x.push_back(1.0);
  roots_x.push_back(0.0);
  roots_x.push_back(0.0);
  roots_y.push_back(0.0);
  roots_y.push_back(0.0);
  roots_y.push_back(-1.0);
  roots_y.push_back(1.0);
  f.nrf(nrf_z4m1,grid,grid,1000,1.0e9,t3d,roots_x,roots_y,min,max);

  hdf_file hf;
  hf.open_or_create("nrf.o2");
  hdf_output(hf,(const table3d &)t3d,"nrf");
  hf.close();

  cout << "roots: ";
  cout.setf(ios::showpos);
  for(size_t i=0;i<roots_x.size();i++) {
    if (i!=0) {
      cout << "       ";
    }
    cout << "(" << roots_x[i] << "," << roots_y[i] << ")" << endl;
  }
  cout.unsetf(ios::showpos);
  cout << "min: ";
  vector_out(cout,min,true);
  cout << "max: ";
  vector_out(cout,max,true);

  t.test_gen(roots_x.size()==min.size(),"array check");

  table3d t3db;
  size_t min2, max2;
  uniform_grid<double> grid2x=uniform_grid_end<double>(-1.7,0.5,599);
  uniform_grid<double> grid2y=uniform_grid_end<double>(-1.1,1.1,599);
  f.itf(itf_mandel,grid2x,grid2y,100,10.0,t3db,min2,max2);
  cout << "min2,max2: " << min2 << " " << max2 << endl;
  
  hf.open_or_create("nrf2.o2");
  hdf_output(hf,(const table3d &)t3db,"nrf2");
  hf.close();

  t.report();
  
  return 0;
}
