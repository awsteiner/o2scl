/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015, Andrew W. Steiner
  
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
#include <o2scl/tensor.h>
#include <o2scl/tensor_grid.h>
#include <o2scl/test_mgr.h>
#if O2SCL_HDF
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
using namespace o2scl_hdf;
#endif

using namespace std;
using namespace o2scl;

double func(double x, double y, double z) {
  return x*x*x-y+z*z;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // ----------------------------------------------------------------
  // Example of a tensor built out of std::vector objects

  // Create a rank three tensor
  tensor_grid<> m3;
  size_t i3[3], j3[3], k3[3];
  i3[0]=4;
  i3[1]=3;
  i3[2]=3;
  m3.resize(3,i3);
  
  // Create and set grid
  std::vector<double> grid;
  size_t j4[3];

  grid.push_back(1.0);
  grid.push_back(2.0);
  grid.push_back(3.0);
  grid.push_back(4.0);
  
  grid.push_back(1.0);
  grid.push_back(2.0);
  grid.push_back(3.0);
  
  grid.push_back(1.0);
  grid.push_back(2.0);
  grid.push_back(3.0);
  m3.set_grid_packed(grid);
  
  // Fill the tensor with data
  for(size_t i=0;i<i3[0];i++) {
    for(size_t j=0;j<i3[1];j++) {
      for(size_t k=0;k<i3[2];k++) {
	double x=m3.get_grid(0,i);
	double y=m3.get_grid(1,j);
	double z=m3.get_grid(2,k);
	j3[0]=i;
	j3[1]=j;
	j3[2]=k;
	m3.set(j3,func(x,y,z));
      }
    }
  }
  
  // Test interpolation between grid points
  cout << "Interpolation: " << endl;
  
  double vals[3]={2.5,2.5,1.5};
  cout << "Exact: " << func(vals[0],vals[1],vals[2])
       << " interpolated: "
       << m3.interp_linear(vals) << endl;
  t.test_rel(func(vals[0],vals[1],vals[2]),
	     m3.interp_linear(vals),2.0e-1,"interp linear");

  vals[0]=2.0;
  vals[1]=3.0;
  vals[2]=1.0;
  cout << "Exact: " << func(vals[0],vals[1],vals[2])
       << " interpolated: "
       << m3.interp_linear(vals) << endl;
  t.test_rel(func(vals[0],vals[1],vals[2]),
	     m3.interp_linear(vals),2.0e-1,"interp linear");
  cout << endl;

  // ----------------------------------------------------------------
  // Example of a tensor built out of ublas vector objects
  // The interpolate() function works easily with this type
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
  
  tensor_grid<ubvector,ubvector_size_t> m3b;

  // Size and set grid
  m3b.resize(3,i3);
  m3b.set_grid_packed(grid);
  
  // Fill the tensor with data
  for(size_t i=0;i<i3[0];i++) {
    for(size_t j=0;j<i3[1];j++) {
      for(size_t k=0;k<i3[2];k++) {
	double x=m3b.get_grid(0,i);
	double y=m3b.get_grid(1,j);
	double z=m3b.get_grid(2,k);
	j3[0]=i;
	j3[1]=j;
	j3[2]=k;
	m3b.set(j3,func(x,y,z));
      }
    }
  }

  // Test interpolation between grid points
  cout << "Interpolation: " << endl;

  vals[0]=2.5;
  vals[1]=2.5;
  vals[2]=1.5;
  cout << "Exact: " << func(vals[0],vals[1],vals[2])
       << " linear: "
       << m3b.interp_linear(vals)
       << " linear(2): "
       << m3b.interpolate(vals) << endl;
  m3b.set_interp_type(itp_cspline);
  cout << " cubic spline: "
       << m3b.interpolate(vals) << endl;
  cout << endl;
  
  t.report();

  return 0;
}
