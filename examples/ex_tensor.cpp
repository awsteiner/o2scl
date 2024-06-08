/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015-2024, Andrew W. Steiner
  
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
// sphinx-example-start
/* Example: ex_tensor.cpp
   -------------------------------------------------------------------
   A simple example for tensors. See "License Information" 
   section of the documentation for license information.
*/
#include <o2scl/tensor_grid.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

double func(double x, double y, double z) {
  return x*x*x-y+z*z;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // Create a rank three tensor
  tensor_grid<> m3;
  vector<size_t> i3(3), j3(3);
  i3[0]=4;
  i3[1]=3;
  i3[2]=3;
  m3.resize(3,i3);
  
  // Create and set grid
  std::vector<double> grid;

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
  
  // ----------------------------------------------------------------
  // Demonstrate linear interpolation
  
  cout << "Interpolation: " << endl;
  
  double vals[3]={2.5,2.5,1.5};
  cout << "Exact: " << func(vals[0],vals[1],vals[2])
       << " interpolated: "
       << m3.interp_linear(vals) << endl;
  t.test_rel(func(vals[0],vals[1],vals[2]),
	     m3.interp_linear(vals),2.0e-1,"interp linear");
  
  vals[0]=3.5;
  vals[1]=2.5;
  vals[2]=1.0;
  cout << "Exact: " << func(vals[0],vals[1],vals[2])
       << " interpolated: "
       << m3.interp_linear(vals) << endl;
  t.test_rel(func(vals[0],vals[1],vals[2]),
	     m3.interp_linear(vals),2.0e-1,"interp linear");
  cout << endl;

  // ----------------------------------------------------------------
  // Demonstrate HDF5 I/O for tensors

  // Output to a file
  hdf_file hf;
  hf.open_or_create("ex_tensor_dat.o2");
  hdf_output(hf,m3,"tens_grid_test");
  hf.close();

  // Now read from that file
  tensor_grid<> m3c;
  hf.open("ex_tensor_dat.o2");
  hdf_input(hf,m3c);
  hf.close();

  // Show that the tensor is the same
  t.test_gen(m3==m3c,"tensor equality");

  // ----------------------------------------------------------------
  // Demonstrate tensor rearrangement. Create a new tensor which
  // interpolates the value 1.5 into the second index, sums over the
  // third index, and interpolates a new grid for the first index.
  // The result will be a rank 1 tensor with 10 entries. 

  cout << "grid(0),interp(1),sum(2):" << endl;
  tensor_grid<> m3r=grid_rearrange_and_copy<tensor_grid<>,double>
    (m3,{ix_grid(0,1.0,4.0,9),ix_interp(1,1.5),ix_sum(2)});
  
  t.test_gen(m3r.get_rank()==1,"rearrange 1");
  t.test_gen(m3r.get_size(0)==10,"rearrange 2");
  
  t.test_gen(m3r.get_rank()==1,"rearrange 1");
  t.test_gen(m3r.get_size(0)==10,"rearrange 2");
  
  for(size_t i=0;i<10;i++) {
    vector<size_t> ix={i};
    cout << i << " " << m3r.get(ix) << " ";
    cout << (func(1.0+((double)i)/3.0,1.5,1.0)+
             func(1.0+((double)i)/3.0,1.5,2.0)+
             func(1.0+((double)i)/3.0,1.5,3.0)) << endl;
    t.test_rel(m3r.get(ix),(func(1.0+((double)i)/3.0,1.5,1.0)+
                            func(1.0+((double)i)/3.0,1.5,2.0)+
                            func(1.0+((double)i)/3.0,1.5,3.0)),
               3.0e-1,"tensor rearrange 1");
  }
  cout << endl;
  
  cout << "interp(1),grid(0),sum(2):" << endl;
  tensor_grid<> m3s=grid_rearrange_and_copy<tensor_grid<>,double>
    (m3,{ix_interp(1,1.5),ix_grid(0,1.0,4.0,9),ix_sum(2)});
  
  for(size_t i=0;i<10;i++) {
    vector<size_t> ix={i};
    cout << i << " " << m3s.get(ix) << " ";
    cout << (func(1.0+((double)i)/3.0,1.5,1.0)+
             func(1.0+((double)i)/3.0,1.5,2.0)+
             func(1.0+((double)i)/3.0,1.5,3.0)) << endl;
    t.test_rel(m3s.get(ix),(func(1.0+((double)i)/3.0,1.5,1.0)+
                            func(1.0+((double)i)/3.0,1.5,2.0)+
                            func(1.0+((double)i)/3.0,1.5,3.0)),
               3.0e-1,"tensor rearrange 2");
  }
  cout << endl;
  
  t.report();

  return 0;
}
