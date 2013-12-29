/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/table3d.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
#if O2SCL_HDF_SVAR
#include <o2scl/hdf_file.h>
using namespace o2scl_hdf;
#endif

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  table3d at;
  at.set_size(2,4);
  at.set_size(3,5);

  // Set the grid size
  size_t nx1, ny1;
  at.get_size(nx1,ny1);
  cout << "nx: " << nx1 << " ny: " << ny1 << endl;
    
  // Add some slices
  at.new_slice("z1");
  at.new_slice("z2");

  at.get_size(nx1,ny1);
  t.test_gen(nx1==3,"size 1");
  t.test_gen(ny1==5,"size 2");
  cout << "nx: " << nx1 << " ny: " << ny1 << endl;
  cout << endl;

  // Add the data
  for(size_t i=0;i<nx1;i++) {
    for(size_t j=0;j<ny1;j++) {
      at.set(i,j,"z1",sqrt(((double)(i+j))));
      at.set(i,j,"z2",sqrt(2.0*((double)(i+j))));
    }
  }
    
  // Definite the grid
  double xv[3]={0,1,2};
  double yv[5]={0,1,2,3,4};
  at.set_xy<double *>("x",3,xv,"y",5,yv);

  // Output the table3d object
  cout.precision(5);

  cout << " z1:        ";
  for(size_t i=0;i<5;i++) {
    cout << yv[i] << " ";
  }
  cout << endl;
  for(size_t i=0;i<nx1;i++) {
    cout << xv[i] << " ";
    for(size_t j=0;j<ny1;j++) {
      cout << at.get(i,j,"z1") << " ";
    }
    cout << endl;
  }
  cout << endl;
    
  cout << " z2:        ";
  for(size_t i=0;i<5;i++) {
    cout << yv[i] << " ";
  }
  cout << endl;
  for(size_t i=0;i<nx1;i++) {
    cout << xv[i] << " ";
    for(size_t j=0;j<ny1;j++) {
      cout << at.get(i,j,"z2") << " ";
    }
    cout << endl;
  }
  cout << endl;

  // Test the interpolation
  cout << "1.5,0.5,z2: " << at.interp(1.5,0.5,"z2") << endl;
  cout << "1.9,0.75,z2: " << at.interp(1.9,0.75,"z2") << endl;
  cout << "2.0,1.0,z2: " << at.interp(2.0,1.0,"z2") << endl;
  cout << "2.1,1.25,z2: " << at.interp(2.1,1.25,"z2") << endl;
  cout << "2.5,1.5,z2: " << at.interp(2.5,1.5,"z2") << endl;
  cout << endl;

  // Test extract functions
  table<> tx, ty;
  at.extract_x(1.0,tx);

  for(size_t i=0;i<tx.get_ncolumns();i++) {
    cout.width(12);
    cout << tx.get_column_name(i) << " ";
  }
  cout << endl;
  for(size_t j=0;j<tx.get_nlines();j++) {
    for(size_t i=0;i<tx.get_ncolumns();i++) {
      cout.width(12);
      cout << tx.get(i,j) << " ";
    }
    cout << endl;
  }
  cout << endl;

  at.extract_y(2.0,ty);

  for(size_t i=0;i<ty.get_ncolumns();i++) {
    cout.width(12);
    cout << ty.get_column_name(i) << " ";
  }
  cout << endl;
  for(size_t j=0;j<ty.get_nlines();j++) {
    for(size_t i=0;i<ty.get_ncolumns();i++) {
      cout.width(12);
      cout << ty.get(i,j) << " ";
    }
    cout << endl;
  }
  cout << endl;

  // Test of table3d I/O
  {
    table3d bt, bt2;

    // Set the grid size
    bt.set_size(3,5);
    size_t nnx1, nny1;
    bt.get_size(nnx1,nny1);
    
    // Add some slices
    bt.new_slice("z1");
    bt.new_slice("z2");
    
    // Add the data
    for(size_t i=0;i<nnx1;i++) {
      for(size_t j=0;j<nny1;j++) {
	bt.set(i,j,"z1",sqrt(((double)(i+j))));
	bt.set(i,j,"z2",sqrt(2.0*((double)(i+j))));
      }
    }
    
#if O2SCL_HDF_SVAR

    hdf_file hf;
    hf.open_or_create("table3d.o2");
    hdf_output(hf,bt,"table3d_test");
    hf.close();
    hf.open("table3d.o2");
    hdf_input(hf,bt2,"table3d_test");
    hf.close();

    t.test_gen(bt.get_nx()==bt2.get_nx(),"nx");
    t.test_gen(bt.get_ny()==bt2.get_ny(),"ny");
    t.test_gen(bt.get_nslices()==bt2.get_nslices(),"cols");
    t.test_gen(bt.get_nconsts()==bt2.get_nconsts(),"cols");

#endif

  }    

  {
    table3d atf;
    ubvector x(4), y(4);
    for(size_t i=0;i<4;i++) {
      x[i]=i;
      y[i]=i;
    }
    atf.set_xy("x",4,x,"y",4,y);
    atf.new_slice("z1");
    atf.new_slice("z2");
    for(size_t i=0;i<4;i++) {
      for(size_t j=0;j<4;j++) {
	atf.set(i,j,"z1",sqrt(((double)(i+j))));
      }
    }
    atf.init_slice("z2",0.0);

    atf.summary(&cout);

    // Test make_slice().
    atf.function_slice("sin(z1)+z2+1.0","z3");
    for(size_t i=0;i<4;i++) {
      for(size_t j=0;j<4;j++) {
	cout << atf.get(i,j,"z3") << " ";
	t.test_rel(atf.get(i,j,"z3"),sin(sqrt((double)(i+j)))+1.0,1.0e-8,
		   "make_slice() results.");
      }
      cout << endl;
    }
    cout << endl;
  }

  t.report();

  return 0;
}


