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

/* Example: ex_contour.cpp
   -------------------------------------------------------------------
   Example for generating contour lines
*/

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/contour.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

// A function defining the three-dimensional surface
// for which we want to compute contour levels
double fun(double x, double y) {
  return 15.0*exp(-pow(x-20.0,2.0)/400.0-pow(y-5.0,2.0)/25.0)
    +40.0*exp(-pow(x-70.0,2.0)/500.0-pow(y-2.0,2.0)/4.0);
}

// A function for outputting the data to cout
int print_data(int nx, int ny, ubvector &x, ubvector &y,
	       ubmatrix &data);

// A function for printing the contour information to a file
int file_out(string prefix, ubvector &x, ubvector &y,
	     ubmatrix &data, vector<contour_line> &conts, 
	     vector<edge_crossings> &xecs,
	     vector<edge_crossings> &yecs);

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  contour co;

  // Initialize the data

  ubvector x(12), y(10);
  ubmatrix data(12,10);
  for(size_t i=0;i<10;i++) {
    y[i]=((double)i);
  }
  for(size_t i=0;i<12;i++) {
    x[i]=((double)i)*((double)i);
  }
  
  for(size_t j=0;j<12;j++) {
    for(size_t k=0;k<10;k++) {
      data(j,k)=fun(x[j],y[k]);
    }
  }
  co.set_data(12,10,x,y,data);
  
  // Print out the data

  print_data(12,10,x,y,data);

  // Set the contour levels

  ubvector levels(7);
  levels[0]=5.0;
  levels[1]=10.0;
  levels[2]=0.002;
  levels[3]=20.0;
  levels[4]=0.2;
  levels[5]=30.0;
  levels[6]=2.0;
  
  co.set_levels(7,levels);

  // Compute the contours

  vector<contour_line> conts;
  co.calc_contours(conts);

  vector<edge_crossings> xecs, yecs;
  co.get_edges(xecs,yecs);

  // Output to a file

  file_out("c1",x,y,data,conts,xecs,yecs);

  // Print the contours to the screen and test to make sure
  // that they match the requested level

  size_t nc=conts.size();
  for(size_t i=0;i<nc;i++) {
    cout << "Contour " << i << " at level " << conts[i].level << ":" << endl;
    size_t cs=conts[i].x.size();
    for(size_t j=0;j<cs;j++) {
      cout << "(" << conts[i].x[j] << ", " << conts[i].y[j] << ") " 
	   << fun(conts[i].x[j],conts[i].y[j]) << endl;
      t.test_rel(fun(conts[i].x[j],conts[i].y[j]),conts[i].level,
		 1.0,"curve");
    }
    cout << endl;
  }

  // Refine the data using cubic spline interpolation

  co.regrid_data(5,5);
  ubvector *x2;
  ubvector *y2;
  ubmatrix *data2;
  size_t sx, sy;
  co.get_data(sx,sy,x2,y2,data2);

  // Recompute the contours

  vector<contour_line> conts2;
  co.calc_contours(conts2);

  vector<edge_crossings> xecs2, yecs2;
  co.get_edges(xecs2,yecs2);

  // Output to a file

  file_out("c2",*x2,*y2,*data2,conts2,xecs2,yecs2);

  t.report();

  return 0;
}
// End of example

int file_out(string prefix, ubvector &x, ubvector &y,
	     ubmatrix &data, vector<contour_line> &conts, 
	     vector<edge_crossings> &xecs,
	     vector<edge_crossings> &yecs) {
  
  hdf_file hf;
  hf.open_or_create("ex_contour.o2");
  hf.setd_vec_copy(prefix+"_x",x);
  hf.setd_vec_copy(prefix+"_y",y);
  hf.setd_mat_copy(prefix+"_data",data);
  hdf_output(hf,conts,prefix+"_cl");
  hdf_output(hf,xecs,prefix+"_xe");
  hdf_output(hf,yecs,prefix+"_ye");
  hf.close();

  return 0;
}
  
int print_data(int nx, int ny, ubvector &x, ubvector &y,
	       ubmatrix &data) {
  int j, k;
  cout.setf(ios::showpos);
  cout.precision(3);
  cout << "              ";
  for(k=0;k<ny;k++) {
    string stx="y";
    stx+=itos(k);
    cout.width(11);
    cout << stx;
  }
  cout << endl;
  cout << "               ";
  for(k=0;k<ny;k++) {
    cout << y[k] << " ";
  }
  cout << endl;
  for(j=0;j<nx;j++) {
    string stx="x";
    stx+=itos(j);
    cout.width(3);
    cout << stx << " ";
    cout << x[j] << " ";
    for(k=0;k<ny;k++) {
      cout << data(j,k) << " ";
    }
    cout << endl;
  }
  cout << endl;
  cout.unsetf(ios::showpos);
  cout.precision(6);
  return 0;
}

