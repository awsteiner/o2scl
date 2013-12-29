/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2012, Andrew W. Steiner
  
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

/* Example: ex_mandel.cpp
   -------------------------------------------------------------------
   Mandelbrot example demonstrating table3d and complex arithmetic
*/

#include <iostream>
#include <o2scl/cx_arith.h>
#include <o2scl/test_mgr.h>
#include <o2scl/table3d.h>
#ifdef O2SCL_HDF_IN_EXAMPLES
#include <o2scl/hdf_file.h>
#endif

using namespace std;
using namespace o2scl;
#ifdef O2SCL_HDF_IN_EXAMPLES
using namespace o2scl_hdf;
#endif

int main(void) {
  
  test_mgr tm;
  tm.set_output_level(2);

  // Create a table3d object
  table3d t;

  // Add parameters
  double delta=0.001, minx=-1.5, maxx=0.8, miny=-1.0, maxy=1.0;
  size_t maxtime=0, limit=100;
  t.add_constant("delta",delta);
  t.add_constant("minx",minx);
  t.add_constant("maxx",maxx);
  t.add_constant("miny",miny);
  t.add_constant("maxy",maxy);

  // Set grid
  ovector ox, oy;
  for(double x=minx;x<=maxx;x+=delta) ox.push_back(x);
  for(double y=miny;y<=maxy;y+=delta) oy.push_back(y);
  t.set_xy("x",ox.size(),ox,"y",oy.size(),oy);

  // Create slice
  t.new_slice("time");

  // Compute escape times
  for(size_t i=0;i<ox.size();i++) {
    for(size_t j=0;j<oy.size();j++) {
      gsl_complex c={{ox[i],oy[j]}};
      gsl_complex z={{0,0}};
      size_t time=0;
      for(size_t k=0;k<limit;k++) {
	// Arithmetic with gsl_complex objects
	z=z*z+c;
	if (abs(z)>10.0) {
	  time=k;
	  k=limit;
	}
      }
      t.set(i,j,"time",time);
      if (time>maxtime) maxtime=time;
    }
  }

  // Maximum escape time for color normalization
  t.add_constant("maxtime",maxtime);
  
  // Output to file if O2scl is compiled with HDF support
#ifdef O2SCL_HDF_IN_EXAMPLES
  hdf_file hf;
  hf.open_or_create("ex_mandel.o2");
  hdf_output(hf,t,"mandel");
  hf.close();
#endif

  tm.test_gen(maxtime==99,"maxtime test");
  tm.report();

  return 0;
}
// End of example
