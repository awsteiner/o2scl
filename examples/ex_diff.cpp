/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/funct.h>
#include <o2scl/table.h>
#include <o2scl/deriv_eqi.h>
#include <o2scl/interp.h>
#include <o2scl/test_mgr.h>
#ifdef O2SCL_HDF
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#endif

/* Example: ex_diff.cpp
   -------------------------------------------------------------------
   A simple example for showing computation of derivatives with
   interpolation from the table class and using the deriv_eqi class.
   This also shows simple I/O with the table.
*/

using namespace std;
using namespace o2scl;
#ifdef O2SCL_HDF
using namespace o2scl_hdf;
#endif

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  
  size_t bign=30, i;
  
  ubvector x(bign), y(bign), dydx(bign), dydx2(bign);

  // We create some data, y=exp(x).
  for(i=0;i<bign;i++) {
    x[i]=((double)i)/3.0;
    y[i]=exp(x[i]);
  }

  // Compute dydx using derivative formulas for equally-spaced abscissas
  deriv_eqi<> ei;
  ei.deriv_vector(bign,1.0/3.0,y,dydx);

  // Compute the derivative using the interpolation from o2scl_interp_vec
  interp_vec<> gi(bign,x,y);
  for(i=0;i<bign;i++) {
    dydx2[i]=gi.deriv(x[i]);
  }

  // Create a data table 
  table<> ta(100);
  // Give the column names
  ta.line_of_names("x y dydx dydx2 err err2");

  cout << "x            y            dydx         dydx'"
       << "        err          err'" << endl;
  for(i=0;i<bign;i++) {
    cout << x[i] << " " << y[i] << " " << dydx[i] << " " << dydx2[i] << " ";
    cout << fabs(y[i]-dydx[i])/y[i] << " " << fabs(y[i]-dydx2[i])/y[i] << endl;
    // Add the data to the table as well...
    double line[6]={x[i],y[i],dydx[i],dydx2[i],
		    fabs(y[i]-dydx[i])/y[i],fabs(y[i]-dydx2[i])/y[i]};
    ta.line_of_data(6,line);
  }
  
#ifdef O2SCL_HDF
  hdf_file hf;
  hf.open_or_create("ex_diff.o2");
  o2scl_hdf::hdf_output(hf,ta,"table");
  hf.close();
#endif

  t.report();
  return 0;
}
// End of example
