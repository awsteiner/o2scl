/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/test_mgr.h>
#include <o2scl/table.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_iv_table.h>

using namespace std;
using namespace o2scl;

int derivs(double x, size_t nv, const ubvector &y, ubvector &dydx) {
  dydx[0]=y[1];
  dydx[1]=-y[0];
  return 0;
}

typedef double arr_t[1];

int derivsa(double x, size_t nv, const arr_t &y, arr_t &dydx) {
  dydx[0]=y[1];
  dydx[1]=-y[0];
  return 0;
}

int main(void) {
  test_mgr t;

  t.set_output_level(1);

  cout.setf(ios::scientific);

#ifdef O2SCL_NEVER_DEFINED

  // Function and ODE types
  ode_funct_fptr<> od(derivs);
  ode_iv_table<> oit;
  ode_funct_fptr<arr_t> oda(derivsa);
  ode_iv_table<ode_funct<arr_t>,arr_t,arr_t,array_alloc<arr_t> > oita;

  // ------------------------------------------------

  {
    const size_t ngrid=100;

    table tab;
    tab.set_nlines(ngrid);
    tab.new_column("x");
    for(size_t i=0;i<ngrid;i++) {
      tab.set("x",i,((double)i)/(ngrid-1)*5.0);
    }
    ubvector ystart(2);
    ystart[0]=0.5;
    ystart[1]=1.0;
    oit.solve_grid_table(2,ystart,tab,"x","y","dydx","yerr",od);
    tab.summary(&cout);
  
    cout.precision(3);
    cout.setf(ios::showpos);
    for(size_t i=0;i<ngrid;i+=9) {
      cout << tab.get("x",i) << " " << tab.get("y0",i) << " " 
	   << tab.get("y1",i) << " " << tab.get("yerr0",i) << " " 
	   << tab.get("yerr1",i) << " " << tab.get("dydx0",i) << " " 
	   << tab.get("dydx1",i) << endl;
    }
    cout.precision(6);
    cout.unsetf(ios::showpos);
  }

  // ------------------------------------------------

  {
    table tab;
    ubvector ystart(2);
    ystart[0]=0.5;
    ystart[1]=1.0;
    size_t n_sol=100;
    oit.solve_store_table(0.0,5.0,0.05,2,ystart,n_sol,tab,
			  "x","y","dydx","yerr",od);
    tab.set_nlines(n_sol);
    tab.summary(&cout);
  
    cout.precision(3);
    cout.setf(ios::showpos);
    for(size_t i=0;i<tab.get_nlines();i++) {
      cout << tab.get("x",i) << " " << tab.get("y0",i) << " " 
	   << tab.get("y1",i) << " " << tab.get("yerr0",i) << " " 
	   << tab.get("yerr1",i) << " " << tab.get("dydx0",i) << " " 
	   << tab.get("dydx1",i) << endl;
    }
    cout.precision(6);
    cout.unsetf(ios::showpos);
  }

  // ------------------------------------------------

#endif

  t.report();
  return 0;
}

