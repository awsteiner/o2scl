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
#include <o2scl/interp_smart.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {

  test_mgr t;
  t.set_output_level(1);

#ifdef O2SCL_NEVER_DEFINED

  ovector x(5), y(5), rx(5), ry(5);
  double xa[5], ya[5], rxa[5], rya[5];

  cout.setf(ios::scientific);

  for(size_t i=0;i<5;i++) {
    x[i]=((double)i);
    y[i]=((double)i*i);
    rx[4-i]=((double)i);
    ry[4-i]=((double)i*i);
    xa[i]=((double)i);
    ya[i]=((double)i*i);
    rxa[4-i]=((double)i);
    rya[4-i]=((double)i*i);
  }

  // ---------------------------------------------------------------
  // Test reversing arrays

  ovector_reverse rx2(x);
  for(size_t i=0;i<5;i++) {
    t.test_rel(rx[i],rx2[i],1.0e-10,"rev1");
  }
  
  array_reverse<5> rxa2(xa);
  for(size_t i=0;i<5;i++) {
    t.test_rel(rxa[i],rxa2[i],1.0e-10,"rev2");
  }

  // ---------------------------------------------------------------
  // Use linear interpolation

  interp_o2scl<ovector_base> gi(0);
  interp_sm giv1(5,x,y,0);
  interp_sm giv2(5,rx,ry,0);

  // ---------------------------------------------------------------
  // Test normal interpolation

  double x0=2.5;
  double y0;

  y0=gi.interp(x0,5,x,y);
  t.test_rel(y0,6.5,1.0e-5,"intp 1a.");
  y0=gi.deriv(x0,5,x,y);
  t.test_rel(y0,5.0,1.0e-5,"intp 1b.");
  y0=gi.deriv2(x0,5,x,y);
  t.test_rel(y0,0.0,1.0e-5,"intp 1c.");
  y0=gi.integ(x0,3.5,5,x,y);
  t.test_rel(y0,9.25,1.0e-5,"intp 1d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array
  
  y0=gi.interp(x0,5,rx,ry);
  t.test_rel(y0,6.5,1.0e-5,"rintp 1a.");
  y0=gi.deriv(x0,5,rx,ry);
  t.test_rel(y0,5.0,1.0e-5,"rintp 1b.");
  y0=gi.deriv2(x0,5,rx,ry);
  t.test_rel(y0,0.0,1.0e-5,"rintp 1c.");
  y0=gi.integ(x0,3.5,5,rx,ry);
  t.test_rel(y0,9.25,1.0e-5,"rintp 1d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with alternative interface

  y0=giv1.interp(x0);
  t.test_rel(y0,6.5,1.0e-5,"intp 2a.");
  y0=giv1.deriv(x0);
  t.test_rel(y0,5.0,1.0e-5,"intp 2b.");
  y0=giv1.deriv2(x0);
  t.test_rel(y0,0.0,1.0e-5,"intp 2c.");
  y0=giv1.integ(x0,3.5);
  t.test_rel(y0,9.25,1.0e-5,"intp 2d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array and 
  // alternate interface
  
  y0=giv2.interp(x0);
  t.test_rel(y0,6.5,1.0e-5,"rintp 2a.");
  y0=giv2.deriv(x0);
  t.test_rel(y0,5.0,1.0e-5,"rintp 2b.");
  y0=giv2.deriv2(x0);
  t.test_rel(y0,0.0,1.0e-5,"rintp 2c.");
  y0=giv2.integ(x0,3.5);
  t.test_rel(y0,9.25,1.0e-5,"rintp 2d.");

  // ---------------------------------------------------------------
  // Use cspline interpolation

  interp_o2scl<ovector_base> hi;
  interp_sm hiv1(5,x,y);
  interp_sm hiv2(5,x,y);

  // ---------------------------------------------------------------
  // Test normal interpolation

  y0=hi.interp(x0,5,x,y);
  t.test_rel(y0,6.232143,1.0e-5,"intp 3a.");
  y0=hi.deriv(x0,5,x,y);
  t.test_rel(y0,4.964286,1.0e-5,"intp 3b.");
  y0=hi.deriv2(x0,5,x,y);
  t.test_rel(y0,2.142857,1.0e-5,"intp 3c.");
  y0=hi.integ(x0,3.5,5,x,y);
  t.test_rel(y0,9.098214,1.0e-5,"intp 3d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array
  
  y0=hi.interp(x0,5,rx,ry);
  t.test_rel(y0,6.232143,1.0e-5,"rintp 3a.");
  y0=hi.deriv(x0,5,rx,ry);
  t.test_rel(y0,4.964286,1.0e-5,"rintp 3b.");
  y0=hi.deriv2(x0,5,rx,ry);
  t.test_rel(y0,2.142857,1.0e-5,"rintp 3c.");
  y0=hi.integ(x0,3.5,5,rx,ry);
  t.test_rel(y0,9.098214,1.0e-5,"rintp 3d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with alternative interface
  
  y0=hiv1.interp(x0);
  t.test_rel(y0,6.232143,1.0e-5,"intp 4a.");
  y0=hiv1.deriv(x0);
  t.test_rel(y0,4.964286,1.0e-5,"intp 4b.");
  y0=hiv1.deriv2(x0);
  t.test_rel(y0,2.142857,1.0e-5,"intp 4c.");
  y0=hiv1.integ(x0,3.5);
  t.test_rel(y0,9.098214,1.0e-5,"intp 4d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array and 
  // alternate interface
  
  y0=hiv2.interp(x0);
  t.test_rel(y0,6.232143,1.0e-5,"rintp 4a.");
  y0=hiv2.deriv(x0);
  t.test_rel(y0,4.964286,1.0e-5,"rintp 4b.");
  y0=hiv2.deriv2(x0);
  t.test_rel(y0,2.142857,1.0e-5,"rintp 4c.");
  y0=hiv2.integ(x0,3.5);
  t.test_rel(y0,9.098214,1.0e-5,"rintp 4d.");

  // ---------------------------------------------------------------
  // Now test with C-style arrays

  // ---------------------------------------------------------------
  // Use linear interpolation

  interp_array<5> gia(0);
  interp_sma<double[5]> giv1a(5,xa,ya,0);
  interp_sma<double[5]> giv2a(5,rxa,rya,0);

  // ---------------------------------------------------------------
  // Test normal interpolation

  x0=2.5;

  y0=gia.interp(x0,5,xa,ya);
  t.test_rel(y0,6.5,1.0e-5,"intp 1a.");
  y0=gia.deriv(x0,5,xa,ya);
  t.test_rel(y0,5.0,1.0e-5,"intp 1b.");
  y0=gia.deriv2(x0,5,xa,ya);
  t.test_rel(y0,0.0,1.0e-5,"intp 1c.");
  y0=gia.integ(x0,3.5,5,xa,ya);
  t.test_rel(y0,9.25,1.0e-5,"intp 1d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array
  
  y0=gia.interp(x0,5,rxa,rya);
  t.test_rel(y0,6.5,1.0e-5,"rintp 1a.");
  y0=gia.deriv(x0,5,rxa,rya);
  t.test_rel(y0,5.0,1.0e-5,"rintp 1b.");
  y0=gia.deriv2(x0,5,rxa,rya);
  t.test_rel(y0,0.0,1.0e-5,"rintp 1c.");
  y0=gia.integ(x0,3.5,5,rxa,rya);
  t.test_rel(y0,9.25,1.0e-5,"rintp 1d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with alternative interface

  y0=giv1a.interp(x0);
  t.test_rel(y0,6.5,1.0e-5,"intp 2a.");
  y0=giv1a.deriv(x0);
  t.test_rel(y0,5.0,1.0e-5,"intp 2b.");
  y0=giv1a.deriv2(x0);
  t.test_rel(y0,0.0,1.0e-5,"intp 2c.");
  y0=giv1a.integ(x0,3.5);
  t.test_rel(y0,9.25,1.0e-5,"intp 2d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array and 
  // alternate interface
  
  y0=giv2a.interp(x0);
  t.test_rel(y0,6.5,1.0e-5,"rintp 2a.");
  y0=giv2a.deriv(x0);
  t.test_rel(y0,5.0,1.0e-5,"rintp 2b.");
  y0=giv2a.deriv2(x0);
  t.test_rel(y0,0.0,1.0e-5,"rintp 2c.");
  y0=giv2a.integ(x0,3.5);
  t.test_rel(y0,9.25,1.0e-5,"rintp 2d.");

  // ---------------------------------------------------------------
  // Use cspline interpolation

  interp_array<5> hia;
  interp_sma<double[5]> hiv1a(5,xa,ya);
  interp_sma<double[5]> hiv2a(5,rxa,rya);

  // ---------------------------------------------------------------
  // Test normal interpolation

  y0=hia.interp(x0,5,xa,ya);
  t.test_rel(y0,6.232143,1.0e-5,"intp 3a.");
  y0=hia.deriv(x0,5,xa,ya);
  t.test_rel(y0,4.964286,1.0e-5,"intp 3b.");
  y0=hia.deriv2(x0,5,xa,ya);
  t.test_rel(y0,2.142857,1.0e-5,"intp 3c.");
  y0=hia.integ(x0,3.5,5,xa,ya);
  t.test_rel(y0,9.098214,1.0e-5,"intp 3d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array
  
  y0=hia.interp(x0,5,rxa,rya);
  t.test_rel(y0,6.232143,1.0e-5,"rintp 3a.");
  y0=hia.deriv(x0,5,rxa,rya);
  t.test_rel(y0,4.964286,1.0e-5,"rintp 3b.");
  y0=hia.deriv2(x0,5,rxa,rya);
  t.test_rel(y0,2.142857,1.0e-5,"rintp 3c.");
  y0=hia.integ(x0,3.5,5,rxa,rya);
  t.test_rel(y0,9.098214,1.0e-5,"rintp 3d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with alternative interface
  
  y0=hiv1a.interp(x0);
  t.test_rel(y0,6.232143,1.0e-5,"intp 4a.");
  y0=hiv1a.deriv(x0);
  t.test_rel(y0,4.964286,1.0e-5,"intp 4b.");
  y0=hiv1a.deriv2(x0);
  t.test_rel(y0,2.142857,1.0e-5,"intp 4c.");
  y0=hiv1a.integ(x0,3.5);
  t.test_rel(y0,9.098214,1.0e-5,"intp 4d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array and 
  // alternate interface
  
  y0=hiv2a.interp(x0);
  t.test_rel(y0,6.232143,1.0e-5,"rintp 4a.");
  y0=hiv2a.deriv(x0);
  t.test_rel(y0,4.964286,1.0e-5,"rintp 4b.");
  y0=hiv2a.deriv2(x0);
  t.test_rel(y0,2.142857,1.0e-5,"rintp 4c.");
  y0=hiv2a.integ(x0,3.5);
  t.test_rel(y0,9.098214,1.0e-5,"rintp 4d.");

  // ---------------------------------------------------------------

  // Test a malformed array
  err_hnd=&alt_err_hnd;
  alt_err_hnd.set_mode(0);
  
  cout << "--------------------------------------------" << endl;
  
  ovector mx(15), my(15);
  interp_o2scl<ovector_base> mi;

  err_hnd->reset();
  
  // Increasing and flat
  
  for(int i=0;i<15;i++) {
    if (i<5 || i>10) mx[i]=0.0;
    else mx[i]=((double)i-4);
    my[i]=mx[i]*mx[i];
  }
  
  cout << mi.interp(2.0,10,mx,my) << endl;
  cout << err_hnd->get_str() << endl;
  
  // Decreasing and flat
  
  for(int i=0;i<15;i++) {
    if (i<5 || i>10) mx[i]=0.0;
    else mx[i]=((double)(10-i));
    my[i]=mx[i]*mx[i];
  }
  
  cout << mi.interp(2.0,10,mx,my) << endl;
  cout << err_hnd->get_str() << endl;
  
  // Both

  for(int i=0;i<15;i++) {
    mx[i]=fabs(4.5-i);
    my[i]=mx[i]*mx[i];
  }
  
  cout << mi.interp(2.0,10,mx,my) << endl;
  cout << mi.interp(6.0,10,mx,my) << endl;

  // zero

  for(int i=0;i<15;i++) {
    mx[i]=0.0;
    my[i]=0.0;
  }

  interp_sm sizero(15,mx,my);
  cout << sizero.interp(1.0) << endl;
  cout << endl;

  // Test vector_find_level(), etc.
  {
    ovector xdat(100), ydat(100);
    for(size_t i=0;i<100;i++) {
      xdat[i]=-4+((double)i)/12.5;
      ydat[i]=exp(-xdat[i]*xdat[i])*(sin(xdat[i]-1.4)+1.0);
    }
    // Ensure endpoints are equal
    ydat[0]=0.0;
    ydat[99]=0.0;

    cout << "Testing find_level():" << endl;
    ovector locs;
    double tlev;

    vector_find_level<ovector,ovector>(0.3,100,xdat,ydat,locs);
    t.test_gen(locs.size()==0,"find level 1");
    cout << 0.3 << " " << locs.size() << " " << locs << endl;
    
    vector_find_level<ovector,ovector>(0.2,100,xdat,ydat,locs);
    t.test_gen(locs.size()==2,"find level 2");
    cout << 0.2 << " " << locs.size() << " " << locs << endl;

    vector_find_level<ovector,ovector>(0.08,100,xdat,ydat,locs);
    t.test_gen(locs.size()==4,"find level 3");
    cout << 0.08 << " " << locs.size() << " " << locs << endl;
    cout << endl;

    cout << "Testing vector_invert_enclosed_sum():" << endl;
    interp_o2scl<ovector_base> si;
    double total=si.integ(-3.95,3.95,100,xdat,ydat);
    cout << total << endl;
    double lev;

    vector_invert_enclosed_sum(total*0.1,100,xdat,ydat,lev);
    vector_find_level<ovector,ovector>(lev,100,xdat,ydat,locs);
    t.test_gen(locs.size()==2,"vector_endpoints 1");
    cout << total*0.1 << " " << locs.size() << " " << locs << endl;

    vector_invert_enclosed_sum(total*0.3,100,xdat,ydat,lev);
    vector_find_level<ovector,ovector>(lev,100,xdat,ydat,locs);
    t.test_gen(locs.size()==2,"vector_endpoints 2");
    cout << total*0.3 << " " << locs.size() << " " << locs << endl;

    vector_invert_enclosed_sum(total*0.5,100,xdat,ydat,lev);
    vector_find_level<ovector,ovector>(lev,100,xdat,ydat,locs);
    t.test_gen(locs.size()==2,"vector_endpoints 3");
    cout <<total*0.5 << " " <<  locs.size() << " " << locs << endl;

    vector_invert_enclosed_sum(total*0.7,100,xdat,ydat,lev);
    vector_find_level<ovector,ovector>(lev,100,xdat,ydat,locs);
    t.test_gen(locs.size()==4,"vector_endpoints 4");
    cout << total*0.7 << " " << locs.size() << " " << locs << endl;

    vector_invert_enclosed_sum(total*0.9,100,xdat,ydat,lev);
    vector_find_level<ovector,ovector>(lev,100,xdat,ydat,locs);
    t.test_gen(locs.size()==4,"vector_endpoints 5");
    cout << total*0.9 << " " << locs.size() << " " << locs << endl;
  }

#endif

  t.report();

  return 0;
}
