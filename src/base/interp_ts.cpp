/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
#include <o2scl/interp.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // ---------------------------------------------------------------
  // Create test data

  ubvector x(5), y(5), rx(5), ry(5);
  double xa[5], ya[5], rxa[5], rya[5];
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
  // Use linear interpolation

  interp<> gi(itp_linear);
  interp_vec<> giv1(5,x,y,itp_linear);

  // ---------------------------------------------------------------
  // Test normal interpolation

  double x0=2.5;
  double y0;
  
  y0=gi.eval(x0,5,x,y);
  t.test_rel(y0,6.5,1.0e-5,"intp 1a.");
  y0=gi.deriv(x0,5,x,y);
  t.test_rel(y0,5.0,1.0e-5,"intp 1b.");
  y0=gi.deriv2(x0,5,x,y);
  t.test_rel(y0,0.0,1.0e-5,"intp 1c.");
  y0=gi.integ(x0,3.5,5,x,y);
  t.test_rel(y0,9.25,1.0e-5,"intp 1d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array
  
  y0=gi.eval(x0,5,rx,ry);
  t.test_rel(y0,6.5,1.0e-5,"rintp 1a.");
  y0=gi.deriv(x0,5,rx,ry);
  t.test_rel(y0,5.0,1.0e-5,"rintp 1b.");
  y0=gi.deriv2(x0,5,rx,ry);
  t.test_rel(y0,0.0,1.0e-5,"rintp 1c.");
  y0=gi.integ(x0,3.5,5,rx,ry);
  t.test_rel(y0,9.25,1.0e-5,"rintp 1d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with alternative interface

  y0=giv1.eval(x0);
  t.test_rel(y0,6.5,1.0e-5,"intp 2a.");
  y0=giv1.deriv(x0);
  t.test_rel(y0,5.0,1.0e-5,"intp 2b.");
  y0=giv1.deriv2(x0);
  t.test_rel(y0,0.0,1.0e-5,"intp 2c.");
  y0=giv1.integ(x0,3.5);
  t.test_rel(y0,9.25,1.0e-5,"intp 2d.");

  // ---------------------------------------------------------------
  // Use cspline interpolation

  gi.set_type(itp_cspline);
  interp_vec<> giv3(5,x,y,itp_cspline);
  interp_vec<> giv4(5,rx,ry,itp_cspline);

  // ---------------------------------------------------------------
  // Test normal interpolation

  y0=gi.eval(x0,5,x,y);
  t.test_rel(y0,6.232143,1.0e-5,"intp 3a.");
  y0=gi.deriv(x0,5,x,y);
  t.test_rel(y0,4.964286,1.0e-5,"intp 3b.");
  y0=gi.deriv2(x0,5,x,y);
  t.test_rel(y0,2.142857,1.0e-5,"intp 3c.");
  y0=gi.integ(x0,3.5,5,x,y);
  t.test_rel(y0,9.098214,1.0e-5,"intp 3d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array
  
  y0=gi.eval(x0,5,rx,ry);
  t.test_rel(y0,6.232143,1.0e-5,"rintp 3a.");
  y0=gi.deriv(x0,5,rx,ry);
  t.test_rel(y0,4.964286,1.0e-5,"rintp 3b.");
  y0=gi.deriv2(x0,5,rx,ry);
  t.test_rel(y0,2.142857,1.0e-5,"rintp 3c.");
  y0=gi.integ(x0,3.5,5,rx,ry);
  t.test_rel(y0,9.098214,1.0e-5,"rintp 3d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with alternative interface
  
  y0=giv3.eval(x0);
  t.test_rel(y0,6.232143,1.0e-5,"intp 4a.");
  y0=giv3.deriv(x0);
  t.test_rel(y0,4.964286,1.0e-5,"intp 4b.");
  y0=giv3.deriv2(x0);
  t.test_rel(y0,2.142857,1.0e-5,"intp 4c.");
  y0=giv3.integ(x0,3.5);
  t.test_rel(y0,9.098214,1.0e-5,"intp 4d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array and 
  // alternate interface
  
  y0=giv4.eval(x0);
  t.test_rel(y0,6.232143,1.0e-5,"rintp 4a.");
  y0=giv4.deriv(x0);
  t.test_rel(y0,4.964286,1.0e-5,"rintp 4b.");
  y0=giv4.deriv2(x0);
  t.test_rel(y0,2.142857,1.0e-5,"rintp 4c.");
  y0=giv4.integ(x0,3.5);
  t.test_rel(y0,9.098214,1.0e-5,"rintp 4d.");

  // ---------------------------------------------------------------
  // Testing with C-style arrays

  // ---------------------------------------------------------------
  // Use linear interpolation

  interp_array<5> gia(itp_linear);
  interp_array_vec<double[5]> giv1a(5,xa,ya,itp_linear);

  // ---------------------------------------------------------------
  // Test normal interpolation
  
  x0=2.5;
  
  y0=gia.eval(x0,5,xa,ya);
  t.test_rel(y0,6.5,1.0e-5,"intp 1a.");
  y0=gia.deriv(x0,5,xa,ya);
  t.test_rel(y0,5.0,1.0e-5,"intp 1b.");
  y0=gia.deriv2(x0,5,xa,ya);
  t.test_rel(y0,0.0,1.0e-5,"intp 1c.");
  y0=gia.integ(x0,3.5,5,xa,ya);
  t.test_rel(y0,9.25,1.0e-5,"intp 1d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array
  
  y0=gia.eval(x0,5,rxa,rya);
  t.test_rel(y0,6.5,1.0e-5,"rintp 1a.");
  y0=gia.deriv(x0,5,rxa,rya);
  t.test_rel(y0,5.0,1.0e-5,"rintp 1b.");
  y0=gia.deriv2(x0,5,rxa,rya);
  t.test_rel(y0,0.0,1.0e-5,"rintp 1c.");
  y0=gia.integ(x0,3.5,5,rxa,rya);
  t.test_rel(y0,9.25,1.0e-5,"rintp 1d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with alternative interface

  y0=giv1a.eval(x0);
  t.test_rel(y0,6.5,1.0e-5,"intp 2a.");
  y0=giv1a.deriv(x0);
  t.test_rel(y0,5.0,1.0e-5,"intp 2b.");
  y0=giv1a.deriv2(x0);
  t.test_rel(y0,0.0,1.0e-5,"intp 2c.");
  y0=giv1a.integ(x0,3.5);
  t.test_rel(y0,9.25,1.0e-5,"intp 2d.");

  // ---------------------------------------------------------------
  // Use cspline interpolation
  
  gia.set_type(itp_cspline);
  interp_array_vec<double[5]> giv3a(5,xa,ya,itp_cspline);
  interp_array_vec<double[5]> giv4a(5,rxa,rya,itp_cspline);
  
  // ---------------------------------------------------------------
  // Test normal interpolation

  y0=gia.eval(x0,5,xa,ya);
  t.test_rel(y0,6.232143,1.0e-5,"intp 3a.");
  y0=gia.deriv(x0,5,xa,ya);
  t.test_rel(y0,4.964286,1.0e-5,"intp 3b.");
  y0=gia.deriv2(x0,5,xa,ya);
  t.test_rel(y0,2.142857,1.0e-5,"intp 3c.");
  y0=gia.integ(x0,3.5,5,xa,ya);
  t.test_rel(y0,9.098214,1.0e-5,"intp 3d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array
  
  y0=gia.eval(x0,5,rxa,rya);
  t.test_rel(y0,6.232143,1.0e-5,"rintp 3a.");
  y0=gia.deriv(x0,5,rxa,rya);
  t.test_rel(y0,4.964286,1.0e-5,"rintp 3b.");
  y0=gia.deriv2(x0,5,rxa,rya);
  t.test_rel(y0,2.142857,1.0e-5,"rintp 3c.");
  y0=gia.integ(x0,3.5,5,rxa,rya);
  t.test_rel(y0,9.098214,1.0e-5,"rintp 3d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with alternative interface
  
  y0=giv3a.eval(x0);
  t.test_rel(y0,6.232143,1.0e-5,"intp 4a.");
  y0=giv3a.deriv(x0);
  t.test_rel(y0,4.964286,1.0e-5,"intp 4b.");
  y0=giv3a.deriv2(x0);
  t.test_rel(y0,2.142857,1.0e-5,"intp 4c.");
  y0=giv3a.integ(x0,3.5);
  t.test_rel(y0,9.098214,1.0e-5,"intp 4d.");

  // ---------------------------------------------------------------
  // Test normal interpolation with a reversed array and 
  // alternate interface
  
  y0=giv4a.eval(x0);
  t.test_rel(y0,6.232143,1.0e-5,"rintp 4a.");
  y0=giv4a.deriv(x0);
  t.test_rel(y0,4.964286,1.0e-5,"rintp 4b.");
  y0=giv4a.deriv2(x0);
  t.test_rel(y0,2.142857,1.0e-5,"rintp 4c.");
  y0=giv4a.integ(x0,3.5);
  t.test_rel(y0,9.098214,1.0e-5,"rintp 4d.");

  // ---------------------------------------------------------------
  // Test akima

  gia.set_type(itp_akima);
  y0=gia.eval(x0,5,xa,ya);
  t.test_rel(y0,6.25,1.0e-7,"intp 3a.");
  y0=gia.deriv(x0,5,xa,ya);
  t.test_rel(y0,5.0,1.0e-7,"intp 3b.");
  y0=gia.deriv2(x0,5,xa,ya);
  t.test_rel(y0,2.0,1.0e-7,"intp 3c.");
  y0=gia.integ(x0,3.5,5,xa,ya);
  t.test_rel(y0,9.0833333,1.0e-7,"intp 3d.");

  // ---------------------------------------------------------------
  // Compare manual reversal in interp with new find()
  // in base_interp classes

  {
    interp<> oic(itp_linear);
    interp_linear<ubvector> cic;
    cic.set(5,rx,ry);
    double yc;

    for(double x0c=sqrt(2.0);x0c<=10.0;x0c+=2.0) {
      yc=cic.eval(x0c);
      //cout << yc << " " << oic.eval(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.eval(x0c,5,rx,ry),1.0e-12,"interp");
      yc=cic.deriv(x0c);
      //cout << yc << " " << oic.deriv(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv(x0c,5,rx,ry),1.0e-12,"deriv");
      yc=cic.deriv2(x0c);
      //cout << yc << " " << oic.deriv2(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv2(x0c,5,rx,ry),1.0e-12,"deriv2");
      yc=cic.integ(x0c,x0c-1.0);
      //cout << yc << " " << oic.integ(x0c,x0c-1.0,5,rx,ry) << endl;
      t.test_rel(yc,oic.integ(x0c,x0c-1.0,5,rx,ry),1.0e-12,"integ l");
    }
  }

  {
    interp<> oic(itp_cspline);
    interp_cspline<ubvector> cic;
    cic.set(5,rx,ry);
    double yc;

    for(double x0c=sqrt(2.0);x0c<=10.0;x0c+=2.0) {
      yc=cic.eval(x0c);
      //cout << yc << " " << oic.eval(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.eval(x0c,5,rx,ry),1.0e-12,"interp");
      yc=cic.deriv(x0c);
      //cout << yc << " " << oic.deriv(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv(x0c,5,rx,ry),1.0e-12,"deriv");
      yc=cic.deriv2(x0c);
      //cout << yc << " " << oic.deriv2(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv2(x0c,5,rx,ry),1.0e-12,"deriv2");
      yc=cic.integ(x0c,x0c-1.0);
      //cout << yc << " " << oic.integ(x0c,x0c-1.0,5,rx,ry) << endl;
      t.test_rel(yc,oic.integ(x0c,x0c-1.0,5,rx,ry),1.0e-12,"integ c");
    }
  }

  {
    interp<> oic(itp_cspline_peri);
    interp_cspline_peri<ubvector> cic;
    cic.set(5,rx,ry);
    double yc;

    for(double x0c=sqrt(2.0);x0c<=10.0;x0c+=2.0) {
      yc=cic.eval(x0c);
      //cout << yc << " " << oic.eval(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.eval(x0c,5,rx,ry),1.0e-12,"interp");
      yc=cic.deriv(x0c);
      //cout << yc << " " << oic.deriv(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv(x0c,5,rx,ry),1.0e-12,"deriv");
      yc=cic.deriv2(x0c);
      //cout << yc << " " << oic.deriv2(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv2(x0c,5,rx,ry),1.0e-12,"deriv2");
      yc=cic.integ(x0c,x0c-1.0);
      //cout << yc << " " << oic.integ(x0c,x0c-1.0,5,rx,ry) << endl;
      t.test_rel(yc,oic.integ(x0c,x0c-1.0,5,rx,ry),1.0e-12,"integ c");
    }
  }

  {
    interp<> oic(itp_akima);
    interp_akima<ubvector> cic;
    cic.set(5,rx,ry);
    double yc;

    for(double x0c=sqrt(2.0);x0c<=10.0;x0c+=2.0) {
      yc=cic.eval(x0c);
      //cout << yc << " " << oic.eval(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.eval(x0c,5,rx,ry),1.0e-12,"interp");
      yc=cic.deriv(x0c);
      //cout << yc << " " << oic.deriv(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv(x0c,5,rx,ry),1.0e-12,"deriv");
      yc=cic.deriv2(x0c);
      //cout << yc << " " << oic.deriv2(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv2(x0c,5,rx,ry),1.0e-12,"deriv2");
      yc=cic.integ(x0c,x0c-1.0);
      //cout << yc << " " << oic.integ(x0c,x0c-1.0,5,rx,ry) << endl;
      t.test_rel(yc,oic.integ(x0c,x0c-1.0,5,rx,ry),1.0e-12,"integ c");
    }
  }

  {
    interp<> oic(itp_akima_peri);
    interp_akima_peri<ubvector> cic;
    cic.set(5,rx,ry);
    double yc;

    for(double x0c=sqrt(2.0);x0c<=10.0;x0c+=2.0) {
      yc=cic.eval(x0c);
      //cout << yc << " " << oic.eval(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.eval(x0c,5,rx,ry),1.0e-12,"interp");
      yc=cic.deriv(x0c);
      //cout << yc << " " << oic.deriv(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv(x0c,5,rx,ry),1.0e-12,"deriv");
      yc=cic.deriv2(x0c);
      //cout << yc << " " << oic.deriv2(x0c,5,rx,ry) << " ";
      t.test_rel(yc,oic.deriv2(x0c,5,rx,ry),1.0e-12,"deriv2");
      yc=cic.integ(x0c,x0c-1.0);
      //cout << yc << " " << oic.integ(x0c,x0c-1.0,5,rx,ry) << endl;
      t.test_rel(yc,oic.integ(x0c,x0c-1.0,5,rx,ry),1.0e-12,"integ c");
    }
  }

  // ---------------------------------------------------------------
  // Test monotonicity-preserving interpolation

  // First, with a sparse dataset to show that monotonicity works. This
  // is challenging because the underlying function is not actually
  // monotonic.
  {
    static const size_t N=21;
    double a=2.01;

    bool debug=false;
    if (debug) cout.precision(4);
    
    // Test data 
    ubvector vx(N), vy(N);
    for(size_t i=0;i<N;i++) {
      vx[i]=((double)i);
      vy[i]=(sin(vx[i]/2.0)+vx[i]/a);
      if (i>0) t.test_gen(vy[i]>vy[i-1],"monotonic 1");
      if (debug) cout << i << " " << vx[i] << " " << vy[i] << endl;
    }
    if (debug) cout << endl;

    interp<> io(itp_monotonic);
    interp_vec<> iov(N,vx,vy,itp_monotonic);

    double last=-1.0;
    for(double x=0.0;x<=7.0;x+=0.02) {
      double exact=sin(x/2.0)+x/a;
      double deriv=cos(x/2.0)/2.0+1.0/a;
      double deriv2=-sin(x/2.0)/4.0;
      double interpa=io.eval(x,N,vx,vy);
      double interpb=iov.eval(x);
      double deriva=io.deriv(x,N,vx,vy);
      double derivb=iov.deriv(x);
      double deriv2a=io.deriv2(x,N,vx,vy);
      double deriv2b=iov.deriv2(x);
      if (debug) {
	cout << x << " " << exact << " " << interpa << " " 
	     << deriv << " " << deriva << " "
	     << deriv2 << " " << deriv2a << endl;
      }
      t.test_gen(interpa>last,"monotonic 2");
      t.test_rel(interpa,exact,4.0e-2,"interp");
      last=interpa;
    }

    if (debug) cout.precision(6);
  }

  // Second with a dense data set to show that the values are correct.
  // The underlying function is now monotonic.
  {
    static const size_t N=201;
    double a=1.99;

    bool debug=false;
    if (debug) cout.precision(4);
    
    // Test data 
    ubvector vx(N), vy(N);
    for(size_t i=0;i<N;i++) {
      vx[i]=((double)i)/10.0;
      vy[i]=(sin(vx[i]/2.0)+vx[i]/a);
      if (i>0) t.test_gen(vy[i]>vy[i-1],"monotonic 1");
      if (debug) cout << i << " " << vx[i] << " " << vy[i] << endl;
    }
    if (debug) cout << endl;

    interp<> io(itp_monotonic);
    interp_vec<> iov(N,vx,vy,itp_monotonic);

    double last=-1.0;
    for(double x=0.0;x<=7.0;x+=0.005) {
      double exact=sin(x/2.0)+x/a;
      double deriv=cos(x/2.0)/2.0+1.0/a;
      double deriv2=-sin(x/2.0)/4.0;
      double interpa=io.eval(x,N,vx,vy);
      double interpb=iov.eval(x);
      double deriva=io.deriv(x,N,vx,vy);
      double derivb=iov.deriv(x);
      double deriv2a=io.deriv2(x,N,vx,vy);
      double deriv2b=iov.deriv2(x);
      if (debug) {
	cout << x << " " << exact << " " << interpa << " " 
	     << deriv << " " << deriva << " "
	     << deriv2 << " " << deriv2a << endl;
      }
      t.test_gen(interpa>last,"monotonic 2");
      t.test_gen(interpa==interpb,"types");
      t.test_rel(interpa,exact,3.0e-4,"interp");
      t.test_gen(deriva==derivb,"types d");
      t.test_rel(deriva,deriv,1.0e-1,"deriv");
      t.test_gen(deriv2a==deriv2b,"types d2");
      t.test_rel(deriv2a,deriv2,2.0e2,"deriv2");
      for(double x2=0.0;x2<=7.0001;x2+=1.0) {
	if (fabs(x-x2)>1.0e-5) {
	  double intega=io.integ(x,x2,N,vx,vy);
	  double integb=iov.integ(x,x2);
	  double integ=(-x*x+x2*x2+4.0*a*cos(x/2.0)-4.0*a*cos(x2/2.0))/2.0/a;
	  t.test_gen(intega==integb,"types integ");
	  t.test_rel(intega,integ,2.0e-1,"integ");
	}
      }
      last=interpa;
    }

    if (debug) cout.precision(6);
  }

  // Third, same dense data but now decreasing
  {
    static const size_t N=201;
    double a=1.99;

    cout.precision(4);
    
    // Test data 
    ubvector vx(N), vy(N);
    for(size_t i=0;i<N;i++) {
      vx[i]=20.0-((double)i)/10.0;
      vy[i]=(sin(vx[i]/2.0)+vx[i]/a);
      if (i>0) t.test_gen(vy[i]<vy[i-1],"monotonic 1");
      //cout << i << " " << vx[i] << " " << vy[i] << endl;
    }
    //cout << endl;

    interp<> io(itp_monotonic);
    interp_vec<> iov(N,vx,vy,itp_monotonic);

    double last=-1.0;
    for(double x=0.0;x<=7.0;x+=0.005) {
      double exact=sin(x/2.0)+x/a;
      double deriv=cos(x/2.0)/2.0+1.0/a;
      double deriv2=-sin(x/2.0)/4.0;
      double interpa=io.eval(x,N,vx,vy);
      double interpb=iov.eval(x);
      double deriva=io.deriv(x,N,vx,vy);
      double derivb=iov.deriv(x);
      double deriv2a=io.deriv2(x,N,vx,vy);
      double deriv2b=iov.deriv2(x);
      //cout << x << " " << exact << " " << interpa << " " 
      //<< deriv << " " << deriva << " "
      //<< deriv2 << " " << deriv2a << endl;
      t.test_gen(interpa>last,"monotonic 2");
      t.test_gen(interpa==interpb,"types");
      t.test_rel(interpa,exact,3.0e-4,"interp");
      t.test_gen(deriva==derivb,"types d");
      t.test_rel(deriva,deriv,1.0e-1,"deriv");
      t.test_gen(deriv2a==deriv2b,"types d2");
      t.test_rel(deriv2a,deriv2,2.0e2,"deriv2");
      /*
	for(double x2=0.0;x2<=7.0001;x2+=1.0) {
	if (fabs(x-x2)>1.0e-5) {
	double intega=io.integ(x,x2,N,vx,vy);
	double integb=iov.integ(x,x2);
	double integ=(-x*x+x2*x2+4.0*a*cos(x/2.0)-4.0*a*cos(x2/2.0))/2.0/a;
	t.test_gen(intega==integb,"types integ pass2");
	t.test_rel(intega,integ,2.0e-1,"integ pass2");
	}
	}
      */
      last=interpa;
    }

    cout.precision(6);
  }

  // Four, non-monotonic data set
  {
    static const size_t N=201;

    cout.precision(3);
    
    // Test data 
    ubvector vx(N), vy(N);
    for(size_t i=0;i<N;i++) {
      vx[i]=((double)i)/10.0;
      vy[i]=sin(vx[i]/2.0);
      //cout << i << " " << vx[i] << " " << vy[i] << endl;
    }
    //cout << endl;

    interp<> io(itp_monotonic);
    interp_vec<> iov(N,vx,vy,itp_monotonic);

    double last=-1.0;
    for(double x=0.0;x<=7.0;x+=0.005) {
      double exact=sin(x/2.0);
      double deriv=cos(x/2.0)/2.0;
      double deriv2=-sin(x/2.0)/4.0;
      double interpa=io.eval(x,N,vx,vy);
      double interpb=iov.eval(x);
      double deriva=io.deriv(x,N,vx,vy);
      double derivb=iov.deriv(x);
      double deriv2a=io.deriv2(x,N,vx,vy);
      double deriv2b=iov.deriv2(x);
      //cout << x << " " << exact << " " << interpa << " " 
      //<< deriv << " " << deriva << " "
      //<< deriv2 << " " << deriv2a << endl;
      t.test_gen(interpa==interpb,"types");
      t.test_rel(interpa,exact,4.0e-3,"interp");
      t.test_gen(deriva==derivb,"types d");
      t.test_rel(deriva,deriv,7.0,"deriv");
      t.test_gen(deriv2a==deriv2b,"types d2");
      t.test_rel(deriv2a,deriv2,2.0e2,"deriv2");
      for(double x2=0.0;x2<=7.0001;x2+=1.0) {
	if (fabs(x-x2)>1.0e-5) {
	  double intega=io.integ(x,x2,N,vx,vy);
	  double integb=iov.integ(x,x2);
	  double integ=2.0*(cos(x/2.0)-cos(x2/2.0));
	  t.test_gen(intega==integb,"types integ pass3");
	  t.test_rel(intega,integ,5.0,"integ pass3");
	}
      }
      last=interpa;
    }

    cout.precision(6);
  }
  
  if (false) {
    // Try a data set which is flat over large regions,
    // which emphasizes the distinctions between the 
    // interpolation types

    double vx[13]={0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0};
    double vy[13]={0.0,0.1,0.5,0.9,1.0,1.0,0.9,0.5,0.1,0.0,0.0,0.1,0.5};

    interp_vec<double *> io1(13,vx,vy,itp_linear);
    interp_vec<double *> io2(13,vx,vy,itp_cspline);
    interp_vec<double *> io3(13,vx,vy,itp_cspline_peri);
    interp_vec<double *> io4(13,vx,vy,itp_akima);
    interp_vec<double *> io5(13,vx,vy,itp_akima_peri);
    interp_vec<double *> io6(13,vx,vy,itp_monotonic);

    cout.setf(ios::showpos);
    cout << "x i1 i2 i3 i4 i5 i6 ";
    cout << "d1 d2 d3 d4 d5 d6 ";
    cout << "e1 e2 e3 e4 e5 e6" << endl;
    for(double x=0.0;x<12.001;x+=0.2) {
      cout << x << " " << io1(x) << " " << io2(x) << " "
	   << io3(x) << " " << io4(x) << " " << io5(x) << " " 
	   << io6(x) << " ";
      cout << io1.deriv(x) << " " << io2.deriv(x) << " "
	   << io3.deriv(x) << " " << io4.deriv(x) << " " 
	   << io5.deriv(x) << " " << io6.deriv(x) << " ";
      cout << io1.deriv2(x) << " " << io2.deriv2(x) << " "
	   << io3.deriv2(x) << " " << io4.deriv2(x) << " " 
	   << io5.deriv2(x) << " " << io6.deriv2(x) << endl;
    }
    cout.unsetf(ios::showpos);

  }

  // ---------------------------------------------------------------
  // Test Steffen interpolation

  // First, with a sparse dataset to show that monotonicity works. This
  // is challenging because the underlying function is not actually
  // monotonic.
  {
    static const size_t N=21;
    double a=2.01;

    bool debug=false;
    if (debug) cout.precision(4);
    
    // Test data 
    ubvector vx(N), vy(N);
    for(size_t i=0;i<N;i++) {
      vx[i]=((double)i);
      vy[i]=(sin(vx[i]/2.0)+vx[i]/a);
      if (i>0) t.test_gen(vy[i]>vy[i-1],"steffen 1");
      if (debug) cout << i << " " << vx[i] << " " << vy[i] << endl;
    }
    if (debug) cout << endl;

    interp<> io(itp_steffen);
    interp_vec<> iov(N,vx,vy,itp_steffen);

    double last=-1.0;
    for(double x=0.0;x<=7.0;x+=0.02) {
      double exact=sin(x/2.0)+x/a;
      double deriv=cos(x/2.0)/2.0+1.0/a;
      double deriv2=-sin(x/2.0)/4.0;
      double interpa=io.eval(x,N,vx,vy);
      double interpb=iov.eval(x);
      double deriva=io.deriv(x,N,vx,vy);
      double derivb=iov.deriv(x);
      double deriv2a=io.deriv2(x,N,vx,vy);
      double deriv2b=iov.deriv2(x);
      double integ=-2.0*cos(x/2.0)+x*x/2.0/a+2.0;
      double intega=io.integ(0.0,x,N,vx,vy);
      double integb=iov.integ(0.0,x);
      if (debug) {
	cout.precision(4);
	cout << x << " " << exact << " " << interpa << " " 
	     << deriv << " " << deriva << " "
	     << deriv2 << " " << deriv2a << " "
	     << integ << " " << intega << " " << integb << endl;
	cout.precision(6);
      }
      t.test_gen(interpa>last,"steffen 2");
      t.test_gen(interpa==interpb,"types");
      t.test_rel(interpa,exact,4.0e-2,"interp");
      t.test_gen(deriva==derivb,"types d");
      t.test_rel(deriva,deriv,100.0,"deriv");
      t.test_gen(deriv2a==deriv2b,"types d2");
      t.test_rel(deriv2a,deriv2,50.0,"deriv2");
      t.test_gen(intega==integb,"types i");
      t.test_rel(intega,integ,4.0e-1,"integ2");
      last=interpa;
    }

    if (debug) cout.precision(6);
  }

  t.report();

  return 0;
}
