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
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/interp2_neigh.h>
#include <o2scl/interp2_planar.h>
#include <o2scl/rng_gsl.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double ft(double x, double y, double z) {
  return 3.0-2.0*x*x+7.0*y*z-5.0*z*x;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

    rng_gsl rg;

    // Construct the data
  ubvector x(8), y(8), dp(8);
  
  x[0]=1.04; y[0]=0.02; 
  x[1]=0.03; y[1]=1.01; 
  x[2]=0.81; y[2]=0.23; 
  x[3]=0.03; y[3]=0.83; 
  x[4]=0.03; y[4]=0.99; 
  x[5]=0.82; y[5]=0.84; 
  x[6]=0.03; y[6]=0.24; 
  x[7]=0.03; y[7]=1.02; 
  
  for(size_t i=0;i<8;i++) {
    dp[i]=1.0-pow(x[i]-0.5,2.0)-pow(y[i]-0.5,2.0);
  }

  // Reformat data into std::vector objects
  std::vector<ubvector> dat(3);
  dat[0]=x;
  dat[1]=y;
  dat[2]=dp;

  // Specify the data in the interpolation objects
  interp2_neigh<ubvector> i2n;
  interp2_planar<ubvector> i2p;
  interpm_idw<ubvector> imi;

  imi.set_data(2,1,8,dat);
  i2n.set_data(8,x,y,dp);
  i2p.set_data(8,x,y,dp);

  // Temporary storage
  double val, err;

  cout << "Interpolate at a point and compare the three methods:" << endl;
  ubvector point(2);
  point[0]=0.4;
  point[1]=0.5;
  imi.eval_err(point,val,err);
  cout << imi.eval(point) << " " << val << " " << err << " ";
  cout << i2n.eval(0.4,0.5) << " ";
  cout << i2p.eval(0.4,0.5) << endl;
  t.test_rel(imi.eval(point),i2n.eval(0.4,0.5),8.0e-2,"imi vs. i2n 1");
  t.test_rel(imi.eval(point),i2p.eval(0.4,0.5),4.0e-2,"imi vs. i2p 1");
  cout << endl;

  cout << "Interpolate at another point and compare the three methods:"
       << endl;
  point[0]=0.03;
  point[1]=1.0;
  imi.eval_err(point,val,err);
  cout << imi.eval(point) << " " << val << " " << err << " ";
  cout << i2n.eval(0.03,1.0) << " ";
  cout << i2p.eval(0.03,1.0) << endl;
  t.test_rel(imi.eval(point),i2n.eval(0.03,1.0),4.0e-2,"imi vs. i2n 2");
  t.test_rel(imi.eval(point),i2p.eval(0.03,1.0),1.0e-2,"imi vs. i2p 2");
  cout << endl;

  // Show how to swap a pointer instead
  std::vector<double> x2, y2, dp2;
  o2scl::vector_copy(x,x2);
  o2scl::vector_copy(y,y2);
  o2scl::vector_copy(dp,dp2);
  interpm_idw<double *> imi2;

  std::vector<double *> dat2(3);
  dat2[0]=&(x2[0]);
  dat2[1]=&(y2[0]);
  dat2[2]=&(dp2[0]);
  imi2.set_data(2,1,8,dat2);

  cout << "Same interpolation as above, but with pointers for storage:"
       << endl;
  imi.eval_err(point,val,err);
  cout << imi.eval(point) << " " << val << " " << err << endl;
  cout << endl;
  
  cout << "Show that interpolation gets better with more points." << endl;
  for(size_t N=10;N<1000000;N*=10) {
    // Create a random data set
    interpm_idw<std::vector<double> > imi3;
    std::vector<double> x3, y3, z3, f3;
    double scale=10.0;
    for(size_t i=0;i<N;i++) {
      x3.push_back(0.2+(2.0*rg.random()-1.0)/scale);
      y3.push_back(0.2+(2.0*rg.random()-1.0)/scale);
      z3.push_back(0.2+(2.0*rg.random()-1.0)/scale);
      f3.push_back(ft(x3[i],y3[i],z3[i]));
    }

    std::vector<double> p3={0.2,0.2,0.2};
    std::vector<std::vector<double> > dat3(4);
    std::vector<double> derivs(3), errs(3);
    double f;
    dat3[0]=x3;
    dat3[1]=y3;
    dat3[2]=z3;
    dat3[3]=f3;
    //imi3.verbose=1;
    imi3.set_data(3,1,N,dat3);
    imi3.eval_err(p3,val,err);
    cout.width(6);
    cout << N << " " << val << " " << err << " " << fabs(val-3.0) << endl;
  }
  cout << endl;

  cout << "Show that partial derivatives get better with more points."
       << endl;
  for(size_t N=10;N<1000000;N*=10) {
    // Create a random data set
    interpm_idw<std::vector<double> > imi3;
    std::vector<double> x3, y3, z3, f3;
    double scale=10.0;
    x3.push_back(0.2);
    y3.push_back(0.2);
    z3.push_back(0.2);
    f3.push_back(ft(0.2,0.2,0.2));
    for(size_t i=0;i<N;i++) {
      x3.push_back(0.2+(2.0*rg.random()-1.0)/scale);
      y3.push_back(0.2+(2.0*rg.random()-1.0)/scale);
      z3.push_back(0.2+(2.0*rg.random()-1.0)/scale);
      f3.push_back(ft(x3[i+1],y3[i+1],z3[i+1]));
    }

    std::vector<double> p3={0.2,0.2,0.2};
    std::vector<std::vector<double> > dat3(4);
    std::vector<double> derivs(3), errs(3);
    double f;
    dat3[0]=x3;
    dat3[1]=y3;
    dat3[2]=z3;
    dat3[3]=f3;
    imi3.verbose=1;
    imi3.set_data(3,1,N,dat3);
    cout.width(6);
    cout << N << endl;
    imi3.f_derivs_err(0,0,derivs,errs);
    cout << "\t" << -1.8 << " " << 1.4 << " " << 0.4 << endl;
    cout << "\t" << derivs[0] << " " << derivs[1] << " " << derivs[2] << endl;
    cout << "\t" << errs[0] << " " << errs[1] << " " << errs[2] << endl;
    cout << endl;
  }

  t.report();
  return 0;
}

