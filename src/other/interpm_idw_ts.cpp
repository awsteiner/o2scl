/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/interp2_neigh.h>
#include <o2scl/rng.h>
#include <o2scl/table.h>

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

  rng<> rg;
  
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
  std::vector<ubvector> dat_x(2);
  dat_x[0]=x;
  dat_x[1]=y;
  std::vector<ubvector> dat_y(1);
  dat_y[0]=dp;
  matrix_view_vec_vec_trans<ubvector> mv3_x(dat_x);
  matrix_view_vec_vec_trans<ubvector> mv3_y(dat_y);

  // Try a table representation
  table<> tab;
  tab.line_of_names("x y z");
  tab.line_of_data(2,vector<double>({1.04,0.02}));
  tab.line_of_data(2,vector<double>({0.03,1.01}));
  tab.line_of_data(2,vector<double>({0.81,0.23}));
  tab.line_of_data(2,vector<double>({0.03,0.83}));
  tab.line_of_data(2,vector<double>({0.03,0.99}));
  tab.line_of_data(2,vector<double>({0.82,0.84}));
  tab.line_of_data(2,vector<double>({0.03,0.24}));
  tab.line_of_data(2,vector<double>({0.03,1.02}));
  for(size_t i=0;i<8;i++) {
    tab.set("z",i,1.0-pow(x[i]-0.5,2.0)-pow(y[i]-0.5,2.0));
  }
  const_matrix_view_table<> cmvt_x(tab,{"x","y"});
  const_matrix_view_table<> cmvt_y(tab,{"z"});

  // Specify the data in the interpolation objects
  interp2_neigh<ubvector> i2n;
  interpm_idw<ubvector,
              matrix_view_vec_vec_trans<ubvector>,
              matrix_view_vec_vec_trans<ubvector>> imi;
  interpm_idw<> imi2;

  i2n.set_data(8,x,y,dp);
  imi.set_data(2,1,8,mv3_x,mv3_y);
  imi2.set_data(2,1,8,cmvt_x,cmvt_y);

  // Temporary storage
  double val, err;

  cout << "Interpolate at a point and compare the three methods:" << endl;
  ubvector point(2);
  point[0]=0.4;
  point[1]=0.5;
  double extrap;
  imi.eval_one_unc_tl(point,val,err,extrap);
  cout << imi.eval_one_tl(point) << " " << val << " " << err << " ";
  cout << i2n.eval(0.4,0.5) << " ";
  t.test_rel(imi.eval_one_tl(point),i2n.eval(0.4,0.5),8.0e-2,
             "imi vs. i2n 1");
  cout << endl;

  cout << "Test the extrapolation factor." << endl;
  for(point[0]=0.5;point[0]<20.0;point[0]*=2.0) {
    point[1]=point[0];
    imi.eval_one_unc_tl(point,val,err,extrap);
    cout << point[0] << " " << extrap << endl;
  }
  cout << endl;
  
  cout << "Interpolate at another point and compare the three methods:"
       << endl;
  point[0]=0.03;
  point[1]=1.0;
  imi.eval_one_unc_tl(point,val,err,extrap);
  cout << imi.eval_one_tl(point) << " " << val << " " << err << " ";
  cout << i2n.eval(0.03,1.0) << " ";
  t.test_rel(imi.eval_one_tl(point),
             i2n.eval(0.03,1.0),4.0e-2,"imi vs. i2n 2");
  cout << endl;

  cout << "Show that interpolation gets better with more points." << endl;
  for(size_t N=10;N<1000000;N*=10) {
    // Create a large data set
    interpm_idw<ubvector, matrix_view_vec_vec_trans<vector<double> >,
                matrix_view_vec_vec_trans<vector<double> > > imi3;
    
    std::vector<double> x3, y3, z3, f3;
    double scale=10.0;
    for(size_t i=0;i<N;i++) {
      double r1=fabs(sin(((double)(10*(i+1))))/1.04)+0.01;
      double r2=fabs(sin(((double)(10*(i+101))))/1.04)+0.01;
      double r3=fabs(sin(((double)(10*(i+201))))/1.04)+0.01;
      x3.push_back(0.2+(2.0*r1-1.0)/scale);
      y3.push_back(0.2+(2.0*r2-1.0)/scale);
      z3.push_back(0.2+(2.0*r3-1.0)/scale);
      f3.push_back(ft(x3[i],y3[i],z3[i]));
    }

    std::vector<double> p3={0.2,0.2,0.2};
    std::vector<std::vector<double> > dat3_x(3);
    std::vector<std::vector<double> > dat3_y(1);
    std::vector<double> derivs(3), errs(3);
    double f;
    dat3_x[0]=x3;
    dat3_x[1]=y3;
    dat3_x[2]=z3;
    dat3_y[0]=f3;
    //imi3.verbose=1;
    matrix_view_vec_vec_trans<vector<double> > mv3b_x(dat3_x);
    matrix_view_vec_vec_trans<vector<double> > mv3b_y(dat3_y);
    imi3.set_data(3,1,N,mv3b_x,mv3b_y);
    imi3.eval_one_unc_tl(p3,val,err,extrap);
    cout.width(6);
    cout << N << " " << val << " " << err << " " << fabs(val-3.0) << endl;
    if (N==1000000) {
      t.test_rel(val,3.0,10.0*err,"interp");
    }
  }
  cout << endl;

  cout << "Show that partial derivatives get better with more points."
       << endl;
  for(size_t N=10;N<1000000;N*=10) {
    // Create a random data set
    interpm_idw<ubvector, matrix_view_vec_vec_trans<vector<double> > ,
                matrix_view_vec_vec_trans<vector<double> > > imi3;
    std::vector<double> x3, y3, z3, f3;
    double scale=10.0;
    x3.push_back(0.2);
    y3.push_back(0.2);
    z3.push_back(0.2);
    f3.push_back(ft(0.2,0.2,0.2));
    for(size_t i=0;i<N;i++) {
      double r1=fabs(sin(((double)(10*(i+1))))/1.04)+0.01;
      double r2=fabs(sin(((double)(10*(i+101))))/1.04)+0.01;
      double r3=fabs(sin(((double)(10*(i+201))))/1.04)+0.01;
      x3.push_back(0.2+(2.0*r1-1.0)/scale);
      y3.push_back(0.2+(2.0*r2-1.0)/scale);
      z3.push_back(0.2+(2.0*r3-1.0)/scale);
      f3.push_back(ft(x3[i+1],y3[i+1],z3[i+1]));
    }

    std::vector<double> p3={0.2,0.2,0.2};
    std::vector<std::vector<double> > dat3_x(3);
    std::vector<std::vector<double> > dat3_y(1);
    std::vector<double> derivs(3), errs(3);
    double f;
    dat3_x[0]=x3;
    dat3_x[1]=y3;
    dat3_x[2]=z3;
    dat3_y[0]=f3;
    imi3.verbose=2;
    matrix_view_vec_vec_trans<vector<double> > mv3b_x(dat3_x);
    matrix_view_vec_vec_trans<vector<double> > mv3b_y(dat3_y);
    imi3.set_data(3,1,N,mv3b_x,mv3b_y);
    cout.width(6);
    imi3.derivs_err(0,0,derivs,errs);
    cout << N << endl;
    if (N==1000000) {
      t.test_rel(derivs[0],-1.8,errs[0]*10.0,"derivs 1");
      t.test_rel(derivs[1],1.4,errs[1]*10.0,"derivs 2");
      t.test_rel(derivs[2],0.4,errs[2]*10.0,"derivs 3");
    }
    cout << "\t" << -1.8 << " " << 1.4 << " " << 0.4 << endl;
    cout << "\t" << derivs[0] << " " << derivs[1] << " " << derivs[2] << endl;
    cout << "\t" << errs[0] << " " << errs[1] << " " << errs[2] << endl;
    cout << endl;
  }
    
  t.report();
  return 0;
}

