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
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/interp2_neigh.h>
#include <o2scl/interp2_planar.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);
  
  ubvector x(8), y(8), dp(8);
  
  interp2_neigh<ubvector> i2n;
  interp2_planar<ubvector> i2p;
  interpm_idw<ubvector> imi;

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

  // Reformat data into a std::vector
  std::vector<ubvector> dat(3);
  dat[0]=x;
  dat[1]=y;
  dat[2]=dp;

  imi.set_data(2,8,dat);
  i2n.set_data(8,x,y,dp);
  i2p.set_data(8,x,y,dp);

  ubvector point(2);
  point[0]=0.4;
  point[1]=0.5;
  cout << imi.eval(point) << " ";
  cout << i2n.eval(0.4,0.5) << " ";
  cout << i2p.eval(0.4,0.5) << endl;
  t.test_rel(imi.eval(point),i2n.eval(0.4,0.5),8.0e-2,"imi vs. i2n 1");
  t.test_rel(imi.eval(point),i2p.eval(0.4,0.5),4.0e-2,"imi vs. i2p 1");

  point[0]=0.03;
  point[1]=1.0;
  cout << imi.eval(point) << " ";
  cout << i2n.eval(0.03,1.0) << " ";
  cout << i2p.eval(0.03,1.0) << endl;
  t.test_rel(imi.eval(point),i2n.eval(0.03,1.0),4.0e-2,"imi vs. i2n 2");
  t.test_rel(imi.eval(point),i2p.eval(0.03,1.0),1.0e-2,"imi vs. i2p 2");

  t.report();
  return 0;
}

