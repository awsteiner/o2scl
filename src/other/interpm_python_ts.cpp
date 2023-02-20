/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/interp2_neigh.h>
#include <o2scl/interp2_planar.h>
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

#ifdef O2SCL_PYTHON

  if (true) {
  
    tensor<> tin, tout;
    vector<size_t> in_size={8,2}, out_size={8,1};
    tin.resize(2,in_size);
    tout.resize(2,out_size);
    for(size_t j=0;j<8;j++) {
      vector<size_t> ix={j,0};
      tin.get(ix)=x[j];
      tout.get(ix)=dp[j];
      ix[1]++;
      tin.get(ix)=y[j];
    }
    interpm_python ip("o2sclpy","set_data_str","eval",2,8,1,
                      tin,tout,"verbose=2","interpm_sklearn_gpr",2);
  }
    
#endif
    
  t.report();
  return 0;
}

