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
#include <o2scl/interpm_python.h>
#include <o2scl/rng.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

double ft(double x, double y, double z) {
  return 3.0-2.0*x*x+7.0*y*z-5.0*z*x;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  // Construct the data
  static const size_t N=50;
  ubvector x(N), y(N), dp(N), dp2(N);

  for(size_t i=0;i<N;i++) {
    x[i]=((double)i)/((double)N);
    y[i]=fabs(sin(1.0e8*i));
    dp[i]=1.0-pow(x[i]-0.5,2.0)-pow(y[i]-0.5,2.0);
    dp2[i]=2.0-dp[i]*dp[i]+dp[i];
  }

  table<> tab;
  tab.line_of_names("x y dp dp2");
  for(size_t i=0;i<N;i++) {
    vector<double> line={x[i],y[i],dp[i],dp2[i]};
    tab.line_of_data(line.size(),line);
  }

  hdf_file hf;
  hf.open_or_create("interpm_python_data.o2");
  hdf_output(hf,tab,"tab");
  hf.close();
  
#ifdef O2SCL_PYTHON

  if (true) {
    
    tensor<> tin, tout;
    vector<size_t> in_size={N,2}, out_size={N,1};
    tin.resize(2,in_size);
    tout.resize(2,out_size);
    for(size_t j=0;j<N;j++) {
      vector<size_t> ix;
      ix={j,0};
      tin.get(ix)=x[j];
      tout.get(ix)=dp[j];
      ix={j,1};
      tin.get(ix)=y[j];
    }
    
    interpm_python ip("o2sclpy","set_data_str","eval",2,N,1,
                      tin,tout,"verbose=2","interpm_sklearn_gp",2);
    
    std::vector<double> ex(2), ey(1);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval(ex,ey);
    cout << ey[0] << endl;
    cout << 1.0 << endl;
    t.test_rel(ey[0],1.0,0.5,"sklearn gp 1");
    
    cout << endl;
  }
    
  if (true) {
    
    tensor<> tin, tout;
    vector<size_t> in_size={N,2}, out_size={N,2};
    tin.resize(2,in_size);
    tout.resize(2,out_size);
    for(size_t j=0;j<N;j++) {
      vector<size_t> ix;
      ix={j,0};
      tin.get(ix)=x[j];
      tout.get(ix)=dp[j];
      ix={j,1};
      tin.get(ix)=y[j];
      tout.get(ix)=dp2[j];
    }
    
    interpm_python ip("o2sclpy","set_data_str","eval",2,N,2,
                      tin,tout,"verbose=2","interpm_sklearn_gp",2);
    
    std::vector<double> ex(2), ey(2);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval(ex,ey);
    cout << ey[0] << " " << ey[1] << endl;
    cout << 1.0 << " " << 2.0 << endl;
    t.test_rel(ey[0],1.0,0.5,"sklearn gp 2");
    t.test_rel(ey[1],2.0,0.5,"sklearn gp 3");
    
    cout << endl;
  }
    
  if (true) {
    
    tensor<> tin, tout;
    vector<size_t> in_size={N,2}, out_size={N,1};
    tin.resize(2,in_size);
    tout.resize(2,out_size);
    for(size_t j=0;j<N;j++) {
      vector<size_t> ix;
      ix={j,0};
      tin.get(ix)=x[j];
      tout.get(ix)=dp[j];
      ix={j,1};
      tin.get(ix)=y[j];
    }
    
    interpm_python ip("o2sclpy","set_data_str","eval",2,N,1,
                      tin,tout,"verbose=2","interpm_tf_dnn",2);
    
    std::vector<double> ex(2), ey(1);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval(ex,ey);
    cout << ey[0] << endl;
    cout << 1.0 << endl;
    t.test_rel(ey[0],1.0,0.5,"tf_dnn 1");
    
    cout << endl;
  }
    
  if (true) {
    
    tensor<> tin, tout;
    vector<size_t> in_size={N,2}, out_size={N,2};
    tin.resize(2,in_size);
    tout.resize(2,out_size);
    for(size_t j=0;j<N;j++) {
      vector<size_t> ix;
      ix={j,0};
      tin.get(ix)=x[j];
      tout.get(ix)=dp[j];
      ix={j,1};
      tin.get(ix)=y[j];
      tout.get(ix)=dp2[j];
    }
    
    interpm_python ip("o2sclpy","set_data_str","eval",2,N,2,
                      tin,tout,"verbose=2","interpm_tf_dnn",2);
    
    std::vector<double> ex(2), ey(2);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval(ex,ey);
    cout << ey[0] << " " << ey[1] << endl;
    cout << 1.0 << " " << 2.0 << endl;
    t.test_rel(ey[0],1.0,0.5,"tf_dnn 2");
    t.test_rel(ey[1],2.0,0.5,"tf_dnn 3");
    
    cout << endl;
  }
    
#endif
    
  t.report();
  return 0;
}

