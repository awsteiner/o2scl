/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
#include <o2scl/interpm_python.h>
#include <o2scl/rng.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

double f(double x, double y) {
  return (sin(x*10)+2.0*tan(y))/5.0+0.14;
}

double f2(double x, double y) {
  double fv=f(x,y);
  return 2.0-fv*fv+fv;
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

#ifdef O2SCL_SET_PYTHON
  
  // Construct the data
  static const size_t N=100;
  ubvector x(N), y(N), dp(N), dp2(N);

  for(size_t i=0;i<N;i++) {
    x[i]=((double)i)/((double)N);
    y[i]=fabs(sin(1.0e8*i));
    dp[i]=f(x[i],y[i]);
    dp2[i]=f2(x[i],y[i]);
  }

  table<> tab;
  tab.line_of_names("x y dp dp2");
  for(size_t i=0;i<N;i++) {
    vector<double> line={x[i],y[i],dp[i],dp2[i]};
    tab.line_of_data(line.size(),line);
  }

  hdf_file hf2;
  hf2.open_or_create("interpm_python_data.o2");
  hdf_output(hf2,tab,"tab");
  hf2.close();

  table3d t3d;
  uniform_grid<double> ugx=uniform_grid_end<double>(0,1,99);
  uniform_grid<double> ugy=uniform_grid_end<double>(0,1,99);
  t3d.set_xy("x",ugx,"y",ugy);
  t3d.line_of_names("exact gp exact2 gp2 dnn dnn2");
  
  for(size_t i=0;i<100;i++) {
    for(size_t j=0;j<100;j++) {
      t3d.set(i,j,"exact",f(t3d.get_grid_x(i),t3d.get_grid_y(j)));
      t3d.set(i,j,"exact2",f2(t3d.get_grid_x(i),t3d.get_grid_y(j)));
    }
  }
    
  if (true) {

    // Sklearn Gaussian process, n_out=1
    
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

    interpm_python ip("interpm_sklearn_gp","verbose=3",1);
                      
    ip.set_data_tensor(2,1,N,tin,tout);
    
    std::vector<double> ex(2), ey(1), eyp(1);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval_std_vec(ex,ey);
    cout << ey[0] << endl;
    cout << f(0.5,0.5) << endl;
    t.test_rel(ey[0],f(ex[0],ex[1]),0.1,"sklearn gp 1");

    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval_unc_std_vec(ex,ey,eyp);
    cout << ey[0] << endl;
    cout << f(0.5,0.5) << endl;
    t.test_rel(ey[0],f(ex[0],ex[1]),0.1,"sklearn gp 2");
    t.test_abs(eyp[0],0.0,1.0e-5,"sklearn gp 3");
    
    cout << endl;
  }
    
  if (true) {

    // Sklearn Gaussian process, n_out=1, different kernel
    
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
    
    interpm_python ip("interpm_sklearn_gp",
                      ((std::string)"verbose=3,kernel=ConstantKernel")+
                      "(1.0,constant_value_bounds=\"fixed\")*"+
                      "RBF(1.0,length_scale_bounds=\"fixed\")",
                      1);
    ip.set_data_tensor(2,1,N,tin,tout);
    
    std::vector<double> ex(2), ey(1), eyp(1);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval_std_vec(ex,ey);
    cout << ey[0] << endl;
    cout << f(0.5,0.5) << endl;
    t.test_rel(ey[0],f(ex[0],ex[1]),0.1,"sklearn gp 1");

    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval_unc_std_vec(ex,ey,eyp);
    cout << ey[0] << endl;
    cout << f(0.5,0.5) << endl;
    t.test_rel(ey[0],f(ex[0],ex[1]),0.1,"sklearn gp 2");
    t.test_abs(eyp[0],0.0,1.0e-5,"sklearn gp 3");
    
    cout << endl;
  }
    
  if (true) {
    
    // Sklearn Gaussian process, n_out=2
    
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
    
    interpm_python ip("interpm_sklearn_gp","verbose=1",1);
                      
    ip.set_data_tensor(2,2,N,tin,tout);
    
    std::vector<double> ex(2), ey(2);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval_std_vec(ex,ey);
    cout << ey[0] << " " << ey[1] << endl;
    cout << f(0.5,0.5) << " " << f2(0.5,0.5) << endl;
    t.test_rel(ey[0],f(ex[0],ex[1]),0.1,"sklearn gp 1");
    t.test_rel(ey[1],f2(ex[0],ex[1]),0.1,"sklearn gp 2");
    
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        ex[0]=t3d.get_grid_x(i);
        ex[1]=t3d.get_grid_y(j);
        ip.eval_std_vec(ex,ey);
        t3d.set(i,j,"gp",ey[0]);
        t3d.set(i,j,"gp2",ey[1]);
      }
    }
    
    cout << endl;
  }
    
  if (true) {

    // TF dnn n_out=1
    
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

    // AWS, 5/8/24: I had to set verbose=0 because of the tensorflow
    // progress bars which cause problems in docker terminals
    
    interpm_python ip("interpm_tf_dnn",
                      ((std::string)"verbose=0,")+
                      "test_size=0.15,batch_size=10",0);
    ip.set_data_tensor(2,1,N,tin,tout);
    
    std::vector<double> ex(2), ey(1);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval_std_vec(ex,ey);
    cout << ey[0] << endl;
    cout << f(0.5,0.5) << endl;
    t.test_rel(ey[0],f(0.5,0.5),20.0,"tf_dnn 1");
    
    cout << endl;
  }
  
  if (false) {
    
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
    
    // AWS, 5/8/24: I had to set verbose=0 because of the tensorflow
    // progress bars which cause problems in docker terminals
    interpm_python ip("interpm_tf_dnn","verbose=0",0);
                      
    ip.set_data_tensor(2,2,N,tin,tout);
    
    std::vector<double> ex(2), ey(2);
    ex[0]=0.5;
    ex[1]=0.5;
    ip.eval_std_vec(ex,ey);
    cout << ey[0] << " " << ey[1] << endl;
    cout << f(0.5,0.5) << " " << f2(0.5,0.5) << endl;
    t.test_rel(ey[0],f(ex[0],ex[1]),20.0,"tf_dnn 2");
    t.test_rel(ey[1],f2(ex[0],ex[1]),20.0,"tf_dnn 3");

    if (false) {
      for(size_t i=0;i<100;i++) {
        std::cout << "i: " << i << std::endl;
        for(size_t j=0;j<100;j++) {
          ex[0]=t3d.get_grid_x(i);
          ex[1]=t3d.get_grid_y(j);
          ip.eval_std_vec(ex,ey);
          t3d.set(i,j,"dnn",ey[0]);
          t3d.set(i,j,"dnn2",ey[1]);
        }
      }
    }
    
    cout << endl;
  }
    
  hdf_file hf;
  hf.open_or_create("interpm_python_data.o2");
  hdf_output(hf,tab,"tab");
  hdf_output(hf,t3d,"t3d");
  hf.close();
  
#endif
    
  t.report();
  return 0;
}

