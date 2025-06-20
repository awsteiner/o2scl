/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2025, Andrew W. Steiner
  
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
#include <o2scl/classify_python.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

int f(double x, double y) {
  return 3*((sin(x*10)+2.0*tan(y))/5.0+0.14);
}

int f2(double x, double y) {
  int fv=f(x,y);
  return 2-2*fv*fv+fv+2.0*x;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

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
  hf2.open_or_create("classify_python_data.o2");
  hdf_output(hf2,tab,"tab");
  hf2.close();

  table3d t3d;
  uniform_grid<double> ugx=uniform_grid_end<double>(0,1,99);
  uniform_grid<double> ugy=uniform_grid_end<double>(0,1,99);
  t3d.set_xy("x",ugx,"y",ugy);
  t3d.line_of_names("exact exact2 mlpc dtc gnb");
  
  for(size_t i=0;i<100;i++) {
    for(size_t j=0;j<100;j++) {
      t3d.set(i,j,"exact",f(t3d.get_grid_x(i),t3d.get_grid_y(j)));
      t3d.set(i,j,"exact2",f2(t3d.get_grid_x(i),t3d.get_grid_y(j)));
    }
  }
    
  if (true) {

    // Sklearn DTC, n_out=1
    
    tensor<> tin;
    tensor<int> tout;
    vector<size_t> in_size={N,2}, out_size={N,1};
    tin.resize(2,in_size);
    tout.resize(2,out_size);
    for(size_t j=0;j<N;j++) {
      vector<size_t> ix;
      ix={j,0};
      tin.get(ix)=x[j];
      tout.get(ix)=(int)dp[j];
      ix={j,1};
      tin.get(ix)=y[j];
    }

    classify_python<> cp("classify_sklearn_dtc",
                         ((std::string)"verbose=1"),1);
    classify_python<> cp2("classify_sklearn_dtc",
                          ((std::string)"verbose=1"),1);
                         
    cp.set_data_tensor(2,1,N,tin,tout);
    
    std::vector<double> ex(2);
    std::vector<int> ey(1);
    cout << "x y z_exact z_intp" << endl;
    for(double dx=0.1;dx<1.01;dx+=0.1) {
      ex[0]=dx;
      ex[1]=dx;
      cp.eval_std_vec(ex,ey);
      cout << ex[0] << " " << ex[1] << " ";
      cout << f(ex[0],ex[1]) << " ";
      cout << ey[0] << endl;
      t.test_gen(abs(ey[0]-f(ex[0],ex[1]))<=1,"sklearn dtc 1");
      if (dx<0.1001) {
        // Test save and load
        cp.save("classify_python_dtc.o2","cp");
        cp2.load("classify_python_dtc.o2","cp");
        cp2.eval_std_vec(ex,ey);
        cout << ex[0] << " " << ex[1] << " ";
        cout << f(ex[0],ex[1]) << " ";
        cout << ey[0] << endl;
        t.test_gen(abs(ey[0]-f(ex[0],ex[1]))<=1,"sklearn dtc 1b");
      }
    }

    for(size_t i=0;i<N;i++) {
      for(size_t j=0;j<N;j++) {
        ex[0]=t3d.get_grid_x(i);
        ex[1]=t3d.get_grid_y(j);
        cp.eval_std_vec(ex,ey);
        t3d.set(i,j,"dtc",ey[0]);
      }
    }
    
    cout << endl;
  }
    
  if (true) {

    // Sklearn GNB, n_out=1
    
    tensor<> tin;
    tensor<int> tout;
    vector<size_t> in_size={N,2}, out_size={N,1};
    tin.resize(2,in_size);
    tout.resize(2,out_size);
    for(size_t j=0;j<N;j++) {
      vector<size_t> ix;
      ix={j,0};
      tin.get(ix)=x[j];
      tout.get(ix)=(int)dp[j];
      ix={j,1};
      tin.get(ix)=y[j];
    }

    classify_python<> cp("classify_sklearn_gnb",
                         ((std::string)"verbose=1"),1);
    classify_python<> cp2("classify_sklearn_gnb",
                          ((std::string)"verbose=1"),1);
                         
    cp.set_data_tensor(2,1,N,tin,tout);
    
    std::vector<double> ex(2);
    std::vector<int> ey(1);
    cout << "x y z_exact z_intp" << endl;
    for(double dx=0.1;dx<1.01;dx+=0.1) {
      ex[0]=dx;
      ex[1]=dx;
      cp.eval_std_vec(ex,ey);
      cout << ex[0] << " " << ex[1] << " ";
      cout << f(ex[0],ex[1]) << " ";
      cout << ey[0] << endl;
      t.test_gen(abs(ey[0]-f(ex[0],ex[1]))<=1,"sklearn gnb 1");
      if (dx<0.1001) {
        // Test save and load
        cp.save("classify_python_gnb.o2","cp");
        cp2.load("classify_python_gnb.o2","cp");
        cp2.eval_std_vec(ex,ey);
        cout << ex[0] << " " << ex[1] << " ";
        cout << f(ex[0],ex[1]) << " ";
        cout << ey[0] << endl;
        t.test_gen(abs(ey[0]-f(ex[0],ex[1]))<=1,"sklearn gnb 1b");
      }
    }

    for(size_t i=0;i<N;i++) {
      for(size_t j=0;j<N;j++) {
        ex[0]=t3d.get_grid_x(i);
        ex[1]=t3d.get_grid_y(j);
        cp.eval_std_vec(ex,ey);
        t3d.set(i,j,"gnb",ey[0]);
      }
    }
    
    cout << endl;
  }
    
  if (true) {

    // Sklearn MLPC, n_out=1
    
    tensor<> tin;
    tensor<int> tout;
    vector<size_t> in_size={N,2}, out_size={N,1};
    tin.resize(2,in_size);
    tout.resize(2,out_size);
    for(size_t j=0;j<N;j++) {
      vector<size_t> ix;
      ix={j,0};
      tin.get(ix)=x[j];
      tout.get(ix)=(int)dp[j];
      ix={j,1};
      tin.get(ix)=y[j];
    }

    classify_python<> cp("classify_sklearn_mlpc",
                         ((std::string)"hlayers=[100,100],activation=")+
                         "relu,verbose=1,max_iter=2000",1);
    classify_python<> cp2("classify_sklearn_mlpc",
                          ((std::string)"hlayers=[100,100],activation=")+
                          "relu,verbose=1,max_iter=2000",1);
    cp.set_data_tensor(2,1,N,tin,tout);
    
    std::vector<double> ex(2);
    std::vector<int> ey(1);
    cout << "x y z_exact z_intp" << endl;
    for(double dx=0.1;dx<1.01;dx+=0.1) {
      ex[0]=dx;
      ex[1]=dx;
      cp.eval_std_vec(ex,ey);
      cout << ex[0] << " " << ex[1] << " ";
      cout << f(ex[0],ex[1]) << " ";
      cout << ey[0] << endl;
      t.test_gen(abs(ey[0]-f(ex[0],ex[1]))<=1,"sklearn mlpc 1");
      if (dx<0.1001) {
        // Test save and load
        cp.save("classify_python_mlpc.o2","cp");
        cp2.load("classify_python_mlpc.o2","cp");
        cp2.eval_std_vec(ex,ey);
        cout << ex[0] << " " << ex[1] << " ";
        cout << f(ex[0],ex[1]) << " ";
        cout << ey[0] << endl;
        t.test_gen(abs(ey[0]-f(ex[0],ex[1]))<=1,"sklearn dtc 1b");
      }
    }

    for(size_t i=0;i<N;i++) {
      for(size_t j=0;j<N;j++) {
        ex[0]=t3d.get_grid_x(i);
        ex[1]=t3d.get_grid_y(j);
        cp.eval_std_vec(ex,ey);
        t3d.set(i,j,"mlpc",ey[0]);
      }
    }
    
    cout << endl;
  }
    
  hdf_file hf;
  hf.open_or_create("classify_python_data.o2");
  hdf_output(hf,tab,"tab");
  hdf_output(hf,t3d,"t3d");
  hf.close();
  
#endif
    
  t.report();
  return 0;
}

