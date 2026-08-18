/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2026, Andrew W. Steiner
  
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
#include <o2scl/rng.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/interpm_rbf_libtorch.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

double f(double x, double y) {
  return (sin(x*10)+2.0*tan(y)+cos(x*y))/5.0+0.14;
}

double dfdx(double x, double y) {
  return (10.0*cos(x*10)-y*sin(x*y))/5.0;
}

double dfdy(double x, double y) {
  return (2.0/pow(cos(y),2.0)-x*sin(x*y))/5.0;
}

double d2fdx2(double x, double y) {
  return (-100.0*sin(x*10)-y*y*cos(x*y))/5.0;
}

double d2fdxy(double x, double y) {
  return (-sin(x*y)-x*y*cos(x*y))/5.0;
}

double d2fdy2(double x, double y) {
  return (4.0*sin(y)/pow(cos(y),3.0)-x*x*cos(x*y))/5.0;
}

double g(double x, double y) {
  return (sin(x*5.0)+2.0*tan(y)+sin(x*y))/5.0+0.24;
}

double dgdx(double x, double y) {
  return (5.0*cos(x*5.0)+y*cos(x*y))/5.0;
}

double dgdy(double x, double y) {
  return (2.0/cos(y)/cos(y)+x*cos(x*y))/5.0;
}

double d2gdx2(double x, double y) {
  return (-25.0*sin(5.0*x)-y*y*sin(x*y))/5.0;
}

double d2gdxy(double x, double y) {
  return (cos(x*y)-x*y*sin(x*y))/5.0;
}

double d2gdy2(double x, double y) {
  return (4.0*tan(y)/cos(y)/cos(y)-x*x*sin(x*y))/5.0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

#ifdef O2SCL_SET_LIBTORCH
  
  // Construct the data
  static const size_t N=400;
  
  ubvector x(N), y(N), dp(N), dp2(N);

  for(size_t i=0;i<N;i++) {
    x[i]=((double)i)/((double)N);
    y[i]=fabs(sin(1.0e8*i));
    dp[i]=f(x[i],y[i]);
    dp2[i]=g(x[i],y[i]);
  }

  table<> tab;
  tab.line_of_names("x y dp dp2");
  for(size_t i=0;i<N;i++) {
    vector<double> line={x[i],y[i],dp[i],dp2[i]};
    tab.line_of_data(line.size(),line);
  }

  hdf_file hf2;
  hf2.open_or_create("interpm_rbf_libtorch_data.o2");
  hdf_output(hf2,tab,"tab");
  hf2.close();

  if (true) {

    // Sklearn Gaussian process, n_out=1
    
    tensor2<> tin, tout;
    vector<size_t> in_size={N,2}, out_size={N,1};
    tin.resize(N,2);
    tout.resize(N,2);
    for(size_t j=0;j<N;j++) {
      tin.get(j,0)=x[j];
      tin.get(j,1)=y[j];
      tout.get(j,0)=dp[j];
      tout.get(j,1)=dp2[j];
    }

    interpm_rbf_libtorch<> ip;

    ip.test_size=0.2;
    ip.set_data_tensor(2,2,N,tin,tout);

    t.test_abs(ip.val_loss,0.0,2.0e-2,"val loss");
    t.test_abs(ip.train_loss,0.0,2.0e-2,"train loss");

    ubvector point(2);
    point[0]=0.5;
    point[1]=0.5;
    ubvector result(2);
    ip.eval(point,result);
    cout << "f: " << result[0] << " " << f(0.5,0.5) << endl;
    cout << "g: " << result[1] << " " << g(0.5,0.5) << endl;
    ip.deriv(point,result,0);
    cout << "dfdx: " << result[0] << " " << dfdx(0.5,0.5) << endl;
    cout << "dgdx: " << result[1] << " " << dgdx(0.5,0.5) << endl;
    ip.deriv(point,result,1);
    cout << "dfdy: " << result[0] << " " << dfdy(0.5,0.5) << endl;
    cout << "dgdy: " << result[1] << " " << dgdy(0.5,0.5) << endl;
    ip.deriv2(point,result,0,0);
    cout << "d2fdx2: " << result[0] << " " << d2fdx2(0.5,0.5) << endl;
    cout << "d2gdx2: " << result[1] << " " << d2gdx2(0.5,0.5) << endl;
    ip.deriv2(point,result,0,1);
    cout << "d2fdxy: " << result[0] << " " << d2fdxy(0.5,0.5) << endl;
    cout << "d2gdxy: " << result[1] << " " << d2gdxy(0.5,0.5) << endl;
    ip.deriv2(point,result,1,0);
    cout << "d2fdxy: " << result[0] << " " << d2fdxy(0.5,0.5) << endl;
    cout << "d2gdxy: " << result[1] << " " << d2gdxy(0.5,0.5) << endl;
    ip.deriv2(point,result,1,1);
    cout << "d2fdy2: " << result[0] << " " << d2fdy2(0.5,0.5) << endl;
    cout << "d2gdy2: " << result[1] << " " << d2gdy2(0.5,0.5) << endl;
    
    cout << endl;

    table3d t3d;
    uniform_grid<double> ugx=uniform_grid_end<double>(0,1,99);
    uniform_grid<double> ugy=uniform_grid_end<double>(0,1,99);
    t3d.set_xy("x",ugx,"y",ugy);
    t3d.line_of_names("f g lt_f lt_g dfdx dfdy dgdx dgdy");
    t3d.line_of_names("lt_dfdx lt_dfdy lt_dgdx lt_dgdy");
    t3d.line_of_names("d2fdx2 d2fdxy d2fdy2");
    t3d.line_of_names("d2gdx2 d2gdxy d2gdy2");
    t3d.line_of_names("lt_d2fdx2 lt_d2fdxy lt_d2fdy2");
    t3d.line_of_names("lt_d2gdx2 lt_d2gdxy lt_d2gdy2");
    
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        t3d.set(i,j,"f",f(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"g",g(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"dfdx",dfdx(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"dfdy",dfdy(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"dgdx",dgdx(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"dgdy",dgdy(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"d2fdx2",d2fdx2(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"d2fdxy",d2fdxy(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"d2fdy2",d2fdy2(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"d2gdx2",d2gdx2(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"d2gdxy",d2gdxy(t3d.get_grid_x(i),t3d.get_grid_y(j)));
        t3d.set(i,j,"d2gdy2",d2gdy2(t3d.get_grid_x(i),t3d.get_grid_y(j)));
      }
    }

    tin.resize(100*100,2);
    tout.resize(100*100,2);
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        tin.get(i*100+j,0)=t3d.get_grid_x(i);
        tin.get(i*100+j,1)=t3d.get_grid_y(j);
      }
    }
    ip.eval_list_tensor(tin,tout);
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        t3d.set(i,j,"lt_f",tout.get(i*100+j,0));
        t3d.set(i,j,"lt_g",tout.get(i*100+j,1));
      }
    }
    ip.deriv_list_tensor(tin,tout,0);
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        t3d.set(i,j,"lt_dfdx",tout.get(i*100+j,0));
        t3d.set(i,j,"lt_dgdx",tout.get(i*100+j,1));
      }
    }
    ip.deriv_list_tensor(tin,tout,1);
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        t3d.set(i,j,"lt_dfdy",tout.get(i*100+j,0));
        t3d.set(i,j,"lt_dgdy",tout.get(i*100+j,1));
      }
    }
    ip.deriv2_list_tensor(tin,tout,0,0);
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        t3d.set(i,j,"lt_d2fdx2",tout.get(i*100+j,0));
        t3d.set(i,j,"lt_d2gdx2",tout.get(i*100+j,1));
      }
    }
    ip.deriv2_list_tensor(tin,tout,0,1);
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        t3d.set(i,j,"lt_d2fdxy",tout.get(i*100+j,0));
        t3d.set(i,j,"lt_d2gdxy",tout.get(i*100+j,1));
      }
    }
    ip.deriv2_list_tensor(tin,tout,1,1);
    for(size_t i=0;i<100;i++) {
      for(size_t j=0;j<100;j++) {
        t3d.set(i,j,"lt_d2fdy2",tout.get(i*100+j,0));
        t3d.set(i,j,"lt_d2gdy2",tout.get(i*100+j,1));
      }
    }

    hdf_file hf;
    t3d.deriv_x("f","num_dfdx");
    t3d.deriv_y("f","num_dfdy");
    t3d.deriv_x("g","num_dgdx");
    t3d.deriv_y("g","num_dgdy");
    t3d.deriv_x("num_dfdx","num_d2fdx2");
    t3d.deriv_y("num_dfdx","num_d2fdxy");
    t3d.deriv_y("num_dfdy","num_d2fdy2");
    t3d.deriv_x("num_dgdx","num_d2gdx2");
    t3d.deriv_y("num_dgdx","num_d2gdxy");
    t3d.deriv_y("num_dgdy","num_d2gdy2");
    hf.open_or_create("interpm_rbf_libtorch.o2");
    hdf_output(hf,t3d,"t3d");
    hf.close();

    cout << "dfdx: " << t3d.interp(0.5,0.5,"lt_dfdx") << endl;
    cout << "dfdy: " << t3d.interp(0.5,0.5,"lt_dfdy") << endl;
    cout << "dgdx: " << t3d.interp(0.5,0.5,"lt_dgdx") << endl;
    cout << "dgdy: " << t3d.interp(0.5,0.5,"lt_dgdy") << endl;
    cout << "d2fdx2: " << t3d.interp(0.5,0.5,"lt_d2fdx2") << endl;
    cout << "d2fdxy: " << t3d.interp(0.5,0.5,"lt_d2fdxy") << endl;
    cout << "d2fdy2: " << t3d.interp(0.5,0.5,"lt_d2fdy2") << endl;
    cout << "d2gdx2: " << t3d.interp(0.5,0.5,"lt_d2gdx2") << endl;
    cout << "d2gdxy: " << t3d.interp(0.5,0.5,"lt_d2gdxy") << endl;
    cout << "d2gdy2: " << t3d.interp(0.5,0.5,"lt_d2gdy2") << endl;

  }
  
#endif
    
  t.report();
  return 0;
}

